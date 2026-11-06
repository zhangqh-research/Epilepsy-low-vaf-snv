"""Fragment evidence and the six BAM-supported criteria."""
import argparse
import bisect
import csv
import gzip
import statistics
from pathlib import Path
from collections import Counter, defaultdict
import pysam
MASK = 3844

def load(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'xt', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        for row in rows:
            w.writerow({k: 'NA' if v is None else v for k, v in row.items()})


def ratio(n, d):
    return n / d if d else None


def median(a):
    return statistics.median(a) if a else None


def observation(read, q):
    seq = read.query_sequence
    if seq is None or q is None or not 0 <= q < len(seq):
        return None
    qualities = read.query_qualities
    bq = qualities[q] if qualities is not None and qualities[q] != 255 else None
    mq = read.mapping_quality if read.mapping_quality < 255 else None
    direction = None
    if read.is_paired and read.is_read1 and not read.is_read2:
        direction = 'reverse' if read.is_reverse else 'forward'
    elif (read.is_paired and read.is_read2 and not read.is_read1 and not read.mate_is_unmapped
          and read.next_reference_id == read.reference_id):
        direction = 'reverse' if read.mate_is_reverse else 'forward'
    return (seq[q].upper(), bq, mq, min(q + 1, len(seq) - q), direction,
            int(read.is_read1), read.reference_start, read.flag)


def representative_key(o):
    return (o[1] if o[1] is not None else -1,
            o[2] if o[2] is not None else -1, o[5], -o[6], -o[7])


def summarize(fragments, ref, alt, prefix):
    bases = Counter()
    conflicts = no_acgt = overlap = ambiguous_direction = 0
    mq_values, bq_values, ends = [], [], []
    directions = Counter()
    joint = mq_missing = bq_missing = 0
    for obs in fragments.values():
        assert obs
        overlap += int(len(obs) > 1)
        acgt = [o for o in obs if o[0] in 'ACGT']
        if not acgt:
            no_acgt += 1
            continue
        if len({o[0] for o in acgt}) > 1:
            conflicts += 1
            continue
        rep = max(acgt, key=representative_key)
        base, bq, mq, end = rep[:4]
        bases[base] += 1
        if base != alt:
            continue
        if mq is None:
            mq_missing += 1
        else:
            mq_values.append(mq)
        if bq is None:
            bq_missing += 1
        else:
            bq_values.append(bq)
        joint += int(mq is not None and mq >= 40 and bq is not None and bq >= 30)
        ends.append(end)
        ds = {o[4] for o in acgt if o[4] is not None}
        direction = next(iter(ds)) if len(ds) == 1 else 'unknown'
        ambiguous_direction += int(len(ds) > 1)
        directions[direction] += 1
    a, dp = bases[alt], sum(bases.values())
    d = {'fragment_DP': dp, 'ALT_fragments': a, 'REF_fragments': bases[ref],
         'OTHER_fragments': dp - a - bases[ref], 'fragment_VAF': ratio(a, dp),
         'templates_seen': len(fragments), 'overlapping_templates': overlap,
         'conflicting_fragments': conflicts, 'nonACGT_only_fragments': no_acgt}
    for base in 'ACGT':
        d['base_' + base + '_fragments'] = bases[base]
    for name, vals, c1, c2, missing in [('MAPQ', mq_values, 20, 40, mq_missing),
                                      ('baseQ', bq_values, 20, 30, bq_missing)]:
        d['ALT_' + name + '_median'] = median(vals)
        d['ALT_' + name + '_missing_fragments'] = missing
        for label, n in [('low', sum(v < c1 for v in vals)),
                         ('mid', sum(c1 <= v < c2 for v in vals)),
                         ('high', sum(v >= c2 for v in vals))]:
            d['ALT_' + name + '_' + label + '_fragments'] = n
            d['ALT_' + name + '_' + label + '_fraction'] = ratio(n, a)
    d['ALT_joint_high_fragments'] = joint
    d['ALT_joint_high_fraction'] = ratio(joint, a)
    for direction in ['forward', 'reverse', 'unknown']:
        d['ALT_template_' + direction + '_fragments'] = directions[direction]
    d['ALT_direction_conflict_fragments'] = ambiguous_direction
    d['ALT_minor_direction_fraction'] = ratio(min(directions['forward'], directions['reverse']), a)
    d['ALT_end_distance_median'] = median(ends)
    d['ALT_end5_fragments'] = sum(e <= 5 for e in ends)
    d['ALT_end5_fraction'] = ratio(d['ALT_end5_fragments'], a)
    assert dp + conflicts + no_acgt == len(fragments)
    return {prefix + '_' + k: v for k, v in d.items()}


def query_cigar(read, positions):
    rp, qp = read.reference_start, 0
    for op, n in read.cigartuples or []:
        if op in (0, 7, 8):
            i = bisect.bisect_left(positions, rp)
            while i < len(positions) and positions[i] < rp + n:
                yield positions[i], qp + positions[i] - rp
                i += 1
            rp += n
            qp += n
        elif op in (1, 4):
            qp += n
        elif op in (2, 3):
            rp += n
        elif op not in (5, 6):
            raise ValueError(('unknown CIGAR', op))


def count_batch(bam, events, prefix, direct=False):
    bychr = defaultdict(lambda: defaultdict(list))
    for i, r in enumerate(events):
        bychr[r['chrom']][int(r['pos']) - 1].append(i)
    results = [None] * len(events)
    for chrom, lookup in bychr.items():
        blocks = []
        for p in sorted(lookup):
            if not direct and blocks and p - blocks[-1][-1] <= 500 and p - blocks[-1][0] <= 10000:
                blocks[-1].append(p)
            else:
                blocks.append([p])
        for positions in blocks:
            states = {p: defaultdict(list) for p in positions}
            for read in bam.fetch(chrom, positions[0], positions[-1] + 1):
                if read.flag & MASK or read.query_sequence is None:
                    continue
                key = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
                if direct:
                    mapped = ((rp, q) for q, rp in read.get_aligned_pairs(matches_only=False)
                              if rp in states and q is not None)
                else:
                    mapped = query_cigar(read, positions)
                for p, q in mapped:
                    o = observation(read, q)
                    if o is not None:
                        states[p][key].append(o)
            for p in positions:
                for i in lookup[p]:
                    results[i] = summarize(states[p], events[i]['ref'], events[i]['alt'], prefix)
    assert all(r is not None for r in results)
    return results


def label(row, plan):
    a, ldp, ndp = row['lesion_ALT_fragments'], row['lesion_fragment_DP'], row['normal_fragment_DP']
    minor = row['lesion_ALT_minor_direction_fraction']
    unknown = row['lesion_ALT_template_unknown_fragments']
    direction = minor >= .1 if minor is not None else None
    if direction is False and unknown:
        direction = None
    l_v, n_v = row['lesion_fragment_VAF'], row['normal_fragment_VAF']
    conditions = {
        'lesion_coverage': ldp >= plan['lesion_fragment_DP_min'],
        'ALT_support': a >= plan['lesion_ALT_fragment_min'],
        'joint_quality': row['lesion_ALT_joint_high_fraction'] >= .8 if a else None,
        'template_direction': direction,
        'read_position': row['lesion_ALT_end5_fraction'] <= .2 if a else None,
        'normal_background': (row['normal_ALT_fragments'] <= 1 and l_v is not None and l_v > n_v) if ndp >= 30 else None}
    # Boundaries use the locked thresholds; comparisons remain integer exact when possible.
    if a:
        conditions['joint_quality'] = 5 * row['lesion_ALT_joint_high_fragments'] >= 4 * a
        conditions['read_position'] = 5 * row['lesion_ALT_end5_fragments'] <= a
        if direction is not None:
            conditions['template_direction'] = (10 * min(row['lesion_ALT_template_forward_fragments'],
                                                        row['lesion_ALT_template_reverse_fragments']) >= a)
    if ndp >= 30 and ldp:
        conditions['normal_background'] = row['normal_ALT_fragments'] <= 1 and a * ndp > row['normal_ALT_fragments'] * ldp
    if all(v is True for v in conditions.values()):
        status = 'BAM-supported'
    elif any(v is False for v in conditions.values()):
        status = 'not_supported_by_prespecified_rule'
    else:
        status = 'not_fully_evaluable'
    for k, v in conditions.items():
        row['criterion_' + k] = None if v is None else int(v)
    row['normal_background_evaluable'] = int(ndp >= 30)
    row['BAM_supported'] = int(status == 'BAM-supported')
    row['support_status'] = status
    row['failed_criteria'] = ';'.join(k for k, v in conditions.items() if v is False)
    row['unevaluable_criteria'] = ';'.join(k for k, v in conditions.items() if v is None)
    return row


def validate(row):
    for p in ['lesion', 'normal']:
        dp, a = row[p + '_fragment_DP'], row[p + '_ALT_fragments']
        assert 0 <= a <= dp
        assert dp == row[p + '_REF_fragments'] + a + row[p + '_OTHER_fragments']
        assert dp == sum(row[p + '_base_' + b + '_fragments'] for b in 'ACGT')
        assert dp + row[p + '_conflicting_fragments'] + row[p + '_nonACGT_only_fragments'] == row[p + '_templates_seen']
        v = row[p + '_fragment_VAF']
        assert (dp == 0 and v is None) or (dp and abs(v - a / dp) < 1e-12)
        for q in ['MAPQ', 'baseQ']:
            assert a == sum(row[p + '_ALT_' + q + '_' + c + '_fragments'] for c in ['low', 'mid', 'high', 'missing'])
        assert a == sum(row[p + '_ALT_template_' + c + '_fragments'] for c in ['forward', 'reverse', 'unknown'])
        assert 0 <= row[p + '_ALT_joint_high_fragments'] <= min(row[p + '_ALT_MAPQ_high_fragments'], row[p + '_ALT_baseQ_high_fragments'])
        assert 0 <= row[p + '_ALT_end5_fragments'] <= a


def main():
    parser = argparse.ArgumentParser(description='Count fragment evidence at candidate SNVs.')
    parser.add_argument('--events', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    events, pairs = load(args.events), load(args.manifest)
    assert len(pairs) == len({r['pair_id'] for r in pairs}) == 83
    bypair = defaultdict(list)
    for e in events:
        bypair[e['pair_id']].append(e)
    assert set(bypair) <= {r['pair_id'] for r in pairs}
    rows = []
    for pair in pairs:
        subset = bypair[pair['pair_id']]
        if not subset:
            continue
        with pysam.AlignmentFile(pair['lesion_bam'], 'rb', index_filename=pair['lesion_bai']) as bam:
            lesion = count_batch(bam, subset, 'lesion')
        with pysam.AlignmentFile(pair['normal_bam'], 'rb', index_filename=pair['normal_bai']) as bam:
            normal = count_batch(bam, subset, 'normal')
        for e, l, n in zip(subset, lesion, normal):
            row = label({**e, **l, **n}, {'lesion_fragment_DP_min': 100, 'lesion_ALT_fragment_min': 5})
            validate(row)
            rows.append(row)
    rows.sort(key=lambda r: r['event_id'])
    assert len(rows) == len(events)
    write(args.output, rows)

if __name__ == '__main__':
    main()
