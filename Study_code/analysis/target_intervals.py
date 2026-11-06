"""Exact WES target containment after allele normalization."""
import bisect
from collections import defaultdict

class TargetBed:
    def __init__(self, path):
        intervals = defaultdict(list)
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                chrom, start, end = line.split()[:3]
                start, end = int(start), int(end)
                assert 0 <= start < end
                intervals[chrom].append((start, end))
        self.merged = {}
        self.starts = {}
        for chrom, values in intervals.items():
            out = []
            for start, end in sorted(values):
                if out and start <= out[-1][1]:
                    out[-1] = (out[-1][0], max(end, out[-1][1]))
                else:
                    out.append((start, end))
            self.merged[chrom] = out
            self.starts[chrom] = [a for a, b in out]

    def contains(self, chrom, start, end):
        if start < 0 or end <= start:
            return False
        i = bisect.bisect_right(self.starts.get(chrom, []), start) - 1
        return i >= 0 and end <= self.merged[chrom][i][1]

    def classify(self, chrom, pos, ref, alt):
        start = int(pos) - 1
        # Remove all shared sequence, so VCF padding is not treated as altered DNA.
        while ref and alt and ref[0] == alt[0]:
            start += 1
            ref, alt = ref[1:], alt[1:]
        while ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        if not ref and not alt:
            raise ValueError('No-change allele')
        if not ref:
            left, right, kind = start - 1, start + 1, 'insertion_two_flanks'
        else:
            left, right = start, start + len(ref)
            kind = 'snv_base' if len(ref) == len(alt) == 1 else 'changed_reference_span'
        return self.contains(chrom, left, right), left, right, kind
