"""Pre-spike-in depth and target-ALT eligibility checks."""
import pysam
WORKER_LESION = None
WORKER_NORMAL = None
min_lesion_raw_dp = 80
min_lesion_hq_dp = 80
min_normal_raw_dp = 20
min_normal_hq_dp = 20

def init_worker(
    lesion_path,
    lesion_index,
    normal_path,
    normal_index,
):
    global WORKER_LESION
    global WORKER_NORMAL

    WORKER_LESION = pysam.AlignmentFile(
        lesion_path,
        "rb",
        index_filename=lesion_index,
    )

    WORKER_NORMAL = pysam.AlignmentFile(
        normal_path,
        "rb",
        index_filename=normal_index,
    )


def count_bases(bam, chrom, pos, ref, alt):
    raw_ref = raw_alt = raw_other = 0
    hq_ref = hq_alt = hq_other = 0

    for column in bam.pileup(
        chrom,
        pos - 1,
        pos,
        truncate=True,
        stepper="all",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=1000000,
        ignore_overlaps=True,
    ):
        if column.reference_pos != pos - 1:
            continue

        for pileupread in column.pileups:
            aln = pileupread.alignment

            if (
                aln.is_unmapped
                or aln.is_qcfail
                or aln.is_secondary
                or aln.is_supplementary
                or pileupread.is_del
                or pileupread.is_refskip
                or pileupread.query_position is None
                or aln.query_sequence is None
            ):
                continue

            qpos = pileupread.query_position
            base = aln.query_sequence[qpos].upper()

            if base not in {"A", "C", "G", "T"}:
                continue

            if base == ref:
                raw_ref += 1
            elif base == alt:
                raw_alt += 1
            else:
                raw_other += 1

            if (
                aln.is_duplicate
                or aln.mapping_quality < 20
                or aln.query_qualities is None
                or qpos >= len(aln.query_qualities)
                or aln.query_qualities[qpos] < 20
            ):
                continue

            if base == ref:
                hq_ref += 1
            elif base == alt:
                hq_alt += 1
            else:
                hq_other += 1

        break

    return {
        "raw_dp": raw_ref + raw_alt + raw_other,
        "raw_ref": raw_ref,
        "raw_alt": raw_alt,
        "raw_other": raw_other,
        "hq_dp": hq_ref + hq_alt + hq_other,
        "hq_ref": hq_ref,
        "hq_alt": hq_alt,
        "hq_other": hq_other,
    }


def screen(candidate):
    rank, chrom, pos, ref, alt = candidate

    lesion = count_bases(
        WORKER_LESION,
        chrom,
        pos,
        ref,
        alt,
    )

    normal = count_bases(
        WORKER_NORMAL,
        chrom,
        pos,
        ref,
        alt,
    )

    eligible = (
        lesion["raw_dp"] >= min_lesion_raw_dp
        and lesion["hq_dp"] >= min_lesion_hq_dp
        and normal["raw_dp"] >= min_normal_raw_dp
        and normal["hq_dp"] >= min_normal_hq_dp
        and lesion["raw_alt"] == 0
        and lesion["hq_alt"] == 0
        and normal["raw_alt"] == 0
        and normal["hq_alt"] == 0
    )

    return {
        "rank": rank,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "lesion": lesion,
        "normal": normal,
        "eligible": eligible,
    }
