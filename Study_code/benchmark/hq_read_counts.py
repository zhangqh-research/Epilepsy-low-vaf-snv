"""HQ read recount used for the first spike-in truth set."""
MIN_MAPQ = 20
MIN_BASEQ = 20

def count_bases(bam, chrom, pos, ref, alt):
    ref_count = alt_count = other_count = 0

    for pileup_column in bam.pileup(
        chrom, pos - 1, pos,
        truncate=True,
        stepper="nofilter",
        max_depth=1000000,
        min_base_quality=0,
        min_mapping_quality=0,
        ignore_overlaps=False,
        ignore_orphans=False,
        compute_baq=False,
        redo_baq=False,
    ):
        if pileup_column.reference_pos != pos - 1:
            continue

        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment

            if (
                aln.is_unmapped
                or aln.is_secondary
                or aln.is_supplementary
                or aln.is_duplicate
                or aln.is_qcfail
                or aln.mapping_quality < MIN_MAPQ
                or pileup_read.is_del
                or pileup_read.is_refskip
            ):
                continue

            query_pos = pileup_read.query_position
            if query_pos is None or aln.query_sequence is None or aln.query_qualities is None:
                continue
            if query_pos >= len(aln.query_sequence) or query_pos >= len(aln.query_qualities):
                continue
            if aln.query_qualities[query_pos] < MIN_BASEQ:
                continue

            base = aln.query_sequence[query_pos].upper()
            if base == ref:
                ref_count += 1
            elif base == alt:
                alt_count += 1
            else:
                other_count += 1

    return ref_count + alt_count + other_count, ref_count, alt_count, other_count
