"""Caller-specific filtering of biallelic target ALT records."""
import gzip
import math
from collections import Counter
FREEBAYES_QUAL_MIN = 1.0
FREEBAYES_NORMAL_VAF_MAX = 0.02

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def safe_float(x):
    try:
        if x is None or x == "" or x == ".":
            return None
        y = float(x)
        if math.isnan(y):
            return None
        return y
    except Exception:
        return None


def safe_int(x):
    try:
        if x is None or x == "" or x == ".":
            return None
        return int(float(x))
    except Exception:
        return None


def split_multi_value(value):
    if value is None or value == "" or value == ".":
        return []
    return str(value).split(",")


def get_alt_index_value(value, alt_index, cast="float"):
    vals = split_multi_value(value)
    if not vals:
        return None
    if len(vals) == 1:
        v = vals[0]
    elif alt_index < len(vals):
        v = vals[alt_index]
    else:
        return None

    if cast == "int":
        return safe_int(v)
    return safe_float(v)


def parse_info(info):
    d = {}
    if info in ("", "."):
        return d
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def parse_format(format_keys, sample_value):
    keys = format_keys.split(":") if format_keys else []
    vals = sample_value.split(":") if sample_value else []
    return dict(zip(keys, vals))


def find_sample_index_strict(sample_names, sample_id):
    hits = []
    for i, s in enumerate(sample_names):
        if s == sample_id:
            hits.append((i, s, "exact"))
        elif s.endswith("_" + sample_id):
            hits.append((i, s, "suffix_underscore"))

    if len(hits) == 1:
        return hits[0][0], hits[0][1], hits[0][2]

    return None, "", f"non_unique_or_missing:{hits}"


def variant_class(ref, alt):
    if alt is None or alt == "" or alt == ".":
        return "other"
    if alt.startswith("<") or alt == "*":
        return "other"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) != len(alt):
        return "short_indel"
    return "other"


def pass_filter(filter_value):
    return filter_value == "PASS"


def freebayes_normal_vaf(format_keys, samples, sample_names, normal_id, alt_index):
    normal_idx, matched_name, match_mode = find_sample_index_strict(sample_names, normal_id)

    if normal_idx is None:
        return None, None, None, "", "", "normal_sample_not_uniquely_matched"

    normal_fmt = parse_format(format_keys, samples[normal_idx])

    dp = safe_int(normal_fmt.get("DP"))
    alt_count = None

    if "AO" in normal_fmt:
        alt_count = get_alt_index_value(normal_fmt.get("AO"), alt_index, cast="int")

    if alt_count is None and "AD" in normal_fmt:
        ad_vals = split_multi_value(normal_fmt.get("AD"))
        if len(ad_vals) > alt_index + 1:
            alt_count = safe_int(ad_vals[alt_index + 1])

    if dp is None:
        return None, alt_count, None, matched_name, match_mode, "normal_dp_missing"

    if dp <= 0:
        return dp, alt_count if alt_count is not None else 0, 0.0, matched_name, match_mode, ""

    if alt_count is None:
        return dp, None, None, matched_name, match_mode, "normal_alt_missing"

    return dp, alt_count, alt_count / dp, matched_name, match_mode, ""


def freebayes_alt_support_pass(info, alt_index):
    required = ["SAF", "SAR", "RPL", "RPR"]
    values = {}
    missing = []

    for key in required:
        v = get_alt_index_value(info.get(key), alt_index, cast="float")
        values[key] = v
        if v is None:
            missing.append(key)

    if missing:
        return False, values, "missing_" + "_".join(missing)

    failed = [k for k, v in values.items() if v <= 0]
    if failed:
        return False, values, "zero_" + "_".join(failed)

    return True, values, ""


def parse_vcf_for_caller(pair_row, caller, vcf_path, forced_scope=None):
    pair_id = pair_row["pair_id"]
    lesion_id = pair_row["lesion_id"]
    normal_id = pair_row["normal_id"]

    sample_names = []
    counts = Counter()
    filter_summary = Counter()
    candidates = []

    input_records = 0
    input_alleles = 0
    freebayes_normal_match_name = ""
    freebayes_normal_match_mode = ""

    with open_text(vcf_path) as f:
        for line in f:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]

                if caller == "freebayes":
                    n_idx, n_name, n_mode = find_sample_index_strict(sample_names, normal_id)
                    l_idx, l_name, l_mode = find_sample_index_strict(sample_names, lesion_id)

                    if n_idx is None:
                        raise RuntimeError(
                            f"FreeBayes normal sample column cannot be uniquely matched: "
                            f"pair={pair_id}; normal_id={normal_id}; sample_names={sample_names}; detail={n_mode}"
                        )

                    if l_idx is None:
                        raise RuntimeError(
                            f"FreeBayes lesion sample column cannot be uniquely matched: "
                            f"pair={pair_id}; lesion_id={lesion_id}; sample_names={sample_names}; detail={l_mode}"
                        )

                    freebayes_normal_match_name = n_name
                    freebayes_normal_match_mode = n_mode

                continue

            if line.startswith("#"):
                continue

            input_records += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                filter_summary["malformed_record"] += 1
                continue

            chrom, pos, vid, ref, alt_field, qual_s, filter_value, info_s = parts[:8]
            fmt = parts[8] if len(parts) > 8 else ""
            samples = parts[9:] if len(parts) > 9 else []
            info = parse_info(info_s)
            qual = safe_float(qual_s)
            alts = alt_field.split(",") if alt_field not in ("", ".") else []

            for alt_index, alt in enumerate(alts):
                input_alleles += 1
                vclass = variant_class(ref, alt)

                if caller == "muse" and vclass != "snv":
                    filter_summary["muse_non_snv_record"] += 1
                    continue

                if forced_scope == "snv" and vclass != "snv":
                    filter_summary["scope_non_snv"] += 1
                    continue

                if forced_scope == "short_indel" and vclass != "short_indel":
                    filter_summary["scope_non_short_indel"] += 1
                    continue

                if vclass not in {"snv", "short_indel"}:
                    filter_summary["other_variant_class"] += 1
                    counts["other"] += 1
                    continue

                keep = True
                fail_reason = ""
                normal_dp = ""
                normal_alt = ""
                normal_vaf = ""
                fb_support = {"SAF": "", "SAR": "", "RPL": "", "RPR": ""}

                if caller in {"mutect2", "strelka2"}:
                    if not pass_filter(filter_value):
                        keep = False
                        fail_reason = "filter_not_PASS"

                elif caller == "muse":
                    if filter_value not in {"PASS", "Tier1", "Tier2", "Tier3", "Tier4"}:
                        keep = False
                        fail_reason = "filter_not_PASS_Tier1_Tier4"

                elif caller == "freebayes":
                    if qual is None or qual < FREEBAYES_QUAL_MIN:
                        keep = False
                        fail_reason = "qual_lt_1"

                    if keep:
                        support_ok, fb_support, support_reason = freebayes_alt_support_pass(info, alt_index)
                        if not support_ok:
                            keep = False
                            fail_reason = support_reason

                    normal_dp_calc, normal_alt_calc, normal_vaf_calc, matched_name, match_mode, normal_reason = freebayes_normal_vaf(
                        fmt, samples, sample_names, normal_id, alt_index
                    )

                    if matched_name:
                        freebayes_normal_match_name = matched_name
                    if match_mode:
                        freebayes_normal_match_mode = match_mode

                    normal_dp = normal_dp_calc
                    normal_alt = normal_alt_calc
                    normal_vaf = normal_vaf_calc

                    if normal_dp_calc is not None and normal_dp_calc <= 0:
                        filter_summary["audit_normal_dp_nonpositive"] += 1

                    if keep and normal_reason:
                        keep = False
                        fail_reason = normal_reason

                    if keep and normal_vaf is not None and normal_vaf > FREEBAYES_NORMAL_VAF_MAX:
                        keep = False
                        fail_reason = "normal_vaf_gt_2pct"

                else:
                    raise RuntimeError(f"Unknown caller: {caller}")

                if not keep:
                    filter_summary[fail_reason] += 1
                    continue

                counts[vclass] += 1
                counts["total"] += 1
                filter_summary["kept_" + vclass] += 1

                candidates.append({
                    "pair_id": pair_id,
                    "lesion_id": lesion_id,
                    "normal_id": normal_id,
                    "caller": caller,
                    "variant_class": vclass,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "variant_key": f"{chrom}:{pos}:{ref}:{alt}",
                    "qual": qual_s,
                    "filter": filter_value,
                    "normal_dp": normal_dp if caller == "freebayes" else "",
                    "normal_alt": normal_alt if caller == "freebayes" else "",
                    "normal_vaf": f"{normal_vaf:.6f}" if caller == "freebayes" and normal_vaf is not None else "",
                    "freebayes_normal_match_name": freebayes_normal_match_name if caller == "freebayes" else "",
                    "freebayes_normal_match_mode": freebayes_normal_match_mode if caller == "freebayes" else "",
                    "freebayes_SAF": fb_support.get("SAF", "") if caller == "freebayes" else "",
                    "freebayes_SAR": fb_support.get("SAR", "") if caller == "freebayes" else "",
                    "freebayes_RPL": fb_support.get("RPL", "") if caller == "freebayes" else "",
                    "freebayes_RPR": fb_support.get("RPR", "") if caller == "freebayes" else "",
                    "vcf": vcf_path,
                    "_parts": parts,
                    "_alt_index": alt_index,
                    "_sample_names": sample_names,
                    "_record_index": input_records,
                })

    metrics = {
        "pair_id": pair_id,
        "lesion_id": lesion_id,
        "normal_id": normal_id,
        "caller": caller,
        "variant_type_scope": forced_scope if forced_scope is not None else "snv_indel",
        "vcf": vcf_path,
        "input_records": input_records,
        "input_alleles": input_alleles,
        "snv_count": counts["snv"],
        "short_indel_count": counts["short_indel"],
        "other_count": counts["other"],
        "total_count": counts["total"],
        "freebayes_normal_match_name": freebayes_normal_match_name,
        "freebayes_normal_match_mode": freebayes_normal_match_mode,
    }

    for k, v in sorted(filter_summary.items()):
        metrics["filter_" + k] = v

    return metrics, candidates
