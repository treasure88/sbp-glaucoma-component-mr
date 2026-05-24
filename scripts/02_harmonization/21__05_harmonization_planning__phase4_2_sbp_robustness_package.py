#!/usr/bin/env python3
import csv
import gzip
import math
import os
import random
import time
from datetime import datetime
from html import escape

OUTDIR = "../../22_sbp_robustness"
FIGDIR = os.path.join(OUTDIR, "figures")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

SBP_INPUTS = {
    "GBS_nonIOPcomponent": "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz",
    "GBS_IOPcomponent": "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz",
}

METHOD_OUT = os.path.join(OUTDIR, "phase4_2_sbp_method_consistency.tsv")
RADIAL_SUMMARY_OUT = os.path.join(OUTDIR, "phase4_2_sbp_radial_outlier_summary.tsv")
RADIAL_SNPS_OUT = os.path.join(OUTDIR, "phase4_2_sbp_radial_outlier_snps.tsv")
OUTLIER_CORRECTED_OUT = os.path.join(OUTDIR, "phase4_2_sbp_outlier_corrected_results.tsv")
CHR_OUT = os.path.join(OUTDIR, "phase4_2_sbp_leave_one_chromosome_out.tsv")
INFLUENCE_OUT = os.path.join(OUTDIR, "phase4_2_sbp_single_snp_influence_summary.tsv")
CONTRAST_OUT = os.path.join(OUTDIR, "phase4_2_sbp_component_contrast_by_method.tsv")
SUMMARY_OUT = os.path.join(OUTDIR, "phase4_2_sbp_robustness_summary.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase4_2_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase4_2_runtime_log.tsv")

FIG_METHOD = os.path.join(FIGDIR, "phase4_2_sbp_method_consistency.svg")
FIG_CHR = os.path.join(FIGDIR, "phase4_2_sbp_leave_one_chromosome_out.svg")

R_VALUES = [-0.25, 0.0, 0.25, 0.5]
BOOT_N = 500
RANDOM_SEED = 20260517

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def parse_float(x):
    if x is None:
        return None
    x = str(x).strip()
    if x in ("", ".", "NA"):
        return None
    try:
        return float(x)
    except Exception:
        return None

def fmt(x, digits=8):
    if x is None:
        return "NA"
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    try:
        if math.isnan(x) or math.isinf(x):
            return "NA"
    except Exception:
        pass
    return f"{x:.{digits}g}"

def normal_p(z):
    if z is None:
        return None
    return math.erfc(abs(z) / math.sqrt(2.0))

def direction(beta):
    if beta is None:
        return "NA"
    if beta > 0:
        return "positive"
    if beta < 0:
        return "negative"
    return "zero"

def choose_col(header, candidates):
    hset = set(header)
    for c in candidates:
        if c in hset:
            return c
    return None

def read_mr_input(path, outcome_suffix):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames

        col_snp = choose_col(header, ["SNP", "rsid", "variant"])
        col_chr = choose_col(header, ["chr", "chromosome", "CHR"])
        col_bx = choose_col(header, ["beta_exposure", "beta.exposure", "exposure_beta", "beta_exp"])
        col_sx = choose_col(header, ["se_exposure", "se.exposure", "exposure_se", "se_exp"])
        col_by = choose_col(header, ["beta_outcome_harmonized", "beta_outcome", "beta.outcome", "outcome_beta", "beta_out"])
        col_sy = choose_col(header, ["se_outcome", "se.outcome", "outcome_se", "se_out"])
        col_eaf = choose_col(header, ["eaf_exposure", "eaf", "effect_allele_frequency"])
        col_f = choose_col(header, ["F_stat", "f_stat", "F", "f_statistic"])
        col_action = choose_col(header, ["harmonization_action", "action"])

        required = {
            "SNP": col_snp,
            "chr": col_chr,
            "beta_exposure": col_bx,
            "se_exposure": col_sx,
            "beta_outcome": col_by,
            "se_outcome": col_sy,
        }

        missing = [k for k, v in required.items() if v is None]
        if missing:
            raise RuntimeError(f"Missing required columns in {path}: {missing}; header={header}")

        for r in reader:
            snp = r.get(col_snp, "")
            chrom = r.get(col_chr, "")
            bx = parse_float(r.get(col_bx))
            sx = parse_float(r.get(col_sx))
            by = parse_float(r.get(col_by))
            sy = parse_float(r.get(col_sy))

            if not snp or bx is None or sx is None or by is None or sy is None:
                continue
            if bx == 0 or sy <= 0:
                continue

            rows.append({
                "outcome_suffix": outcome_suffix,
                "SNP": snp,
                "chr": chrom,
                "beta_exposure": bx,
                "se_exposure": sx,
                "beta_outcome": by,
                "se_outcome": sy,
                "eaf_exposure": r.get(col_eaf, "NA") if col_eaf else "NA",
                "F_stat": r.get(col_f, "NA") if col_f else "NA",
                "harmonization_action": r.get(col_action, "NA") if col_action else "NA",
            })

    return rows

def ivw_no_intercept(rows, random_effects=False):
    sum_wxx = 0.0
    sum_wxy = 0.0

    for r in rows:
        bx = r["beta_exposure"]
        by = r["beta_outcome"]
        sy = r["se_outcome"]
        w = 1.0 / (sy * sy)
        sum_wxx += w * bx * bx
        sum_wxy += w * bx * by

    if sum_wxx <= 0:
        return None

    beta = sum_wxy / sum_wxx

    q = 0.0
    for r in rows:
        bx = r["beta_exposure"]
        by = r["beta_outcome"]
        sy = r["se_outcome"]
        w = 1.0 / (sy * sy)
        resid = by - beta * bx
        q += w * resid * resid

    df = max(len(rows) - 1, 1)
    phi = max(q / df, 1.0) if random_effects else 1.0
    se = math.sqrt(phi / sum_wxx)
    z = beta / se if se > 0 else None
    p = normal_p(z)

    return {
        "beta": beta,
        "se": se,
        "pval": p,
        "Q": q,
        "Q_df": df,
        "phi": phi,
        "n": len(rows),
    }

def weighted_quantile(values, weights, q):
    pairs = sorted(zip(values, weights), key=lambda x: x[0])
    total = sum(w for _, w in pairs)
    if total <= 0:
        return None
    target = q * total
    c = 0.0
    for v, w in pairs:
        c += w
        if c >= target:
            return v
    return pairs[-1][0]

def weighted_median(rows):
    ratios = []
    weights = []
    for r in rows:
        bx = r["beta_exposure"]
        by = r["beta_outcome"]
        sy = r["se_outcome"]
        if bx == 0 or sy <= 0:
            continue
        ratios.append(by / bx)
        weights.append((bx * bx) / (sy * sy))

    beta = weighted_quantile(ratios, weights, 0.5)
    se = bootstrap_weighted_median(rows, BOOT_N)
    z = beta / se if beta is not None and se is not None and se > 0 else None
    p = normal_p(z)

    return {
        "beta": beta,
        "se": se,
        "pval": p,
        "n": len(rows),
    }

def bootstrap_weighted_median(rows, n_boot):
    if len(rows) < 5:
        return None
    estimates = []
    n = len(rows)

    for _ in range(n_boot):
        sample = [rows[random.randrange(n)] for __ in range(n)]
        ratios = []
        weights = []
        for r in sample:
            bx = r["beta_exposure"]
            by = r["beta_outcome"]
            sy = r["se_outcome"]
            if bx == 0 or sy <= 0:
                continue
            ratios.append(by / bx)
            weights.append((bx * bx) / (sy * sy))
        est = weighted_quantile(ratios, weights, 0.5)
        if est is not None:
            estimates.append(est)

    if len(estimates) < 10:
        return None

    mean = sum(estimates) / len(estimates)
    var = sum((x - mean) ** 2 for x in estimates) / max(len(estimates) - 1, 1)
    return math.sqrt(var)

def weighted_mode(rows, bins=120):
    ratios = []
    weights = []
    for r in rows:
        bx = r["beta_exposure"]
        by = r["beta_outcome"]
        sy = r["se_outcome"]
        if bx == 0 or sy <= 0:
            continue
        ratios.append(by / bx)
        weights.append((bx * bx) / (sy * sy))

    if not ratios:
        return {"beta": None, "se": None, "pval": None, "n": len(rows)}

    lo = weighted_quantile(ratios, weights, 0.02)
    hi = weighted_quantile(ratios, weights, 0.98)
    if lo is None or hi is None or hi <= lo:
        return {"beta": None, "se": None, "pval": None, "n": len(rows)}

    bin_w = [0.0] * bins
    bin_vals = [[] for _ in range(bins)]
    width = (hi - lo) / bins

    for v, w in zip(ratios, weights):
        if v < lo or v > hi:
            continue
        idx = int((v - lo) / width)
        if idx >= bins:
            idx = bins - 1
        bin_w[idx] += w
        bin_vals[idx].append((v, w))

    best = max(range(bins), key=lambda i: bin_w[i])
    vals = bin_vals[best]
    if not vals:
        return {"beta": None, "se": None, "pval": None, "n": len(rows)}

    beta = sum(v * w for v, w in vals) / sum(w for _, w in vals)
    return {"beta": beta, "se": None, "pval": None, "n": len(rows)}

def weighted_regression_with_intercept(rows):
    sw = sx = sy = sxx = sxy = 0.0

    for r in rows:
        x = r["beta_exposure"]
        y = r["beta_outcome"]
        w = 1.0 / (r["se_outcome"] ** 2)
        sw += w
        sx += w * x
        sy += w * y
        sxx += w * x * x
        sxy += w * x * y

    det = sw * sxx - sx * sx
    if det <= 0:
        return None

    intercept = (sxx * sy - sx * sxy) / det
    slope = (sw * sxy - sx * sy) / det

    q = 0.0
    for r in rows:
        x = r["beta_exposure"]
        y = r["beta_outcome"]
        w = 1.0 / (r["se_outcome"] ** 2)
        resid = y - intercept - slope * x
        q += w * resid * resid

    df = max(len(rows) - 2, 1)
    phi = max(q / df, 1.0)

    var_intercept = phi * sxx / det
    var_slope = phi * sw / det

    se_intercept = math.sqrt(var_intercept)
    se_slope = math.sqrt(var_slope)

    p_slope = normal_p(slope / se_slope)
    p_intercept = normal_p(intercept / se_intercept)

    return {
        "beta": slope,
        "se": se_slope,
        "pval": p_slope,
        "egger_intercept": intercept,
        "egger_intercept_se": se_intercept,
        "egger_intercept_pval": p_intercept,
        "Q": q,
        "Q_df": df,
        "phi": phi,
        "n": len(rows),
    }

def radial_outlier_check(rows, ivw_result, threshold=3.0):
    beta = ivw_result["beta"]
    phi = ivw_result.get("phi", 1.0)
    phi = max(phi, 1.0)

    outliers = []
    stats = []

    for r in rows:
        bx = r["beta_exposure"]
        by = r["beta_outcome"]
        sy = r["se_outcome"]
        resid = by - beta * bx
        std_resid = resid / (sy * math.sqrt(phi))
        abs_std = abs(std_resid)

        item = dict(r)
        item["radial_residual"] = resid
        item["radial_std_residual"] = std_resid
        item["radial_abs_std_residual"] = abs_std
        item["radial_outlier_flag"] = "YES" if abs_std > threshold else "NO"
        stats.append(item)

        if abs_std > threshold:
            outliers.append(item)

    return outliers, stats

def leave_one_chromosome_out(rows):
    chroms = sorted(set(r["chr"] for r in rows), key=lambda x: (str(x).isdigit() == False, str(x)))
    out = []
    full = ivw_no_intercept(rows, random_effects=True)

    for chrom in chroms:
        subset = [r for r in rows if r["chr"] != chrom]
        if len(subset) < 3:
            continue
        est = ivw_no_intercept(subset, random_effects=True)
        if est is None:
            continue
        out.append({
            "left_out_chr": chrom,
            "n_remaining": est["n"],
            "beta": est["beta"],
            "se": est["se"],
            "pval": est["pval"],
            "direction": direction(est["beta"]),
            "delta_from_full_beta": est["beta"] - full["beta"],
            "full_beta": full["beta"],
        })

    return out

def single_snp_influence(rows):
    full = ivw_no_intercept(rows, random_effects=True)
    out = []

    for r in rows:
        subset = [x for x in rows if x["SNP"] != r["SNP"]]
        if len(subset) < 3:
            continue
        est = ivw_no_intercept(subset, random_effects=True)
        if est is None:
            continue
        out.append({
            "left_out_SNP": r["SNP"],
            "chr": r["chr"],
            "beta_leave_one_out": est["beta"],
            "se_leave_one_out": est["se"],
            "pval_leave_one_out": est["pval"],
            "delta_from_full_beta": est["beta"] - full["beta"],
            "abs_delta_from_full_beta": abs(est["beta"] - full["beta"]),
            "direction_leave_one_out": direction(est["beta"]),
            "full_beta": full["beta"],
        })

    out.sort(key=lambda x: x["abs_delta_from_full_beta"], reverse=True)
    return out

def contrast(beta_iop, se_iop, beta_non, se_non, r):
    beta_diff = beta_iop - beta_non
    var = se_iop ** 2 + se_non ** 2 - 2 * r * se_iop * se_non
    if var <= 0:
        return None
    se = math.sqrt(var)
    z = beta_diff / se
    p = normal_p(z)
    lo = beta_diff - 1.96 * se
    hi = beta_diff + 1.96 * se
    return beta_diff, se, z, p, lo, hi

def make_method_results(all_data):
    method_rows = []

    for outcome, rows in all_data.items():
        ivw_fixed = ivw_no_intercept(rows, random_effects=False)
        ivw_random = ivw_no_intercept(rows, random_effects=True)
        wmed = weighted_median(rows)
        wmode = weighted_mode(rows)
        egger = weighted_regression_with_intercept(rows)

        methods = [
            ("IVW_fixed", ivw_fixed),
            ("IVW_random", ivw_random),
            ("weighted_median", wmed),
            ("weighted_mode", wmode),
            ("MR_Egger_slope", egger),
        ]

        for method, res in methods:
            if res is None:
                continue

            row = {
                "exposure_id": "SBP",
                "outcome_suffix": outcome,
                "method": method,
                "n_instruments": res.get("n", len(rows)),
                "beta": res.get("beta"),
                "se": res.get("se"),
                "pval": res.get("pval"),
                "direction": direction(res.get("beta")),
                "Q": res.get("Q", None),
                "Q_df": res.get("Q_df", None),
                "phi": res.get("phi", None),
                "egger_intercept": res.get("egger_intercept", None),
                "egger_intercept_se": res.get("egger_intercept_se", None),
                "egger_intercept_pval": res.get("egger_intercept_pval", None),
                "note": "pure_python_implementation_no_external_packages",
            }
            method_rows.append(row)

    return method_rows

def write_tsv(path, rows, cols):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: fmt(r.get(c)) for c in cols})

def svg_forest(path, rows, title, beta_col="beta", se_col="se"):
    points = []
    for r in rows:
        beta = parse_float(r.get(beta_col))
        se = parse_float(r.get(se_col))
        if beta is None or se is None:
            continue
        label = f"{r.get('outcome_suffix')} | {r.get('method')}"
        points.append({
            "label": label,
            "beta": beta,
            "lo": beta - 1.96 * se,
            "hi": beta + 1.96 * se,
            "p": r.get("pval", "NA"),
        })

    points = list(reversed(points))
    width = 1050
    height = max(360, 90 + 36 * len(points))
    left = 360
    right = width - 80
    top = 75
    row_h = 36

    vals = [0.0]
    for p in points:
        vals.extend([p["lo"], p["hi"]])
    xmin = min(vals)
    xmax = max(vals)
    pad = (xmax - xmin) * 0.15 if xmax > xmin else 0.1
    xmin -= pad
    xmax += pad

    def xscale(v):
        return left + (v - xmin) / (xmax - xmin) * (right - left)

    zero_x = xscale(0)

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="30" y="30" font-family="Arial" font-size="20" font-weight="bold">{escape(title)}</text>')
    svg.append('<text x="30" y="54" font-family="Arial" font-size="13">Beta estimates with 95% CI where SE is available.</text>')
    svg.append(f'<line x1="{zero_x:.1f}" y1="62" x2="{zero_x:.1f}" y2="{height-50}" stroke="black" stroke-dasharray="4,4" stroke-width="1"/>')

    for i, p in enumerate(points):
        y = top + i * row_h
        xlo = xscale(p["lo"])
        xhi = xscale(p["hi"])
        xb = xscale(p["beta"])
        label = f'{p["label"]} | p={p["p"]}'
        svg.append(f'<text x="30" y="{y+5}" font-family="Arial" font-size="11">{escape(label)}</text>')
        svg.append(f'<line x1="{xlo:.1f}" y1="{y}" x2="{xhi:.1f}" y2="{y}" stroke="black" stroke-width="1.3"/>')
        svg.append(f'<circle cx="{xb:.1f}" cy="{y}" r="4" fill="black"/>')

    axis_y = height - 35
    svg.append(f'<line x1="{left}" y1="{axis_y}" x2="{right}" y2="{axis_y}" stroke="black"/>')
    svg.append(f'<text x="{left}" y="{axis_y+18}" font-family="Arial" font-size="11">{xmin:.3g}</text>')
    svg.append(f'<text x="{zero_x-10}" y="{axis_y+18}" font-family="Arial" font-size="11">0</text>')
    svg.append(f'<text x="{right-40}" y="{axis_y+18}" font-family="Arial" font-size="11">{xmax:.3g}</text>')
    svg.append('</svg>')

    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(svg))

def svg_chr(path, rows):
    points = []
    for r in rows:
        beta = parse_float(r.get("beta"))
        se = parse_float(r.get("se"))
        if beta is None or se is None:
            continue
        points.append({
            "label": f"{r.get('outcome_suffix')} | leave chr {r.get('left_out_chr')}",
            "beta": beta,
            "lo": beta - 1.96 * se,
            "hi": beta + 1.96 * se,
            "p": r.get("pval", "NA"),
        })

    points = points[:60]
    width = 1150
    height = max(500, 80 + 22 * len(points))
    left = 400
    right = width - 70
    top = 70
    row_h = 22

    vals = [0.0]
    for p in points:
        vals.extend([p["lo"], p["hi"]])
    xmin = min(vals)
    xmax = max(vals)
    pad = (xmax - xmin) * 0.15 if xmax > xmin else 0.1
    xmin -= pad
    xmax += pad

    def xscale(v):
        return left + (v - xmin) / (xmax - xmin) * (right - left)

    zero_x = xscale(0)

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="30" y="30" font-family="Arial" font-size="20" font-weight="bold">SBP leave-one-chromosome-out IVW random-effects</text>')
    svg.append(f'<line x1="{zero_x:.1f}" y1="52" x2="{zero_x:.1f}" y2="{height-45}" stroke="black" stroke-dasharray="4,4" stroke-width="1"/>')

    for i, p in enumerate(points):
        y = top + i * row_h
        xlo = xscale(p["lo"])
        xhi = xscale(p["hi"])
        xb = xscale(p["beta"])
        svg.append(f'<text x="30" y="{y+4}" font-family="Arial" font-size="10">{escape(p["label"])}</text>')
        svg.append(f'<line x1="{xlo:.1f}" y1="{y}" x2="{xhi:.1f}" y2="{y}" stroke="black" stroke-width="1"/>')
        svg.append(f'<circle cx="{xb:.1f}" cy="{y}" r="3" fill="black"/>')

    svg.append('</svg>')

    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(svg))

def main():
    print("===== Phase 4.2 SBP-focused robustness package =====")
    print(f"Start: {now_str()}")

    random.seed(RANDOM_SEED)
    t0 = time.time()

    all_data = {}
    for outcome, path in SBP_INPUTS.items():
        all_data[outcome] = read_mr_input(path, outcome)

    method_rows = make_method_results(all_data)

    method_cols = [
        "exposure_id", "outcome_suffix", "method", "n_instruments",
        "beta", "se", "pval", "direction",
        "Q", "Q_df", "phi",
        "egger_intercept", "egger_intercept_se", "egger_intercept_pval",
        "note"
    ]
    write_tsv(METHOD_OUT, method_rows, method_cols)

    radial_summary = []
    radial_snps = []
    outlier_corrected_rows = []
    chr_rows = []
    influence_rows = []

    for outcome, rows in all_data.items():
        full_random = ivw_no_intercept(rows, random_effects=True)
        outliers, radial_stats = radial_outlier_check(rows, full_random, threshold=3.0)

        radial_summary.append({
            "exposure_id": "SBP",
            "outcome_suffix": outcome,
            "n_instruments": len(rows),
            "n_radial_outliers_abs_std_resid_gt3": len(outliers),
            "max_abs_std_residual": max([x["radial_abs_std_residual"] for x in radial_stats]) if radial_stats else None,
            "threshold": 3.0,
            "note": "Radial-style outlier screen using IVW random-effects residuals; not a replacement for formal RadialMR package",
        })

        for x in outliers:
            radial_snps.append({
                "exposure_id": "SBP",
                "outcome_suffix": outcome,
                "SNP": x["SNP"],
                "chr": x["chr"],
                "beta_exposure": x["beta_exposure"],
                "beta_outcome": x["beta_outcome"],
                "se_outcome": x["se_outcome"],
                "radial_residual": x["radial_residual"],
                "radial_std_residual": x["radial_std_residual"],
                "radial_abs_std_residual": x["radial_abs_std_residual"],
                "radial_outlier_flag": x["radial_outlier_flag"],
            })

        corrected = [r for r in rows if r["SNP"] not in set(x["SNP"] for x in outliers)]
        if len(corrected) >= 3:
            corr_est = ivw_no_intercept(corrected, random_effects=True)
        else:
            corr_est = None

        outlier_corrected_rows.append({
            "exposure_id": "SBP",
            "outcome_suffix": outcome,
            "method": "IVW_random_after_radial_outlier_exclusion",
            "n_original": len(rows),
            "n_outliers_removed": len(outliers),
            "n_remaining": len(corrected),
            "beta": corr_est["beta"] if corr_est else None,
            "se": corr_est["se"] if corr_est else None,
            "pval": corr_est["pval"] if corr_est else None,
            "direction": direction(corr_est["beta"]) if corr_est else "NA",
            "phi": corr_est["phi"] if corr_est else None,
            "note": "Outlier correction based on abs standardized residual > 3",
        })

        for r in leave_one_chromosome_out(rows):
            r["exposure_id"] = "SBP"
            r["outcome_suffix"] = outcome
            chr_rows.append(r)

        infl = single_snp_influence(rows)
        for r in infl[:25]:
            r["exposure_id"] = "SBP"
            r["outcome_suffix"] = outcome
            influence_rows.append(r)

    write_tsv(RADIAL_SUMMARY_OUT, radial_summary, [
        "exposure_id", "outcome_suffix", "n_instruments",
        "n_radial_outliers_abs_std_resid_gt3", "max_abs_std_residual",
        "threshold", "note"
    ])

    write_tsv(RADIAL_SNPS_OUT, radial_snps, [
        "exposure_id", "outcome_suffix", "SNP", "chr",
        "beta_exposure", "beta_outcome", "se_outcome",
        "radial_residual", "radial_std_residual",
        "radial_abs_std_residual", "radial_outlier_flag"
    ])

    write_tsv(OUTLIER_CORRECTED_OUT, outlier_corrected_rows, [
        "exposure_id", "outcome_suffix", "method", "n_original",
        "n_outliers_removed", "n_remaining", "beta", "se", "pval",
        "direction", "phi", "note"
    ])

    write_tsv(CHR_OUT, chr_rows, [
        "exposure_id", "outcome_suffix", "left_out_chr", "n_remaining",
        "beta", "se", "pval", "direction", "delta_from_full_beta", "full_beta"
    ])

    write_tsv(INFLUENCE_OUT, influence_rows, [
        "exposure_id", "outcome_suffix", "left_out_SNP", "chr",
        "beta_leave_one_out", "se_leave_one_out", "pval_leave_one_out",
        "delta_from_full_beta", "abs_delta_from_full_beta",
        "direction_leave_one_out", "full_beta"
    ])

    method_lookup = {}
    for r in method_rows:
        method_lookup[(r["outcome_suffix"], r["method"])] = r

    contrast_rows = []
    methods_for_contrast = ["IVW_fixed", "IVW_random", "weighted_median", "MR_Egger_slope"]

    for method in methods_for_contrast:
        non = method_lookup.get(("GBS_nonIOPcomponent", method))
        iop = method_lookup.get(("GBS_IOPcomponent", method))
        if not non or not iop:
            continue

        beta_non = parse_float(non["beta"])
        se_non = parse_float(non["se"])
        beta_iop = parse_float(iop["beta"])
        se_iop = parse_float(iop["se"])

        if None in (beta_non, se_non, beta_iop, se_iop):
            continue

        for r in R_VALUES:
            res = contrast(beta_iop, se_iop, beta_non, se_non, r)
            if res is None:
                continue
            beta_diff, se_diff, z, p, lo, hi = res
            contrast_rows.append({
                "exposure_id": "SBP",
                "method": method,
                "assumed_component_correlation_r": r,
                "beta_nonIOP": beta_non,
                "se_nonIOP": se_non,
                "direction_nonIOP": direction(beta_non),
                "beta_IOP": beta_iop,
                "se_IOP": se_iop,
                "direction_IOP": direction(beta_iop),
                "beta_difference_IOP_minus_nonIOP": beta_diff,
                "se_difference": se_diff,
                "z_contrast": z,
                "p_contrast": p,
                "ci_lower_difference": lo,
                "ci_upper_difference": hi,
                "contrast_direction": "IOP_MORE_POSITIVE_THAN_NONIOP" if beta_diff > 0 else "IOP_MORE_NEGATIVE_THAN_NONIOP",
            })

    write_tsv(CONTRAST_OUT, contrast_rows, [
        "exposure_id", "method", "assumed_component_correlation_r",
        "beta_nonIOP", "se_nonIOP", "direction_nonIOP",
        "beta_IOP", "se_IOP", "direction_IOP",
        "beta_difference_IOP_minus_nonIOP", "se_difference",
        "z_contrast", "p_contrast", "ci_lower_difference", "ci_upper_difference",
        "contrast_direction"
    ])

    non_methods = [r for r in method_rows if r["outcome_suffix"] == "GBS_nonIOPcomponent"]
    iop_methods = [r for r in method_rows if r["outcome_suffix"] == "GBS_IOPcomponent"]

    non_direction_negative = sum(1 for r in non_methods if r["direction"] == "negative")
    iop_direction_positive = sum(1 for r in iop_methods if r["direction"] == "positive")

    contrast_r0 = [r for r in contrast_rows if str(r["assumed_component_correlation_r"]) == "0.0" or str(r["assumed_component_correlation_r"]) == "0"]
    contrast_positive = sum(1 for r in contrast_r0 if r["contrast_direction"] == "IOP_MORE_POSITIVE_THAN_NONIOP")
    contrast_nominal = sum(1 for r in contrast_r0 if parse_float(r["p_contrast"]) is not None and parse_float(r["p_contrast"]) < 0.05)

    n_outliers_total = sum(parse_float(r["n_radial_outliers_abs_std_resid_gt3"]) or 0 for r in radial_summary)

    summary_rows = [
        {
            "item": "nonIOP_method_direction_consistency",
            "value": f"{non_direction_negative}/{len(non_methods)} negative",
            "interpretation": "Higher consistency supports SBP negative association with nonIOP component",
        },
        {
            "item": "IOP_method_direction_consistency",
            "value": f"{iop_direction_positive}/{len(iop_methods)} positive",
            "interpretation": "Higher consistency supports SBP positive association with IOP component",
        },
        {
            "item": "component_contrast_direction_consistency_r0",
            "value": f"{contrast_positive}/{len(contrast_r0)} IOP_MORE_POSITIVE_THAN_NONIOP",
            "interpretation": "Checks whether method-specific contrasts support directional divergence under r=0",
        },
        {
            "item": "component_contrast_nominal_count_r0",
            "value": f"{contrast_nominal}/{len(contrast_r0)} p<0.05",
            "interpretation": "Method-specific contrast significance under r=0",
        },
        {
            "item": "radial_style_outlier_count_total",
            "value": str(int(n_outliers_total)),
            "interpretation": "Outliers from abs standardized residual > 3 screen; formal MR-PRESSO/RadialMR package can be added later",
        },
        {
            "item": "core_package_status",
            "value": "COMPLETED",
            "interpretation": "Pure-Python SBP robustness package completed without external R packages",
        },
    ]

    write_tsv(SUMMARY_OUT, summary_rows, ["item", "value", "interpretation"])

    svg_forest(FIG_METHOD, method_rows, "SBP method consistency across GBS components")
    svg_chr(FIG_CHR, chr_rows)

    elapsed = time.time() - t0

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write("Phase 4.2 SBP-focused robustness package\tPASSED_CORE_PACKAGE\tPure-Python SBP robustness outputs generated\tIncludes method consistency, radial-style outliers, outlier-corrected IVW, leave-one-chromosome-out, single-SNP influence, and method-specific component contrasts\n")
        f.write(f"nonIOP_direction_consistency\tINFO\t{non_direction_negative}/{len(non_methods)} negative\tAcross available methods\n")
        f.write(f"IOP_direction_consistency\tINFO\t{iop_direction_positive}/{len(iop_methods)} positive\tAcross available methods\n")
        f.write(f"component_contrast_direction_consistency_r0\tINFO\t{contrast_positive}/{len(contrast_r0)} IOP_MORE_POSITIVE_THAN_NONIOP\tMethod-specific contrasts under assumed r=0\n")
        f.write(f"component_contrast_nominal_count_r0\tINFO\t{contrast_nominal}/{len(contrast_r0)} p<0.05\tMethod-specific contrasts under assumed r=0\n")
        f.write(f"radial_style_outlier_count_total\tINFO\t{int(n_outliers_total)}\tAbs standardized residual > 3 across two SBP outcome datasets\n")
        f.write("MR_PRESSO_status\tPLANNED_OPTIONAL\tNot run in this pure-Python core package\tRun formal MR-PRESSO in R if package is available\n")
        f.write("MR_RAPS_status\tPLANNED_OPTIONAL\tNot run in this pure-Python core package\tRun formal mr.raps in R if package is available\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 4.2 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstart_time\tend_time\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 4.2\tNA\t{now_str()}\t{elapsed:.6f}\t{elapsed:.1f}s\n")

    print("===== Phase 4.2 completed =====")
    print(f"End: {now_str()}")
    print(f"Elapsed: {elapsed:.1f}s")
    print(f"Wrote: {METHOD_OUT}")
    print(f"Wrote: {RADIAL_SUMMARY_OUT}")
    print(f"Wrote: {OUTLIER_CORRECTED_OUT}")
    print(f"Wrote: {CHR_OUT}")
    print(f"Wrote: {INFLUENCE_OUT}")
    print(f"Wrote: {CONTRAST_OUT}")
    print(f"Wrote: {SUMMARY_OUT}")
    print(f"Wrote: {STATUS_OUT}")
    print(f"Wrote: {FIG_METHOD}")
    print(f"Wrote: {FIG_CHR}")

if __name__ == "__main__":
    main()
