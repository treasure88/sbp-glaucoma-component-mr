#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

INDIR = "../../32_poag_external_validation_inputs/mr_input"
OUTDIR = "../../33_poag_external_mr_results"
os.makedirs(OUTDIR, exist_ok=True)

PAIRWISE = [
    ("SBP", "POAG", os.path.join(INDIR, "SBP__POAG.external_mr_input.tsv.gz")),
    ("ART_STIFFNESS", "POAG", os.path.join(INDIR, "ART_STIFFNESS__POAG.external_mr_input.tsv.gz")),
]

RESULTS_OUT = os.path.join(OUTDIR, "phase5_8E_POAG_external_mr_results.tsv")
HET_OUT = os.path.join(OUTDIR, "phase5_8E_POAG_external_heterogeneity.tsv")
EGGER_OUT = os.path.join(OUTDIR, "phase5_8E_POAG_external_egger_intercept.tsv")
SINGLE_OUT = os.path.join(OUTDIR, "phase5_8E_POAG_external_single_snp_ratios.tsv.gz")
LOO_OUT = os.path.join(OUTDIR, "phase5_8E_POAG_external_leave_one_out.tsv.gz")
STATUS_OUT = os.path.join(OUTDIR, "phase5_8E_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase5_8E_runtime_log.tsv")

def is_missing(x):
    return x is None or str(x).strip() in ("", ".", "NA", "NaN", "nan")

def fnum(x):
    try:
        if is_missing(x):
            return None
        return float(x)
    except Exception:
        return None

def normal_p(z):
    if z is None:
        return "NA"
    return f"{math.erfc(abs(z) / math.sqrt(2)):.12g}"

def direction(beta):
    b = fnum(beta)
    if b is None:
        return "NA"
    if b > 0:
        return "positive"
    if b < 0:
        return "negative"
    return "zero"

def chi2_p_approx(q, df):
    # Wilson-Hilferty approximation, adequate for QC flagging without scipy
    if q is None or df is None or df <= 0:
        return "NA"
    if q <= 0:
        return "1"
    z = ((q / df) ** (1.0 / 3.0) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    return normal_p(z)

def odds_ratio(beta):
    b = fnum(beta)
    if b is None:
        return "NA"
    try:
        return f"{math.exp(b):.12g}"
    except OverflowError:
        return "NA"

def odds_ratio_ci(beta, se):
    b = fnum(beta)
    s = fnum(se)
    if b is None or s is None:
        return ("NA", "NA")
    lo = b - 1.96 * s
    hi = b + 1.96 * s
    try:
        return (f"{math.exp(lo):.12g}", f"{math.exp(hi):.12g}")
    except OverflowError:
        return ("NA", "NA")

def read_mr_input(path):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            bx = fnum(r.get("beta_exposure"))
            by = fnum(r.get("beta_outcome"))
            sx = fnum(r.get("se_exposure"))
            sy = fnum(r.get("se_outcome"))
            if bx is None or by is None or sx is None or sy is None:
                continue
            if bx == 0 or sy <= 0:
                continue
            r["_bx"] = bx
            r["_by"] = by
            r["_sx"] = sx
            r["_sy"] = sy
            r["_w"] = 1.0 / (sy * sy)
            rows.append(r)
    return rows

def ivw(rows, random_effects=False):
    n = len(rows)
    if n < 1:
        return None
    swxx = sum(r["_w"] * r["_bx"] * r["_bx"] for r in rows)
    swxy = sum(r["_w"] * r["_bx"] * r["_by"] for r in rows)
    if swxx <= 0:
        return None

    beta = swxy / swxx
    q = sum(r["_w"] * (r["_by"] - beta * r["_bx"]) ** 2 for r in rows)
    df = n - 1
    phi = max(q / df, 1.0) if df > 0 else 1.0
    se = math.sqrt(1.0 / swxx)
    if random_effects:
        se *= math.sqrt(phi)

    z = beta / se if se > 0 else None
    return {
        "beta": beta,
        "se": se,
        "pval": normal_p(z),
        "Q": q,
        "Q_df": df,
        "Q_pval": chi2_p_approx(q, df),
        "phi": phi,
    }

def weighted_median(rows):
    vals = []
    for r in rows:
        bx = r["_bx"]
        by = r["_by"]
        sy = r["_sy"]
        ratio = by / bx
        ratio_se = abs(sy / bx)
        if ratio_se <= 0:
            continue
        w = 1.0 / (ratio_se * ratio_se)
        vals.append((ratio, w))
    if not vals:
        return None

    vals.sort(key=lambda x: x[0])
    total_w = sum(w for _, w in vals)
    c = 0.0
    for ratio, w in vals:
        c += w
        if c >= 0.5 * total_w:
            return ratio
    return vals[-1][0]

def egger(rows):
    n = len(rows)
    if n < 3:
        return None

    sw = sum(r["_w"] for r in rows)
    swx = sum(r["_w"] * r["_bx"] for r in rows)
    swy = sum(r["_w"] * r["_by"] for r in rows)
    swxx = sum(r["_w"] * r["_bx"] * r["_bx"] for r in rows)
    swxy = sum(r["_w"] * r["_bx"] * r["_by"] for r in rows)

    det = sw * swxx - swx * swx
    if det == 0:
        return None

    intercept = (swxx * swy - swx * swxy) / det
    slope = (sw * swxy - swx * swy) / det

    q = sum(r["_w"] * (r["_by"] - intercept - slope * r["_bx"]) ** 2 for r in rows)
    df = n - 2
    phi = max(q / df, 1.0) if df > 0 else 1.0

    se_intercept = math.sqrt(phi * swxx / det)
    se_slope = math.sqrt(phi * sw / det)

    return {
        "slope": slope,
        "slope_se": se_slope,
        "slope_pval": normal_p(slope / se_slope if se_slope > 0 else None),
        "intercept": intercept,
        "intercept_se": se_intercept,
        "intercept_pval": normal_p(intercept / se_intercept if se_intercept > 0 else None),
        "Q": q,
        "Q_df": df,
        "Q_pval": chi2_p_approx(q, df),
        "phi": phi,
    }

def bh_fdr(pvals):
    valid = []
    for i, p in enumerate(pvals):
        x = fnum(p)
        if x is not None:
            valid.append((i, x))
    m = len(valid)
    q = ["NA"] * len(pvals)
    if m == 0:
        return q

    valid_sorted = sorted(valid, key=lambda x: x[1])
    prev = 1.0
    for rank_from_end, (i, p) in enumerate(reversed(valid_sorted), start=1):
        rank = m - rank_from_end + 1
        val = min(prev, p * m / rank)
        prev = val
        q[i] = f"{val:.12g}"
    return q

def single_snp_ratios(exposure_id, outcome_id, rows):
    out = []
    for r in rows:
        ratio = r["_by"] / r["_bx"]
        se = abs(r["_sy"] / r["_bx"])
        p = normal_p(ratio / se if se > 0 else None)
        lo, hi = odds_ratio_ci(ratio, se)
        out.append({
            "exposure_id": exposure_id,
            "outcome_id": outcome_id,
            "SNP": r.get("SNP", "NA"),
            "beta_exposure": r.get("beta_exposure", "NA"),
            "se_exposure": r.get("se_exposure", "NA"),
            "beta_outcome_log_odds": r.get("beta_outcome", "NA"),
            "se_outcome": r.get("se_outcome", "NA"),
            "wald_ratio_log_odds": f"{ratio:.12g}",
            "wald_se": f"{se:.12g}",
            "wald_pval": p,
            "wald_or": odds_ratio(ratio),
            "wald_or_ci_lower": lo,
            "wald_or_ci_upper": hi,
            "direction": direction(ratio),
        })
    return out

def loo_results(exposure_id, outcome_id, rows):
    out = []
    n = len(rows)
    if n < 3:
        return out

    for i, r in enumerate(rows):
        sub = rows[:i] + rows[i+1:]
        res = ivw(sub, random_effects=True)
        if res is None:
            continue
        lo, hi = odds_ratio_ci(res["beta"], res["se"])
        out.append({
            "exposure_id": exposure_id,
            "outcome_id": outcome_id,
            "left_out_SNP": r.get("SNP", "NA"),
            "n_remaining": str(len(sub)),
            "beta_log_odds": f"{res['beta']:.12g}",
            "se": f"{res['se']:.12g}",
            "pval": res["pval"],
            "or": odds_ratio(res["beta"]),
            "or_ci_lower": lo,
            "or_ci_upper": hi,
            "direction": direction(res["beta"]),
            "Q": f"{res['Q']:.12g}",
            "Q_pval": res["Q_pval"],
        })
    return out

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def write_gz_tsv(path, fieldnames, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start_all = time.time()

    result_rows = []
    het_rows = []
    egger_rows = []
    single_rows = []
    loo_rows = []
    runtime_rows = []

    for exposure_id, outcome_id, path in PAIRWISE:
        t0 = time.time()
        print(f"Running POAG external MR: {exposure_id} -> {outcome_id}", flush=True)

        rows = read_mr_input(path)
        n = len(rows)

        fixed = ivw(rows, random_effects=False)
        random = ivw(rows, random_effects=True)
        wm = weighted_median(rows)
        eg = egger(rows)

        for method_name, res, phi_value, note in [
            ("IVW_fixed_effects", fixed, "1", "External POAG MR; fixed effects"),
            ("IVW_multiplicative_random_effects", random, None, "Primary external POAG MR method"),
        ]:
            if res is None:
                continue
            lo, hi = odds_ratio_ci(res["beta"], res["se"])
            result_rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "outcome_scale": "log_odds",
                "method": method_name,
                "n_instruments": str(n),
                "beta_log_odds": f"{res['beta']:.12g}",
                "se": f"{res['se']:.12g}",
                "pval": res["pval"],
                "qval_bh_ivw_random": "NA",
                "or": odds_ratio(res["beta"]),
                "or_ci_lower": lo,
                "or_ci_upper": hi,
                "direction": direction(res["beta"]),
                "Q": f"{res['Q']:.12g}",
                "Q_df": str(res["Q_df"]),
                "Q_pval": res["Q_pval"],
                "phi": phi_value if phi_value is not None else f"{res['phi']:.12g}",
                "note": note,
            })

        if random is not None:
            het_rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "n_instruments": str(n),
                "Q": f"{random['Q']:.12g}",
                "Q_df": str(random["Q_df"]),
                "Q_pval": random["Q_pval"],
                "heterogeneity_flag": "YES" if fnum(random["Q_pval"]) is not None and fnum(random["Q_pval"]) < 0.05 else "NO",
                "phi": f"{random['phi']:.12g}",
            })

        if wm is not None:
            result_rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "outcome_scale": "log_odds",
                "method": "weighted_median_ratio_descriptive",
                "n_instruments": str(n),
                "beta_log_odds": f"{wm:.12g}",
                "se": "NA",
                "pval": "NA",
                "qval_bh_ivw_random": "NA",
                "or": odds_ratio(wm),
                "or_ci_lower": "NA",
                "or_ci_upper": "NA",
                "direction": direction(wm),
                "Q": "NA",
                "Q_df": "NA",
                "Q_pval": "NA",
                "phi": "NA",
                "note": "Descriptive weighted median of SNP ratio estimates",
            })

        if eg is not None:
            lo, hi = odds_ratio_ci(eg["slope"], eg["slope_se"])
            result_rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "outcome_scale": "log_odds",
                "method": "MR_Egger_slope",
                "n_instruments": str(n),
                "beta_log_odds": f"{eg['slope']:.12g}",
                "se": f"{eg['slope_se']:.12g}",
                "pval": eg["slope_pval"],
                "qval_bh_ivw_random": "NA",
                "or": odds_ratio(eg["slope"]),
                "or_ci_lower": lo,
                "or_ci_upper": hi,
                "direction": direction(eg["slope"]),
                "Q": f"{eg['Q']:.12g}",
                "Q_df": str(eg["Q_df"]),
                "Q_pval": eg["Q_pval"],
                "phi": f"{eg['phi']:.12g}",
                "note": "MR-Egger slope; low power when instrument count is small",
            })

            egger_rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "n_instruments": str(n),
                "egger_intercept": f"{eg['intercept']:.12g}",
                "egger_intercept_se": f"{eg['intercept_se']:.12g}",
                "egger_intercept_pval": eg["intercept_pval"],
                "egger_slope": f"{eg['slope']:.12g}",
                "egger_slope_se": f"{eg['slope_se']:.12g}",
                "egger_slope_pval": eg["slope_pval"],
                "pleiotropy_flag": "YES" if fnum(eg["intercept_pval"]) is not None and fnum(eg["intercept_pval"]) < 0.05 else "NO",
                "note": "MR-Egger intercept test",
            })

        single_rows.extend(single_snp_ratios(exposure_id, outcome_id, rows))
        loo_rows.extend(loo_results(exposure_id, outcome_id, rows))

        elapsed = time.time() - t0
        runtime_rows.append({
            "exposure_id": exposure_id,
            "outcome_id": outcome_id,
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_human": f"{elapsed:.1f}s",
        })

    idx = [i for i, r in enumerate(result_rows) if r["method"] == "IVW_multiplicative_random_effects"]
    qs = bh_fdr([result_rows[i]["pval"] for i in idx])
    for j, i in enumerate(idx):
        result_rows[i]["qval_bh_ivw_random"] = qs[j]

    write_tsv(
        RESULTS_OUT,
        [
            "exposure_id", "outcome_id", "outcome_scale", "method",
            "n_instruments", "beta_log_odds", "se", "pval", "qval_bh_ivw_random",
            "or", "or_ci_lower", "or_ci_upper", "direction",
            "Q", "Q_df", "Q_pval", "phi", "note",
        ],
        result_rows,
    )

    write_tsv(
        HET_OUT,
        ["exposure_id", "outcome_id", "n_instruments", "Q", "Q_df", "Q_pval", "heterogeneity_flag", "phi"],
        het_rows,
    )

    write_tsv(
        EGGER_OUT,
        [
            "exposure_id", "outcome_id", "n_instruments",
            "egger_intercept", "egger_intercept_se", "egger_intercept_pval",
            "egger_slope", "egger_slope_se", "egger_slope_pval",
            "pleiotropy_flag", "note",
        ],
        egger_rows,
    )

    write_gz_tsv(
        SINGLE_OUT,
        [
            "exposure_id", "outcome_id", "SNP",
            "beta_exposure", "se_exposure",
            "beta_outcome_log_odds", "se_outcome",
            "wald_ratio_log_odds", "wald_se", "wald_pval",
            "wald_or", "wald_or_ci_lower", "wald_or_ci_upper", "direction",
        ],
        single_rows,
    )

    write_gz_tsv(
        LOO_OUT,
        [
            "exposure_id", "outcome_id", "left_out_SNP", "n_remaining",
            "beta_log_odds", "se", "pval",
            "or", "or_ci_lower", "or_ci_upper",
            "direction", "Q", "Q_pval",
        ],
        loo_rows,
    )

    write_tsv(
        RUNTIME_OUT,
        ["exposure_id", "outcome_id", "elapsed_seconds", "elapsed_human"],
        runtime_rows,
    )

    ivw_random = [r for r in result_rows if r["method"] == "IVW_multiplicative_random_effects"]
    nominal = sum(1 for r in ivw_random if fnum(r["pval"]) is not None and fnum(r["pval"]) < 0.05)
    fdr = sum(1 for r in ivw_random if fnum(r["qval_bh_ivw_random"]) is not None and fnum(r["qval_bh_ivw_random"]) < 0.05)

    elapsed_all = time.time() - start_all

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.8E POAG external MR\tPASSED\tivw_random_tests={len(ivw_random)};nominal={nominal};fdr={fdr}\tSBP and ART_STIFFNESS tested against POAG\n")
        f.write("outcome\tINFO\tPOAG\tExternal clinical glaucoma validation outcome\n")
        f.write("outcome_scale\tDOCUMENTED\tlog_odds\tMR beta is log odds; OR is exp(beta)\n")
        f.write("primary_external_exposure\tINFO\tSBP\tMain clinical endpoint validation exposure\n")
        f.write("secondary_external_exposure\tINFO\tART_STIFFNESS\tExploratory low-power exposure with 3 instruments\n")
        f.write(f"runtime\tINFO\t{elapsed_all:.3f}s\tPhase 5.8E runtime\n")

    print("===== Phase 5.8E completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", RESULTS_OUT)
    print("Wrote:", HET_OUT)
    print("Wrote:", EGGER_OUT)

if __name__ == "__main__":
    main()
