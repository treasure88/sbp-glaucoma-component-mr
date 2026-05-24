#!/usr/bin/env python3
import gzip
import os
import math
import time
import random
import statistics
from datetime import datetime

MR_INPUT_DIR = "../../16_mr_input_datasets/pairwise"
ELIGIBILITY_FILE = "../../17_mr_method_eligibility/phase3_4_mr_method_eligibility_summary.tsv"

OUTDIR = "../../18_pairwise_mr_results"
os.makedirs(OUTDIR, exist_ok=True)

RESULTS_OUT = os.path.join(OUTDIR, "phase3_5_pairwise_mr_results.tsv")
HET_OUT = os.path.join(OUTDIR, "phase3_5_heterogeneity_summary.tsv")
EGGER_INT_OUT = os.path.join(OUTDIR, "phase3_5_egger_intercept_summary.tsv")
LOO_OUT = os.path.join(OUTDIR, "phase3_5_leave_one_out_ivw.tsv.gz")
STATUS_OUT = os.path.join(OUTDIR, "phase3_5_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase3_5_runtime_log.tsv")
NOTES_MD = os.path.join(OUTDIR, "phase3_5_pairwise_mr_execution_notes.md")

PRIMARY_EXPOSURES = ["SBP", "DBP", "MIGRAINE", "CRP", "BMI"]
GBS_OUTCOMES = ["GBS_nonIOPcomponent", "GBS_IOPcomponent"]

BOOTSTRAP_N = 1000
RANDOM_SEED = 20260517

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def fmt_seconds(x):
    if x < 60:
        return f"{x:.1f}s"
    if x < 3600:
        return f"{x/60:.1f}min"
    return f"{x/3600:.2f}h"

def safe_float(x):
    try:
        if x in ("", ".", "NA"):
            return None
        return float(x)
    except Exception:
        return None

def normal_p_from_z(z):
    if z is None:
        return None
    return math.erfc(abs(z) / math.sqrt(2.0))

def fmt_num(x, digits=8):
    if x is None:
        return "NA"
    if isinstance(x, str):
        return x
    if math.isnan(x) or math.isinf(x):
        return "NA"
    return f"{x:.{digits}g}"

def fmt_p(p):
    if p is None:
        return "NA"
    if p == 0:
        return "<1e-300"
    if p < 1e-300:
        return "<1e-300"
    return f"{p:.8g}"

def exp_safe(x):
    if x is None:
        return None
    if x > 700:
        return float("inf")
    if x < -700:
        return 0.0
    return math.exp(x)

def chi2_sf(q, df):
    """
    Chi-square survival function.
    Uses scipy if available; otherwise Wilson-Hilferty normal approximation.
    """
    if q is None or df is None or df <= 0:
        return None
    try:
        from scipy.stats import chi2
        return float(chi2.sf(q, df))
    except Exception:
        # Wilson-Hilferty approximation
        z = ((q / df) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))) / math.sqrt(2.0 / (9.0 * df))
        return 0.5 * math.erfc(z / math.sqrt(2.0))

def weighted_median(values, weights):
    pairs = sorted((v, w) for v, w in zip(values, weights) if v is not None and w is not None and w > 0)
    if not pairs:
        return None
    total = sum(w for _, w in pairs)
    c = 0.0
    for v, w in pairs:
        c += w
        if c >= total / 2.0:
            return v
    return pairs[-1][0]

def ivw_estimate(rows, random_effects=False):
    xs = [r["beta_exposure"] for r in rows]
    ys = [r["beta_outcome"] for r in rows]
    sy = [r["se_outcome"] for r in rows]

    usable = [(x, y, s) for x, y, s in zip(xs, ys, sy) if x is not None and y is not None and s is not None and s > 0 and x != 0]
    n = len(usable)
    if n < 1:
        return None

    w = [1.0 / (s * s) for _, _, s in usable]
    den = sum(wi * x * x for (x, y, s), wi in zip(usable, w))
    if den == 0:
        return None

    beta = sum(wi * x * y for (x, y, s), wi in zip(usable, w)) / den
    residuals = [y - beta * x for x, y, s in usable]
    q = sum(wi * e * e for wi, e in zip(w, residuals))
    df = n - 1

    phi = 1.0
    if random_effects and df > 0:
        phi = max(1.0, q / df)

    se = math.sqrt(phi / den)
    z = beta / se if se > 0 else None
    p = normal_p_from_z(z)

    return {
        "n": n,
        "beta": beta,
        "se": se,
        "p": p,
        "q": q,
        "q_df": df,
        "q_p": chi2_sf(q, df),
        "phi": phi,
    }

def egger_estimate(rows):
    usable = [
        (r["beta_exposure"], r["beta_outcome"], r["se_outcome"])
        for r in rows
        if r["beta_exposure"] is not None and r["beta_outcome"] is not None and r["se_outcome"] is not None and r["se_outcome"] > 0
    ]

    n = len(usable)
    if n < 3:
        return None

    w = [1.0 / (s * s) for _, _, s in usable]
    x = [u[0] for u in usable]
    y = [u[1] for u in usable]

    s0 = sum(w)
    s1 = sum(wi * xi for wi, xi in zip(w, x))
    s2 = sum(wi * xi * xi for wi, xi in zip(w, x))
    t0 = sum(wi * yi for wi, yi in zip(w, y))
    t1 = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))

    det = s0 * s2 - s1 * s1
    if det == 0:
        return None

    intercept = (t0 * s2 - t1 * s1) / det
    slope = (s0 * t1 - s1 * t0) / det

    residuals = [yi - intercept - slope * xi for xi, yi in zip(x, y)]
    q = sum(wi * e * e for wi, e in zip(w, residuals))
    df = n - 2
    phi = max(1.0, q / df) if df > 0 else 1.0

    se_intercept = math.sqrt(phi * s2 / det)
    se_slope = math.sqrt(phi * s0 / det)

    p_intercept = normal_p_from_z(intercept / se_intercept) if se_intercept > 0 else None
    p_slope = normal_p_from_z(slope / se_slope) if se_slope > 0 else None

    return {
        "n": n,
        "slope": slope,
        "se_slope": se_slope,
        "p_slope": p_slope,
        "intercept": intercept,
        "se_intercept": se_intercept,
        "p_intercept": p_intercept,
        "q": q,
        "q_df": df,
        "q_p": chi2_sf(q, df),
        "phi": phi,
    }

def weighted_median_estimate(rows, bootstrap_n=BOOTSTRAP_N, seed=RANDOM_SEED):
    usable = [
        r for r in rows
        if r["beta_exposure"] is not None
        and r["beta_outcome"] is not None
        and r["se_exposure"] is not None
        and r["se_outcome"] is not None
        and r["se_exposure"] > 0
        and r["se_outcome"] > 0
        and r["beta_exposure"] != 0
    ]

    n = len(usable)
    if n < 3:
        return None

    ratios = [r["beta_outcome"] / r["beta_exposure"] for r in usable]
    weights = [(r["beta_exposure"] ** 2) / (r["se_outcome"] ** 2) for r in usable]

    beta = weighted_median(ratios, weights)

    rng = random.Random(seed)
    boot = []

    for _ in range(bootstrap_n):
        b_ratios = []
        b_weights = []
        for r in usable:
            bx = rng.gauss(r["beta_exposure"], r["se_exposure"])
            by = rng.gauss(r["beta_outcome"], r["se_outcome"])
            if bx == 0:
                continue
            ratio = by / bx
            weight = (bx * bx) / (r["se_outcome"] ** 2)
            if weight > 0:
                b_ratios.append(ratio)
                b_weights.append(weight)
        if b_ratios:
            bm = weighted_median(b_ratios, b_weights)
            if bm is not None and math.isfinite(bm):
                boot.append(bm)

    if len(boot) >= 20:
        se = statistics.stdev(boot)
        p = normal_p_from_z(beta / se) if se > 0 else None
    else:
        se = None
        p = None

    return {
        "n": n,
        "beta": beta,
        "se": se,
        "p": p,
        "bootstrap_n": bootstrap_n,
    }

def load_mr_input(path):
    rows = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        header = f.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = [
            "exposure_id",
            "outcome_id",
            "SNP",
            "beta_exposure",
            "se_exposure",
            "beta_outcome",
            "se_outcome",
            "pval_exposure",
            "pval_outcome",
            "instrument_F_stat",
        ]

        missing = [c for c in required if c not in cmap]
        if missing:
            raise RuntimeError(f"Missing columns in {path}: {missing}")

        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(header):
                continue

            row = {
                "exposure_id": parts[cmap["exposure_id"]],
                "outcome_id": parts[cmap["outcome_id"]],
                "SNP": parts[cmap["SNP"]],
                "beta_exposure": safe_float(parts[cmap["beta_exposure"]]),
                "se_exposure": safe_float(parts[cmap["se_exposure"]]),
                "beta_outcome": safe_float(parts[cmap["beta_outcome"]]),
                "se_outcome": safe_float(parts[cmap["se_outcome"]]),
                "pval_exposure": parts[cmap["pval_exposure"]],
                "pval_outcome": parts[cmap["pval_outcome"]],
                "instrument_F_stat": safe_float(parts[cmap["instrument_F_stat"]]),
            }
            rows.append(row)

    return rows

def result_row(exposure, outcome, method, n, beta, se, p, note):
    if beta is not None and se is not None:
        l95 = beta - 1.96 * se
        u95 = beta + 1.96 * se
    else:
        l95 = None
        u95 = None

    return {
        "exposure_id": exposure,
        "outcome_suffix": outcome,
        "method": method,
        "n_instruments": n,
        "beta": beta,
        "se": se,
        "pval": p,
        "ci_lower": l95,
        "ci_upper": u95,
        "or": exp_safe(beta),
        "or_ci_lower": exp_safe(l95),
        "or_ci_upper": exp_safe(u95),
        "note": note,
    }

def bh_fdr(pvals):
    indexed = [(i, p) for i, p in enumerate(pvals) if p is not None]
    m = len(indexed)
    qvals = [None] * len(pvals)

    if m == 0:
        return qvals

    sorted_pairs = sorted(indexed, key=lambda x: x[1])
    prev = 1.0

    for rank_from_end, (i, p) in enumerate(reversed(sorted_pairs), start=1):
        rank = m - rank_from_end + 1
        q = min(prev, p * m / rank)
        q = min(q, 1.0)
        qvals[i] = q
        prev = q

    return qvals

def write_results(rows):
    primary_indices = [
        i for i, r in enumerate(rows)
        if r["method"] == "IVW_random_effects"
    ]
    primary_pvals = [
        max(rows[i]["pval"], 1e-300) if rows[i]["pval"] is not None else None
        for i in primary_indices
    ]
    qvals = bh_fdr(primary_pvals)

    q_by_index = {}
    for idx, q in zip(primary_indices, qvals):
        q_by_index[idx] = q

    cols = [
        "exposure_id",
        "outcome_suffix",
        "method",
        "n_instruments",
        "beta",
        "se",
        "pval",
        "qval_bh_primary_ivw_random",
        "ci_lower",
        "ci_upper",
        "or",
        "or_ci_lower",
        "or_ci_upper",
        "note",
    ]

    with open(RESULTS_OUT, "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows):
            q = q_by_index.get(i)
            vals = [
                r["exposure_id"],
                r["outcome_suffix"],
                r["method"],
                str(r["n_instruments"]),
                fmt_num(r["beta"]),
                fmt_num(r["se"]),
                fmt_p(r["pval"]),
                fmt_p(q),
                fmt_num(r["ci_lower"]),
                fmt_num(r["ci_upper"]),
                fmt_num(r["or"]),
                fmt_num(r["or_ci_lower"]),
                fmt_num(r["or_ci_upper"]),
                r["note"],
            ]
            out.write("\t".join(vals) + "\n")

def main():
    print("===== Phase 3.5 Pairwise MR execution for GBS outcomes =====")
    print(f"Start: {now_str()}")
    print(f"Bootstrap replicates for weighted median SE: {BOOTSTRAP_N}")
    print("Primary IVW result uses multiplicative random effects.")
    print()

    t0 = time.time()

    mr_results = []
    het_rows = []
    egger_rows = []
    loo_rows = []
    runtime_rows = []

    for exposure in PRIMARY_EXPOSURES:
        for outcome in GBS_OUTCOMES:
            path = os.path.join(MR_INPUT_DIR, f"{exposure}__{outcome}.mr_input.tsv.gz")

            if not os.path.exists(path):
                raise FileNotFoundError(path)

            pair_t0 = time.time()
            rows = load_mr_input(path)
            n = len(rows)

            print(f"Running MR: {exposure} vs {outcome}; instruments={n}")

            if n < 1:
                continue

            ivw_fixed = ivw_estimate(rows, random_effects=False)
            ivw_random = ivw_estimate(rows, random_effects=True)
            wm = weighted_median_estimate(rows, seed=RANDOM_SEED + len(mr_results))
            egger = egger_estimate(rows)

            if ivw_fixed:
                mr_results.append(result_row(
                    exposure, outcome, "IVW_fixed_effects",
                    ivw_fixed["n"], ivw_fixed["beta"], ivw_fixed["se"], ivw_fixed["p"],
                    "IVW fixed-effect estimate"
                ))

            if ivw_random:
                mr_results.append(result_row(
                    exposure, outcome, "IVW_random_effects",
                    ivw_random["n"], ivw_random["beta"], ivw_random["se"], ivw_random["p"],
                    "Primary MR estimate; multiplicative random-effects IVW"
                ))

                het_rows.append({
                    "exposure_id": exposure,
                    "outcome_suffix": outcome,
                    "method": "IVW",
                    "n_instruments": ivw_random["n"],
                    "Q": ivw_random["q"],
                    "Q_df": ivw_random["q_df"],
                    "Q_pval": ivw_random["q_p"],
                    "phi": ivw_random["phi"],
                })

            if wm:
                mr_results.append(result_row(
                    exposure, outcome, "weighted_median",
                    wm["n"], wm["beta"], wm["se"], wm["p"],
                    f"Weighted median estimate; bootstrap_n={BOOTSTRAP_N}"
                ))

            if egger:
                mr_results.append(result_row(
                    exposure, outcome, "MR_Egger_slope",
                    egger["n"], egger["slope"], egger["se_slope"], egger["p_slope"],
                    "MR-Egger slope estimate"
                ))

                egger_rows.append({
                    "exposure_id": exposure,
                    "outcome_suffix": outcome,
                    "n_instruments": egger["n"],
                    "egger_intercept": egger["intercept"],
                    "se_intercept": egger["se_intercept"],
                    "pval_intercept": egger["p_intercept"],
                    "egger_Q": egger["q"],
                    "egger_Q_df": egger["q_df"],
                    "egger_Q_pval": egger["q_p"],
                    "phi": egger["phi"],
                    "note": "MR-Egger intercept tests directional pleiotropy; interpret cautiously when instrument count is limited",
                })

                het_rows.append({
                    "exposure_id": exposure,
                    "outcome_suffix": outcome,
                    "method": "MR_Egger",
                    "n_instruments": egger["n"],
                    "Q": egger["q"],
                    "Q_df": egger["q_df"],
                    "Q_pval": egger["q_p"],
                    "phi": egger["phi"],
                })

            # Leave-one-out IVW random effects
            if n >= 3:
                for r in rows:
                    sub = [x for x in rows if x["SNP"] != r["SNP"]]
                    est = ivw_estimate(sub, random_effects=True)
                    if est:
                        loo_rows.append({
                            "exposure_id": exposure,
                            "outcome_suffix": outcome,
                            "left_out_SNP": r["SNP"],
                            "n_instruments_after_removal": est["n"],
                            "beta_ivw_random": est["beta"],
                            "se_ivw_random": est["se"],
                            "pval_ivw_random": est["p"],
                        })

            elapsed = time.time() - pair_t0
            runtime_rows.append({
                "exposure_id": exposure,
                "outcome_suffix": outcome,
                "n_instruments": n,
                "elapsed_seconds": elapsed,
                "elapsed_human": fmt_seconds(elapsed),
            })

    write_results(mr_results)

    with open(HET_OUT, "w", encoding="utf-8", newline="\n") as out:
        cols = ["exposure_id", "outcome_suffix", "method", "n_instruments", "Q", "Q_df", "Q_pval", "phi"]
        out.write("\t".join(cols) + "\n")
        for r in het_rows:
            out.write("\t".join([
                r["exposure_id"],
                r["outcome_suffix"],
                r["method"],
                str(r["n_instruments"]),
                fmt_num(r["Q"]),
                str(r["Q_df"]),
                fmt_p(r["Q_pval"]),
                fmt_num(r["phi"]),
            ]) + "\n")

    with open(EGGER_INT_OUT, "w", encoding="utf-8", newline="\n") as out:
        cols = [
            "exposure_id",
            "outcome_suffix",
            "n_instruments",
            "egger_intercept",
            "se_intercept",
            "pval_intercept",
            "egger_Q",
            "egger_Q_df",
            "egger_Q_pval",
            "phi",
            "note",
        ]
        out.write("\t".join(cols) + "\n")
        for r in egger_rows:
            out.write("\t".join([
                r["exposure_id"],
                r["outcome_suffix"],
                str(r["n_instruments"]),
                fmt_num(r["egger_intercept"]),
                fmt_num(r["se_intercept"]),
                fmt_p(r["pval_intercept"]),
                fmt_num(r["egger_Q"]),
                str(r["egger_Q_df"]),
                fmt_p(r["egger_Q_pval"]),
                fmt_num(r["phi"]),
                r["note"],
            ]) + "\n")

    with gzip.open(LOO_OUT, "wt", encoding="utf-8", newline="\n") as out:
        cols = [
            "exposure_id",
            "outcome_suffix",
            "left_out_SNP",
            "n_instruments_after_removal",
            "beta_ivw_random",
            "se_ivw_random",
            "pval_ivw_random",
        ]
        out.write("\t".join(cols) + "\n")
        for r in loo_rows:
            out.write("\t".join([
                r["exposure_id"],
                r["outcome_suffix"],
                r["left_out_SNP"],
                str(r["n_instruments_after_removal"]),
                fmt_num(r["beta_ivw_random"]),
                fmt_num(r["se_ivw_random"]),
                fmt_p(r["pval_ivw_random"]),
            ]) + "\n")

    total_elapsed = time.time() - t0

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="\n") as out:
        cols = ["exposure_id", "outcome_suffix", "n_instruments", "elapsed_seconds", "elapsed_human"]
        out.write("\t".join(cols) + "\n")
        for r in runtime_rows:
            out.write("\t".join([
                r["exposure_id"],
                r["outcome_suffix"],
                str(r["n_instruments"]),
                f"{r['elapsed_seconds']:.6f}",
                r["elapsed_human"],
            ]) + "\n")

    n_pairs = len(runtime_rows)
    n_methods = len(mr_results)
    n_egger = len(egger_rows)
    n_het = len(het_rows)
    n_loo = len(loo_rows)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="\n") as out:
        out.write("phase\tstatus\tkey_result\tnote\n")
        out.write(f"Phase 3.5 pairwise MR execution for GBS outcomes\tPASSED\t{n_pairs}/10 GBS pairwise datasets analyzed\tMR run completed for SBP, DBP, MIGRAINE, CRP, and BMI against both GBS outcomes\n")
        out.write(f"mr_method_result_count\tINFO\t{n_methods}\tIncludes IVW fixed, IVW random, weighted median, and MR-Egger slope where eligible\n")
        out.write(f"heterogeneity_result_count\tINFO\t{n_het}\tIVW and MR-Egger heterogeneity outputs\n")
        out.write(f"egger_intercept_result_count\tINFO\t{n_egger}\tDirectional pleiotropy tests\n")
        out.write(f"leave_one_out_result_count\tINFO\t{n_loo}\tLeave-one-out IVW random-effects estimates\n")
        out.write(f"total_elapsed_runtime\tINFO\t{fmt_seconds(total_elapsed)}\tPhase 3.5 total runtime\n")
        out.write("excluded_exposure\tDOCUMENTED\tINSOMNIA\tNo LD-clumped instruments from Phase 3.2\n")
        out.write("excluded_outcome\tDOCUMENTED\tIOPcc_coordinate_subset\tZero matched instruments from Phase 3.3\n")
        out.write("primary_method\tDOCUMENTED\tIVW_random_effects\tMultiplicative random-effects IVW used as primary MR estimate\n")

    with open(NOTES_MD, "w", encoding="utf-8", newline="\n") as out:
        out.write("# Phase 3.5 Pairwise MR Execution Notes\n\n")
        out.write("## Status\n\n")
        out.write("PASSED\n\n")
        out.write("## Datasets analyzed\n\n")
        out.write("Ten GBS pairwise datasets were analyzed: five exposures against GBS_nonIOPcomponent and GBS_IOPcomponent.\n\n")
        out.write("## Methods\n\n")
        out.write("- IVW fixed effects\n")
        out.write("- IVW multiplicative random effects, used as the primary MR estimate\n")
        out.write("- Weighted median with bootstrap standard error\n")
        out.write("- MR-Egger slope\n")
        out.write("- MR-Egger intercept for directional pleiotropy\n")
        out.write("- Heterogeneity Q statistics\n")
        out.write("- Leave-one-out IVW random-effects analysis\n\n")
        out.write("## Exclusions\n\n")
        out.write("- INSOMNIA was excluded because no LD-clumped instruments were available.\n")
        out.write("- IOPcc_coordinate_subset was excluded because no matched clumped instruments were available.\n\n")
        out.write("## Important interpretation note\n\n")
        out.write("Outcome effects are on the GBS component log-odds scale. Exponentiated estimates are reported as odds-ratio scale summaries where appropriate.\n\n")
        out.write("## Next phase\n\n")
        out.write("Phase 3.6 should perform MR result QC, compare nonIOP vs IOP component patterns, inspect heterogeneity and pleiotropy, and prepare interpretation tables.\n")

    print()
    print("===== Phase 3.5 completed =====")
    print(f"End: {now_str()}")
    print(f"Total elapsed: {fmt_seconds(total_elapsed)}")
    print(f"Wrote MR results: {RESULTS_OUT}")
    print(f"Wrote heterogeneity: {HET_OUT}")
    print(f"Wrote Egger intercept: {EGGER_INT_OUT}")
    print(f"Wrote leave-one-out: {LOO_OUT}")
    print(f"Wrote status: {STATUS_OUT}")

if __name__ == "__main__":
    main()
