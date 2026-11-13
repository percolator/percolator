#!/usr/bin/env python3
"""Check the calibration of PEPs and q-values in a Percolator PSM or peptide
result file for a synthetic input made with make_synthetic_pin.py
--mark-correct.

For bins of PEP it compares the mean PEP of the target PSMs with the fraction
of them that are incorrect, and for some q-value thresholds it compares the
threshold with the fraction of incorrect targets accepted. It then reports
the mean PEP of the accepted targets (about the threshold for calibrated
PEPs), the PEP bins below 0.001, the Brier score and a Hosmer-Lemeshow
statistic along the score axis: the targets with PEP below --hl-max-pep,
ordered by score, are cut into bins of at least --hl-bin targets and at
least 5 expected incorrect targets (sum of PEPs), and for each bin
z = (incorrect - sum PEP) / sqrt(sum PEP (1 - PEP)). chi2/df is the mean of
z^2, about 1 for calibrated PEPs. Without the minimum expected count, a
single incorrect target among PEPs of 1e-6 would dominate the statistic.

A peptide result file is checked through the PSM ids of its peptides.

  pep_calibration.py [--hl-bin N] [--hl-max-pep P] psms.txt
"""

import argparse
import sys

EDGES = [0.0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0000001]
LOW_EDGES = [0.0, 0.0001, 0.001]
THRESHOLDS = (0.001, 0.01, 0.05, 0.1)
MIN_EXPECTED = 5.0


def read_marked(path):
    """(score, q, PEP, correct) of the marked target rows, in file order."""
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        q_col = header.index("q-value")
        pep_col = header.index("posterior_error_prob")
        score_col = header.index("score") if "score" in header else None
        for line in f:
            fields = line.rstrip("\n").split("\t")
            psm_id = fields[0]
            if psm_id.startswith("target_c_"):
                correct = True
            elif psm_id.startswith("target_i_"):
                correct = False
            else:
                continue
            score = float(fields[score_col]) if score_col is not None else 0.0
            rows.append((score, float(fields[q_col]), float(fields[pep_col]), correct))
    return rows, score_col is not None


def hosmer_lemeshow(rows, bin_size, max_pep, min_expected=MIN_EXPECTED):
    """chi2/df and number of bins along the score axis. A bin is closed when
    it has bin_size targets and min_expected expected incorrect targets, so
    that bins of tiny PEPs, where a single incorrect target would dominate
    the statistic, are merged; the incomplete last bin is dropped."""
    sel = sorted((r for r in rows if r[2] < max_pep), key=lambda r: -r[0])
    chi2, bins = 0.0, 0
    n = expected = observed = var = 0
    for r in sel:
        n += 1
        expected += r[2]
        observed += 0 if r[3] else 1
        var += r[2] * (1.0 - r[2])
        if n >= bin_size and expected >= min_expected:
            if var > 0.0:
                chi2 += (observed - expected) ** 2 / var
                bins += 1
            n = expected = observed = var = 0
    return (chi2 / bins if bins else float("nan")), bins


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result")
    ap.add_argument("--hl-bin", type=int, default=2000,
                    help="targets per Hosmer-Lemeshow bin (default 2000)")
    ap.add_argument("--hl-max-pep", type=float, default=0.5,
                    help="only targets with a PEP below this enter the Hosmer-Lemeshow "
                         "statistic (default 0.5)")
    args = ap.parse_args()
    if args.hl_bin < 1:
        ap.error("--hl-bin must be at least 1")
    rows, has_score = read_marked(args.result)
    if not rows:
        sys.exit("no marked target PSMs found")

    print("PEP bin            n    mean PEP  incorrect")
    weighted_error = 0.0
    for lo, hi in zip(EDGES, EDGES[1:]):
        sel = [r for r in rows if lo <= r[2] < hi]
        if not sel:
            continue
        mean_pep = sum(r[2] for r in sel) / len(sel)
        incorrect = sum(1 for r in sel if not r[3]) / len(sel)
        weighted_error += len(sel) * abs(mean_pep - incorrect)
        print("[%.3f, %.3f) %8d   %.4f    %.4f" % (lo, min(hi, 1.0), len(sel), mean_pep, incorrect))
    print("mean |PEP - observed| weighted by bin size: %.4f" % (weighted_error / len(rows)))

    print("q threshold  accepted  observed FDR")
    for t in THRESHOLDS:
        acc = [r for r in rows if r[1] < t]
        fdr = sum(1 for r in acc if not r[3]) / len(acc) if acc else 0.0
        print("%10.3f %9d   %.4f" % (t, len(acc), fdr))

    print("q threshold  mean PEP of the accepted targets")
    for t in THRESHOLDS:
        acc = [r for r in rows if r[1] < t]
        mean_pep = sum(r[2] for r in acc) / len(acc) if acc else 0.0
        print("%10.3f   %.5f" % (t, mean_pep))

    print("PEP bin              n    mean PEP    incorrect")
    for lo, hi in zip(LOW_EDGES, LOW_EDGES[1:]):
        sel = [r for r in rows if lo <= r[2] < hi]
        if not sel:
            continue
        mean_pep = sum(r[2] for r in sel) / len(sel)
        incorrect = sum(1 for r in sel if not r[3]) / len(sel)
        print("[%.4f, %.4f) %8d   %.3e   %.3e" % (lo, hi, len(sel), mean_pep, incorrect))

    brier = sum((r[2] - (0.0 if r[3] else 1.0)) ** 2 for r in rows) / len(rows)
    print("Brier score: %.6f" % brier)
    if has_score:
        chi2_df, bins = hosmer_lemeshow(rows, args.hl_bin, args.hl_max_pep)
        print("Hosmer-Lemeshow chi2/df along the score: %.3f (%d bins of %d targets with PEP < %g)"
              % (chi2_df, bins, args.hl_bin, args.hl_max_pep))


if __name__ == "__main__":
    main()
