#!/usr/bin/env python3
"""Compare the PEPs of Percolator result files with the true PEPs that
make_synthetic_pin.py wrote for the same input: with --truth for PSM result
files, matched by PSM id, and with --peptide-truth for peptide result files,
matched by peptide. (With --truth, a peptide file is compared through the
PSM ids of its peptides, i.e. with the PSM-level true PEP of its best PSM,
which is lower than the peptide's at the top of the list.)

For bands of target ranks (by score, best first) it prints the median of
log10(PEP / true PEP), the median of its absolute value, and the mean PEP
and mean true PEP of the band; then the mean |PEP - true PEP| over all
targets and the Brier score of the true PEPs, the lowest Brier score any
PEP estimate can reach on average. PEPs of 0 are counted as 1e-300 in the
log ratio.

  true_pep.py truth.txt result.txt [result2.txt ...]
"""

import math
import sys

BANDS = ((1, 10), (11, 100), (101, 1000), (1001, 10000), (10001, 100000))
TINY = 1e-300


def read_truth(path):
    """truth by PSM id (--truth) or by peptide (--peptide-truth), and
    whether it is by peptide"""
    truth = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        by_peptide = "Peptide" in header
        key_col = header.index("Peptide" if by_peptide else "SpecId")
        correct_col, pep_col = header.index("correct"), header.index("true_pep")
        label_col = None if by_peptide else header.index("Label")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if label_col is not None and fields[label_col] != "1":
                continue
            truth[fields[key_col]] = (fields[correct_col] == "1", float(fields[pep_col]))
    return truth, by_peptide


def read_targets(path, truth, by_peptide):
    """(score, PEP, true PEP, correct) of the target rows found in truth,
    best score first; the first row of a repeated key is used."""
    rows, seen = [], set()
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        score_col = header.index("score")
        pep_col = header.index("posterior_error_prob")
        key_col = header.index("peptide") if by_peptide else 0
        for line in f:
            fields = line.rstrip("\n").split("\t")
            key = fields[key_col]
            if key in seen or key not in truth:
                continue
            seen.add(key)
            correct, true_pep = truth[key]
            rows.append((float(fields[score_col]), float(fields[pep_col]), true_pep, correct))
    rows.sort(key=lambda r: -r[0])  # stable: file order within ties
    return rows


def median(values):
    v = sorted(values)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    truth, by_peptide = read_truth(sys.argv[1])
    for path in sys.argv[2:]:
        rows = read_targets(path, truth, by_peptide)
        if not rows:
            print("%s: no targets of the truth file found" % path)
            continue
        print("%s: %d targets" % (path, len(rows)))
        print("  ranks           median log10(PEP/true)  median |log10|  mean PEP   mean true PEP")
        for lo, hi in BANDS:
            band = rows[lo - 1:hi]
            if not band:
                break
            ratios = [math.log10(max(r[1], TINY) / max(r[2], TINY)) for r in band]
            print("  %6d-%-8d %+10.2f              %8.2f      %.3e  %.3e"
                  % (lo, lo + len(band) - 1, median(ratios), median(abs(x) for x in ratios),
                     sum(r[1] for r in band) / len(band), sum(r[2] for r in band) / len(band)))
        zeros = sum(1 for r in rows if r[1] <= 0.0)
        mae = sum(abs(r[1] - r[2]) for r in rows) / len(rows)
        brier_true = sum((r[2] - (0.0 if r[3] else 1.0)) ** 2 for r in rows) / len(rows)
        print("  mean |PEP - true PEP| = %.5f, Brier score of the true PEPs = %.6f%s"
              % (mae, brier_true, ", %d PEPs are 0" % zeros if zeros else ""))


if __name__ == "__main__":
    main()
