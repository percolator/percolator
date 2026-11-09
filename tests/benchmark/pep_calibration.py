#!/usr/bin/env python3
"""Check the calibration of PEPs and q-values in a Percolator PSM result file
for a synthetic input made with make_synthetic_pin.py --mark-correct.

For bins of PEP it compares the mean PEP of the target PSMs with the fraction
of them that are incorrect, and for some q-value thresholds it compares the
threshold with the fraction of incorrect targets accepted.

  pep_calibration.py psms.txt
"""

import sys

EDGES = [0.0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0000001]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    rows = []
    with open(sys.argv[1]) as f:
        header = f.readline().rstrip("\n").split("\t")
        q_col = header.index("q-value")
        pep_col = header.index("posterior_error_prob")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            psm_id = fields[0]
            if psm_id.startswith("target_c_"):
                correct = True
            elif psm_id.startswith("target_i_"):
                correct = False
            else:
                continue
            rows.append((float(fields[q_col]), float(fields[pep_col]), correct))
    if not rows:
        sys.exit("no marked target PSMs found")

    print("PEP bin            n    mean PEP  incorrect")
    weighted_error = 0.0
    for lo, hi in zip(EDGES, EDGES[1:]):
        sel = [r for r in rows if lo <= r[1] < hi]
        if not sel:
            continue
        mean_pep = sum(r[1] for r in sel) / len(sel)
        incorrect = sum(1 for r in sel if not r[2]) / len(sel)
        weighted_error += len(sel) * abs(mean_pep - incorrect)
        print("[%.3f, %.3f) %8d   %.4f    %.4f" % (lo, min(hi, 1.0), len(sel), mean_pep, incorrect))
    print("mean |PEP - observed| weighted by bin size: %.4f" % (weighted_error / len(rows)))

    print("q threshold  accepted  observed FDR")
    for t in (0.001, 0.01, 0.05, 0.1):
        acc = [r for r in rows if r[0] < t]
        fdr = sum(1 for r in acc if not r[2]) / len(acc) if acc else 0.0
        print("%10.3f %9d   %.4f" % (t, len(acc), fdr))


if __name__ == "__main__":
    main()
