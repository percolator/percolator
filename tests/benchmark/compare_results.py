#!/usr/bin/env python3
"""Compare two Percolator result files (PSM, peptide or protein level).

Exits with status 1 if the files differ beyond the given tolerances.

  compare_results.py --exact base.txt new.txt
  compare_results.py base.txt new.txt   # tolerance mode

Tolerance mode joins rows on the first column (PSMId or protein id) and
reports the number of targets below the q-value threshold, the Jaccard index
of the accepted sets, the largest q-value and PEP differences and the Spearman
correlation of the scores.
"""

import argparse
import math
import sys


def read_table(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = {}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            rows[fields[0]] = fields
    return header, rows


def column(header, *names):
    for name in names:
        if name in header:
            return header.index(name)
    raise SystemExit("none of the columns %s found in header" % (names,))


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            result[order[k]] = (i + j) / 2.0
        i = j + 1
    return result


def spearman(a, b):
    """Spearman rank correlation; NaN (undefined) for fewer than two values or
    when either list is constant."""
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    if n < 2:
        return math.nan
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    va = sum((x - ma) ** 2 for x in ra)
    vb = sum((y - mb) ** 2 for y in rb)
    return cov / (va * vb) ** 0.5 if va > 0 and vb > 0 else math.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base")
    ap.add_argument("new")
    ap.add_argument("--exact", action="store_true", help="require byte-identical files")
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--max-count-change", type=float, default=0.002,
                    help="largest relative change of the number of accepted targets")
    ap.add_argument("--min-jaccard", type=float, default=0.995)
    ap.add_argument("--max-dq", type=float, default=1e-3)
    ap.add_argument("--max-dpep", type=float, default=None)
    args = ap.parse_args()

    if args.exact:
        with open(args.base, "rb") as f1, open(args.new, "rb") as f2:
            same = f1.read() == f2.read()
        print("identical" if same else "DIFFERENT", args.base, args.new)
        sys.exit(0 if same else 1)

    hb, base = read_table(args.base)
    hn, new = read_table(args.new)
    qb, qn = column(hb, "q-value"), column(hn, "q-value")
    pb, pn = column(hb, "posterior_error_prob"), column(hn, "posterior_error_prob")
    sb = column(hb, "score") if "score" in hb else None
    sn = column(hn, "score") if "score" in hn else None

    acc_b = {k for k, r in base.items() if float(r[qb]) < args.fdr}
    acc_n = {k for k, r in new.items() if float(r[qn]) < args.fdr}
    union = acc_b | acc_n
    jaccard = len(acc_b & acc_n) / len(union) if union else 1.0
    count_change = abs(len(acc_n) - len(acc_b)) / max(1, len(acc_b))
    common = sorted(base.keys() & new.keys())
    max_dq = max((abs(float(base[k][qb]) - float(new[k][qn])) for k in common), default=0.0)
    max_dpep = max((abs(float(base[k][pb]) - float(new[k][pn])) for k in common), default=0.0)
    rho = 1.0
    if sb is not None and sn is not None:
        rho = spearman([float(base[k][sb]) for k in common], [float(new[k][sn]) for k in common])

    print("rows base=%d new=%d common=%d" % (len(base), len(new), len(common)))
    print("accepted (q<%g) base=%d new=%d rel.change=%.5f jaccard=%.5f" %
          (args.fdr, len(acc_b), len(acc_n), count_change, jaccard))
    print("max|dq|=%.3g max|dPEP|=%.3g spearman(score)=%.6f" % (max_dq, max_dpep, rho))

    ok = (count_change <= args.max_count_change and jaccard >= args.min_jaccard
          and max_dq <= args.max_dq and (args.max_dpep is None or max_dpep <= args.max_dpep))
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
