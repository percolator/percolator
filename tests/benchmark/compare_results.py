#!/usr/bin/env python3
"""Compare two Percolator result files (PSM, peptide or protein level).

Exits with status 1 if the files differ beyond the given tolerances.

  compare_results.py --exact base.txt new.txt
  compare_results.py base.txt new.txt   # tolerance mode

Tolerance mode joins the rows of the two files and reports the number of
targets below the q-value threshold, the Jaccard index of the accepted sets,
the largest q-value and PEP differences and the Spearman correlation of the
scores. The rows are joined on (--key auto):
  - the peptide, if every row of both files has a different peptide (peptide
    level: the PSMId names the best PSM of a peptide, which may differ);
  - otherwise the PSMId and the peptide (PSM level: the input may repeat a
    SpecId for different peptides);
  - the first column if there is no peptide column (protein level).
A PSM-level file in which every peptide occurs once looks like a peptide-level
file. If such files are joined on the peptide and some peptides have a
different PSMId in the two files, a warning gives their number; use --key psm
for PSM-level files. If a file repeats a key, the row with the highest score
is used (the first of them on ties, or the first row without a score column),
and the number of repeated keys is reported on stderr.
"""

import argparse
import math
import sys


def read_rows(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [(line_nr, line.rstrip("\n").split("\t"))
                for line_nr, line in enumerate(f, 2) if line.strip()]
    return header, rows


def key_function(mode, header_base, rows_base, header_new, rows_new):
    """Return the name of the join key and a function from a header to a key
    function of a row."""
    def peptide_unique(header, rows):
        ix = header.index("peptide")
        peptides = set()
        for _, fields in rows:
            if fields[ix] in peptides:
                return False
            peptides.add(fields[ix])
        return True
    has_peptide = "peptide" in header_base and "peptide" in header_new
    if mode == "auto":
        if not has_peptide:
            mode = "first"
        elif peptide_unique(header_base, rows_base) and peptide_unique(header_new, rows_new):
            mode = "peptide"
            warn_if_psmids_differ(header_base, rows_base, header_new, rows_new)
        else:
            mode = "psm"
    if mode != "first" and not has_peptide:
        raise SystemExit("--key %s needs a peptide column in both files" % mode)
    if mode == "first":
        return mode, lambda header: (lambda fields: fields[0])
    if mode == "peptide":
        return mode, lambda header: (lambda fields, ix=header.index("peptide"): fields[ix])
    return mode, lambda header: (lambda fields, ix=header.index("peptide"): (fields[0], fields[ix]))


def warn_if_psmids_differ(header_base, rows_base, header_new, rows_new):
    """Called when every peptide occurs once in both files. Peptide-level and
    PSM-level files then look alike; the level matters only for peptides
    whose PSMId differs between the files, so name their number."""
    def psmid_of(header, rows):
        ix = header.index("peptide")
        return {fields[ix]: fields[0] for _, fields in rows}
    base, new = psmid_of(header_base, rows_base), psmid_of(header_new, rows_new)
    common = base.keys() & new.keys()
    moved = sorted(p for p in common if base[p] != new[p])
    if moved:
        print("warning: every peptide occurs once in both files, so the rows are "
              "joined on the peptide, as for peptide-level results; %d of %d "
              "common peptides have a different PSMId, e.g. %r: %r and %r. If "
              "these are PSM-level results, use --key psm."
              % (len(moved), len(common), moved[0], base[moved[0]], new[moved[0]]),
              file=sys.stderr)


def index_rows(path, rows, key, score_ix=None):
    """Map each key to its row with the highest score (the first of them on
    ties, or the first row if there is no score column)."""
    def score(fields):
        try:
            value = float(fields[score_ix])
        except (IndexError, ValueError):
            return -math.inf
        return value if not math.isnan(value) else -math.inf
    table, line_of = {}, {}
    repeated, example = 0, None
    for line_nr, fields in rows:
        k = key(fields)
        if k in table:
            repeated += 1
            if example is None:
                example = (k, line_of[k], line_nr)
            if score_ix is not None and score(fields) > score(table[k]):
                table[k] = fields
                line_of[k] = line_nr
            continue
        table[k] = fields
        line_of[k] = line_nr
    if repeated:
        print("warning: %s: %d rows repeat an earlier key; the %s row of each key "
              "is used, e.g. %r on lines %d and %d"
              % ((path, repeated, "best-scoring" if score_ix is not None else "first")
                 + example), file=sys.stderr)
    return table


def same_bytes(path1, path2, block_size=1 << 20):
    with open(path1, "rb") as f1, open(path2, "rb") as f2:
        while True:
            b1, b2 = f1.read(block_size), f2.read(block_size)
            if b1 != b2:
                return False
            if not b1:
                return True


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


def max_abs_diff(base, new, common, cb, cn):
    """Largest |base - new| of columns cb and cn over the common rows; inf if
    any difference is not finite (NaN or inf), 0 if there are no rows."""
    diffs = (abs(float(base[k][cb]) - float(new[k][cn])) for k in common)
    return max((d if math.isfinite(d) else math.inf for d in diffs), default=0.0)


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
    ap.add_argument("--key", choices=["auto", "first", "psm", "peptide"], default="auto",
                    help="join rows on the first column, on PSMId and peptide, or on the "
                         "peptide; auto chooses as described above")
    args = ap.parse_args()

    if args.exact:
        same = same_bytes(args.base, args.new)
        print("identical" if same else "DIFFERENT", args.base, args.new)
        sys.exit(0 if same else 1)

    hb, rows_b = read_rows(args.base)
    hn, rows_n = read_rows(args.new)
    key_name, key_for = key_function(args.key, hb, rows_b, hn, rows_n)
    base = index_rows(args.base, rows_b, key_for(hb), hb.index("score") if "score" in hb else None)
    new = index_rows(args.new, rows_n, key_for(hn), hn.index("score") if "score" in hn else None)
    del rows_b, rows_n
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
    max_dq = max_abs_diff(base, new, common, qb, qn)
    max_dpep = max_abs_diff(base, new, common, pb, pn)
    rho = math.nan  # undefined without a score column (e.g. protein level)
    if sb is not None and sn is not None:
        rho = spearman([float(base[k][sb]) for k in common], [float(new[k][sn]) for k in common])

    print("rows base=%d new=%d common=%d key=%s" % (len(base), len(new), len(common), key_name))
    print("accepted (q<%g) base=%d new=%d rel.change=%.5f jaccard=%.5f" %
          (args.fdr, len(acc_b), len(acc_n), count_change, jaccard))
    print("max|dq|=%.3g max|dPEP|=%.3g spearman(score)=%.6f" % (max_dq, max_dpep, rho))

    ok = (count_change <= args.max_count_change and jaccard >= args.min_jaccard
          and max_dq <= args.max_dq and (args.max_dpep is None or max_dpep <= args.max_dpep))
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
