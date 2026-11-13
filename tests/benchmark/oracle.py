#!/usr/bin/env python3
"""Best possible numbers of correct targets for a synthetic input made with
make_synthetic_pin.py --truth, and the numbers a Percolator ranking reaches.

The oracle counts, along a ranking of the targets, the largest number of
correct targets in a top list whose true FDP (fraction of incorrect targets)
is at most 1% and 5%. Ranking the targets by their true PEP is the best
possible ranking (the Bayes oracle); each result file is ranked by its
scores. Unlike accepted counts at a q-value threshold, these numbers do not
depend on the FDR estimate, only on the ranking.

  oracle.py truth.txt [result.txt ...]
"""

import sys

LEVELS = (0.01, 0.05)


def oracle_counts(correct_flags):
    best = {t: 0 for t in LEVELS}
    incorrect = 0
    for k, correct in enumerate(correct_flags, 1):
        if not correct:
            incorrect += 1
        for t in LEVELS:
            if incorrect <= t * k:
                best[t] = max(best[t], k - incorrect)
    return best


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    truth = {}
    with open(sys.argv[1]) as f:
        header = f.readline().rstrip("\n").split("\t")
        id_col, label_col, correct_col, pep_col = (
            header.index(c) for c in ("SpecId", "Label", "correct", "true_pep"))
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[label_col] == "1":
                truth[fields[id_col]] = (fields[correct_col] == "1", float(fields[pep_col]))
    ranked = sorted(truth.values(), key=lambda v: v[1])
    best = oracle_counts(c for c, _ in ranked)
    print("%s: Bayes oracle: %d correct targets at true FDP <= 1%%, %d at <= 5%% "
          "(%d targets, %d correct)" % (sys.argv[1], best[0.01], best[0.05], len(ranked),
                                        sum(1 for c, _ in ranked if c)))
    for path in sys.argv[2:]:
        rows, seen = [], set()
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            score_col = header.index("score")
            for line in f:
                fields = line.rstrip("\n").split("\t")
                key = fields[0]
                if key in truth and key not in seen:
                    seen.add(key)
                    rows.append((float(fields[score_col]), truth[key][0]))
        rows.sort(key=lambda r: -r[0])
        best = oracle_counts(c for _, c in rows)
        print("%s: %d correct targets at true FDP <= 1%%, %d at <= 5%% (%d targets)"
              % (path, best[0.01], best[0.05], len(rows)))


if __name__ == "__main__":
    main()
