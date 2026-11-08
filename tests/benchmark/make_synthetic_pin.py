#!/usr/bin/env python3
"""Write a synthetic Percolator input (pin) file for speed and regression tests.

Uses only the Python standard library. The output is deterministic for a given
seed. There is no DefaultDirection row, so Percolator searches for an initial
direction as it would on real data.

Each scan has a correct target match with probability --correct-fraction. The
features of a correct match are drawn around feature-specific shifts (scaled
by --shift-scale), those of an incorrect match or of a decoy around 0.
Correct matches are to peptides from one half of the target peptides (the
peptides present in the sample), incorrect target matches to the other half,
and decoys to as many decoy peptides, so that a target peptide is correct
exactly when its PSMs are, and incorrect target peptides and decoy peptides
are alike, as peptide-level target-decoy competition assumes.

Search modes:
  separate      one target and one decoy PSM per scan (mix-max path); the
                target is correct with probability --correct-fraction,
                otherwise it is an incorrect match drawn like the decoy
  native        like separate, but a scan with a correct match (a native
                spectrum) also has an incorrect target candidate, and the
                better-scoring of the two is reported, labelled by the
                winner, as a search engine reports its best target match
  concatenated  one PSM per scan (target-decoy competition): a correct target
                with probability --correct-fraction, otherwise an incorrect
                match that is a target or a decoy with probability 1/2 each,
                as target-decoy competition assumes

--legacy reproduces the files of earlier versions (concatenated and separate
modes). Their concatenated mode drew a decoy with probability 1/2 and made
--correct-fraction of the targets correct; incorrect targets were then 35%
of the PSMs and decoys 50%, so the decoys overestimated the incorrect
targets by a factor of 1/0.7 and the observed FDR at q < 0.01 was about
0.007. All matches drew from all peptides, so a target peptide's best PSM
could be incorrect while the peptide had correct PSMs, and peptides with
only incorrect PSMs were fewer than decoy peptides. Use it only to
regenerate earlier files.

--truth writes the ground truth of every PSM: whether it is correct, and its
true PEP, the probability that a target with its features is incorrect under
the generating model (1 for decoys). Ranking the targets by true PEP is the
best possible ranking. For a peptide, the true PEP of its best PSM is close to the peptide's true
PEP at the top of the list, where the other PSMs of a correct peptide almost
surely score lower.

Example:
  make_synthetic_pin.py -n 1000000 -f 30 --mode concatenated -o 1M.pin
"""

import argparse
import math
import random
import sys

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def make_vocabulary(rng, size, prefix):
    peptides = []
    for _ in range(size):
        length = rng.randint(7, 30)
        seq = "".join(rng.choice(AMINO_ACIDS) for _ in range(length - 1))
        seq += rng.choice("KR")
        if rng.random() < 0.3:
            pos = rng.randrange(length)
            seq = seq[:pos + 1] + "[15.99]" + seq[pos + 1:]
        peptides.append((rng.choice("KR-") + "." + seq + "." + rng.choice(AMINO_ACIDS + "-"), length))
    proteins = ["%sprot_%d" % (prefix, i) for i in range(max(1, size // 20))]
    return peptides, proteins


def zipf_index(rng, n):
    # Approximate Zipf(1) draw over [0, n)
    return min(n - 1, int(math.exp(rng.random() * math.log(n + 1))) - 1)


def norm_cdf(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


class TruthModel:
    """Exact posteriors under the generating model.

    A correct match has features x_k = N(shift_k, 1) + x_(k-1) / 2, an
    incorrect match or a decoy the same with shift 0. The log likelihood
    ratio of correct against incorrect is therefore
      llr = sum_k shift_k * (x_k - x_(k-1) / 2) - sum_k shift_k^2 / 2,
    which is N(S/2, S) for a correct and N(-S/2, S) for an incorrect match,
    with S = sum_k shift_k^2. It is computed from the values as written, i.e.
    as Percolator reads them."""

    def __init__(self, shifts, correct_fraction, mode):
        self.shifts = shifts
        self.half = sum(s * s for s in shifts) / 2.0
        self.sd = math.sqrt(2.0 * self.half)
        self.c = correct_fraction
        self.mode = mode
        if mode == "concatenated":
            # among the targets: correct c, incorrect (1 - c) / 2
            correct, incorrect = self.c, (1.0 - self.c) / 2.0
        else:  # separate and legacy: a target is correct with probability c
            correct, incorrect = self.c, 1.0 - self.c
        if correct <= 0.0:
            self.log_prior_odds = -math.inf
        elif incorrect <= 0.0:
            self.log_prior_odds = math.inf
        else:
            self.log_prior_odds = math.log(correct / incorrect)

    def llr(self, feats):
        s, prev = 0.0, 0.0
        for shift, text in zip(self.shifts, feats):
            x = float(text)
            s += shift * (x - 0.5 * prev)
            prev = x
        return s - self.half

    def true_pep(self, llr):
        """Probability that a target with this log likelihood ratio is
        incorrect."""
        if self.mode == "native":
            # incorrect: a foreign spectrum, or a native one whose correct
            # candidate scored lower; correct: a native spectrum whose
            # incorrect candidate scored lower
            if self.sd > 0.0:
                cdf_correct = norm_cdf((llr - self.half) / self.sd)
                cdf_incorrect = norm_cdf((llr + self.half) / self.sd)
            else:
                cdf_correct = cdf_incorrect = 0.5
            incorrect = (1.0 - self.c) + self.c * cdf_correct
            if self.c <= 0.0 or cdf_incorrect <= 0.0:
                return 1.0
            if incorrect <= 0.0:
                return 0.0
            z = math.log(self.c * cdf_incorrect / incorrect) + llr
        else:
            z = self.log_prior_odds + llr
        if z > 700.0:
            return 0.0
        if z < -700.0:
            return 1.0
        return 1.0 / (1.0 + math.exp(z))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--num-psms", type=int, required=True,
                    help="number of PSMs; must be even with --mode separate or native")
    ap.add_argument("-f", "--num-features", type=int, default=30,
                    help="total number of feature columns (at least 10)")
    ap.add_argument("--mode", choices=["separate", "native", "concatenated"],
                    default="concatenated")
    ap.add_argument("--legacy", action="store_true",
                    help="draw as earlier versions did, which does not follow the "
                         "target-decoy competition assumption; only to regenerate earlier "
                         "files (not with --mode native)")
    ap.add_argument("--correct-fraction", type=float, default=0.3,
                    help="fraction of the scans that have a correct target match "
                         "(with --legacy --mode concatenated: fraction of the target PSMs "
                         "that are correct)")
    ap.add_argument("--shift-scale", type=float, default=1.0,
                    help="scale of the feature shifts of correct matches; smaller values "
                         "give data that are harder to separate (default 1)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--decoy-prefix", default="decoy_")
    ap.add_argument("--mark-correct", action="store_true",
                    help="mark correct and incorrect targets in SpecId (target_c_..., "
                         "target_i_...) for checking PEP and q-value calibration")
    ap.add_argument("--truth", metavar="FILE",
                    help="write the SpecId, label, correctness and true PEP of every PSM "
                         "to FILE")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()

    if args.num_features < 10:
        ap.error("--num-features must be at least 10")
    if args.mode in ("separate", "native") and args.num_psms % 2:
        ap.error("--num-psms must be even with --mode %s" % args.mode)
    if args.legacy and args.mode == "native":
        ap.error("--legacy cannot be used with --mode native")
    if not 0.0 <= args.correct_fraction <= 1.0:
        ap.error("--correct-fraction must be between 0 and 1")
    rng = random.Random(args.seed)
    num_cont = args.num_features - 7  # continuous; the rest are charge, enzyme and length
    shifts = [args.shift_scale * rng.uniform(0.0, 2.0) if i % 3 else 0.0 for i in range(num_cont)]
    model = TruthModel(shifts, args.correct_fraction,
                       "legacy" if args.legacy else args.mode)

    vocab_size = max(100, args.num_psms // 4)
    target_peps, target_prots = make_vocabulary(rng, vocab_size, "")
    decoy_peps, decoy_prots = make_vocabulary(rng, vocab_size, args.decoy_prefix)

    header = ["SpecId", "Label", "ScanNr", "ExpMass"]
    header += ["feat%d" % i for i in range(num_cont)]
    header += ["Charge2", "Charge3", "Charge4", "enzN", "enzC", "enzInt", "PepLen", "Peptide", "Proteins"]

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write("\t".join(header) + "\n")
    truth = open(args.truth, "w") if args.truth else None
    if truth:
        truth.write("SpecId\tLabel\tcorrect\ttrue_pep\n")

    def draw_features(correct):
        feats = []
        prev = 0.0
        for i in range(num_cont):
            value = rng.gauss(shifts[i] if correct else 0.0, 1.0) + 0.5 * prev
            feats.append("%.5f" % value)
            prev = value
        return feats

    def write_psm(scan, exp_mass, charge, is_decoy, correct=None, feats=None):
        if correct is None:
            correct = (not is_decoy) and rng.random() < args.correct_fraction
        peps, prots = (decoy_peps, decoy_prots) if is_decoy else (target_peps, target_prots)
        if args.legacy:
            index = zipf_index(rng, len(peps))
        else:
            # present target peptides, and decoy peptides, from the first
            # half; absent target peptides from the second half
            half = len(peps) // 2
            index = zipf_index(rng, half)
            if not is_decoy and not correct:
                index += half
        peptide, length = peps[index]
        if feats is None:
            feats = draw_features(correct)
        true_pep = 1.0 if is_decoy or not truth else model.true_pep(model.llr(feats))
        feats = feats + ["1" if charge == c else "0" for c in (2, 3, 4)]
        feats += [str(int(rng.random() < 0.9)), str(int(rng.random() < 0.9)), str(rng.choice((0, 0, 0, 1, 2)))]
        feats.append(str(length))
        num_prots = rng.choice((1, 1, 1, 2, 3))
        proteins = [prots[zipf_index(rng, len(prots))] for _ in range(num_prots)]
        label = "-1" if is_decoy else "1"
        kind = "decoy" if is_decoy else "target"
        if args.mark_correct and not is_decoy:
            kind += "_c" if correct else "_i"
        spec_id = "%s_%d_%d_%s" % (kind, scan, charge, label)
        out.write("\t".join([spec_id, label, str(scan), "%.4f" % exp_mass] + feats + [peptide] + proteins) + "\n")
        if truth:
            truth.write("%s\t%s\t%d\t%.6g\n" % (spec_id, label, int(correct), true_pep))

    scan = 0
    written = 0
    while written < args.num_psms:
        scan += 1
        exp_mass = rng.uniform(600.0, 4000.0)
        charge = rng.choice((2, 2, 3, 3, 4))
        if args.mode == "separate":
            write_psm(scan, exp_mass, charge, False)
            write_psm(scan, exp_mass, charge, True)
            written += 2
        elif args.mode == "native":
            if rng.random() < args.correct_fraction:
                # a native spectrum: the better of its correct and an
                # incorrect target candidate is reported
                correct_feats = draw_features(True)
                incorrect_feats = draw_features(False)
                correct = model.llr(correct_feats) > model.llr(incorrect_feats)
                feats = correct_feats if correct else incorrect_feats
            else:
                correct, feats = False, draw_features(False)
            write_psm(scan, exp_mass, charge, False, correct, feats)
            write_psm(scan, exp_mass, charge, True, False)
            written += 2
        elif args.legacy:
            write_psm(scan, exp_mass, charge, rng.random() < 0.5)
            written += 1
        else:
            # target-decoy competition: an incorrect match is a target or a
            # decoy with equal probability, a correct match is a target
            if rng.random() < args.correct_fraction:
                write_psm(scan, exp_mass, charge, False, True)
            else:
                write_psm(scan, exp_mass, charge, rng.random() < 0.5, False)
            written += 1
    if out is not sys.stdout:
        out.close()
    if truth:
        truth.close()


if __name__ == "__main__":
    main()
