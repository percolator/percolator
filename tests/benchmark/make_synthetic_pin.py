#!/usr/bin/env python3
"""Write a synthetic Percolator input (pin) file for speed and regression tests.

Uses only the Python standard library. The output is deterministic for a given
seed. There is no DefaultDirection row, so Percolator searches for an initial
direction as it would on real data.

Search modes:
  separate      one target and one decoy PSM per scan (mix-max path)
  concatenated  one PSM per scan, target or decoy (target-decoy competition)

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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--num-psms", type=int, required=True)
    ap.add_argument("-f", "--num-features", type=int, default=30,
                    help="total number of feature columns (at least 10)")
    ap.add_argument("--mode", choices=["separate", "concatenated"], default="concatenated")
    ap.add_argument("--correct-fraction", type=float, default=0.3,
                    help="fraction of target PSMs that are correct")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--decoy-prefix", default="decoy_")
    ap.add_argument("--mark-correct", action="store_true",
                    help="mark correct and incorrect targets in SpecId (target_c_..., "
                         "target_i_...) for checking PEP and q-value calibration")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()

    if args.num_features < 10:
        ap.error("--num-features must be at least 10")
    if args.mode == "separate" and args.num_psms % 2:
        ap.error("--num-psms must be even with --mode separate")
    rng = random.Random(args.seed)
    num_cont = args.num_features - 7  # continuous; the rest are charge, enzyme and length
    shifts = [rng.uniform(0.0, 2.0) if i % 3 else 0.0 for i in range(num_cont)]

    vocab_size = max(100, args.num_psms // 4)
    target_peps, target_prots = make_vocabulary(rng, vocab_size, "")
    decoy_peps, decoy_prots = make_vocabulary(rng, vocab_size, args.decoy_prefix)

    header = ["SpecId", "Label", "ScanNr", "ExpMass"]
    header += ["feat%d" % i for i in range(num_cont)]
    header += ["Charge2", "Charge3", "Charge4", "enzN", "enzC", "enzInt", "PepLen", "Peptide", "Proteins"]

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write("\t".join(header) + "\n")

    def write_psm(scan, exp_mass, charge, is_decoy):
        correct = (not is_decoy) and rng.random() < args.correct_fraction
        peps, prots = (decoy_peps, decoy_prots) if is_decoy else (target_peps, target_prots)
        peptide, length = peps[zipf_index(rng, len(peps))]
        feats = []
        prev = 0.0
        for i in range(num_cont):
            value = rng.gauss(shifts[i] if correct else 0.0, 1.0) + 0.5 * prev
            feats.append("%.5f" % value)
            prev = value
        feats += ["1" if charge == c else "0" for c in (2, 3, 4)]
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
        else:
            write_psm(scan, exp_mass, charge, rng.random() < 0.5)
            written += 1
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
