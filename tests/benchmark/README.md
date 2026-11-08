# Benchmark helpers

Scripts for timing Percolator on large inputs and for checking that a code
change leaves the results unchanged. They need only the Python standard library.

- `make_synthetic_pin.py` writes a deterministic synthetic pin file of any size
  (an even size for separate searches, which pair a target and a decoy per
  spectrum), for separate or concatenated target-decoy searches.
  - In concatenated mode an incorrect match is a target or a decoy with
    probability 1/2 each, as target-decoy competition assumes, so the decoys
    estimate the incorrect targets without bias.
  - Correct and incorrect target matches are to different peptides, so a
    peptide is correct exactly when its PSMs are, and incorrect target
    peptides and decoy peptides are alike.
  - `--mode native` is a separate search in which a spectrum with a correct
    match also has an incorrect target candidate and reports the better one.
  - `--shift-scale` below 1 makes the data harder to separate.
  - `--truth FILE` writes whether each PSM is correct and its true PEP under
    the generating model; it does not change the pin file.
  - `--legacy` reproduces files of earlier versions (e.g. a `1M_marked.pin`
    made with `-n 1000000 --mark-correct`). In their concatenated mode the
    decoys are 50% of the PSMs and the incorrect targets 35%, so Percolator's
    FDR and PEPs look conservative by a factor of 0.7, and correct and
    incorrect PSMs share peptides, which biases peptide-level checks; do not
    use them to judge PEPs or FDR control.
- `compare_results.py` compares two result files, either byte for byte
  (`--exact`) or within tolerances on the accepted targets and the q-values.
  The largest PEP difference is reported, but only checked when `--max-dpep`
  is given. Peptide-level rows are matched on the peptide, PSM-level rows on
  the PSMId and the peptide, and protein-level rows on the first column
  (`--key` overrides the choice); when a file repeats a key, its first
  (best-scoring) row is used and a warning gives the count.
- `pep_calibration.py` checks PEP and q-value calibration of a PSM result
  file for an input made with `make_synthetic_pin.py --mark-correct`, which
  marks correct and incorrect targets in the PSM ids.

Per-stage wall clock and CPU times are printed to stderr with `-v 3`
(lines starting with `Stage`).

Example:

    python3 make_synthetic_pin.py -n 1000000 --mode concatenated -o 1M.pin
    percolator -S 1 -v 3 -m new_psms.txt 1M.pin
    python3 compare_results.py --exact old_psms.txt new_psms.txt
