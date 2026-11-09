# Benchmark helpers

Scripts for timing Percolator on large inputs and for checking that a code
change leaves the results unchanged. They need only the Python standard library.

- `make_synthetic_pin.py` writes a deterministic synthetic pin file of any size
  (an even size for separate searches, which pair a target and a decoy per
  spectrum), for separate or concatenated target-decoy searches.
- `compare_results.py` compares two result files, either byte for byte
  (`--exact`) or within tolerances on the accepted targets and the q-values.
  The largest PEP difference is reported, but only checked when `--max-dpep`
  is given. Peptide-level rows are matched on the peptide, PSM-level rows on
  the PSMId and the peptide, and protein-level rows on the first column
  (`--key` overrides the choice). A PSM-level file in which every peptide
  occurs once looks like a peptide-level file; it is joined on the peptide,
  with a warning if some peptides changed their PSMId, so use `--key psm`
  for such files. When a file repeats a key, its highest-scoring row is used
  and a warning gives the count.
- `pep_calibration.py` checks PEP and q-value calibration of a PSM result
  file for an input made with `make_synthetic_pin.py --mark-correct`, which
  marks correct and incorrect targets in the PSM ids.

Per-stage wall clock and CPU times are printed to stderr with `-v 3`
(lines starting with `Stage`).

Example:

    python3 make_synthetic_pin.py -n 1000000 --mode concatenated -o 1M.pin
    percolator -S 1 -v 3 -m new_psms.txt 1M.pin
    python3 compare_results.py --exact old_psms.txt new_psms.txt
