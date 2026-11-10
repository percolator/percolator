# CI scripts

The scripts in this directory are run by `.github/workflows/build.yml` and
can be run by hand on a build tree. They need bash, awk and Python 3 (for
`tests/benchmark/make_synthetic_pin.py`), nothing else.

| Script | What it checks |
|---|---|
| `option_matrix.sh <build-dir> <out-dir> [percolator options]` | Runs percolator and qvality over a matrix of option combinations; every run must exit with status 0 and print no sanitizer report. |
| `check_golden.sh [--update \| --report \| --strict] <build-dir> <out-dir>` | Runs the option matrix on inputs of 8,000 PSMs and compares the SHA-256 hash of every output file with `option_matrix.sha256`. |
| `check_determinism.sh <build-dir> <out-dir>` | Runs the option matrix with 1 and with 4 threads; the outputs must be byte-identical. |
| `count_warnings.sh <build-log> <budget>` | Counts the distinct compiler warnings of a build log. |
| `check_symbol_versions.sh` | The glibc symbol versions of the portable Linux binaries. |

## The manifest of output hashes

`option_matrix.sha256` records the SHA-256 hash of each of the output files
of `option_matrix.sh` at `MATRIX_PSMS=8000` (result files, weights, standard
output; not stderr) and of its generated inputs. A change that is meant to
keep the results byte for byte is checked with

    bash tests/ci/check_golden.sh <build-dir> <out-dir>

which takes about 15 s and prints `291 files, all as in the manifest`, or
the files that differ (`DIFF`), are no longer written (`MISSING`) or are not
in the manifest (`NEW`). The check is strict (exit status 1) on x86-64 Linux,
where the manifest was made; elsewhere, and with `--report`, differences are
listed without failing, since another processor architecture or math library
may round a last digit differently. A difference in `inputs/` means that the
input generator, not percolator, gives other files.

A commit that changes results on purpose, or that adds, removes or renames
runs in `option_matrix.sh`, writes the manifest again in the same commit:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j 4
    bash tests/ci/check_golden.sh --update build build/golden
    git diff --stat tests/ci/option_matrix.sha256

Use a Release build at the default flags on x86-64 Linux (the `build-ubuntu`
job compares with it), and name in the commit message which runs change and
why. `git diff tests/ci/option_matrix.sha256` shows exactly which files
changed: a result-preserving change must leave it empty.

The manifest holds for inputs of 8,000 PSMs and the default compiler flags.
Larger inputs and other flags are deliberately not covered: at 1M PSMs a
build with FMA instructions was found to differ from the default
build in the last digit of a PEP.
