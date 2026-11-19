# Conda package

The bioconda package is built by `admin/builders/conda_build.sh`. It links
dynamically against OpenMP of the conda environment (`llvm-openmp` on macOS,
the compiler's `libgomp` on Linux) and installs `percolator` and `qvality`
into `$PREFIX`. Percolator needs no other library.

`meta.yaml` and `build.sh` in this directory are a test recipe that follows
the bioconda recipe but builds the checked-out source.
`.github/workflows/conda.yml` builds and tests it with the conda-forge
compilers and pinnings on linux-64, osx-arm64 (`macos-15`) and osx-64
(`macos-15-intel`). To run it locally:

```bash
conda install -c conda-forge conda-build conda-forge-pinning
conda build admin/conda -c conda-forge -m "$CONDA_PREFIX/conda_build_config.yaml"
```

On macOS, also pass `-m` with a file that sets `CONDA_BUILD_SYSROOT` to
the output of `xcrun --show-sdk-path`, as the workflow does.

## Updating bioconda after a release

1. Check that the conda workflow passes on the release commit.
2. Tag the release (`rel-<major>-<minor>[-<patch>]`, e.g. `rel-4-00`) and
   create the GitHub release.
3. Compute the checksum of the source archive that bioconda downloads:

   ```bash
   curl -sL https://github.com/percolator/percolator/archive/refs/tags/rel-4-00.tar.gz | sha256sum
   ```

4. In a fork of bioconda-recipes, update `recipes/percolator`: the version
   (`4.0.0` for `rel-4-00`), the source URL and sha256, and reset the build
   number to 0. Copy `build.sh` from this directory, and keep the
   requirements and tests in step with `meta.yaml` here.
5. Open a pull request against bioconda-recipes. The bioconda CI builds all
   platforms listed in the recipe.

### Changes to the bioconda recipe for 3.9.1

Compared with the recipe for 3.9:

- Remove `skip: True  # [osx]` and add `osx-arm64` to
  `extra: additional-platforms`.
- Add `git` to the build requirements (Eigen is downloaded during the
  configure step).
- Host requirements: keep `libboost-devel`, add `llvm-openmp  # [osx]`, and
  remove `sqlite`, `bzip2` and `zlib`, which Percolator does not use. Run
  requirements: `llvm-openmp  # [osx]`; `libboost` comes from the run
  exports of `libboost-devel`.
- Replace `build.sh` with the one in this directory.

### Changes to the bioconda recipe for 4.0.0

Compared with the recipe for 3.9.1 (apply the changes for 3.9.1 first if
bioconda still has 3.9):

- Remove `libboost-devel` from the host requirements: Percolator uses
  `std::filesystem` instead of Boost. The run requirement `libboost` came
  from its run exports and goes with it.
- Remove `git` from the build requirements: Eigen is no longer used, and
  nothing is downloaded during the configure step.
- The requirements are then those of `meta.yaml` here: the C++ compiler,
  `stdlib('c')`, `cmake` and `make` to build, and `llvm-openmp  # [osx]` as
  host and run requirement.
- `std::filesystem` needs macOS 10.15 or later; the conda-forge pinnings set
  11.0 for both macOS platforms (November 2026), so the recipe needs no
  setting of its own. CMake keeps the deployment target that conda sets
  (`MACOSX_DEPLOYMENT_TARGET`); only without one does it build for macOS
  12.0.

### Changes to the bioconda recipe for 4.0.1

None besides the version (`4.0.1` for the tag `rel-4-00-01`), the source
URL and sha256, and the build number 0: the requirements and `build.sh`
are those of 4.0.0.
