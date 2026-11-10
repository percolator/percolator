#!/usr/bin/env python3
"""Checks of the binaries-only wheel (pyproject.toml).

  check_wheel.py --wheel <file.whl>...
      The contents of each wheel: the launcher package with bin/percolator
      and bin/qvality (.exe on Windows, with vcomp140.dll beside them), the
      two console commands, a platform tag that is not tied to a Python
      version, and no test binary.

  check_wheel.py --installed <pin-file>
      The installed wheel: runs the percolator and qvality commands found
      on the PATH, which must be the launchers of the wheel, on <pin-file>
      (data/percolator/tab/percolatorTab) and on the scores it gives, and
      checks the exit status and that results were written. cibuildwheel
      runs this as its test command on every platform.

Python standard library only. Exits with status 1 at the first failure.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

PACKAGE = "percolator_launcher"


def fail(message):
    sys.exit("check_wheel: FAILED: " + message)


def check_wheel(path):
    name = os.path.basename(path)
    parts = name[:-len(".whl")].split("-")
    if len(parts) != 5:
        fail("%s: unexpected file name" % name)
    python_tag, abi_tag, platform_tag = parts[2:]
    if (python_tag, abi_tag) != ("py3", "none") or platform_tag == "any":
        fail("%s: the tags must be py3-none-<platform>" % name)
    windows = platform_tag.startswith("win")
    suffix = ".exe" if windows else ""
    with zipfile.ZipFile(path) as wheel:
        names = wheel.namelist()
        entry_points = [n for n in names
                        if n.endswith(".dist-info/entry_points.txt")]
        if len(entry_points) != 1:
            fail("%s: no entry_points.txt" % name)
        scripts = wheel.read(entry_points[0]).decode("utf-8")
        modes = dict((info.filename, info.external_attr >> 16)
                     for info in wheel.infolist())
    required = [PACKAGE + "/__init__.py"]
    for program in ("percolator", "qvality"):
        required.append("%s/bin/%s%s" % (PACKAGE, program, suffix))
        if "%s = %s:%s" % (program, PACKAGE, program) not in scripts:
            fail("%s: no console command %s" % (name, program))
    if windows:
        # OpenMP of MSVC has no static library; the DLL must be beside the
        # executables
        required.append(PACKAGE + "/bin/vcomp140.dll")
    for member in required:
        if member not in names:
            fail("%s: %s is missing" % (name, member))
    if not windows:
        for member in required[1:]:
            if not modes[member] & 0o111:
                fail("%s: %s is not executable" % (name, member))
    for member in names:
        if "gtest" in member.lower():
            fail("%s: holds the test binary %s" % (name, member))
    print("check_wheel: %s ok (%d files)" % (name, len(names)))
    for member in sorted(names):
        print("  " + member)


def run(command, **options):
    print("+ " + " ".join(command))
    sys.stdout.flush()
    status = subprocess.call(command, **options)
    if status != 0:
        fail("%s ended with status %d" % (command[0], status))


def data_lines(path):
    with open(path) as handle:
        return [line for line in handle][1:]


def check_installed(pin):
    import percolator_launcher
    package_dir = os.path.dirname(os.path.abspath(percolator_launcher.__file__))
    bin_dir = os.path.join(package_dir, "bin")
    print("check_wheel: bundled files: " + ", ".join(sorted(os.listdir(bin_dir))))
    for program in ("percolator", "qvality"):
        percolator_launcher.executable(program)  # raises if it is missing
        if shutil.which(program) is None:
            fail("the command %s is not on the PATH" % program)
    if os.name == "nt" and not os.path.isfile(os.path.join(bin_dir, "vcomp140.dll")):
        fail("vcomp140.dll is not beside the executables")
    work = tempfile.mkdtemp()
    try:
        psms = os.path.join(work, "psms.txt")
        decoys = os.path.join(work, "decoy_psms.txt")
        peptides = os.path.join(work, "peptides.txt")
        run(["percolator", "--version"])
        run(["percolator", "-S", "1", "-m", psms, "-M", decoys, "-r", peptides, pin])
        for path in (psms, decoys, peptides):
            if len(data_lines(path)) < 100:
                fail("%s holds fewer than 100 results" % path)
        # qvality on the scores percolator gave
        scores = []
        for path in (psms, decoys):
            scores.append(path + ".scores")
            with open(scores[-1], "w") as out:
                for line in data_lines(path):
                    out.write(line.split("\t")[1] + "\n")
        qvalues = os.path.join(work, "qvality.txt")
        with open(qvalues, "w") as out:
            run(["qvality"] + scores, stdout=out)
        if len(data_lines(qvalues)) < 100:
            fail("qvality wrote fewer than 100 lines")
        # the launcher passes the exit status on
        with open(os.devnull, "w") as null:
            status = subprocess.call(["percolator", os.path.join(work, "missing.pin")],
                                     stdout=null, stderr=null)
        if status == 0:
            fail("percolator on a missing input ended with status 0")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("check_wheel: the installed commands work")


def main(argv):
    if len(argv) >= 3 and argv[1] == "--wheel":
        for path in argv[2:]:
            check_wheel(path)
    elif len(argv) == 3 and argv[1] == "--installed":
        check_installed(argv[2])
    else:
        sys.stderr.write(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main(sys.argv)
