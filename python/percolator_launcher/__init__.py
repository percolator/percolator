"""Starts the percolator and qvality executables that the wheel bundles.

The wheel holds the binaries of the CMake build in the bin directory of this
package (see pyproject.toml). The console commands "percolator" and
"qvality" call the functions of the same names below. There is no Python
API: the functions pass the command line on and end with the exit status of
the program.
"""

import os
import subprocess
import sys

__all__ = ["executable", "percolator", "qvality"]


def executable(name):
    """Returns the path of the bundled program (percolator or qvality)."""
    suffix = ".exe" if os.name == "nt" else ""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin",
                        name + suffix)
    if not os.path.isfile(path):
        raise FileNotFoundError("the bundled program %s is missing" % path)
    return path


def _run(name):
    exe = executable(name)
    args = [exe] + sys.argv[1:]
    if os.name != "nt":
        # replaces this process, so signals and the exit status are the
        # program's own
        os.execv(exe, args)
    # Windows has no exec that keeps the console waiting: start the program
    # and wait for it. Ctrl-C reaches the child as well, so keep waiting for
    # it to end instead of leaving it behind. vcomp140.dll (OpenMP) is found
    # in the directory of the executable.
    process = subprocess.Popen(args)
    while True:
        try:
            sys.exit(process.wait())
        except KeyboardInterrupt:
            pass


def percolator():
    _run("percolator")


def qvality():
    _run("qvality")
