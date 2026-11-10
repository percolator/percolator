#!/usr/bin/env python3
"""Hashes of the outputs of tests/ci/option_matrix.sh, to compare platforms.

  matrix_hashes.py hash <out-dir> > hashes.txt
      One line "<sha256>  <file>" for every file below <out-dir>, sorted,
      without the stderr.txt files (they hold timings, the version and the
      command line). Line ends are hashed as "\\n", so that Windows gives
      the hashes of the other platforms. The inputs are listed too
      (inputs/...): where they differ, the Python that generated them does,
      not Percolator.

  matrix_hashes.py compare [--lines N] <reference> <other>...
      Each argument is a hash list, or a directory that holds hashes.txt
      and, optionally, the output files beside it. Prints for every <other>
      the files whose hashes differ from the reference's and the files that
      only one of the two has. Where both directories hold a differing
      file, its first N differing lines are printed as well (default 5).
      Differences are a report, not an error: the exit status is 0 unless
      an argument cannot be read. Compilers that contract a*b+c into a
      fused multiply-add (the default on arm64) round differently, which
      changes a printed digit now and then.

Python standard library only.
"""

import hashlib
import os
import sys


def file_hash(path):
    with open(path, "rb") as handle:
        data = handle.read()
    return hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def hash_tree(root):
    lines = []
    for directory, _, files in os.walk(root):
        for name in files:
            path = os.path.join(directory, name)
            relative = os.path.relpath(path, root).replace(os.sep, "/")
            # hashes.txt: this list itself, when it is written into <out-dir>
            if name == "stderr.txt" or relative == "hashes.txt":
                continue
            lines.append((relative, file_hash(path)))
    return sorted(lines)


def read_list(argument):
    """Returns the name, the hashes by file and the directory of the outputs
    (None for a plain hash list)."""
    directory = None
    path = argument
    if os.path.isdir(argument):
        directory = argument
        path = os.path.join(argument, "hashes.txt")
    hashes = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                digest, name = line.rstrip("\r\n").split("  ", 1)
                hashes[name] = digest
    return argument, hashes, directory


def read_lines(path):
    with open(path, "rb") as handle:
        return handle.read().replace(b"\r\n", b"\n").split(b"\n")


def print_line_differences(reference_file, other_file, limit):
    reference = read_lines(reference_file)
    other = read_lines(other_file)
    if len(reference) != len(other):
        print("      %d lines against %d" % (len(reference), len(other)))
    differing = [i for i in range(min(len(reference), len(other)))
                 if reference[i] != other[i]]
    print("      %d of %d lines differ" % (len(differing), len(reference)))
    for i in differing[:limit]:
        print("      line %d" % (i + 1))
        print("        < " + reference[i].decode("utf-8", "replace"))
        print("        > " + other[i].decode("utf-8", "replace"))


def compare(arguments, limit):
    reference_name, reference, reference_dir = read_list(arguments[0])
    print("reference: %s (%d files)" % (reference_name, len(reference)))
    for argument in arguments[1:]:
        name, hashes, directory = read_list(argument)
        differing = sorted(f for f in reference
                           if f in hashes and hashes[f] != reference[f])
        missing = sorted(f for f in reference if f not in hashes)
        extra = sorted(f for f in hashes if f not in reference)
        print("%s: %d files, %d differ, %d missing, %d not in the reference"
              % (name, len(hashes), len(differing), len(missing), len(extra)))
        for f in differing:
            print("  DIFF %s" % f)
            if reference_dir and directory:
                pair = [os.path.join(d, *f.split("/"))
                        for d in (reference_dir, directory)]
                if all(os.path.isfile(p) for p in pair):
                    print_line_differences(pair[0], pair[1], limit)
        for f in missing:
            print("  MISSING %s" % f)
        for f in extra:
            print("  EXTRA %s" % f)


def main(argv):
    if len(argv) == 3 and argv[1] == "hash":
        if not os.path.isdir(argv[2]):
            sys.exit("error: %s is not a directory" % argv[2])
        for name, digest in hash_tree(argv[2]):
            sys.stdout.write("%s  %s\n" % (digest, name))
        return
    if len(argv) >= 3 and argv[1] == "compare":
        arguments = argv[2:]
        limit = 5
        if arguments[0] == "--lines":
            limit = int(arguments[1])
            arguments = arguments[2:]
        if arguments:
            compare(arguments, limit)
            return
    sys.stderr.write(__doc__)
    sys.exit(2)


if __name__ == "__main__":
    main(sys.argv)
