#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import argparse
from collections import defaultdict

def merge_groups(in_dir, middle_name):
    pattern = re.compile(r"^(.*?)_mpi\d+_(.*?)\.out$")
    files = glob.glob(os.path.join(in_dir, "*.out"))

    # Group files by (prefix, suffix)
    groups = defaultdict(list)

    for f in files:
        fname = os.path.basename(f)
        m = pattern.match(fname)
        if m:
            prefix = m.group(1)
            suffix = m.group(2)
            groups[(prefix, suffix)].append(f)

    if not groups:
        raise RuntimeError("No files matching <prefix>_mpi<number>_<suffix>.out")

    for (prefix, suffix), flist in groups.items():
        flist = sorted(flist)

        out_file = os.path.join(in_dir, f"{prefix}_{middle_name}_{suffix}.out")

        # Remove existing merged file
        if os.path.exists(out_file):
            print(f"Removing old file: {out_file}")
            os.remove(out_file)

        # Merge
        with open(out_file, "w") as fout:
            for fname in flist:
                with open(fname) as fin:
                    fout.write(fin.read())

        print(f"Merged {len(flist)} files → {out_file}")

        # Delete originals
        for f in flist:
            print(f"Removing {f}")
            os.remove(f)

    print("All merges complete.")


def main():
    parser = argparse.ArgumentParser(description="Merge all <prefix>_mpiX_<suffix>.out groups")
    parser.add_argument("in_dir", help="Directory containing .out files")
    parser.add_argument("middle_name", help="Middle part of merged filename")
    args = parser.parse_args()

    merge_groups(args.in_dir, args.middle_name)


if __name__ == "__main__":
    main()
