#!/usr/bin/env python3

import argparse
import os
import shutil
from pathlib import Path


DEFAULT_SHIFT = 150.0 + 300.0 / 512.0


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Copy an initial-condition folder and shift the x coordinate in "
            "its ic_part file."
        )
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("nfw"),
        help="Source IC directory to copy. Default: nfw",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path("nfw_shifted"),
        help="Destination directory to create. Default: nfw_shifted",
    )
    parser.add_argument(
        "--ic-name",
        default="ic_part",
        help="Particle IC filename inside source/dest. Default: ic_part",
    )
    parser.add_argument(
        "--x-shift",
        type=float,
        default=DEFAULT_SHIFT,
        help=f"Value added to the x column. Default: {DEFAULT_SHIFT}",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace destination directory if it already exists.",
    )
    return parser.parse_args()


def copy_source_tree(source, dest, overwrite=False):
    if not source.is_dir():
        raise FileNotFoundError(f"Source directory does not exist: {source}")

    if dest.exists():
        if not overwrite:
            raise FileExistsError(
                f"Destination already exists: {dest}. Use --overwrite to replace it."
            )
        shutil.rmtree(dest)

    shutil.copytree(source, dest)


def shift_ic_part(ic_file, x_shift):
    tmp_file = ic_file.with_name(f"{ic_file.name}.tmp")

    with ic_file.open("r") as src, tmp_file.open("w") as dst:
        for line_number, line in enumerate(src, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                dst.write(line)
                continue

            columns = stripped.split()
            if len(columns) < 1:
                dst.write(line)
                continue

            try:
                columns[0] = f"{float(columns[0]) + x_shift:.10g}"
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse x coordinate on line {line_number} of {ic_file}"
                ) from exc

            dst.write(" ".join(columns) + "\n")

    os.replace(tmp_file, ic_file)


def main():
    args = parse_args()
    copy_source_tree(args.source, args.dest, overwrite=args.overwrite)
    shift_ic_part(args.dest / args.ic_name, args.x_shift)

    print(f"Created {args.dest}")
    print(f"Shifted {args.dest / args.ic_name} x column by {args.x_shift:.10g}")


if __name__ == "__main__":
    main()
