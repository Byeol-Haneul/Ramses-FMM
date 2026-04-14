#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

LEVELMIN_RE = re.compile(r"^(\s*levelmin\s*=\s*)([-+]?\d+)(\s*(?:!.*)?)$", re.IGNORECASE)
LEVELMAX_RE = re.compile(r"^(\s*levelmax\s*=\s*)([-+]?\d+)(\s*(?:!.*)?)$", re.IGNORECASE)


def update_level_bounds(text: str, level: int) -> str:
    lines = text.splitlines(keepends=True)
    found_min = False
    found_max = False

    for i, line in enumerate(lines):
        m = LEVELMIN_RE.match(line.rstrip("\n"))
        if m:
            lines[i] = f"{m.group(1)}{level}{m.group(3)}\n"
            found_min = True
            continue

        m = LEVELMAX_RE.match(line.rstrip("\n"))
        if m:
            lines[i] = f"{m.group(1)}{level}{m.group(3)}\n"
            found_max = True

    if not found_min or not found_max:
        raise ValueError("Could not find both levelmin and levelmax in template namelist.")
    return "".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate level-specific spheres namelists by copying a template and "
            "setting levelmin=levelmax=L."
        )
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=Path("namelist/spheres.nml"),
        help="Template namelist path (default: namelist/spheres.nml).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("namelist/generated_spheres"),
        help="Directory for generated namelists.",
    )
    parser.add_argument(
        "--level-start",
        type=int,
        default=6,
        help="First level to generate (inclusive).",
    )
    parser.add_argument(
        "--level-end",
        type=int,
        default=10,
        help="Last level to generate (inclusive).",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="spheres_l",
        help="Output file prefix (default: spheres_l).",
    )
    args = parser.parse_args()

    if args.level_end < args.level_start:
        raise ValueError("--level-end must be >= --level-start.")

    template = args.template.resolve()
    if not template.exists():
        raise FileNotFoundError(f"Template namelist not found: {template}")

    template_text = template.read_text(encoding="utf-8")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for level in range(args.level_start, args.level_end + 1):
        out_path = output_dir / f"{args.prefix}{level:02d}.nml"
        out_path.write_text(update_level_bounds(template_text, level), encoding="utf-8")
        written.append(out_path)

    print(f"Wrote {len(written)} namelists to {output_dir}")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
