#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


def update_block_value(text: str, block_name: str, key: str, value_expr: str) -> str:
    pattern = re.compile(rf"(?is)(&\s*{re.escape(block_name)}\b)(.*?)(/)")
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"Could not find &{block_name} block.")

    start, body, end = match.group(1), match.group(2), match.group(3)
    key_pattern = re.compile(rf"(?im)^\s*{re.escape(key)}\s*=.*$")
    new_line = f"{key}={value_expr}"

    if key_pattern.search(body):
        new_body = key_pattern.sub(new_line, body)
    else:
        if body and not body.endswith("\n"):
            body += "\n"
        new_body = body + f"{new_line}\n"

    return text[: match.start()] + start + new_body + end + text[match.end() :]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a halo smoke-test namelist with nstepmax=1 for analysis pipeline checks."
    )
    parser.add_argument("--input", type=Path, required=True, help="Source namelist.")
    parser.add_argument("--output", type=Path, required=True, help="Output namelist path.")
    parser.add_argument(
        "--foutput",
        type=int,
        default=1,
        help="Output cadence in coarse steps (default: 1).",
    )
    parser.add_argument(
        "--ncontrol",
        type=int,
        default=1,
        help="Print cadence in coarse steps (default: 1).",
    )
    args = parser.parse_args()

    src = args.input.resolve()
    dst = args.output.resolve()
    text = src.read_text()

    text = update_block_value(text, "RUN_PARAMS", "nstepmax", "1")
    text = update_block_value(text, "RUN_PARAMS", "ncontrol", str(args.ncontrol))
    text = update_block_value(text, "OUTPUT_PARAMS", "foutput", str(args.foutput))

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text)
    print(f"Wrote {dst}")


if __name__ == "__main__":
    main()

