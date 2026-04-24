#!/usr/bin/env python3

import argparse
import csv
import os
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage a subset of SWAM-g input FASTQs from a manifest TSV.")
    parser.add_argument("--manifest", required=True, help="Manifest TSV with sample, r1, r2 columns.")
    parser.add_argument("--dest", required=True, help="Destination directory for the subset.")
    parser.add_argument(
        "--mode",
        choices=("symlink", "copy"),
        default="symlink",
        help="How to stage the files into the destination directory.",
    )
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    dest_dir = Path(args.dest)
    dest_dir.mkdir(parents=True, exist_ok=True)

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for key in ("r1", "r2"):
                src = row.get(key, "").strip()
                if not src:
                    continue

                src_path = Path(src)
                if not src_path.exists():
                    raise FileNotFoundError(f"Manifest entry does not exist: {src_path}")

                dest_path = dest_dir / src_path.name
                if dest_path.exists() or dest_path.is_symlink():
                    dest_path.unlink()

                if args.mode == "copy":
                    shutil.copy2(src_path, dest_path)
                else:
                    os.symlink(src_path.resolve(), dest_path)


if __name__ == "__main__":
    main()
