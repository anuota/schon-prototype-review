#!/usr/bin/env python3
# SPDX-License-Identifier: EUPL-1.2
"""
Compare the contents of two Bruker .d folders (recursively).

It will:
- list files only in folder A
- list files only in folder B
- list files common to both (by relative path)
- optionally show file size differences for common files
"""

from pathlib import Path


def collect_files(root: Path) -> dict[str, Path]:
    """
    Recursively collect all files under `root` and return
    a mapping: relative_path_str -> absolute Path.
    """
    files = {}
    for p in root.rglob("*"):
        if p.is_file():
            rel = p.relative_to(root)
            # normalize to POSIX-style string for easier comparison
            files[str(rel)] = p
    return files


def compare_d_folders(a: Path, b: Path, show_sizes: bool = True) -> None:
    if not a.is_dir():
        raise NotADirectoryError(f"A is not a directory: {a}")
    if not b.is_dir():
        raise NotADirectoryError(f"B is not a directory: {b}")

    print(f"Comparing .d folders:\n  A = {a}\n  B = {b}\n")

    files_a = collect_files(a)
    files_b = collect_files(b)

    set_a = set(files_a.keys())
    set_b = set(files_b.keys())

    only_a = sorted(set_a - set_b)
    only_b = sorted(set_b - set_a)
    common = sorted(set_a & set_b)

    print(f"Total files in A: {len(set_a)}")
    print(f"Total files in B: {len(set_b)}")
    print(f"Common files:     {len(common)}")
    print(f"Only in A:        {len(only_a)}")
    print(f"Only in B:        {len(only_b)}")
    print()

    # Files only in A
    if only_a:
        print("=== Files only in A ===")
        for rel in only_a:
            print(f"  {rel}")
        print()

    # Files only in B
    if only_b:
        print("=== Files only in B ===")
        for rel in only_b:
            print(f"  {rel}")
        print()

    # Common files
    if common:
        print("=== Files common to A and B (same relative path) ===")
        if show_sizes:
            print("(size in bytes: A_size -> B_size)")
        for rel in common:
            if show_sizes:
                size_a = files_a[rel].stat().st_size
                size_b = files_b[rel].stat().st_size
                same = "==" if size_a == size_b else "!="
                print(f"  {rel} : {size_a} {same} {size_b}")
            else:
                print(f"  {rel}")


if __name__ == "__main__":
    # Adjust these paths as needed.
    # Example for your specific folders:
    a = Path("/Users/anya/Coding/CoreMS/tests/tests_data/ftms/ESI_NEG_SRFA.d")
    b = Path("/Users/anya/Coding/ML_FT-ICR-MS/Measurement_d-files/calibrated_d-files/ESI_neg_G015548-0_100_200sc_000001.d")

    compare_d_folders(a, b)