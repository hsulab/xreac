"""Extract Gai2016 Pt/O/H parameters from jp6b01064_si_001.pdf, pp. S7–S10."""

import argparse
from pathlib import Path
import re
import subprocess


def extract(pdf):
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
    rows = text.splitlines()
    start = next(i for i, row in enumerate(rows) if row.startswith("Reactive MD-force field:"))
    rows = [r.strip() for r in rows[start:] if r.strip() and not re.fullmatch(r"\s*S\d+\s*", r)]
    # PDF extraction joins two adjacent fixed-width fields in the dummy X row.
    # Restore their separator; do not change either value.
    rows = [r.replace("5.00009999.9999", "5.0000 9999.9999") for r in rows]
    end = next(i for i, r in enumerate(rows) if "Nr of hydrogen bonds" in r)
    return "\n".join(rows[: end + 1 + int(rows[end].split()[0])]) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(extract(args.pdf))
