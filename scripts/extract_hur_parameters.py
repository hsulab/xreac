"""Extract Hur et al. (2021) ReaxFF tables from their supplementary PDF.

Usage: python scripts/extract_hur_parameters.py supplement.pdf output.ff
Requires the Poppler pdftotext executable. Only PDF layout is removed; numerical
entries, ordering, repeated records, and dummy atom type X are preserved.
"""

import argparse
from pathlib import Path
import re
import subprocess

HEADER = (
    "Hur et al., RSC Adv. 11, 29298-29307 (2021), DOI:10.1039/D1RA04397H; "
    "CC-BY-NC-3.0; extracted from ESI pp. S10-S12; see NOTICE"
)


def extract(pdf):
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
    rows = text.splitlines()
    start = next(i for i, row in enumerate(rows) if re.match(r"\s*39\s*!\s*Number of general parameters", row))
    rows = [row.strip() for row in rows[start:] if row.strip() and not re.fullmatch(r"\s*S\d+\s*", row)]
    end = next(i for i, row in enumerate(rows) if "Nr of hydrogen bonds" in row)
    count = int(rows[end].split()[0])
    return HEADER + "\n" + "\n".join(rows[: end + 1 + count]) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(extract(args.pdf))


if __name__ == "__main__":
    main()
