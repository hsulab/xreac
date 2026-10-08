"""Extract the complete Kim2013 force field from la4006983_si_002.pdf."""

import argparse
from pathlib import Path
import subprocess


def extract(pdf):
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
    rows = [row.strip() for row in text.splitlines() if row.strip()]
    if not rows[0].startswith("Reactive MD-force field:"):
        raise ValueError("Expected the second Kim2013 supplement")
    rows[0] += " | Kim et al. 2013, doi:10.1021/la4006983.s002; CC BY-NC 4.0; PDF whitespace extraction only"
    return "\n".join(rows) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(extract(args.pdf))
