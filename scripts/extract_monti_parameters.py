"""Extract the complete Monti2012 force field from jp2121593_si_001.pdf."""

import argparse
from pathlib import Path
import re
import subprocess


def extract(pdf):
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
    text = text[text.index("Reactive MD-force field:") :]
    rows = [row.strip() for row in text.splitlines() if row.strip() and not re.fullmatch(r"\s*S\d+\s*", row)]
    rows[0] += (
        " | Monti et al. 2012; doi:10.1021/jp2121593.s001; CC BY-NC 4.0; PDF whitespace/page-number removal only"
    )
    return "\n".join(rows) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(extract(args.pdf))
