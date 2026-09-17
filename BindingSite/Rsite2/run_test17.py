"""Run Rsite2 on the canonical Test17 benchmark."""

import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent
os.chdir(HERE)

from method import get_performence


if __name__ == "__main__":
    get_performence(
        str(REVIEW / "data" / "Rsite2" / "SS_NDS"),
        str(REVIEW / "data" / "common" / "test17.fasta"),
        str(REVIEW / "results" / "Rsite2" / "test17"),
    )
