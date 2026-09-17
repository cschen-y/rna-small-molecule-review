"""Run RNet on the canonical Test17 benchmark."""

import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
os.chdir(HERE)

from method import run_with_seeds


if __name__ == "__main__":
    run_with_seeds()
