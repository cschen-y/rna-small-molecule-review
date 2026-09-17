"""Run Rsite on the canonical Test17 benchmark."""

import os
import runpy
from pathlib import Path


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
runpy.run_path(str(HERE / "rsite.py"), run_name="__main__")
