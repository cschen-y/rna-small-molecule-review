"""Run RLBind on the canonical Test17 benchmark."""

import os
import runpy
from pathlib import Path


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
runpy.run_path(str(HERE / "train_test17.py"), run_name="__main__")
