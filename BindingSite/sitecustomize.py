"""Expose project-vendored dependencies to the default Python interpreter."""

import sys
from pathlib import Path


VENDOR = Path(__file__).resolve().parent / "vendor"
if VENDOR.is_dir() and str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
