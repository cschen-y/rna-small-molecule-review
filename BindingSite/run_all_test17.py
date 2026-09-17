import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


REVIEW = Path(__file__).resolve().parent
ENVIRONMENTS = {
    "rlbind": "core",
    "rnet": "core",
    "rbind": "core",
    "rsite": "core",
    "rsite2": "core",
    "rnasite": "core",
    "rnabind": "MVRBind",
    "multimodrlbp": "Mul",
    "mvrbind": "MVRBind",
}
FOLDERS = {
    "rlbind": "RLBind",
    "rnet": "RNet",
    "rbind": "RBind",
    "rsite": "Rsite",
    "rsite2": "Rsite2",
    "rnasite": "RNAsite",
    "rnabind": "RNABind",
    "multimodrlbp": "MultiModRLBP",
    "mvrbind": "MVRBind",
}


def find_conda(explicit):
    candidates = [
        explicit,
        os.environ.get("CONDA_EXE"),
        shutil.which("conda"),
        str(Path(sys.executable).resolve().parent / "Scripts" / "conda.exe"),
        str(Path(sys.executable).resolve().parent.parent / "Scripts" / "conda.exe"),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    raise FileNotFoundError("Conda was not found. Set CONDA_EXE or pass --conda /path/to/conda.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--methods", nargs="+", choices=[*FOLDERS, "all"], default=["all"])
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument("--conda")
    args = parser.parse_args()
    conda = find_conda(args.conda)
    selected = list(FOLDERS) if "all" in args.methods else args.methods

    failures = []
    for method in selected:
        script = REVIEW / FOLDERS[method] / "run_test17.py"
        command = [
            str(conda),
            "run",
            "--no-capture-output",
            "-n",
            ENVIRONMENTS[method],
            "python",
            str(script),
        ]
        print(f"\n========== {FOLDERS[method]} / Test17 [{ENVIRONMENTS[method]}] ==========", flush=True)
        completed = subprocess.run(command, cwd=str(script.parent))
        if completed.returncode:
            failures.append((method, completed.returncode))
            if not args.continue_on_error:
                raise SystemExit(completed.returncode)

    if failures:
        print("Failed methods: " + ", ".join(f"{name}({code})" for name, code in failures))
        raise SystemExit(1)
    print("\nAll selected Test17 runs completed.")


if __name__ == "__main__":
    main()
