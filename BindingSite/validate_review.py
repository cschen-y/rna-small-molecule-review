"""Validate the standardized Test17 workspace without training models."""

import ast
import json
import pickle
from pathlib import Path


ROOT = Path(__file__).resolve().parent
METHODS = ["RLBind", "RNet", "RBind", "Rsite", "Rsite2", "RNAsite", "RNABind", "MultiModRLBP", "MVRBind"]


def main():
    for method in METHODS:
        script = ROOT / method / "run_test17.py"
        if not script.exists():
            raise FileNotFoundError(script)

    for script in ROOT.rglob("*.py"):
        ast.parse(script.read_text(encoding="utf-8"), filename=str(script))

    fasta_files = list((ROOT / "data" / "common" / "test17_fastas").glob("*.fasta"))
    train_files = list((ROOT / "data" / "common" / "train60_fastas").glob("*.fasta"))
    with (ROOT / "data" / "common" / "all_label" / "label" / "test17_labels.pkl").open("rb") as file:
        labels = pickle.load(file)
    if len(fasta_files) != 17 or len(train_files) != 60 or len(labels) != 17 or sum(map(len, labels)) != 583:
        raise ValueError("Canonical Train60/Test17 data validation failed")

    available_results = []
    for method in METHODS:
        result = ROOT / "results" / method / "test17_results.json"
        if result.exists():
            json.loads(result.read_text(encoding="utf-8"))
            available_results.append(method)
    print("Layout valid: 9 method runners, Train60=60, Test17=17/583 nt")
    print("JSON results: " + ", ".join(available_results))


if __name__ == "__main__":
    main()
