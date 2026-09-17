"""Collect all available Test17 JSON results into one CSV table."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
METHODS = ["RLBind", "RNet", "RBind", "Rsite", "Rsite2", "RNAsite", "RNABind", "MultiModRLBP", "MVRBind"]
METRICS = ["Accuracy", "Precision", "Recall", "F1", "MCC", "AUC", "AUPR", "BACC"]
PROTOCOLS = {
    "RLBind": "1-epoch smoke run",
    "RNet": "full configured run (5 seeds)",
    "RBind": "deterministic full evaluation",
    "Rsite": "deterministic full evaluation",
    "Rsite2": "deterministic full evaluation",
    "RNAsite": "full configured run (5 seeds)",
    "RNABind": "1-epoch smoke run (one-hot variant)",
    "MultiModRLBP": "1-epoch smoke run (5 seeds)",
    "MVRBind": "5 pretrained checkpoints",
}


def main():
    rows = []
    for method in METHODS:
        path = ROOT / "results" / method / "test17_results.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        source = payload.get("summary", payload.get("metrics", {}))
        row = {
            "method": payload.get("method", method),
            "protocol": PROTOCOLS[method],
            "epochs": payload.get("epochs", ""),
        }
        for metric in METRICS:
            value = source.get(metric, "")
            row[metric] = value.get("mean", "") if isinstance(value, dict) else value
        rows.append(row)

    output = ROOT / "results" / "test17_summary.csv"
    with output.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=["method", "protocol", "epochs", *METRICS])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Collected {len(rows)} methods into {output}")


if __name__ == "__main__":
    main()
