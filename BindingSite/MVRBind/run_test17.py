"""Evaluate the five MVRBind checkpoints on Test17 (Test18 minus 6EZ0A)."""

import json
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    f1_score, matthews_corrcoef, precision_score, recall_score, roc_auc_score,
)
from torch_geometric.loader import DataLoader

from data_process.test_18_dataset import Test18Dataset
from model import MVRBind


HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent
DATA = REVIEW / "data" / "MVRBind"
SEEDS = [134, 221, 333, 432, 511]


def evaluate(seed, loader):
    model = MVRBind(136)
    checkpoint = DATA / "model_parameters" / f"model_seed_{seed}.pt"
    model.load_state_dict(torch.load(checkpoint, map_location="cpu"))
    model.eval()
    labels, scores = [], []
    with torch.no_grad():
        for batch in loader:
            scores.append(model(batch).cpu().numpy())
            labels.append(batch.y.cpu().numpy())
    y_true = np.concatenate(labels).astype(int)
    y_score = np.concatenate(scores)
    y_pred = (y_score >= 0.5).astype(int)
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred, zero_division=0),
        "Recall": recall_score(y_true, y_pred, zero_division=0),
        "F1": f1_score(y_true, y_pred, zero_division=0),
        "MCC": matthews_corrcoef(y_true, y_pred),
        "AUC": roc_auc_score(y_true, y_score),
        "AUPR": average_precision_score(y_true, y_score),
        "BACC": balanced_accuracy_score(y_true, y_pred),
    }


def main():
    test18 = Test18Dataset(str(DATA / "pt"))
    test17 = [test18[index] for index in range(17)]
    nucleotide_count = sum(int(graph.y.numel()) for graph in test17)
    if len(test18) != 18 or len(test17) != 17 or nucleotide_count != 583:
        raise ValueError(
            f"Expected Test18=18 and Test17=17/583 nt, got {len(test18)} and {len(test17)}/{nucleotide_count}"
        )
    loader = DataLoader(test17, batch_size=17, shuffle=False)
    all_metrics = {key: [] for key in ("Accuracy", "Precision", "Recall", "F1", "MCC", "AUC", "AUPR", "BACC")}
    for seed in SEEDS:
        metrics = evaluate(seed, loader)
        for key, value in metrics.items():
            all_metrics[key].append(float(value))
        print(f"seed={seed} " + " ".join(f"{key}={value:.4f}" for key, value in metrics.items()))

    payload = {
        "method": "MVRBind", "dataset": "Test17", "rna_count": 17,
        "nucleotide_count": 583, "seeds": SEEDS, "per_seed": all_metrics,
        "summary": {
            key: {"mean": float(np.mean(values)), "std": float(np.std(values))}
            for key, values in all_metrics.items()
        },
    }
    output = REVIEW / "results" / "MVRBind" / "test17_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
