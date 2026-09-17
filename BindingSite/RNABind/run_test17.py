"""Train the RNABind EGNN (one-hot variant) on Train60 and test on Test17."""

import argparse
import json
import pickle
import random
import sys
from pathlib import Path

import numpy as np
import torch
from Bio.PDB import PDBParser
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    f1_score, matthews_corrcoef, precision_score, recall_score, roc_auc_score,
)
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent
COMMON = REVIEW / "data" / "common"
sys.path.insert(0, str(REVIEW))

from RNABind.models.model import BindingSiteModel


RNA_NAMES = {"A", "U", "C", "G", "DA", "DU", "DC", "DG", "PSU", "CBV", "5BU", "UMS", "CSL", "CCC", "GTP", "GDP", "A23", "U37", "IU"}
BASE_NAMES = {"DA": "A", "DU": "U", "DC": "C", "DG": "G"}


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def residue_coordinate(residue):
    if "C3'" in residue:
        return residue["C3'"].coord
    atoms = [atom.coord for atom in residue if atom.name[0] in {"C", "N", "O", "P"}]
    if not atoms:
        raise ValueError(f"No usable atom in residue {residue}")
    return np.mean(atoms, axis=0)


def make_graph(pdb_path, chain_id, labels, pdb_id):
    structure = PDBParser(QUIET=True).get_structure(pdb_id, str(pdb_path))
    chain = structure[0][chain_id]
    residues = [r for r in chain if r.resname.replace(" ", "") in RNA_NAMES]
    if len(residues) != len(labels):
        raise ValueError(f"{pdb_id}: PDB residues={len(residues)}, labels={len(labels)}")

    coords = torch.tensor(np.asarray([residue_coordinate(r) for r in residues]), dtype=torch.float32)
    x = torch.zeros((len(residues), 4), dtype=torch.float32)
    alphabet = {"A": 0, "U": 1, "C": 2, "G": 3}
    for index, residue in enumerate(residues):
        name = BASE_NAMES.get(residue.resname.strip(), residue.resname.strip())
        x[index, alphabet.get(name, 0)] = 1.0

    nodes = torch.arange(len(residues), dtype=torch.long)
    row = nodes.repeat_interleave(len(residues))
    col = nodes.repeat(len(residues))
    keep = row != col
    edge_index = torch.stack([row[keep], col[keep]], dim=0)
    distances = torch.linalg.vector_norm(coords[edge_index[0]] - coords[edge_index[1]], dim=1)
    centers = torch.linspace(0.0, 20.0, 16)
    edge_attr = torch.exp(-((distances[:, None] - centers[None, :]) / 1.5) ** 2)

    return Data(
        x=x, coord=coords, edge_index=edge_index, edge_attr=edge_attr,
        binding_site=torch.tensor(labels, dtype=torch.float32), pdb_id=pdb_id,
    )


def load_split(fasta_dir, label_path):
    fasta_files = sorted(Path(fasta_dir).glob("*.fasta"))
    with open(label_path, "rb") as file:
        labels = pickle.load(file)
    if len(fasta_files) != len(labels):
        raise ValueError(f"FASTA/label count mismatch: {len(fasta_files)} != {len(labels)}")
    graphs = []
    for fasta, one_labels in zip(fasta_files, labels):
        pdb_id = fasta.stem.upper()
        graphs.append(make_graph(COMMON / "pdbFiles" / f"{pdb_id[:4]}.pdb", pdb_id[4], one_labels, pdb_id))
    return graphs


def evaluate(model, loader, device):
    model.eval()
    labels, scores = [], []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            scores.append(model(batch).view(-1).cpu().numpy())
            labels.append(batch.binding_site.view(-1).cpu().numpy())
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


def run(args):
    train_graphs = load_split(COMMON / "train60_fastas", COMMON / "all_label" / "label" / "train_label.pkl")
    test_graphs = load_split(COMMON / "test17_fastas", COMMON / "all_label" / "label" / "test17_labels.pkl")
    if len(test_graphs) != 17 or sum(graph.num_nodes for graph in test_graphs) != 583:
        raise ValueError("RNABind Test17 must contain 17 RNAs and 583 nucleotides")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    all_metrics = {key: [] for key in ("Accuracy", "Precision", "Recall", "F1", "MCC", "AUC", "AUPR", "BACC")}
    for seed in args.seeds:
        set_seed(seed)
        order = np.random.RandomState(42).permutation(len(train_graphs))
        split = int(0.9 * len(order))
        train_set = [train_graphs[i] for i in order[:split]]
        validation_set = [train_graphs[i] for i in order[split:]]
        train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True)
        validation_loader = DataLoader(validation_set, batch_size=args.batch_size)
        test_loader = DataLoader(test_graphs, batch_size=args.batch_size)

        model = BindingSiteModel(embedding_type="onehot", in_node_nf=128, in_edge_nf=16).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
        criterion = torch.nn.BCELoss()
        best_mcc, best_state = -2.0, None
        for epoch in range(args.epochs):
            model.train()
            for batch in train_loader:
                batch = batch.to(device)
                optimizer.zero_grad()
                loss = criterion(model(batch), batch.binding_site.view(-1, 1))
                loss.backward()
                optimizer.step()
            validation_metrics = evaluate(model, validation_loader, device)
            if validation_metrics["MCC"] > best_mcc:
                best_mcc = validation_metrics["MCC"]
                best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            print(f"seed={seed} epoch={epoch + 1}/{args.epochs} val_MCC={validation_metrics['MCC']:.4f}")

        model.load_state_dict(best_state)
        model.to(device)
        metrics = evaluate(model, test_loader, device)
        for key, value in metrics.items():
            all_metrics[key].append(float(value))
        print(f"seed={seed} Test17 " + " ".join(f"{key}={value:.4f}" for key, value in metrics.items()))

    payload = {
        "method": "RNABind-onehot", "embedding": "onehot", "dataset": "Test17", "rna_count": 17,
        "epochs": args.epochs,
        "nucleotide_count": 583, "seeds": args.seeds,
        "per_seed": all_metrics,
        "summary": {key: {"mean": float(np.mean(values)), "std": float(np.std(values))} for key, values in all_metrics.items()},
    }
    output = REVIEW / "results" / "RNABind" / "test17_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Saved {output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--seeds", type=int, nargs="+", default=[8124, 27045, 58392, 17765, 44322])
    run(parser.parse_args())


if __name__ == "__main__":
    main()
