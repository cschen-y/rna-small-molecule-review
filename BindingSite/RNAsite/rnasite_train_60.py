import numpy as np
import os
import json
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, matthews_corrcoef
)
from get_asa_rnasite import get_new_asa
from get_ln_rnasite import get_ln
from get_connection_ad import calculate_contact_matrix
from get_topology import get_node_topology
from sklearn.metrics import average_precision_score
from sklearn.metrics import balanced_accuracy_score

def read_clustal(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.strip() == "" or line.startswith("CLUSTAL") or line.startswith(" "):
            continue
        parts = line.split()
        if len(parts) == 2:
            seq_id, seq = parts
            sequences.setdefault(seq_id, "")
            sequences[seq_id] += seq
    return list(sequences.values()), sequences

def henikoff_weights(msa):
    n_seq = len(msa)
    n_pos = len(msa[0])
    weights = np.zeros(n_seq)
    for j in range(n_pos):
        column = [seq[j].replace("T", "U") for seq in msa]
        freqs = Counter(column)
        num_types = len(freqs)
        for i, nucleotide in enumerate(column):
            weights[i] += (1 / num_types) * freqs[nucleotide]
    return weights / np.sum(weights)

def conservation_scores(msa, weights):
    n_pos = len(msa[0])
    scores = []
    for j in range(n_pos):
        column = [seq[j].replace("T", "U") for seq in msa]
        weighted_counts = {nuc: 0 for nuc in ['A', 'C', 'G', 'U', '-']}
        for i, nucleotide in enumerate(column):
            if nucleotide in weighted_counts:
                weighted_counts[nucleotide] += weights[i]
        scores.append(weighted_counts)
    return scores

def sliding_window_context(features, window_size):
    num_positions = len(features)
    num_features = len(features[0])
    context_features = []
    for j in range(num_positions):
        feature_vector = []
        for offset in range(-window_size, window_size + 1):
            pos = j + offset
            if 0 <= pos < num_positions:
                feature_vector.extend(features[pos])
            else:
                feature_vector.extend([0] * num_features)
        context_features.append(feature_vector)
    return context_features

def encode_features(msa, scores, window_size, msa_name_seq):
    index_gap = [i for i, ch in enumerate(next(v for k, v in msa_name_seq.items() if not k.startswith("seq"))) if ch == '-']
    raw_features = [[score.get(nuc, 0) for nuc in ['A', 'C', 'G', 'U', '-']] for score in scores]
    context_features = sliding_window_context(raw_features, window_size)
    return np.array([context_features[i] for i in range(len(context_features)) if i not in index_gap], dtype=object)

def process_multiple_clustal_files(directory, window_size):
    all_features, lengths = [], []
    for file_name in sorted(os.listdir(directory)):
        if file_name.endswith(".aln"):
            msa, msa_name_seq = read_clustal(os.path.join(directory, file_name))
            weights = henikoff_weights(msa)
            scores = conservation_scores(msa, weights)
            features = encode_features(msa, scores, window_size, msa_name_seq)
            lengths.append(len(features))
            all_features.append(features)
    return np.vstack(all_features), lengths

def get_msa_huang(directory, pdb_list ,window_size):
    all_features, lengths = [], []
    for file_name in pdb_list:
        file_path = os.path.join(directory, f"{file_name[1:]}_msa.npy")
        data = np.load(file_path)
        lengths.append(len(data))
        features = sliding_window_context(data, window_size)
        all_features.append(features)
    return np.vstack(all_features), lengths

def symmetric_sliding_window(features, half_window_size, padding_mode="zero"):
    N, D = features.shape
    window_size = 2 * half_window_size + 1
    results = []
    for i in range(N):
        start_idx = max(0, i - half_window_size)
        end_idx = min(N, i + half_window_size + 1)
        window_data = features[start_idx:end_idx]
        if padding_mode == "extend":
            left_pad = features[0:1].repeat(max(0, half_window_size - i), axis=0) if start_idx == 0 else np.empty((0, D))
            right_pad = features[-1:].repeat(max(0, i + half_window_size + 1 - N), axis=0) if end_idx == N else np.empty((0, D))
        elif padding_mode == "zero":
            left_pad = np.zeros((max(0, half_window_size - i), D))
            right_pad = np.zeros((max(0, i + half_window_size + 1 - N), D))
        else:
            raise ValueError("padding_mode must be 'extend' or 'zero'")
        padded_window = np.vstack([left_pad, window_data, right_pad])
        results.append(padded_window.flatten())
    return np.array(results)

def get_cl_dg(pdb_path, fasta_path, output_dir):
    result = []
    for pdb in sorted(os.listdir(fasta_path)):
        adj = calculate_contact_matrix(f'{pdb_path}/{pdb[:4]}.pdb', pdb[4], 8.0, "RNABind", output_dir)
        topo = get_node_topology(adj)
        with open(f"{fasta_path}/{pdb}", 'r', encoding='utf-8') as f:
            lines = f.readlines()
        if len(topo) != len(lines[1].strip()):
            print(pdb[:4])
        result.append(np.array(topo))
    return result

def get_cl_dg_huang(pdb_path, pdb_list, output_dir):
    result = []
    for pdb in pdb_list:
        pdb = pdb[1:]
        adj = calculate_contact_matrix(f'{pdb_path}/{pdb}.pdb', pdb[4], 8.0, "huang", output_dir)
        topo = get_node_topology(adj)
        result.append(np.array(topo))
    return result

def merge_feature(cl_dg, asa, ln, msa=None):
    all_features = []
    for i in range(len(cl_dg)):
        if msa is not None:
            merged = np.hstack((cl_dg[i], asa[i], ln[i], msa[i].reshape(-1, 1)))
        else:
            merged = np.hstack((cl_dg[i], asa[i], ln[i]))
            
        all_features.append(symmetric_sliding_window(merged, window_size))
    return np.vstack(all_features), [len(f) for f in cl_dg]

def split_feature(array, split_lengths):
    result, start = [], 0
    for length in split_lengths:
        result.append(array[start:start + length])
        start += length
    return result

def get_asa_own(file_path,pdb_list):
    all_asa = []
    for pdb in pdb_list:
        rna_asa_path = os.path.join(file_path, f"{pdb[1:]}.npy")
        rna_asa = np.load(rna_asa_path)[:,-1]
        all_asa.append(np.array(rna_asa).reshape(-1, 1))
    return all_asa


train_fasta = "../data/RNAsite/data_v2/RB_Train60_label2.fasta"
test18_fasta = "../data/common/test17.fasta"
train_msa_path = "../data/RNAsite/data_v2/selected_RB_msa/"
test_msa_path = "../data/RNAsite/test17_msa"
train_pdb_path = "../data/RNAsite/data_v2/normalized_pdb"
test_pdb_path = "../data/common/pdbFiles"
train_asa_path = "../data/RNAsite/data_v2/All_SASA_new"
test_asa_path = "../data/RNAsite/train60_test18_asa"
test_fasta_dir = "../data/common/test17_fastas"
train_ad_output_path = "../data/RNAsite/data_v2/ad"
test_ad_output_path = "../data/RNAsite/test17_ad"
window_size = 12




def read_fasta_label(fasta_file):
    seq_list = []
    label_list = []

    with open(fasta_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    records = [lines[i:i + 3] for i in range(0, len(lines), 3)]

    for rec in records:
        seq_list.append(rec[0])
        label_list.extend(list(map(int, rec[2].split(","))))

    return seq_list, label_list


train_list, train_labels = read_fasta_label(train_fasta)
test_list, test_labels = read_fasta_label(test18_fasta)

if len(test_list) != 17 or len(test_labels) != 583:
    raise ValueError(
        f"Test17 must contain 17 RNAs and 583 nucleotides; "
        f"got {len(test_list)} and {len(test_labels)}"
    )

train_labels = np.array(train_labels)
test_labels = np.array(test_labels)




train_asa = get_asa_own(train_asa_path, train_list)
test_asa = get_new_asa(test_asa_path, test_fasta_dir)

train_cl_dg = get_cl_dg_huang(train_pdb_path, train_list, train_ad_output_path)
test_cl_dg = get_cl_dg(test_pdb_path, test_fasta_dir, test_ad_output_path)

train_msa, train_lengths = get_msa_huang(train_msa_path, train_list, window_size)
test_msa, test_lengths = process_multiple_clustal_files(test_msa_path, window_size)

train_ln = get_ln(train_pdb_path, train_list)
test_ln = get_ln(test_pdb_path, test_list)

if sum(test_lengths) != 583:
    raise ValueError(f"Test17 MSA features must contain 583 nucleotides, got {sum(test_lengths)}")




num_runs = 5
seeds = [0, 1, 2, 3, 4]

metrics = {
    "Accuracy": [],
    "Precision": [],
    "Recall": [],
    "F1": [],
    "AUC": [],
    "AUPR": [],
    "BACC": [],
    "MCC": []
}

for run, seed in enumerate(seeds):
    print(f"\n===== Run {run + 1} / {num_runs} (seed={seed}) =====")

    
    model_msa = RandomForestClassifier(
        n_estimators=200,
        random_state=seed,
        n_jobs=-1
    )

    model_msa.fit(train_msa, train_labels)
    train_msa_f = model_msa.predict_proba(train_msa)[:, 1]
    test_msa_f = model_msa.predict_proba(test_msa)[:, 1]

    train_msa_p = split_feature(train_msa_f, train_lengths)
    test_msa_p = split_feature(test_msa_f, test_lengths)

    
    train_X_all, _ = merge_feature(train_cl_dg, train_asa, train_ln, train_msa_p)
    test_X_all, _ = merge_feature(test_cl_dg, test_asa, test_ln, test_msa_p)

    model_all = RandomForestClassifier(
        n_estimators=200,
        random_state=seed,
        n_jobs=-1
    )

    model_all.fit(train_X_all, train_labels)
    y_pred_proba = model_all.predict_proba(test_X_all)[:, 1]
    y_pred = (y_pred_proba >= 0.5).astype(int)

    
    metrics["Accuracy"].append(accuracy_score(test_labels, y_pred))
    metrics["Precision"].append(precision_score(test_labels, y_pred))
    metrics["Recall"].append(recall_score(test_labels, y_pred))
    metrics["F1"].append(f1_score(test_labels, y_pred))
    metrics["AUC"].append(roc_auc_score(test_labels, y_pred_proba))
    metrics["AUPR"].append(average_precision_score(test_labels, y_pred_proba))
    metrics["BACC"].append(balanced_accuracy_score(test_labels, y_pred))
    metrics["MCC"].append(matthews_corrcoef(test_labels, y_pred))




print("\n===== Final Results (Mean ± Std over 5 runs) =====")
for k, v in metrics.items():
    mean = np.mean(v)
    std = np.std(v)
    print(f"{k}: {mean:.4f} ± {std:.4f}")

output_path = "../results/RNAsite/test17_results.json"
os.makedirs(os.path.dirname(output_path), exist_ok=True)
payload = {
    "method": "RNAsite",
    "dataset": "Test17",
    "rna_count": 17,
    "nucleotide_count": 583,
    "seeds": seeds,
    "per_seed": {key: [float(value) for value in values] for key, values in metrics.items()},
    "summary": {
        key: {"mean": float(np.mean(values)), "std": float(np.std(values))}
        for key, values in metrics.items()
    },
}
with open(output_path, "w", encoding="utf-8") as f_out:
    json.dump(payload, f_out, indent=2)
print(f"Test17 results saved to: {output_path}")

