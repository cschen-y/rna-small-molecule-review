import os
import re
import json
import pickle
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, matthews_corrcoef, average_precision_score,
    balanced_accuracy_score
)

def process_resullt(file_path):
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    predicted_sites = []
    
    
    matches = re.findall(r'Site#\d+:\s*(\d+(?:,\d+)*)', content)
    
    for match in matches:
        sites = match.split(',')
        for site in sites:
            predicted_sites.append(int(site))
    
    
    return predicted_sites

def get_label(file_path):
    with open(file_path, 'rb') as file:  
        test_labels = pickle.load(file)
        temp = []
        for i in test_labels:
            temp.extend(i)
        test_labels = temp
    return test_labels




































def get_performence(result_file_path, fasta_file, output_prefix="../../results/rsite2_test17"):
    all_predict_labels = []
    all_ground_truth = []

    all_predict_rna_list = []
    all_ground_truth_rna_list = []

    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            fasta_name = lines[i][1:].strip()  
            seq_line = lines[i + 1].strip()
            label_line = lines[i + 2].strip()

            rna_result_path = os.path.join(result_file_path, f"{fasta_name[:5]}_ss_nds_gd_extrema_2-2.txt")
            if not os.path.exists(rna_result_path):
                raise FileNotFoundError(f"Missing Rsite2 output for {fasta_name}: {rna_result_path}")
            rna_predict_label_index = process_resullt(rna_result_path)

            rna_length = len(seq_line)
            rna_predict_label = []
            for j in range(rna_length):
                if j + 1 in rna_predict_label_index:
                    rna_predict_label.append(1)
                else:
                    rna_predict_label.append(0)

            ground_truth = [int(x.strip()) for x in label_line.split(',')]
            assert len(rna_predict_label) == len(ground_truth), f"Length mismatch for {fasta_name}"

            all_predict_labels.extend(rna_predict_label)
            all_ground_truth.extend(ground_truth)
            all_predict_rna_list.append(rna_predict_label)
            all_ground_truth_rna_list.append(ground_truth)

            i += 3
        else:
            i += 1

    if len(all_ground_truth_rna_list) != 17:
        raise ValueError(f"Test17 requires 17 RNA records, got {len(all_ground_truth_rna_list)}")
    if len(all_ground_truth) != 583:
        raise ValueError(f"Test17 must contain 583 nucleotides, got {len(all_ground_truth)}")

    
    accuracy = accuracy_score(all_ground_truth, all_predict_labels)
    precision = precision_score(all_ground_truth, all_predict_labels, average="binary")
    recall = recall_score(all_ground_truth, all_predict_labels, average="binary")
    f1 = f1_score(all_ground_truth, all_predict_labels, average="binary")
    auc = roc_auc_score(all_ground_truth, all_predict_labels)
    mcc = matthews_corrcoef(all_ground_truth, all_predict_labels)
    aupr = average_precision_score(all_ground_truth, all_predict_labels)
    bacc = balanced_accuracy_score(all_ground_truth, all_predict_labels)

    
    print(f"准确率 (Accuracy): {accuracy:.3f}")
    print(f"精确率 (Precision): {precision:.3f}")
    print(f"召回率 (Recall): {recall:.3f}")
    print(f"F1 分数: {f1:.3f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUC: {auc:.3f}")
    print(f"AUPR: {aupr:.3f}")
    print(f"BACC: {bacc:.3f}")

    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_prefix}_ground_truth.pkl", 'wb') as file:
        pickle.dump(all_ground_truth_rna_list, file)
    with open(f"{output_prefix}_predict.pkl", 'wb') as file:
        pickle.dump(all_predict_rna_list, file)

    metrics = {
        "Accuracy": accuracy, "Precision": precision, "Recall": recall,
        "F1": f1, "MCC": mcc, "AUC": auc, "AUPR": aupr, "BACC": bacc,
    }
    payload = {
        "method": "Rsite2", "dataset": "Test17",
        "rna_count": 17, "nucleotide_count": len(all_ground_truth),
        "metrics": {key: float(value) for key, value in metrics.items()},
    }
    with open(f"{output_prefix}_results.json", "w", encoding="utf-8") as file:
        json.dump(payload, file, indent=2)
    print(f"Test17 results saved to: {output_prefix}_results.json")
    return metrics


if __name__ == "__main__":
    get_performence("PS/SS_NDS", "../rnasite/test17.fasta")
