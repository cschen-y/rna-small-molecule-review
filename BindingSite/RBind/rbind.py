import numpy as np
import os
import json
import pickle
from get_connection_ad import calculate_contact_matrix
from get_topology import get_node_topology
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, matthews_corrcoef, average_precision_score,
    balanced_accuracy_score, confusion_matrix
)

def predict(closeness, degree):
    
    closeness_avg = np.mean(closeness)
    degree_avg = np.mean(degree)
    closeness_std = np.std(closeness)
    degree_std = np.std(degree)

    
    closeness_cutoff = closeness_avg + closeness_std
    degree_cutoff = degree_avg + degree_std

    
    result = [
        1 if c > closeness_cutoff and d > degree_cutoff else 0
        for c, d in zip(closeness, degree)
    ]
    return result


import os
import numpy as np
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, matthews_corrcoef

def rbind_huang_test9(fasta_label_file_path, pdb_file_path):
    all_prediction = []
    all_labels = []

    with open(fasta_label_file_path, 'r') as f:
        lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            pdb_id = lines[i][1:].strip()  
            sequence = lines[i+1].strip()
            labels_line = lines[i+2].strip()
            labels = list(map(int, labels_line.split(',')))

            
            pdb_filename = f'{pdb_file_path}/{pdb_id[:4].lower()}.pdb'
            chain_id = pdb_id[4]  

            adj_matrix = calculate_contact_matrix(pdb_filename, chain_id=chain_id, cutoff=8.0, mode="rep", output_dir="data/huang_new_ad")
            topology = get_node_topology(adj_matrix)
            topology = np.array(topology)
            degree = topology[:, 0]
            closeness = topology[:, 1]

            rna_predict = predict(closeness, degree)

            print(f'Processing: {pdb_id}')
            print(f'Prediction Length: {len(rna_predict)}, Label Length: {len(labels)}')

            all_prediction.extend(rna_predict)
            all_labels.extend(labels)

            i += 3
        else:
            i += 1

    
    accuracy = accuracy_score(all_labels, all_prediction)
    precision = precision_score(all_labels, all_prediction)
    recall = recall_score(all_labels, all_prediction)
    f1 = f1_score(all_labels, all_prediction)
    auc = roc_auc_score(all_labels, all_prediction)
    mcc = matthews_corrcoef(all_labels, all_prediction)
    aupr = average_precision_score(all_labels, all_prediction)
    bacc = balanced_accuracy_score(all_labels, all_prediction)

    
    print(f"准确率 (Accuracy): {accuracy:.3f}")
    print(f"精确率 (Precision): {precision:.3f}")
    print(f"召回率 (Recall): {recall:.3f}")
    print(f"F1 分数: {f1:.3f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUC: {auc:.3f}")
    print(f"AUPR: {aupr:.3f}")
    print(f"BACC: {bacc:.3f}")

def rbind_all_pdb(fasta_label_file_path, pdb_file_path):
    all_prediction = []
    all_labels = []

    with open(fasta_label_file_path, 'r') as f:
        lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            pdb_id = lines[i][1:].strip()  
            sequence = lines[i+1].strip()
            labels_line = lines[i+2].strip()
            labels = list(map(int, labels_line.split(',')))

            pdb_filename = f'{pdb_file_path}/{pdb_id[:5].upper()}.pdb'
            
            chain_id = pdb_id[4]  

            adj_matrix = calculate_contact_matrix(pdb_filename, chain_id=chain_id, cutoff=8.0, mode="rep", output_dir="data/huang_new_ad")
            topology = get_node_topology(adj_matrix)
            topology = np.array(topology)
            degree = topology[:, 0]
            closeness = topology[:, 1]

            rna_predict = predict(closeness, degree)

            print(f'Processing: {pdb_id}')
            print(f'Prediction Length: {len(rna_predict)}, Label Length: {len(labels)}')

            all_prediction.extend(rna_predict)
            all_labels.extend(labels)

            i += 3
        else:
            i += 1

    
    accuracy = accuracy_score(all_labels, all_prediction)
    precision = precision_score(all_labels, all_prediction)
    recall = recall_score(all_labels, all_prediction)
    f1 = f1_score(all_labels, all_prediction)
    auc = roc_auc_score(all_labels, all_prediction)
    mcc = matthews_corrcoef(all_labels, all_prediction)
    aupr = average_precision_score(all_labels, all_prediction)
    bacc = balanced_accuracy_score(all_labels, all_prediction)

    
    print(f"准确率 (Accuracy): {accuracy:.3f}")
    print(f"精确率 (Precision): {precision:.3f}")
    print(f"召回率 (Recall): {recall:.3f}")
    print(f"F1 分数: {f1:.3f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUC: {auc:.3f}")
    print(f"AUPR: {aupr:.3f}")
    print(f"BACC: {bacc:.3f}")

def read_labels(label_file_path):
    labels_result = []
    with open(label_file_path, 'rb') as f:
        labels = pickle.load(f)
    for i in labels:
        labels_result.extend(i)
    return labels_result

def rbind(fasta_file_path, pdb_file_path, ground_true_path, output_path="../results/RBind/test17_results.json"):
    all_prediction = []
    all_labels = []
    ground_true_all_labels = []

    with open(ground_true_path, 'rb') as f:
        ground_true_all_labels = pickle.load(f)

    fasta_files = sorted(name for name in os.listdir(fasta_file_path) if name.lower().endswith(".fasta"))
    if len(fasta_files) != 17 or len(ground_true_all_labels) != 17:
        raise ValueError(
            f"Test17 requires 17 FASTA files and 17 label arrays; got "
            f"{len(fasta_files)} and {len(ground_true_all_labels)}"
        )

    for pdb_id, one_true_label in zip(fasta_files, ground_true_all_labels):
        pdb_filename = f'{pdb_file_path}/{pdb_id[:4].upper()}.pdb'
        
        chain_id = pdb_id[4]  
        adj_matrix = calculate_contact_matrix(pdb_filename, chain_id=chain_id, cutoff=8.0, mode="rep", output_dir="data/test17_ad")
        topology = get_node_topology(adj_matrix)
        topology = np.array(topology)
        degree = topology[:, 0]
        closeness = topology[:, 1]

        rna_predict = predict(closeness, degree)

        print(f'Processing: {pdb_id}')
        print(f'Prediction Length: {len(rna_predict)}')
        if len(rna_predict) != len(one_true_label):
            raise ValueError(
                f"Length mismatch for {pdb_id}: prediction={len(rna_predict)}, "
                f"label={len(one_true_label)}"
            )

        all_prediction.extend(rna_predict)
        all_labels.extend(one_true_label)
    
    accuracy = accuracy_score(all_labels, all_prediction)
    precision = precision_score(all_labels, all_prediction)
    recall = recall_score(all_labels, all_prediction)
    f1 = f1_score(all_labels, all_prediction)
    auc = roc_auc_score(all_labels, all_prediction)
    mcc = matthews_corrcoef(all_labels, all_prediction)
    aupr = average_precision_score(all_labels, all_prediction)
    bacc = balanced_accuracy_score(all_labels, all_prediction)
    if len(all_labels) != 583:
        raise ValueError(f"Test17 must contain 583 nucleotides, got {len(all_labels)}")

    metrics = {
        "Accuracy": accuracy, "Precision": precision, "Recall": recall,
        "F1": f1, "MCC": mcc, "AUC": auc, "AUPR": aupr, "BACC": bacc,
    }

    
    print(f"准确率 (Accuracy): {accuracy:.4f}")
    print(f"精确率 (Precision): {precision:.4f}")
    print(f"召回率 (Recall): {recall:.4f}")
    print(f"F1 分数: {f1:.4f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUC: {auc:.4f}")
    print(f"AUPR: {aupr:.4f}")
    print(f"BACC: {bacc:.4f}")

    payload = {
        "method": "RBind", "dataset": "Test17",
        "rna_count": 17, "nucleotide_count": len(all_labels),
        "metrics": {key: float(value) for key, value in metrics.items()},
    }
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f_out:
        json.dump(payload, f_out, indent=2)
    print(f"Test17 results saved to: {output_path}")
    return metrics


def rbind_conformation():
    
    all_auc_data1 = []
    all_auc_data2 = []
    all_auc_data3 = []
    all_auc_combined = []

    
    with open("rnasite/label/rnasite_conformation.pkl", "rb") as f:
        data1, data2, data3 = pickle.load(f)

    
    combined_ground_truth = []
    combined_predictions = []

    
    print("Evaluating data1:")
    for features, label in zip(data1[0], data1[1]):
        features = features[:, -10:-8]
        rna_predict = predict(features[:, 0], features[:, 1])

        
        accuracy = accuracy_score(label, rna_predict)
        precision = precision_score(label, rna_predict, average="binary")  
        recall = recall_score(label, rna_predict, average="binary")
        f1 = f1_score(label, rna_predict, average="binary")
        auc = roc_auc_score(label, rna_predict)
        mcc = matthews_corrcoef(label, rna_predict)

        
        print(
            f"Accuracy: {accuracy:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, F1: {f1:.3f}, MCC: {mcc:.3f}, AUC: {auc:.3f}")

        all_auc_data1.append(auc)
        combined_ground_truth.extend(label)
        combined_predictions.extend(rna_predict)

    
    print("\nEvaluating data2:")
    for features, label in zip(data2[0], data2[1]):
        features = features[:, -10:-8]
        rna_predict = predict(features[:, 0], features[:, 1])

        
        accuracy = accuracy_score(label, rna_predict)
        precision = precision_score(label, rna_predict, average="binary")
        recall = recall_score(label, rna_predict, average="binary")
        f1 = f1_score(label, rna_predict, average="binary")
        auc = roc_auc_score(label, rna_predict)
        mcc = matthews_corrcoef(label, rna_predict)

        
        print(
            f"Accuracy: {accuracy:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, F1: {f1:.3f}, MCC: {mcc:.3f}, AUC: {auc:.3f}")

        all_auc_data2.append(auc)
        combined_ground_truth.extend(label)
        combined_predictions.extend(rna_predict)

    
    print("\nEvaluating data3:")
    for features, label in zip(data3[0], data3[1]):
        features = features[:, -10:-8]
        rna_predict = predict(features[:, 0], features[:, 1])

        
        accuracy = accuracy_score(label, rna_predict)
        precision = precision_score(label, rna_predict, average="binary")
        recall = recall_score(label, rna_predict, average="binary")
        f1 = f1_score(label, rna_predict, average="binary")
        auc = roc_auc_score(label, rna_predict)
        mcc = matthews_corrcoef(label, rna_predict)

        
        print(
            f"Accuracy: {accuracy:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, F1: {f1:.3f}, MCC: {mcc:.3f}, AUC: {auc:.3f}")

        all_auc_data3.append(auc)
        combined_ground_truth.extend(label)
        combined_predictions.extend(rna_predict)

    
    print("\nEvaluating combined data (all datasets):")
    accuracy_combined = accuracy_score(combined_ground_truth, combined_predictions)
    precision_combined = precision_score(combined_ground_truth, combined_predictions, average="binary")
    recall_combined = recall_score(combined_ground_truth, combined_predictions, average="binary")
    f1_combined = f1_score(combined_ground_truth, combined_predictions, average="binary")
    auc_combined = roc_auc_score(combined_ground_truth, combined_predictions)
    mcc_combined = matthews_corrcoef(combined_ground_truth, combined_predictions)

    
    print(
        f"Combined Accuracy: {accuracy_combined:.3f}, Combined Precision: {precision_combined:.3f}, Combined Recall: {recall_combined:.3f}")
    print(f"Combined F1: {f1_combined:.3f}, Combined MCC: {mcc_combined:.3f}, Combined AUC: {auc_combined:.3f}")

    
    all_auc_combined.extend(all_auc_data1 + all_auc_data2 + all_auc_data3)

    
    print(f"AUC for data1: {np.mean(all_auc_data1):.3f}")
    print(f"AUC for data2: {np.mean(all_auc_data2):.3f}")
    print(f"AUC for data3: {np.mean(all_auc_data3):.3f}")
    print(f"Average AUC for combined data: {np.mean(all_auc_combined):.3f}")

if __name__ == '__main__':
    
    pdb_file_path = "data/pdbFiles"

    
    
    

    

    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    

    rbind("./data/test17_fastas","./data/pdbFiles", "./all_label/label/test17_labels.pkl")
    
    





