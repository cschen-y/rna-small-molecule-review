import os
import json

import numpy as np
import pickle
from scipy.spatial.distance import pdist, squareform
from scipy.ndimage import gaussian_filter1d
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
    classification_report,
    roc_auc_score,
    roc_curve,
    matthews_corrcoef,
    balanced_accuracy_score,
    average_precision_score
)
from Bio.PDB import PDBParser

def euclidean_distance(coords):
    
    return squareform(pdist(coords))


def compute_distance_curve(coords):
    
    distances = euclidean_distance(coords)
    return np.sum(distances, axis=1)


def smooth_distance_curve(distance_curve, window_size=2):
    
    return gaussian_filter1d(distance_curve, window_size)


def find_local_extrema(smoothed_curve):
    
    local_max = (np.diff(np.sign(np.diff(smoothed_curve))) < 0).nonzero()[0] + 1
    local_min = (np.diff(np.sign(np.diff(smoothed_curve))) > 0).nonzero()[0] + 1
    return local_max, local_min


def check_start_end_functionality(distance_curve, threshold=50):
    
    mean_distance = np.mean(distance_curve)
    start_deviation = abs(distance_curve[0] - mean_distance) / mean_distance * 100
    end_deviation = abs(distance_curve[-1] - mean_distance) / mean_distance * 100

    start_is_functional = start_deviation > threshold
    end_is_functional = end_deviation > threshold

    return start_is_functional, end_is_functional


def identify_functional_sites(distance_curve, smoothed_curve, local_max, local_min, threshold=50):
    
    functional_sites = [0] * len(distance_curve)

    
    start_functional, end_functional = check_start_end_functionality(distance_curve, threshold)
    if start_functional:
        functional_sites[0] = 1
    if end_functional:
        functional_sites[-1] = 1

    
    for idx in np.concatenate((local_max, local_min)):
        functional_sites[idx] = 1

    return functional_sites


def get_single_rna_coordinate_con(chain_id, pdb_file_path):
    residue_id = ['A', 'U', 'C', 'G', 'DA', 'DU', 'DC', 'DG', 'PSU', 'CBV', '5BU', 'UMS', 'CSL', 'CCC', 'GTP', 'GDP',
                  'A23', 'U37', 'IU']
    rna_coordinate = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("MY", pdb_file_path)
    all_coordinate = []
    for model in structure:
        rna_coordinate = []
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    count = 0
                    x, y, z = 0, 0, 0
                    if residue.get_resname().replace(" ", "") in residue_id:
                        for atom in residue:
                            x_atom, y_atom, z_atom = atom.get_coord()
                            
                            x = x + x_atom
                            y = y + y_atom
                            z = z + z_atom
                            count += 1
                        rna_coordinate.append([x / count, y / count, z / count])
        all_coordinate.append(rna_coordinate)
    return all_coordinate

def get_single_rna_coordinate_no_con(chain_id, pdb_file_path):
    residue_id = ['A', 'U', 'C', 'G', 'DA', 'DU', 'DC', 'DG', 'PSU', 'CBV', '5BU', 'UMS', 'CSL', 'CCC', 'GTP', 'GDP',
                  'A23', 'U37', 'IU']
    rna_coordinate = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("MY", pdb_file_path)
    model = structure[0]
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                count = 0
                x, y, z = 0, 0, 0
                if residue.get_resname().replace(" ", "") in residue_id:
                    for atom in residue:
                        x_atom, y_atom, z_atom = atom.get_coord()
                        
                        x = x + x_atom
                        y = y + y_atom
                        z = z + z_atom
                        count += 1
                    rna_coordinate.append([x / count, y / count, z / count])
    return rna_coordinate

def rsite_conformation(pdb_file_path, count):
    all_predict_labels = []
    all_coord = get_single_rna_coordinate_con("A",pdb_file_path)
    for coords in all_coord:
        
        distance_curve = compute_distance_curve(coords)
        
        smoothed_curve = smooth_distance_curve(distance_curve, window_size=2)
        
        local_max, local_min = find_local_extrema(smoothed_curve)
        
        functional_sites = identify_functional_sites(distance_curve, smoothed_curve, local_max, local_min)
        print("Functional Sites (0 or 1):", functional_sites)
        all_predict_labels.extend(functional_sites)

    with open("./all_label/label/label_Tapo.pkl", 'rb') as file:
        test_labels = pickle.load(file)

    test_labels  = [test_labels[count] for i in range(len(all_coord))]
    test_labels = np.array(test_labels)
    test_labels = np.hstack(test_labels)

    accuracy = accuracy_score(test_labels, all_predict_labels)
    precision = precision_score(test_labels, all_predict_labels, average="binary")  
    recall = recall_score(test_labels, all_predict_labels, average="binary")
    f1 = f1_score(test_labels, all_predict_labels, average="binary")
    auc = roc_auc_score(test_labels, all_predict_labels)
    mcc = matthews_corrcoef(test_labels, all_predict_labels)
    print(f"准确率 (Accuracy): {accuracy:.3f}")
    print(f"精确率 (Precision): {precision:.3f}")
    print(f"召回率 (Recall): {recall:.3f}")
    print(f"F1 分数: {f1:.3f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUC: {auc:.3f}")

    return test_labels, np.array(all_predict_labels)


from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score


def read_labels(label_file_path):
    with open(label_file_path, 'rb') as f:
        labels = pickle.load(f)
    if len(labels) != 17:
        raise ValueError(f"Test17 requires 17 label arrays, got {len(labels)}")
    return labels


def rsite(fasta_file_path, pdb_file_path, label_file_path, output_path="../results/Rsite/test17_results.json"):
    all_predict_labels = []
    fasta_list = sorted(name for name in os.listdir(fasta_file_path) if name.lower().endswith(".fasta"))
    if len(fasta_list) != 17:
        raise ValueError(f"Test17 requires 17 FASTA files, got {len(fasta_list)}")

    for fasta_name in fasta_list:
        all_coord = get_single_rna_coordinate_no_con(fasta_name[4], f"{pdb_file_path}/{fasta_name[:4]}.pdb")

        
        distance_curve = compute_distance_curve(all_coord)

        
        smoothed_curve = smooth_distance_curve(distance_curve, window_size=2)

        
        local_max, local_min = find_local_extrema(smoothed_curve)

        
        functional_sites = identify_functional_sites(distance_curve, smoothed_curve, local_max, local_min)

        print(f"Predicted Functional Sites for {fasta_name}: {functional_sites}")
        all_predict_labels.append(functional_sites)

    
    true_labels = read_labels(label_file_path)

    
    all_predict_labels = np.concatenate(all_predict_labels)
    true_labels = np.concatenate(true_labels)
    if len(true_labels) != 583:
        raise ValueError(f"Test17 must contain 583 nucleotides, got {len(true_labels)}")
    if len(all_predict_labels) != len(true_labels):
        raise ValueError(
            f"Prediction/label length mismatch: {len(all_predict_labels)} != {len(true_labels)}"
        )

    accuracy = accuracy_score(true_labels, all_predict_labels)
    precision = precision_score(true_labels, all_predict_labels)
    recall = recall_score(true_labels, all_predict_labels)
    f1 = f1_score(true_labels, all_predict_labels)
    auc = roc_auc_score(true_labels, all_predict_labels)  
    mcc = matthews_corrcoef(true_labels, all_predict_labels)  
    aupr = average_precision_score(true_labels, all_predict_labels)
    bacc = balanced_accuracy_score(true_labels, all_predict_labels)

    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"AUC: {auc:.4f}")
    print(f"MCC: {mcc:.4f}")
    print(f"AUPR: {aupr:.4f}")
    print(f"BACC: {bacc:.4f}")

    metrics = {
        "Accuracy": accuracy, "Precision": precision, "Recall": recall,
        "F1": f1, "MCC": mcc, "AUC": auc, "AUPR": aupr, "BACC": bacc,
    }
    payload = {
        "method": "Rsite", "dataset": "Test17",
        "rna_count": 17, "nucleotide_count": len(true_labels),
        "metrics": {key: float(value) for key, value in metrics.items()},
    }
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f_out:
        json.dump(payload, f_out, indent=2)
    print(f"Test17 results saved to: {output_path}")
    return metrics


    
    

if __name__ == "__main__":
    rsite("./data/test17_fastas", "./data/pdbFiles", "./all_label/label/test17_labels.pkl")






















