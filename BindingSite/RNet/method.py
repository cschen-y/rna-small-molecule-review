import os
import json
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, matthews_corrcoef, roc_auc_score
import os
import re
import pickle
from Bio.PDB import PDBParser
import networkx as nx

def sliding_window_context(features, window_size):
    num_positions = len(features)
    num_features_per_position = len(features[0])
    context_features = []

    for j in range(num_positions):
        feature_vector = []
        for offset in range(-window_size, window_size + 1):
            pos = j + offset
            if 0 <= pos < num_positions:
                feature_vector.extend(features[pos])
            else:
                feature_vector.extend([0] * num_features_per_position)
        context_features.append(feature_vector)
    return context_features

def get_node_topology(adj_matrix):
    """
    计算给定邻接矩阵的图的节点拓扑属性。

    参数：
        adj_matrix (numpy.ndarray): 图的邻接矩阵。

    返回：
        list: 每个节点的属性列表 [度, 接近中心性, 邻居平均度, 中介中心性, 离心率]。
    """
    
    G = nx.from_numpy_matrix(adj_matrix)

    
    DG = dict(G.degree())

    
    NC = {
        node: np.mean([DG[neighbor] for neighbor in G.neighbors(node)]) if DG[node] > 0 else 0
        for node in G.nodes()
    }

    
    BC = nx.betweenness_centrality(G)  
    CL = nx.closeness_centrality(G)  

    
    try:
        EC = nx.eccentricity(G)  
    except nx.NetworkXError:
        
        EC = {}
        for component in nx.connected_components(G):
            subgraph = G.subgraph(component)
            sub_eccentricity = nx.eccentricity(subgraph)
            EC.update(sub_eccentricity)

    
    all_properties = []
    for node in G.nodes():
        node_properties = [DG[node], NC[node], BC[node], CL[node], EC[node]]
        all_properties.append(node_properties)

    return all_properties

def calculate_contact_matrix(input_file, chain_id, cutoff, mode, output_dir):
    
    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    output_filename = f"{base_filename}_{chain_id}_{mode}.npy"
    output_path = os.path.join(output_dir, output_filename)

    
    if os.path.exists(output_path):
        return np.load(output_path)

    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', input_file)
    residue_id = ['A', 'U', 'C', 'G', 'DA', 'DU', 'DC', 'DG', 'PSU', 'CBV', '5BU', 'UMS', 'CSL', 'CCC', 'GTP', 'GDP',
                  'A23', 'U37', 'IU']

    atoms = []
    model = structure[0]
    count = 0
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                if residue.get_resname().replace(" ", "") in residue_id:
                    count += 1
                    for atom in residue:
                        atoms.append((atom.get_serial_number(),
                                      residue.get_id()[1],
                                      atom.get_coord(),
                                      chain.get_id()))

    number_of_atoms = len(atoms)

    atom_serial_numbers = np.array([atom[0] for atom in atoms])
    amino_acid_numbers = np.array([atom[1] for atom in atoms])
    atom_positions = np.array([atom[2] for atom in atoms])
    chain_ids = [atom[3] for atom in atoms]

    number_of_amino_acids = 1
    revised_amino_acid_numbers = np.ones(number_of_atoms, dtype=int)
    for i in range(1, number_of_atoms):
        if abs(amino_acid_numbers[i] - amino_acid_numbers[i - 1]) > 0:
            number_of_amino_acids += 1
            revised_amino_acid_numbers[i:] = number_of_amino_acids

    contact_matrix = np.zeros((count, count), dtype=int)

    for i in range(number_of_atoms):
        for j in range(number_of_atoms):
            if abs(revised_amino_acid_numbers[i] - revised_amino_acid_numbers[j]) <= 1:
                continue
            distance = np.linalg.norm(atom_positions[i] - atom_positions[j])
            if distance <= cutoff:
                contact_matrix[revised_amino_acid_numbers[i] - 1, revised_amino_acid_numbers[j] - 1] = 1

    
    os.makedirs(output_dir, exist_ok=True)
    np.save(output_path, contact_matrix)

    return contact_matrix

def get_data(train_file_path, test_file_path, mode, train_label, test_label, window_size):
    all_train_data = []
    all_test_data = []
    
    train_data_list_path = sorted(os.listdir(train_file_path))
    test_data_list_path = sorted(os.listdir(test_file_path))

    
    for rna_train_data_path in train_data_list_path:
        train_data_list = os.path.join(train_file_path, rna_train_data_path)
        train_data = pd.read_csv(train_data_list, header=None)
        for i in range(train_data.shape[0]):
            all_train_data.append(list(train_data.iloc[i, 1:-2]))
    with open(train_label, 'rb') as f:
        train_labels = pickle.load(f)
        train_labels = [item for sublist in train_labels for item in sublist]

    
    for rna_test_data_path in test_data_list_path:
        test_data_list = os.path.join(test_file_path, rna_test_data_path)
        test_data = pd.read_csv(test_data_list, header=None)
        for i in range(test_data.shape[0]):
            all_test_data.append(list(test_data.iloc[i, 1:-2]))
    with open(test_label, 'rb') as f:
        test_labels = pickle.load(f)
        test_labels = [item for sublist in test_labels for item in sublist]
    if mode == 'apo':
        all_test_data = []
        for pdb in sorted(os.listdir("../data/apo_fastas")):
            adj_matrix = calculate_contact_matrix(f'../data/apo_pdb/{pdb[:4]}.pdb', chain_id=pdb[4], cutoff=8.0,
                                                  mode="rep",
                                                  output_dir="./ad")
            topology = get_node_topology(adj_matrix)
            topology = np.array(topology)
            topology = sliding_window_context(np.array(topology), window_size)
            all_test_data.append(topology)
        all_test_data = np.vstack(all_test_data)
        with open(test_label, 'rb') as f:
            test_labels = pickle.load(f)
            test_labels = [item for sublist in test_labels for item in sublist]
    return np.array(all_train_data),np.array(train_labels), np.array(all_test_data),np.array(test_labels)


def get_data_predict(pdb_file_path,train_fastas_path,test_fastas_path,train_label_path,test_label_path, window_size):
    train_fastas = os.listdir(train_fastas_path)
    test_fastas = os.listdir(test_fastas_path)
    train_fastas = sorted(train_fastas)
    test_fastas = sorted(test_fastas)
    all_train_data  = []
    all_test_data = []
    for fasta_name in train_fastas:
        
        pdb_name = f"{fasta_name[:4]}.pdb"
        pdb_file = os.path.join(pdb_file_path, pdb_name)
        adj_matrix = calculate_contact_matrix(pdb_file, chain_id=fasta_name[4], cutoff=8.0,
                                              mode="test17",
                                              output_dir="./ad")
        topology = get_node_topology(adj_matrix)
        topology = np.array(topology)
        topology = sliding_window_context(np.array(topology), window_size)
        all_train_data.extend(topology)
    for fasta_name in test_fastas:
        
        pdb_name = f"{fasta_name[:4]}.pdb"
        pdb_file = os.path.join(pdb_file_path, pdb_name)
        adj_matrix = calculate_contact_matrix(pdb_file, chain_id=fasta_name[4], cutoff=8.0,
                                              mode="test17",
                                              output_dir="./ad")
        topology = get_node_topology(adj_matrix)
        topology = np.array(topology)
        topology = sliding_window_context(np.array(topology), window_size)
        all_test_data.extend(topology)
        with open(train_label_path, 'rb') as f:
            train_labels = pickle.load(f)
            train_labels = [item for sublist in train_labels for item in sublist]
        with open(test_label_path, 'rb') as f:
            test_labels = pickle.load(f)
            test_labels = [item for sublist in test_labels for item in sublist]
    return np.array(all_train_data),np.array(train_labels), np.array(all_test_data),np.array(test_labels)

def rnet(rf_train_data_file_path,rf_test_data_file_path,lgb_train_data_file_path, lgb_test_data_file_path,xgb_train_data_file_path,xgb_test_data_file_path,mode,test_labels,window_size,seed):

    rf_data_train, rf_train_labels, rf_data_test, rf_test_labels = get_data(rf_train_data_file_path,rf_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[0])
    lgb_data_train, lgb_train_labels, lgb_data_test, lgb_test_labels= get_data(lgb_train_data_file_path, lgb_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[1])
    xgb_data_train, xgb_train_labels, xgb_data_test, xgb_test_labels = get_data(xgb_train_data_file_path,xgb_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[2])


    rf_model = RandomForestClassifier(n_estimators=200, random_state=seed)
    lgb_model = LGBMClassifier(max_bin=60, num_leaves=20, random_state=seed)
    xgb_model = XGBClassifier(learning_rate=0.07, gamma=1, random_state=seed)

    
    rf_model.fit(rf_data_train, rf_train_labels)
    lgb_model.fit(lgb_data_train, lgb_train_labels)
    xgb_model.fit(xgb_data_train, xgb_train_labels)

    rf_pred = rf_model.predict(rf_data_test)
    lgb_pred = lgb_model.predict(lgb_data_test)
    xgb_pred = xgb_model.predict(xgb_data_test)

    predictions = np.vstack([rf_pred, lgb_pred, xgb_pred]).T

    voting_pred = []

    
    rf_proba = rf_model.predict_proba(rf_data_test)[:, 1]  
    lgb_proba = lgb_model.predict_proba(lgb_data_test)[:, 1]
    xgb_proba = xgb_model.predict_proba(xgb_data_test)[:, 1]

    
    
    
    

    
    
    voting_proba = (rf_proba + lgb_proba + xgb_proba) / 3

    for i in range(len(predictions)):
        
        vote = np.bincount(predictions[i]).argmax()
        voting_pred.append(vote)

    
    
    
    
    auc = roc_auc_score(rf_test_labels, voting_proba)
    accuracy = accuracy_score(rf_test_labels, voting_pred)
    precision = precision_score(rf_test_labels, voting_pred)
    recall = recall_score(rf_test_labels, voting_pred)
    f1 = f1_score(rf_test_labels, voting_pred)
    mcc = matthews_corrcoef(rf_test_labels, voting_pred)

    return rf_test_labels, voting_pred, voting_proba


def rnet_conformation(rf_train_data_file_path,rf_test_data_file_path,lgb_train_data_file_path, lgb_test_data_file_path,xgb_train_data_file_path,xgb_test_data_file_path,mode,test_labels,window_size):
    rf_data_train, rf_train_labels, rf_data_test, rf_test_labels = get_data(rf_train_data_file_path,rf_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[0])
    lgb_data_train, lgb_train_labels, lgb_data_test, lgb_test_labels= get_data(lgb_train_data_file_path, lgb_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[1])
    xgb_data_train, xgb_train_labels, xgb_data_test, xgb_test_labels = get_data(xgb_train_data_file_path,xgb_test_data_file_path,mode,"../all_label/label/train_label.pkl",test_labels,window_size[2])

    
    rf_model = RandomForestClassifier(n_estimators=200, random_state=12)
    lgb_model = LGBMClassifier(max_bin=60, num_leaves=20, random_state=12)
    xgb_model = XGBClassifier(learning_rate=0.07, gamma=1, random_state=12)

    
    rf_model.fit(rf_data_train, rf_train_labels)
    lgb_model.fit(lgb_data_train, lgb_train_labels)
    xgb_model.fit(xgb_data_train, xgb_train_labels)

    
    all_auc_data1 = []
    all_auc_data2 = []
    all_auc_data3 = []

    
    with open("../rnasite/label/rnasite_conformation_net.pkl", "rb") as f:
        data1, data2, data3 = pickle.load(f)

    
    combined_ground_truth = []
    combined_probabilities = []  

    
    print("Evaluating data1:")
    for features, label in zip(data1[0], data1[1]):
        features_rf = sliding_window_context(features[:, -10:-5], window_size[0])
        features_lgb = sliding_window_context(features[:, -10:-5], window_size[1])
        features_xgb = sliding_window_context(features[:, -10:-5], window_size[2])

        rf_pred = rf_model.predict(features_rf)
        lgb_pred = lgb_model.predict(features_lgb)
        xgb_pred = xgb_model.predict(features_xgb)

        predictions = np.vstack([rf_pred, lgb_pred, xgb_pred]).T

        
        rf_proba = rf_model.predict_proba(features_rf)[:, 1]
        lgb_proba = lgb_model.predict_proba(features_lgb)[:, 1]
        xgb_proba = xgb_model.predict_proba(features_xgb)[:, 1]

        
        voting_proba = (rf_proba + lgb_proba + xgb_proba) / 3

        voting_pred = []
        for i in range(len(predictions)):
            vote = np.bincount(predictions[i]).argmax()  
            voting_pred.append(vote)

        
        auc = roc_auc_score(label, voting_proba)
        accuracy = accuracy_score(label, voting_pred)
        precision = precision_score(label, voting_pred)
        recall = recall_score(label, voting_pred)
        f1 = f1_score(label, voting_pred)
        mcc = matthews_corrcoef(label, voting_pred)

        all_auc_data1.append(auc)
        combined_ground_truth.extend(label)
        combined_probabilities.extend(voting_proba)  

        
        print(f"Accuracy: {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1 Score: {f1:.4f}")
        print(f"MCC: {mcc:.4f}")
        print(f"AUC (Voting): {auc:.4f}")

    
    print("\nEvaluating data2:")
    for features, label in zip(data2[0], data2[1]):
        features_rf = sliding_window_context(features[:, -10:-5], window_size[0])
        features_lgb = sliding_window_context(features[:, -10:-5], window_size[1])
        features_xgb = sliding_window_context(features[:, -10:-5], window_size[2])

        rf_pred = rf_model.predict(features_rf)
        lgb_pred = lgb_model.predict(features_lgb)
        xgb_pred = xgb_model.predict(features_xgb)

        predictions = np.vstack([rf_pred, lgb_pred, xgb_pred]).T

        
        rf_proba = rf_model.predict_proba(features_rf)[:, 1]
        lgb_proba = lgb_model.predict_proba(features_lgb)[:, 1]
        xgb_proba = xgb_model.predict_proba(features_xgb)[:, 1]

        
        voting_proba = (rf_proba + lgb_proba + xgb_proba) / 3

        voting_pred = []
        for i in range(len(predictions)):
            vote = np.bincount(predictions[i]).argmax()  
            voting_pred.append(vote)

        
        auc = roc_auc_score(label, voting_proba)
        accuracy = accuracy_score(label, voting_pred)
        precision = precision_score(label, voting_pred)
        recall = recall_score(label, voting_pred)
        f1 = f1_score(label, voting_pred)
        mcc = matthews_corrcoef(label, voting_pred)

        all_auc_data2.append(auc)
        combined_ground_truth.extend(label)
        combined_probabilities.extend(voting_proba)

        
        print(f"Accuracy: {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1 Score: {f1:.4f}")
        print(f"MCC: {mcc:.4f}")
        print(f"AUC (Voting): {auc:.4f}")

    
    print("\nEvaluating data3:")
    for features, label in zip(data3[0], data3[1]):
        features_rf = sliding_window_context(features[:, -10:-5], window_size[0])
        features_lgb = sliding_window_context(features[:, -10:-5], window_size[1])
        features_xgb = sliding_window_context(features[:, -10:-5], window_size[2])

        rf_pred = rf_model.predict(features_rf)
        lgb_pred = lgb_model.predict(features_lgb)
        xgb_pred = xgb_model.predict(features_xgb)

        predictions = np.vstack([rf_pred, lgb_pred, xgb_pred]).T

        
        rf_proba = rf_model.predict_proba(features_rf)[:, 1]
        lgb_proba = lgb_model.predict_proba(features_lgb)[:, 1]
        xgb_proba = xgb_model.predict_proba(features_xgb)[:, 1]

        
        voting_proba = (rf_proba + lgb_proba + xgb_proba) / 3

        voting_pred = []
        for i in range(len(predictions)):
            vote = np.bincount(predictions[i]).argmax()  
            voting_pred.append(vote)

        
        auc = roc_auc_score(label, voting_proba)
        accuracy = accuracy_score(label, voting_pred)
        precision = precision_score(label, voting_pred)
        recall = recall_score(label, voting_pred)
        f1 = f1_score(label, voting_pred)
        mcc = matthews_corrcoef(label, voting_pred)

        all_auc_data3.append(auc)
        combined_ground_truth.extend(label)
        combined_probabilities.extend(voting_proba)

        
        print(f"Accuracy: {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1 Score: {f1:.4f}")
        print(f"MCC: {mcc:.4f}")
        print(f"AUC (Voting): {auc:.4f}")

    threshold = 0.5
    
    print("\nEvaluating combined data (all datasets):")
    accuracy_combined = accuracy_score(combined_ground_truth, np.array(combined_probabilities) > threshold)
    precision_combined = precision_score(combined_ground_truth, np.array(combined_probabilities) > threshold)
    recall_combined = recall_score(combined_ground_truth, np.array(combined_probabilities) > threshold)
    f1_combined = f1_score(combined_ground_truth, np.array(combined_probabilities) > threshold)
    auc_combined = roc_auc_score(combined_ground_truth, combined_probabilities)
    mcc_combined = matthews_corrcoef(combined_ground_truth, np.array(combined_probabilities) > threshold)

    
    print(
        f"Combined Accuracy: {accuracy_combined:.3f}, Combined Precision: {precision_combined:.3f}, Combined Recall: {recall_combined:.3f}")
    print(f"Combined F1: {f1_combined:.3f}, Combined MCC: {mcc_combined:.3f}, Combined AUC: {auc_combined:.3f}")

    
    print(f"AUC for data1: {np.mean(all_auc_data1):.3f}")
    print(f"AUC for data2: {np.mean(all_auc_data2):.3f}")
    print(f"AUC for data3: {np.mean(all_auc_data3):.3f}")


def rnet_predict(pdb_file_path,train_fastas_path,test_fastas_path,train_label_path,test_label_path, window_size, random_state):
    rf_data_train, rf_train_labels, rf_data_test, rf_test_labels = get_data_predict(pdb_file_path,train_fastas_path,test_fastas_path,train_label_path,test_label_path, window_size[0])
    lgb_data_train, lgb_train_labels, lgb_data_test, lgb_test_labels = get_data_predict(pdb_file_path,train_fastas_path,test_fastas_path,train_label_path,test_label_path, window_size[1])
    xgb_data_train, xgb_train_labels, xgb_data_test, xgb_test_labels = get_data_predict(pdb_file_path,train_fastas_path,test_fastas_path,train_label_path,test_label_path, window_size[2])

    rf_model = RandomForestClassifier(n_estimators=200, random_state=random_state)
    lgb_model = LGBMClassifier(max_bin=60, num_leaves=20, random_state=random_state)
    xgb_model = XGBClassifier(learning_rate=0.07, gamma=1, random_state=random_state)

    
    rf_model.fit(rf_data_train, rf_train_labels)
    lgb_model.fit(lgb_data_train, lgb_train_labels)
    xgb_model.fit(xgb_data_train, xgb_train_labels)

    rf_pred = rf_model.predict(rf_data_test)
    lgb_pred = lgb_model.predict(lgb_data_test)
    xgb_pred = xgb_model.predict(xgb_data_test)

    predictions = np.vstack([rf_pred, lgb_pred, xgb_pred]).T

    voting_pred = []

    
    rf_proba = rf_model.predict_proba(rf_data_test)[:, 1]  
    lgb_proba = lgb_model.predict_proba(lgb_data_test)[:, 1]
    xgb_proba = xgb_model.predict_proba(xgb_data_test)[:, 1]

    
    
    
    

    
    
    voting_proba = (rf_proba + lgb_proba + xgb_proba) / 3

    for i in range(len(predictions)):
        
        vote = np.bincount(predictions[i]).argmax()
        voting_pred.append(vote)

    
    
    
    
    auc = roc_auc_score(rf_test_labels, voting_proba)
    accuracy = accuracy_score(rf_test_labels, voting_pred)
    precision = precision_score(rf_test_labels, voting_pred)
    recall = recall_score(rf_test_labels, voting_pred)
    f1 = f1_score(rf_test_labels, voting_pred)
    mcc = matthews_corrcoef(rf_test_labels, voting_pred)

    return rf_test_labels, voting_pred, voting_proba

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    




















def run_with_seeds(num_runs=5, output_path="../results/RNet/test17_results.json"):
    import random
    import numpy as np
    from sklearn.metrics import (
        accuracy_score, precision_score, recall_score,
        f1_score, matthews_corrcoef, roc_auc_score,
        average_precision_score, confusion_matrix
    )

    
    seeds = [1, 2, 3, 4, 5][:num_runs]
    print(f"\nRandom Seeds: {seeds}\n")

    
    def bacc_score(y_true, y_pred):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        sensitivity = tp / (tp + fn + 1e-8)
        specificity = tn / (tn + fp + 1e-8)
        return 0.5 * (sensitivity + specificity)

    
    metrics_all = {
        'Accuracy': [], 'Precision': [], 'Recall': [], 'F1': [],
        'MCC': [], 'AUC': [], 'AUPR': [], 'BACC': []
    }

    
    for seed in seeds:
        print(f"\n================ SEED {seed} ================\n")
        np.random.seed(seed)
        random.seed(seed)

        
        y_true, y_pred, y_prob = rnet("./data/test17/train_dataset_w3","./data/test17/test_dataset_w3","./data/test17/train_dataset_w7","./data/test17/test_dataset_w7","./data/test17/train_dataset_w8","./data/test17/test_dataset_w8","reproduction","../all_label/label/test17_labels.pkl",[3,7,8],seed)

        if len(y_true) != 583:
            raise ValueError(f"Test17 must contain 583 nucleotides, got {len(y_true)}")


        
        acc = accuracy_score(y_true, y_pred)
        pre = precision_score(y_true, y_pred)
        rec = recall_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)
        mcc = matthews_corrcoef(y_true, y_pred)
        auc = roc_auc_score(y_true, y_prob)
        aupr = average_precision_score(y_true, y_prob)
        bacc = bacc_score(y_true, y_pred)

        
        metrics_all['Accuracy'].append(acc)
        metrics_all['Precision'].append(pre)
        metrics_all['Recall'].append(rec)
        metrics_all['F1'].append(f1)
        metrics_all['MCC'].append(mcc)
        metrics_all['AUC'].append(auc)
        metrics_all['AUPR'].append(aupr)
        metrics_all['BACC'].append(bacc)

        print(f"Seed {seed}: "
              f"ACC={acc:.4f}, PRE={pre:.4f}, REC={rec:.4f}, "
              f"F1={f1:.4f}, MCC={mcc:.4f}, AUC={auc:.4f}, "
              f"AUPR={aupr:.4f}, BACC={bacc:.4f}")

    
    print("\n============== Final Results (mean ± std) ==============")
    for k, v in metrics_all.items():
        print(f"{k:8s}: {np.mean(v):.4f} ± {np.std(v):.4f}")

    summary = {
        key: {"mean": float(np.mean(values)), "std": float(np.std(values))}
        for key, values in metrics_all.items()
    }
    payload = {
        "method": "RNet",
        "dataset": "Test17",
        "rna_count": 17,
        "nucleotide_count": 583,
        "seeds": seeds,
        "per_seed": {key: [float(value) for value in values] for key, values in metrics_all.items()},
        "summary": summary,
    }
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f_out:
        json.dump(payload, f_out, indent=2)
    print(f"Test17 results saved to: {output_path}")

    return metrics_all


if __name__ == "__main__":
    run_with_seeds()
