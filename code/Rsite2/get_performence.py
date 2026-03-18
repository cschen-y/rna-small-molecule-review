import os
import re
import pickle
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
    balanced_accuracy_score
)


def process_resullt(file_path):
    # 读取文件内容
    with open(file_path, 'r') as file:
        content = file.read()
    # 提取"Predicted sites number"后的站点信息
    predicted_sites = []
    # 使用正则表达式来查找站点
    # 假设"Predicted sites"后面的行包含了站点位置，按格式提取
    matches = re.findall(r'Site#\d+:\s*(\d+(?:,\d+)*)', content)
    # 处理每个站点位置，将其转换为整数并存入列表
    for match in matches:
        sites = match.split(',')
        for site in sites:
            predicted_sites.append(int(site))
    # 输出结果
    # print(predicted_sites)
    return predicted_sites

def get_performence(result_file_path, fastas_file_path):
    fasta_files = sorted(os.listdir(fastas_file_path))
    all_predict_labels = []
    for fasta_name in fasta_files:
        # rna_result_2_2_path = f"{result_file_path}/{fasta_name[:5]}_ss_ndc_gd_extrema_2-2.txt"
        rna_result_1_5_path = f"{result_file_path}/{fasta_name[:5]}_ss_ndS_gd_extrema_1-5.txt"
        rna_predict_label_index = process_resullt(rna_result_1_5_path)
        rna_predict_label = []
        fasta_rna = os.path.join(fastas_file_path, fasta_name)
        with open(fasta_rna, 'r') as file:
            lines = file.readlines()
            rna_length = len(lines[1].strip())
        for i in range(rna_length):
            if i+1 in rna_predict_label_index:
                rna_predict_label.append(1)
            else:
                rna_predict_label.append(0)
        all_predict_labels.extend(rna_predict_label)
    return all_predict_labels

def read_labels(label_file_path):
    labels_result = []
    with open(label_file_path, 'rb') as f:
        labels = pickle.load(f)
    for i in labels:
        labels_result.extend(i)
    return labels_result

def get_result(result_file_path, fastas_file_path, ground_true_labels):
    all_predict_labels = get_performence(result_file_path,
                    fastas_file_path)
    true_labels = read_labels(ground_true_labels)
    accuracy = accuracy_score(true_labels, all_predict_labels)
    precision = precision_score(true_labels, all_predict_labels)
    recall = recall_score(true_labels, all_predict_labels)
    f1 = f1_score(true_labels, all_predict_labels)
    auc = roc_auc_score(true_labels, all_predict_labels)  # AUC
    mcc = matthews_corrcoef(true_labels, all_predict_labels)  # MCC
    bacc = balanced_accuracy_score(true_labels, all_predict_labels)

    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"AUC: {auc:.4f}")
    print(f"MCC: {mcc:.4f}")
    print(f"bacc: {bacc:.4f}")


get_result(r"D:\Desktop\rna\code\RLBind\v2\reproduction\Rsite2\PS\SS_NDS", r"D:\Desktop\rna\code\RLBind\v2\reproduction\data\test17_fastas", "./test17_labels.pkl")
