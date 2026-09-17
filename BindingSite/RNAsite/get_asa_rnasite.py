import os
import pandas as pd
import numpy as np

def restore_features_with_zeros(original, modified, modified_features):
    """
    在被删除的位置补零，同时对齐剩余的特征。

    参数：
        original: 原始序列 (字符串)
        modified: 删除字符后的序列 (字符串)
        modified_features: 删除字符后的特征列表 (二维数组，每个字符对应多个特征)

    返回：
        restored_features: 补齐后的特征列表 (二维数组，与原始序列长度一致)
    """
    deleted_positions = []
    i, j = 0, 0
    while i < len(original) and j < len(modified):
        if original[i] == modified[j]:
            i += 1
            j += 1
        else:
            deleted_positions.append(i)
            i += 1
    
    while i < len(original):
        deleted_positions.append(i)
        i += 1

    
    restored_features = []
    j = 0  
    feature_dim = len(modified_features.iloc[0]) if isinstance(modified_features, pd.DataFrame) else len(
        modified_features[0])
    for i in range(len(original)):
        if i in deleted_positions:
            
            restored_features.append([0] * feature_dim)
        else:
            
            restored_features.append(
                modified_features.iloc[j] if isinstance(modified_features, pd.DataFrame) else modified_features[j])
            j += 1

    
    return restored_features



def get_new_asa(asa_file_path="./data/train60_test18_asa",fasta_file_path = './data/test18_fastas'):
    fastas_file = sorted(os.listdir(fasta_file_path))
    all_asa = []
    
    for fasta_name in fastas_file:
        with open(f"{fasta_file_path}/{fasta_name}", 'r') as f:
            lines = f.readlines()
            o_seq = lines[1].strip()
            rna_length = len(o_seq)
        
        df = pd.read_csv(f'{asa_file_path}/{fasta_name[:5]}.csv')
        
        filtered_df = df[df['Chain'].isin([fasta_name[4]])]
        final_result = filtered_df[['Phob.A.2', 'Phil.A.2', 'SASA.A.2']]
        asa_residues_seq = ''.join(filtered_df['ResidNe'])
        
        if final_result.shape[0] != rna_length:
            final_result = restore_features_with_zeros(o_seq,asa_residues_seq,final_result)
            
            
        final_result = np.array(final_result)
        final_result = final_result[:,-1].reshape(-1,1)
        
        all_asa.append(final_result)
        
    
    return all_asa

