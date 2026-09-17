import os
import sys
import math
import time
import json
import pickle
import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt
from sklearn.metrics import precision_recall_curve
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable
import torch.utils.data.sampler as sampler
import torch.optim as optim
from sklearn.model_selection import KFold

from RL_bind import RLBind
from dataset import DataSet
from evaluation import compute_roc, compute_mcc, micro_score, compute_performance
from utils import *
from utils import DefinedConfig

defconstant = DefinedConfig()
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix,\
    matthews_corrcoef, auc, precision_recall_curve

device = torch.device("cuda:0")


class AverageMeter(object):
    """
    Computes and stores the average and current value
    Copied from: https://github.com/pytorch/examples/blob/master/imagenet/main.py
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


def params_init(m):
    if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
        nn.init.normal_(m.weight.data, mean=0, std=min(1.0 / math.sqrt(m.weight.data.shape[-1]), 0.1))
        nn.init.constant_(m.bias, 0)


def t_epoch(model, load, optimizer, epoch, epochs, trainnum=None):
    print_freq = 10
    global threads
    model.train()
    losses = AverageMeter()

    for batch_indx, (global_feats, local_feats, labels, nucle_idx, rna_id) in enumerate(load):  
        with torch.no_grad():
            if torch.cuda.is_available():
                global_vert = torch.autograd.Variable(global_feats.cuda().float())
                local_vert = torch.autograd.Variable(local_feats.cuda().float())
                label_vert = torch.autograd.Variable(labels.cuda().float())
            else:
                global_vert = torch.autograd.Variable(global_feats.float())
                local_vert = torch.autograd.Variable(local_feats.float())
                label_vert = torch.autograd.Variable(labels.float())

        batch_size = global_feats.size(0)
        output = model(global_vert, local_vert)
        output = torch.cat(tuple(output), 0)
        loss = torch.nn.functional.binary_cross_entropy(output, label_vert).cuda()
        pred_value = output.ge(threads)
        pred_value = pred_value + 0
        MiP, MiR, MiF, PNum, RNum = micro_score(pred_value.data.cpu().numpy(), label_vert.data.cpu().numpy())

        losses.update(loss.item(), batch_size)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch_indx % print_freq == 0:
            performance = '\t'.join([
                'Epoch: %d/%d' % (epoch, epochs),
                'Iter: %d/%d' % (batch_indx, len(load)),
                'Loss: %0.4f' % (losses.avg),
                'MiP: %0.4f' % (MiP),
                'MiR: %0.4f' % (MiR),
                'MiF: %0.4f' % (MiF)
            ])
    return losses.avg


def v_epoch(model, load, is_test=True, validnum=None):
    print_freq = 1
    global threads
    
    losses = AverageMeter()
    model.eval()

    v_labels = []
    v_predicts = []
    nucle_nums = []
    rna_names = []

    for batch_indx, (global_feats, local_feats, labels, nucle_id, rna_id) in enumerate(load):
        with torch.no_grad():
            if torch.cuda.is_available():
                global_vert = torch.autograd.Variable(global_feats.cuda().float())
                local_vert = torch.autograd.Variable(local_feats.cuda().float())
                label_vert = torch.autograd.Variable(labels.cuda().float())
            else:
                globaa_vert = torch.autograd.Variable(global_feats.float())
                local_vert = torch.autograd.Variable(local_feats.float())
                label_vert = torch.autograd.Variable(labels.float())

        batch_size = global_feats.size(0)
        output = model(global_vert, local_vert)
        output = torch.cat(tuple(output), 0)

        loss = torch.nn.functional.binary_cross_entropy(output, label_vert).cuda()
        losses.update(loss.item(), batch_size)

        if batch_indx % print_freq == 0:
            performance = '\t'.join([
                'Test' if is_test else 'Valid',
                'Iter: [%d/%d]' % (batch_indx + 1, len(load)),
                'Loss %0.4f' % (losses.avg),
            ])

        v_labels.append(labels.numpy())
        v_predicts.append(output.data.cpu().numpy())
        nucle_nums.append(nucle_id.numpy())
        rna_names.append(rna_id)

    v_labels = np.concatenate(v_labels, axis=0)
    print(v_labels.size)
    v_predicts = np.concatenate(v_predicts, axis=0)
    nucle_nums = np.concatenate(nucle_nums, axis=0)
    rna_names = np.concatenate(rna_names, axis=0)

    auc = compute_roc(v_predicts, v_labels)
    p_max, r_max, t_max, predictions_max = compute_performance(v_predicts, v_labels)
    mcc = compute_mcc(predictions_max, v_labels)
    return losses.avg, p_max, r_max, auc, t_max, mcc, predictions_max, v_labels, v_predicts, nucle_nums, rna_names


def compute_metrics(v_predicts, v_labels):
    """
    Compute multiple evaluation metrics for binary classification.

    Parameters:
    v_predicts (array-like): Predicted probabilities or binary class predictions (0 or 1).
    v_labels (array-like): True labels (0 or 1).

    Returns:
    dict: A dictionary containing various evaluation metrics.
    """
    print(v_predicts)
    for th in range(1, 10):
        
        binary_predicts = (v_predicts >= th * 0.1).astype(int)
        
        metrics = {}
        metrics['AUC'] = roc_auc_score(v_labels, v_predicts)  
        metrics['Accuracy'] = accuracy_score(v_labels, binary_predicts)
        metrics['Precision'] = precision_score(v_labels, binary_predicts)
        metrics['Recall'] = recall_score(v_labels, binary_predicts)
        metrics['F1-score'] = f1_score(v_labels, binary_predicts)
        metrics['MCC'] = matthews_corrcoef(v_labels, binary_predicts)
        metrics['Confusion Matrix'] = confusion_matrix(v_labels, binary_predicts)
        print(metrics)
        


def test_epoch(model, load, is_test=None):
    print_freq = 1
    global threads
    losses = AverageMeter()
    model.eval()

    v_labels = []
    v_predicts = []
    nucle_nums = []
    rna_names = []
    for batch_indx, (global_feats, local_feats, labels, nucle_id, rna_id) in enumerate(load):
        with torch.no_grad():
            if torch.cuda.is_available():
                global_vert = torch.autograd.Variable(global_feats.cuda().float())
                local_vert = torch.autograd.Variable(local_feats.cuda().float())
                label_vert = torch.autograd.Variable(labels.cuda().float())
            else:
                globaa_vert = torch.autograd.Variable(global_feats.float())
                local_vert = torch.autograd.Variable(local_feats.float())
                label_vert = torch.autograd.Variable(labels.float())
        batch_size = global_feats.size(0)
        output = model(global_vert, local_vert)
        output = torch.cat(tuple(output), 0)
        loss = torch.nn.functional.binary_cross_entropy(output, label_vert).cuda()
        losses.update(loss.item(), batch_size)

        if batch_indx % print_freq == 0:
            performance = '\t'.join([
                'Test19' if is_test else 'T19',
                'Iter: [%d/%d]' % (batch_indx + 1, len(load)),
                'Loss %0.4f' % (losses.avg),
            ])
        v_labels.append(labels.numpy())
        v_predicts.append(output.data.cpu().numpy())
        nucle_nums.append(nucle_id.numpy())
        rna_names.append(rna_id)
    v_labels = np.concatenate(v_labels, axis=0)
    v_predicts = np.concatenate(v_predicts, axis=0)
    nucle_nums = np.concatenate(nucle_nums, axis=0)
    rna_names = np.concatenate(rna_names, axis=0)

    precision, recall, thresholds = precision_recall_curve(v_labels, v_predicts)
    aupr = auc(recall, precision)
    auc1 = compute_roc(v_predicts, v_labels)
    p_max, r_max, t_max, predictions_max = compute_performance(v_predicts, v_labels)
    mcc = compute_mcc(predictions_max, v_labels)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    

    return losses.avg, p_max, r_max, auc1, t_max, mcc, predictions_max, v_labels, v_predicts, nucle_nums, rna_names, aupr


def train(model, train_dataset_all, test19, save=None, batch_size=32, train_number=6, epochs=30, train_files=None,
          test19_files=None):
    cutoff_seq_len = defconstant.cutoff_seq_len
    global threads
    global splite_rate
    with open(train_files, 'rb') as f_rna:
        rna_list = pickle.load(f_rna)
    
    with open(test19_files, 'rb') as f_test19:
        test19_list = pickle.load(f_test19)
    nucleo_samples = len(rna_list)
    split_id = int(splite_rate * nucleo_samples)
    valid_id = nucleo_samples - split_id
    np.random.shuffle(rna_list)
    
    
    
    
    
    
    
    
    
    from sklearn.model_selection import train_test_split
    from torch.utils.data import DataLoader, Subset

    
    indices = list(range(len(train_dataset_all)))  
    train_idx, valid_idx = train_test_split(indices, test_size=0.2, random_state=11)

    
    train_dataset = Subset(train_dataset_all, train_idx)
    validation_dataset = Subset(train_dataset_all, valid_idx)

    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True,
                              pin_memory=torch.cuda.is_available(), num_workers=6, drop_last=False)

    valid_loader = DataLoader(validation_dataset, batch_size=batch_size, shuffle=False,
                              pin_memory=torch.cuda.is_available(), num_workers=6, drop_last=False)

    
    test19_samples = sampler.SubsetRandomSampler(test19_list)

    test19_loader = torch.utils.data.DataLoader(test19, batch_size=batch_size, sampler=test19_samples,
                                                pin_memory=(torch.cuda.is_available()), num_workers=6, drop_last=False)
    if torch.cuda.is_available():
        model = model.cuda()
    
    model_wrapper = model
    optimizer = torch.optim.Adam(model_wrapper.parameters(), lr=0.0001)
    best_p_max = 0.0
    best_r_max = 0.0
    best_auc = 0.0
    
    train_losses = []

    for epoch in range(epochs):
        
        t_loss = t_epoch(model=model_wrapper, load=train_loader, optimizer=optimizer, epoch=epoch, epochs=epochs,
                         trainnum=split_id)
        train_losses.append(t_loss)  

        
        v_loss, p_max, r_max, auc, t_max, mcc, v_predictions, v_labels, v_predicts, nucle_nums, rna_names = v_epoch(
            model=model_wrapper, load=valid_loader, is_test=(not valid_loader), validnum=valid_id)

        
        test19_loss, p_max19, r_max19, auc19, t_max19, mcc19, predictions19, v_labels19, v_predicts19, nucle_nums19, rna_names19, aupr = test_epoch(
            model=model_wrapper, load=test19_loader, is_test=test19_loader)

        
        

        
        if auc > best_auc:
            print(f"best:{epoch + 1}")
            best_auc = auc
            threadhold = t_max
            torch.save(model.state_dict(), os.path.join(save, 'best_model.dat'))

    
    
    
    
    
    
    
    








































if __name__ == '__main__':
    import argparse
    import random
    from scipy import stats

    parser = argparse.ArgumentParser(description="Train RLBind on Train60 and evaluate Test17")
    parser.add_argument("--epochs", type=int, default=60)
    parser.add_argument("--seeds", type=int, nargs="+", default=[8124, 27045, 58392, 17765, 44322])
    args = parser.parse_args()

    def set_seed(seed):
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    def bacc_score(y_true, y_pred):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        sensitivity = tp / (tp + fn + 1e-8)
        specificity = tn / (tn + fp + 1e-8)
        return 0.5 * (sensitivity + specificity)

    
    batch_size = defconstant.batch_size
    splite_rate = defconstant.splite_rate
    epochs = args.epochs
    class_nums = defconstant.class_nums
    threads = 0.577
    path_dir = '../results/RLBind'
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)

    
    
    
    train_global_files = ["../data/RLBind/data_T60.pkl"]
    train_global_labels = ["../data/RLBind/label_T60.pkl"]
    train_local_files = ["../data/RLBind/mot11_T60.pkl"]
    train_all_nucleotide_files = "../data/RLBind/train_all.pkl"
    train_index_files = "../data/RLBind/train_index.pkl"

    test19_global_files = ["../data/RLBind/delete_6ez0/delete_data/data_T18.pkl"]
    test19_global_labels = ["../data/RLBind/delete_6ez0/delete_data/label_T18.pkl"]
    test19_local_files = ["../data/RLBind/delete_6ez0/delete_data/mot11_T18.pkl"]
    test19_all_nucleotide_files = "../data/RLBind/delete_6ez0/delete_data/test18_all.pkl"
    test19_index_files = "../data/RLBind/delete_6ez0/delete_data/test18_index.pkl"

    
    
    
    
    

    
    seeds = args.seeds
    metrics_all = {
        'Accuracy': [], 'Precision': [], 'Recall': [], 'F1': [],
        'MCC': [], 'AUC': [], 'AUPR': [], 'BACC': []
    }
    evaluated_nucleotide_count = None

    for seed in seeds:
        print(f"\n================ SEED {seed} ================\n")
        set_seed(seed)

        
        train_data = DataSet(train_global_files, train_local_files, train_global_labels, train_all_nucleotide_files)
        test19_data = DataSet(test19_global_files, test19_local_files, test19_global_labels, test19_all_nucleotide_files)

        model = RLBind()
        model.apply(params_init)

        
        train(model, train_data, test19_data, path_dir, batch_size, 0, epochs, train_index_files, test19_index_files)

        
        model.load_state_dict(torch.load(os.path.join(path_dir, 'best_model.dat')))
        model = model.cuda() if torch.cuda.is_available() else model

        
        with open(test19_index_files, 'rb') as f_test19:
            test19_list = pickle.load(f_test19)
        test19_samples = sampler.SubsetRandomSampler(test19_list)
        test19_loader = torch.utils.data.DataLoader(
            test19_data, batch_size=batch_size, sampler=test19_samples,
            pin_memory=torch.cuda.is_available(), num_workers=6, drop_last=False
        )

        
        loss, p_max, r_max, auc_value, t_max, mcc, pred_bin, y_true, y_prob, _, _, aupr = test_epoch(
            model, test19_loader, is_test=True
        )
        if evaluated_nucleotide_count is None:
            evaluated_nucleotide_count = len(y_true)
            if evaluated_nucleotide_count != 583:
                print(
                    f"WARNING: RLBind's supplied Test17 index evaluates "
                    f"{evaluated_nucleotide_count}/583 nucleotides."
                )
        elif len(y_true) != evaluated_nucleotide_count:
            raise ValueError("Inconsistent evaluated nucleotide count across seeds")

        y_bin = (y_prob >= t_max).astype(int)
        acc = accuracy_score(y_true, y_bin)
        pre = precision_score(y_true, y_bin)
        rec = recall_score(y_true, y_bin)
        f1 = f1_score(y_true, y_bin)
        bacc = bacc_score(y_true, y_bin)

        
        metrics_all['Accuracy'].append(acc)
        metrics_all['Precision'].append(pre)
        metrics_all['Recall'].append(rec)
        metrics_all['F1'].append(f1)
        metrics_all['MCC'].append(mcc)
        metrics_all['AUC'].append(auc_value)
        metrics_all['AUPR'].append(aupr)
        metrics_all['BACC'].append(bacc)

        print(f"[SEED {seed}] Acc={acc:.4f}, Pre={pre:.4f}, Rec={rec:.4f}, F1={f1:.4f}, "
              f"MCC={mcc:.4f}, AUC={auc_value:.4f}, AUPR={aupr:.4f}, BACC={bacc:.4f}")

    
    print("\n================ FINAL RESULTS ================\n")
    for key, values in metrics_all.items():
        mean = np.mean(values)
        ci = 0.0 if len(values) < 2 else stats.sem(values) * stats.t.ppf((1 + 0.95) / 2, len(values) - 1)
        print(f"{key:10s}: {mean:.4f} ± {ci:.4f} (95% CI)")

    
    result_file = os.path.join(path_dir, "rlbind_test17_multi_seed_results.txt")
    with open(result_file, "w") as f_out:
        def log(msg):
            print(msg)
            f_out.write(msg + "\n")

        log("\n================ FINAL RESULTS ================\n")

        for key, values in metrics_all.items():
            mean = np.mean(values)
            ci = 0.0 if len(values) < 2 else stats.sem(values) * stats.t.ppf((1 + 0.95) / 2, len(values) - 1)
            log(f"{key:10s}: {mean:.4f} ± {ci:.4f} (95% CI)")

        log("\n================ ALL SEED DETAILS ================\n")
        for i, seed in enumerate(seeds):
            log(f"[SEED {seed}] Acc={metrics_all['Accuracy'][i]:.4f}, "
                f"Pre={metrics_all['Precision'][i]:.4f}, "
                f"Rec={metrics_all['Recall'][i]:.4f}, "
                f"F1={metrics_all['F1'][i]:.4f}, "
                f"MCC={metrics_all['MCC'][i]:.4f}, "
                f"AUC={metrics_all['AUC'][i]:.4f}, "
                f"AUPR={metrics_all['AUPR'][i]:.4f}, "
                f"BACC={metrics_all['BACC'][i]:.4f}")

    print(f"\nMulti-seed results saved to: {result_file}\n")

    json_file = os.path.join(path_dir, "test17_results.json")
    payload = {
        "method": "RLBind", "dataset": "Test17", "rna_count": 17,
        "canonical_nucleotide_count": 583,
        "evaluated_nucleotide_count": evaluated_nucleotide_count,
        "epochs": epochs, "seeds": seeds,
        "per_seed": {key: [float(value) for value in values] for key, values in metrics_all.items()},
        "summary": {
            key: {"mean": float(np.mean(values)), "std": float(np.std(values))}
            for key, values in metrics_all.items()
        },
    }
    with open(json_file, "w", encoding="utf-8") as f_json:
        json.dump(payload, f_json, indent=2)
    print(f"JSON results saved to: {json_file}")
