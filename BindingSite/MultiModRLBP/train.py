import os
import sys
import argparse
import json
import torch
import random
import numpy as np

from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, matthews_corrcoef,
    roc_auc_score, average_precision_score,
    balanced_accuracy_score
)




script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '..', '..'))
sys.path.append('RnaBert')

from Model import model, learn
from RGCN.data_loading import graphloader
from RnaBert.MLM_SFP import get_config, TRAIN, BertModel, BertForMaskedLM




node_target = ['binding_small-molecule']
node_features = [
    'nt_code', 'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
    'epsilon_zeta', 'chi', 'C5prime_xyz', 'P_xyz', 'ssZp', 'Dp',
    'splay_angle', 'splay_distance', 'splay_ratio',
    'eta', 'theta', 'eta_prime', 'theta_prime',
    'eta_base', 'theta_base',
    'v0', 'v1', 'v2', 'v3', 'v4',
    'amplitude', 'phase_angle',
    'suiteness', 'filter_rmsd',
    'puckering', 'sugar_class', 'bin',
    'TotalAsa', 'PolarAsa', 'ApolarAsa'
]

device = "cuda:0" if torch.cuda.is_available() else "cpu"

seeds = [8124, 27045, 58392, 17765, 44322]




def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


@torch.no_grad()
def evaluate_all_metrics(model, loader, threshold=0.573):
    model.eval()
    device = model.current_device

    y_true_all = []
    y_score_all = []

    for batch_idx, (graph, indxs, features, seqs, seqlens, len_graphs, chain_idx) in enumerate(loader):

        graph = graph.to(device)
        features = features.to(device)

        indxs = torch.as_tensor(indxs, device=device)
        chain_idx = torch.as_tensor(chain_idx, device=device)

        seqs = torch.from_numpy(seqs.astype(np.int64)).to(device)
        len_graphs = torch.as_tensor(len_graphs, device=device)

        
        labels = graph.ndata["target"]
        labels = torch.index_select(labels, 0, indxs)

        
        outputs = model(
            graph,
            features,
            indxs,
            seqs,
            seqlens,
            len_graphs,
            chain_idx
        )

        
        pos = 0
        idx = 0
        filtered_labels = []

        for seqlen in seqlens:
            if idx in chain_idx:
                filtered_labels.append(labels[pos:pos + seqlen])
            pos += seqlen
            idx += 1

        filtered_labels = torch.cat(filtered_labels, dim=0)

        
        y_true_all.append(filtered_labels.view(-1).cpu().numpy())
        y_score_all.append(outputs.view(-1).cpu().numpy())

    
    y_true = np.concatenate(y_true_all)
    y_score = np.concatenate(y_score_all)

    
    y_score = np.nan_to_num(y_score, nan=0.0, posinf=1.0, neginf=0.0)

    
    y_pred = (y_score >= threshold).astype(int)

    
    acc = accuracy_score(y_true, y_pred)
    pre = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_true, y_pred)

    auc = roc_auc_score(y_true, y_score)
    aupr = average_precision_score(y_true, y_score)
    bacc = balanced_accuracy_score(y_true, y_pred)

    return acc, pre, rec, f1, mcc, auc, aupr, bacc





if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--train_batch", type=int, default=5)
    parser.add_argument("--test_batch", type=int, default=17)
    parser.add_argument("--result_path", type=str, default="../results/MultiModRLBP/test17_results.json")
    args = parser.parse_args()

    
    
    
    with open("data/baseline/train60.txt") as f:
        train_graphs = [l.strip().lower() for l in f if len(l.strip()) > 3]

    with open("data/baseline/test17.txt") as f:
        test_graphs = [l.strip().lower() for l in f if len(l.strip()) > 3]
    if len(test_graphs) != 17:
        raise ValueError(f"Test17 requires 17 RNA chains, got {len(test_graphs)}")

    
    
    
    metrics_all = {
        "Accuracy": [],
        "Precision": [],
        "Recall": [],
        "F1": [],
        "MCC": [],
        "AUC": [],
        "AUPR": [],
        "BACC": []
    }

    
    
    
    for seed in seeds:
        print(f"\n========== Seed {seed} ==========")
        set_seed(seed)

        
        train_dataset = graphloader.SupervisedDataset(
            data_path="data/myData_asa",
            hashing_path="data/2.5DGraph/iguana/all_graphs_annot_hash.p",
            node_features=node_features,
            node_target=node_target,
            all_graphs=train_graphs
        )

        train_loader = graphloader.GraphLoader(
            dataset=train_dataset,
            split=False,
            batch_size=args.train_batch,
            num_workers=0
        ).get_data()

        test_dataset = graphloader.SupervisedDataset(
            data_path="data/myData_asa",
            hashing_path="data/2.5DGraph/iguana/all_graphs_annot_hash.p",
            node_features=node_features,
            node_target=node_target,
            all_graphs=test_graphs
        )
        test_dataset.setNorm(train_dataset.getNorm())

        test_loader = graphloader.GraphLoader(
            dataset=test_dataset,
            split=False,
            batch_size=args.test_batch,
            num_workers=0
        ).get_data()

        
        config = get_config("RNA_bert_config.json")
        config.hidden_size = config.num_attention_heads * config.multiple

        train_helper = TRAIN(config, device)
        rna_bert = BertForMaskedLM(config, BertModel(config))
        rna_bert = train_helper.model_device(rna_bert)
        rna_bert.load_state_dict(torch.load("bert_mul_2.pth", map_location=device))

        for p in rna_bert.parameters():
            p.requires_grad = False

        
        embedder = model.RGATEmbedder(
            infeatures_dim=train_dataset.input_dim + 2,
            dims=[64, 64]
        )

        classifier = model.RGATClassifier(
            rgat_embedder=embedder,
            rbert_embedder=rna_bert,
            conv_output=False,
            return_loss=False,
            classif_dims=[train_dataset.output_dim]
        )
        classifier.deactivate_loss()
        classifier.to(device)

        optimizer = torch.optim.Adam(
            filter(lambda p: p.requires_grad, classifier.parameters()),
            lr=args.lr
        )

        routine = learn.LearningRoutine(
            num_epochs=args.epochs,
            device=device,
            test18_loader=test_loader,
            save_path=None
        )

        learn.train_supervised(
            model=classifier,
            optimizer=optimizer,
            train_loader=train_loader,
            learning_routine=routine,
        )

        
        acc, pre, rec, f1, mcc, auc, aupr, bacc = evaluate_all_metrics(
            classifier, test_loader
        )

        metrics_all['Accuracy'].append(acc)
        metrics_all['Precision'].append(pre)
        metrics_all['Recall'].append(rec)
        metrics_all['F1'].append(f1)
        metrics_all['MCC'].append(mcc)
        metrics_all['AUC'].append(auc)
        metrics_all['AUPR'].append(aupr)
        metrics_all['BACC'].append(bacc)

        print(
            f"Seed {seed}: "
            f"ACC={acc:.4f}, PRE={pre:.4f}, REC={rec:.4f}, "
            f"F1={f1:.4f}, MCC={mcc:.4f}, AUC={auc:.4f}, "
            f"AUPR={aupr:.4f}, BACC={bacc:.4f}"
        )

    
    
    
    print("\n============== Final Results (mean ± std) ==============")
    for k, v in metrics_all.items():
        print(f"{k:8s}: {np.mean(v):.4f} ± {np.std(v):.4f}")

    payload = {
        "method": "MultiModRLBP", "dataset": "Test17", "rna_count": 17,
        "epochs": args.epochs,
        "seeds": seeds,
        "per_seed": {key: [float(value) for value in values] for key, values in metrics_all.items()},
        "summary": {
            key: {"mean": float(np.mean(values)), "std": float(np.std(values))}
            for key, values in metrics_all.items()
        },
    }
    os.makedirs(os.path.dirname(args.result_path), exist_ok=True)
    with open(args.result_path, "w", encoding="utf-8") as f_out:
        json.dump(payload, f_out, indent=2)
    print(f"Test17 results saved to: {args.result_path}")
