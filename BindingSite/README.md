# Reproducing RNA–Small-Molecule Binding-Site Methods

This directory consolidates nine RNA–small-molecule binding-site prediction methods: **RLBind, RNet, RBind, Rsite, Rsite2, RNAsite, RNABind, MultiModRLBP, and MVRBind**. A standardized `run_test17.py` entry point is provided for each method.

> Last verified: 2026-09-17. Public repositories, websites, and dependency packages may change after this date.

## 1. Test17 Data

All methods use the same Test17 set. It contains 17 RNA chains and 583 nucleotides and was obtained by removing `6EZ0A` from the public Test18 set. Labels and evaluation order for every method follow this version.

Shared data are stored in `data/common/`:

- `test17.fasta`: combined Test17 FASTA file;
- `test17_fastas/`: 17 FASTA files named by PDB ID and chain;
- `train60_fastas/`: FASTA files for the Train60 training set;
- `pdbFiles/`: standardized RNA structures;
- `all_label/label/train_label.pkl`: Train60 labels;
- `all_label/label/test17_labels.pkl`: Test17 labels.

Method-specific features are stored in `data/<method>/` and are not shared across methods. The `all_label` entry and the `data` or label entries under RBind, Rsite, RNet, and MultiModRLBP are Windows junctions that point to the physical data in this directory. Recreate these links if the package is copied through a ZIP archive or a file system that does not preserve junctions.

## 2. Environments

The nine methods require three Conda environments because their Python, Biopython, PyTorch, and CUDA constraints are incompatible.

| Environment | Methods | Python | PyTorch | CUDA |
|---|---|---:|---:|---:|
| `core` | RLBind, RNet, RBind, Rsite, Rsite2, RNAsite | 3.9.21 | 2.1.0 | 11.8 |
| `MVRBind` | MVRBind, RNABind | 3.10.18 | 2.2.0 | 12.1 |
| `Mul` | MultiModRLBP | 3.7.12 | 1.13.1 | 11.6 |

Create the environments from the `BindingSite` directory:

```bash
conda env create -f environments/core.yml
conda env create -f environments/mvrbind.yml
conda env create -f environments/multimodrlbp.yml
```

MultiModRLBP includes an `alignment_C` extension built for Windows CPython 3.7. Rebuild it in the `Mul` environment on Linux or when the included binary is incompatible:

```bash
conda activate Mul
cd MultiModRLBP/RnaBert
python setup.py build_ext --inplace
cd ../..
```

See `ENVIRONMENTS.md` for more detailed version and GPU information.

## 3. Unified Execution and Validation

Run the following commands from the `BindingSite` directory.

Run all methods:

```bash
python run_all_test17.py --methods all --continue-on-error
```

Run selected methods only:

```bash
python run_all_test17.py --methods rbind rsite rsite2 mvrbind
```

The unified runner automatically selects the `core`, `MVRBind`, or `Mul` environment for each method. If Conda is not on `PATH`, set `CONDA_EXE` or pass `--conda /path/to/conda`.

Check directory structure, Python syntax, Train60/Test17 counts, label lengths, and result JSON files:

```bash
python validate_review.py
```

Aggregate generated results:

```bash
python collect_results.py
```

The summary table is written to `results/test17_summary.csv`. Complete output for each method is stored in `results/<method>/test17_results.json`.

## 4. RLBind

- Paper: [RLBind: a deep learning method to predict RNA–ligand binding sites](https://doi.org/10.1093/bib/bbac486)
- Original source code and data: [KailiWang1/RLBind](https://github.com/KailiWang1/RLBind)
- Local code: `RLBind/`
- Local data: `data/RLBind/`
- Environment: `core`

RLBind combines global RNA sequence features with local neighborhood features. The local Test17 entry point retrains the model on Train60 using 60 epochs and five random seeds by default.

```bash
conda run --no-capture-output -n core python RLBind/run_test17.py
```

Run a one-epoch pipeline check:

```bash
conda run --no-capture-output -n core python RLBind/run_test17.py --epochs 1 --seeds 8124
```

## 5. RNet

- Paper: [RNet: a network strategy to predict RNA binding preferences](https://doi.org/10.1093/bib/bbad482)
- Source-code and data page reported by the authors: [RNetsite](http://zhaoserver.com.cn/RNet/RNet.html)
- Local code: `RNet/`
- Local data: `data/RNet/`
- Environment: `core`

RNetsite in RNet extracts local and global network features from an RNA three-dimensional contact network and predicts nucleotide-level binding sites with conventional machine-learning models. The local entry point runs five configured seeds.

```bash
conda run --no-capture-output -n core python RNet/run_test17.py
```

## 6. RBind

- Paper: [RBind: computational network method to predict RNA binding sites](https://doi.org/10.1093/bioinformatics/bty345)
- Original code and data page: [RBind](https://zhaolab.com.cn/RBind)
- Local code: `RBind/`
- Local data: `data/common/`
- Environment: `core`

RBind represents nucleotides as nodes in a contact network and predicts binding sites using topological features and statistical thresholds. The Test17 entry point performs deterministic evaluation.

```bash
conda run --no-capture-output -n core python RBind/run_test17.py
```

## 7. Rsite

- Paper: [Rsite: a computational method to identify the functional sites of noncoding RNAs](https://doi.org/10.1038/srep09179)
- Original source code and example data: [Rsite official page](https://www.cuilab.cn/rsite)
- Local code: `Rsite/`
- Local data: `data/common/`
- Environment: `core`

Rsite calculates nucleotide-to-centroid distance profiles from RNA three-dimensional coordinates and identifies candidate functional sites from extrema in the smoothed profile.

```bash
conda run --no-capture-output -n core python Rsite/run_test17.py
```

## 8. Rsite2

- Paper: [Rsite2: an efficient computational method to predict the functional sites of noncoding RNAs](https://doi.org/10.1038/srep19016)
- Original source code and data: [Rsite official page](https://www.cuilab.cn/rsite)
- Local code: `Rsite2/`
- Local data: `data/Rsite2/`
- Environment: `core`

Rsite2 replaces the three-dimensional structural requirement of Rsite with two-dimensional coordinates and distance profiles derived from RNA secondary structure. The local Test17 entry point reads precomputed `SS_NDS` files and reports the standardized metrics.

```bash
conda run --no-capture-output -n core python Rsite2/run_test17.py
```

## 9. RNAsite

- Paper: [Recognition of small molecule–RNA binding sites using RNA sequence and structure](https://doi.org/10.1093/bioinformatics/btaa1092)
- Official web server: [RNAsite](https://yanglab.qd.sdu.edu.cn/RNAsite/)
- Local code: `RNAsite/`
- Local data: `data/RNAsite/` and `data/common/`
- Environment: `core`

RNAsite combines sequence conservation, three-dimensional network topology, Laplacian norm, and SASA features and makes predictions with a random forest. The local entry point uses preprocessed Train60/Test17 features and runs five configured seeds.

```bash
conda run --no-capture-output -n core python RNAsite/run_test17.py
```

## 10. RNABind

- Paper: [Identifying RNA-small Molecule Binding Sites Using Geometric Deep Learning with Language Models](https://doi.org/10.1016/j.jmb.2025.169010)
- Original source code and data: [jaminzzz/RNABind](https://github.com/jaminzzz/RNABind)
- Local code: `RNABind/`
- Local data: `data/common/`
- Environment: `MVRBind`

RNABind uses RNA structural graphs and an EGNN. The original paper also evaluates several RNA language-model embeddings. This local package does not include large-model weights such as ERNIE-RNA. The Test17 script therefore explicitly uses the one-hot input variant of the official RNABind EGNN architecture and does not download other models implicitly. Training uses 30 epochs and five random seeds by default.

```bash
conda run --no-capture-output -n MVRBind python RNABind/run_test17.py
```

Run a one-epoch pipeline check:

```bash
conda run --no-capture-output -n MVRBind python RNABind/run_test17.py --epochs 1 --seeds 8124
```

## 11. MultiModRLBP

- Paper: [MultiModRLBP: A Deep Learning Approach for Multi-Modal RNA-Small Molecule Ligand Binding Sites Prediction](https://doi.org/10.1109/JBHI.2024.3400521)
- Original source code and data: [lennylv/MultiModRLBP](https://github.com/lennylv/MultiModRLBP)
- Local code: `MultiModRLBP/`
- Local data: `data/MultiModRLBP/`
- Environment: `Mul`

MultiModRLBP combines nucleotide-level three-dimensional features, an RNA relation graph, and RNABert sequence representations. The local data include the preprocessed features and model files required for Test17. Training uses 100 epochs by default.

```bash
conda run --no-capture-output -n Mul python MultiModRLBP/run_test17.py
```

Run a one-epoch pipeline check:

```bash
conda run --no-capture-output -n Mul python MultiModRLBP/run_test17.py --epochs 1
```

## 12. MVRBind

- Paper: [MVRBind: multi-view learning for RNA-small molecule binding site prediction](https://doi.org/10.1093/bib/bbaf489)
- Original source code and data: [cschen-y/MVRBind](https://github.com/cschen-y/MVRBind)
- Local code: `MVRBind/`
- Local data: `data/MVRBind/`
- Environment: `MVRBind`

MVRBind models RNA primary, secondary, and tertiary structures as multiple views and fuses node representations across spatial scales. The local Test17 entry point reads preprocessed graphs from `data/MVRBind/pt/`, removes `6EZ0A` from Test18, and evaluates five pretrained checkpoints in `data/MVRBind/model_parameters/`.

```bash
conda run --no-capture-output -n MVRBind python MVRBind/run_test17.py
```
