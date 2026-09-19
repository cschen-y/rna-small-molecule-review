# Reproducing RNA–Small-Molecule Binding-Preference Methods

This document covers RNA–small-molecule binding-preference and binding-affinity methods, including **BioLLMNet, DeepRSMA, RSAPred, RLaffinity, RLASIF, SPRank, RNAmigos, R-BIND, RNALigands, and ZHMol-RLinter**. It provides links to the papers and source code, data sources, and execution steps. When executable source code is currently unavailable, only the specific reason that prevents reproduction is described.

> Last verified: 2026-09-17. Public repositories and download links may change after this date.

## 1. BioLLMNet

- Paper: [BioLLMNet: a multimodal framework for RNA-centric interaction prediction using large language model embeddings](https://doi.org/10.1093/bib/bbaf549)
- Source code: [abrarrahmanabir/BioLLMNet](https://github.com/abrarrahmanabir/BioLLMNet)
- Data: [preprocessed data provided by the authors](https://drive.google.com/drive/folders/1qDX5u_5BgptB0Ah5o4-uEF91oSV4-7ci?usp=sharing)
- Main entry point: `rna_molecule_run.py`
- Evaluation entry point: `rna_molecule_eval.py`

The author-provided data already contain RiNALMo RNA representations, MoleBERT small-molecule representations, and labels. The training script does not regenerate embeddings from raw RNA sequences or SMILES strings. The RNA–small-molecule CSV file must contain `Compound`, `Protein`, and `Label` columns. RNA-embedding and drug-embedding pickle files are also required.

### Run

```bash
git clone https://github.com/abrarrahmanabir/BioLLMNet.git
cd BioLLMNet

conda create -n biollmnet python=3.10 -y
conda activate biollmnet
pip install torch pandas numpy scikit-learn

# Download the RNA_Molecule data from the Google Drive link above and extract it.
# At the end of rna_molecule_run.py, edit train_files and explicitly specify:
# 1. the RNA-embedding pickle file
# 2. the drug-embedding pickle file
# 3. the interaction CSV file
python rna_molecule_run.py
```

The script performs 10-fold cross-validation and writes results and weights under `final_results/`. For evaluation, edit `CODE_DIR` in `rna_molecule_eval.py` so that the validation data and model filenames match the local outputs:

```bash
python rna_molecule_eval.py
```

Before running, replace the absolute Windows paths from the authors' machine at the end of `rna_molecule_run.py` with local paths. Specify the three input files explicitly rather than relying on directory ordering. Before running the evaluation script, also change `CODE_DIR` in `rna_molecule_eval.py` to the actual result directory.

## 2. DeepRSMA

- Paper: [DeepRSMA: a deep learning model for RNA–small molecule activity prediction](https://doi.org/10.1093/bioinformatics/btae678)
- Source code: [Hhhzj-7/DeepRSMA](https://github.com/Hhhzj-7/DeepRSMA)
- Data: `data/` in the repository
- Cross-validation entry point: `main_cv.py`
- Blind-test entry point: `main_blind.py`
- Independent-test entry point: `main_independent.py`

The repository includes raw data, precomputed RNA representations, RNA contact information, and preprocessing scripts. Generating RNA features from scratch additionally requires [RNA-FM](https://github.com/ml4bio/RNA-FM) and [SPOT-RNA-2D](https://github.com/jaswindersingh2/SPOT-RNA-2D). To reproduce the reported experiments, the representations and contact data already provided in the repository can be used directly.

### Run

```bash
git clone https://github.com/Hhhzj-7/DeepRSMA.git
cd DeepRSMA
conda env create -f environment.yml

# Read the name on the first line of environment.yml and activate that environment.
conda activate <environment-name>

python main_cv.py
python main_blind.py
python main_independent.py
```

Save outputs from the three tasks separately. Also record the random seed, GPU model, CUDA and PyTorch versions, and commit SHA. The experiments in the paper use multiple GPUs. Training on a single GPU can take substantially longer.

## 3. RSAPred

- Paper: [RSAPred: a machine learning approach for predicting RNA–small molecule binding affinity](https://doi.org/10.1093/bib/bbae002)
- Source code: [Sowmya-R-Krishnan/RSAPred](https://github.com/Sowmya-R-Krishnan/RSAPred)
- Web server and data download: [RSAPred prediction page](https://web.iitm.ac.in/bioinfo2/RSAPred/Predict.html)
- Original RNA–small-molecule interaction database: [R-SIM](https://web.iitm.ac.in/bioinfo2/R_SIM/)

RSAPred combines conventional descriptors with linear-regression models. The repository contains sample data and scripts for each processing stage. Complete training and test data should be obtained from the authors' website or R-SIM. The repository README requires Python 3.8 or later. Major dependencies include pandas, NumPy, SciPy, scikit-learn, sklearn-genetic, Mordred, Open Babel, and ViennaRNA.

### Install

```bash
git clone https://github.com/Sowmya-R-Krishnan/RSAPred.git
cd RSAPred
conda create -n rsapred python=3.8 -y
conda activate rsapred
pip install -r requirements.txt
```

Open Babel and ViennaRNA are generally easier to install from conda-forge. If `pip install -r requirements.txt` fails on either package, use:

```bash
conda install -c conda-forge openbabel viennarna -y
pip install -r requirements.txt
```

### Data preprocessing

Run the following commands from the corresponding data-preprocessing directories in the repository. Confirm the directory names against the README in the cloned version.

```bash
python calc_kmer_composition_features_v1.py "./data/rna_data.csv" "./sample_output/"
python calc_ppseudoDNC_features_v1.py "./data/rna_data.csv" "./data/Physicochemical_indices_RNA.csv" "./sample_output/"
python calc_triplet_str_composition_features_v1.py "./data/rna_data.csv" "./sample_output/"
python calc_pseudo_structure_composition_v1.py "./data/rna_data.csv" "./sample_output/"

obabel -ismi ./data/mol_data.smi -osdf -O ./sample_output/mol_data.sdf --gen3d -h
python -m mordred ./sample_output/mol_data.sdf -t sdf -o ./sample_output/Mol_features_v1.csv -3

python combine_RNA_features_v1.py \
  "./sample_output/Mononucleotide_v1.out" \
  "./sample_output/Dinucleotide_v1.out" \
  "./sample_output/Trinucleotide_v1.out" \
  "./sample_output/Tetranucleotide_v1.out" \
  "./sample_output/pPseudoDNC_features_v1.out" \
  "./sample_output/Triplet_features_v1.out" \
  "./sample_output/Pseudo_structure_status_composition_v1.out" \
  "./sample_output/RNA_features_v1.csv"

python create_dataset_v1.py \
  "./sample_output/RNA_features_v1.csv" \
  "./sample_output/Mol_features_v1.csv" \
  "./data/sample_data.csv" \
  "./sample_output/Final_sample_dataset_v1.csv"
```

### Feature selection, cross-validation, and external testing

```bash
# <nfeat> is the target number of features, and <output_path> is the output directory.
python -u ffs_final.py \
  "./data/Final_sample_dataset_v1.csv" \
  "./data/Pairwise_featcorr_0.8.pkl" \
  3 <nfeat> "all" <output_path>

python -u perform_10_fold_CV.py \
  "./data/Final_sample_dataset_v1.csv" \
  "best_models.log" 3 <output_path>

python test_on_regression_dataset.py \
  "./data/QSAR_dataset_final.csv" \
  "../Feature_selection/data/Final_sample_dataset_v1.csv" \
  "best_models.log" \
  "test_results.csv"
```

README files in the repository's subdirectories also provide commands for RFECV, genetic algorithms, leave-one-out cross-validation, and external classification tests. Record the RNA subtype, number of features, and best-model log used for each run. The repository README also states that the software is for academic use only. Contact the authors before commercial use or redistribution, even if the repository contains a separate license file.

## 4. RLaffinity

- Paper: [RLaffinity: a deep learning model for RNA–ligand binding affinity prediction](https://doi.org/10.1093/bioinformatics/btae155)
- Source code: [SaisaiSun/RLaffinity](https://github.com/SaisaiSun/RLaffinity)
- Main directory: `3dcnn_lba/`
- Data: nucleic-acid–ligand complexes from PDBbind NL2020; the repository provides cleaned labels, training/validation/test lists, and selected model outputs

RLaffinity requires three-dimensional receptor and ligand structures. The repository contains `data/train_list.txt`, `data/val_list.txt`, `data/test_list.txt`, `data/input_label/pdbbind_NL_cleaned.csv`, and trained weights. The original PDBbind structures must still be obtained under the PDBbind access terms.

### Run

```bash
git clone https://github.com/SaisaiSun/RLaffinity.git
cd RLaffinity/3dcnn_lba

conda create -n rlaffinity python=3.9 -y
conda activate rlaffinity
pip install torch numpy pandas scipy tqdm biopython

# Use these commands to verify argument names and input formats in the current commit.
python process_pdbbind.py --help
python prepare_lmdb.py --help
python train.py --help

# 1. Preprocess the three-dimensional structures.
python process_pdbbind.py <receptor_dir> <ligand_dir> --out_dir <processed_dir>

# 2. Generate the LMDB used for contrastive pretraining.
python prepare_lmdb.py <processed_dir> <pretrain_lmdb>
python trainstage1.py

# 3. Generate the supervised dataset.
python prepare_lmdb.py <processed_dir> <supervised_lmdb> -s \
  --train_txt data/train_list.txt \
  --val_txt data/val_list.txt \
  --test_txt data/test_list.txt \
  --score_path data/input_label/pdbbind_NL_cleaned.csv

# 4. Train or evaluate the model.
python train.py --data_dir data --mode test --output_dir output_train
```

For a new complex, process the structures and generate a test LMDB with the same procedure, then run:

```bash
python test.py --data_dir data --output_dir output_test
```

The original PDBbind structures must be obtained under the applicable access terms. A GitHub repository should document the download and preprocessing procedure rather than redistribute restricted structures. The official repository does not pin dependency versions. Export the working environment immediately after the first successful run:

```bash
conda env export --no-builds > environment.yml
```

## 5. RLASIF

- Paper: [RLASIF: RNA–ligand affinity prediction based on surface interaction fingerprints](https://doi.org/10.1016/j.compbiolchem.2025.108367)
- Currently discoverable repository: [ZUSTSTTLAB/RLASIF](https://github.com/ZUSTSTTLAB/RLASIF)
- Data: RNA–ligand structures derived from PDBbind NL2020; the comparison in the review uses an affinity subset of 95 complexes

### Why the method cannot currently be run

The core `RLASIF` entry in the public repository is only a Git submodule pointer, but the repository does not provide the submodule URL in `.gitmodules`. The actual model source code therefore cannot be retrieved. The remaining repository content consists mainly of macOS metadata and does not include a usable README, environment file, preprocessing scripts, training entry point, data split, or pretrained weights. A verifiable execution command cannot be constructed without a complete repository or a source archive tied to a fixed commit from the authors.

## 6. SPRank

- Paper: [SPRank: a knowledge-based scoring function for RNA–ligand complexes](https://doi.org/10.1021/acs.jctc.4c00681)
- Full text and supporting information: [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12960052/)
- Source-code URL reported in the paper: `https://github.com/Vfold-RNA/SPRank`

### Why the method cannot currently be run

The reported source-code URL, `https://github.com/Vfold-RNA/SPRank`, returned 404 on the verification date. The standalone SPRank program, feature and statistical-potential files, random-forest weights, dependency versions, and command-line instructions are therefore unavailable. The supporting information identifies the training and test sets but cannot replace the missing program and model files. Although the paper describes generating candidate poses with rDock or AutoDock Vina and scoring them with SPRank, that description alone is insufficient to construct a verifiable command. The authors would need to restore the repository or provide a source archive.

## 7. RNAmigos

- Paper: [Augmented base pairing networks encode RNA-small molecule binding preferences](https://doi.org/10.1093/nar/gkaa583)
- Source code: [cgoliver/RNAmigos](https://github.com/cgoliver/RNAmigos)
- Paper data: [Zenodo 8338267](https://zenodo.org/records/8338267)
- Training entry point: `learning/main.py`
- Custom-structure inference entry point: `inference.py`

RNAmigos represents a known RNA binding pocket as a graph with canonical and noncanonical base-pair types and predicts the MACCS fingerprint of a candidate ligand. An input structure must contain only previously identified pocket residues. The program does not locate the binding site itself.

### Environment and paper data

```bash
git clone https://github.com/cgoliver/RNAmigos.git
cd RNAmigos

conda env create -f environment.yml
conda activate rnamigos_minimal

cd data
tar -xzvf pockets_nx_symmetric_orig.tar.gz
mkdir -p annotated
mv pockets_nx_symmetric_orig annotated/
cd ..
```

The repository's `environment.yml` pins Python 3.6, PyTorch 1.5.1, and DGL 0.4.3. Use the authors' environment first because the old code is incompatible with parts of the current DGL and PyTorch APIs. Zenodo provides the cleaned training and validation data and decoy sets used in the paper. Running `make_nice.py` from that archive produces `rnamigos1_dataset.csv`.

### Train

```bash
python learning/main.py \
  -da pockets_nx_symmetric_orig \
  -n rnamigos_reproduction

# Display all training arguments.
python learning/main.py -h
```

Models and logs are saved in the run directory specified with `-n`.

### Infer on a custom RNA pocket

```bash
mkdir -p data/my_pdbs data/my_graphs

# Place a .cif file containing only the target binding-pocket residues in data/my_pdbs/.
# Follow the inference.py example to generate a graph with rnaglib's fr3d_to_graph and run the model.
python inference.py
```

The output is a vector of predicted MACCS-fingerprint probabilities. To screen a small-molecule library, calculate MACCS fingerprints for the candidate molecules separately and rank them by similarity as illustrated in the repository.

## 8. R-BIND

- Paper: [R-BIND: An Interactive Database for Exploring and Developing RNA-Targeted Chemical Probes](https://doi.org/10.1021/acschembio.9b00631)
- Updated paper: [R-BIND 2.0](https://doi.org/10.1021/acschembio.2c00224)
- Database and online search: [R-BIND](https://rbind.chem.duke.edu/)
- R-BIND 2.0 data tables: [full paper and Supporting Information](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015/)

### Why the UNK96 comparison in the table cannot currently be run

R-BIND is a database and an online cheminformatics search platform rather than a standalone prediction package with a public command-line entry point. The R-BIND paper describes the database content, 20 chemical descriptors, and the online nearest-neighbor search. The Supporting Information for R-BIND 2.0 provides `RBIND_v2.0_A.xlsx` and `RBIND_v2.0_B.xlsx`, but the website back-end source code, a fixed-version feature-calculation script, and an environment specification are not public.

The UNK96 table in the review uses a ligand-ranking workflow specific to that test set. The public materials do not provide the UNK96 input files, preprocessing map, complete candidate library, or executable scoring script used for that comparison. The reported values therefore cannot be reconstructed directly from the R-BIND website or data tables. Rewriting a similarity-ranking algorithm from the paper would be a reimplementation, not execution of the authors' code.

## 9. RNALigands

- Paper: [RNALigands: a database and web server for RNA-ligand interactions](https://doi.org/10.1261/rna.078889.121)
- Source code: [SaisaiSun/RNALigands](https://github.com/SaisaiSun/RNALigands)
- Data: motif–ligand data files, substitution matrices, and examples under `Package/` in the repository
- Entry point: `Package/run.pl`

RNALigands extracts hairpin, internal, bulge, and multibranch-loop motifs from an RNA sequence or dot-bracket secondary structure. It then searches for similar motifs and associated ligands in data derived from PDB, R-BIND, and miRBase. The code uses Perl, the ViennaRNA command-line tools, and Unix file commands and should be run on Linux or WSL.

### Environment and execution

```bash
git clone https://github.com/SaisaiSun/RNALigands.git
cd RNALigands/Package

conda create -n rnaligands -c conda-forge perl viennarna -y
conda activate rnaligands
chmod +x *.pl

# First replace the hard-coded path in run.pl:
# /var/www/rnaligands/ViennaRNA/bin/RNAfold
# with the RNAfold executable in the current environment.
which RNAfold

# Use a FASTA file and let RNAfold generate the secondary structure before motif searching.
perl run.pl -f example/1ddy_A.fasta

# If a dot-bracket secondary structure is already available, use -s directly.
perl run.pl -s example/1ddy_A_dot.txt
```

Run the command from the `Package/` directory because `run.pl` locates the other Perl scripts and database files relative to the current working directory. Output is written to the directory containing the input example. Direct execution in Windows PowerShell fails because the script depends on the Unix `cp` command, Unix path behavior, and executable permissions. Use WSL or Linux.

## 10. ZHMol-RLinter

- Paper: [A Machine Learning Method for RNA-Small Molecule Binding Preference Prediction](https://doi.org/10.1021/acs.jcim.4c01324)
- Data from the original paper: [ACS Supporting Information](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01324)
- Inference code and models released later: [Zenodo 17157778](https://doi.org/10.5281/zenodo.17157778)
- Zenodo file: `ZHMol-RLinter_on_TAR.tar`

The Supporting Information for the original paper provides data tables and experimental results for the RNA–small-molecule database, RL98, UNK96, and PC40, but no complete training repository was released. In a later study, the authors published TAR-inhibitor inference examples, random-forest `.mat` models, and selected feature-generation scripts for ZHMol-RLinter.

### Run the public TAR inference example

```bash
curl -L \
  https://zenodo.org/api/records/17157778/files/ZHMol-RLinter_on_TAR.tar/content \
  -o ZHMol-RLinter_on_TAR.tar

tar -xf ZHMol-RLinter_on_TAR.tar
cd ZHMol-RLinter_on_TAR/example_TAR_110FA

matlab -batch "test_TAR"
```

Other prepared example directories include `example_TAR_115FA`, `example_TAR_AM6538`, `example_TAR_DB00594`, and `example_TAR_F07#13`. Each directory contains a feature table, a random-forest model, and `test_TAR.m`. Predictions are written to `predict_result/scores_motif.xlsx` in the corresponding directory. Running these examples requires MATLAB and the Statistics and Machine Learning Toolbox.

### Prepare features for a new input

The workflow described in the archive README is:

1. Predict RNA secondary structure with [MXfold2](https://github.com/mxfold/mxfold2) and extract loop motifs from a PDB structure.
2. Calculate the Laplacian norm with `feature preparation/Laplacian Norm calculation/ln.pl`.
3. Calculate the physicochemical environment with `feature preparation/Physicochemical environment/PE_feature.py`. This script requires NumPy and Open Babel, and its absolute input and output paths must be edited before use.
4. Run `feature preparation/Network topology/Network_T.m` in MATLAB.
5. Extract the loop-motif pocket with [GHECOM](https://pdbj.org/ghecom/).
6. Generate a MACCS fingerprint with `feature preparation/small molecule fingerprint/fingerprint.py`. This script requires RDKit, and its input and output paths must be edited before use.
7. Merge the motif and ligand features into the 188-dimensional `feature.xlsx` format shown in the examples, then run the corresponding MATLAB test script.

### Why the original training and UNK96 test cannot be fully reconstructed

The Zenodo archive contains only TAR inference examples and feature-preparation scripts. It does not include the random-forest training program or the complete `test_program/`, `test_UNK96_1.m`, and input directories referenced in the README. The Supporting Information for the original paper provides data tables and results rather than the missing programs. The TAR examples can therefore be run, but the currently public files cannot train the model from scratch or strictly reconstruct the RL98 and UNK96 experiments reported in the paper.
