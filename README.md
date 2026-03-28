# Machine Learning for RNA–Small Molecule Interactions: Methodologies, Evaluation, and Future Direction

📄 **Paper**: *Machine Learning for RNA–Small Molecule Interactions: Methodologies, Evaluation, and Future Directions*

---

## 🔬 Overview

RNA–small molecule interactions play a central role in drug discovery, especially for targeting non-coding RNAs and structured RNA elements.  

This repository accompanies our survey and benchmark study, providing:

- A **systematic review** of machine learning and deep learning methods
- A **task-oriented taxonomy** of RNA–ligand modeling
- A **critical analysis of widely used datasets (Test9, Test18)**
- A **clean benchmark (Test17)** for reliable evaluation

---

## 🧠 Task Taxonomy

We categorize RNA–small molecule modeling into three core tasks:

| Task | Description |
|------|------------|
| Binding Site Prediction | Identify ligand-binding residues |
| Binding Preference Prediction | Predict binding affinity / interaction strength |
| Binding Pose Prediction | Predict 3D binding conformation |

---

# 📚 Literature Overview

# Binding Site Prediction

## Statistical Methods

* [Rsite](https://scholar.google.com/scholar?q=Rsite+RNA+binding+site+zeng+2015)
* [Rsite2](https://scholar.google.com/scholar?q=Rsite2+RNA+binding+site+zeng+2016)
* [RBind](https://scholar.google.com/scholar?q=RBind+RNA+binding+site+wang+2018)
* [MetalionRNA](https://scholar.google.com/scholar?q=MetalionRNA+philips+2012)

## Traditional Machine Learning

* [RNAsite](https://scholar.google.com/scholar?q=RNAsite+su+2021+RNA+binding+site)
* [RNet](https://scholar.google.com/scholar?q=RNet+RNA+binding+site+liu+2024)
* [NABS](https://scholar.google.com/scholar?q=NABS+jiang+2022+nucleic+acid+binding+site)

## Deep Learning

* [RLBind](https://scholar.google.com/scholar?q=RLBind+RNA+binding+site+wang+2023)
* [MultiModRLBP](https://scholar.google.com/scholar?q=MultiModRLBP+wang+2024)
* [ZeSTa](https://scholar.google.com/scholar?q=ZeSTa+RNA+binding+site+gao+2024)
* [BiteNet_N](https://scholar.google.com/scholar?q=BiteNet+RNA+binding+site+kozlovskii+2021)
* [RLBSIF](https://schalar.google.com/scholar?q=RLBSIF+RNA+binding+site+sang+2025)
* [RLsite](https://scholar.google.com/scholar?q=RLsite+RNA+binding+site+zou+2025)

## Graph & Geometric Deep Learning / Foundation Model-based

* [CoBRA](https://scholar.google.com/scholar?q=CoBRA+RNA+binding+site+jang+2026)
* [GATRsite](https://scholar.google.com/scholar?q=GATRsite+sun+2025+RNA)
* [RLsite](https://scholar.google.com/scholar?q=RLsite+ERNIE-RNA+sun+2025)
* [RNABind](https://scholar.google.com/scholar?q=RNABind+zhu+2025)
* [MVRBind](https://scholar.google.com/scholar?q=MVRBind+chen+2025)

---

# Binding Preference Prediction

## Classification-based Methods

* [ZHMol-RLinter](https://scholar.google.com/scholar?q=ZHMol-RLinter+zhuo+2024)
* [DistRMI](https://scholar.google.com/scholar?q=DistRMI+liu+2025)

## Affinity-based Regression

* [RSAPred](https://scholar.google.com/scholar?q=RSAPred+krishnan+2024)
* [SPrank](https://scholar.google.com/scholar?q=SPrank+zhou+2024)
* [DeepRSMA](https://scholar.google.com/scholar?q=DeepRSMA+huang+2024)
* [RLaffinity](https://scholar.google.com/scholar?q=RLaffinity+sun+2024)
* [RLASIF](https://scholar.google.com/scholar?q=RLASIF+xia+2025)
* [EMMPTNet](https://scholar.google.com/scholar?q=EMMPTNet+wu+2025)
* [DLRNA-BERTa](https://scholar.google.com/scholar?q=DLRNA-BERTa+lobascio+2025)
* [BioLLMNet](https://scholar.google.com/scholar?q=BioLLMNet+abir+2025)

## Score-based Regression

* [RNAmigos](https://scholar.google.com/scholar?q=RNAmigos+oliver+2020)
* [RNAmigos2](https://scholar.google.com/scholar?q=RNAmigos2+carvajal+2025)
* [RNAsmol](https://scholar.google.com/scholar?q=RNAsmol+ma+2025)
* [SMRTnet](https://scholar.google.com/scholar?q=SMRTnet+fei+2026)
* [DeepRNA-DTI](https://scholar.google.com/scholar?q=DeepRNA-DTI+bae+2025)
* [GerNA-Bind](https://scholar.google.com/scholar?q=GerNA-Bind+xia+2025)

---

# Binding Pose Prediction

## Scoring / Rescoring Methods

* [LigandRNA](https://scholar.google.com/scholar?q=LigandRNA+philips+2013)
* [RNAPosers](https://scholar.google.com/scholar?q=RNAPosers+chhabra+2020)
* [AnnapuRNA](https://scholar.google.com/scholar?q=AnnapuRNA+stefaniak+2021)
* [RmsdXNA](https://scholar.google.com/scholar?q=RmsdXNA+tan+2024)

## Docking Methods

* [AutoDock](https://autodock.scripps.edu/)
* [Dock6](https://dock.compbio.ucsf.edu/)
* [rDock](https://rdock.github.io/)
* [RLDOCK](https://scholar.google.com/scholar?q=RLDOCK+sun+2020)
* [NLDock](https://scholar.google.com/scholar?q=NLDock+feng+2021)

## Pocket Detection Methods

* [fPocket](https://github.com/Discngine/fpocket)
* [fPocketR](https://scholar.google.com/scholar?q=fPocketR+veenbaas+2025)

---


# 📊 Dataset Analysis

## ⚠️ Existing Problems

### Test18

- Contains **non-real ligand case (6EZ0)**
- Modified nucleotide treated as ligand

---

### Test9

- Severe redundancy with training data:

| Metric | Value |
|--------|------|
| Sequence identity | up to 100% |
| TM-score | up to 0.987 |

---

## ✅ Our Solution: Test17

We construct a **clean benchmark dataset**:

- Remove redundant samples  
- Remove erroneous annotations  
- Ensure strict separation from training data  

---

# 🧾 Key Insights

* Dataset quality critically affects model evaluation
* Many reported results are **overestimated due to leakage**
* Proper data splitting is essential for generalization

---

# 🚀 Future Directions

* RNA-specific foundation models
* Interpretable deep learning
* Data-centric AI for RNA
* RNA-targeted drug generation

---

# 📖 Citation

If you find this work useful, please consider citing:

```bibtex
```

---


# 📬 Contact

For questions about the paper, please contact:

```
leideng@csu.edu.cn
wumin@a-star.edu.sg
```
