# RNA–小分子 Binding Preference 方法复现说明

本文档整理 RNA–小分子 binding preference / affinity 方法，包括 **BioLLMNet、DeepRSMA、RSAPred、RLaffinity、RLASIF、SPRank、RNAmigos、R-BIND、RNALigands 和 ZHMol-RLinter**。内容包括论文与源代码入口、数据来源和运行步骤；对当前缺少可执行源码的方法，仅说明阻塞原因。

> 核对日期：2026-09-17。公开仓库和下载链接可能后续发生变化。

## 1. BioLLMNet

- 论文：[BioLLMNet: a multimodal framework for RNA-centric interaction prediction using large language model embeddings](https://doi.org/10.1093/bib/bbaf549)
- 源码：[abrarrahmanabir/BioLLMNet](https://github.com/abrarrahmanabir/BioLLMNet)
- 数据：[作者提供的预处理数据](https://drive.google.com/drive/folders/1qDX5u_5BgptB0Ah5o4-uEF91oSV4-7ci?usp=sharing)
- 主要入口：`rna_molecule_run.py`
- 评估入口：`rna_molecule_eval.py`

作者提供的数据已经包含 RiNALMo RNA 表示、MoleBERT 小分子表示和标签，训练脚本本身不负责从原始 RNA/SMILES 重新生成 embedding。RNA–小分子 CSV 需要包含 `Compound`、`Protein` 和 `Label` 列，同时需要 RNA embedding 与 drug embedding 的 pickle 文件。

### 运行

```bash
git clone https://github.com/abrarrahmanabir/BioLLMNet.git
cd BioLLMNet

conda create -n biollmnet python=3.10 -y
conda activate biollmnet
pip install torch pandas numpy scikit-learn

# 从上面的 Google Drive 下载 RNA_Molecule 数据并解压。
# 修改 rna_molecule_run.py 底部的 train_files，明确填写：
# 1. RNA embedding pickle
# 2. drug embedding pickle
# 3. interaction CSV
python rna_molecule_run.py
```

脚本执行 10 折交叉验证，并在 `final_results/` 下写出结果和权重。评估时还需要修改 `rna_molecule_eval.py` 中的 `CODE_DIR`，使验证数据和模型文件名与本地输出一致：

```bash
python rna_molecule_eval.py
```

运行前必须把 `rna_molecule_run.py` 末尾作者机器上的 Windows 绝对路径替换成本地路径。三个输入文件应逐一明确指定，不要依赖目录中文件的排序。运行评估脚本前，还要把 `rna_molecule_eval.py` 中的 `CODE_DIR` 改为实际结果目录。

## 2. DeepRSMA

- 论文：[DeepRSMA: a deep learning model for RNA–small molecule activity prediction](https://doi.org/10.1093/bioinformatics/btae678)
- 源码：[Hhhzj-7/DeepRSMA](https://github.com/Hhhzj-7/DeepRSMA)
- 数据：仓库中的 `data/`
- 交叉验证入口：`main_cv.py`
- 盲测入口：`main_blind.py`
- 独立测试入口：`main_independent.py`

仓库已包含原始数据、预计算 RNA 表示、RNA contact 信息和预处理脚本。若要从头生成 RNA 特征，需要额外使用 [RNA-FM](https://github.com/ml4bio/RNA-FM) 和 [SPOT-RNA-2D](https://github.com/jaswindersingh2/SPOT-RNA-2D)。仅复现实验时，可优先使用仓库中已提供的表示和 contact 数据。

### 运行

```bash
git clone https://github.com/Hhhzj-7/DeepRSMA.git
cd DeepRSMA
conda env create -f environment.yml

# 查看 environment.yml 第一行的 name，并激活该环境。
conda activate <environment-name>

python main_cv.py
python main_blind.py
python main_independent.py
```

建议把三种任务的输出分别保存，并额外记录随机种子、GPU 型号、CUDA/PyTorch 版本和 commit SHA。论文实验使用多 GPU；若改成单 GPU，需要预期训练时间明显增加。

## 3. RSAPred

- 论文：[RSAPred: a machine learning approach for predicting RNA–small molecule binding affinity](https://doi.org/10.1093/bib/bbae002)
- 源码：[Sowmya-R-Krishnan/RSAPred](https://github.com/Sowmya-R-Krishnan/RSAPred)
- Web 服务及数据下载：[RSAPred prediction page](https://web.iitm.ac.in/bioinfo2/RSAPred/Predict.html)
- 原始 RNA–小分子相互作用数据库：[R-SIM](https://web.iitm.ac.in/bioinfo2/R_SIM/)

RSAPred 是传统描述符与线性回归方法。仓库包含样例数据和各阶段脚本；完整训练/测试数据应从作者网站或 R-SIM 获取。README 要求 Python 3.8+，主要依赖包括 pandas、NumPy、SciPy、scikit-learn、sklearn-genetic、Mordred、Open Babel 和 ViennaRNA。

### 安装

```bash
git clone https://github.com/Sowmya-R-Krishnan/RSAPred.git
cd RSAPred
conda create -n rsapred python=3.8 -y
conda activate rsapred
pip install -r requirements.txt
```

Open Babel 和 ViennaRNA 通常更适合通过 conda-forge 安装；若 `pip install -r requirements.txt` 在这两个包处失败，可使用：

```bash
conda install -c conda-forge openbabel viennarna -y
pip install -r requirements.txt
```

### 数据预处理

以下命令应在仓库对应的 data preprocessing 目录执行；路径名称以克隆后的实际 README 为准。

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

### 特征选择、交叉验证和外部测试

```bash
# <nfeat> 为目标特征数，<output_path> 为输出目录。
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

仓库各子目录的 README 还提供 RFECV、遗传算法、LOO-CV 和分类外部测试命令。复现时应严格记录所使用的 RNA subtype、特征数和最佳模型日志。仓库 README 同时写有“仅限学术使用”的说明；商业使用或再次分发前应联系作者确认，即使仓库中另有许可证文件也不要忽略该限制。

## 4. RLaffinity

- 论文：[RLaffinity: a deep learning model for RNA–ligand binding affinity prediction](https://doi.org/10.1093/bioinformatics/btae155)
- 源码：[SaisaiSun/RLaffinity](https://github.com/SaisaiSun/RLaffinity)
- 主要目录：`3dcnn_lba/`
- 数据：PDBbind NL2020 中的核酸–配体复合物；仓库提供清洗标签、训练/验证/测试列表和部分模型输出

RLaffinity 需要受体和配体的三维结构。仓库中可见 `data/train_list.txt`、`data/val_list.txt`、`data/test_list.txt`、`data/input_label/pdbbind_NL_cleaned.csv` 以及训练得到的权重，但原始 PDBbind 结构仍需按 PDBbind 的授权方式自行取得。

### 运行

```bash
git clone https://github.com/SaisaiSun/RLaffinity.git
cd RLaffinity/3dcnn_lba

conda create -n rlaffinity python=3.9 -y
conda activate rlaffinity
pip install torch numpy pandas scipy tqdm biopython

# 运行以下命令核对当前 commit 的参数名和输入格式。
python process_pdbbind.py --help
python prepare_lmdb.py --help
python train.py --help

# 1. 预处理三维结构。
python process_pdbbind.py <receptor_dir> <ligand_dir> --out_dir <processed_dir>

# 2. 生成对比预训练使用的 LMDB。
python prepare_lmdb.py <processed_dir> <pretrain_lmdb>
python trainstage1.py

# 3. 生成有监督数据。
python prepare_lmdb.py <processed_dir> <supervised_lmdb> -s \
  --train_txt data/train_list.txt \
  --val_txt data/val_list.txt \
  --test_txt data/test_list.txt \
  --score_path data/input_label/pdbbind_NL_cleaned.csv

# 4. 训练/评估。
python train.py --data_dir data --mode test --output_dir output_train
```

对新复合物预测时，先以相同步骤处理结构并生成测试 LMDB，再执行：

```bash
python test.py --data_dir data --output_dir output_test
```

PDBbind 原始结构需要按其授权方式自行取得。GitHub 仓库中只记录下载与预处理步骤，不应直接提交受限的原始结构。由于官方仓库没有锁定依赖版本，首次成功运行后应立即导出实际环境：

```bash
conda env export --no-builds > environment.yml
```

## 5. RLASIF

- 论文：[RLASIF: RNA–ligand affinity prediction based on surface interaction fingerprints](https://doi.org/10.1016/j.compbiolchem.2025.108367)
- 当前可找到的仓库：[ZUSTSTTLAB/RLASIF](https://github.com/ZUSTSTTLAB/RLASIF)
- 数据：论文使用 PDBbind NL2020 派生的 RNA–配体结构数据；综述的统一比较使用 95 个复合物的 affinity 子集

### 无法运行的原因

当前公开仓库中的核心 `RLASIF` 条目只是 Git 子模块指针，但仓库没有提供 `.gitmodules` 中的子模块 URL，因此无法取得真正的模型源码。仓库其余内容主要是 macOS 元数据，也没有可用的 README、环境文件、数据预处理脚本、训练入口、数据划分或预训练权重。缺少这些文件时无法构造可验证的运行命令，需要作者补充完整仓库或固定 commit 的源码归档。

## 6. SPRank

- 论文：[SPRank: a knowledge-based scoring function for RNA–ligand complexes](https://doi.org/10.1021/acs.jctc.4c00681)
- 论文全文与补充材料：[PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12960052/)
- 论文声明的源码地址：`https://github.com/Vfold-RNA/SPRank`

### 无法运行的原因

论文声明的源码地址 `https://github.com/Vfold-RNA/SPRank` 在核对日期返回 404，因而无法取得 SPRank standalone 程序、特征与统计势文件、随机森林权重、依赖版本或命令行说明。论文补充材料可以确认训练集和测试集组成，但不能代替缺失的程序与模型文件。虽然论文描述了使用 rDock/AutoDock Vina 生成候选构象并用 SPRank 评分的流程，但仅凭方法描述无法写出可验证的运行命令，需要作者恢复仓库或提供源码归档。

## 7. RNAmigos

- 论文：[Augmented base pairing networks encode RNA-small molecule binding preferences](https://doi.org/10.1093/nar/gkaa583)
- 源码：[cgoliver/RNAmigos](https://github.com/cgoliver/RNAmigos)
- 论文数据：[Zenodo 8338267](https://zenodo.org/records/8338267)
- 训练入口：`learning/main.py`
- 自定义结构推理入口：`inference.py`

RNAmigos 将已知 RNA binding pocket 表示成带有 canonical 和 non-canonical base-pair 类型的图，并预测候选配体的 MACCS fingerprint。输入结构必须只包含已经确定的 pocket residues；该程序本身不负责寻找 binding site。

### 环境与论文数据

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

仓库中的 `environment.yml` 固定了旧版 Python 3.6、PyTorch 1.5.1 和 DGL 0.4.3。为避免旧代码与新版 DGL/PyTorch API 不兼容，应先使用作者环境，不要直接升级核心依赖。Zenodo 提供清理后的论文训练/验证数据和 decoy 集合；运行其中的 `make_nice.py` 可以生成 `rnamigos1_dataset.csv`。

### 训练

```bash
python learning/main.py \
  -da pockets_nx_symmetric_orig \
  -n rnamigos_reproduction

# 查看全部训练参数。
python learning/main.py -h
```

模型和日志保存在 `-n` 指定的运行目录中。

### 对自定义 RNA pocket 推理

```bash
mkdir -p data/my_pdbs data/my_graphs

# 将仅包含目标 binding pocket residues 的 .cif 文件放入 data/my_pdbs/。
# 按 inference.py 中的示例用 rnaglib 的 fr3d_to_graph 生成图并执行模型。
python inference.py
```

输出是预测的 MACCS fingerprint 概率。使用候选小分子库进行筛选时，需要另外计算候选分子的 MACCS fingerprint，再按照仓库示例进行相似度排序。

## 8. R-BIND

- 论文：[R-BIND: An Interactive Database for Exploring and Developing RNA-Targeted Chemical Probes](https://doi.org/10.1021/acschembio.9b00631)
- 更新版论文：[R-BIND 2.0](https://doi.org/10.1021/acschembio.2c00224)
- 数据库与在线检索：[R-BIND](https://rbind.chem.duke.edu/)
- R-BIND 2.0 数据表：[论文全文和 Supporting Information](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015/)

### 无法运行表中 UNK96 比较的原因

R-BIND 本身是数据库和在线 cheminformatics 检索平台，不是带有公开命令行入口的独立预测软件。R-BIND 论文公开了数据库内容、20 个化学描述符以及在线 nearest-neighbor search 的原理，R-BIND 2.0 的 Supporting Information 也提供 `RBIND_v2.0_A.xlsx` 和 `RBIND_v2.0_B.xlsx`，但没有公开网站后端源码、固定版本的特征计算脚本或环境文件。

图中 UNK96 表格使用的是针对该测试集的 ligand ranking 流程。当前公开材料没有提供该比较所用的 UNK96 输入文件、预处理映射、完整候选库和可执行评分脚本，因此不能从 R-BIND 网站或数据表直接重建表中的数值。仅根据论文描述重新编写相似度算法属于重新实现，而不是运行原作者代码。

## 9. RNALigands

- 论文：[RNALigands: a database and web server for RNA-ligand interactions](https://doi.org/10.1261/rna.078889.121)
- 源码：[SaisaiSun/RNALigands](https://github.com/SaisaiSun/RNALigands)
- 数据：仓库 `Package/` 中的 motif-ligand 数据文件、替换矩阵和示例
- 入口：`Package/run.pl`

RNALigands 从 RNA 序列或 dot-bracket 二级结构中提取 hairpin、internal、bulge 和 multibranch loop motif，然后在 PDB、R-BIND 和 miRBase 派生的数据中搜索相似 motif 及其配体。代码使用 Perl、ViennaRNA 命令行工具和 Unix 文件命令，建议在 Linux 或 WSL 中运行。

### 环境与运行

```bash
git clone https://github.com/SaisaiSun/RNALigands.git
cd RNALigands/Package

conda create -n rnaligands -c conda-forge perl viennarna -y
conda activate rnaligands
chmod +x *.pl

# 必须先把 run.pl 中硬编码的
# /var/www/rnaligands/ViennaRNA/bin/RNAfold
# 改为当前环境中的 RNAfold。
which RNAfold

# 使用 FASTA，让 RNAfold 生成二级结构后执行 motif 搜索。
perl run.pl -f example/1ddy_A.fasta

# 如果已经有 dot-bracket 二级结构，可直接使用 -s。
perl run.pl -s example/1ddy_A_dot.txt
```

必须从 `Package/` 目录执行命令，因为 `run.pl` 使用当前工作目录定位其余 Perl 脚本和数据库文件。输出写入输入示例所在目录。若直接在 Windows PowerShell 中运行，脚本中的 Unix `cp`、路径和可执行权限处理会失败，因此应使用 WSL/Linux。

## 10. ZHMol-RLinter

- 论文：[A Machine Learning Method for RNA-Small Molecule Binding Preference Prediction](https://doi.org/10.1021/acs.jcim.4c01324)
- 原论文数据：[ACS Supporting Information](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01324)
- 后续公开的推理代码与模型：[Zenodo 17157778](https://doi.org/10.5281/zenodo.17157778)
- Zenodo 文件：`ZHMol-RLinter_on_TAR.tar`

原论文的 Supporting Information 提供 RNA-small molecule database、RL98、UNK96 和 PC40 的数据表及实验结果，但没有发布完整训练仓库。作者在后续研究中公开了 ZHMol-RLinter 的 TAR inhibitor inference 示例、随机森林 `.mat` 模型和部分特征生成脚本。

### 运行公开的 TAR 推理示例

```bash
curl -L \
  https://zenodo.org/api/records/17157778/files/ZHMol-RLinter_on_TAR.tar/content \
  -o ZHMol-RLinter_on_TAR.tar

tar -xf ZHMol-RLinter_on_TAR.tar
cd ZHMol-RLinter_on_TAR/example_TAR_110FA

matlab -batch "test_TAR"
```

其他已准备好的示例目录包括 `example_TAR_115FA`、`example_TAR_AM6538`、`example_TAR_DB00594` 和 `example_TAR_F07#13`。每个目录已经包含 feature 表、随机森林模型和 `test_TAR.m`；预测结果写入该目录的 `predict_result/scores_motif.xlsx`。运行这些示例需要 MATLAB 及 Statistics and Machine Learning Toolbox。

### 为新输入准备特征

归档 README 给出的流程是：

1. 用 [MXfold2](https://github.com/mxfold/mxfold2) 预测 RNA 二级结构，并从 PDB 中提取 loop motif。
2. 用 `feature preparation/Laplacian Norm calculation/ln.pl` 计算 Laplacian norm。
3. 用 `feature preparation/Physicochemical environment/PE_feature.py` 计算 physicochemical environment；该脚本需要 NumPy 和 Open Babel，并且必须先修改其中的输入、输出绝对路径。
4. 在 MATLAB 中运行 `feature preparation/Network topology/Network_T.m`。
5. 用 [GHECOM](https://pdbj.org/ghecom/) 提取 loop motif pocket。
6. 用 `feature preparation/small molecule fingerprint/fingerprint.py` 生成 MACCS fingerprint；该脚本需要 RDKit，并且必须修改输入、输出路径。
7. 按示例把 motif 与 ligand 特征合并为 188 维 `feature.xlsx`，再运行相应的 MATLAB 测试脚本。

### 无法完整重建原论文训练和 UNK96 测试的原因

Zenodo 归档只包含 TAR 推理示例和特征准备脚本，没有随机森林训练程序，也没有 README 中提到的完整 `test_program/`、`test_UNK96_1.m` 及其输入目录。原论文 Supporting Information 提供的是数据表和结果，而不是这些缺失的程序。因此可以运行归档中的 TAR 示例，但无法用当前公开文件从头训练模型或严格重建论文中的 RL98/UNK96 实验。

## 11. 运行时统一记录的元数据

每个方法至少应记录：

- 原论文 DOI；
- 官方仓库 URL 和实际使用的 commit SHA；
- 数据集名称、下载日期、许可证/访问限制和原始文件校验值；
- 训练、验证、测试划分文件；
- Python、CUDA、PyTorch、RDKit/Open Babel、ViennaRNA 等版本；
- 随机种子、训练轮数、批大小和硬件；
- 原始命令、标准输出、错误日志、模型权重和最终指标；
- 对官方代码所做的每一处修改。

建议每个项目至少提供 `environment.yml`（或锁定版本的 `requirements.txt`）、`run_train.*`、`run_test.*`、`data/README.md` 和 `results/README.md`。数据受许可证限制时，`data/README.md` 只保存来源、下载步骤和校验方式，不提交原始数据。

## 12. 与综述统一评测的关系

这些方法并非使用完全相同的原始输入：

- **序列/表示类**：RSAPred、DeepRSMA、BioLLMNet，综述中以 R-SIM 派生任务进行统一比较。
- **三维结构类**：RLaffinity、RLASIF、SPRank，综述中以 PDBbind NL2020 派生的 95 个 RNA–配体复合物 affinity 集合进行比较。
- **binding preference / ligand ranking 类**：RNAmigos、R-BIND、RNALigands 和 ZHMol-RLinter，图中的比较使用 UNK96；前三者按 top-10 ranking 判定，ZHMol-RLinter 按是否正确分类为 binding 判定。

因此，复现“原论文结果”和复现“综述中的统一测试结果”应作为两个独立任务。前者严格遵循各论文数据划分，后者必须额外建立一致的数据 ID、标签单位、去重规则和 train/validation/test 划分，不能直接比较各仓库默认输出。
