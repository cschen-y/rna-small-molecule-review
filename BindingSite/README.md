# RNA–小分子 Binding Site 方法复现说明

本目录统一整理 **RLBind、RNet、RBind、Rsite、Rsite2、RNAsite、RNABind、MultiModRLBP 和 MVRBind** 九个 RNA–小分子 binding-site 预测方法，并为每个方法提供统一的 `run_test17.py`。

> 核对日期：2026-09-17。公开仓库、网站和依赖包可能后续变化。

## 1. Test17 数据

本目录使用统一的 Test17：从公开 Test18 中删除 `6EZ0A` 后得到 17 条 RNA 链，共 583 个核苷酸。所有方法的标签和评估顺序以这一版本为准。

共用数据位于 `data/common/`：

- `test17.fasta`：合并的 Test17 FASTA；
- `test17_fastas/`：17 个按 PDB ID 和链命名的 FASTA；
- `train60_fastas/`：Train60 训练集 FASTA；
- `pdbFiles/`：统一使用的 RNA 结构；
- `all_label/label/train_label.pkl`：Train60 标签；
- `all_label/label/test17_labels.pkl`：Test17 标签。

方法专用特征位于 `data/<method>/`，不与其他方法混用。`all_label` 以及 RBind、Rsite、RNet、MultiModRLBP 目录中的 `data`/标签入口是 Windows junction，已指向本目录下的实际数据。若通过 ZIP 打包或在不保留 junction 的文件系统上复制，需要重新建立这些链接。

## 2. 环境

九个方法需要三个 Conda 环境，原因是 Python、Biopython、PyTorch 和 CUDA 版本约束不兼容。

| 环境 | 方法 | Python | PyTorch | CUDA |
|---|---|---:|---:|---:|
| `core` | RLBind、RNet、RBind、Rsite、Rsite2、RNAsite | 3.9.21 | 2.1.0 | 11.8 |
| `MVRBind` | MVRBind、RNABind | 3.10.18 | 2.2.0 | 12.1 |
| `Mul` | MultiModRLBP | 3.7.12 | 1.13.1 | 11.6 |

在 `BindingSite` 目录下创建环境：

```bash
conda env create -f environments/core.yml
conda env create -f environments/mvrbind.yml
conda env create -f environments/multimodrlbp.yml
```

MultiModRLBP 包含 Windows CPython 3.7 的 `alignment_C` 扩展。Linux 或二进制不兼容时，在 `Mul` 环境中重新编译：

```bash
conda activate Mul
cd MultiModRLBP/RnaBert
python setup.py build_ext --inplace
cd ../..
```

更完整的版本和 GPU 说明见 `ENVIRONMENTS.md`。

## 3. 统一运行与校验

以下命令必须在 `BindingSite` 目录下执行。

运行全部方法：

```bash
python run_all_test17.py --methods all --continue-on-error
```

只运行指定方法：

```bash
python run_all_test17.py --methods rbind rsite rsite2 mvrbind
```

总运行脚本会根据方法自动选择 `core`、`MVRBind` 或 `Mul` 环境。Conda 不在 `PATH` 时，设置 `CONDA_EXE` 或传入 `--conda /path/to/conda`。

检查目录、Python 语法、Train60/Test17 数量、标签长度和结果 JSON：

```bash
python validate_review.py
```

汇总已生成的结果：

```bash
python collect_results.py
```

汇总表写入 `results/test17_summary.csv`，每个方法的完整输出位于 `results/<method>/test17_results.json`。

## 4. RLBind

- 论文：[RLBind: a deep learning method to predict RNA–ligand binding sites](https://doi.org/10.1093/bib/bbac486)
- 原始源码与数据：[KailiWang1/RLBind](https://github.com/KailiWang1/RLBind)
- 本地代码：`RLBind/`
- 本地数据：`data/RLBind/`
- 环境：`core`

RLBind 融合全局 RNA 序列特征与局部邻域特征。本地 Test17 入口使用 Train60 重新训练，默认 60 epochs 和 5 个随机种子。

```bash
conda run --no-capture-output -n core python RLBind/run_test17.py
```

只做单轮管线检查：

```bash
conda run --no-capture-output -n core python RLBind/run_test17.py --epochs 1 --seeds 8124
```

## 5. RNet

- 论文：[RNet: a network strategy to predict RNA binding preferences](https://doi.org/10.1093/bib/bbad482)
- 作者声明的源码与数据页：[RNetsite](http://zhaoserver.com.cn/RNet/RNet.html)
- 本地代码：`RNet/`
- 本地数据：`data/RNet/`
- 环境：`core`

RNet 中的 RNetsite 从 RNA 三维接触网络提取局部和全局网络特征，再使用传统机器学习模型预测核苷酸级 binding site。本地入口运行 5 个配置种子。

```bash
conda run --no-capture-output -n core python RNet/run_test17.py
```

## 6. RBind

- 论文：[RBind: computational network method to predict RNA binding sites](https://doi.org/10.1093/bioinformatics/bty345)
- 原作者代码与数据页：[RBind](https://zhaolab.com.cn/RBind)
- 本地代码：`RBind/`
- 本地数据：`data/common/`
- 环境：`core`

RBind 将核苷酸表示为接触网络节点，基于拓扑特征和统计阈值预测 binding site。Test17 入口是确定性评估。

```bash
conda run --no-capture-output -n core python RBind/run_test17.py
```

## 7. Rsite

- 论文：[Rsite: a computational method to identify the functional sites of noncoding RNAs](https://doi.org/10.1038/srep09179)
- 原始源码和示例数据：[Rsite official page](https://www.cuilab.cn/rsite)
- 本地代码：`Rsite/`
- 本地数据：`data/common/`
- 环境：`core`

Rsite 根据 RNA 三维坐标计算核苷酸到质心的距离曲线，用平滑后的极值点识别候选功能位点。

```bash
conda run --no-capture-output -n core python Rsite/run_test17.py
```

## 8. Rsite2

- 论文：[Rsite2: an efficient computational method to predict the functional sites of noncoding RNAs](https://doi.org/10.1038/srep19016)
- 原始源码与数据：[Rsite official page](https://www.cuilab.cn/rsite)
- 本地代码：`Rsite2/`
- 本地数据：`data/Rsite2/`
- 环境：`core`

Rsite2 用 RNA 二级结构的二维坐标和距离曲线代替 Rsite 对三维结构的要求。本地 Test17 入口读取预计算的 `SS_NDS` 并输出统一指标。

```bash
conda run --no-capture-output -n core python Rsite2/run_test17.py
```

## 9. RNAsite

- 论文：[Recognition of small molecule–RNA binding sites using RNA sequence and structure](https://doi.org/10.1093/bioinformatics/btaa1092)
- 官方 Web 服务：[RNAsite](https://yanglab.qd.sdu.edu.cn/RNAsite/)
- 本地代码：`RNAsite/`
- 本地数据：`data/RNAsite/` 和 `data/common/`
- 环境：`core`

RNAsite 结合序列保守性、三维网络拓扑、Laplacian norm 和 SASA 特征，并通过随机森林预测。本地入口使用 Train60/Test17 预处理特征并运行 5 个配置种子。

```bash
conda run --no-capture-output -n core python RNAsite/run_test17.py
```

## 10. RNABind

- 论文：[Identifying RNA-small Molecule Binding Sites Using Geometric Deep Learning with Language Models](https://doi.org/10.1016/j.jmb.2025.169010)
- 原始源码与数据：[jaminzzz/RNABind](https://github.com/jaminzzz/RNABind)
- 本地代码：`RNABind/`
- 本地数据：`data/common/`
- 环境：`MVRBind`

RNABind 是基于 RNA 结构图和 EGNN 的方法，原论文还评估了多种 RNA language-model embedding。当前本地整理包不含 ERNIE-RNA 等大模型权重，因此 Test17 脚本明确使用官方 RNABind EGNN 架构的 one-hot 输入变体，不会隐式下载其他模型。默认训练 30 epochs 和 5 个随机种子。

```bash
conda run --no-capture-output -n MVRBind python RNABind/run_test17.py
```

单轮管线检查：

```bash
conda run --no-capture-output -n MVRBind python RNABind/run_test17.py --epochs 1 --seeds 8124
```

## 11. MultiModRLBP

- 论文：[MultiModRLBP: A Deep Learning Approach for Multi-Modal RNA-Small Molecule Ligand Binding Sites Prediction](https://doi.org/10.1109/JBHI.2024.3400521)
- 原始源码与数据：[lennylv/MultiModRLBP](https://github.com/lennylv/MultiModRLBP)
- 本地代码：`MultiModRLBP/`
- 本地数据：`data/MultiModRLBP/`
- 环境：`Mul`

MultiModRLBP 融合核苷酸级三维特征、RNA 关系图和 RNABert 序列表示。本地数据已包含 Test17 所需预处理特征和模型文件；默认训练 100 epochs。

```bash
conda run --no-capture-output -n Mul python MultiModRLBP/run_test17.py
```

单轮管线检查：

```bash
conda run --no-capture-output -n Mul python MultiModRLBP/run_test17.py --epochs 1
```

## 12. MVRBind

- 论文：[MVRBind: multi-view learning for RNA-small molecule binding site prediction](https://doi.org/10.1093/bib/bbaf489)
- 原始源码与数据：[cschen-y/MVRBind](https://github.com/cschen-y/MVRBind)
- 本地代码：`MVRBind/`
- 本地数据：`data/MVRBind/`
- 环境：`MVRBind`

MVRBind 将 RNA 一级、二级和三级结构建模为多视图，并在多个空间尺度上融合节点表示。本地 Test17 入口读取 `data/MVRBind/pt/` 的预处理图，删除 Test18 中的 `6EZ0A`，并评估 `data/MVRBind/model_parameters/` 中的 5 个预训练 checkpoint。

```bash
conda run --no-capture-output -n MVRBind python MVRBind/run_test17.py
```

## 13. 结果解读

所有方法统一输出 Accuracy、Precision、Recall、F1、MCC、AUC、AUPR 和 BACC。比较前必须同时核对结果 JSON 中的训练轮数、随机种子、模型变体和是否使用预训练 checkpoint。

已保存的 RLBind、RNABind 和 MultiModRLBP 单轮结果只用于证明管线能完整执行，不是论文配置。正式重跑时不传入 `--epochs 1`，使用各脚本的默认训练轮数。RNet 和 RNAsite 使用 5 个配置种子；MVRBind 评估 5 个预训练 checkpoint；RBind、Rsite 和 Rsite2 为确定性评估。

## 14. GitHub 复现记录

提交到 GitHub 时，至少记录：

- 原论文 DOI、原始源码 URL、commit SHA 和下载日期；
- 环境 YAML、操作系统、GPU、NVIDIA driver 和 CUDA 版本；
- Train60/Test17 的 PDB ID、链 ID、去重规则与原始文件校验值；
- 随机种子、训练轮数、批大小、学习率和预训练权重校验值；
- 完整命令、标准输出、错误日志、原始 JSON 和汇总 CSV；
- 相对原作者代码的所有修改。

数据、权重或第三方代码受原许可限制时，只提交来源、下载脚本和校验方式，不应默认转发原文件。
