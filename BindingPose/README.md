# RNA–小分子 Binding Pose 方法复现说明

本文档整理图中的 RNA–小分子 binding pose 方法。其中 **LigandRNA、DrugScoreRNA、RmsdXNA 和 SPRank** 是对候选构象进行打分或重排的方法；**NLDock、AutoDock 4、rDock 和 DOCK 6** 是生成并排序对接构象的框架。两类方法的输入和评估流程不同，不能把打分程序当成完整对接程序使用。

> 核对日期：2026-09-17。公开仓库、下载页和许可条款可能后续变化。

## 1. 公共数据与评估规则

图中的 Yan、Ruiz、Chen 和 Philips 是四套不同来源的 RNA–ligand pose benchmark，不是某个对接软件自带的通用训练集。复现时应优先使用原论文补充材料中的 PDB ID、结构清理规则、对接参数和成功判定标准。

- 四套 benchmark 的统一比较与对接参数：[NLDock 论文及 Supporting Information](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00341)
- 四套数据的列表与重排评估：[RmsdXNA 论文及补充材料](https://doi.org/10.1093/bib/bbae166)
- Philips set 的来源和 LigandRNA/Dock6 候选构象流程：[LigandRNA 论文](https://doi.org/10.1261/rna.039834.113)
- 原始共晶结构：[RCSB Protein Data Bank](https://www.rcsb.org/)

补充材料主要提供 PDB ID、参数和结果表，不等于已准备好的受体、配体和全部 decoy pose 文件。若需要严格重建表中数值，必须保持下列条件一致：

1. 受体链、配体、金属离子、结构水和质子化状态；
2. 原子类型、部分电荷、可旋转键和网格/口袋范围；
3. 软件版本、随机种子、每个配体的搜索次数和输出 pose 数；
4. 对称原子处理方式和 heavy-atom RMSD 计算程序；
5. 论文定义的 top-1 success 阈值，通常为排名第一的 pose 与晶体配体的 heavy-atom RMSD 不超过 2 Å。

## 2. LigandRNA

- 论文：[LigandRNA: computational predictor of RNA–ligand interactions](https://doi.org/10.1261/rna.039834.113)
- 全文与补充材料：[PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC3860260/)
- 论文声明的 Web 服务：`http://ligandrna.genesilico.pl/`
- 论文声明的镜像：`http://ligandrna.biol.amu.edu.pl/`

### 无法运行的原因

LigandRNA 只发布为 Web 服务，当前找不到公开源码、独立可执行文件、统计势参数文件或容器镜像。核对日期访问主站返回 HTTP 503，镜像站也无法建立有效响应。论文描述了以 RNA PDB 和候选配体 MOL2 为输入的排序流程，但补充材料不包含可代替服务器的程序和势函数文件，因此不能写出可验证的本地运行命令。

## 3. DrugScoreRNA

- 论文：[DrugScoreRNA—Knowledge-Based Scoring Function To Predict RNA−Ligand Interactions](https://doi.org/10.1021/ci700134p)
- 作者课题组软件页：[Computational Pharmaceutical Chemistry Lab – Software](https://cpclab.uni-duesseldorf.de/index.php/Software)
- 数据说明：论文 Supporting Information 的 Table S1/S2 列出了构建势函数和对接验证使用的 PDB 结构

### 无法运行的原因

作者的官方软件页仍列出 DrugScoreRNA，但没有提供下载链接、源码仓库、可执行文件或运行手册。论文公开了方法和 PDB ID，却没有发布从 670 个核酸复合物推导的距离相关 pair-potential 参数、原子类型映射和 AutoDock 网格转换程序。缺少这些核心文件时，无法用原方法对 pose 打分或重建论文实验。

## 4. RmsdXNA

- 论文：[RmsdXNA: RMSD prediction of nucleic acid–ligand docking poses using machine-learning method](https://doi.org/10.1093/bib/bbae166)
- 全文与补充材料：[PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11063749/)
- 论文声明的源码地址：`https://github.com/laiheng001/RmsdXNA`

### 无法运行的原因

论文声明的 GitHub 仓库在核对日期返回 404，也没有找到可替代的官方归档。论文补充材料可用于核对 benchmark 和 PDB 列表，但不包含距离特征生成脚本、原子类型规则、训练好的 XGBoost 模型、锁定依赖或推理入口。只根据论文重写一个 XGBoost 回归器不是运行原作者代码，无法作为严格复现。

## 5. SPRank

- 论文：[SPRank: a knowledge-based scoring function for RNA–ligand complexes](https://doi.org/10.1021/acs.jctc.4c00681)
- 全文与补充材料：[PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12960052/)
- 论文声明的源码地址：`https://github.com/Vfold-RNA/SPRank`

### 无法运行的原因

论文声明的 GitHub 仓库在核对日期返回 404。补充材料可以核对训练/测试数据组成和结果，但没有 standalone 程序、统计势文件、特征生成器、随机森林权重、环境文件或命令行说明。尽管论文描述了使用 rDock/AutoDock Vina 生成候选构象再由 SPRank 重排的流程，但缺少模型与势函数后不能执行原方法。

## 6. NLDock

- 论文：[NLDock: a Fast Nucleic Acid–Ligand Docking Algorithm for Modeling RNA/DNA–Ligand Complexes](https://doi.org/10.1021/acs.jcim.1c00341)
- 官方下载页：[NLDock v1.0](http://huanglab.phys.hust.edu.cn/software/NLDock/)
- 数据与参数：论文 Supporting Information 提供 NLDock、AutoDock、rDock 和 DOCK 6 的对接参数与绑定位点定义

NLDock 不是公开 Git 仓库。官方页面要求填写姓名、电子邮箱和机构，并接受非商业使用条款后才能下载。这些信息必须由使用者本人提交，不应在 GitHub 中转存获取到的受限二进制包。

### 获取与运行

1. 在官方下载页完成登记和许可确认，下载 NLDock v1.0 的 Linux 软件包。
2. 解压后先运行软件包自带的示例，核对可执行文件、辅助工具和数据文件是否完整。
3. 按软件包手册准备 RNA/DNA 受体 PDB 和配体构象 MOL2。局部对接需先用包内工具生成 binding-site sphere points；全局对接使用包内对应模式。
4. 使用下载包 README/手册中与 v1.0 对应的命令行执行对接，然后用包内排序工具处理输出 pose。
5. 复现表中数据时，用论文 Supporting Information 的 binding-site 定义和参数替换默认示例参数。

官方页在下载前不公开软件包手册，因此本文档不猜测二进制文件名和命令行参数。下载后应记录包版本、文件校验值和手册中的官方命令。

## 7. AutoDock 4

- 论文：[AutoDock4 and AutoDockTools4: Automated Docking with Selective Receptor Flexibility](https://doi.org/10.1002/jcc.21256)
- 源码：[ccsb-scripps/AutoDock4](https://github.com/ccsb-scripps/AutoDock4)
- AutoGrid 源码：[ccsb-scripps/AutoGrid](https://github.com/ccsb-scripps/AutoGrid)
- 官方下载：[AutoDock 4](https://autodock.scripps.edu/download-autodock4/)
- 手册：[AutoDock 4.2.6 User Guide](https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf)
- benchmark 数据：按 NLDock Supporting Information 的 PDB ID 从 RCSB PDB 下载，并使用其 AutoDock 参数

AutoDock 4 需要 AutoGrid4 生成受体能量网格，AutoDock4 执行搜索，AutoDockTools/MGLTools 负责生成 PDBQT、GPF 和 DPF。下面是通用的本地对接流程；严格复现表中结果时，还必须把 GPF/DPF 改为论文补充材料中的 RNA 参数。

### 安装

从官方下载页获取 AutoDock 4.2.6、AutoGrid 4.2.6 和 AutoDockTools/MGLTools。如需从源码构建，分别克隆 AutoDock4 和 AutoGrid 仓库，按各仓库 `INSTALL` 文件编译。安装后先确认以下命令可用：

```bash
autodock4 -h
autogrid4 -h
```

### 准备与对接

`pythonsh` 和 `Utilities24` 由 MGLTools/AutoDockTools 提供：

```bash
pythonsh Utilities24/prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
pythonsh Utilities24/prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt
pythonsh Utilities24/prepare_gpf4.py -l ligand.pdbqt -r receptor.pdbqt -o receptor.gpf
pythonsh Utilities24/prepare_dpf4.py -l ligand.pdbqt -r receptor.pdbqt -o ligand_receptor.dpf

autogrid4 -p receptor.gpf -l receptor.glg
autodock4 -p ligand_receptor.dpf -l ligand_receptor.dlg
```

运行 AutoGrid 前必须核对 `receptor.gpf` 中的网格中心、尺寸和原子类型；运行 AutoDock 前必须核对 `ligand_receptor.dpf` 的搜索次数、随机种子和算法参数。RNA 的磷酸根电荷、金属离子和结构水不应在不核对原论文的情况下自动删除。最终 pose 和聚类结果从 `ligand_receptor.dlg` 中读取或由 AutoDockTools 导出。

## 8. rDock

- 论文：[rDock: a fast, versatile and open source program for docking ligands to proteins and nucleic acids](https://doi.org/10.1371/journal.pcbi.1003571)
- 源码：[CBDD/rDock](https://github.com/CBDD/rDock)
- 文档：[rDock Documentation](https://rdock.github.io/documentation/)
- 三步对接教程：[Docking in 3 steps](https://rdock.github.io/docking-in-3-steps/)
- benchmark 数据：按 NLDock Supporting Information 的 PDB ID 从 RCSB PDB 下载，并使用其 rDock 参数

### 安装

建议在 Linux 下从源码构建，并使用当前仓库自带的测试验证二进制文件：

```bash
sudo apt update
sudo apt install -y make git libpopt0 libpopt-dev g++
git clone https://github.com/CBDD/rDock.git
cd rDock
make -j 4
make test

PREFIX=/path/to/rdock-install make install
export PATH=/path/to/rdock-install/bin:$PATH
export LD_LIBRARY_PATH=/path/to/rdock-install/lib:$LD_LIBRARY_PATH
export RBT_ROOT=/path/to/rdock-install/rDock
```

### 准备与对接

rDock 需要 RNA 受体 MOL2、配体 SDF、一个描述受体和口袋的 `receptor.prm`，以及对接协议 `dock.prm`。`receptor.prm` 的完整模板见官方三步教程。

```bash
rbcavity -was -d -r receptor.prm
rbdock -i ligands.sd -o rdock_output -r receptor.prm -p dock.prm -n 50
```

`rbcavity` 生成 `.as` 口袋文件和可视化网格；应先在 PyMOL 等程序中确认口袋位置。`rbdock` 输出 `rdock_output.sd`，其中包含每个 pose 的分数字段。复现原论文时应使用固定 commit，并通过 `-s` 记录随机种子；表中旧版 rDock/RiboDock 参数与当前仓库默认参数可能不同，以 NLDock Supporting Information 为准。

## 9. DOCK 6

- RNA 对接论文：[DOCK 6: combining techniques to model RNA–small molecule complexes](https://doi.org/10.1261/rna.1563609)
- 源码：[docking-org/dock6](https://github.com/docking-org/dock6)
- 官方网站：[UCSF DOCK 6](https://dock.docking.org/DOCK_6/index.htm)
- 手册：[DOCK 6.13 User Manual](https://dock.docking.org/DOCK_6/dock6_manual.htm)
- 教程：[DOCK 6 Tutorials](https://dock.docking.org/DOCK_6/tutorials/index.htm)
- benchmark 数据：按 NLDock Supporting Information 的 PDB ID 从 RCSB PDB 下载，并使用其 DOCK 6 参数

### 安装

DOCK 6 是 Unix/Linux 软件。当前公开仓库为 DOCK 6.13.1；复现旧论文时应在 GitHub Releases 选择相应版本，不要默认将最新版结果与表中数值等同。

```bash
git clone https://github.com/docking-org/dock6.git
cd dock6/install
./configure gnu
make install
make test
```

### 准备与对接

DOCK 6 的标准局部对接流程是：准备带电荷的 RNA 受体和配体 MOL2，为受体生成分子表面，使用 `sphgen` 生成 spheres，用 `sphere_selector` 保留绑定口袋附近的 spheres，然后通过 `showbox` 和 `grid` 生成对接网格。具体输入文件模板见官方 flexible ligand tutorial。

```bash
dms receptor.mol2 -n -w 1.4 -v -o receptor.ms
sphgen -i INSPH -o OUTSPH
sphere_selector receptor.sph reference_ligand.mol2 10.0
showbox < showbox.in
grid -i grid.in -o grid.out
dock6 -i dock.in -o dock.out
```

`dock.in`、`grid.in`、`INSPH` 和 `showbox.in` 不能通用，必须由当前受体、参考配体和论文参数生成。最终 pose 通常以 MOL2 输出；应使用与论文一致的 heavy-atom RMSD 程序进行排名第一 pose 的成功判定。

## 10. 严格复现时的文件记录

对可运行的对接框架，每个方法至少应保存：

- 论文 DOI、源码 URL、commit SHA/软件版本和获取日期；
- PDB ID、原始下载文件校验值与结构清理日志；
- 受体、参考配体、候选配体和口袋定义文件；
- 原子类型、电荷、质子化、金属/水处理和可旋转键记录；
- 完整配置文件、命令、标准输出、错误日志和随机种子；
- 所有输出 pose、原始分数、重排结果和 RMSD 程序版本。

表中的 scoring-function 比较应对同一批候选 pose 重排；docking-framework 比较应分别记录每个框架自己生成的 pose 集。如果候选 pose、口袋或成功阈值不同，就不能将得到的 success rate 直接与图中数值比较。
