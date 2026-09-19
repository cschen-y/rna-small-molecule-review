# Reproducing RNA–Small-Molecule Binding-Pose Methods

This document covers the RNA–small-molecule binding-pose methods shown in the review. **LigandRNA, DrugScoreRNA, RmsdXNA, and SPRank** score or rerank candidate poses. **NLDock, AutoDock 4, rDock, and DOCK 6** are docking frameworks that generate and rank poses. These two method classes use different inputs and evaluation workflows. A scoring program should not be treated as a complete docking program.

> Last verified: 2026-09-17. Public repositories, download pages, and licensing terms may change after this date.

## 1. Shared Data and Evaluation Protocol

Yan, Ruiz, Chen, and Philips are four RNA–ligand pose benchmarks from different sources. They are not generic training sets distributed with a particular docking package. Reproduction should use the PDB IDs, structure-cleaning rules, docking parameters, and success criteria reported in the supporting information of the original papers.

- Harmonized comparison and docking parameters for all four benchmarks: [NLDock paper and Supporting Information](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00341)
- Dataset lists and reranking evaluation for all four benchmarks: [RmsdXNA paper and supplementary material](https://doi.org/10.1093/bib/bbae166)
- Origin of the Philips set and the LigandRNA/DOCK 6 candidate-pose workflow: [LigandRNA paper](https://doi.org/10.1261/rna.039834.113)
- Original cocrystal structures: [RCSB Protein Data Bank](https://www.rcsb.org/)

The supplementary materials mainly provide PDB IDs, parameters, and result tables. They are not collections of prepared receptors, ligands, and complete decoy-pose files. A strict reconstruction of the tabulated values requires consistency in all of the following:

1. receptor chains, ligands, metal ions, crystallographic water molecules, and protonation states;
2. atom types, partial charges, rotatable bonds, and grid or pocket boundaries;
3. software versions, random seeds, number of searches per ligand, and number of output poses;
4. treatment of symmetry-equivalent atoms and the program used to calculate heavy-atom RMSD;
5. the top-1 success threshold defined in the paper, which is usually a heavy-atom RMSD of no more than 2 Å between the highest-ranked pose and the crystallographic ligand.

## 2. LigandRNA

- Paper: [LigandRNA: computational predictor of RNA–ligand interactions](https://doi.org/10.1261/rna.039834.113)
- Full text and supplementary material: [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC3860260/)
- Web server reported in the paper: `http://ligandrna.genesilico.pl/`
- Mirror reported in the paper: `http://ligandrna.biol.amu.edu.pl/`

### Why the method cannot currently be run

LigandRNA was released only as a web service. No public source code, standalone executable, statistical-potential parameter file, or container image could be found. On the verification date, the primary server returned HTTP 503 and the mirror did not provide a valid response. The paper describes a ranking workflow that accepts an RNA PDB structure and candidate ligands in MOL2 format, but the supplementary material does not contain a program or potential files that can replace the server. A verifiable local execution command therefore cannot be provided.

## 3. DrugScoreRNA

- Paper: [DrugScoreRNA—Knowledge-Based Scoring Function To Predict RNA−Ligand Interactions](https://doi.org/10.1021/ci700134p)
- Software page of the authors' laboratory: [Computational Pharmaceutical Chemistry Lab – Software](https://cpclab.uni-duesseldorf.de/index.php/Software)
- Data description: Tables S1 and S2 in the Supporting Information list the PDB structures used to derive the potential and validate docking

### Why the method cannot currently be run

The official software page still lists DrugScoreRNA but provides no download link, source-code repository, executable, or user manual. The paper describes the method and reports PDB IDs, but it does not release the distance-dependent pair-potential parameters derived from 670 nucleic-acid complexes, the atom-type mapping, or the AutoDock grid-conversion program. Without these core files, the original method cannot score poses or reproduce the reported experiments.

## 4. RmsdXNA

- Paper: [RmsdXNA: RMSD prediction of nucleic acid–ligand docking poses using machine-learning method](https://doi.org/10.1093/bib/bbae166)
- Full text and supplementary material: [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11063749/)
- Source-code URL reported in the paper: `https://github.com/laiheng001/RmsdXNA`

### Why the method cannot currently be run

The GitHub repository reported in the paper returned 404 on the verification date, and no alternative official archive was found. The supplementary material can be used to verify the benchmark and PDB lists, but it does not contain the distance-feature generation scripts, atom-typing rules, trained XGBoost model, pinned dependencies, or inference entry point. Reimplementing an XGBoost regressor from the paper is not equivalent to running the authors' code and would not constitute a strict reproduction.

## 5. SPRank

- Paper: [SPRank: a knowledge-based scoring function for RNA–ligand complexes](https://doi.org/10.1021/acs.jctc.4c00681)
- Full text and supporting information: [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12960052/)
- Source-code URL reported in the paper: `https://github.com/Vfold-RNA/SPRank`

### Why the method cannot currently be run

The GitHub repository reported in the paper returned 404 on the verification date. The supplementary material identifies the training and test sets and reports results, but it does not provide a standalone program, statistical-potential files, a feature generator, random-forest weights, an environment file, or command-line instructions. The paper describes generating candidate poses with rDock or AutoDock Vina and reranking them with SPRank. The original method cannot be executed without the missing model and potential files.

## 6. NLDock

- Paper: [NLDock: a Fast Nucleic Acid–Ligand Docking Algorithm for Modeling RNA/DNA–Ligand Complexes](https://doi.org/10.1021/acs.jcim.1c00341)
- Official download page: [NLDock v1.0](http://huanglab.phys.hust.edu.cn/software/NLDock/)
- Data and parameters: the Supporting Information provides docking parameters and binding-site definitions for NLDock, AutoDock, rDock, and DOCK 6

NLDock is not distributed through a public Git repository. Its official page requires the user to provide a name, email address, and affiliation and to accept noncommercial-use terms before downloading. Each user must submit this information personally. The restricted binary package should not be redistributed through GitHub.

### Obtain and run

1. Complete the registration and license confirmation on the official download page and download the NLDock v1.0 Linux package.
2. After extraction, run the included example first and confirm that the executable, auxiliary tools, and data files are complete.
3. Prepare an RNA or DNA receptor in PDB format and ligand conformations in MOL2 format according to the package manual. For local docking, first generate binding-site sphere points with the tools in the package. Use the corresponding package mode for global docking.
4. Run docking with the v1.0 command documented in the downloaded README or manual, then process the output poses with the package's ranking tool.
5. To reproduce the table, replace the example parameters with the binding-site definitions and parameters from the paper's Supporting Information.

The official page does not expose the package manual before download. This document therefore does not guess executable names or command-line arguments. After downloading, record the package version, file checksums, and official commands from the manual.

## 7. AutoDock 4

- Paper: [AutoDock4 and AutoDockTools4: Automated Docking with Selective Receptor Flexibility](https://doi.org/10.1002/jcc.21256)
- Source code: [ccsb-scripps/AutoDock4](https://github.com/ccsb-scripps/AutoDock4)
- AutoGrid source code: [ccsb-scripps/AutoGrid](https://github.com/ccsb-scripps/AutoGrid)
- Official download: [AutoDock 4](https://autodock.scripps.edu/download-autodock4/)
- Manual: [AutoDock 4.2.6 User Guide](https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf)
- Benchmark data: download structures from RCSB PDB using the PDB IDs in the NLDock Supporting Information and use the AutoDock parameters reported there

AutoDock 4 uses AutoGrid4 to generate receptor energy grids and AutoDock4 to perform the search. AutoDockTools or MGLTools prepares PDBQT, GPF, and DPF files. The commands below show a general local docking workflow. To reproduce the table strictly, the GPF and DPF files must use the RNA-specific parameters from the paper's supplementary material.

### Install

Obtain AutoDock 4.2.6, AutoGrid 4.2.6, and AutoDockTools or MGLTools from the official download page. To build from source, clone the AutoDock4 and AutoGrid repositories separately and follow the `INSTALL` file in each repository. Confirm that the commands are available after installation:

```bash
autodock4 -h
autogrid4 -h
```

### Prepare and dock

MGLTools or AutoDockTools provides `pythonsh` and `Utilities24`:

```bash
pythonsh Utilities24/prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt
pythonsh Utilities24/prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt
pythonsh Utilities24/prepare_gpf4.py -l ligand.pdbqt -r receptor.pdbqt -o receptor.gpf
pythonsh Utilities24/prepare_dpf4.py -l ligand.pdbqt -r receptor.pdbqt -o ligand_receptor.dpf

autogrid4 -p receptor.gpf -l receptor.glg
autodock4 -p ligand_receptor.dpf -l ligand_receptor.dlg
```

Before running AutoGrid, verify the grid center, grid dimensions, and atom types in `receptor.gpf`. Before running AutoDock, verify the number of searches, random seed, and algorithm parameters in `ligand_receptor.dpf`. Do not remove RNA phosphate charges, metal ions, or crystallographic water molecules automatically without checking the protocol in the original paper. Read final poses and clustering results from `ligand_receptor.dlg` or export them with AutoDockTools.

## 8. rDock

- Paper: [rDock: a fast, versatile and open source program for docking ligands to proteins and nucleic acids](https://doi.org/10.1371/journal.pcbi.1003571)
- Source code: [CBDD/rDock](https://github.com/CBDD/rDock)
- Documentation: [rDock Documentation](https://rdock.github.io/documentation/)
- Tutorial: [Docking in 3 steps](https://rdock.github.io/docking-in-3-steps/)
- Benchmark data: download structures from RCSB PDB using the PDB IDs in the NLDock Supporting Information and use the rDock parameters reported there

### Install

Building from source on Linux is recommended. Validate the resulting binaries with the tests included in the current repository:

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

### Prepare and dock

rDock requires an RNA receptor in MOL2 format, ligands in SDF format, a `receptor.prm` file describing the receptor and pocket, and a `dock.prm` docking protocol. A complete `receptor.prm` template is available in the official three-step tutorial.

```bash
rbcavity -was -d -r receptor.prm
rbdock -i ligands.sd -o rdock_output -r receptor.prm -p dock.prm -n 50
```

`rbcavity` produces an `.as` pocket file and a visualization grid. Inspect the pocket location in PyMOL or another molecular viewer before docking. `rbdock` writes `rdock_output.sd`, which contains score fields for each pose. To reproduce an original study, use a fixed commit and record the random seed with `-s`. The older rDock or RiboDock parameters used in the table may differ from current repository defaults. Follow the NLDock Supporting Information for the benchmark parameters.

## 9. DOCK 6

- RNA docking paper: [DOCK 6: combining techniques to model RNA–small molecule complexes](https://doi.org/10.1261/rna.1563609)
- Source code: [docking-org/dock6](https://github.com/docking-org/dock6)
- Official website: [UCSF DOCK 6](https://dock.docking.org/DOCK_6/index.htm)
- Manual: [DOCK 6.13 User Manual](https://dock.docking.org/DOCK_6/dock6_manual.htm)
- Tutorials: [DOCK 6 Tutorials](https://dock.docking.org/DOCK_6/tutorials/index.htm)
- Benchmark data: download structures from RCSB PDB using the PDB IDs in the NLDock Supporting Information and use the DOCK 6 parameters reported there

### Install

DOCK 6 runs on Unix or Linux. The current public repository contains DOCK 6.13.1. When reproducing an older paper, select the corresponding version from GitHub Releases rather than treating results from the latest version as equivalent to the tabulated values.

```bash
git clone https://github.com/docking-org/dock6.git
cd dock6/install
./configure gnu
make install
make test
```

### Prepare and dock

A standard local-docking workflow in DOCK 6 prepares charged RNA receptor and ligand MOL2 files, generates a molecular surface for the receptor, creates spheres with `sphgen`, retains spheres near the binding pocket with `sphere_selector`, and generates docking grids with `showbox` and `grid`. Input-file templates are available in the official flexible-ligand tutorial.

```bash
dms receptor.mol2 -n -w 1.4 -v -o receptor.ms
sphgen -i INSPH -o OUTSPH
sphere_selector receptor.sph reference_ligand.mol2 10.0
showbox < showbox.in
grid -i grid.in -o grid.out
dock6 -i dock.in -o dock.out
```

The `dock.in`, `grid.in`, `INSPH`, and `showbox.in` files are receptor-specific and cannot be reused unchanged. Generate them from the current receptor, reference ligand, and paper-specific parameters. Final poses are commonly written in MOL2 format. Evaluate the highest-ranked pose with the same heavy-atom RMSD program used in the paper.
