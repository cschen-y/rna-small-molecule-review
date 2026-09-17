# Machine Learning for RNA–Small Molecule Interactions

This repository accompanies the survey **“Machine Learning for RNA–Small Molecule Interactions: Methodologies, Evaluation, and Future Directions.”** It organizes reproducibility materials around three tasks: binding-site prediction, binding-preference prediction, and binding-pose prediction.

The scope differs by task. The binding-site results were generated in this study on the annotation-corrected Test17 benchmark and are accompanied by executable code, environments, inputs, raw outputs, and summary files. The binding-preference and binding-pose tables primarily summarize results reported in the original publications; their folders provide provenance, official code and data locations, and reproduction instructions where the required artifacts are publicly available. A literature-reported value is not presented as a rerun result.

## Repository contents

| Directory | Contents | Documentation |
|---|---|---|
| `BindingSite/` | Test17 benchmark package for nine methods, three Conda environments, per-method runners, unified runner, validator, raw JSON results, and summary CSV | [BindingSite/README.md](BindingSite/README.md) |
| `BindingPreference/` | Papers, source-code and data links, input requirements, execution instructions, and explicit blockers for methods whose complete artifacts are unavailable | [BindingPreference/README.md](BindingPreference/README.md) |
| `BindingPose/` | Benchmark provenance, pose-success protocol, scoring/rescoring methods, docking frameworks, installation steps, and explicit access or artifact limitations | [BindingPose/README.md](BindingPose/README.md) |

Redistributable third-party source code and data are included only when the applicable license or author terms permit redistribution. Otherwise, the documentation points to the official source and explains how to obtain the required artifact. Users remain responsible for complying with the licenses of the original methods and datasets.

## Task taxonomy

| Task | Prediction level | Typical outputs |
|---|---|---|
| Binding-site prediction | Nucleotide/residue | Binding probability or binary label for each RNA residue |
| Binding-preference prediction | RNA–ligand pair | Binding class, affinity, or ranking score |
| Binding-pose prediction | Three-dimensional complex | Generated or ranked ligand conformations in an RNA pocket |

Success rate is protocol-specific. In the UNK96 comparison, success means that the true ligand is retrieved within the top 10 for ranking methods, whereas for ZHMol-RLinter it means that a true binding pair is classified as binding; these values are not directly comparable. In the pose tables, success means that the highest-ranked ligand pose has RMSD ≤ 2.0 Å from the experimental pose.

## Binding-site benchmark

### Test17 definition

Test17 is the original Test18 collection with PDB entry `6EZ0A` excluded because U37 is an unnatural nucleotide incorporated into the RNA chain rather than a separate small-molecule ligand. Test17 contains 17 RNA chains and 583 nucleotides.

Test17 is an **annotation-corrected** version of Test18. The separate Test9 analysis found sequence identity up to 100% and RNA-align TM-score up to 0.987 relative to Train60. Only four Test9 RNA chains satisfy both the sequence- and structure-based independence criteria, and three of these (`1NEM`, `1Q8N`, and `2TOB`) are also present in Test18. Thus, only `1Y26` is independent of both Train60 and Test18. Because a one-sample set cannot support a stable benchmark comparison, the quantitative binding-site evaluation uses Test17.

### Included methods

The package contains Test17 entry points for:

- Rsite, Rsite2, and RBind;
- RNAsite and RNet;
- RLBind, RNABind, MultiModRLBP, and MVRBind.

Each method has its own `run_test17.py`. Shared benchmark files are under `BindingSite/data/common/`, while method-specific features are kept under `BindingSite/data/<method>/`.

### Environments

The nine methods require three environments because their Python, PyTorch, CUDA, and compiled-extension requirements are not mutually compatible.

| Environment | Methods | Python | PyTorch | CUDA toolkit |
|---|---|---:|---:|---:|
| `core` | RLBind, RNet, RBind, Rsite, Rsite2, RNAsite | 3.9.21 | 2.1.0 | 11.8 |
| `MVRBind` | MVRBind, RNABind | 3.10.18 | 2.2.0 | 12.1 |
| `Mul` | MultiModRLBP | 3.7.12 | 1.13.1 | 11.6 |

Create the environments from the `BindingSite` directory:

```bash
cd BindingSite
conda env create -f environments/core.yml
conda env create -f environments/mvrbind.yml
conda env create -f environments/multimodrlbp.yml
```

See `BindingSite/ENVIRONMENTS.md` for operating-system, GPU, CUDA, and extension-building details.

### Validate and run

All commands below are executed from `BindingSite/`.

Validate the directory layout, Python syntax, Train60/Test17 counts, label lengths, and result schema:

```bash
python validate_review.py
```

Run all nine methods. The unified runner selects the required Conda environment for each method:

```bash
python run_all_test17.py --methods all --continue-on-error
```

Run selected methods only:

```bash
python run_all_test17.py --methods rbind rsite rsite2 mvrbind
```

If Conda is not on `PATH`, set `CONDA_EXE` or pass `--conda /path/to/conda`. Full training can require a compatible NVIDIA GPU and may take substantially longer than a pipeline check. Method-specific commands, default epochs, seeds, checkpoints, and model variants are documented in `BindingSite/README.md`.

Collect all available result files:

```bash
python collect_results.py
```

Outputs are written to:

- `BindingSite/results/<method>/test17_results.json` for complete per-method results;
- `BindingSite/results/test17_summary.csv` for the consolidated table.

The JSON metadata should be checked before comparing methods, particularly the number of epochs, random seeds, model variant, and whether a pretrained checkpoint was used. One-epoch outputs are pipeline checks and should not be interpreted as paper-configuration results.

## Binding-preference reproduction guide

`BindingPreference/README.md` covers the methods represented in the review tables:

- affinity-oriented methods: BioLLMNet, DeepRSMA, RSAPred, RLaffinity, RLASIF, and SPRank;
- ranking or classification methods: RNAmigos, R-BIND, RNALigands, and ZHMol-RLinter.

For each method, the guide records the original paper, source-code location, dataset source, expected inputs, environment or dependency information, and execution steps. When a complete run cannot be reconstructed, the guide names the specific missing component, such as unavailable training code, model weights, data, a broken repository, or a distribution restriction.

The affinity comparison based on PDBbind NL2020 may require separately licensed PDBbind structures. These files are not redistributed here; use the acquisition and preprocessing instructions in the task README.

## Binding-pose reproduction guide

`BindingPose/README.md` covers:

- scoring or rescoring methods: LigandRNA, DrugScoreRNA, RmsdXNA, and SPRank;
- docking frameworks: NLDock, AutoDock 4, rDock, and DOCK 6.

The Yan, Ruiz, Chen, and Philips benchmarks have different source protocols. Strict reconstruction requires the original PDB list, receptor and ligand preparation, pocket definition, charge and protonation settings, software version, sampling parameters, random seeds, candidate poses, and RMSD implementation. The task README links the primary papers and supporting information and provides executable installation or docking steps where the software is publicly obtainable.

Scoring-function comparisons rerank a common candidate-pose set, whereas docking-framework comparisons include pose generation. These experiment types should not be conflated.

## Reproducibility records

For a publication-quality rerun, retain:

- the original paper DOI, source URL, commit SHA or software version, and acquisition date;
- environment lock files, operating system, GPU, driver, CUDA, and compiler information;
- input identifiers, downloaded-file checksums, structure-cleaning logs, and method-specific features;
- complete commands, configuration files, seeds, hyperparameters, checkpoints, standard output, and error logs;
- raw predictions, generated poses where applicable, evaluation outputs, and consolidated tables;
- a record of all changes made relative to the original implementation.

## Dataset observations

- Test18 contains the problematic `6EZ0A` entry described above.
- Test9 has substantial sequence and structural overlap with Train60. After accounting for overlap with both Train60 and Test18, only `1Y26` remains independent, which is insufficient for a statistically meaningful multi-method benchmark.
- Dataset curation rules, train–test separation, and task-specific metric definitions should be reported alongside benchmark results.

## Citation

If this repository supports your work, please cite the accompanying article:

> *Machine Learning for RNA–Small Molecule Interactions: Methodologies, Evaluation, and Future Directions.*

Please also cite each original method and dataset that you use. The task-specific READMEs provide links to the corresponding primary publications.

## Repository URL

<https://github.com/cschen-y/rna-small-molecule-review>
