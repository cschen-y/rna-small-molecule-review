# Reproducible environments

The nine Test17 methods require three environments because their Python, Biopython, PyTorch and CUDA constraints are incompatible.

| Environment | Methods | Python | PyTorch | CUDA | Key constraint |
|---|---|---:|---:|---:|---|
| `core` | RLBind, RNet, RBind, Rsite, Rsite2, RNAsite | 3.9.21 | 2.1.0 | 11.8 | Compatible with the six legacy/core methods |
| `MVRBind` | MVRBind, RNABind | 3.10.18 | 2.2.0 | 12.1 | PyG 2.6.1 and CUDA 12.1 extension wheels |
| `Mul` | MultiModRLBP | 3.7.12 | 1.13.1 | 11.6 | `Bio.Alphabet` requires Biopython 1.77; bundled native module targets CPython 3.7 |

Create all environments from the repository root:

```bash
conda env create -f environments/core.yml
conda env create -f environments/mvrbind.yml
conda env create -f environments/multimodrlbp.yml
```

MultiModRLBP includes a Windows CPython 3.7 binary for `alignment_C`. On Linux, or when the binary is incompatible, rebuild it after activating `Mul`:

```bash
conda activate Mul
cd MultiModRLBP/RnaBert
python setup.py build_ext --inplace
cd ../..
```

Run every method from the repository root:

```bash
python run_all_test17.py --methods all --continue-on-error
```

The runner uses `conda run` and the three environment names above. Set `CONDA_EXE` or pass `--conda /path/to/conda` when Conda is not on `PATH`.

The pinned files describe the NVIDIA GPU configuration used for the recorded results. A compatible NVIDIA driver is required. CPU-only reproduction requires replacing the PyTorch/CUDA and compiled PyG/DGL packages with matching CPU builds and will be slower.
