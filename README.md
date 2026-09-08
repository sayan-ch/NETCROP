# NETCROP replication code

This repository provides two versions of the code for:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2026).
> “Network Cross-Validation and Model Selection via Subsampling.”
> arXiv:2504.06903.

## Choose a version

### [NETCROP_simple](NETCROP_simple/README.md) — recommended for quick verification

This version installs the released [netOP](https://github.com/sayan-ch/netOP) R package and uses its exported network generators, estimators, NETCROP, NCV, ECV, simulation, and plotting functions. It does not use an RStudio project or `renv`.

The scripts deliberately omit seeds. They preserve the paper's scenarios and analysis workflow, but their random networks and numerical results will not be bit-for-bit identical to the paper because netOP handles randomness differently from the original helper code.

Start here if you want to verify the methods quickly or read the current, package-based implementation.

### [NETCROP_renv](NETCROP_renv/README.md) — exact seeded replication environment

This is the original seeded paper-code version, preserved unchanged. It uses an RStudio project, a locked `renv` environment, and locally compiled R/C++ helpers. Restoring the environment and compiling its dependencies can take substantially longer, but this is the version intended for exact replication of the archived workflow.

## Download the repository

### Full clone

```bash
git clone https://github.com/sayan-ch/NETCROP.git
cd NETCROP
```

### Pull only the simple version

Git cannot clone a subdirectory as an independent repository, but sparse checkout downloads and checks out only the selected folder plus the root documentation:

```bash
git clone --filter=blob:none --sparse https://github.com/sayan-ch/NETCROP.git
cd NETCROP
git sparse-checkout set NETCROP_simple
```

Then follow [NETCROP_simple/README.md](NETCROP_simple/README.md).

### Pull only the renv version

```bash
git clone --filter=blob:none --sparse https://github.com/sayan-ch/NETCROP.git
cd NETCROP
git sparse-checkout set NETCROP_renv
```

Then follow [NETCROP_renv/README.md](NETCROP_renv/README.md).

To add the other version later:

```bash
git sparse-checkout set NETCROP_simple NETCROP_renv
```

## Before running full experiments

Run the quick tests first. Full experiments include networks with 10,000 nodes, repeated model fits, and competing cross-validation procedures; they can require substantial memory and hours of computation.

The supplementary small-network benchmark is named **Figure S3** throughout the new simple version's code and filenames. The frozen renv snapshot is preserved byte-for-byte, including any historical output-name spelling.
