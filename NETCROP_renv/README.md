# NETCROP: Network Cross-Validation with Overlapping Partitions

This repository contains the code used for the paper:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2025).
> “Network Cross-Validation and Model Selection via Subsampling.”
> arXiv:2504.06903.

## Start here

### 1. Install the prerequisites and clone the repo

Install:

- R 4.6 (the version recorded in `renv.lock`, however the codes should run in R 4.4 or newer)
- RStudio, recommended
- A C++ compiler:
  - macOS: Xcode Command Line Tools (`xcode-select --install` in Terminal)
  - Windows: the Rtools version matching your R installation
  - Linux: `build-essential`, or the equivalent compiler toolchain
  
Clone:

- Open a terminal in a target directory and run the following:

```bash
git clone https://github.com/sayan-ch/NETCROP.git
```

### 2. Open the project correctly

Open `NETCROP.Rproj` in RStudio. Do not open an individual R script first.

In the R console, verify the working directory:

```r
getwd()
file.exists("NETCROP.Rproj")
file.exists("renv.lock")
```

Both `file.exists()` calls must return `TRUE`. If either returns `FALSE`, close
RStudio and open `NETCROP.Rproj` directly.

### 3. Restore the package environment

Run exactly this command in the R console:

```r
source("first_time_renv_setup.R")
```

The setup script:

1. checks that R is at the repository root;
2. installs `renv` from CRAN if necessary;
3. restores the exact packages recorded in `renv.lock`;
4. reports dependency inconsistencies; and
5. compiles and loads the general C++ helpers as a smoke test.

If package installation is requested, answer yes. A successful run ends with:

```text
NETCROP environment restored and helper smoke test passed.
```

You normally run this setup only once per fresh clone or after deleting the
project library.

### 4. Run a small test first

The small tests use `n = 500`, `K` or `d = 3`, a maximum candidate value of 5,
`p.test = 0.1`, and two simulations.

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/xx_small_network_test.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table2/xx_small_network_test.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table3/xx_small_network_test.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table4_realdata/xx_simulated_small_network_test.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Figure2/xx_small_network_test.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/FigureS3/xx_quicker_test.R"))
```

Run one command at a time. The Table 4 small script is a synthetic DCBM smoke
test for the same model-selection routines; it does not reproduce a real-data
result.

## Running the paper cases

Full cases can be expensive. Networks with `n = 10000` may require substantial
RAM and can run for hours depending on hardware and the selected method.

### Table 1: SBM and DCBM

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case1_sbmK5.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case2_sbmK20.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case3_dcbmK10.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case4_dcbmK20.R"))
```

### Table 2: RDPG

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table2/case1_zeta0_75.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table2/case2_zeta0_70.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table2/case3_zeta0_65.R"))
```

### Table 3: latent space models

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table3/case1_d2_a0.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table3/case2_d2_a1.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table3/case3_d5_a0.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table3/case4_d5_a1.R"))
```

### Table 4: real networks

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table4_realdata/DBLP_analysis.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table4_realdata/Twitch_analysis.R"))
```

### Figures

```r
source(file.path(here::here(), "NETCROP_PAPER_CODES/Figure2/Figure2_partune_rsc.R"))
source(file.path(here::here(), "NETCROP_PAPER_CODES/FigureS3/FigureS3_small_networks.R"))
```

## Existing results: replace, resume, or archive

Every simulation script has its own directories:

```text
<table-or-figure>/output/<script-name>/
<table-or-figure>/logs/<script-name>/
```

When an interactive R session finds existing results, it asks whether to:

1. **Replace** — delete that script's existing results and restart at
   simulation 1.
2. **Resume** — remove duplicate or incomplete trailing rows and continue from
   the first incomplete simulation. NETCROP, NCV, and ECV resume independently.
3. **Archive** — rename the script's output and log directories with a timestamp
   and start new directories.

Cancelling the menu stops the simulation without guessing.

### Batch jobs

Batch jobs cannot answer an interactive menu. Set the environment variable
`NETCROP_OUTPUT_ACTION` to `replace`, `resume`, or `archive`.

If it is not set, batch mode defaults to `replace` and restarts from simulation
1.

From a shell:

```sh
NETCROP_OUTPUT_ACTION=resume Rscript -e 'source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case1_sbmK5.R")'
```

From R before sourcing a script:

```r
Sys.setenv(NETCROP_OUTPUT_ACTION = "resume")
source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case1_sbmK5.R"))
```

To return to interactive menus in the same R session:

```r
Sys.unsetenv("NETCROP_OUTPUT_ACTION")
```

## NCV and ECV on large networks

Before every NCV or ECV loop with `n > 1000`, an interactive session warns that
the method can consume substantial time and memory. If you decline, the method
is skipped and the console displays a runnable, clickable RStudio command for
the relevant `xx_small_network_test.R` script.

For batch jobs, large NCV and ECV runs are skipped unless explicitly enabled:

```sh
NETCROP_RUN_LARGE_CV=yes NETCROP_OUTPUT_ACTION=resume \
Rscript -e 'source(file.path(here::here(), "NETCROP_PAPER_CODES/Table1/case1_sbmK5.R")'
```

Accepted true values are `yes`, `true`, and `1`; accepted false values are
`no`, `false`, and `0`.

## Repository structure

```text
NETCROP/
├── .Rprofile
├── .gitignore
├── NETCROP.Rproj
├── README.md
├── first_time_renv_setup.R
├── renv.lock
├── renv/
│   ├── activate.R
│   └── settings.json
└── NETCROP_PAPER_CODES/
    ├── Table1/
    │   ├── case1_sbmK5.R
    │   ├── case2_sbmK20.R
    │   ├── case3_dcbmK10.R
    │   ├── case4_dcbmK20.R
    │   └── xx_small_network_test.R
    ├── Table2/
    │   ├── case1_zeta0_75.R
    │   ├── case2_zeta0_70.R
    │   ├── case3_zeta0_65.R
    │   └── xx_small_network_test.R
    ├── Table3/
    │   ├── case1_d2_a0.R
    │   ├── case2_d2_a1.R
    │   ├── case3_d5_a0.R
    │   ├── case4_d5_a1.R
    │   └── xx_small_network_test.R
    ├── Table4_realdata/
    │   ├── DBLP/
    │   ├── Twitch/
    │   ├── DBLP_analysis.R
    │   ├── Twitch_analysis.R
    │   └── xx_simulated_small_network_test.R
    ├── Figure2/
    │   ├── Figure2_partune_rsc.R
    │   └── xx_small_network_test.R
    ├── FigureS3/
    │   └── FigureS3_small_networks.R
    │   └── xx_quicker_test.R
    └── helpers/
        ├── General_helpers.R
        ├── General_helpers.cpp
        ├── SBM_DCBM_helpers.R
        ├── RDPG_helpers.R
        ├── LSM_helpers.R
        ├── LSM_PGD_Cpp.cpp
        └── PARTUNE_RSC_helpers.R
```

Generated `output/`, `logs/`, local package libraries, RStudio state, and macOS
metadata are ignored by Git.

## Main functions

- `netcrop_blockmodel()` — NETCROP for SBM and DCBM.
- `netcrop_rdpg()` — NETCROP for random dot product graphs.
- `netcrop_lsm()` — NETCROP for latent space models.
- `netcrop.tune.regsp()` — NETCROP tuning for regularized spectral clustering.
- `NCV.stability.BM()` and `ECV.stability.BM()` — block-model comparison
  procedures.
- `ECV.stability.RDPG()` — RDPG comparison procedure.

General AUC, outer-addition, universal singular-value thresholding, and
Procrustes operations are implemented locally in `General_helpers.R` and
`General_helpers.cpp`. The external packages previously used for those
operations are not required.

## Troubleshooting

### `source()` says a file does not exist

You are probably not at the repository root. Reopen `NETCROP.Rproj` and run:

```r
getwd()
list.files()
```

### Package restore fails

Check your internet connection, restart R, and run:

```r
source("first_time_renv_setup.R")
```

Do not install individual project packages manually unless the restore error
specifically instructs you to do so.

### C++ compilation fails

Confirm that the platform compiler listed under “Install the prerequisites” is
installed. Then restart R and rerun `first_time_renv_setup.R`.

### A process runs out of memory

Stop the full script and run its `xx_small_network_test.R` file. For large NCV
or ECV runs, close other memory-intensive applications and use a machine with
adequate RAM.

### Check the restored environment

```r
renv::status()
```

The project library is intentionally not committed. A fresh clone reconstructs
it from `renv.lock` by running `source("first_time_renv_setup.R")`.

## Citation

```bibtex
@misc{chakrabarty2025network,
  title = {Network Cross-Validation and Model Selection via Subsampling},
  author = {Chakrabarty, Sayan and Sengupta, Srijan and Chen, Yuguo},
  year = {2025},
  eprint = {2504.06903},
  archivePrefix = {arXiv},
  primaryClass = {stat.ME},
  url = {https://arxiv.org/abs/2504.06903}
}
```
