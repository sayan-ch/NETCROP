# NETCROP simple replication with netOP

This version reproduces the paper scenarios with the released `netOP` R package. It has no `.Rproj`, `.Rprofile`, or `renv` environment, and it does not source the original helper implementations.

Use this version for quick verification. For the original seeded environment intended for exact replication, see [../NETCROP_renv/README.md](../NETCROP_renv/README.md).

## What is and is not reproduced

The scripts preserve the network models, parameter grids, candidate models or dimensions, repetitions, comparisons, summaries, and figures from the paper code. They deliberately omit every seed argument. Consequently, reruns generate new networks and will not reproduce the paper's exact random realizations or numerical values.

Simulation results are saved twice:

- a package-native RDS file containing the `netOP::run_simulations()` records;
- a derived CSV file suitable for inspection and plotting.

The scripts print simulation progress, method names, intermediate selections, elapsed time, and failures in a style similar to the original scripts.

## Download the `NETCROP_simple` directory

Run the following commands in a terminal. Replace `<HOME_DIR>` with a directory of your choice.

```bash
cd <HOME_DIR>
git clone --filter=blob:none --sparse https://github.com/sayan-ch/NETCROP.git
cd NETCROP
git sparse-checkout set NETCROP_simple
```

## Install netOP

Start an R session and run:

```{r}
install.packages("netOP")
```

## Verify the installation

Start a fresh R session and run:

```r
library(netOP)
packageVersion("netOP")
stopifnot(packageVersion("netOP") == "0.1.1")

A <- generate_sbm(
  n = 60,
  K = 3,
  alpha = 0.35,
  beta = 0.08,
  representation = "sparse",
  ncores = 1
)
stopifnot(identical(dim(A), c(60L, 60L)))
```

The version check must report `0.1.1` and the final command must finish without an error.

## Run the quick tests first

Open R or RStudio and change the working directory to `NETCROP_simple`.

```r
setwd("<HOME_DIR>/NETCROP/NETCROP_simple")
getwd()
file.exists("README.md")
file.exists("NETCROP_PAPER_CODES")
```

Then run one quick script at a time. Alternatively, open the .R files manually in RStudio and click the Source button.

```r
source("NETCROP_PAPER_CODES/Table1/xx_small_network_test.R")
source("NETCROP_PAPER_CODES/Table2/xx_small_network_test.R")
source("NETCROP_PAPER_CODES/Table3/xx_small_network_test.R")
source("NETCROP_PAPER_CODES/Table4_realdata/xx_simulated_small_network_test.R")
source("NETCROP_PAPER_CODES/Figure2/xx_small_network_test.R")
source("NETCROP_PAPER_CODES/FigureS3/xx_quicker_test.R")
```

All six use 10 simulations. The first five use `n = 500`; the Figure S3 quick test uses `n = 100` and `n = 200`, matching its original quick scenario.

## Run the paper cases

### Table 1: SBM and DCBM

```r
source("NETCROP_PAPER_CODES/Table1/case1_sbmK5.R")
source("NETCROP_PAPER_CODES/Table1/case2_sbmK20.R")
source("NETCROP_PAPER_CODES/Table1/case3_dcbmK10.R")
source("NETCROP_PAPER_CODES/Table1/case4_dcbmK20.R")
```

### Table 2: RDPG

```r
source("NETCROP_PAPER_CODES/Table2/case1_zeta0_75.R")
source("NETCROP_PAPER_CODES/Table2/case2_zeta0_70.R")
source("NETCROP_PAPER_CODES/Table2/case3_zeta0_65.R")
```

### Table 3: latent-space models

```r
source("NETCROP_PAPER_CODES/Table3/case1_d2_a0.R")
source("NETCROP_PAPER_CODES/Table3/case2_d2_a1.R")
source("NETCROP_PAPER_CODES/Table3/case3_d5_a0.R")
source("NETCROP_PAPER_CODES/Table3/case4_d5_a1.R")
```

### Table 4: real networks

```r
source("NETCROP_PAPER_CODES/Table4_realdata/DBLP_analysis.R")
source("NETCROP_PAPER_CODES/Table4_realdata/Twitch_analysis.R")
```

### Figures

```r
source("NETCROP_PAPER_CODES/Figure2/Figure2_partune_rsc.R")
source("NETCROP_PAPER_CODES/FigureS3/FigureS3_small_networks.R")
```

Full cases can take hours and require substantial RAM, especially the `n = 10000` cases and stabilized NCV/ECV comparisons.

## Results, replacement, and resuming

Each script keeps its generated files beneath its own `output/<script-name>/` directory. The RDS file is the authoritative netOP simulation record; the CSV is regenerated from successful and failed records for convenient review.

At the start of a run, choose one of the actions supported by `run_simulations()`:

- `replace`: start again and replace the active result file;
- `resume`: retain completed records and retry missing or failed simulations;
- `archive`: timestamp the existing result file and begin a new run.

Because `run_simulations()` writes its durable RDS after the current simulation call returns, an interrupted R process can lose work from that in-progress call. This differs from the original scripts' incremental CSV workflow.

## Cores and progress

Most scripts use one outer simulation loop and allow netOP methods to use a configurable number of cores. Figure S3 intentionally follows the legacy layout: simulations run over half of the detected cores while each individual NETCROP, NCV, or ECV call uses one core.

Every script prints the current simulation, method, repetition count, selected model/dimension/regularizer, elapsed time, and recorded errors. Since the scripts are unseeded, repeated runs can legitimately make different selections.

For the full Table 1, Table 2, and DBLP workflows, all NETCROP simulations run, save, and summarize first. The scripts then warn that NCV/ECV can be slow, memory intensive, and capable of crashing a personal R session. An interactive R session asks for confirmation; `Rscript` prints the warning and proceeds automatically. Quick tests print the same warning and proceed without prompting.

Each simulation prints its bound result rows immediately before returning. Final simulated-data summaries are grouped by true model, algorithm, and repetition count (`R`) and report the modal `best_model`, its frequency, accuracy, and mean absolute deviation of the selected K or dimension from the truth.

## Troubleshooting

- **Wrong netOP version:** restart R after installing, run `find.package("netOP")`, and confirm `packageVersion("netOP") == "0.1.1"`.
- **Compilation failure:** install the matching compiler toolchain, restart R, and retry the pinned source installation.
- **Binary cannot be installed:** confirm that the asset matches both the R major/minor series and CPU architecture; otherwise install from source.
- **Files cannot be found:** set the working directory to `NETCROP_simple`, not to an individual Table/Figure directory.
- **Out of memory:** run a quick test first, use one core, close other large applications, and avoid retaining large intermediates.
- **Slow NCV/ECV:** these competing methods are expected to be much slower on large networks. Validate the quick cases before starting a full experiment.
- **Failed simulation records:** inspect the `error` field in the RDS-derived CSV, fix the underlying problem, and rerun with the resume action.

## Directory tree

```text
NETCROP_simple/
├── .gitignore
├── README.md
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
    └── FigureS3/
        ├── FigureS3_small_networks.R
        └── xx_quicker_test.R
```
