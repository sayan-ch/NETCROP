# NETCROP: Network Cross-Validation with Overlapping Partitions

This repository contains the implementations of NETCROP, a method for
network cross-validation with overlapping partitions. The codes in the
"NETCROP_PAPER_CODES" can be used to replicate the numerical results in
the paper "Network Cross-Validation and Model Selection via Subsampling"
[1] by Sayan Chakrabarty, Srijan Sengupta and Yuguo Chen
(<https://arxiv.org/abs/2504.06903>).

## Organization

The repository is organized as follows:

### `NETCROP_PAPER_CODES/`:

#### `Table1/`:

-   `case1_sbmK5.R`: SBM with $n = 10000$ and $K = 5$.

-   `case2_sbmK20.R`: SBM with $n = 10000$ and $K = 20$.

-   `case3_dcbmK10.R`: DCBM with $n = 10000$ and $K = 10$.

-   `case4_dcbmK20.R`: DCBM with $n = 10000$ and $K = 20$.

#### `Table2/`:

-   `case1_zeta0_75.R`: RDPG with $n = 10000$, $d = 5$ and
    $\zeta = 0.75$.

-   `case2_zeta0_70.R`: RDPG with $n = 10000$, $d = 5$ and
    $\zeta = 0.70$.

-   `case3_zeta0_65.R`: RDPG with $n = 10000$, $d = 5$ and
    $\zeta = 0.65$.

#### `Table3/`:

-   `case1_d2_a0.R`: LSM with $n = 1000$, $d = 2$ and $\alpha = 0$.

-   `case2_d2_a1.R`: LSM with $n = 1000$, $d = 2$ and $\alpha = 1$.

-   `case3_d5_a0.R`: LSM with $n = 1000$, $d = 2$ and $\alpha = 2$.

-   `case4_d5_a1.R`: LSM with $n = 1000$, $d = 2$ and $\alpha = 3$.

#### `Table4_realdata`:

-   `DBLP/`: Contains DBLP 4-area network in .csv and .rds formats.

-   `DBLP_analysis.R`: NETCROP, NCV and ECV on DBLP author-conference
    network.

-   `Twitch/`: Contains Twitch social network in .csv and .rds formats.
    It also contains the full Twitch network and the cleanup codes to
    obtain the exact subnetwork used in the paper.

-   `Twitch_analysis.R`: NETCROP on Twitch social network.

#### `Figure2/`:

-   `Figure2_partune_rsc.R`: Parameter tuning for regularized spectral
    clustering on DCBM with $n = 10000$ and $K = 5$.

#### `FigureS3/`:

-   `FigureS3_small_networks.R`: NETCROP, NCV and ECV on small networks
    with $n \in \{200, 300, 400, 500\}$ for SBM and DCBM with $K = 3$.

#### `helpers/`:

-   `General_helpers.R`: Contains general helper functions used in all
    examples, including parameter selection for NETCROP.

-   `SBM_DCBM_helper.R`: Contains the functions

    -   **`netcrop_blockmodel`: Function to perform NETCROP for block
        models (SBM and DCBM).**

    -   `SBM.gen`, `DCBM.gen`: Generates SBM and DCBM

    -   `best.perm.label.match`: MatchGreedy algorithm for label
        matching

    -   `fast.SBM.est`, `NCV.SBM.est`: Functions to estimate SBM
        parameters using the full adjacency matrix and the same for NCV
        [2] that uses rectangular submatrices, respectively.

    -   `fast.DCBM.est`, `NCV.DCBM.est`: Functions to estimate DCBM
        parameters by profile likelihood using the full adjacency matrix
        and the same for NCV that uses rectangular submatrices,
        respectively.

    -   `eigen.DCBM.est`, `NCV.eigen.DCBM.est`: Functions to estimate
        DCBM parameters by spectral method using the full adjacency
        matrix and the same for NCV that uses rectangular submatrices,
        respectively, using spectral clustering.

    -   `ECV.stability.BM`: Wrapper for ECV [3] for blockmodels (taken
        from `randnet` package [4]) adding the choice of loss function
        and stability (repetition).

    -   `NCV.stability.BM`: Wrapper for NCV for blockmodels (taken from
        `randnet` package) adding the choice of loss function and
        stability (repetition).

-   `RDPG_helpers.R`: Contains the functions

    -   **`netcrop_rdpg`: Function to perform NETCROP for RDPG.**

    -   `RDPG.gen`: Generates RDPG

    -   `ECV.stability.RDPG`: Wrapper for ECV for RDPG (taken from
        `randnet` package) adding the choice of loss function and
        stability (repetition).

-   `LSM_helpers.R`: Contains the functions

    -   **`netcrop_lsm`: Function to perform NETCROP for LSM.**

    -   `LSM.gen`: Generates latent space model [5, 6]

    -   `pgd.LSM`: Function to estimate LSM parameters by projected
        gradient descent (PGD, [6]) using the full adjacency matrix.
        Uses Rcpp for speed and is in "LSM_PGD_Cpp.cpp".

    -   `ECV.stability.LSM`: Wrapper for ECV for LSM (taken from
        `randnet` package) adding the choice of loss function and
        stability (repetition).

-   `PARTUNE_RSC_helpers.R`: Contains the functions

    -   **`netcrop_tune_regsp`: Function to perform NETCROP for
        selecting the best tuning parameter for regularized spectral
        clusrtering.**

    -   `DKest`: Davis-Kahan estimator for the tuning parameter [7].

## References:

[1] Chakrabarty, S., Sengupta, S., and Chen, Y. (2025), "Network
Cross-Validation and Model Selection via Subsampling". arXiv preprint
arXiv:2504.06903.

[2] Chen, K. and Lei, J. (2018), “Network Cross-Validation for
Determining the Number of Communities in Network Data,” Journal of the
American Statistical Association, 113, 241–251.

[3] Li, T., Levina, E., and Zhu, J. (2020), “Network Cross-Validation by
Edge Sampling,” Biometrika, 107, 257–276.

[4] Li, T., Levina, E., Zhu, J., and Le, C. M. (2023), “randnet: Random
Network Model Estimation, Selection and Parameter Tuning,” CRAN, R
Package Version 0.7, <https://CRAN.R-project.org/package=randnet>.

[5] Hoff, P. D., Raftery, A. E., and Handcock, M. S. (2002), “Latent
Space Approaches to Social Network Analysis,” Journal of the American
Statistical Association, 97, 1090–1098.

[6] Ma, Z., Ma, Z., and Yuan, H. (2020), “Universal Latent Space Model
Fitting for Large Networks with Edge Covariates,” Journal of Machine
Learning Research, 21, 1–67.

[7] Joseph, A. and Yu, B. (2016), “Impact of Regularization on Spectral
Clustering,” The Annals of Statistics, 44, 1765–1791.
