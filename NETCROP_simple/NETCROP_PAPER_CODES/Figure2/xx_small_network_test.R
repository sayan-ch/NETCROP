#!/usr/bin/env Rscript

# Fast Figure 2 algorithm check using only netOP and its imported packages.

suppressPackageStartupMessages(library(netOP))
if (packageVersion("netOP") != "0.1.1") {
  stop("This script requires netOP v0.1.1.")
}

script_dir <- function() {
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(...) NULL)
  if (!is.null(source_file)) return(dirname(normalizePath(source_file)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
  }
  normalizePath(getwd())
}

root <- script_dir()
output_dir <- file.path(root, "output", "xx_small_network_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ncores <- 1L
nsim <- 10L
n <- 500L
K <- 3L
params <- netcrop_param_select(test_prop = 0.1, n = n, o_range = 0)

one_simulation <- function(simulation) {
  psi_pool <- stats::rbeta(300L, 1, 4)
  psi <- sample(psi_pool, n, replace = TRUE)
  A <- generate_dcbm(
    n = n, K = K, alpha = 0.3, beta = 1 / 3, psi = psi,
    degree_scale = "max_by_community", representation = "sparse",
    ncores = ncores
  )
  netcrop_fit <- netcrop_tune_regularizer(
    A, K = K, tau_candidates = seq(0, 2, 0.1), use_dcbm = TRUE,
    num_subnetworks = params$num_subnetworks[[1L]],
    overlap_size = params$overlap_size[[1L]], nrep = 5L,
    use_laplacian = TRUE, dcbm_est_method = "plugin", losses = "sse",
    ncores = ncores, retain_intermediates = "minimal"
  )
  dk_fit <- dkest_tune_regularizer(
    A, K = K, tau_candidates = seq(0, 0.1, 0.01), use_laplacian = TRUE,
    use_dcbm = TRUE, dcbm_est_method = "plugin", ncores = ncores,
    retain_intermediates = "minimal"
  )
  cat(sprintf("[%d/%d] NETCROP tau %.3f; DKEST tau %.3f\n",
              simulation, nsim,
              netcrop_fit$overall_best$tau_hat[[1L]], dk_fit$tau_hat))
  print(netcrop_fit$overall_best, row.names = FALSE)
  print(dk_fit$overall_best, row.names = FALSE)
  rows <- list(data.frame(
      simulation = simulation,
      netcrop_tau = netcrop_fit$overall_best$tau_hat[[1L]],
      dkest_tau = dk_fit$tau_hat,
      netcrop_seconds = netcrop_fit$timing[["total"]],
      dkest_seconds = dk_fit$timing[["total"]]
    ))
  summary <- do.call(rbind, rows)
  print(summary, row.names = FALSE)
  list(summary = summary, netcrop = netcrop_fit, dkest = dk_fit)
}

records <- run_simulations(
  one_simulation = one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "small_figure2_results.rds"),
  action = "archive", show_progress = TRUE
)
stopifnot(all(vapply(records, function(x) isTRUE(x$success), logical(1))))
summary_table <- do.call(rbind, lapply(records, function(x) x$result$summary))
csv_file <- file.path(output_dir, "small_figure2_results.csv")
if (file.exists(csv_file)) {
  archived_csv <- file.path(
    output_dir,
    paste0("small_figure2_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  )
  file.rename(csv_file, archived_csv)
}
utils::write.csv(summary_table, csv_file, row.names = FALSE)
cat("Small Figure 2 test completed:", normalizePath(output_dir), "\n")
