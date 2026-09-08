#!/usr/bin/env Rscript

# Figure 2: NETCROP and Davis--Kahan regularizer selection for DCBM spectral
# clustering. Requires netOP 0.1.1. Random seeds are intentionally omitted.

suppressPackageStartupMessages(library(netOP))
if (packageVersion("netOP") != "0.1.1") stop("This script requires netOP v0.1.1.")

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
output_dir <- file.path(root, "output", "Figure2_partune_rsc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
detected <- parallel::detectCores()
ncores <- if (is.na(detected)) 1L else max(1L, floor(detected / 2))
nsim <- 100L
n <- 10000L
K <- 5L
action <- Sys.getenv("NETCROP_ACTION", "resume")
if (!action %in% c("replace", "resume", "archive")) {
  stop("NETCROP_ACTION must be replace, resume, or archive.")
}

one_simulation <- function(simulation) {
  psi_pool <- stats::rbeta(300L, 1, 4)
  psi <- sample(psi_pool, n, replace = TRUE)
  A <- generate_dcbm(
    n = n, K = K, alpha = 0.3, beta = 1 / 3, psi = psi,
    degree_scale = "max_by_community", representation = "sparse",
    ncores = ncores
  )
  truth <- get_generator_parameters(A)$g_true
  average_degree <- mean(Matrix::rowSums(A))

  netcrop_seconds <- system.time({
    netcrop_fit <- netcrop_tune_regularizer(
      A = A, K = K, tau_candidates = seq(0, 2, 0.1),
      use_dcbm = TRUE, num_subnetworks = 3L, overlap_size = 8002L,
      nrep = 5L, use_laplacian = TRUE, dcbm_est_method = "plugin",
      losses = "sse", ncores = ncores, retain_intermediates = "minimal"
    )
  })[["elapsed"]]
  dk_seconds <- system.time({
    dk_fit <- dkest_tune_regularizer(
      A = A, K = K, tau_candidates = seq(0, 0.1, 0.01),
      use_laplacian = TRUE, use_dcbm = TRUE,
      dcbm_est_method = "plugin", ncores = ncores,
      retain_intermediates = "minimal"
    )
  })[["elapsed"]]

  cat(sprintf(
    "[%d/%d] average degree %.2f; NETCROP tau %.3f (%.2fs); DKEST tau %.3f (%.2fs)\n",
    simulation, nsim, average_degree,
    netcrop_fit$overall_best$tau_hat[[1L]], netcrop_seconds,
    dk_fit$tau_hat, dk_seconds
  ))
  print(netcrop_fit$overall_best, row.names = FALSE)
  print(dk_fit$overall_best, row.names = FALSE)
  rows <- list(data.frame(
      simulation = simulation, n = n, K = K,
      average_degree = average_degree,
      netcrop_tau = netcrop_fit$overall_best$tau_hat[[1L]],
      dkest_tau = dk_fit$tau_hat,
      netcrop_seconds = netcrop_seconds, dkest_seconds = dk_seconds
    ))
  summary <- do.call(rbind, rows)
  print(summary, row.names = FALSE)
  list(A = A, truth = truth, netcrop = netcrop_fit, dkest = dk_fit,
       summary = summary)
}

records <- run_simulations(
  one_simulation = one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "Figure2_results.rds"),
  action = action, show_progress = TRUE
)
successful <- records[vapply(records, function(x) isTRUE(x$success), logical(1))]
if (!length(successful)) stop("No successful Figure 2 simulations are available.")
values <- lapply(successful, `[[`, "result")
summary_table <- do.call(rbind, lapply(values, `[[`, "summary"))
utils::write.csv(summary_table, file.path(output_dir, "Figure2_results.csv"),
                 row.names = FALSE)

# netOP's plotting helper recreates the oracle, unregularized, DKEST, and
# NETCROP-selected clustering-accuracy comparison from the fitted objects.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Figure rendering uses netOP::oracle_plotter() and requires ggplot2.")
}
figure <- oracle_plotter(
  A = lapply(values, `[[`, "A"),
  g_true = lapply(values, `[[`, "truth"), K = K,
  netcrop_outcomes = lapply(values, `[[`, "netcrop"),
  dkest_outcomes = lapply(values, `[[`, "dkest"),
  losses = "sse", engines = "spectral_cluster", ncores = ncores
)
ggplot2::ggsave(file.path(output_dir, "Figure2.pdf"), figure,
                width = 8, height = 4.5, dpi = 600, device = "pdf")
accuracy_data <- attr(figure, "accuracy_data")
plot_data <- attr(figure, "plot_data")
utils::write.csv(accuracy_data,
                 file.path(output_dir, "Figure2_accuracy_by_simulation.csv"),
                 row.names = FALSE)
utils::write.csv(plot_data, file.path(output_dir, "Figure2_plot_data.csv"),
                 row.names = FALSE)
saveRDS(list(plot = figure, accuracy_data = accuracy_data,
             plot_data = plot_data, metadata = attr(figure, "metadata")),
        file.path(output_dir, "Figure2_plot.rds"))
print(plot_data, row.names = FALSE)
