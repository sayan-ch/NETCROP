#!/usr/bin/env Rscript

# Fast synthetic check for the three block-model selectors used in Table 4.
# It intentionally uses only netOP and packages imported by netOP.

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
output_dir <- file.path(root, "output", "xx_simulated_small_network_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ncores <- 1L
nsim <- 10L
n <- 500L
K <- 3L
max_K <- 5L
params <- netcrop_param_select(test_prop = 0.1, n = n, o_range = 0)

one_simulation <- function(simulation) {
  psi_pool <- 1 / stats::rbeta(300L, 4, 1)
  psi <- sample(psi_pool, n, replace = TRUE)
  A <- generate_dcbm(
    n = n, K = K, alpha = 0.5, beta = 0.3, psi = psi,
    degree_scale = "none", average_degree = 30,
    representation = "sparse", ncores = ncores
  )

  fits <- list()
  for (nrep in c(1L, 5L)) {
    fits[[paste0("NETCROP_R", nrep)]] <- netcrop_blockmodel(
        A, K_candidates = seq_len(max_K),
        num_subnetworks = params$num_subnetworks[[1L]],
        overlap_size = params$overlap_size[[1L]], nrep = nrep,
        losses = c("sse", "auc_as_loss"), ncores = ncores,
        dcbm_est_options = list(estimate = list(method = "plugin")),
        retain_intermediates = "minimal"
      )
  }
  for (nrep in c(1L, 20L)) {
    fits[[paste0("NCV_R", nrep)]] <- ncv_stability_blockmodel(
        A, max_K = max_K, cv = 3L, nrep = nrep, dc_est = "plugin",
        losses = c("sse", "auc_as_loss"), ncores = ncores,
        retain_intermediates = "minimal"
      )
    fits[[paste0("ECV_R", nrep)]] <- ecv_stability_blockmodel(
        A, max_K = max_K, cv = 3L, nrep = nrep,
        losses = c("sse", "auc_as_loss"), ncores = ncores,
        retain_intermediates = "minimal"
      )
  }

  rows <- do.call(rbind, lapply(names(fits), function(algorithm) {
    fit <- fits[[algorithm]]
    best <- fit$overall_best$best_model[fit$overall_best$loss == "sse"]
    cat(sprintf("[%d/%d] %s selected %s\n",
                simulation, nsim, algorithm, best))
    print(fit$overall_best, row.names = FALSE)
    name_parts <- strsplit(algorithm, "_R", fixed = TRUE)[[1L]]
    data.frame(simulation = simulation, algorithm = name_parts[[1L]],
               nrep = as.integer(name_parts[[2L]]), best_model = best,
               elapsed_seconds = fit$timing[["total"]])
  }))
  list(summary = rows, fits = fits)
}

records <- run_simulations(
  one_simulation = one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "small_table4_results.rds"),
  action = "archive", show_progress = TRUE
)
stopifnot(all(vapply(records, function(x) isTRUE(x$success), logical(1))))
summary_table <- do.call(rbind, lapply(records, function(x) x$result$summary))
csv_file <- file.path(output_dir, "small_table4_results.csv")
if (file.exists(csv_file)) {
  archived_csv <- file.path(
    output_dir,
    paste0("small_table4_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  )
  file.rename(csv_file, archived_csv)
}
utils::write.csv(summary_table, csv_file, row.names = FALSE)
cat("Small Table 4 test completed:", normalizePath(output_dir), "\n")
