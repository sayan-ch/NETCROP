#!/usr/bin/env Rscript

# Table 4: block-model selection on the DBLP conference network.
# Requires netOP 0.1.1. Results are stored as native R objects and as a
# convenient flat CSV table. Random seeds are intentionally not supplied.

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
output_dir <- file.path(root, "output", "DBLP_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

detected <- parallel::detectCores()
ncores <- if (is.na(detected)) 1L else max(1L, floor(detected / 2))
nsim <- 100L
action <- Sys.getenv("NETCROP_ACTION", "resume")
if (!action %in% c("replace", "resume", "archive")) {
  stop("NETCROP_ACTION must be replace, resume, or archive.")
}

edge_list <- utils::read.csv(file.path(root, "DBLP", "DBLP_conf_edge_list.csv"))
n <- max(edge_list$i, edge_list$j)
A <- Matrix::sparseMatrix(
  i = edge_list$i, j = edge_list$j, x = 1,
  dims = c(n, n)
)
params <- netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)
num_subnetworks <- params$num_subnetworks[[1L]]
overlap_size <- params$overlap_size[[1L]]

chosen_auc <- function(fit, selected_model) {
  bits <- strsplit(selected_model, "-", fixed = TRUE)[[1L]]
  selected <- fit$cv_all_loss[
    fit$cv_all_loss$loss == "auc_as_loss" &
      fit$cv_all_loss$model == bits[[1L]] &
      fit$cv_all_loss$K == as.integer(bits[[2L]]),
    "loss_value"
  ]
  1 - mean(selected, na.rm = TRUE)
}

one_simulation <- function(simulation) {
  rows <- list()
  fits <- list(NETCROP = list(), NCV = list(), ECV = list())
  row_id <- 0L

  for (nrep in c(1L, 5L, 10L, 20L)) {
    elapsed <- system.time({
      fit <- netcrop_blockmodel(
        A = A,
        K_candidates = 1:10,
        num_subnetworks = num_subnetworks,
        overlap_size = overlap_size,
        nrep = nrep,
        losses = c("sse", "auc_as_loss"),
        model_candidates = c("SBM", "DCBM"),
        dcbm_est_options = list(estimate = list(method = "plugin")),
        ncores = ncores,
        laplacian = FALSE,
        regularize_tau = 0,
        retain_intermediates = "minimal"
      )
    })[["elapsed"]]
    best_sse <- fit$overall_best$best_model[fit$overall_best$loss == "sse"]
    best_auc <- fit$overall_best$best_model[fit$overall_best$loss == "auc_as_loss"]
    cat(sprintf("[%d/%d] NETCROP R=%d: SSE=%s, AUC=%s (%.2fs)\n",
                simulation, nsim, nrep, best_sse, best_auc, elapsed))
    print(fit$overall_best, row.names = FALSE)
    row_id <- row_id + 1L
    rows[[row_id]] <- data.frame(
      simulation = simulation, data = "DBLP-Conf", algorithm = "NETCROP",
      num_subnetworks = num_subnetworks, overlap_size = overlap_size,
      nrep = nrep, best_sse = best_sse, best_auc = best_auc,
      test_auc_at_best_sse = chosen_auc(fit, best_sse),
      test_auc_at_best_auc = chosen_auc(fit, best_auc),
      elapsed_seconds = elapsed
    )
    fits$NETCROP[[as.character(nrep)]] <- fit
  }

  for (algorithm in c("NCV", "ECV")) {
    for (nrep in c(1L, 20L)) {
      elapsed <- system.time({
        fit <- if (algorithm == "NCV") {
          ncv_stability_blockmodel(
            A = A, max_K = 10L, cv = 3L, nrep = nrep,
            dc_est = "plugin", tau = 0, use_laplacian = FALSE,
            losses = c("sse", "auc_as_loss"), ncores = ncores,
            retain_intermediates = "minimal"
          )
        } else {
          ecv_stability_blockmodel(
            A = A, max_K = 10L, train_proportion = 0.9,
            cv = 3L, nrep = nrep, tau = 0,
            losses = c("sse", "auc_as_loss"), ncores = ncores,
            retain_intermediates = "minimal"
          )
        }
      })[["elapsed"]]
      best_sse <- fit$overall_best$best_model[fit$overall_best$loss == "sse"]
      best_auc <- fit$overall_best$best_model[fit$overall_best$loss == "auc_as_loss"]
      cat(sprintf("[%d/%d] %s R=%d: SSE=%s, AUC=%s (%.2fs)\n",
                  simulation, nsim, algorithm, nrep, best_sse, best_auc, elapsed))
      print(fit$overall_best, row.names = FALSE)
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        simulation = simulation, data = "DBLP-Conf", algorithm = algorithm,
        num_subnetworks = NA_integer_, overlap_size = NA_integer_,
        nrep = nrep, best_sse = best_sse, best_auc = best_auc,
        test_auc_at_best_sse = chosen_auc(fit, best_sse),
        test_auc_at_best_auc = chosen_auc(fit, best_auc),
        elapsed_seconds = elapsed
      )
      fits[[algorithm]][[as.character(nrep)]] <- fit
    }
  }

  list(summary = do.call(rbind, rows), fits = fits)
}

results_file <- file.path(output_dir, "DBLP_results.rds")
records <- run_simulations(
  one_simulation = one_simulation,
  nsim = nsim,
  results_file = results_file,
  action = action,
  show_progress = TRUE
)

successful <- records[vapply(records, function(x) isTRUE(x$success), logical(1))]
summary_table <- do.call(rbind, lapply(successful, function(x) x$result$summary))
utils::write.csv(summary_table, file.path(output_dir, "DBLP_results.csv"),
                 row.names = FALSE)
cat("\nTable 4 DBLP summary:\n")
print(aggregate(
  cbind(accuracy = 100 * (summary_table$best_sse == "DCBM-4"),
        test_auc = summary_table$test_auc_at_best_sse,
        elapsed_seconds = summary_table$elapsed_seconds),
  summary_table[c("algorithm", "nrep")], mean
), row.names = FALSE)
