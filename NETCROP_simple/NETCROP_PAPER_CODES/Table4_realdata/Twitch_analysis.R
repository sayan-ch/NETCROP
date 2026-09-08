#!/usr/bin/env Rscript

# Table 4: NETCROP block-model selection on the Twitch network.
# Requires netOP 0.1.1. Random seeds are intentionally not supplied.

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
output_dir <- file.path(root, "output", "Twitch_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
detected <- parallel::detectCores()
ncores <- if (is.na(detected)) 1L else max(1L, floor(detected / 2))
nsim <- 100L
action <- Sys.getenv("NETCROP_ACTION", "resume")
if (!action %in% c("replace", "resume", "archive")) {
  stop("NETCROP_ACTION must be replace, resume, or archive.")
}

edge_list <- utils::read.csv(file.path(root, "Twitch", "Twitch_edge_list.csv"))
n <- max(edge_list$i, edge_list$j)
A <- Matrix::sparseMatrix(
  i = edge_list$i, j = edge_list$j, x = 1,
  dims = c(n, n)
)
params <- netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)

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
  elapsed <- system.time({
    fit <- netcrop_blockmodel(
      A = A, K_candidates = 1:30,
      num_subnetworks = params$num_subnetworks[[1L]],
      overlap_size = params$overlap_size[[1L]], nrep = 1L,
      losses = c("sse", "auc_as_loss"),
      model_candidates = c("SBM", "DCBM"), ncores = ncores,
      laplacian = TRUE, regularize_tau = 0,
      sbm_est_options = list(spectral_cluster = list(
        cluster_engine = "clara",
        cluster_options = list(samples = 50L, sampsize = 200L)
      )),
      dcbm_est_options = list(
        spectral_cluster = list(
          cluster_engine = "clara",
          cluster_options = list(samples = 50L, sampsize = 200L)
        ),
        estimate = list(method = "plugin")
      ),
      retain_intermediates = "minimal"
    )
  })[["elapsed"]]
  best_model <- fit$overall_best$best_model[fit$overall_best$loss == "sse"]
  cat(sprintf("[%d/%d] NETCROP: SSE-selected model=%s (%.2fs)\n",
              simulation, nsim, best_model, elapsed))
  print(fit$overall_best[fit$overall_best$loss == "sse", ], row.names = FALSE)
  rows <- list(data.frame(
      simulation = simulation, data = "Twitch", algorithm = "NETCROP",
      num_subnetworks = params$num_subnetworks[[1L]],
      overlap_size = params$overlap_size[[1L]], nrep = 1L,
      best_model = best_model,
      test_auc = chosen_auc(fit, best_model),
      elapsed_seconds = elapsed
    ))
  data <- do.call(rbind, rows)
  print(data, row.names = FALSE)
  list(summary = data, fit = fit)
}

records <- run_simulations(
  one_simulation = one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "Twitch_results.rds"),
  action = action, show_progress = TRUE
)
successful <- records[vapply(records, function(x) isTRUE(x$success), logical(1))]
summary_table <- do.call(rbind, lapply(successful, function(x) x$result$summary))
utils::write.csv(summary_table, file.path(output_dir, "Twitch_results.csv"),
                 row.names = FALSE)
cat("\nTable 4 Twitch summary:\n")
groups<-split(summary_table,summary_table$best_model)
out<-do.call(rbind,lapply(groups,function(x)data.frame(algorithm=x$algorithm[1L],R=x$nrep[1L],best_model=x$best_model[1L],selected_count=nrow(x),selected_percent=100*nrow(x)/nrow(summary_table),test_auc=mean(x$test_auc,na.rm=TRUE))))
print(out[order(-out$selected_percent),],row.names=FALSE)
