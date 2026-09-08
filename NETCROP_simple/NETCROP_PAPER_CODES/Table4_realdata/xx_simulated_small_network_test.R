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

phase <- "NETCROP"
one_simulation <- function(simulation) {
  psi_pool <- 1 / stats::rbeta(300L, 4, 1)
  psi <- sample(psi_pool, n, replace = TRUE)
  A <- generate_dcbm(
    n = n, K = K, alpha = 0.5, beta = 0.3, psi = psi,
    degree_scale = "none", average_degree = 30,
    representation = "sparse", ncores = ncores
  )

  fits <- list()
  if (phase == "NETCROP") for (nrep in c(1L, 5L)) {
    fits[[paste0("NETCROP_R", nrep)]] <- netcrop_blockmodel(
        A, K_candidates = seq_len(max_K),
        num_subnetworks = params$num_subnetworks[[1L]],
        overlap_size = params$overlap_size[[1L]], nrep = nrep,
        losses = c("sse", "auc_as_loss"), ncores = ncores,
        dcbm_est_options = list(estimate = list(method = "plugin")),
        retain_intermediates = "minimal"
      )
  }
  if (phase == "comparison") for (nrep in c(1L, 20L)) {
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

  rows <- lapply(names(fits), function(algorithm) {
    fit <- fits[[algorithm]]
    best <- fit$overall_best$best_model[fit$overall_best$loss == "sse"]
    cat(sprintf("[%d/%d] %s selected %s\n",
                simulation, nsim, algorithm, best))
    print(fit$overall_best, row.names = FALSE)
    name_parts <- strsplit(algorithm, "_R", fixed = TRUE)[[1L]]
    data.frame(simulation = simulation, model = "DCBM", K = K, algorithm = name_parts[[1L]],
               nrep = as.integer(name_parts[[2L]]), best_model = best,
               elapsed_seconds = fit$timing[["total"]])
  })
  data <- do.call(rbind, rows)
  print(data, row.names = FALSE)
  list(summary = data, fits = fits)
}

summarize_results<-function(results){results$selected_K<-as.integer(sub(".*-","",results$best_model));groups<-split(results,interaction(results$model,results$algorithm,results$nrep,drop=TRUE));out<-do.call(rbind,lapply(groups,function(x){counts<-sort(table(x$best_model),decreasing=TRUE);data.frame(model=x$model[1L],algorithm=x$algorithm[1L],R=x$nrep[1L],best_model=names(counts)[1L],selected_count=as.integer(counts[1L]),selected_percent=100*as.integer(counts[1L])/nrow(x),accuracy=100*mean(x$best_model==paste0(x$model,"-",x$K)),MAD=mean(abs(x$selected_K-x$K)))}));print(out,row.names=FALSE);invisible(out)}
records <- run_simulations(
  one_simulation = one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "small_table4_results_netcrop.rds"),
  action = "resume", show_progress = TRUE
)
stopifnot(all(vapply(records, function(x) isTRUE(x$success), logical(1))))
tables<-lapply(records,function(x)x$result$summary);cat("\nNETCROP quick-test summary:\n");summarize_results(do.call(rbind,tables))
message("\nWARNING: NCV and ECV are very slow and memory consuming and may crash R on a personal computer. The quick test will continue automatically.")
phase<-"comparison";comparison<-run_simulations(one_simulation=one_simulation,nsim=nsim,results_file=file.path(output_dir,"small_table4_results_comparison.rds"),action="resume",show_progress=TRUE)
stopifnot(all(vapply(comparison,function(x)isTRUE(x$success),logical(1))))
tables<-c(tables,lapply(comparison,function(x)x$result$summary));summary_table<-do.call(rbind,tables)
csv_file <- file.path(output_dir, "small_table4_results.csv")
if (file.exists(csv_file)) {
  archived_csv <- file.path(
    output_dir,
    paste0("small_table4_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  )
  file.rename(csv_file, archived_csv)
}
utils::write.csv(summary_table, csv_file, row.names = FALSE)
cat("\nFinal Table 4 quick-test summary:\n");summarize_results(summary_table)
cat("Small Table 4 test completed:", normalizePath(output_dir), "\n")
