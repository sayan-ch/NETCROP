# Fast integration test for the Table 1 workflow (10 unseeded simulations).
sourced_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (!is.null(sourced_file)) sourced_file else if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else NULL
script_dir <- if (!is.null(script_file)) dirname(normalizePath(script_file)) else getwd()
output_dir <- file.path(script_dir, "output", "xx_small_network_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(netOP))
if (utils::packageVersion("netOP") != "0.1.1") stop("These replication scripts require netOP 0.1.1.", call. = FALSE)

nsim <- 10L; n <- 500L; K <- 3L; max_K <- 5L; ncores <- 1L
partition <- netcrop_param_select(test_prop = 0.1, n = n, o_range = 0)
phase <- "NETCROP"
one_simulation <- function(simulation) {
  model <- "SBM"
  cat(sprintf("\nTable 1 quick test | simulation %d/%d | %s\n", simulation, nsim, model))
  P_block <- matrix(0.03, K, K); diag(P_block) <- 0.1
  A <- generate_sbm(n = n, K = K, P_block = P_block, ncores = ncores)
  calls <- list()
  if (phase == "NETCROP") for (r in c(1L, 5L)) calls[[paste0("NETCROP_R", r)]] <- local({rr <- r; function() netcrop_blockmodel(A, K_candidates = seq_len(max_K),
      num_subnetworks = partition$num_subnetworks[1L], overlap_size = partition$overlap_size[1L],
      nrep = rr, losses = "sse", dcbm_est_options = list(estimate = list(method = "plugin")),
      ncores = ncores, verbose = FALSE, retain_intermediates = "minimal")})
  if (phase == "comparison") for (r in c(1L, 20L)) {
    calls[[paste0("NCV_R", r)]] <- local({rr <- r; function() ncv_stability_blockmodel(A, max_K = max_K, cv = 3L, nrep = rr,
      dc_est = "plugin", losses = "sse", ncores = ncores, verbose = FALSE, retain_intermediates = "minimal")})
    calls[[paste0("ECV_R", r)]] <- local({rr <- r; function() ecv_stability_blockmodel(A, max_K = max_K,
      train_proportion = 0.9, cv = 3L, nrep = rr, losses = "sse", ncores = ncores,
      verbose = FALSE, retain_intermediates = "minimal")})
  }
  rows <- list(); fits <- list()
  for (method in names(calls)) {
    cat("  ", method, " ... ", sep = ""); started <- proc.time()[["elapsed"]]
    fit <- calls[[method]](); elapsed <- proc.time()[["elapsed"]] - started
    best <- fit$overall_best$best_model[fit$overall_best$loss == "sse"][1L]; cat(best, "\n")
    print(fit$overall_best, row.names = FALSE)
    pieces <- strsplit(method, "_R", fixed = TRUE)[[1L]]
    rows[[method]] <- data.frame(simulation, model, algorithm = pieces[1L], nrep = as.integer(pieces[2L]),
      K = K, best_model = best, run_time = elapsed)
    fits[[method]] <- fit
  }
  data <- do.call(rbind, rows)
  print(data, row.names = FALSE)
  list(data = data, fits = fits)
}
summarize_results <- function(results) {
  results$selected_K <- as.integer(sub(".*-", "", results$best_model))
  groups <- split(results, interaction(results$model, results$algorithm, results$nrep, drop=TRUE))
  out <- do.call(rbind, lapply(groups, function(x) { counts <- sort(table(x$best_model), decreasing=TRUE); data.frame(
    model=x$model[1L], algorithm=x$algorithm[1L], R=x$nrep[1L], best_model=names(counts)[1L],
    selected_count=as.integer(counts[1L]), selected_percent=100*as.integer(counts[1L])/nrow(x),
    accuracy=100*mean(x$best_model==paste0(x$model,"-",x$K)), MAD=mean(abs(x$selected_K-x$K))) }))
  print(out, row.names=FALSE); invisible(out)
}
records_netcrop <- run_simulations(one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "Table1_quick_test_netcrop.rds"), action = "resume",
  show_progress = TRUE, continue_on_error = FALSE)
tab_netcrop <- do.call(rbind, lapply(records_netcrop, function(x) x$result$data))
cat("\nNETCROP quick-test summary:\n"); summarize_results(tab_netcrop)
message("\nWARNING: NCV and ECV are very slow and memory consuming and may crash R on a personal computer. The quick test will continue automatically.")
phase <- "comparison"
records_comparison <- run_simulations(one_simulation, nsim=nsim,
  results_file=file.path(output_dir,"Table1_quick_test_comparison.rds"), action="resume",
  show_progress=TRUE, continue_on_error=FALSE)
tab <- rbind(tab_netcrop, do.call(rbind, lapply(records_comparison, function(x) x$result$data)))
utils::write.csv(tab, file.path(output_dir, "Table1_quick_test.csv"), row.names = FALSE)
cat("\nFinal Table 1 quick-test summary:\n"); summarize_results(tab)
