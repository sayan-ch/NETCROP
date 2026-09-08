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
one_simulation <- function(simulation) {
  model <- "SBM"
  cat(sprintf("\nTable 1 quick test | simulation %d/%d | %s\n", simulation, nsim, model))
  P_block <- matrix(0.03, K, K); diag(P_block) <- 0.1
  A <- generate_sbm(n = n, K = K, P_block = P_block, ncores = ncores)
  calls <- list()
  for (r in c(1L, 5L)) calls[[paste0("NETCROP_R", r)]] <- local({rr <- r; function() netcrop_blockmodel(A, K_candidates = seq_len(max_K),
      num_subnetworks = partition$num_subnetworks[1L], overlap_size = partition$overlap_size[1L],
      nrep = rr, losses = "sse", dcbm_est_options = list(estimate = list(method = "plugin")),
      ncores = ncores, verbose = FALSE, retain_intermediates = "minimal")})
  for (r in c(1L, 20L)) {
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
    best <- fit$overall_best$best_model[1L]; cat(best, "\n")
    pieces <- strsplit(method, "_R", fixed = TRUE)[[1L]]
    rows[[method]] <- data.frame(simulation, model, algorithm = pieces[1L], nrep = as.integer(pieces[2L]),
      best_model = best, run_time = elapsed)
    fits[[method]] <- fit
  }
  list(data = do.call(rbind, rows), fits = fits)
}
records <- run_simulations(one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "Table1_quick_test.rds"), action = "replace",
  show_progress = TRUE, continue_on_error = FALSE)
tab <- do.call(rbind, lapply(records, function(x) x$result$data))
utils::write.csv(tab, file.path(output_dir, "Table1_quick_test.csv"), row.names = FALSE)
print(aggregate(best_model ~ model + algorithm, tab, function(x) paste(names(sort(table(x), decreasing = TRUE))[1L], sprintf("(%d/%d)", max(table(x)), length(x)))))
