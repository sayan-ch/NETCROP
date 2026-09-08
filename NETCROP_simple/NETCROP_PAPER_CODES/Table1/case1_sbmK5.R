# Table 1, case 1: SBM, n = 10000, K = 5.
sourced_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (!is.null(sourced_file)) sourced_file else if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else NULL
script_dir <- if (!is.null(script_file)) dirname(normalizePath(script_file)) else getwd()
output_dir <- file.path(script_dir, "output", "case1_sbmK5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(netOP))
if (utils::packageVersion("netOP") != "0.1.1") stop("These replication scripts require netOP 0.1.1.", call. = FALSE)

ncores <- max(1L, floor(parallel::detectCores() / 2), na.rm = TRUE)
nsim <- 100L
n <- 10000L; K <- 5L; ratio <- 0.3; alpha <- 0.1; max_K <- 10L
P_block <- matrix(alpha * ratio, K, K); diag(P_block) <- alpha
partition <- netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)

phase <- "NETCROP"
one_simulation <- function(simulation) {
  cat(sprintf("\nTable 1 case 1 | simulation %d/%d | generating SBM\n", simulation, nsim))
  A <- generate_sbm(n = n, K = K, community_probabilities = rep(1 / K, K),
                    P_block = P_block, ncores = ncores)
  lambda <- mean(Matrix::rowSums(A)); rows <- list(); fits <- list(); id <- 0L
  run_method <- function(method, reps) {
    cat(sprintf("  %s R=%d ... ", method, reps)); started <- proc.time()[["elapsed"]]
    fit <- switch(method,
      NETCROP = netcrop_blockmodel(A, K_candidates = seq_len(max_K),
        num_subnetworks = partition$num_subnetworks[1L], overlap_size = partition$overlap_size[1L],
        nrep = reps, losses = "sse", model_candidates = c("SBM", "DCBM"),
        dcbm_est_options = list(estimate = list(method = "plugin")), ncores = ncores,
        retain_intermediates = "minimal", verbose = FALSE),
      NCV = ncv_stability_blockmodel(A, max_K = max_K, cv = 3L, nrep = reps,
        dc_est = "plugin", losses = "sse", ncores = ncores, retain_intermediates = "minimal", verbose = FALSE),
      ECV = ecv_stability_blockmodel(A, max_K = max_K, train_proportion = 0.9,
        cv = 3L, nrep = reps, losses = "sse", ncores = ncores,
        retain_intermediates = "minimal", verbose = FALSE))
    elapsed <- proc.time()[["elapsed"]] - started
    best <- fit$overall_best$best_model[fit$overall_best$loss == "sse"][1L]
    cat(sprintf("selected %s (%.3f s)\n", best, elapsed))
    print(fit$overall_best, row.names = FALSE)
    id <<- id + 1L
    rows[[id]] <<- data.frame(simulation, model = "SBM", n, K, out_in_ratio = ratio,
      alpha, average_degree = lambda, max_K, loss = "sse", algorithm = method,
      num_subnetworks = if (method == "NETCROP") partition$num_subnetworks[1L] else NA_integer_,
      overlap_size = if (method == "NETCROP") partition$overlap_size[1L] else NA_integer_,
      nrep = reps, best_model = best, run_time = elapsed)
    fits[[paste(method, reps, sep = "_R")]] <<- fit
  }
  if (phase == "NETCROP") {
    for (r in c(1L, 5L)) run_method("NETCROP", r)
  } else {
    for (r in c(1L, 20L)) run_method("NCV", r)
    for (r in c(1L, 20L)) run_method("ECV", r)
  }
  data <- do.call(rbind, rows)
  print(data, row.names = FALSE)
  list(data = data, fits = fits)
}

summarize_results <- function(results) {
  results$selected_K <- as.integer(sub(".*-", "", results$best_model))
  groups <- split(results, interaction(results$model, results$algorithm, results$nrep, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(x) {
    counts <- sort(table(x$best_model), decreasing = TRUE)
    data.frame(model = x$model[1L], algorithm = x$algorithm[1L], R = x$nrep[1L],
      best_model = names(counts)[1L], selected_count = as.integer(counts[1L]),
      selected_percent = 100 * as.integer(counts[1L]) / nrow(x),
      accuracy = 100 * mean(x$best_model == paste0(x$model, "-", x$K)),
      MAD = mean(abs(x$selected_K - x$K)))
  }))
  print(out, row.names = FALSE)
  invisible(out)
}

records_netcrop <- run_simulations(one_simulation, nsim = nsim,
  results_file = file.path(output_dir, "case1_sbmK5_netcrop.rds"), action = "resume",
  show_progress = TRUE, continue_on_error = TRUE)
tables <- lapply(records_netcrop, function(x) if (isTRUE(x$success)) x$result$data else NULL)
tables <- Filter(Negate(is.null), tables)
if (length(tables)) summarize_results(do.call(rbind, tables))

warning_text <- paste("NCV and ECV are very slow and memory consuming on this",
  "network and may crash R on a personal computer.")
message("\nWARNING: ", warning_text)
proceed <- TRUE
if (interactive()) proceed <- tolower(trimws(readline("Proceed with NCV and ECV? [y/N]: "))) %in% c("y", "yes")
if (proceed) {
  phase <- "comparison"
  records_comparison <- run_simulations(one_simulation, nsim = nsim,
    results_file = file.path(output_dir, "case1_sbmK5_comparison.rds"), action = "resume",
    show_progress = TRUE, continue_on_error = TRUE)
  comparison <- lapply(records_comparison, function(x) if (isTRUE(x$success)) x$result$data else NULL)
  tables <- c(tables, Filter(Negate(is.null), comparison))
} else message("NCV and ECV were skipped; NETCROP results remain saved.")
if (length(tables)) {
  results <- do.call(rbind, tables)
  utils::write.csv(results, file.path(output_dir, "case1_sbmK5.csv"), row.names = FALSE)
  cat("\nFinal Table 1 case 1 summary:\n")
  summarize_results(results)
}
