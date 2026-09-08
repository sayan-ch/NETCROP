# Quick Figure S3 test (n = 100, 200; 10 simulations): runtime, R-managed memory, and selection accuracy on small networks.
sourced_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (!is.null(sourced_file)) sourced_file else if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else NULL
script_dir <- if (!is.null(script_file)) dirname(normalizePath(script_file)) else getwd()
output_dir <- file.path(script_dir, "output", "xx_quicker_FigureS3_small_networks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(netOP))
if (utils::packageVersion("netOP") != "0.1.1") stop("These replication scripts require netOP 0.1.1.", call. = FALSE)

nsim <- 10L
n_values <- c(100L, 200L)
models <- c("SBM", "DCBM")
ncores_outer <- max(1L, floor(parallel::detectCores() / 2), na.rm = TRUE)

measure_call <- function(call) {
  if (requireNamespace("peakRAM", quietly = TRUE)) {
    return(netOP::measure_peak_ram(call))
  }
  invisible(gc(reset = TRUE))
  elapsed <- system.time(result <- call())[["elapsed"]]
  memory <- sum(gc()[, 7L])
  list(result = result, metrics = list(
    elapsed_seconds = elapsed,
    peak_ram_used_mib = memory,
    memory_method = "base R gc high-water mark"
  ))
}

run_configuration <- function(n, model, configuration) {
  K <- 3L; ratio <- 0.3; density <- 0.25
  partition <- netcrop_param_select(test_prop = 0.1, n = n, o_range = 0)
  one_simulation <- function(simulation) {
    cat(sprintf("\nFigure S3 config %d | n=%d %s | simulation %d/%d | generate\n",
                configuration, n, model, simulation, nsim))
    if (model == "SBM") {
      P_block <- matrix(ratio, K, K); diag(P_block) <- 1
      A <- generate_sbm(n = n, K = K, P_block = P_block,
        average_degree = n * density, average_degree_method = "naive", ncores = 1L)
    } else {
      psi_pool <- 1 / stats::rbeta(300L, 4, 1)
      A <- generate_dcbm(n = n, K = K, alpha = 1, beta = ratio,
        psi = sample(psi_pool, n, replace = TRUE), average_degree = n * density,
        average_degree_method = "naive", ncores = 1L)
    }
    calls <- list()
    for (r in c(1L, 3L, 5L)) calls[[paste0("NETCROP_R", r)]] <- local({
      rr <- r; function() netcrop_blockmodel(A, K_candidates = 1:5,
        num_subnetworks = partition$num_subnetworks[1L], overlap_size = partition$overlap_size[1L],
        nrep = rr, losses = "sse", model_candidates = c("SBM", "DCBM"),
        dcbm_est_options = list(estimate = list(method = "plugin")), ncores = 1L,
        verbose = FALSE, retain_intermediates = "minimal")
    })
    for (r in c(1L, 20L)) {
      calls[[paste0("NCV_R", r)]] <- local({rr <- r; function() ncv_stability_blockmodel(
        A, max_K = 5L, cv = 3L, nrep = rr, dc_est = "plugin", losses = "sse",
        ncores = 1L, verbose = FALSE, retain_intermediates = "minimal")})
      calls[[paste0("ECV_R", r)]] <- local({rr <- r; function() ecv_stability_blockmodel(
        A, max_K = 5L, train_proportion = 0.9, cv = 3L, nrep = rr, losses = "sse",
        ncores = 1L, verbose = FALSE, retain_intermediates = "minimal")})
    }
    rows <- list(); fits <- list()
    for (label in names(calls)) {
      pieces <- strsplit(label, "_R", fixed = TRUE)[[1L]]
      cat(sprintf("  simulation %d: %s R=%s ... ", simulation, pieces[1L], pieces[2L]))
      attempt <- tryCatch({
        measured <- measure_call(calls[[label]])
        list(ok = TRUE, fit = measured$result,
             elapsed = measured$metrics$elapsed_seconds,
             memory = measured$metrics$peak_ram_used_mib)
      }, error = function(e) list(ok = FALSE, error = conditionMessage(e), elapsed = -1, memory = -1))
      elapsed <- attempt$elapsed; memory <- attempt$memory
      best <- if (attempt$ok) attempt$fit$overall_best$best_model[1L] else NA_character_
      cat(sprintf("%s (%.3f s, %.1f MiB)\n", ifelse(is.na(best), "FAILED", best), elapsed, memory))
      rows[[label]] <- data.frame(looper = configuration, simulation, model, n, K,
        community_balance = "balanced", beta = ratio, density, test_prop = if (pieces[1L] == "NETCROP") 0.1 else NA_real_,
        o_range = if (pieces[1L] == "NETCROP") 0 else NA_real_, algorithm = pieces[1L],
        num_subnetworks = if (pieces[1L] == "NETCROP") partition$num_subnetworks[1L] else NA_integer_,
        overlap_size = if (pieces[1L] == "NETCROP") partition$overlap_size[1L] else NA_integer_,
        nrep = as.integer(pieces[2L]), best_model = best, run_time = if (attempt$ok) elapsed else -1,
        ram_MiB = if (attempt$ok) memory else -1, error = if (attempt$ok) NA_character_ else attempt$error)
      if (attempt$ok) fits[[label]] <- attempt$fit
    }
    data <- do.call(rbind, rows)
    print(data, row.names = FALSE)
    list(data = data, fits = fits)
  }
  file <- file.path(output_dir, sprintf("FigureS3_config%02d_%s_n%d.rds", configuration, model, n))
  run_simulations(one_simulation, nsim = nsim, use_parallel_simulations = TRUE,
    ncores_outer = ncores_outer, results_file = file, action = "replace",
    continue_on_error = TRUE, show_progress = TRUE)
}

configuration <- 0L; all_records <- list()
for (model in models) for (n in n_values) {
  configuration <- configuration + 1L
  cat(sprintf("\n===== Figure S3 configuration %d/%d: %s n=%d =====\n", configuration, length(models) * length(n_values), model, n))
  all_records[[configuration]] <- run_configuration(n, model, configuration)
}
rows <- unlist(lapply(all_records, function(records) lapply(records, function(x) if (isTRUE(x$success)) x$result$data else NULL)), recursive = FALSE)
rows <- Filter(Negate(is.null), rows); results <- do.call(rbind, rows)
utils::write.csv(results, file.path(output_dir, "FigureS3_small_networks.csv"), row.names = FALSE)

valid <- results[is.finite(results$run_time) & results$run_time >= 0, ]
valid$correct <- valid$best_model == paste0(valid$model, "-", valid$K)
summary_mean <- aggregate(cbind(correct, run_time, ram_MiB) ~ model + n + algorithm + nrep, valid, mean)
summary_n <- aggregate(simulation ~ model + n + algorithm + nrep, valid, length)
names(summary_n)[5L] <- "nsim"
summary_table <- merge(summary_mean, summary_n, by = c("model", "n", "algorithm", "nrep"))
summary_table$accuracy_percent <- 100 * summary_table$correct
utils::write.csv(summary_table, file.path(output_dir, "FigureS3_summary.csv"), row.names = FALSE)

valid$selected_K <- as.integer(sub(".*-", "", valid$best_model))
selection_groups <- split(valid, interaction(valid$model, valid$algorithm, valid$nrep, drop=TRUE))
selection_summary <- do.call(rbind, lapply(selection_groups, function(x) {
  counts <- sort(table(x$best_model), decreasing=TRUE)
  data.frame(model=x$model[1L], algorithm=x$algorithm[1L], R=x$nrep[1L],
    best_model=names(counts)[1L], selected_count=as.integer(counts[1L]),
    selected_percent=100*as.integer(counts[1L])/nrow(x),
    accuracy=100*mean(x$selected_K==x$K), MAD=mean(abs(x$selected_K-x$K)))
}))
cat("\nFinal Figure S3 selection summary:\n")
print(selection_summary, row.names=FALSE)

png(file.path(output_dir, "FigureS3_small_plot.png"), width = 2400, height = 1500, res = 300)
old <- par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
colors <- c(NCV = "#0072B2", ECV = "#D55E00", NETCROP = "#009E73")
for (model in models) for (measure in c("accuracy_percent", "run_time", "ram_MiB")) {
  block <- summary_table[summary_table$model == model, ]; y <- if (measure == "run_time") log(block[[measure]]) else block[[measure]]
  plot(range(n_values), range(y, finite = TRUE), type = "n", xlab = "Number of nodes (n)",
       ylab = switch(measure, accuracy_percent = "Accuracy (%)", run_time = "log(Mean runtime (sec.))", ram_MiB = "Mean R-managed RAM (MiB)"), main = paste0(model, "-3"))
  for (key in unique(paste(block$algorithm, block$nrep))) {
    take <- paste(block$algorithm, block$nrep) == key; ord <- order(block$n[take]); alg <- block$algorithm[take][1L]
    lines(block$n[take][ord], y[take][ord], type = "b", col = colors[alg], lty = 1L + (block$nrep[take][1L] %% 5L), pch = match(alg, names(colors)))
  }
}
dev.off(); par(old)
cat("Figure S3 outputs written to", normalizePath(output_dir), "\n")
