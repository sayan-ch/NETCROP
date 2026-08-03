library(Matrix)
library(parallel)

################################################################################
.netcrop_cpp_file <- here::here(
  "NETCROP_PAPER_CODES",
  "helpers",
  "General_helpers.cpp"
)

.netcrop_native_functions <- c(
  "netcrop_auc_cpp",
  "netcrop_outer_add_cpp",
  "netcrop_procrustes_translated_cpp"
)

if (!all(vapply(
  .netcrop_native_functions,
  exists,
  logical(1),
  mode = "function",
  inherits = TRUE
))) {
  Rcpp::sourceCpp(.netcrop_cpp_file)
}

rm(.netcrop_cpp_file, .netcrop_native_functions)

netcrop_nmi <- function(g1, g2) {
  cross.table <- table(g1, g2)
  entropy::mi.empirical(as.matrix(cross.table)) /
    entropy::entropy.empirical(as.vector(cross.table))
}

# Universal singular value thresholding with irlba explicitly namespaced.
netcrop_usvt <- function(A) {
  n <- nrow(A)
  K <- ceiling(n^(1 / 3))
  decomposition <- irlba::irlba(A, nv = K, nu = K)
  estimate <- decomposition$u %*% (t(decomposition$v) * decomposition$d)
  estimate <- pmin(estimate, 1)
  estimate <- pmax(estimate, 0)
  estimate
}

# Rank-based binary AUC with average-rank handling for tied predictions.
netcrop_auc <- function(group, predictions) {
  netcrop_auc_cpp(group, predictions)
}

# The original outer-add routine returns length(y) rows by length(x) columns.
# Only addition is used by NETCROP, so unsupported operators fail explicitly.
netcrop_outer <- function(x, y, operator = "+") {
  if (!identical(operator, "+")) {
    stop("netcrop_outer() currently supports only operator = \"+\".", call. = FALSE)
  }
  netcrop_outer_add_cpp(x, y)
}

# Hybrid Procrustes implementation. Rotation-only alignment uses BLAS-backed
# base matrix operations, which benchmarked slightly faster than Eigen. The
# translated path uses the native implementation, which benchmarked about four
# times faster at 100,000-by-20 and avoids an N-by-N centering matrix.
netcrop_procrustes <- function(X, Xstar, translate = FALSE, dilate = FALSE,
                               sumsq = FALSE) {
  if (!is.logical(translate) || length(translate) != 1L || is.na(translate)) {
    stop("'translate' must be a single logical indicator", call. = FALSE)
  }
  if (!is.logical(dilate) || length(dilate) != 1L || is.na(dilate)) {
    stop("'dilate' must be a single logical indicator", call. = FALSE)
  }
  if (!is.logical(sumsq) || length(sumsq) != 1L || is.na(sumsq)) {
    stop("'sumsq' must be a single logical indicator", call. = FALSE)
  }
  if ((N <- nrow(X)) != nrow(Xstar)) {
    stop("X and Xstar do not have the same number of rows", call. = FALSE)
  }
  P <- ncol(X)
  P2 <- ncol(Xstar)
  if (P != P2) {
    if (P < P2) {
      warning(
        "X padded out to same number of columns as Xstar\n",
        call. = FALSE,
        immediate. = TRUE
      )
      X <- cbind(X, matrix(0, nrow = N, ncol = P2 - P))
      P <- P2
    } else {
      stop("X cannot have more columns than Xstar", call. = FALSE)
    }
  }
  if (P2 == 0L) {
    stop("Xstar must contain at least one column", call. = FALSE)
  }
  if (anyNA(Xstar) || anyNA(X)) {
    stop("X and Xstar are not allowed to contain missing values", call. = FALSE)
  }

  if (translate && !dilate) {
    native <- netcrop_procrustes_translated_cpp(X, Xstar)
    X.new <- unname(native$X.new)
    R <- if (P2 == 1L) 1L else native$R
    translation <- as.numeric(native$t)
    scale <- 1
  } else {
    if (translate) {
      mean_x <- colMeans(X)
      mean_xstar <- colMeans(Xstar)
      cross <- crossprod(Xstar, X) - N * tcrossprod(mean_xstar, mean_x)
    } else {
      cross <- crossprod(Xstar, X)
    }

    if (P2 == 1L) {
      R <- 1
    } else {
      decomposition <- svd(cross)
      R <- tcrossprod(decomposition$v, decomposition$u)
    }

    if (dilate) {
      denominator <- if (translate) {
        sum(X^2) - N * sum(colMeans(X)^2)
      } else {
        sum(X^2)
      }
      scale <- sum(cross * R) / denominator
    } else {
      scale <- 1
    }

    scaled_rotated <- scale * X %*% R
    translation <- if (translate) colMeans(Xstar - scaled_rotated) else 0
    X.new <- unname(if (translate) {
      sweep(scaled_rotated, 2L, translation, `+`)
    } else {
      scaled_rotated
    })
  }

  result <- list(X.new = X.new, R = R)
  if (translate) result$t <- translation
  if (dilate) result$d <- scale
  if (sumsq) result$ss <- sum((X[, seq_len(P2), drop = FALSE] - X.new)^2)
  result
}

################################################################################
# Simulation workflow helpers

netcrop_output_action <- function(output_dir, log_dir = NULL) {
  output_exists <- dir.exists(output_dir) &&
    length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L
  log_exists <- !is.null(log_dir) && dir.exists(log_dir) &&
    length(list.files(log_dir, all.files = TRUE, no.. = TRUE)) > 0L

  requested <- tolower(trimws(Sys.getenv("NETCROP_OUTPUT_ACTION", unset = "")))
  valid <- c("replace", "resume", "archive")
  if (nzchar(requested) && !requested %in% valid) {
    stop(
      "NETCROP_OUTPUT_ACTION must be replace, resume, or archive.",
      call. = FALSE
    )
  }

  if (!output_exists && !log_exists) {
    action <- if (nzchar(requested)) requested else "replace"
  } else if (nzchar(requested)) {
    action <- requested
  } else if (interactive()) {
    choice <- utils::menu(
      c(
        "Replace existing results and restart",
        "Resume from the last complete simulation",
        "Archive existing results and restart"
      ),
      title = paste0(
        "Existing results were found in ", normalizePath(output_dir, mustWork = FALSE),
        ". Choose how to continue:"
      )
    )
    if (choice == 0L) {
      stop("No output action selected; simulation cancelled.", call. = FALSE)
    }
    action <- valid[choice]
  } else {
    action <- "replace"
    message(
      "NETCROP_OUTPUT_ACTION was not set; batch mode will replace existing ",
      "results and restart."
    )
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (action == "replace") {
    if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)
    if (!is.null(log_dir) && dir.exists(log_dir)) {
      unlink(log_dir, recursive = TRUE, force = TRUE)
    }
  } else if (action == "archive") {
    archive_one <- function(path) {
      if (!dir.exists(path)) return(invisible(NULL))
      archived <- paste0(path, "_", timestamp)
      suffix <- 1L
      while (file.exists(archived)) {
        archived <- paste0(path, "_", timestamp, "_", suffix)
        suffix <- suffix + 1L
      }
      if (!file.rename(path, archived)) {
        stop("Could not archive directory: ", path, call. = FALSE)
      }
      message("Archived ", path, " to ", archived)
      invisible(archived)
    }
    archive_one(output_dir)
    if (!is.null(log_dir)) archive_one(log_dir)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  list(action = action, output_dir = output_dir, log_dir = log_dir)
}

netcrop_resume_csv <- function(file, nsim, action, expected_rows = 1L,
                               key_columns = "nsim",
                               simulation_column = "nsim") {
  all_simulations <- seq_len(nsim)
  if (!identical(action, "resume") || !file.exists(file)) {
    return(all_simulations)
  }

  existing <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
  if (!simulation_column %in% names(existing)) {
    stop("Cannot resume because ", file, " has no ", simulation_column,
         " column.", call. = FALSE)
  }

  available_keys <- intersect(key_columns, names(existing))
  if (length(available_keys) > 0L) {
    existing <- existing[!duplicated(existing[available_keys]), , drop = FALSE]
  }

  simulation <- suppressWarnings(as.integer(existing[[simulation_column]]))
  if (anyNA(simulation)) {
    stop("Cannot resume because the simulation column in ", file,
         " contains missing or non-integer values.", call. = FALSE)
  }

  counts <- table(factor(simulation, levels = all_simulations))
  complete <- as.integer(counts) >= as.integer(expected_rows)
  first_incomplete <- which(!complete)[1L]
  start <- if (is.na(first_incomplete)) nsim + 1L else first_incomplete

  keep <- simulation < start
  normalized <- existing[keep, , drop = FALSE]
  if (nrow(normalized) != nrow(existing) || file.exists(file)) {
    readr::write_csv(normalized, file)
  }

  if (start > nsim) {
    message("All ", nsim, " simulations are already complete in ", file, ".")
    return(integer())
  }

  message("Resuming ", file, " at simulation ", start, " of ", nsim, ".")
  seq.int(start, nsim)
}

netcrop_resume_rds <- function(file, nsim, action) {
  if (!identical(action, "resume") || !file.exists(file)) {
    return(list(results = vector("list", nsim), simulations = seq_len(nsim)))
  }

  results <- readRDS(file)
  if (!is.list(results)) {
    stop("Cannot resume because ", file, " does not contain a list.",
         call. = FALSE)
  }
  length(results) <- nsim
  complete <- !vapply(results, is.null, logical(1))
  first_incomplete <- which(!complete)[1L]
  start <- if (is.na(first_incomplete)) nsim + 1L else first_incomplete
  if (start <= nsim) results[start:nsim] <- vector("list", nsim - start + 1L)

  simulations <- if (start > nsim) integer() else seq.int(start, nsim)
  if (length(simulations) == 0L) {
    message("All ", nsim, " simulations are already complete in ", file, ".")
  } else {
    message("Resuming ", file, " at simulation ", start, " of ", nsim, ".")
  }
  list(results = results, simulations = simulations)
}

netcrop_status <- function(simulation, nsim, method, elapsed,
                           estimate = NULL, repetition = NULL) {
  fields <- c(
    sprintf("Simulation %d/%d", simulation, nsim),
    paste0("method=", method)
  )
  if (!is.null(repetition)) fields <- c(fields, paste0("R=", repetition))
  fields <- c(fields, sprintf("elapsed=%.3f s", as.numeric(elapsed)))
  if (!is.null(estimate)) {
    fields <- c(fields, paste0("estimate=", paste(estimate, collapse = ", ")))
  }
  cat(paste(fields, collapse = " | "), "\n")
}

netcrop_confirm_large_cv <- function(n, method, small_script) {
  if (n <= 1000) return(TRUE)

  command <- sprintf("source(%s)", encodeString(normalizePath(
    small_script, mustWork = FALSE
  ), quote = '"'))
  warning(
    method, " with n > 1000 is very time- and memory-intensive and may ",
    "cause memory overload.",
    call. = FALSE,
    immediate. = TRUE
  )
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_text("For a small-network test, run {.run {command}}")
  } else {
    message("For a small-network test, run: ", command)
  }

  requested <- tolower(trimws(Sys.getenv("NETCROP_RUN_LARGE_CV", unset = "")))
  yes <- c("1", "true", "yes", "y")
  no <- c("0", "false", "no", "n")
  if (nzchar(requested) && !requested %in% c(yes, no)) {
    stop(
      "NETCROP_RUN_LARGE_CV must be yes/true/1 or no/false/0.",
      call. = FALSE
    )
  }
  if (requested %in% yes) return(TRUE)
  if (requested %in% no) return(FALSE)

  if (!interactive()) {
    message(
      "NETCROP_RUN_LARGE_CV was not set; batch mode will skip ", method, "."
    )
    return(FALSE)
  }

  answer <- readline(sprintf("Proceed with %s? [y/N]: ", method))
  tolower(trimws(answer)) %in% yes
}
################################################################################
## General helpers: usually used in all functions=
modal <- function(x) { # computes modal value in a vector
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

l2 <- function(x, y) { #squared l2 distance
  sum((x - y)^2)
}

bin.dev <- function(x, y) { # binomial deviance loss: y: prob, x: binary response
  y <- pmax(y, 1e-5)
  y <- pmin(y, 1 - 1e-5)

  tmp <- -x * log(y) - (1 - x) * log(1 - y)

  return(sum(tmp, rm.na = T))
}

AUC <- function(AA, PP) {
  - netcrop_auc(group = as(AA, "numeric"), predictions = as(PP, "numeric"))
}

softplus <- function(x) pmax(x,0) + log1p(exp(-abs(x)))

sigmoid  <- function(x) 1/(1+exp(-x))
################################################################################
# NETCROP parameter selection
netcrop_param <- function(p.test = 0.02, n = NULL,
                          o.range = c(0, 0.8)){
  over.upper <- (1 - sqrt(p.test))
  over.lower <- 1 - sqrt(2*p.test)

  o.range[o.range == 1] <- 0.8

  over.prop <- over.lower + (over.upper - over.lower) * o.range

  s.range <- ceiling((1 - over.prop)^2 / ((1 - over.prop)^2 - p.test))

  if(is.null(n))
    return(data.frame(p.test = p.test, p_o = over.prop, s = s.range))

  over.range <- ceiling(n * over.prop)
  m.range <- floor((n - over.range) / s.range)

  over.extra <- n - over.range - s.range * m.range
  over.range <- over.range + over.extra

  return(data.frame(p.test = p.test, p_o = over.prop, n = n,
                    s = s.range, o = over.range))
}
