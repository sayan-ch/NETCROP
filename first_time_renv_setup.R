# Run this script once from the NETCROP repository root.

required_files <- c("NETCROP.Rproj", "renv.lock", "renv/activate.R")
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Run this script from the NETCROP repository root. Missing: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv from CRAN...")
  install.packages("renv")
}

lockfile <- renv::lockfile_read("renv.lock")
locked_r <- lockfile$R$Version
running_r <- paste(R.version$major, R.version$minor, sep = ".")
if (!is.null(locked_r) && !identical(locked_r, running_r)) {
  warning(
    "renv.lock records R ", locked_r, " but this session uses R ", running_r,
    ". Package compilation may differ across R versions.",
    call. = FALSE,
    immediate. = TRUE
  )
}

message("Restoring the project library from renv.lock...")
tryCatch(
  renv::restore(prompt = interactive()),
  error = function(error) {
    stop(
      "renv::restore() failed: ", conditionMessage(error),
      "\nCheck your internet connection and compiler installation, then run ",
      "source(\"first_time_renv_setup.R\") again.",
      call. = FALSE
    )
  }
)

message("Checking the restored environment...")
renv::status()

tryCatch(
  {
    source(file.path("NETCROP_PAPER_CODES", "helpers", "General_helpers.R"))
    stopifnot(
      exists("netcrop_auc", mode = "function"),
      exists("netcrop_outer", mode = "function"),
      exists("netcrop_procrustes", mode = "function")
    )
  },
  error = function(error) {
    stop(
      "Packages restored, but the compiled helper smoke test failed: ",
      conditionMessage(error),
      "\nInstall the platform compiler described in README.md and rerun this script.",
      call. = FALSE
    )
  }
)

message("NETCROP environment restored and helper smoke test passed.")
message("Next: source a folder's xx_small_network_test.R script.")
