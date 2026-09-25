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

compiler_help <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
  paste0(
    "\n\nOn macOS, if this error involves compilation, gfortran, an SDK, ",
    "or linker failures, see the \"Compilation fails on macOS\" section ",
    "of README.md. Apple Silicon systems may require Homebrew GCC and ",
    "~/.R/Makevars configuration."
  )
} else {
  "\nCheck the platform compiler instructions in README.md."
}

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

minimum_r_version <- numeric_version("4.4.0")

if (getRversion() < minimum_r_version) {
  stop(
    "NETCROP requires R 4.4.0 or newer. ",
    "You are using R ", getRversion(), ". ",
    "Please update R before running this setup script.",
    call. = FALSE
  )
}

recommended_packages <- c("Matrix", "cluster", "lattice")
missing_recommended <- recommended_packages[
  !vapply(recommended_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_recommended) > 0L) {
  stop(
    "Your R installation is missing recommended packages: ",
    paste(missing_recommended, collapse = ", "),
    "\nInstall the complete R distribution, then rerun this script.",
    call. = FALSE
  )
}

message("Restoring the project library from renv.lock...")
tryCatch(
  renv::restore(prompt = interactive()),
  error = function(error) {
    stop(
      "renv::restore() failed: ", conditionMessage(error),
      "\nCheck your internet connection, resolve any compiler problem, then run ",
      "source(\"first_time_renv_setup.R\") again.",
      compiler_help,
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
      compiler_help,
      call. = FALSE
    )
  }
)

message("NETCROP environment restored and helper smoke test passed.")
message("Next: source a folder's xx_small_network_test.R script.")
