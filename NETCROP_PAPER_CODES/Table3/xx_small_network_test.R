setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table3")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))
source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "LSM_helpers.R"))

RUN_NAME <- "xx_small_network_test"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores)) 1L else min(2L, detected.cores)
nsim <- 2L
n <- 500L
d <- 3L
max.d <- 5L
alpha <- 0
loss.use <- "l2"

p.test <- 0.1
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o
R <- c(1L, 5L)

nc.file <- file.path(OUTPUT_DIR, "small_netcrop.csv")
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R", "s", "o")
)

for (sim in nc.simulations) {
  net <- LSM.gen(n = n, d = d, alpha = alpha, ncore = ncore,
                 seed = 100 + sim)
  for (R.use in R) {
    timing <- system.time({
      result <- netcrop_lsm(
        A = net$A, d.cand = max.d, s = s, o = o, R = R.use,
        loss = loss.use, ncore = ncore, seed = 500 + sim,
        step.size = 0.3, niter = 100, trace = 0
      )
    })
    netcrop_status(sim, nsim, "NETCROP", timing[3],
                   result$`d.hat.each.rep (l2)`, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, d = d, max_d = max.d, alpha = alpha,
      p_test = p.test, s = s, o = o, R = R.use,
      d_hat = result$l2.model, run_time = timing[3]
    ), nc.file, append = file.exists(nc.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_netcrop_sim", sim, "_R", R.use, ".rds")))
  }
}

message("Small Table 3 test completed. Results: ", normalizePath(OUTPUT_DIR))

################################################################################
nc.all <- readr::read_csv(nc.file, show_col_types = FALSE)


print(nc.all |> dplyr::group_by(s, o, R) |>
        dplyr::summarize(
          nsim = dplyr::n(),
          mean.time = mean(run_time),
          accu = 100 * mean(d_hat == d)
        ))
