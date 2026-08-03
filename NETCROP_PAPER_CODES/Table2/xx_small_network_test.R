setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table2")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))
source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "RDPG_helpers.R"))

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
xi <- 0.75
loss.use <- "l2"

p.test <- 0.1
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

nc.file <- file.path(OUTPUT_DIR, "small_netcrop.csv")
ecv.file <- file.path(OUTPUT_DIR, "small_ecv.csv")

R <- c(1L, 5L)
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R", "s", "o")
)
for (sim in nc.simulations) {
  net <- RDPG.gen(n = n, d = d, rho = xi, ncore = ncore,
                  seed = 200 + sim)
  for (R.use in R) {
    timing <- system.time({
      result <- netcrop_rdpg(
        A = net$A, d.cand = max.d, s = s, o = o, R = R.use,
        loss = loss.use, ncore = ncore, seed = 500 + sim
      )
    })
    netcrop_status(sim, nsim, "NETCROP", timing[3],
                   result$`d.hat.each.rep (l2)`, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, d = d, max_d = max.d, p_test = p.test,
      s = s, o = o, R = R.use, d_hat = result$l2.model,
      run_time = timing[3]
    ), nc.file, append = file.exists(nc.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_netcrop_sim", sim, "_R", R.use, ".rds")))
  }
}

R <- c(1L, 20L)
ecv.simulations <- netcrop_resume_csv(
  ecv.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R")
)
for (sim in ecv.simulations) {
  net <- RDPG.gen(n = n, d = d, rho = xi, ncore = ncore,
                  seed = 200 + sim)
  for (R.use in R) {
    timing <- system.time({
      result <- ECV.stability.RDPG(
        A = net$A, max.K = max.d, cv = 3, R = R.use, train.p = 0.9,
        loss = loss.use, ncore = ncore, seed = 2 + sim * 100
      )
    })
    netcrop_status(sim, nsim, "ECV", timing[3],
                   result$best.l2.each.rep, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, d = d, max_d = max.d, R = R.use,
      d_hat = result$best.l2.stable, run_time = timing[3]
    ), ecv.file, append = file.exists(ecv.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_ecv_sim", sim, "_R", R.use, ".rds")))
  }
}

message("Small Table 2 test completed. Results: ", normalizePath(OUTPUT_DIR))
