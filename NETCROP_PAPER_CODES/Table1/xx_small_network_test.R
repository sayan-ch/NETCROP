setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table1")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))
source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

RUN_NAME <- "xx_small_network_test"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores)) 1L else min(2L, detected.cores) # it can be set to anything else
nsim <- 2L
n <- 500L
K <- 3L
max.K <- 5L
out.in.ratio <- 0.3
alpha <- 0.1
B <- alpha * (diag(1 - out.in.ratio, K) + out.in.ratio)
loss.use <- "l2"
model <- "SBM"

p.test <- 0.1
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

nc.file <- file.path(OUTPUT_DIR, "small_netcrop.csv")
ncv.file <- file.path(OUTPUT_DIR, "small_ncv.csv")
ecv.file <- file.path(OUTPUT_DIR, "small_ecv.csv")

R <- c(1L, 5L)
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R", "s", "o")
)
for (sim in nc.simulations) {
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)
  for (R.use in R) {
    timing <- system.time({
      result <- netcrop_blockmodel(
        A = net$A, K.CAND = seq_len(max.K), s = s, o = o, R = R.use,
        tau = 0, laplace = FALSE, dc.est = 2, loss = loss.use,
        mod.cand = c("SBM", "DCBM"), ncore = ncore, seed = 2 + sim * 100,
        rngR = TRUE # recommended for small networks for better clara agreements
      )
    })
    estimate <- result$`Mod.K.hat.each.rep (l2)`
    netcrop_status(sim, nsim, "NETCROP", timing[3], estimate, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, K = K, max_K = max.K, p_test = p.test,
      s = s, o = o, R = R.use, best_model = result$l2.model,
      run_time = timing[3]
    ), nc.file, append = file.exists(nc.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_netcrop_sim", sim, "_R", R.use, ".rds")))
  }
}

R <- c(1L, 20L)
ncv.simulations <- netcrop_resume_csv(
  ncv.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R")
)
for (sim in ncv.simulations) {
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)
  for (R.use in R) {
    timing <- system.time({
      result <- NCV.stability.BM(
        A = net$A, max.K = max.K, cv = 3, R = R.use, tau = 0,
        laplace = FALSE, dc.est = 2, loss = loss.use,
        ncore = ncore, seed = 2 + sim * 100
      )
    })
    netcrop_status(sim, nsim, "NCV", timing[3],
                   result$best.l2.each.rep, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, K = K, max_K = max.K, R = R.use,
      best_model = result$best.l2.stable, run_time = timing[3]
    ), ncv.file, append = file.exists(ncv.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_ncv_sim", sim, "_R", R.use, ".rds")))
  }
}

ecv.simulations <- netcrop_resume_csv(
  ecv.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R")
)
for (sim in ecv.simulations) {
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)
  for (R.use in R) {
    timing <- system.time({
      result <- ECV.stability.BM(
        A = net$A, max.K = max.K, train.p = 0.9, cv = 3, R = R.use,
        tau = 0, dc.est = 2, loss = loss.use,
        ncore = ncore, seed = 2 + sim * 100
      )
    })
    netcrop_status(sim, nsim, "ECV", timing[3],
                   result$best.l2.each.rep, R.use)
    readr::write_csv(tibble::tibble(
      nsim = sim, n = n, K = K, max_K = max.K, R = R.use,
      best_model = result$best.l2.stable, run_time = timing[3]
    ), ecv.file, append = file.exists(ecv.file))
    saveRDS(result, file.path(LOG_DIR,
      paste0("small_ecv_sim", sim, "_R", R.use, ".rds")))
  }
}

message("Small Table 1 test completed. Results: ", normalizePath(OUTPUT_DIR))

################################################################################
################################################################################
## Summarizing all the results
library(dplyr)

true_model <- paste0(model, "-", K)

nc.all <- readr::read_csv(nc.file, show_col_types = FALSE)

print(nc.all |> dplyr::group_by(s, o, R) |>
        dplyr::summarize(
          method = "NETCROP",
          nsim = dplyr::n(),
          mean.time = mean(run_time),
          accu = 100 * mean(best_model == true_model)
        ))

if (file.exists(ncv.file)) {
  ncv.all <- readr::read_csv(ncv.file, show_col_types = FALSE)
  print(ncv.all |> dplyr::group_by(R) |>
          dplyr::summarize(
            method = "NCV",
            nsim = dplyr::n(),
            mean.time = mean(run_time),
            accu = 100 * mean(best_model == true_model)
          ))
}


if (file.exists(ecv.file)) {
  ecv.all <- readr::read_csv(ecv.file, show_col_types = FALSE)
  print(ecv.all |> dplyr::group_by(R) |>
          dplyr::summarize(
            method = "ECV",
            nsim = dplyr::n(),
            mean.time = mean(run_time),
            accu = 100 * mean(best_model == true_model)
          ))
}
