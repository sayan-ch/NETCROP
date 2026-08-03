setwd(file.path(here::here("NETCROP_PAPER_CODES", "Figure2")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))
source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "PARTUNE_RSC_helpers.R"))

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
K <- 3L
K.max <- 5L
beta <- 1 / 3
rho <- 0.3

p.test <- 0.1
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

results.file <- file.path(OUTPUT_DIR, "small_figure2_results.rds")
resume.state <- netcrop_resume_rds(results.file, nsim, OUTPUT_ACTION)
results <- resume.state$results

for (sim in resume.state$simulations) {
  net <- DCBM.gen(n = n, K = K, beta = beta, rho = rho,
                  ncore = ncore, seed = 100 + sim)

  netcrop.time <- system.time({
    netcrop.result <- netcrop.tune.regsp(
      A = net$A, K = K, tau.cand = seq(0, 2, 0.1), DCBM = TRUE,
      s = s, o = o, R = 5, laplace = TRUE, dc.est = 2,
      loss = "l2", true.g = net$member,
      ncore = ncore, seed = 200 + 10 * sim
    )
  })[3]
  netcrop_status(sim, nsim, "NETCROP", netcrop.time,
                 netcrop.result$croissant.all.accu["l2"])

  dk.time <- system.time({
    dk.result <- DKest(
      A = net$A, K = K, true.g = net$member,
      tau.cand = seq(0, 0.1, 0.01), laplace = TRUE,
      DCBM = TRUE, DC.est = 2, ncore = ncore
    )
  })[3]
  best.dk <- dk.result[which.min(dk.result[, "DK.stat"]), "tau.cand"]
  netcrop_status(sim, nsim, "Davis-Kahan", dk.time, best.dk)

  results[[sim]] <- list(
    net = net, nc.out = netcrop.result, nc.time = netcrop.time,
    dk.out = dk.result, dk.time = dk.time,
    n = n, K = K, K.max = K.max, p.test = p.test
  )
  saveRDS(results[[sim]], file.path(LOG_DIR,
    paste0("small_figure2_sim", sim, ".rds")))
  saveRDS(results, results.file)
}

message("Small Figure 2 test completed. Results: ", normalizePath(OUTPUT_DIR))
