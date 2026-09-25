setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table1")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

RUN_NAME <- "case2_sbmK20"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################

version <- 1
detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores) | .Platform$OS.type == "windows" ) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 100L

# Case 2 of Table 1:
# SBM with n = 10000 nodes and K = 20 communities
# out.in.ratio = 1/3 and alpha = 0.3
n <- 10^4
K <- 20
out.in.ratio <- 1/3
alpha <- 0.3
B <- alpha * ( diag(1-out.in.ratio, K) + out.in.ratio )
model <- "SBM"

PI <- rep(1/K, K)

################################################################################
################################################################################
################################################################################
################################################################################

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing lowest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.K <- 30
loss.use <- c("l2")

nc.file <- file.path(OUTPUT_DIR, paste0("case2_netcrop_v", version, ".csv"))
ncv.file <- file.path(OUTPUT_DIR, paste0("case2_ncv_v", version, ".csv"))
ecv.file <- file.path(OUTPUT_DIR, paste0("case2_ecv_v", version, ".csv"))
small.script <- file.path(here::here("NETCROP_PAPER_CODES", "Table1"),
                          "xx_small_network_test.R")

## Run NETCROP
all.nc <- list()
count <- 1
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION,
  expected_rows = length(R) * length(s) * length(o),
  key_columns = c("nsim", "R", "s", "o", "loss_function")
)
for(sim in nc.simulations){
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    for(ss in seq_along(s)){
      for(oo in seq_along(o)){
        s.use <- s[ss]
        o.use <- o[oo]
        R.use <- R[RR]

        time.nc <- system.time({
          out.nc <- netcrop_blockmodel(A = net$A, K.CAND = 1:max.K,
                                       s = s.use, o = o.use, R = R.use,
                                       tau = 0, laplace = F, dc.est = 2,
                                       loss = loss.use, mod.cand = c('SBM', 'DCBM'),
                                       ncore = ncore, seed = 2 + sim*100)
        })

        gc()

        netcrop_status(sim, nsim, "NETCROP", time.nc[3],
                       out.nc$`Mod.K.hat.each.rep (l2)`, R.use)

        nc.tab <- tibble::tibble(
          nsim = sim,
          model = model,
          n = n, K = K,
          out_in_ratio = out.in.ratio, alpha = alpha, lambda = lambda,
          max_K = max.K, loss_function = loss.use,
          s = s.use, o = o.use, R = R.use,
          best_model = do.call('c', out.nc[paste0(loss.use, ".model")]),
          run_time = time.nc[3],
          user_time = time.nc[1],
          system_time = time.nc[2]
        ) |>
          dplyr::mutate(
            model_hat = stringr::str_extract(best_model, "^[^-]+"),
            K_hat = as.integer(stringr::str_extract(best_model, "(?<=-)\\d+")),
            .before = run_time
          )

        readr::write_csv(
          nc.tab, file = nc.file, append = file.exists(nc.file)
        )

        all.nc[[count]] <- list(time = time.nc, out = out.nc)

        saveRDS(
          all.nc[[count]], file = file.path(
            LOG_DIR,
            paste0("case2_netcrop_v", version, "_sim", sim, "_R", R.use, ".rds")
          )
        )

        count <- count + 1
      }
    }
  }
}

################################################################################
################################################################################
################################################################################
################################################################################

## Run NCV
all.ncv <- list()
count <- 1
R <- c(1, 20)
run.ncv <- netcrop_confirm_large_cv(n, "NCV", small.script)
ncv.simulations <- netcrop_resume_csv(
  ncv.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R", "loss_function")
)
if (run.ncv) for(sim in ncv.simulations){
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    R.use <- R[RR]

    time.ncv <- system.time({
      out.ncv <- NCV.stability.BM(A = net$A, max.K = max.K, cv = 3, R = R.use,
                                  tau = 0, laplace = F, dc.est = 2,
                                  loss = loss.use,
                                  ncore = ncore, seed = 2 + sim*100)
    })

    gc()

    netcrop_status(sim, nsim, "NCV", time.ncv[3],
                   out.ncv$best.l2.each.rep, R.use)

    ncv.tab <- tibble::tibble(
      nsim = sim,
      model = model,
      n = n, K = K,
      out_in_ratio = out.in.ratio, alpha = alpha, lambda = lambda,
      max_K = max.K, loss_function = loss.use,
      cv = 3, R = R.use,
      best_model = do.call('c', out.ncv[paste0("best.", loss.use, ".stable")]),
      run_time = time.ncv[3],
      user_time = time.ncv[1],
      system_time = time.ncv[2]
    ) |>
      dplyr::mutate(
        model_hat = stringr::str_extract(best_model, "^[^-]+"),
        K_hat = as.integer(stringr::str_extract(best_model, "(?<=-)\\d+")),
        .before = run_time
      )

    readr::write_csv(
      ncv.tab, file = ncv.file, append = file.exists(ncv.file)
    )

    all.ncv[[count]] <- list(time = time.ncv, out = out.ncv)

    saveRDS(
      all.ncv[[count]], file = file.path(
        LOG_DIR,
        paste0("case2_ncv_v", version, "_sim", sim, "_R", R.use, ".rds")
      )
    )

    count <- count + 1
  }
}

################################################################################
################################################################################
################################################################################
################################################################################

## Run ECV
all.ecv <- list()
count <- 1
R <- c(1, 20)
run.ecv <- netcrop_confirm_large_cv(n, "ECV", small.script)
ecv.simulations <- netcrop_resume_csv(
  ecv.file, nsim, OUTPUT_ACTION, expected_rows = length(R),
  key_columns = c("nsim", "R", "loss_function")
)
if (run.ecv) for(sim in ecv.simulations){
  net <- SBM.gen(n = n, K = K, B = B, ncore = ncore, seed = sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    R.use <- R[RR]

    time.ecv <- system.time({
      out.ecv <- ECV.stability.BM(A = net$A, max.K = max.K, cv = 3,
                                  train.p = 0.9, R = R.use,
                                  tau = 0, dc.est = 2,
                                  loss = loss.use,
                                  ncore = ncore, seed = 2 + sim*100)
    })

    gc()

    netcrop_status(sim, nsim, "ECV", time.ecv[3],
                   out.ecv$best.l2.each.rep, R.use)

    ecv.tab <- tibble::tibble(
      nsim = sim,
      model = model,
      n = n, K = K,
      out_in_ratio = out.in.ratio, alpha = alpha, lambda = lambda,
      max_K = max.K, loss_function = loss.use,
      cv = 3, R = R.use,
      best_model = do.call('c', out.ecv[paste0("best.", loss.use, ".stable")]),
      run_time = time.ecv[3],
      user_time = time.ecv[1],
      system_time = time.ecv[2]
    ) |>
      dplyr::mutate(
        model_hat = stringr::str_extract(best_model, "^[^-]+"),
        K_hat = as.integer(stringr::str_extract(best_model, "(?<=-)\\d+")),
        .before = run_time
      )

    readr::write_csv(
      ecv.tab, file = ecv.file, append = file.exists(ecv.file)
    )

    all.ecv[[count]] <- list(time = time.ecv, out = out.ecv)

    saveRDS(
      all.ecv[[count]], file = file.path(
        LOG_DIR,
        paste0("case2_ecv_v", version, "_sim", sim, "_R", R.use, ".rds")
      )
    )

    count <- count + 1
  }
}

################################################################################
################################################################################
################################################################################
library(dplyr)

true_model <- paste0(model, "-", K)

nc.all <- readr::read_csv(nc.file, show_col_types = FALSE)

print(nc.all |> dplyr::group_by(s, o, R) |>
       dplyr::summarize(
         nsim = dplyr::n(),
         avg.deg = mean(lambda),
         mean.time = mean(run_time),
         accu = 100 * mean(best_model == true_model),
         mad = mean(abs(K_hat - K)),
         mean.Khat = mean(K_hat),
         sb.count = 100*mean(model_hat == "SBM"),
         dc.count = 100*mean(model_hat == "DCBM")
       ))

if (file.exists(ncv.file)) {
  ncv.all <- readr::read_csv(ncv.file, show_col_types = FALSE)
  print(ncv.all |> dplyr::group_by(R) |>
  dplyr::summarize(
    nsim = dplyr::n(),
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(best_model == true_model),
    mad = mean(abs(K_hat - K)),
    mean.Khat = mean(K_hat),
    sb.count = 100*mean(model_hat == "SBM"),
    dc.count = 100*mean(model_hat == "DCBM"),
    mad = mean(abs(K_hat - K))
  ))
}


if (file.exists(ecv.file)) {
  ecv.all <- readr::read_csv(ecv.file, show_col_types = FALSE)
  print(ecv.all |> dplyr::group_by(R) |>
  dplyr::summarize(
    nsim = dplyr::n(),
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(best_model == true_model),
    mad = mean(abs(K_hat - K)),
    mean.Khat = mean(K_hat),
    sb.count = 100*mean(model_hat == "SBM"),
    dc.count = 100*mean(model_hat == "DCBM")
  ))
}













