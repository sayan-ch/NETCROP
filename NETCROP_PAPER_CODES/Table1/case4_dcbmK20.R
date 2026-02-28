setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table1")))

if(!dir.exists(c("output")))
  dir.create("output")

if(!dir.exists("logs"))
  dir.create("logs")


HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

################################################################################

version <- 1
ncore <- 40 # set the number of available processors to parallelize
nsim <- 100

# Case 4 of Table 1:
# DCBM with n = 10000 nodes and K = 20 communities
# out.in.ratio = 1/3 and alpha = 3
n <- 10^4
K <- 20
out.in.ratio <- 1/3
alpha <- 3
B <- alpha * (diag(1-out.in.ratio, K) + out.in.ratio)

PI <- rep(1/K, K)

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing lowest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.K <- 30
loss.use <- c("l2")

## Run NETCROP
all.nc <- list()
count <- 1


for(sim in 1:nsim){
  net <- DCBM.gen(n = n, K = K, avg.deg = 660, beta = out.in.ratio,
                  ncore = ncore, seed = 200 + sim)

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
                                       loss = loss.use,
                                       mod.cand = c('SBM', 'DCBM'),
                                       ncore = ncore, 
                                       seed = 500 + sim*20*RR)
        })

        gc()
        
        cat("Sim ", sim, "::", "Time:", time.nc[3], "::", 
            paste(out.nc$`Mod.K.hat.each.rep (l2)`,
                  collapse = ", "), "\n")

        nc.tab <- data.table::data.table(
          nsim = sim,
          model = dplyr::if_else(all(psi == 1), "SBM", "DCBM"), n = n, K = K,
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
          nc.tab, file = file.path(paste0("output/case4_netcrop_v", version, ".csv")),
          append = file.exists(paste0("output/case4_netcrop_v", version, ".csv"))
        )

        all.nc[[count]] <- list(time = time.nc, out = out.nc)

        saveRDS(
          all.nc[[count]], file = file.path(paste0("logs/case4_netcrop_v", version, ".rds"))
        )

        count <- count + 1
      }
    }
  }
}



## Run NCV
all.ncv <- list()
count <- 1
R <- c(1, 20)
for(sim in 1:nsim){
  net <- DCBM.gen(n = n, K = K, avg.deg = 660, beta = out.in.ratio,
                  ncore = ncore, seed = 200 + sim)

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

    cat("Sim ", sim, "::", "Time:", time.ncv[3], "::", 
        paste(out.ncv$best.l2.each.rep, collapse = ", "), "\n")

    ncv.tab <- data.table::data.table(
      nsim = sim,
      model = dplyr::if_else(all(psi == 1), "SBM", "DCBM"), n = n, K = K,
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
      ncv.tab, file = file.path(paste0("output/case4_ncv_v", version, ".csv")),
      append = file.exists(paste0("output/case4_ncv_v", version, ".csv"))
    )

    all.ncv[[count]] <- list(time = time.ncv, out = out.ncv)

    saveRDS(
      all.ncv[[count]], file = file.path(paste0("logs/case4_ncv_v", version, ".rds"))
    )

    count <- count + 1
  }
}

## Run ECV
all.ecv <- list()
count <- 1
R <- c(1, 20)
for(sim in 1:nsim){
  set.seed(100 + sim)
  g <- sample(1:K, size = n, replace = T, prob = PI)

  psi.pool <- rbeta(300, 1, 4)
  psi <- rep(1, n)
  wh.dc <- sample(1:n, 0.2 * n, F)
  psi[wh.dc] <- sample(psi.pool, length(wh.dc), T)

  net <- DCBM.gen(n = n, K = K, avg.deg = 660, beta = out.in.ratio,
                  ncore = ncore, seed = 200 + sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    R.use <- R[RR]

    time.ecv <- system.time({
      out.ecv <- ECV.stability.BM(A = net$A, max.K = max.K, train.p = 0.9,
                                  cv = 3, R = R.use,
                                  tau = 0, dc.est = 2,
                                  loss = loss.use,
                                  ncore = ncore, seed = 2 + sim*100)
    })

    gc()

    cat("Sim ", sim, "::", "Time:", time.ecv[3], "::", 
        paste(out.ecv$best.l2.each.rep, collapse = ", "), "\n")

    ecv.tab <- data.table::data.table(
      nsim = sim,
      model = dplyr::if_else(all(psi == 1), "SBM", "DCBM"), n = n, K = K,
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
      ecv.tab, file = file.path(paste0("output/case4_ecv_v", version, ".csv")),
      append = file.exists(paste0("output/case4_ecv_v", version, ".csv"))
    )

    all.ecv[[count]] <- list(time = time.ecv, out = out.ecv)

    saveRDS(
      all.ecv[[count]], file = file.path(paste0("logs/case4_ecv_v", version, ".rds"))
    )

    count <- count + 1
  }
}

################################################################################
library(tidyverse)

nc.all <- readr::read_csv(file.path(paste0("output/case4_netcrop_v", version, ".csv")))

View(nc.all |> dplyr::group_by(s, o, R) |>
       dplyr::summarize(
         nsim = dplyr::n(),,
         avg.deg = mean(lambda),
         mean.time = mean(run_time),
         accu = 100 * mean(best_model == "DCBM-20"),
         mad = mean(abs(K_hat - K)),
         mean.Khat = mean(K_hat),
         sb.count = 100*mean(model_hat == "SBM"),
         dc.count = 100*mean(model_hat == "DCBM")
       ))

ncv.all <- readr::read_csv(file.path(paste0("output/case4_ncv_v", version, ".csv")))

ncv.all |> dplyr::group_by(R) |>
  dplyr::summarize(
    nsim = dplyr::n(),,
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(best_model == "DCBM-20"),
    mad = mean(abs(K_hat - K)),
    mean.Khat = mean(K_hat),
    sb.count = 100*mean(model_hat == "SBM"),
    dc.count = 100*mean(model_hat == "DCBM"),
    mad = mean(abs(K_hat - K))
  )


ecv.all <- readr::read_csv(file.path(paste0("output/case4_ecv_v", version, ".csv")))

ecv.all |> dplyr::group_by(R) |>
  dplyr::summarize(
    nsim = dplyr::n(),,
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(best_model == "DCBM-20"),
    mad = mean(abs(K_hat - K)),
    mean.Khat = mean(K_hat),
    sb.count = 100*mean(model_hat == "SBM"),
    dc.count = 100*mean(model_hat == "DCBM")
  )






