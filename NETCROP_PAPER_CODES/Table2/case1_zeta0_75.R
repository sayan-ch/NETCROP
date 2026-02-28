setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table2")))

if(!dir.exists(c("output")))
  dir.create("output")

if(!dir.exists("logs"))
  dir.create("logs")


HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "RDPG_helpers.R"))

################################################################################

version <- 1
ncore <- 40 # set the number of available processors to parallelize
nsim <- 100

# Case 1 of Table 2:
# RDPG with n = 10000 nodes and d = 10 dimension
# P = xi * XX^T / max(XX^T), xi = 0.75
n <- 10^4
d <- 10
xi <- 0.75

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing smallest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.d <- 20
loss.use <- c("l2")

## Run NETCROP
all.nc <- list()
count <- 1

for(sim in 1:nsim){
  net <- RDPG.gen(n = n, d = d, X = NULL, rho = xi,
                  ncore = ncore, seed = 200 + sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    for(ss in seq_along(s)){
      for(oo in seq_along(o)){
        s.use <- s[ss]
        o.use <- o[oo]
        R.use <- R[RR]

        time.nc <- system.time({
          out.nc <- netcrop_rdpg(A = net$A, d.cand = max.d,
                                       s = s.use, o = o.use, R = R.use,
                                       loss = loss.use,
                                       ncore = ncore, seed = 500 + sim*20*RR)
        })

        gc()
        cat("Sim ", sim, "::", "Time:", time.nc[3], ":: d_hat : ",
            paste(out.nc$`d.hat.each.rep (l2)`, collapse = ", "), "\n")

        nc.tab <- data.table::data.table(
          nsim = sim,
          model = "RDPG", n = n, d = d,
          xi = xi, lambda = lambda,
          max_d = max.d, loss_function = loss.use,
          s = s.use, o = o.use, R = R.use,
          d_hat = out.nc$l2.model,
          run_time = time.nc[3],
          user_time = time.nc[1],
          system_time = time.nc[2]
        )

        readr::write_csv(
          nc.tab, file = file.path(paste0("output/case1_netcrop_v", version, ".csv")),
          append = file.exists(paste0("output/case1_netcrop_v", version, ".csv"))
        )

        all.nc[[count]] <- list(time = time.nc, out = out.nc)

        saveRDS(
          all.nc[[count]], file = file.path(paste0("logs/case1_netcrop_v", version, ".rds"))
        )

        count <- count + 1
      }
    }
  }
}


## Run ECV
all.ecv <- list()
count <- 1
R <- c(1, 20)
for(sim in 1:nsim){
  net <- RDPG.gen(n = n, d = d, X = NULL, rho = xi,
                  ncore = ncore, seed = 200 + sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    R.use <- R[RR]

    time.ecv <- system.time({
      out.ecv <- ECV.stability.RDPG(A = net$A, max.K = max.d, cv = 3,
                                    R = R.use, train.p = 0.9,
                                    loss = loss.use,
                                    ncore = ncore, seed = 2 + sim*100)
    })

    gc()

    cat("Sim ", sim, "::", "Time:", time.ecv[3], "::",
        paste(out.ecv$best.l2.each.rep, collapse = ", "), "\n")

    ecv.tab <- data.table::data.table(
      nsim = sim,
      model = "RDPG", n = n, d = d,
      xi = xi, lambda = lambda,
      max_d = max.d, loss_function = loss.use,
      cv = 3, R = R.use,
      d_hat = out.ecv$best.l2.stable,
      run_time = time.ecv[3],
      user_time = time.ecv[1],
      system_time = time.ecv[2]
    )

    readr::write_csv(
      ecv.tab, file = file.path(paste0("output/case1_ecv_v", version, ".csv")),
      append = file.exists(paste0("output/case1_ecv_v", version, ".csv"))
    )

    all.ecv[[count]] <- list(time = time.ecv, out = out.ecv)

    saveRDS(
      all.ecv[[count]], file = file.path(paste0("logs/case1_ecv_v", version, ".rds"))
    )

    count <- count + 1
  }
}

################################################################################

nc.all <- readr::read_csv(file.path(paste0("output/case1_netcrop_v", version, ".csv")))

View(nc.all |> dplyr::group_by(s, o, R) |>
  dplyr::summarize(
    nsim = dplyr::n(),,
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mad = mean(abs(d_hat - d)),
    mean.dhat = mean(d_hat)
  ))

ecv.all <- readr::read_csv(file.path(paste0("output/case1_ecv_v", version, ".csv")))

View(ecv.all |> dplyr::group_by(s, o, R) |>
  dplyr::summarize(
    nsim = dplyr::n(),,
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mad = mean(abs(d_hat - d)),
    mean.dhat = mean(d_hat)
  ))












