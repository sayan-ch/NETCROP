setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table2")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "RDPG_helpers.R"))

RUN_NAME <- "case1_zeta0_75"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################
################################################################################

version <- 1
detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores)) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 100L

# Case 1 of Table 2:
# RDPG with n = 10000 nodes and d = 10 dimension
# P = xi * XX^T / max(XX^T), xi = 0.75
n <- 10^4
d <- 10
xi <- 0.75

################################################################################
################################################################################
################################################################################
################################################################################

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing smallest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.d <- 20
loss.use <- c("l2")

nc.file <- file.path(OUTPUT_DIR, paste0("case1_netcrop_v", version, ".csv"))
ecv.file <- file.path(OUTPUT_DIR, paste0("case1_ecv_v", version, ".csv"))
small.script <- file.path(here::here("NETCROP_PAPER_CODES", "Table2"),
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
        netcrop_status(sim, nsim, "NETCROP", time.nc[3],
                       out.nc$`d.hat.each.rep (l2)`, R.use)

        nc.tab <- tibble::tibble(
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
          nc.tab, file = nc.file, append = file.exists(nc.file)
        )

        all.nc[[count]] <- list(time = time.nc, out = out.nc)

        saveRDS(
          all.nc[[count]], file = file.path(
            LOG_DIR,
            paste0("case1_netcrop_v", version, "_sim", sim, "_R", R.use, ".rds")
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

    netcrop_status(sim, nsim, "ECV", time.ecv[3],
                   out.ecv$best.l2.each.rep, R.use)

    ecv.tab <- tibble::tibble(
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
      ecv.tab, file = ecv.file, append = file.exists(ecv.file)
    )

    all.ecv[[count]] <- list(time = time.ecv, out = out.ecv)

    saveRDS(
      all.ecv[[count]], file = file.path(
        LOG_DIR,
        paste0("case1_ecv_v", version, "_sim", sim, "_R", R.use, ".rds")
      )
    )

    count <- count + 1
  }
}

################################################################################

nc.all <- readr::read_csv(nc.file, show_col_types = FALSE)

print(nc.all |> dplyr::group_by(s, o, R) |>
  dplyr::summarize(
    nsim = dplyr::n(),
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mad = mean(abs(d_hat - d)),
    mean.dhat = mean(d_hat)
  ))

if (file.exists(ecv.file)) {
  ecv.all <- readr::read_csv(ecv.file, show_col_types = FALSE)
  print(ecv.all |> dplyr::group_by(R) |>
  dplyr::summarize(
    nsim = dplyr::n(),
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mad = mean(abs(d_hat - d)),
    mean.dhat = mean(d_hat)
  ))
}











