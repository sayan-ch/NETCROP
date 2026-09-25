setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table3")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "LSM_helpers.R"))

RUN_NAME <- "case1_d2_a0"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################

version <- 1
detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores) | .Platform$OS.type == "windows" ) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 100

# Case 1 of Table 3:
# LSM with n = 1000, d = 2, alpha = 0
n <- 10^3
d <- 2
alpha <- 0

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing smallest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.d <- 5
loss.use <- c("l2")
nc.file <- file.path(OUTPUT_DIR, paste0("case1_netcrop_v", version, ".csv"))

## Run NETCROP
all.nc <- list()
count <- 1
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION,
  expected_rows = length(R) * length(s) * length(o),
  key_columns = c("nsim", "R", "s", "o", "loss_function")
)
for(sim in nc.simulations){
  net <- LSM.gen(n = n, d = d, alpha = alpha, ncore = ncore,
                 seed = 100 + sim)

  lambda <- mean(rowSums(net$A))

  for(RR in seq_along(R)){
    for(ss in seq_along(s)){
      for(oo in seq_along(o)){
        s.use <- s[ss]
        o.use <- o[oo]
        R.use <- R[RR]

        time.nc <- system.time({
          out.nc <- netcrop_lsm(A = net$A, d.cand = max.d,
                                s = s.use, o = o.use, R = R.use,
                                loss = loss.use,
                                ncore = ncore, seed = 500 + sim,
                                step.size = 0.3, niter = 100, trace = 0)
        })

        gc()
        netcrop_status(sim, nsim, "NETCROP", time.nc[3],
                       out.nc$`d.hat.each.rep (l2)`, R.use)

        nc.tab <- tibble::tibble(
          nsim = sim,
          model = "LSM", n = n, d = d,
          alpha = alpha, lambda = lambda,
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
nc.all <- readr::read_csv(nc.file, show_col_types = FALSE)


print(nc.all |> dplyr::group_by(s, o, R) |>
  dplyr::summarize(
    nsim = dplyr::n(),
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mean.dhat = mean(d_hat),
    mad = mean(abs(d_hat - d))
  ))
