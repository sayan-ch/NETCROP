setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table3")))

if(!dir.exists(c("output")))
  dir.create("output")

if(!dir.exists("logs"))
  dir.create("logs")


HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "LSM_helpers.R"))

################################################################################

version <- 1
ncore <- 40 # set the number of available processors to parallelize
nsim <- 100

# Case 3 of Table 3:
# RDPG with n = 1000, d = 5, alpha = 0
n <- 10^3
d <- 5
alpha <- 0

# NETCROP parameters
## Parameter selection
p.test <- 0.02
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0) # choosing smallest o in feasible region

s <- param.out$s
o <- param.out$o
R <- c(1, 5)
max.d <- 10
loss.use <- c("l2")

## Run NETCROP
all.nc <- list()
count <- 1

for(sim in 77:nsim){
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
        cat("Sim ", sim, "::", "Time:", time.nc[3], ":: d_hat : ",
            paste(out.nc$`d.hat.each.rep (l2)`, collapse = ", "), "\n")
        
        nc.tab <- data.table::data.table(
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
          nc.tab, file = file.path(paste0("output/case3_netcrop_v", version, ".csv")),
          append = file.exists(paste0("output/case3_netcrop_v", version, ".csv"))
        )
        
        all.nc[[count]] <- list(time = time.nc, out = out.nc)
        
        saveRDS(
          all.nc[[count]], file = file.path(paste0("logs/case3_netcrop_v", version, ".rds"))
        )
        
        count <- count + 1
      }
    }
  }
}

################################################################################
nc.all <- readr::read_csv(file.path(paste0("output/case3_netcrop_v", version, ".csv")))


View(nc.all |> dplyr::group_by(s, o, R) |>
  dplyr::summarize(
    nsim = dplyr::n(),,
    avg.deg = mean(lambda),
    mean.time = mean(run_time),
    accu = 100 * mean(d_hat == d),
    mean.dhat = mean(d_hat),
    mad = mean(abs(d_hat - d))
  ))

