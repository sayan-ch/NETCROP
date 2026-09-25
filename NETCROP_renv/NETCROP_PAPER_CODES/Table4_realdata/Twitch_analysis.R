setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

RUN_NAME <- "Twitch_analysis"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################
detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores) | .Platform$OS.type == "windows" ) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 100L

library(dplyr)
library(tidyr)
library(readr)

## Read from csv edge list (already symmetric)
twitch.el <- read_csv(file.path("Twitch", "Twitch_edge_list.csv"),
                      show_col_types = FALSE)
n <- max(twitch.el)
twitch.A <- sparseMatrix(i = twitch.el$i, j = twitch.el$j, x = 1,
                         dims = c(n,n))

param.out <- netcrop_param(p.test = 0.02, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

nc.file <- file.path(OUTPUT_DIR, "Twitch_netcrop.csv")
netcrop.R <- 1L
nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION, expected_rows = length(netcrop.R),
  key_columns = c("nsim", "R")
)

for(ii in nc.simulations){
  for(rr in netcrop.R){

    time.nc <- system.time({
      out.nc <- netcrop_blockmodel(
        A = twitch.A,
        K.CAND = 30,
        s = s, o = o,
        R = rr, tau = 0,
        laplace = T, loss = c('l2', 'AUC'),
        mod.cand = c('SBM', 'DCBM'), ncore = ncore,
        seed = 100 + 200*ii,
        # clara parameters
        samples = 50L,
        sampsize = 200L
      )
    })[3]

    netcrop_status(ii, nsim, "NETCROP", time.nc,
                   c(out.nc$l2.model, out.nc$AUC.model), rr)

    test.auc <- out.nc$loss |> group_by(Candidate_Model, Candidate_Value) |>
      summarize(test.auc = mean(AUC), .groups = "drop")

    tab.nc <- data.frame(
      nsim = ii, data = "Twitch", algo = "NETCROP",
      s = s, o = o, R = rr,
      best_l2 = out.nc$l2.model, best_AUC = out.nc$AUC.model,
      test_auc_best_l2 = -mean(test.auc$test.auc[
        which(paste0(test.auc$Candidate_Model, "-", test.auc$Candidate_Value) ==
                out.nc$l2.model)]),
      test_auc_best_auc = -mean(test.auc$test.auc[
        which(paste0(test.auc$Candidate_Model, "-", test.auc$Candidate_Value) ==
                out.nc$AUC.model)]),
      time = time.nc
    )

    write_csv(tab.nc, nc.file, append = file.exists(nc.file))
  }

}






