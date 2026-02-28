setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata")))

if(!dir.exists(c("output")))
  dir.create("output")

if(!dir.exists("logs"))
  dir.create("logs")


HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

################################################################################
ncore <- 80
nsim <- 100

library(tidyverse)

## Read from csv edge list (already symmetric)
twitch.el <- read_csv(file.path("Twitch", "Twitch_edge_list.csv"))
n <- max(twitch.el) 
twitch.A <- sparseMatrix(i = twitch.el$i, j = twitch.el$j, x = 1,
                         dims = c(n,n))

param.out <- netcrop_param(p.test = 0.02, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

for(ii in 1:nsim){
  for(rr in c(1, 5)){
    time.nc <- system.time({
      out.nc <- netcrop_blockmodel(A = twitch.A, K.CAND = 30,
                                   s = s, o = o, R = rr, tau = 0,
                                   laplace = T, loss = c('l2', 'AUC'),
                                   mod.cand = c('SBM', 'DCBM'), ncore = ncore,
                                   seed = 100 + 200*ii)
    })[3]
    
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
    
    write_csv(tab.nc, file.path("output/Twitch_netcrop.csv"),
              append = file.exists(file.path("output/Twitch_netcrop.csv")))
  }
}







