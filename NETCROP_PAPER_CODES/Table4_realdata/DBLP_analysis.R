setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata")))

if(!dir.exists(c("output")))
  dir.create("output")

if(!dir.exists("logs"))
  dir.create("logs")


HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

################################################################################
ncore <- 40
nsim <- 100

library(tidyverse)

## Read from csv edge list (already symmetric)
DBLP.el <- read_csv(file.path("DBLP", "DBLP_conf_edge_list.csv"))
n <- max(DBLP.el) 
DBLP.A <- sparseMatrix(i = DBLP.el$i, j = DBLP.el$j, x = 1,
                         dims = c(n,n))

param.out <- netcrop_param(p.test = 0.02, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

for(ii in 1:nsim){
  for(rr in c(1, 5, 10, 20)){
    time.nc <- system.time({
      out.nc <- netcrop_blockmodel(A = DBLP.A, K.CAND = 10,
                                   s = s, o = o, R = rr, tau = 0,
                                   laplace = F, loss = c('l2', 'AUC'),
                                   mod.cand = c('SBM', 'DCBM'), ncore = ncore,
                                   seed = 100 + 200*ii)
    })[3]
    
    test.auc <- out.nc$loss |> group_by(Candidate_Model, Candidate_Value) |> 
      summarize(test.auc = mean(AUC), .groups = "drop")
    
    tab.nc <- data.frame(
      nsim = ii, data = "DBLP-Conf", algo = "NETCROP",
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
    
    write_csv(tab.nc, file.path("output/DBLP_netcrop.csv"),
              append = file.exists(file.path("output/DBLP_netcrop.csv")))
  }
}

###########################
# Run ECV and NCV
for(ii in 1:nsim){
  for(rr in c(1, 20)){
    time.ncv <- system.time({
      out.ncv <- NCV.stability.BM(A = DBLP.A, max.K = 10, cv = 3, R = rr,
                                   tau = 0, laplace = F, loss = c('l2', 'AUC'),
                                   ncore = ncore, seed = 100 + 200*ii)
    })[3]
    
    best.model.ncv <- do.call('c', out.ncv[paste0("best.l2.stable")])
    model.hat.ncv <- stringr::str_extract(best.model.ncv, "^[^-]+")
    K.hat.ncv <- as.integer(stringr::str_extract(best.model.ncv, "(?<=-)\\d+"))
    
    test.auc.ncv <- -mean(out.ncv$ncv.loss[[paste0(tolower(model.hat.ncv), ".auc")]][
      seq(K.hat.ncv, by = 10, length.out = rr)])
    
    tab.ncv <- data.table::data.table(
      nsim = ii, data = "DBLP", n = n,
      algo = "NCV",
      max_K = 10, loss_function = 'l2',
      cv = 3, R = rr,
      best_model = do.call('c', out.ncv[paste0("best.l2.stable")]),
      test_auc_best_l2 = test.auc.ncv,
      run_time = time.ncv
    ) |>
      dplyr::mutate(
        model_hat = stringr::str_extract(best_model, "^[^-]+"),
        K_hat = as.integer(stringr::str_extract(best_model, "(?<=-)\\d+")),
        .before = run_time
      )
    
    write_csv(tab.ncv, file.path("output/DBLP_ncv.csv"),
              append = file.exists(file.path("output/DBLP_ncv.csv")))
    
    time.ecv <- system.time({
      out.ecv <- ECV.stability.BM(A = DBLP.A, max.K = 10, cv = 3, R = rr,
                                  tau = 0, loss = c('l2', 'AUC'),
                                  ncore = ncore, seed = 100 + 200*ii)
    })[3]
    
    best.model.ecv <- do.call('c', out.ecv[paste0("best.l2.stable")])
    model.hat.ecv <- stringr::str_extract(best.model.ecv, "^[^-]+")
    K.hat.ecv <- as.integer(stringr::str_extract(best.model.ecv, "(?<=-)\\d+"))
    
    if(model.hat.ecv == "SBM"){
      test.auc.ecv <- - mean(sapply(1:rr, function(ll){
        mean(out.ecv$ecv.loss[[ll]]$sbm.auc.mat[, K.hat.ecv])}))
    } else if(model.hat.ecv == "DCBM"){
      test.auc.ecv <- - mean(sapply(1:rr, function(ll){
        mean(out.ecv$ecv.loss[[ll]]$dc.auc.mat[, K.hat.ecv])}))
    }
    
    tab.ecv <- data.table::data.table(
      nsim = ii, data = "DBLP", n = n,
      algo = "ecv",
      max_K = 10, loss_function = 'l2',
      cv = 3, R = rr,
      best_model = do.call('c', out.ecv[paste0("best.l2.stable")]),
      test_auc_best_l2 = test.auc.ecv,
      run_time = time.ecv
    ) |>
      dplyr::mutate(
        model_hat = stringr::str_extract(best_model, "^[^-]+"),
        K_hat = as.integer(stringr::str_extract(best_model, "(?<=-)\\d+")),
        .before = run_time
      )
    
    write_csv(tab.ecv, file.path("output/DBLP_ecv.csv"),
              append = file.exists(file.path("output/DBLP_ecv.csv")))
  }
}


################################################################################
out.nc <- read_csv(file.path("output/DBLP_netcrop.csv"))

View(out.nc |> group_by(s, o, R) |> 
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_l2 == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(time),
         .groups = "drop"))


out.ncv <- read_csv(file.path("output/DBLP_ncv.csv"))

View(out.ncv |> group_by(s, o, R) |> 
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_model == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(run_time),
         .groups = "drop"))


out.ecv <- read_csv(file.path("output/DBLP_ecv.csv"))

View(out.ecv |> group_by(s, o, R) |> 
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_model == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(run_time),
         .groups = "drop"))









