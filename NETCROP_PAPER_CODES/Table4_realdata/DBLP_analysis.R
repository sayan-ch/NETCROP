setwd(file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

RUN_NAME <- "DBLP_analysis"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################
detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores)) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 100L

library(dplyr)
library(readr)

## Read from csv edge list (already symmetric)
DBLP.el <- read_csv(file.path("DBLP", "DBLP_conf_edge_list.csv"))
n <- max(DBLP.el)
DBLP.A <- sparseMatrix(i = DBLP.el$i, j = DBLP.el$j, x = 1,
                         dims = c(n,n))

param.out <- netcrop_param(p.test = 0.02, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

nc.file <- file.path(OUTPUT_DIR, "DBLP_netcrop.csv")
ncv.file <- file.path(OUTPUT_DIR, "DBLP_ncv.csv")
ecv.file <- file.path(OUTPUT_DIR, "DBLP_ecv.csv")
small.script <- file.path(here::here("NETCROP_PAPER_CODES", "Table4_realdata"),
                          "xx_small_network_test.R")
netcrop.R <- c(1, 5, 10, 20)
cv.R <- c(1, 20)

nc.simulations <- netcrop_resume_csv(
  nc.file, nsim, OUTPUT_ACTION, expected_rows = length(netcrop.R),
  key_columns = c("nsim", "R")
)
for(ii in nc.simulations){
  for(rr in netcrop.R){
    time.nc <- system.time({
      out.nc <- netcrop_blockmodel(A = DBLP.A, K.CAND = 10,
                                   s = s, o = o, R = rr, tau = 0,
                                   laplace = F, loss = c('l2', 'AUC'),
                                   mod.cand = c('SBM', 'DCBM'), ncore = ncore,
                                   seed = 100 + 200*ii)
    })[3]

    netcrop_status(ii, nsim, "NETCROP", time.nc,
                   c(out.nc$l2.model, out.nc$AUC.model), rr)

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

    write_csv(tab.nc, nc.file, append = file.exists(nc.file))
  }
}

###########################
# Run NCV
run.ncv <- netcrop_confirm_large_cv(n, "NCV", small.script)
ncv.simulations <- netcrop_resume_csv(
  ncv.file, nsim, OUTPUT_ACTION, expected_rows = length(cv.R),
  key_columns = c("nsim", "R", "loss_function")
)
if (run.ncv) for(ii in ncv.simulations){
  for(rr in cv.R){
    time.ncv <- system.time({
      out.ncv <- NCV.stability.BM(A = DBLP.A, max.K = 10, cv = 3, R = rr,
                                   tau = 0, laplace = F, loss = c('l2', 'AUC'),
                                   ncore = ncore, seed = 100 + 200*ii)
    })[3]

    netcrop_status(ii, nsim, "NCV", time.ncv,
                   out.ncv$best.l2.each.rep, rr)

    best.model.ncv <- do.call('c', out.ncv[paste0("best.l2.stable")])
    model.hat.ncv <- stringr::str_extract(best.model.ncv, "^[^-]+")
    K.hat.ncv <- as.integer(stringr::str_extract(best.model.ncv, "(?<=-)\\d+"))

    test.auc.ncv <- -mean(out.ncv$ncv.loss[[paste0(tolower(model.hat.ncv), ".auc")]][
      seq(K.hat.ncv, by = 10, length.out = rr)])

    tab.ncv <- tibble::tibble(
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

    write_csv(tab.ncv, ncv.file, append = file.exists(ncv.file))
  }
}

###########################
# Run ECV
run.ecv <- netcrop_confirm_large_cv(n, "ECV", small.script)
ecv.simulations <- netcrop_resume_csv(
  ecv.file, nsim, OUTPUT_ACTION, expected_rows = length(cv.R),
  key_columns = c("nsim", "R", "loss_function")
)
if (run.ecv) for(ii in ecv.simulations){
  for(rr in cv.R){
    time.ecv <- system.time({
      out.ecv <- ECV.stability.BM(A = DBLP.A, max.K = 10, cv = 3, R = rr,
                                  tau = 0, loss = c('l2', 'AUC'),
                                  ncore = ncore, seed = 100 + 200*ii)
    })[3]

    netcrop_status(ii, nsim, "ECV", time.ecv,
                   out.ecv$best.l2.each.rep, rr)

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

    tab.ecv <- tibble::tibble(
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

    write_csv(tab.ecv, ecv.file, append = file.exists(ecv.file))
  }
}


################################################################################
out.nc <- read_csv(nc.file, show_col_types = FALSE)

print(out.nc |> group_by(s, o, R) |>
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_l2 == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(time),
         .groups = "drop"))


if (file.exists(ncv.file)) {
  out.ncv <- read_csv(ncv.file, show_col_types = FALSE)
  print(out.ncv |> group_by(R) |>
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_model == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(run_time),
         .groups = "drop"))
}


if (file.exists(ecv.file)) {
  out.ecv <- read_csv(ecv.file, show_col_types = FALSE)
  print(out.ecv |> group_by(R) |>
       summarize(
         nsim = n(),
         accuracy = 100*mean(best_model == "DCBM-4"),
         test_auc = mean(test_auc_best_l2),
         time = mean(run_time),
         .groups = "drop"))
}








