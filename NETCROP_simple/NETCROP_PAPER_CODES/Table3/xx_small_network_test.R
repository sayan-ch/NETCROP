script_dir <- local({
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  frames <- sys.frames(); ofiles <- vapply(frames, function(x) if (is.null(x$ofile)) "" else as.character(x$ofile), character(1))
  path <- if (length(arg)) sub("^--file=", "", arg[[1L]]) else if (any(nzchar(ofiles))) tail(ofiles[nzchar(ofiles)], 1L) else ""
  if (nzchar(path)) dirname(normalizePath(path)) else normalizePath(getwd())
})
if (!requireNamespace("netOP",quietly=TRUE)) stop("Install netOP v0.1.1 first.")
if (packageVersion("netOP")!="0.1.1") stop("This script requires netOP v0.1.1.")
run_name <- "xx_small_network_test"; output_dir <- file.path(script_dir,"output",run_name); dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
results_file <- file.path(output_dir,"small_results.rds"); csv_file <- file.path(output_dir,"small_results.csv")
detected_cores <- parallel::detectCores(); ncores <- if(is.na(detected_cores)) 1L else min(2L,detected_cores)
nsim <- 10L; n <- 500L; d <- 3L; alpha <- 0; d_candidates <- 1:5; losses <- "sse"
partition <- netOP::netcrop_param_select(test_prop=0.1,n=n,o_range=0); num_subnetworks <- partition$num_subnetworks[[1L]]; overlap_size <- partition$overlap_size[[1L]]
one_simulation <- function(simulation) {
  A <- netOP::generate_lsm(n=n,d=d,alpha=alpha,ncores=ncores); lambda <- mean(Matrix::rowSums(A)); rows <- list(); fits <- list()
  for(i in seq_along(c(1L,5L))) {
    nrep <- c(1L,5L)[[i]]
    timing <- system.time(fit <- netOP::netcrop_lsm(A,d_candidates,num_subnetworks,overlap_size,nrep,losses,lsm_options=list(step_size=0.3,niter=100L,trace=FALSE),ncores=ncores,verbose=FALSE))
    selected <- fit$overall_best$d_hat[fit$overall_best$loss==losses]
    message(sprintf("Simulation %d/%d | NETCROP R=%d | per-repetition d_hat=%s | elapsed=%.3f s",simulation,nsim,nrep,paste(fit$best_dimension_cv$d_hat,collapse=","),timing[["elapsed"]]))
    rows[[i]] <- data.frame(simulation,model="LSM",n,d,alpha,lambda,max_d=max(d_candidates),loss_function=losses,s=num_subnetworks,o=overlap_size,R=nrep,d_hat=selected,run_time=timing[["elapsed"]])
    fits[[paste0("netcrop_R",nrep)]] <- fit
  }
  list(summary=do.call(rbind,rows),fits=fits)
}
records <- netOP::run_simulations(one_simulation,nsim=nsim,results_file=results_file,action="resume",show_progress=TRUE,continue_on_error=TRUE)
successful <- Filter(function(x)isTRUE(x$success),records)
if(length(successful)){results <- do.call(rbind,lapply(successful,function(x)x$result$summary));write.csv(results,csv_file,row.names=FALSE);print(aggregate(cbind(lambda,run_time,correct=as.numeric(d_hat==d))~s+o+R,results,function(x)mean(x,na.rm=TRUE)))}
failed <- Filter(function(x)identical(x$success,FALSE),records);if(length(failed))warning(length(failed)," simulation(s) failed; inspect ",results_file)
message("Small Table 3 test completed. Results: ",normalizePath(output_dir))
