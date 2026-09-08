script_dir <- local({
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  frames <- sys.frames(); ofiles <- vapply(frames, function(x) if (is.null(x$ofile)) "" else as.character(x$ofile), character(1))
  path <- if (length(arg)) sub("^--file=", "", arg[[1L]]) else if (any(nzchar(ofiles))) tail(ofiles[nzchar(ofiles)], 1L) else ""
  if (nzchar(path)) dirname(normalizePath(path)) else normalizePath(getwd())
})
if (!requireNamespace("netOP", quietly=TRUE)) stop("Install netOP v0.1.1 first.")
if (packageVersion("netOP") != "0.1.1") stop("This script requires netOP v0.1.1.")

run_name <- "xx_small_network_test"; output_dir <- file.path(script_dir,"output",run_name); dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
results_file <- file.path(output_dir,"small_results.rds"); csv_file <- file.path(output_dir,"small_results.csv")
detected_cores <- parallel::detectCores(); ncores <- if(is.na(detected_cores)) 1L else min(2L,detected_cores)
nsim <- 10L; n <- 500L; d <- 3L; xi <- 1; d_candidates <- 1:5; losses <- "sse"
partition <- netOP::netcrop_param_select(test_prop=0.1,n=n,o_range=0); num_subnetworks <- partition$num_subnetworks[[1L]]; overlap_size <- partition$overlap_size[[1L]]
one_simulation <- function(simulation) {
  A <- netOP::generate_rdpg(n=n,d=d,sparsity_multiplier=xi,ncores=ncores); lambda <- mean(Matrix::rowSums(A)); rows <- list(); fits <- list(); k <- 1L
  for(nrep in c(1L,5L)) {
    timing <- system.time(fit <- netOP::netcrop_rdpg(A,d_candidates,num_subnetworks,overlap_size,nrep,losses,ncores=ncores,verbose=FALSE))
    selected <- fit$overall_best$d_hat[fit$overall_best$loss==losses]
    message(sprintf("Simulation %d/%d | NETCROP R=%d | per-repetition d_hat=%s | elapsed=%.3f s",simulation,nsim,nrep,paste(fit$best_dimension_cv$d_hat,collapse=","),timing[["elapsed"]]))
    rows[[k]] <- data.frame(simulation,algorithm="NETCROP",model="RDPG",n,d,xi,lambda,max_d=max(d_candidates),loss_function=losses,s=num_subnetworks,o=overlap_size,cv=NA_integer_,R=nrep,d_hat=selected,run_time=timing[["elapsed"]])
    fits[[paste0("netcrop_R",nrep)]] <- fit; k <- k+1L
  }
  for(nrep in c(1L,20L)) {
    timing <- system.time(fit <- netOP::ecv_stability_rdpg(A,max(d_candidates),3L,nrep,0.9,losses,ncores=ncores,verbose=FALSE))
    selected <- fit$overall_best$d_hat[fit$overall_best$loss==losses]
    message(sprintf("Simulation %d/%d | ECV R=%d | per-repetition d_hat=%s | elapsed=%.3f s",simulation,nsim,nrep,paste(fit$best_dimension_cv$d_hat,collapse=","),timing[["elapsed"]]))
    rows[[k]] <- data.frame(simulation,algorithm="ECV",model="RDPG",n,d,xi,lambda,max_d=max(d_candidates),loss_function=losses,s=NA_integer_,o=NA_integer_,cv=3L,R=nrep,d_hat=selected,run_time=timing[["elapsed"]])
    fits[[paste0("ecv_R",nrep)]] <- fit; k <- k+1L
  }
  list(summary=do.call(rbind,rows),fits=fits)
}
records <- netOP::run_simulations(one_simulation,nsim=nsim,results_file=results_file,action="resume",show_progress=TRUE,continue_on_error=TRUE)
successful <- Filter(function(x)isTRUE(x$success),records)
if(length(successful)){results <- do.call(rbind,lapply(successful,function(x)x$result$summary));write.csv(results,csv_file,row.names=FALSE);print(aggregate(cbind(lambda,run_time,correct=as.numeric(d_hat==d))~algorithm+R,results,function(x)mean(x,na.rm=TRUE)))}
failed <- Filter(function(x)identical(x$success,FALSE),records);if(length(failed))warning(length(failed)," simulation(s) failed; inspect ",results_file)
message("Small Table 2 test completed. Results: ",normalizePath(output_dir))
