library(Matrix)
library(parallel)

################################################################################
################################################################################
## RDPG generation
# n: number of nodes, d: dimension of latent positions,
# X: latent position matrix (n by d), if not supplied, generated from Unif(0,1) and scaled
# rho: sparsity factor (zeta in paper)
# ncore: number of cores for parallelization, seed: random seed for reproducibility
RDPG.gen <- function(n, d, X = NULL, rho = 1,
                     ncore = 1, seed = NULL){
  on.exit(gc())
  set.seed(seed)

  if(is.null(X)){
    X <- matrix(runif(n*d, 0, 1), nrow = n, ncol = d)
  }
  # P0 <- X %*% t(X)
  P0 <- tcrossprod(X)
  P <- rho * P0 / max(P0)
  diag(P) <- 0

  stor <- do.call('rbind',
                  mclapply(1:(n-1), function(i) {
                    if(!is.null(seed))
                      set.seed(seed + i)

                    tmp <- which(rbinom(n-i, 1, P[i, (i+1):n]) == 1)

                    if(length(tmp) == 0)
                      return(NULL)
                    else
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))

  A <- sparseMatrix(i = stor[,1], j = stor[,2], x = 1, dims = c(n,n))
  A <- A + t(A)

  return(list('A' = A, 'P' = P))
}
################################################################################
################################################################################
## NETCROP for RDPG
# A: adjacency matrix, d.cand: candidate dimensions for latent positions,
# s: number of splits, o: number of overlapping nodes, R: number of repetitions,
# loss: loss functions to evaluate,
# ncore: number of cores for parallelization, seed: random seed for reproducibility
netcrop_rdpg <- function(A, d.cand, s, o, R,
                         loss = c("l2", "bin.dev", "AUC"),
                         ncore = 1, seed = NULL){

  set.seed(seed)

  n <- nrow(A)
  m <- (n-o)/s

  A0 <- A
  rho.hat <- mean(rowSums(A0)) / (n-1)
  A <- A0 / rho.hat

  max.d <- max(d.cand)

  if(length(d.cand) == 1)
    d.cand <- 1:max.d

  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o,
                                              replace = F))

  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  colnames(raw.ind) <- c('r', 's')

  raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    q <- raw.ind[ii, 's']
    r <- raw.ind[ii, 'r']

    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]

    eig <- RSpectra::eigs_sym(A.sonn, k = max.d)

    U <- eig$vectors

    sigma.half <- Diagonal(n = max.d, x = sqrt(abs(eig$values)))

    list(U = U, sigma.half = sigma.half)

  }, mc.cores = min(ncore, nrow(raw.ind)))

  ld <- length(d.cand)
  match.ind <- matrix(nrow = R*s*ld, ncol = 3)
  cc <- 1
  for(r in 1:R)
    for(p in 1:s)
      for(dd in seq_along(d.cand)){
        match.ind[cc, ] <- c(r, p, dd)
        cc <- cc + 1
      }
  colnames(match.ind) <- c('r', 's', 'dd')

  est.out <- mclapply(1:nrow(match.ind), function(ii){
    r <- match.ind[ii, 'r']
    p <- match.ind[ii, 's']
    dd <- match.ind[ii, 'dd']

    iip <- which(raw.ind[, 'r'] == r & raw.ind[, 's'] == p)

    X.sonn <- raw.out[[iip]]$U[, 1:d.cand[dd], drop = F] %*%
      raw.out[[iip]]$sigma.half[1:d.cand[dd], 1:d.cand[dd], drop =F]

    return(X.sonn)
  }, mc.cores = min(ncore, nrow(match.ind)))

  match.out <- mclapply(1:nrow(match.ind), function(ii){
    r <- match.ind[ii, 'r']
    p <- match.ind[ii, 's']
    dd <- match.ind[ii, 'dd']

    if(p == 1) return(est.out[[ii]][-(1:o), , drop = F])

    stand <- which(match.ind[, 'r'] == r & match.ind[, 'dd'] == dd &
                     match.ind[, 's'] == 1)

    proc.mat <- netcrop_procrustes(X = est.out[[ii]][(1:o), , drop = F],
                           Xstar = est.out[[stand]][(1:o), , drop = F])$R

    X.rot <- est.out[[ii]][-(1:o), , drop = F] %*% proc.mat

    return(X.rot)
  }, mc.cores = min(ncore, nrow(match.ind)))

  non.size <- s*(s-1)/2

  non.mat <- matrix(nrow = R*non.size*ld, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(dd in seq_along(d.cand))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, dd, p, q)
          cc <- cc + 1
        }
  colnames(non.mat) <- c('r', 'dd', 'p', 'q')

  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 'r']
    dd <- non.mat[ii, 'dd']
    p <- non.mat[ii, 'p']
    q <- non.mat[ii, 'q']

    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]

    A.non <- A[p.non, q.non]

    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(d.cand[dd])

    ind1 <- which(match.ind[, 'r'] == r & match.ind[, 's'] == p &
                    match.ind[, 'dd'] == d.cand[dd])
    ind2 <- which(match.ind[, 'r'] == r & match.ind[, 's'] == q &
                    match.ind[, 'dd'] == d.cand[dd])

    # P.hat <- match.out[[ind1]] %*% t(match.out[[ind2]]) * rho.hat
    P.hat <- tcrossprod(match.out[[ind1]], match.out[[ind2]]) * rho.hat

    P.hat <- pmax(P.hat, 1e-6)

    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        (do.call(loss[lq], list(A.non, P.hat))) / non.size
    }

    return(L.temp)}, mc.cores = min(ncore, nrow(non.mat)))

  L <- list()

  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:ld, function(dd){
      Reduce('+', L.all[which(non.mat[,'dd'] == dd & non.mat[,'r'] == r)])
    }))

    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(d.cand)
  }

  obj <- tibble::tibble(`Candidate Rank` = d.cand)

  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <-
        as(rbind(L[[r]][loss[lq], ]), "vector")
    }

  obj2 <- list()

  obj2[["Candidate Rank"]] <- d.cand

  for(lq in seq_along(loss)){

    obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      l.rdpg <- d.cand[which.min(L[[r]][loss[lq],])]
    })


    obj2[[paste0(loss[lq], ".model")]] <-
      modal(obj2[[paste0("d.hat.each.rep (", loss[lq], ")")]])
  }

  c(list('loss' = obj), obj2)
}

################################################################################
################################################################################
## ECV for RDPG
iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1){
  if(sparse) A <- Matrix(A,sparse=TRUE)
  avg.p <- mean(as.numeric(A),na.rm=TRUE)
  cap <- 1#kappa*avg.p
  A[which(is.na(A))] <- 0
  A <- A/p.sample
  #svd.new <- svd(A,nu=K,nv=K)
  ##print("begin SVD")
  svd.new <- irlba::irlba(A,nu=Kmax,nv=Kmax)
  ##print("end SVD")
  result <- list()
  for(K in 1:Kmax){
    #print(K)
    if(K==1){
      A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
    }else{
      A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
    }
    A.new.thr <- A.new
    A.new.thr <- pmax(A.new.thr, 0 + tau)
    A.new.thr <- pmin(A.new.thr, cap)

    tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
    result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
  }
  return(result)
}

missing.undirected.Rank.weighted.fast.all <-
  function(holdout.index,A,max.K,soft=FALSE,p.sample=1, loss = loss){
    n <- nrow(A)
    #A.new <- A
    #A.new[holdout.index] <- NA
    edge.index <- which(upper.tri(A))
    edge.n <- length(edge.index)
    A.new <- matrix(0,n,n)
    A.new[upper.tri(A.new)] <- A[edge.index]
    A.new[edge.index[holdout.index]] <- NA
    A.new <- A.new + t(A.new)
    diag(A.new) <- diag(A)
    degrees <- colSums(A.new,na.rm=TRUE)
    no.edge <- 0
    no.edge <- sum(degrees==0)

    Omega <- which(is.na(A.new))
    imputed.A <- list()
    sse <- roc.auc <- dev <- rep(0,max.K)
    SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
    for(k in 1:max.K){
      # print(k)
      tmp.est <- SVD.result[[k]]
      #if(k==1){
      #A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
      #}else{
      #   A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
      #}
      A.approx <- tmp.est$A
      response <- A[Omega]
      predictors <- A.approx[Omega]
      #aa <- AUC::roc(predictions=predictors,labels=factor(response))
      if('AUC' %in% loss)
        roc.auc[k] <- AUC(response, predictors)
      if('l2' %in% loss)
        sse[k] <- mean((response-predictors)^2)

      predictors <- pmax(predictors, 1e-6)
      predictors <- pmin(predictors, 1 - 1e-6)
      if('bin.dev' %in% loss)
        dev[k] <- bin.dev(matrix(response, ncol = k), matrix(predictors, ncol = k))

      imputed.A[[k]] <- A.approx
    }
    return(list(imputed.A=imputed.A,Omega=Omega, roc.auc = roc.auc, sse=sse, dev = dev))
  }

ECV.undirected.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,
                                loss = c('l2', 'bin.dev', 'AUC'),
                                ncore = 1, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)

  n <- nrow(A)
  #edge.index <- 1:n^2
  #edge.n <- length(edge.index)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)

  holdout.index.list <- list()

  holdout.n <- floor(holdout.p*edge.n)

  for(j in 1:B){
    holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
  }
  result <- mclapply(holdout.index.list,
                     missing.undirected.Rank.weighted.fast.all,
                     A=A,max.K=max.K,soft=soft,p.sample=1-holdout.p,
                     loss = loss, mc.cores = ncore)

  sse.mat <- roc.auc.mat <- dev.mat <- matrix(0,nrow=B,ncol=max.K)

  for(b in 1:B){
    roc.auc.mat[b,] <- result[[b]]$roc.auc
    sse.mat[b,] <- result[[b]]$sse
    dev.mat[b,] <- result[[b]]$dev
  }

  auc.seq <- colMeans(roc.auc.mat)
  #auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
  sse.seq <- colMeans(sse.mat)
  #sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
  dev.seq <- colMeans(dev.mat)
  # return(list(sse=sse.seq,sse.sd=sse.sd))
  return(list(rank.sse=which.min(sse.seq), sse=sse.seq,
              rank.dev = which.min(dev.seq), dev = dev.seq,
              rank.auc=which.min(auc.seq),auc=auc.seq
  ))
}

ECV.stability.RDPG <- function(A, max.K, cv=3, R = 1, train.p = 0.9,
                               loss = c('l2', 'bin.dev', 'AUC'),
                               ncore = 40, seed = NULL){
  if(ncore >= R*cv){
    outer.ncore <- R
    inner.ncore <- cv
  }else if(ncore >= min(R, cv) & R >= cv){
    outer.ncore <- R
    inner.ncore <- ceiling(ncore/R)
  }else if(ncore >= min(R, cv) & cv > R){
    outer.ncore <- 1
    inner.ncore <- ceiling(ncore/cv)
  }else{
    outer.ncore <- inner.ncore <- 1
  }

  stab.all <- mclapply(1:R, function(rr){
    ECV.undirected.Rank(A = A, max.K = max.K, B = cv,
                        holdout.p = 1 - train.p, seed = seed,
                        loss = loss,
                        ncore = inner.ncore)
  }, mc.cores = outer.ncore)

  best.l2.each.rep <- sapply(stab.all, function(mm) mm$rank.sse)
  best.bin.dev.each.rep <- sapply(stab.all, function(mm) mm$rank.dev)
  best.auc.each.rep <- sapply(stab.all, function(mm) mm$rank.auc)

  best.l2.stable <- modal(best.l2.each.rep)
  best.bin.dev.stable <- modal(best.bin.dev.each.rep)
  best.auc.stable <- modal(best.auc.each.rep)

  list(ecv.loss = stab.all,
       best.l2.each.rep = best.l2.each.rep,
       best.bin.dev.each.rep = best.bin.dev.each.rep,
       best.auc.each.rep = best.auc.each.rep,
       best.l2.stable = best.l2.stable,
       best.bin.dev.stable = best.bin.dev.stable,
       best.auc.stable = best.auc.stable)

}
