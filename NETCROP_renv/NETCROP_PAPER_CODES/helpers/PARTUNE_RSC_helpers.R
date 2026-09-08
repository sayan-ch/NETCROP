library(Matrix)
library(parallel)
library(cluster)

################################################################################
## General helpers: usually used in all functions
modal <- function(x) { # computes modal value in a vector
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

l2 <- function(x, y) {
  sum((x - y)^2)
}

bin.dev <- function(x, y) { # y: prob, x: binary response
  if(y %in% c(0, 1) & y == x)
    return(0)

  if(y %in% c(0, 1) & y != x)
    return(-log(1e-5))

  tmp <- -x * log(y) - (1 - x) * log(1 - y)

  return(sum(tmp, rm.na = T))
}


AUC <- function(AA, PP) {
  -netcrop_auc(group = as(AA, "numeric"), predictions = as(PP, "numeric"))
}

source(file.path(here::here(), "NETCROP_PAPER_CODES/helpers/General_helpers.R"))
################################################################################
##SBM or DCBM generator
SBM.gen <- function(n, K, B, g = NULL, PI = rep(1 / K, K),
                    avg.deg = NULL, ncore = 1, seed = NULL){
  set.seed(seed)

  if(is.null(g))
    g <- sample.int(K, n, T, PI)

  if(is.null(avg.deg)){
    stor <- do.call('rbind', parallel::mclapply(1:(n - 1), function(i) {
      if(!is.null(seed))
        set.seed(seed + i)
      tmp <- which(rbinom(n - i, 1, B[g[i], g[(i + 1):n]]) == 1)

      if (length(tmp) == 0)
        return(NULL)
      else
        return(cbind(rep(i, length(tmp)), i + tmp))
    }, mc.cores = ncore))

    A <- sparseMatrix(i = c(stor[, 1], stor[, 2]),
                      j = c(stor[, 2], stor[, 1]),
                      x = 1, dims = c(n, n))

    return(list(A = A, g = g, psi = psi))
  }

  psi.mat <- sparseMatrix(
    i = 1:n,
    j = g,
    x = 1,
    dims = c(n, K)
  )
  P0 <- psi.mat %*% B %*% t(psi.mat)

  alpha <-  avg.deg / mean(rowSums(P0))

  P <- P0 * alpha

  P <- pmin(P, 1 - 1e-6)
  P <- pmax(P, 1e-6)
  diag(P) <- 0

  stor <- do.call('rbind', mclapply(1:(n - 1), function(i) {
    if(!is.null(seed))
      set.seed(seed + i)
    tmp <- which(rbinom(n - i, 1, P[i, (i + 1):n]) == 1)

    if (length(tmp) == 0)
      return(NULL)
    else
      return(cbind(rep(i, length(tmp)), i + tmp))
  }, mc.cores = ncore))

  A <- sparseMatrix(i = stor[, 1], j = stor[, 2], x = 1,
                    dims = c(n, n))
  A <- A + t(A)

  return(list(A = A, member = g, alpha = alpha))
}

DCBM.gen <- function(n, K, beta = 0, rho = 1, g = NULL, psi = NULL,
                     avg.deg = NULL,
                     PI = rep(1 / K, K), ncore = 1, seed = NULL) {
  set.seed(seed)

  if(is.null(g))
    g <- sample.int(K, n, T, PI)

  if(is.null(psi)){
    psi.pool <- rbeta(300, 1, 4)
    psi <- sample(psi.pool, n, replace = T)

    for(kk in 1:K){
      psi[g == kk] <- psi[g == kk] / max(psi[g == kk])
    }
  }

  B0 <- rho*(diag(1 - beta, nrow = K) + beta)

  psi.mat <- sparseMatrix(
    i = 1:n,
    j = g,
    x = psi,
    dims = c(n, K)
  )
  P0 <- psi.mat %*% B0 %*% t(psi.mat)

  alpha <- 1

  if(!is.null(avg.deg)){
    alpha <-  avg.deg / mean(rowSums(P0))
  }

  P <- alpha * P0

  P <- pmin(P, 1 - 1e-6)
  P <- pmax(P, 1e-6)
  diag(P) <- 0

  stor <- do.call('rbind', mclapply(1:(n - 1), function(i) {
    if(!is.null(seed))
      set.seed(seed + i)

    tmp <- which(rbinom(n - i, 1, P[i, (i + 1):n]) == 1)

    if (length(tmp) == 0)
      return(NULL)
    else
      return(cbind(rep(i, length(tmp)), i + tmp))
  }, mc.cores = ncore))

  A <- sparseMatrix(i = stor[, 1], j = stor[, 2], x = 1,
                    dims = c(n, n))
  A <- A + t(A)

  return(list(A = A, member = g, psi = psi, alpha = alpha))
}

## MatchGreedy Algorithm
best.perm.label.match <- function(lab, fixed,
                                  n = length(lab), K = max(lab, fixed)){

  if(identical(lab, fixed))
    return(diag(1, K))

  if(K == 1)
    return(matrix(1,1,1))

  if(K == 2){
    if(sum(lab!=fixed) <= n/2)
      return(diag(1,2))
    else
      return(matrix(c(0,1,1,0),2,2,T))
  }

  E <- matrix(0, K, K)

  C.lab <- as(sparseMatrix(i = 1:n, j = lab, dims = c(n, K)), 'dMatrix')
  C.fixed <- as(sparseMatrix(i = 1:n, j = fixed, dims = c(n, K)), 'dMatrix')
  M <- crossprod(C.lab, C.fixed)

  while(max(M) != -1)
  {
    ind <- which(M == max(M), T)[1,]
    E[ind[2],ind[1]] <- 1
    M[ind[1],] <- rep(-1,K)
    M[,ind[2]]  <- rep(-1,K)
  }
  return(E)
}

matched.lab <- function(lab, fixed,
                        n = length(lab), K = max(lab, fixed)){

  E <- best.perm.label.match(lab, fixed, n = n, K = K)

  lmat <- sparseMatrix(i = 1:n, j = lab, dims = c(n,K))

  as.vector(tcrossprod(tcrossprod(lmat, E), rbind(1:K)))
}

## MatchGreedy with Hugarian : not used in simulations
vec.perm.label.match <- function(lab, fixed, n = length(lab), K = max(lab, fixed)) {
  if (identical(lab, fixed))
    return(1:K)

  M <- n - table(x = factor(lab, levels = 1:K),
                 y = factor(fixed, levels = 1:K))

  RcppHungarian::HungarianSolver(M)$pairs[,2]
}

## Estimates B.hat from A and g.hat for SBM
fast.SBM.est <- function(A, g, n = nrow(A), K = max(g)) {
  B <- matrix(0, K, K)
  if (K == 1) {
    B[K, K] <- sum(A) / (n^2 - n)
    return(B)
  }

  G <- lapply(1:K, function(k)
    which(g == k))
  nk <- sapply(G, 'length')

  for (k in 1:K) {
    for (l in k:K) {
      B[k, l] <- B[l, k] <- sum(A[G[[k]], G[[l]]]) / (nk[k] * nk[l])
    }
  }

  diag(B) <- diag(B) * nk / pmax((nk - 1), 1)

  return(B)
}

## Estimates B.hat and Psi.hat from A and g.hat for DCBM. Uses plug-in estimator
fast.DCBM.est <- function(A, g, n = nrow(A), K = max(g), psi.omit = 0) {
  B.sum <- matrix(0, K, K)
  if (K == 1) {
    B.sum[K, K] <- sum(A) + 0.01

    if (psi.omit > 0) {
      psi <- as.numeric(rowSums(A[-(1:psi.omit), ]) / (B.sum[K, K])) # 1e-3 to avoid problems with very sparse networks
      return(list(Bsum = B.sum, psi = psi))
    }

    psi <- as.numeric(rowSums(A) / B.sum[K, K])

    return(list(Bsum = B.sum, psi = psi))
  }

  G <- lapply(1:K, function(k)
    which(g == k))

  for (k in 1:K) {
    for (l in k:K) {
      B.sum[k, l] <- B.sum[l, k] <- sum(A[G[[k]], G[[l]]]) + 0.01
    }
  }

  if (psi.omit > 0) {
    psi <- as.numeric(rowSums(A[-(1:psi.omit), ]) /
                        rowSums(B.sum)[g[-(1:psi.omit)]])
    return(list(Bsum = B.sum, psi = psi))
  }

  psi <- as.numeric(rowSums(A) / (rowSums(B.sum)[g]))

  return(list(Bsum = B.sum, psi = psi))
}

## Estimates DCBM parameters by spectral method as in NCV paper
eigen.DCBM.est <- function(A, g, rownorm = NULL,
                           n = nrow(A), K = max(g), psi.omit = 0) {
  if(is.null(rownorm)){
    U.hat <- irlba::irlba(A, nu = K, nv = K)$v

    psi.hat <- rowSums(U.hat^2)^0.5
  }else{
    psi.hat <- rownorm
  }

  psi.outer <- psi.hat %*% t(psi.hat)

  B.sum <- matrix(0, K, K)
  if (K == 1) {
    B.sum[K, K] <- sum(A) / sum(psi.outer)

    if (psi.omit > 0) {
      return(list(Bsum = B.sum, psi = psi.hat[-(1:psi.omit)]))
    }

    return(list(Bsum = B.sum, psi = psi.hat))
  }

  G <- lapply(1:K, function(k)
    which(g == k))

  for (k in 1:K) {
    for (l in k:K) {
      B.sum[k, l] <- B.sum[l, k] <-
        sum(A[G[[k]], G[[l]]]) / sum(psi.outer[G[[k]], G[[l]]])
    }
  }

  if (psi.omit > 0) {
    return(list(Bsum = B.sum, psi = psi.hat[-(1:psi.omit)]))
  }

  return(list(Bsum = B.sum, psi = psi.hat))
}


################################################################################
pair.NMI.loss <- function(g1, g2, g3, g4){

  pair1 <- expand.grid(g1 = g1, g2 = g2)

  pair11 <- pair1 |> mutate(commu = (g1 - 1) * max(g2) + g2,
                            id = 1:nrow(pair1)) |>
    select(id, commu)

  pair2 <- expand.grid(g3 = g3, g4 = g4)
  pair22 <- pair2 |> mutate(commu = (g3 - 1) * max(g4) + g4,
                            id = 1:nrow(pair2)) |>
    select(id, commu)

  return(-netcrop_nmi(pair11$commu, pair22$commu))
}

pair.hemming.loss <- function(g1, g2, g3, g4){

  pair1 <- expand.grid(g1 = g1, g2 = g2)

  pair11 <- pair1 |> mutate(commu = (g1 - 1) * max(g2) + g2,
                            id = 1:nrow(pair1)) |>
    select(id, commu)

  pair2 <- expand.grid(g3 = g3, g4 = g4)
  pair22 <- pair2 |> mutate(commu = (g3 - 1) * max(g4) + g4,
                            id = 1:nrow(pair2)) |>
    select(id, commu)

  mean(matched.lab(pair11$commu, pair22$commu) != pair22$commu)
}

################################################################################

netcrop.tune.regsp <- function(A, K, tau.cand,
                               DCBM = F,
                               s, o, R,
                               laplace = F,
                               dc.est = 2,
                               loss = c("pair.NMI.loss",
                                        "pair.hemming.loss",
                                        "l2", "bin.dev", "AUC"),
                               true.g = NULL,
                               ncore = 1, seed = NULL){
  set.seed(seed)

  n <- nrow(A)
  m <- (n-o)/s

  L <- list()

  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, replace = F))

  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))

  raw.out <- mclapply(1:nrow(raw.ind), function(ii){

    q <- raw.ind[ii, 2]
    r <- raw.ind[ii, 1]

    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]

    deg2 <- rowSums(A.sonn)
    avg.deg <- mean(deg2)

    out.BM <- list()

    for(tt in seq_along(tau.cand)){
      if(K == 1){
        out.BM[[tt]] <- rep(1, o+m)
        next
      }

      if(tau.cand[tt] > 0){
        A.sonn.tau <- A.sonn + tau.cand[tt]*avg.deg/(o+m)

        d.sonn.tau <- sparseMatrix( i = 1:(o+m), j = 1:(o+m),
                                    x = 1/sqrt(deg2 + tau.cand[tt]*avg.deg))

        L.sonn <- A.sonn.tau + 1 - 1

        if(laplace){
          L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
          L.sonn[is.na(L.sonn)] <- 0
        }

        eig.max <- RSpectra::eigs_sym(L.sonn, K, which = "LM",
                                      opts = list(v0 = rep(1, o+m)))$vectors

        rn.eig <- eig.max

        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max^2))
          rownorm[rownorm == 0] <- 1

          rn.eig <- eig.max/rownorm

          out.BM[[tt]] <- as.integer(clara(rn.eig, K,
                                           metric = "manhattan",
                                           cluster.only = T))
        }else{
          out.BM[[tt]] <- as.integer(kmeans(eig.max, K, nstart = 100,
                                            iter.max = 10^7)$cluster)
        }
      }
      if(tau.cand[tt] == 0){
        out.BM[[tt]] <- vector(length = o+m)

        bad.node <- which(deg2 == 0)
        out.BM[[tt]][bad.node] <- sample(1:K, length(bad.node), replace = T)

        good.node <- which(deg2 > 0)

        L.sonn.good <- A.sonn[good.node, good.node] + 1 - 1
        d.sonn.good <- sparseMatrix( i = 1:length(good.node), j = 1:length(good.node),
                                     x = 1/sqrt(deg2[good.node]))

        if(laplace){
          L.sonn.good <- tcrossprod(crossprod(d.sonn.good, L.sonn.good),
                                    d.sonn.good)
          L.sonn.good[is.na(L.sonn.good)] <- 0
        }

        eig.max.good <- RSpectra::eigs_sym(L.sonn.good, K, which = "LM",
                                           opts = list(v0 = rep(1, length(good.node))))$vectors

        rn.eig.good <- eig.max.good

        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max.good^2))
          rownorm[rownorm == 0] <- 1

          rn.eig.good <- eig.max.good/rownorm

          out.BM[[tt]][good.node] <- as.integer(clara(rn.eig.good, K,
                                                      metric = "manhattan",
                                                      cluster.only = T))
        }else{
          out.BM[[tt]][good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                                       iter.max = 10^7)$cluster)
        }

      }
    }

    return(list('BM' = out.BM))
  },
  mc.cores = min(ncore, nrow(raw.ind)))

  tau.size <- length(tau.cand)

  est.out <- mclapply(1:(tau.size*nrow(raw.ind)), function(ii){

    tt <- ii %% tau.size
    tt <- ifelse(tt == 0, tau.size, tt)

    rot <- ceiling(ii / tau.size)

    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]

    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]

    out.BM.std <- raw.out[raw.ind[,1] == r][[1]]$BM[[tt]]

    out.BM <- raw.out[raw.ind[,1] == r][[q]]$BM[[tt]]

    work.tau <- tau.cand[tt]

    if(K == 1){
      mat.BM <- rep(1, m)

      if(!DCBM){
        B.BM <- fast.SBM.est(A = A.sonn, g = rep(1,o+m))
        mat.BM <- rep(1, m)
        psi.BM <- rep(1, m)

        return(list('gBM' = mat.BM, 'BBM' = B.BM,
                    'psiBM' = psi.BM))
      }

      if(DCBM){
        if(dc.est == 2){
          tmp <- fast.DCBM.est(A = A.sonn, g = rep(1,o+m),
                               ps.omit = o)
        }else{
          tmp <- eigen.DCBM.est(A = A.sonn, g = rep(1,o+m),
                                psi.omit = o)
        }
      }

      B.BM <- tmp$Bsum
      psi.BM <- tmp$psi
      mat.BM <- rep(1, m)

      return(list('gBM' = mat.BM,
                  'BBM' = B.BM,
                  'psiBM' = psi.BM))
    }

    E.BM.kc <- best.perm.label.match(out.BM[1:o],
                                     out.BM.std[1:o],
                                     o, K)

    tmp.BM <- sparseMatrix(i = 1:(o+m), j = out.BM, dims = c((o+m), K))

    mat.BM <- as.vector(tcrossprod(tcrossprod(tmp.BM, E.BM.kc),
                                   rbind(1:K)))
    B.BM <- psi.BM <- NULL

    if(any(loss %in% c("l2", "bin.dev", "AUC"))){
      if(!DCBM){
        B.BM <- fast.SBM.est(A = A.sonn, g = mat.BM)
        mat.BM <- mat.BM[-(1:o)]
        psi.BM <- rep(1, m)

        return(list('gBM' = mat.BM, 'BBM' = B.BM,
                    'psiBM' = psi.BM))
      }

      if(dc.est == 2){
        tmp <- fast.DCBM.est(A = A.sonn, g = mat.BM,
                             psi.omit = o)
      }else{
        tmp <- eigen.DCBM.est(A = A.sonn, g = mat.BM, #o+m, K,
                              psi.omit = o)
      }

      B.BM <- tmp$Bsum
      psi.BM <- tmp$psi
    }

    message(paste0("Est. at s=",q, " finished"))

    return(list('gBM' = mat.BM[(o+1):(o+m)], 'ggBM' = mat.BM,
                'BBM' = B.BM, 'psiBM' = psi.BM))
  },
  mc.cores = min(ncore, tau.size*nrow(raw.ind)))

  g.BM <- list()
  B.BM <- list()
  psi.BM <- list()

  gg.BM <- list()

  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = tau.size), ],
                   rep(1:tau.size, nrow(raw.ind)))

  for(r in 1:R){
    g.BM[[r]] <- list()
    B.BM[[r]] <- list()
    psi.BM[[r]] <- list()
    gg.BM[[r]] <- list()

    for(tt in seq_along(tau.cand)){
      tmp.est <- est.out[which(raw.mat[,3] == tt & raw.mat[,1] == r)]
      g.BM[[r]][[tt]] <- list()
      B.BM[[r]][[tt]] <- 0
      psi.BM[[r]][[tt]] <- list()

      gg.BM[[r]][[tt]] <- vector(length = n)

      gg.BM[[r]][[tt]][over[[r]]] <- tmp.est[[1]]$ggBM[1:o]

      for(q in 1:s){
        g.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$gBM
        B.BM[[r]][[tt]] <- B.BM[[r]][[tt]] +
          tmp.est[[q]]$BBM/s
        psi.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$psiBM

        gg.BM[[r]][[tt]][non.over[[r]][((q-1)*m+1):(q*m)]] <-
          tmp.est[[q]]$ggBM[(o+1):(o+m)]
      }
    }

  }

  non.size <- s*(s-1)/2
  non.mat <- matrix(nrow = R*non.size*tau.size, ncol = 4)
  cc <- 1
  for(r in 1:R)
    for(tt in seq_along(tau.cand))
      for(p in 1:(s-1))
        for(q in (p+1):s){
          non.mat[cc, ] <- c(r, tt, p, q)
          cc <- cc + 1
        }

  ##regSC on the entire network
  deg <- rowSums(A)
  avg.deg <- mean(deg)

  all.out <- mclapply(seq_along(tau.cand), function(tt){
    if(tau.cand[tt] > 0){
      A.tau <- A + tau.cand[tt]*avg.deg/n

      d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                             x = 1/sqrt(deg + tau.cand[tt]*avg.deg))

      L.tau <- as(A.tau, 'dMatrix')

      if(laplace){
        L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
        L.tau[is.na(L.tau)] <- 0
      }

      eig.max <- RSpectra::eigs_sym(L.tau, K, which = "LM",
                                    opts = list(v0 = rep(1, o+m)))$vectors

      if(K == 1){
        out.BM[[tt]] <- rep(1, n)
        next
      }

      rn.eig <- eig.max

      if(DCBM){
        rownorm <- sqrt(rowSums(eig.max^2))
        rownorm[rownorm == 0] <- 1

        rn.eig <- eig.max/rownorm

        out.BM <- as.integer(clara(rn.eig, K,
                                   metric = "manhattan",
                                   cluster.only = T))
        return(out.BM)
      }else{
        out.BM <- as.integer(kmeans(eig.max, K, nstart = 100,
                                    iter.max = 10^7)$cluster)
        return(out.BM)
      }

    }

    if(tau.cand[tt] == 0){
      out.BM <- vector(length = n)

      bad.node <- which(deg == 0)
      out.BM[bad.node] <- sample(1:K, length(bad.node), replace = T)

      good.node <- which(deg > 0)

      L.good <- as(A[good.node, good.node], 'dMatrix')
      d.good <- sparseMatrix( i = 1:length(good.node), j = 1:length(good.node),
                              x = 1/sqrt(deg[good.node]))

      if(laplace){
        L.good <- tcrossprod(crossprod(d.good, L.good), d.good)
        L.good[is.na(L.good)] <- 0
      }

      eig.max.good <- RSpectra::eigs_sym(L.good, K, which = "LM",
                                         opts = list(v0 = rep(1, length(good.node))))$vectors

      rn.eig.good <- eig.max.good

      if(DCBM){
        rownorm <- sqrt(rowSums(eig.max.good^2))
        rownorm[rownorm == 0] <- 1

        rn.eig.good <- eig.max.good/rownorm

        out.BM[good.node] <- as.integer(clara(rn.eig.good, K,
                                              metric = "manhattan",
                                              cluster.only = T))
      }else{
        out.BM[good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                               iter.max = 10^7)$cluster)
      }

      return(out.BM)
    }
  },
  mc.cores = min(ncore, length(tau.cand)))



  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 1]
    tt <- non.mat[ii, 2]
    p <- non.mat[ii, 3]
    q <- non.mat[ii, 4]

    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]

    A.non <- A[p.non, q.non]

    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(tau.cand[tt])

    if(any(loss %in% c("l2", "bin.dev", "AUC"))){
      P.BM <- B.BM[[r]][[tt]][g.BM[[r]][[tt]][[p]],
                              g.BM[[r]][[tt]][[q]] ] *
        tcrossprod(psi.BM[[r]][[tt]][[p]],
                   psi.BM[[r]][[tt]][[q]])

      P.BM <- pmax(P.BM, 1e-6)
      P.BM <- pmin(P.BM, 1 - 1e-6)
    }

    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      if(tmp.nm %in% c("l2", "bin.dev", "AUC")){
        L.temp[tmp.nm, 1] <- L.temp[tmp.nm, 1] +
          do.call(loss[lq], list(as.numeric(A.non), P.BM))/(s*(s-1)*0.5)
      }else{
        L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
          do.call(loss[lq], list(g.BM[[r]][[tt]][[p]],
                                 g.BM[[r]][[tt]][[q]],
                                 all.out[[tt]][p.non],
                                 all.out[[tt]][q.non]))
      }
    }

    return(L.temp)},
    mc.cores = min(ncore, nrow(non.mat)))

  L <- list()

  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:tau.size, function(tt){
      Reduce('+', L.all[which(non.mat[,2] == tt & non.mat[,1] == r)])
    }))

    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(tau.cand)
  }

  obj <- tibble::tibble(`Candidate Tau` = tau.cand)

  for(lq in seq_along(loss))
    for(r in 1:R){
      obj[[paste0(loss[lq], "-Rep=", r)]] <-
        as(rbind(L[[r]][loss[lq], ]), "vector")
    }

  obj2 <- list()

  obj2[["Candidate Tau"]] <- tau.cand

  for(lq in seq_along(loss)){

    obj2[[paste0("tau.hat.each.rep (", loss[lq], ")")]] <- sapply(1:R, function(r){
      tau.BM <- tau.cand[which.min(L[[r]][loss[lq],])]
    })


    obj2[[paste0(loss[lq], ".model")]] <-
      mean(obj2[[paste0("tau.hat.each.rep (", loss[lq], ")")]])
  }

  if(is.null(true.g)){
    return(list('obj' = obj, 'obj2' = obj2))
  }else{

    all.accu <- do.call(c, mclapply(all.out, function(out.g){
      return(mean(matched.lab(out.g, true.g) == true.g))
    }, mc.cores = ncore))

    sonnet.accu <- vector()
    for(tt in seq_along(tau.cand)){
      mode.mat <- matrix(nrow = n, ncol = R)
      for(r in 1:R){
        mode.mat[, r] <- gg.BM[[r]][[tt]]
      }
      vote.mat <- apply(mode.mat, 1, modal)
      sonnet.accu[tt] <- mean(matched.lab(vote.mat, true.g) == true.g)
    }

    croissant.all.accu <- vector()
    sonnet.all.accu <- vector()

    for(lq in seq_along(loss)){

      croissant.all.accu[paste0(loss[lq])] <- c(
        all.accu[which(tau.cand ==
                         obj2[[paste0("tau.hat.each.rep (",
                                      loss[lq], ")")]][1])])

      croissant.all.accu[paste0(loss[lq], ".mode")] <- c(
        all.accu[which(tau.cand ==
                         modal(obj2[[paste0("tau.hat.each.rep (",
                                            loss[lq], ")")]]))])

      sonnet.all.accu[paste0(loss[lq])] <- c(
        sonnet.accu[which(tau.cand ==
                            obj2[[paste0("tau.hat.each.rep (",
                                         loss[lq], ")")]][1])])

      sonnet.all.accu[paste0(loss[lq], ".mode")] <- c(
        sonnet.accu[which(tau.cand ==
                            modal(obj2[[paste0("tau.hat.each.rep (",
                                               loss[lq], ")")]]))])

      tau.mean <- mean(obj2[[paste0("tau.hat.each.rep (",
                                    loss[lq], ")")]])

      A.tau <- A + tau.mean*avg.deg/n

      d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                             x = 1/sqrt(deg + tau.mean*avg.deg))

      L.tau <- as(A.tau, 'dMatrix')

      if(laplace){
        L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
        L.tau[is.na(L.tau)] <- 0
      }

      # eig.max <- irlba::partial_eigen(x = L.tau, n = K,
      # symmetric = T)$vectors
      eig.max <- RSpectra::eigs_sym(L.tau, K, which = "LM",
                                    opts = list(v0 = rep(1, o+m)))$vectors

      if(K == 1){
        out.BM <- rep(1, n)
      }else{

        rn.eig <- eig.max

        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max^2))
          rownorm[rownorm == 0] <- 1

          rn.eig <- eig.max/rownorm
        }

        out.BM <- as.integer(clara(rn.eig, K,
                                   metric = "manhattan",
                                   cluster.only = T))
      }

      croissant.all.accu[paste0(loss[lq], ".mean")] <-
        mean(matched.lab(out.BM, true.g) == true.g)
    }

    return(
      list('obj' = obj, 'obj2' = obj2,
           'all.accu' = all.accu,
           'sonnet.accu' = sonnet.accu,
           'croissant.all.accu' = croissant.all.accu,
           'sonnet.all.accu' = sonnet.all.accu)
    )
  }

}

################################################################################
################################################################################
DKest <- function(A, K, true.g = NULL, tau.cand, laplace = T,
                  DCBM = T, DC.est = 2,
                  ncore = 1){
  n <- nrow(A)

  deg <- rowSums(A)
  avg.deg <- mean(deg)

  DK.out <- mclapply(seq_along(tau.cand), function(tt){

    A.tau <- A + tau.cand[tt]*avg.deg/n

    d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                           x = 1/sqrt(deg + tau.cand[tt]*avg.deg))

    L.tau <- A.tau

    if(laplace){
      L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
      L.tau[is.na(L.tau)] <- 0
    }

    if(K == 1){
      out.comm <- rep(1, n)
    }else{
      if(tau.cand[tt] > 0){

        eig.max <- RSpectra::eigs_sym(L.tau, K, which = "LM",
                                      opts = list(v0 = rep(1, n)))$vectors

        rn.eig <- eig.max

        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max^2))
          rownorm[rownorm == 0] <- 1

          rn.eig <- eig.max/rownorm

          out.comm <- as.integer(pam(rn.eig, K,
                                     metric = "euclidean",
                                     do.swap = F, cluster.only = T,
                                     pamonce = 6))
          # out.comm <- as.integer(clara(rn.eig, K,
          #                              metric = "manhattan",
          #                              cluster.only = T))
        }else{
          # out.comm <- as.integer(kmeans(eig.max, K, nstart = 100,
          #                               iter.max = 10^7)$cluster)
          out.comm <- as.integer(clara(eig.max, K,
                                       metric = "euclidean",
                                       cluster.only = T))
        }
      }

      if(tau.cand[tt] == 0){
        out.comm <- vector(length = n)

        bad.node <- which(deg == 0)
        out.comm[bad.node] <- sample(1:K, length(bad.node), replace = T)

        good.node <- which(deg > 0)

        L.good <- L.tau[good.node, good.node]

        eig.max.good <- RSpectra::eigs_sym(L.good, K, which = "LM",
                                           opts = list(v0 = rep(1, length(good.node))))$vectors

        rn.eig.good <- eig.max.good

        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max.good^2))
          rownorm[rownorm == 0] <- 1

          rn.eig.good <- eig.max.good/rownorm

          out.comm[good.node] <- as.integer(pam(rn.eig.good, K,
                                                metric = "euclidean",
                                                do.swap = F, cluster.only = T,
                                                pamonce = 6))
        }else{
          out.comm[good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                                   iter.max = 10^7)$cluster)
        }
      }
    }

    if(!DCBM){
      B.hat <- fast.SBM.est(A, g = out.comm)

      P.hat <- B.hat[out.comm, out.comm]
    }else if(DC.est == 1){
      par.hat <- fast.DCBM.est(A, g = out.comm)

      P.hat <- par.hat$B[out.comm, out.comm] * tcrossprod(par.hat$psi)
    }else{
      par.hat <- eigen.DCBM.est(A, g = out.comm)

      P.hat <- par.hat$B[out.comm, out.comm] * tcrossprod(par.hat$psi)
    }

    P.hat <- pmin(P.hat, 1)
    P.hat <- pmax(P.hat, 0)
    deg.hat <- rowSums(P.hat)
    P.hat.tau <- P.hat + tau.cand[tt]*mean(deg)/n

    d.tau.hat <- sparseMatrix( i = 1:n, j = 1:n,
                               x = 1/sqrt(deg.hat + tau.cand[tt]*mean(deg)/n))

    L.hat.tau <- P.hat.tau + 1 - 1

    if(laplace){
      L.hat.tau <- tcrossprod(crossprod(d.tau.hat, P.hat.tau), d.tau.hat)
      L.hat.tau[is.na(L.hat.tau)] <- 0
    }

    numerator <- RSpectra::eigs_sym((L.tau - L.hat.tau), K, which = "LM",
                                    opts = list(v0 = rep(1, n)))$values[1]

    denominator <- RSpectra::eigs_sym(L.hat.tau, K, which = "LM",
                                      opts = list(v0 = rep(1, n)))$values[K]

    accuracy <- NULL

    if(!is.null(true.g)){
      accuracy <- mean(matched.lab(out.comm, true.g) == true.g)
    }

    return(c(tau.cand = tau.cand[tt], DK.stat = abs(numerator/denominator),
             accu = accuracy))

  }, mc.cores = ncore)

  DK.out <- do.call('rbind', DK.out)

}
