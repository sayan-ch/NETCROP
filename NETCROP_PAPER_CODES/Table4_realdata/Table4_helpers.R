source(file.path(here::here(), "NETCROP_PAPER_CODES/helpers/General_helpers.R"))



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
  y <- pmax(y, 1e-5)
  y <- pmin(y, 1 - 1e-5)
  tmp <- -x * log(y) - (1 - x) * log(1 - y)

  return(sum(tmp, rm.na = T))
}

auroc <- function(score, bool) {
  n1 <- sum(!bool)
  #n2 <- sum(bool)
  n2 <- length(score) - n1
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

AUC <- function(AA, PP) {
  -netcrop_auc(group = as(AA, "numeric"), predictions = as(PP, "numeric"))
  #-auroc(as(P, 'vector'), as(A, 'vector'))
}

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

DCBM.gen <- function(n, K, avg.deg, beta = 0, g = NULL,
                            PI = rep(1 / K, K), ncore = 1, seed = NULL) {
  set.seed(seed)

  if(is.null(g))
    g <- sample.int(K, n, T, PI)

  psi <- 1 / rbeta(n, 4, 1)

  B0 <- diag(1 - beta, nrow = K) + beta

  psi.mat <- sparseMatrix(
    i = 1:n,
    j = g,
    x = psi,
    dims = c(n, K)
  )
  P0 <- psi.mat %*% B0 %*% t(psi.mat)

  alpha <-  avg.deg / mean(rowSums(P0))

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
best.perm.label.match <- function(lab, fixed, n = length(lab),
                                  K = max(lab, fixed)) {
  if (identical(lab, fixed))
    return(1:K)

  if (K == 2) {
    if (sum(lab != fixed) <= n / 2)
      return(1:2)
    else
      return(2:1)
  }

  E <- rep(0, K)

  M <- table(x = factor(lab, levels = 1:K),
             y = factor(fixed, levels = 1:K))

  while (max(M) != -1)
  {
    ind <- which(M == max(M), T)[1, ]
    E[ind[2]] <- ind[1]
    # E[ind[1]] <- ind[2]
    M[ind[1], ] <- rep(-1, K)
    M[, ind[2]]  <- rep(-1, K)
  }
  return(E)
}

matched.lab <- function(lab, fixed, n = length(lab), K = max(lab, fixed)) {
  E <- best.perm.label.match(lab, fixed, n = n, K = K)

  E[lab]
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

NCV.SBM.est <- function(A, g, n = nrow(A), K = max(g), fold.nodes) {
  B <- matrix(0, K, K)
  if (K == 1) {
    B[K, K] <- sum(A) / (nrow(A) * ncol(A) - nrow(A)^2)
    return(B)
  }

  G.row <- lapply(1:K, function(k)
    which(g[-fold.nodes] == k))
  nk.row <- sapply(G.row, 'length')

  G.col <- lapply(1:K, function(k)
    which(g == k))
  nk.col <- sapply(G.col, 'length')

  for (k in 1:K) {
    for (l in k:K) {
      B[k, l] <- B[l, k] <- sum(A[G.row[[k]], G.col[[l]]]) / (nk.row[k] * nk.col[l])
    }
  }

  diag(B) <- diag(B) * nk.col / pmax((nk.col - 1), 1)

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

NCV.DCBM.est <- function(A, g, n = nrow(A), K = max(g),
                         fold.nodes) {
  B.sum <- matrix(0, K, K)
  if (K == 1) {
    B.sum[K, K] <- sum(A) + 0.01

    psi <- as.numeric(colSums(A) / B.sum[K, K])

    return(list(Bsum = B.sum, psi = psi[fold.nodes]))
  }

  G.row <- lapply(1:K, function(k)
    which(g[-fold.nodes] == k))
  nk.row <- sapply(G.row, 'length')

  G.col <- lapply(1:K, function(k)
    which(g == k))
  nk.col <- sapply(G.col, 'length')

  for (k in 1:K) {
    for (l in k:K) {
      B.sum[k, l] <- B.sum[l, k] <- sum(A[G.row[[k]], G.col[[l]]]) + 0.01
    }
  }

  psi <- as.numeric(colSums(A[, fold.nodes]) / (colSums(B.sum)[g[fold.nodes]]))

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

NCV.eigen.DCBM.est <- function(A, g, rownorm = NULL,
                           n = nrow(A), K = max(g),
                           fold.nodes) {
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

    return(list(Bsum = B.sum, psi = psi.hat[fold.nodes]))
  }

  G.row <- lapply(1:K, function(k)
    which(g[-fold.nodes] == k))
  nk.row <- sapply(G.row, 'length')

  G.col <- lapply(1:K, function(k)
    which(g == k))
  nk.col <- sapply(G.col, 'length')

  for (k in 1:K) {
    for (l in k:K) {
      B.sum[k, l] <- B.sum[l, k] <-
        sum(A[G.row[[k]], G.col[[l]]]) / sum(psi.outer[G.row[[k]], G.col[[l]]])
    }
  }

  return(list(Bsum = B.sum, psi = psi.hat[fold.nodes]))
}

################################################################################
# NETCROP
netcrop_param <- function(p.test = 0.02, n = NULL,
                                      o.range = c(0, 0.95)){
    over.upper <- (1 - sqrt(p.test))
    over.lower <- 1 - sqrt(2*p.test)

    o.range[o.range == 1] <- 0.999

    over.prop <- over.lower + (over.upper - over.lower) * o.range

    s.range <- ceiling((1 - over.prop)^2 / ((1 - over.prop)^2 - p.test))

    if(is.null(n))
      return(data.frame(p.test = p.test, p_o = over.prop, s = s.range))

    over.range <- ceiling(n * over.prop)
    m.range <- floor((n - over.range) / s.range)

    over.extra <- n - over.range - s.range * m.range
    over.range <- over.range + over.extra

    return(data.frame(p.test = p.test, p_o = over.prop, n = n,
                      s = s.range, o = over.range))
}


netcrop_blockmodel <- function(A, K.CAND, s, o, R = 1, tau = 0, laplace = F, dc.est = 2,
                       loss = c("l2", "bin.dev", "AUC"),
                       mod.cand = c('SBM', 'DCBM'), ncore = 1, seed = NULL) {
  set.seed(seed)
  # if (length(K.CAND) == 1)
  #   K.CAND <- 1:K.CAND

  K.max <- max(K.CAND)

  n <- nrow(A)
  m <- (n - o) / s

  L <- list()

  mod <- mod.cand

  over <- lapply(1:R, function(ii){sample.int(n, o, F)})
  non.over <- lapply(1:R, function(ii)
    sample((1:n)[-over[[ii]]], n - o, replace = F))

  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))

  mc.cores <- ifelse(ncore > 1, min(s*R, ncore), 1)

  raw.out <- mclapply(1:nrow(raw.ind), function(ii) {
    q <- raw.ind[ii, 2]
    r <- raw.ind[ii, 1]

    sonn <- c(over[[r]], non.over[[r]][((q - 1) * m + 1):(q * m)])
    A.sonn <- A[sonn, sonn]

    deg <- rowSums(A.sonn)
    avg.deg <- mean(deg)

    L.sonn <- A.sonn + tau * avg.deg / (o + m)
    # d.sonn.tau <- sparseMatrix(
    #   i = 1:(o + m),
    #   j = 1:(o + m),
    #   x = 1 / sqrt(deg + tau * avg.deg)
    # )
    d.sonn.tau <- Diagonal(n = o+m, x = 1 / sqrt(deg + tau * avg.deg))

    # L.sonn <- A.sonn.tau
    if (laplace) {
      # L.sonn <- tcrossprod(crossprod(d.sonn.tau, L.sonn), d.sonn.tau)
      L.sonn <- d.sonn.tau %*% L.sonn %*% d.sonn.tau
      L.sonn[is.na(L.sonn)] <- 0
    }

    # eig.max <- irlba::irlba(L.sonn, nu = K.max, nv = K.max)$v
    eig.max <- RSpectra::eigs_sym(L.sonn, k = K.max)$vectors

    out.SBM <- list()
    out.DCBM <- list()
    rownorm.list <- list()

    for (k.cand in seq_along(K.CAND)) {
      if (K.CAND[k.cand] == 1) {
        out.SBM[[k.cand]] <- out.DCBM[[k.cand]] <- rep(1, o + m)
        next
      }

      work.K <- K.CAND[[k.cand]]

      out.SBM[[k.cand]] <- as.integer(cluster::clara(
        x = eig.max[, 1:work.K],
        k = work.K,
        # samples = 5, sampsize = min(max(n/10, 40 + 2*work.K), n),
        metric = "euclidean",
        cluster.only = T
      ))

      rownorm <- sqrt(rowSums(eig.max[, 1:work.K]^2))
      rownorm[rownorm == 0] <- 1e-6

      rn.eig <- eig.max[, 1:work.K] / rownorm

      out.DCBM[[k.cand]] <- as.integer(cluster::clara(
        x = rn.eig,
        k = work.K,
        # samples = 5, sampsize = min(max(n/10, 40 + 2*work.K), n),
        metric = "manhattan",
        cluster.only = T
      ))

      rownorm.list[[k.cand]] <- rownorm
    }

    return(list('SBM' = out.SBM, 'DCBM' = out.DCBM, 'rownorm' = rownorm.list))
  }, mc.cores = mc.cores)

  K.size <- length(K.CAND)

  mc.cores <- ifelse(ncore > 1, min(K.size * s * R, ncore), 1)

  est.out <- mclapply(1:(K.size * nrow(raw.ind)), function(ii) {
    k.cand <- ii %% K.size
    k.cand <- ifelse(k.cand == 0, K.size, k.cand)

    rot <- ceiling(ii / K.size)

    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]

    sonn <- c(over[[r]], non.over[[r]][((q - 1) * m + 1):(q * m)])
    A.sonn <- A[sonn, sonn]

    out.SBM.std <- raw.out[raw.ind[, 1] == r][[1]]$SBM[[k.cand]]
    out.DCBM.std <- raw.out[raw.ind[, 1] == r][[1]]$DCBM[[k.cand]]

    out.SBM <- raw.out[raw.ind[, 1] == r][[q]]$SBM[[k.cand]]
    out.DCBM <- raw.out[raw.ind[, 1] == r][[q]]$DCBM[[k.cand]]


    work.K <- K.CAND[[k.cand]]

    if (work.K == 1) {
      mat.SBM <- mat.DCBM <- rep(1, m)

      B.SBM <- fast.SBM.est(A.sonn, rep(1, o + m), o + m, 1)
      mat.SBM <- rep(1, m)

      if (dc.est > 1) {
        tmp <- fast.DCBM.est(A.sonn, rep(1, o + m), o + m, 1, o)
      } else{
        tmp <- eigen.DCBM.est(A = A.sonn, g = rep(1, o + m),
                              rownorm = raw.out[raw.ind[, 1] == r][[q]]$rownorm[[1]],
                              n = o + m, K = 1, psi.omit = o)
      }

      B.DCBM <- tmp$Bsum
      psi.DCBM <- tmp$psi
      mat.DCBM <- rep(1, m)

      return(
        list(
          'gSBM' = mat.SBM,
          'BSBM' = B.SBM,
          'gDCBM' = mat.DCBM,
          'BDCBM' = B.DCBM,
          'psiDCBM' = psi.DCBM
        )
      )
    }

    E.SBM.kc <- best.perm.label.match(lab = out.SBM[1:o], fixed = out.SBM.std[1:o],
                                      n = o, K = work.K)

    E.DCBM.kc <- best.perm.label.match(lab = out.DCBM[1:o], fixed = out.DCBM.std[1:o],
                                       n = o, K = work.K)

    mat.SBM <- E.SBM.kc[out.SBM]
    mat.DCBM <- E.DCBM.kc[out.DCBM]

    B.SBM <- fast.SBM.est(A = A.sonn, g = mat.SBM, n = o + m, K = work.K)
    mat.SBM <- mat.SBM[-(1:o)]

    if (dc.est > 1) {
      tmp <- fast.DCBM.est(A.sonn, mat.DCBM, o + m, work.K, o)
    } else{
      tmp <- eigen.DCBM.est(A = A.sonn, g = mat.DCBM,
                            rownorm = raw.out[raw.ind[, 1] == r][[q]]$rownorm[[k.cand]],
                            n = o + m, K = work.K, psi.omit = o)
    }

    B.DCBM <- tmp$Bsum
    psi.DCBM <- tmp$psi
    mat.DCBM <- mat.DCBM[-(1:o)]

    return(
      list(
        'gSBM' = mat.SBM,
        'BSBM' = B.SBM,
        'gDCBM' = mat.DCBM,
        'BDCBM' = B.DCBM,
        'psiDCBM' = psi.DCBM
      )
    )
  }, mc.cores = ncore)

  g.SBM <- list()
  B.SBM <- list()
  g.DCBM <- list()
  B.DCBM <- list()
  psi.DCBM <- list()

  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = K.size), ], rep(1:K.size, nrow(raw.ind)))

  for (r in 1:R) {
    g.SBM[[r]] <- list()
    B.SBM[[r]] <- list()
    g.DCBM[[r]] <- list()
    B.DCBM[[r]] <- list()
    psi.DCBM[[r]] <- list()
    for (k.cand in seq_along(K.CAND)) {
      tmp.est <- est.out[which(raw.mat[, 3] == k.cand & raw.mat[, 1] == r)]
      B.SBM[[r]][[k.cand]] <- 0
      B.DCBM[[r]][[k.cand]] <- 0

      g.SBM[[r]][[k.cand]] <- list()
      g.DCBM[[r]][[k.cand]] <- list()
      psi.DCBM[[r]][[k.cand]] <- list()

      for (q in 1:s) {
        B.SBM[[r]][[k.cand]] <- B.SBM[[r]][[k.cand]] +
          tmp.est[[q]]$BSBM / s
        B.DCBM[[r]][[k.cand]] <- B.DCBM[[r]][[k.cand]] +
          tmp.est[[q]]$BDCBM / s

        g.SBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$gSBM
        g.DCBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$gDCBM
        psi.DCBM[[r]][[k.cand]][[q]] <- tmp.est[[q]]$psiDCBM
      }
    }
  }

  non.size <- s * (s - 1) / 2
  non.mat <- matrix(nrow = R * non.size * K.size, ncol = 4)
  cc <- 1
  for (r in 1:R)
    for (k.cand in seq_along(K.CAND))
      for (p in 1:(s - 1))
        for (q in (p + 1):s) {
          non.mat[cc, ] <- c(r, k.cand, p, q)
          cc <- cc + 1
        }

  L.all <- mclapply(1:nrow(non.mat), function(ii) {
    r <- non.mat[ii, 1]
    k.cand <- non.mat[ii, 2]
    p <- non.mat[ii, 3]
    q <- non.mat[ii, 4]

    p.non <- non.over[[r]][((p - 1) * m + 1):(p * m)]
    q.non <- non.over[[r]][((q - 1) * m + 1):(q * m)]

    A.non <- A[p.non, q.non]

    L.temp <- matrix(0, nrow = 2 * length(loss), ncol = 1)
    row.names(L.temp) <- paste(rep(mod, each = length(loss)), rep(loss, 2), sep = "_")
    colnames(L.temp) <- as.character(K.CAND[k.cand])


    P.SBM <- B.SBM[[r]][[k.cand]][g.SBM[[r]][[k.cand]][[p]], g.SBM[[r]][[k.cand]][[q]]]

    P.DCBM <- B.DCBM[[r]][[k.cand]][g.DCBM[[r]][[k.cand]][[p]], g.DCBM[[r]][[k.cand]][[q]]] *
      tcrossprod(psi.DCBM[[r]][[k.cand]][[p]], psi.DCBM[[r]][[k.cand]][[q]])

    P.DCBM <- pmin(P.DCBM, 1 - 1e-6)
    P.DCBM <- pmax(P.DCBM, 1e-6)
    for (mq in seq_along(mod)) {
      if (mod[mq] == "SBM") {
        for (lq in seq_along(loss)) {
          tmp.nm <- paste(mod[mq], loss[lq], sep = "_")
            L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
            do.call(loss[lq], list(as.numeric(A.non), P.SBM)) / non.size
        }
        next
      }
      else{
        for (lq in seq_along(loss)) {
          tmp.nm <- paste(mod[mq], loss[lq], sep = "_")
          L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
            do.call(loss[lq], list(as.numeric(A.non), P.DCBM)) / non.size
        }
      }
    }

    return(L.temp)
  }, mc.cores = ncore)

  for (r in 1:R) {
    L[[r]] <- do.call('cbind', lapply(1:K.size, function(kk) {
      Reduce('+', L.all[which(non.mat[, 2] == kk &
                                non.mat[, 1] == r)])
    }))

    row.names(L[[r]]) <- paste(rep(mod, each = length(loss)), rep(loss, 2), sep = "_")
    colnames(L[[r]]) <- as.character(K.CAND)
  }

  obj <- tibble::tibble(
    `Repititon` = rep(1:R, each = length(mod)*length(K.CAND)),
    `Candidate_Model` = rep(rep(mod, each = length(K.CAND)), R),
    `Candidate_Value` = rep(K.CAND, length(mod)*R)
  )

  for (lq in seq_along(loss)) {
    column <- paste0(loss[lq])
    obj[[column]] <- NA_real_
    for (r in 1:R) {
      rows <- obj$Repititon == r
      obj[[column]][rows] <-
        c(L[[r]][paste0("SBM_", loss[lq]), , drop = T],
          L[[r]][paste0("DCBM_", loss[lq]), , drop = T])
    }
  }

  obj2 <- list()

  obj2[["Candidate Models"]] <- mod

  obj2[["Candidate Values"]] <- K.CAND

  for (lq in seq_along(loss)) {
    obj2[[paste0("Mod.K.hat.each.rep (", loss[lq], ")")]] <-
      sapply(1:R, function(r) {
        l.sbm <- min(L[[r]][paste0("SBM_", loss[lq]), ])
        l.dcbm <- min(L[[r]][paste0("DCBM_", loss[lq]), ])
        ifelse(l.dcbm < l.sbm,
               paste0("DCBM-", K.CAND[which.min(L[[r]][paste0("DCBM_", loss[lq]), ])]),
               paste0("SBM-", K.CAND[which.min(L[[r]][paste0("SBM_", loss[lq]), ])]))
      })


    obj2[[paste0(loss[lq], ".model")]] <-
      modal(obj2[[paste0("Mod.K.hat.each.rep (", loss[lq], ")")]])
  }

  return(c(list('loss' = obj), obj2))
}

################################################################################
################################################################################
################################################################################
# Code for ECV taken as is from the randnet package. Only changes made are:
# passing loss functions as user choice, creating stabilised version with parallelization and
# formatting the output in manner consistent with NETCROP
holdout.evaluation.fast.all <- function(holdout.index, A, max.K, tau = 0,
                                        dc.est = 2, p.sample = 1,
                                        loss = c('l2', 'bin.dev', 'AUC')) {
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0, n, n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  degrees <- colSums(A.new, na.rm = TRUE)
  no.edge <- 0
  no.edge <- sum(degrees == 0)

  Omega <- which(is.na(A.new))
  non.miss <- which(!is.na(A.new))

  SVD.result <- iter.SVD.core.fast.all(A.new, max.K, p.sample = p.sample)

  dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <-
    block.sq.err <- impute.sq.err <- loglike <- rep(0, max.K)
  sbm.auc <- dc.auc <- rep(0, max.K)

  for (k in 1:max.K) {
    tmp.est <- SVD.result[[k]]
    A.approx <- tmp.est$A.thr

    response <- A[edge.index[holdout.index]]#A[Omega]
    predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]

    trunc.predictors <- predictors
    trunc.predictors <- pmin(trunc.predictors, 1 - 1e-6)
    trunc.predictors <- pmax(trunc.predictors, 1e-6)

    if (k == 1) {
      pb <- (sum(A.new, na.rm = TRUE) + 1) / (sum(!is.na(A.new)) - sum(!is.na(diag(A.new))) +
                                                1)
      if (pb < 1e-6)
        pb <- 1e-6
      if (pb > 1 - 1e-6)
        pb <- 1 - 1e-6
      A.Omega <- A[Omega]
      block.sq.err[k] <- sum((pb - A[Omega])^2)
      if('bin.dev' %in% loss)
        loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))
    }

    if (k == 1) {
      U.approx <- matrix(tmp.est$SVD$v, ncol = k)
    } else{
      U.approx <- tmp.est$SVD$v[, 1:k, drop = F]
      if (tau > 0) {
        A.approx <- A.approx + tau * mean(colSums(A.approx)) / n
        d.approx <- colSums(A.approx)
        L.approx <- diag(1 / sqrt(d.approx)) %*% A.approx %*% diag(1 / sqrt(d.approx))
        A.approx.svd <- irlba::irlba(L.approx, nu = k, nv = k)
        U.approx <- A.approx.svd$v[, 1:k]
      }
    }

    km <- kmeans(
      U.approx,
      centers = k,
      nstart = 30,
      iter.max = 30
    )
    B <- matrix(0, k, k)
    Theta <- matrix(0, n, k)
    for (i in 1:k) {
      for (j in i:k) {
        N.i <- which(km$cluster == i)
        N.j <- which(km$cluster == j)
        if (i != j) {
          B[i, j] <- B[j, i] <- (sum(A.new[N.i, N.j], na.rm = TRUE) + 1) / (sum(!is.na(A.new[N.i, N.j])) +
                                                                              1)
        } else{
          B[i, j] <- B[j, i] <- (sum(A.new[N.i, N.j], na.rm = TRUE) + 1) /
            (sum(!is.na(A.new[N.i, N.j])) - sum(!is.na(diag(A.new[N.i, N.j]))) + 1)
        }

      }
      Theta[N.i, i] <- 1
    }
    P.hat <- Theta %*% B %*% t(Theta)
    diag(P.hat) <- 0
    block.sq.err[k] <- sum((P.hat[Omega] - A[Omega])^2)
    P.hat.Omega <- P.hat[Omega]
    A.Omega <- A[Omega]
    P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
    P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
    if('bin.dev' %in% loss)
      loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    if('AUC' %in% loss)
      sbm.auc[k] <- AUC(A.Omega, P.hat.Omega) ##SC addition

    #### Degree correct model
    V <- U.approx

    ptm <- proc.time()

    if (k == 1) {
      V.norms <- as.numeric(abs(V))
    } else{
      V.norms <- apply(V, 1, function(x)
        sqrt(sum(x^2)))
    }

    iso.index <- which(V.norms == 0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    inv.V.norms <- 1 / V.norms
    inv.V.norms[iso.index] <- 1

    V.normalized <- diag(as.numeric(inv.V.norms)) %*% V

    if (k == 1) {
      if (dc.est > 1) {
        B <- sum(A.new, na.rm = TRUE) + 0.01

        partial.d <- colSums(A.new, na.rm = TRUE)
        partial.gd <- B
        phi <- rep(0, n)
        B.g <- partial.gd
        phi <- as.numeric(partial.d / B.g)
        B <- B / p.sample
        P.hat <- t(t(matrix(B, n, n) * phi) * phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((pb - A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
      P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
      if('bin.dev' %in% loss)
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      if('AUC' %in% loss)
        dc.auc[k] <- AUC(A.Omega, P.hat.Omega)
    } else{
      km <- kmeans(
        V.normalized,
        centers = k,
        nstart = 30,
        iter.max = 30
      )
      if (dc.est > 1) {
        B <- matrix(0, k, k)
        Theta <- matrix(0, n, k)
        for (i in 1:k) {
          for (j in 1:k) {
            N.i <- which(km$cluster == i)
            N.j <- which(km$cluster == j)
            B[i, j] <- sum(A.new[N.i, N.j], na.rm = TRUE) + 0.01
          }
          Theta[N.i, i] <- 1
        }
        Theta <- Matrix(Theta, sparse = TRUE)
        partial.d <- colSums(A.new, na.rm = TRUE)
        partial.gd <- colSums(B)
        phi <- rep(0, n)
        B.g <- Theta %*% partial.gd
        phi <- as.numeric(partial.d / B.g)
        B <- B / p.sample
        tmp.int.mat <- Theta * phi
        P.hat <- as.matrix(tmp.int.mat %*% B %*% t(tmp.int.mat))
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((P.hat[Omega] - A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
      P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
      if('bin.dev' %in% loss)
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      if('AUC' %in% loss)
        dc.auc[k] <- AUC(A.Omega, P.hat.Omega)
    }
  }
  return(
    list(
      impute.sq.err = impute.sq.err,
      block.sq.err = block.sq.err,
      loglike = loglike,
      roc.auc = roc.auc,
      no.edge = no.edge,
      dc.block.sq.err = dc.block.sq.err,
      dc.loglike = dc.loglike,
      bin.dev = bin.dev,
      sbm.auc = sbm.auc,
      dc.auc = dc.auc
    )
  )
}

iter.SVD.core.fast.all <- function(A, Kmax, tol = 1e-5, max.iter = 100,
                                   sparse = T, tau = 0, p.sample = 1) {
  if(sparse)
    A <- Matrix(A, sparse = TRUE)
  avg.p <- mean(as.numeric(A), na.rm = TRUE)
  cap <- 1#kappa*avg.p
  A[which(is.na(A))] <- 0
  A <- A / p.sample
  svd.new <- irlba::irlba(A, nu = Kmax, nv = Kmax)
  result <- list()
  for (K in 1:Kmax) {
    if (K == 1) {
      A.new <- svd.new$d[1] * matrix(svd.new$u[, 1], ncol = 1) %*% t(matrix(svd.new$v[, 1], ncol =
                                                                              1))
    } else{
      A.new <- A.new + svd.new$d[K] * matrix(svd.new$u[, K], ncol = 1) %*% t(matrix(svd.new$v[, K], ncol =
                                                                                      1))
    }
    A.new.thr <- A.new
    A.new.thr <- pmax(A.new.thr, 0 + tau)
    A.new.thr <- pmin(A.new.thr, cap)

    tmp.SVD <- list(u = svd.new$u[, 1:K],
                    v = svd.new$v[, 1:K],
                    d = svd.new$d[1:K])
    result[[K]] <- list(
      SVD = tmp.SVD,
      A = A.new,
      A.thr = A.new.thr
    )
  }

  return(result)
}

ECV.BM <- function (A, max.K, cv = 3, holdout.p = 0.1, tau = 0, dc.est = 2,
                    loss = c('l2', 'bin.dev', 'AUC'), ncore = 1,
                    seed = 100){
  set.seed(seed)
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  holdout.index.list <- list()
  holdout.n <- floor(holdout.p * edge.n)
  for (j in 1:cv) {
    holdout.index.list[[j]] <- sample(x = edge.n, size = holdout.n)
  }

  result <- mclapply(holdout.index.list,
                     holdout.evaluation.fast.all,
                     A = A,
                     max.K = max.K,
                     tau = tau,
                     dc.est = dc.est,
                     p.sample = 1 - holdout.p,
                     loss = loss,
                     mc.cores = min(cv, ncore)
  )

  dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <-
    impute.err.mat <- block.err.mat <- loglike.mat <-
    matrix(0, nrow = cv, ncol = max.K)
  sbm.auc.mat <- dc.auc.mat <- matrix(0, nrow = cv, ncol = max.K)

  no.edge.seq <- rep(0, cv)
  Omega.list <- A.list <- Imputed.A.list <- list()
  for (b in 1:cv) {
    impute.err.mat[b, ] <- result[[b]]$impute.sq.err
    block.err.mat[b, ] <- result[[b]]$block.sq.err
    loglike.mat[b, ] <- result[[b]]$loglike
    roc.auc.mat[b, ] <- result[[b]]$roc.auc
    bin.dev.mat[b, ] <- result[[b]]$bin.dev
    no.edge.seq[b] <- result[[b]]$no.edge
    dc.block.err.mat[b, ] <- result[[b]]$dc.block.sq.err
    dc.loglike.mat[b, ] <- result[[b]]$dc.loglike
    sbm.auc.mat[b, ] <- result[[b]]$sbm.auc
    dc.auc.mat[b, ] <- result[[b]]$dc.auc
  }

  output <- list(
    sbm.l2.mat = block.err.mat,
    sbm.bin.dev.mat = loglike.mat,
    sbm.auc.mat = sbm.auc.mat,
    dc.l2.mat = dc.block.err.mat,
    dcbm.bindev.mat = dc.loglike.mat,
    dc.auc.mat = dc.auc.mat,
    sbm.l2 = colMeans(block.err.mat),
    sbm.bin.dev = colSums(loglike.mat),
    sbm.auc = colMeans(sbm.auc.mat),
    dcbm.l2 = colMeans(dc.block.err.mat),
    dcbm.bin.dev = colSums(dc.loglike.mat),
    dcbm.auc = colMeans(dc.auc.mat)
  )

  if (min(output$sbm.bin.dev) > min(output$dcbm.bin.dev)) {
    bin.dev.model <- paste("DCBM", which.min(output$dcbm.bin.dev), sep = "-")
  }else {
    bin.dev.model <- paste("SBM", which.min(output$sbm.bin.dev), sep = "-")
  }

  if (min(output$sbm.l2) > min(output$dcbm.l2)) {
    l2.model <- paste("DCBM", which.min(output$dcbm.l2), sep = "-")
  }else {
    l2.model <- paste("SBM", which.min(output$sbm.l2), sep = "-")
  }

  if (min(output$sbm.auc) > min(output$dcbm.auc)) {
    auc.model <- paste("DCBM", which.min(output$dcbm.auc), sep = "-")
  }else {
    auc.model <- paste("SBM", which.min(output$sbm.auc), sep = "-")
  }

  output$l2.model <- l2.model
  output$bin.dev.model <- bin.dev.model
  output$auc.model <- auc.model
  return(output)
}

ECV.stability.BM <- function(A, max.K, train.p = 0.9, cv = 3, R = 20,
                             dc.est = 2, tau = 0,
                             loss = c("l2", "bin.dev", "AUC"),
                             ncore = 1, seed = 100){

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
    ECV.BM(A = A, max.K = max.K, cv = cv, holdout.p = 1 - train.p,
           tau = tau, dc.est = dc.est,
           loss = loss, ncore = inner.ncore, seed = seed + 100*rr)
  }, mc.cores = outer.ncore)

  valid <- sapply(stab.all, function(xx) is.list(xx))
  stab.all <- stab.all[valid]

  best.l2.each.rep <- sapply(stab.all, function(mm) mm$l2.model)
  best.bin.dev.each.rep <- sapply(stab.all, function(mm) mm$bin.dev.model)
  best.auc.each.rep <- sapply(stab.all, function(mm) mm$auc.model)

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

################################################################################
## NCV Chen and Lei 2018
NCV.BM <- function(A, max.K, cv = 3, dc.est = 2,
                   tau = 0, laplace = F,
                   loss = c("l2", "bin.dev", "AUC"),
                   ncore = 1, seed = 100){

  if(cv <= 1) stop("CV must be >= 2")
  if(max.K < 2) stop("max.K must be >= 2")

  set.seed(seed)

  mc.cores = ifelse(ncore > 1, min(cv, ncore), 1)

  kcand <- 1:max.K

  n <- nrow(A)
  n.perm <- sample(1:n, size = n, replace = F)

  fold.all <- mclapply(1:cv, function(xx){
    if(!is.null(seed))
      set.seed(seed + 100*cv + 200*xx)

    inc <- floor(n/cv)
    start <- (xx-1) * inc + 1
    end <- ifelse(xx == cv, n, xx * inc)
    fold.nodes <- n.perm[start:end]
    train.nodes <- setdiff(1:n, fold.nodes)
    train.mat <- A[train.nodes, ]
    test.mat <- A[fold.nodes, fold.nodes]

    if(tau > 0){
      train.mat <- train.mat + tau * mean(rowSums(train.mat)) / n # as row deg is 0 to n
    }
    if(laplace){
      row.deg.neg12 <- Diagonal(n = nrow(train.mat),
                                x = 1 / sqrt(rowSums(train.mat)))
      col.deg.neg12 <- Diagonal(n = n,
                                x = 1 / sqrt(colSums(train.mat)))

      train.mat <- tcrossprod(crossprod(row.deg.neg12, train.mat),
                              col.deg.neg12)
    }

    # train.eig <- irlba::irlba(train.mat, nv = max.K)$v
    train.eig <- tryCatch({irlba::irlba(train.mat, nu = max.K, nv = max.K)$v},
                          error = function(e){
                            irlba::irlba(train.mat + 1e-12, nu = max.K, nv = max.K)$v
                          })

    sbm.l2 <- dcbm.l2 <- sbm.bd <- dcbm.bd <-
      sbm.auc <- dcbm.auc <- rep(0, max.K)

    for(kk in 1:max.K){
      if(kk == 1){
        g.sbm <- rep(1, n)
        g.dcbm <- rep(1, n)

      }else{
        VV <- train.eig[, 1:kk]
        VV.rownorm <- sqrt(rowSums(VV^2))
        VV.rownorm[VV.rownorm == 0] <- 1e-6
        VV.norm <- VV / VV.rownorm

        g.sbm <- kmeans(VV, centers = kk,
                        nstart=30, iter.max=30)$cluster

        g.dcbm <- kmeans(VV.norm, centers = kk,
                         nstart=30, iter.max=30)$cluster
      }

      sbm.test.par <- NCV.SBM.est(A = train.mat, g = g.sbm, K = kk,
                                   fold.nodes = fold.nodes)

      if(dc.est > 1){
        dcbm.test.par <- NCV.DCBM.est(A = train.mat, g = g.dcbm,
                                       K = kk, fold.nodes = fold.nodes)
      }else{
        if(kk > 1)
          dcbm.test.par <- NCV.eigen.DCBM.est(A = train.mat, g = g.dcbm,
                                        rownorm = VV.rownorm, fold.nodes = fold.nodes)
        if(kk == 1)
          dcbm.test.par <- NCV.eigen.DCBM.est(A = train.mat, g = g.dcbm,
                                          fold.nodes = fold.nodes)
      }

      sbm.P.hat <- sbm.test.par[g.sbm[fold.nodes], g.sbm[fold.nodes]]

      dcbm.P.hat <- dcbm.test.par$Bsum[g.dcbm[fold.nodes], g.dcbm[fold.nodes]] *
        dcbm.test.par$psi %*% t(dcbm.test.par$psi)

      dcbm.P.hat <- pmax(dcbm.P.hat, 1e-6)
      dcbm.P.hat <- pmin(dcbm.P.hat, 1 - 1e-6)
      if("l2" %in% loss){
        sbm.l2[kk] <- l2(test.mat, sbm.P.hat)
        dcbm.l2[kk] <- l2(test.mat, dcbm.P.hat)
      }

      if("bin.dev" %in% loss){
        sbm.bd[kk] <- bin.dev(test.mat, sbm.P.hat)
        dcbm.bd[kk] <- bin.dev(test.mat, dcbm.P.hat)
      }
      if("AUC" %in% loss){
        sbm.auc[kk] <- AUC(test.mat, sbm.P.hat)
        dcbm.auc[kk] <- AUC(test.mat, dcbm.P.hat)
      }
    }

    rbind(sbm.l2, sbm.bd, sbm.auc, dcbm.l2, dcbm.bd, dcbm.auc)

  }, mc.cores = mc.cores)

  fold.mean <- Reduce('+', fold.all)/cv

  sbm.winner <- apply(fold.mean[1:3, ], 1, which.min)
  dcbm.winner <- apply(fold.mean[4:6, ], 1, which.min)

  if(fold.mean['sbm.l2', sbm.winner['sbm.l2']] <
     fold.mean['dcbm.l2', dcbm.winner['dcbm.l2']]){
    best.l2 <- paste0("SBM-", sbm.winner['sbm.l2'])
  }else{
    best.l2 <- paste0("DCBM-", dcbm.winner['dcbm.l2'])
  }

  if(fold.mean['sbm.bd', sbm.winner['sbm.bd']] <
     fold.mean['dcbm.bd', dcbm.winner['dcbm.bd']]){
    best.bd <- paste0("SBM-", sbm.winner['sbm.bd'])
  }else{
    best.bd <- paste0("DCBM-", dcbm.winner['dcbm.bd'])
  }

  if(fold.mean['sbm.auc', sbm.winner['sbm.auc']] <
     fold.mean['dcbm.auc', dcbm.winner['dcbm.auc']]){
    best.auc <- paste0("SBM-", sbm.winner['sbm.auc'])
  }else{
    best.auc <- paste0("DCBM-", dcbm.winner['dcbm.auc'])
  }

  list(
    ncv.loss = t(fold.mean),
    best.l2 = best.l2, best.bin.dev = best.bd, best.auc = best.auc
  )
}

NCV.stability.BM <- function(A, max.K, cv = 3, R = 20,
                             dc.est = 2, tau = 0, laplace = F,
                             loss = c("l2", "bin.dev", "AUC"),
                             ncore = 1, seed = 100){

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
    NCV.BM(A = A, max.K = max.K, cv = cv,
           dc.est = dc.est, tau = tau, laplace = laplace,
           loss = loss, ncore = inner.ncore, seed = 100*rr)
  }, mc.cores = outer.ncore)

  ncv.loss <- do.call('rbind', lapply(1:R, function(rr){
    as.data.frame(stab.all[[rr]]$ncv.loss) |>
      dplyr::mutate(rep = rr, .before = 1)
  }))

  best.l2.each.rep <- sapply(stab.all, function(mm) mm$best.l2)
  best.bin.dev.each.rep <- sapply(stab.all, function(mm) mm$best.bin.dev)
  best.auc.each.rep <- sapply(stab.all, function(mm) mm$best.auc)

  best.l2.stable <- modal(best.l2.each.rep)
  best.bin.dev.stable <- modal(best.bin.dev.each.rep)
  best.auc.stable <- modal(best.auc.each.rep)

  list(ncv.loss = ncv.loss,
       best.l2.each.rep = best.l2.each.rep,
       best.bin.dev.each.rep = best.bin.dev.each.rep,
       best.auc.each.rep = best.auc.each.rep,
       best.l2.stable = best.l2.stable,
       best.bin.dev.stable = best.bin.dev.stable,
       best.auc.stable = best.auc.stable)

}
