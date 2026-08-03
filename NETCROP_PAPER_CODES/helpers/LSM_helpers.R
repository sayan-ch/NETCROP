library(Rcpp)
library(parallel)
library(Matrix)

################################################################################
## generates latent space model from Ma and Ma 2020 with some modifications
# For Hoff et al 2002 latent space model, supply alpha as a scalar

LSM.gen <- function (n, d, K = 1, alpha = NULL, avg.d = NULL,
                     ncore = 1, seed = NULL) {
  set.seed(seed)
  if(is.null(alpha))
    alpha <- -runif(n, 1, 3) / 2

  mu <- matrix(runif(K * d, -1, 1), K, d)
  idx <- sample(1:K, n, replace = TRUE)
  mu <- mu[idx, ]
  Z <- mu + matrix(truncnorm::rtruncnorm(n * d, -2, 2), n, d)
  Z <- scale(Z, scale = F)
  G <- Z %*% t(Z)
  Z <- Z/sqrt(sqrt(sum(G^2))/n)

  alpha <- alpha/2 - sqrt(rowSums(Z^2))/2

  G <- Z %*% t(Z)
  theta <- netcrop_outer(alpha, alpha, '+') + G
  P <- 1/(1 + exp(-theta))
  if (!is.null(avg.d)) {
    for (i in 1:10) {
      ratio <- mean(rowSums(P))/avg.d
      beta <- -log(ratio)
      alpha <- alpha + beta/2
      theta <- netcrop_outer(alpha, alpha, '+') + G
      P <- 1/(1 + exp(-theta))
    }
  }
  diag(P) <- 0

  stor <- do.call('rbind',
                  mclapply(1:(n-1), function(i) {
                    set.seed(seed + i)

                    tmp <- which(rbinom(n-i, 1, P[i,(i+1):n]) == 1)

                    if(length(tmp) == 0)
                      return(NULL)
                    else
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))

  A <- sparseMatrix(i = stor[,1], j = stor[,2], x = 1, dims = c(n,n))
  A <- A + t(A)

  return(list('A' = A, 'P' = P, 'Z' = Z, 'alpha' = alpha))
}

################################################################################
## LSM estimation via projected gradient descent
# the following Rcpp part is taken from randnet package and
# made some minor changes to format the outputs in a convenient way
Rcpp::sourceCpp(
  file.path(here::here("NETCROP_PAPER_CODES/helpers/LSM_PGD_Cpp.cpp")))

pgd.lsm <- function (A, d, step.size = 0.3, niter = 500, trace = 0){
  N <- nrow(A)
  Jmat <- diag(N) - 1/N

  P.tilde <- netcrop_usvt(A)
  P.tilde <- pmin(P.tilde, 1 - 1e-06)
  P.tilde <- pmax(P.tilde, 1e-06)
  Theta.tilde <- log(P.tilde / (1 - P.tilde))
  alpha_0 <- solve(diag(N, N) + 1, rowSums(Theta.tilde))

  G <- crossprod(Jmat, tcrossprod((Theta.tilde - netcrop_outer(alpha_0, alpha_0, "+")),
                                  Jmat))

  eig <- RSpectra::eigs_sym(A = G, k = d)

  eig$values <- pmax(eig$values, 0)
  Z_0 <- t(t(eig$vectors[, 1:d, drop = F]) * sqrt(eig$values[1:d]))

  step.size.z <- step.size/norm(Z_0, "2")^2
  step.size.alpha <- step.size/(2 * N)

  res <- LSM_PGD_Cpp(as.matrix(A), Z_0, alpha_0, step.size.z, step.size.alpha,
                     niter, trace = (trace > 0))
  return(res)
}

################################################################################
## NETCROP for LSM
# A: adjacency matrix, d.cand: candidate ranks,
# s: number of splits, o: number of overlapping nodes, R: number of repetitions,
# loss: loss function to be used for model selection,
# ncore: number of cores to be used for parallelization, seed: random seed for reproducibility, step.size and niter: parameters for projected gradient descent
netcrop_lsm <- function(A, d.cand, s, o, R,
                        loss = c("l2", "bin.dev", "AUC"),
                        ncore = 1, seed = NULL,
                        step.size = 0.3, niter = 500, trace = 0){
  set.seed(seed)
  n <- nrow(A)
  m <- (n-o)/s

  if(length(d.cand) == 1) d.cand <- 1:d.cand

  dmax <- max(d.cand)

  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o,
                                              replace = F))

  ld <- length(d.cand)
  raw.ind <- matrix(nrow = R*s*ld, ncol = 3)
  cc <- 1
  for(r in 1:R)
    for(dd in seq_along(d.cand))
      for(q in 1:s){
        raw.ind[cc, ] <- c(r, dd, q)
        cc <- cc + 1
      }
  colnames(raw.ind) <- c('r', 'dd', 'q')

  system.time(raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    q <- raw.ind[ii, 'q']
    r <- raw.ind[ii, 'r']
    dd <- raw.ind[ii, 'dd']

    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]

    out.lat <- pgd.lsm(A = A.sonn, d = d.cand[dd],
                       step.size = step.size, niter = niter, trace = trace)
    Z.sonn <- out.lat$Z
    beta.sonn <- out.lat$alpha

    return(list(Z.hat = Z.sonn, beta.hat = beta.sonn))

  },mc.cores = min(nrow(raw.ind), ncore)))

  match.out <- mclapply(1:nrow(raw.ind), function(ii){
    r <- raw.ind[ii, 'r']
    q <- raw.ind[ii, 'q']
    dd <- raw.ind[ii, 'dd']

    if(q == 1) return(list('Z.rot' = raw.out[[ii]]$Z.hat[-(1:o), , drop = F],
                           'beta.hat' = raw.out[[ii]]$beta.hat[-(1:o)]))

    stand <- which(raw.ind[,'r'] == r & raw.ind[,'q'] == 1 &
                     raw.ind[,'dd'] == dd)

    proc.par <- netcrop_procrustes(X = raw.out[[ii]]$Z.hat[(1:o), , drop = F],
                                   Xstar = raw.out[[stand]]$Z.hat[(1:o), , drop = F],
                                   translate = T,
                                   dilate = F)

    Z.rot0 <- raw.out[[ii]]$Z.hat[-(1:o), , drop = F] %*% proc.par$R
    Z.rot <- Z.rot0 + matrix(proc.par$t, nrow = m, ncol = d.cand[dd], byrow = T)

    # Z.rot <- scale(Z.rot, scale = F)

    t.mat <- matrix(proc.par$t)

    beta.rot <- raw.out[[ii]]$beta.hat[-(1:o)] +
      0.5 * sum(t.mat^2) + Z.rot0 %*% t.mat

    return(list('Z.rot' = Z.rot,
                'beta.hat' = beta.rot))
  }, mc.cores = min(nrow(raw.ind), ncore))

  non.size <- s*(s-1)/2
  ld <- length(d.cand)
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

  L.all <- lapply(1:nrow(non.mat), function(ii){
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

    ind1 <- which(raw.ind[, 'r'] == r & raw.ind[, 'q'] == p &
                    raw.ind[,'dd'] == dd)
    ind2 <- which(raw.ind[, 'r'] == r & raw.ind[, 'q'] == q &
                    raw.ind[,'dd'] == d.cand[dd])

    Z1.hat <- match.out[[ind1]]$Z.rot
    Z2.hat <- match.out[[ind2]]$Z.rot

    beta1.hat <- match.out[[ind1]]$beta.hat
    beta2.hat <- match.out[[ind2]]$beta.hat

    theta.hat <- netcrop_outer(beta1.hat, beta2.hat, '+') +
      tcrossprod(Z1.hat, Z2.hat)
    P.hat <- 1 / (1 + exp(-theta.hat))

    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, ] +
        ( do.call(loss[lq], list(A.non, P.hat)) +
            2*(o + m + dd) ) / non.size
    }

    return(L.temp)})

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

  return(c(list('loss' = obj),
           obj2))
}
