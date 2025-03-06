source("~/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT/all_base.R")
setwd("~/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT")

################################################################################
pair.NMI.loss <- function(g1, g2, g3, g4){
  
  # g3p <- matched.lab(lab = g3, fixed = g1)
  # g4p <- matched.lab(lab = g4, fixed = g2)
  
  pair1 <- expand.grid(g1 = g1, g2 = g2)
  
  pair11 <- pair1 %>% mutate(commu = (g1 - 1) * max(g2) + g2, 
                             id = 1:nrow(pair1)) %>%
    select(id, commu)
  
  pair2 <- expand.grid(g3 = g3, g4 = g4)
  pair22 <- pair2 %>% mutate(commu = (g3 - 1) * max(g4) + g4,
                             id = 1:nrow(pair2)) %>%
    select(id, commu)
  
  return(-randnet::NMI(pair11$commu, pair22$commu))
  # return(-randnet::NMI(pair22$commu, matched.lab(pair11$commu, pair22$commu)))
}

pair.hemming.loss <- function(g1, g2, g3, g4){
  
  # g3p <- matched.lab(lab = g3, fixed = g1)
  # g4p <- matched.lab(lab = g4, fixed = g2)
  
  pair1 <- expand.grid(g1 = g1, g2 = g2)
  
  pair11 <- pair1 %>% mutate(commu = (g1 - 1) * max(g2) + g2, 
                             id = 1:nrow(pair1)) %>%
    select(id, commu)
  
  pair2 <- expand.grid(g3 = g3, g4 = g4)
  pair22 <- pair2 %>% mutate(commu = (g3 - 1) * max(g4) + g4,
                             id = 1:nrow(pair2)) %>%
    select(id, commu)
  
  # mean(pair11$commu != pair22$commu)
  
  mean(matched.lab(pair11$commu, pair22$commu) != pair22$commu)
}

A = net$A; K = 5; tau.cand = c(0, 0.5, 1); DCBM = T; 
s = 2; o = 4000; R = 3; laplace = T; loss = c("pair.NMI.loss");
true.g = net$member; ncore = 8

# tune parameter
croissant.tune.regsp <- function(A, K, tau.cand,
                                 DCBM = F,
                                 s, o, R,
                                 laplace = F,
                                 loss = c("pair.NMI.loss"),
                                 true.g = NULL,
                                 ncore = 1){
  n <- nrow(A)
  m <- (n-o)/s
  
  # deg <- rowSums(A)
  # avg.deg <- mean(deg)
  
  L <- list()
  
  over <- lapply(1:R, function(ii) sample.int(n, o, F))
  non.over <- lapply(1:R, function(ii) sample((1:n)[-over[[ii]]], n-o, replace = F))
  
  raw.ind <- cbind(rep(1:R, each = s), rep(1:s, R))
  
  raw.out <- mclapply(1:nrow(raw.ind), function(ii){
    
    q <- raw.ind[ii, 2]
    r <- raw.ind[ii, 1]
    
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    # deg2 <- deg[sonn]
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
      
      # L.sonn <- as(A.sonn.tau, 'dMatrix')
      L.sonn <- A.sonn.tau + 1 - 1
      
      if(laplace){
        L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
        L.sonn[is.na(L.sonn)] <- 0
      }
      
      # eig.max <- irlba::partial_eigen(x = L.sonn, n = K,
      #                                 symmetric = T)$vectors
      eig.max <- RSpectra::eigs_sym(L.sonn, K, which = "LM",
                                    opts = list(v0 = rep(1, o+m)))$vectors
      
      rn.eig <- eig.max
      
      if(DCBM){
        rownorm <- sqrt(rowSums(eig.max^2))
        rownorm[rownorm == 0] <- 1
        
        rn.eig <- eig.max/rownorm
        
        out.BM[[tt]] <- as.integer(pam(rn.eig, K,
                                       metric = "euclidean",
                                       do.swap = F, cluster.only = T,
                                       pamonce = 6))
      }else{
        out.BM[[tt]] <- as.integer(kmeans(eig.max, K, nstart = 100,
                                          iter.max = 10^7)$cluster)
      }
      
      
      
      # out.BM[[tt]] <- as.integer(kmeans(rn.eig, K, nstart = 100,
      #                                   iter.max = 10^7)$cluster)
      }
      if(tau.cand[tt] == 0){
        out.BM[[tt]] <- vector(length = o+m)
        # deg2 <- deg[sonn]
        
        bad.node <- which(deg2 == 0)
        out.BM[[tt]][bad.node] <- sample(1:K, length(bad.node), replace = T)
        
        good.node <- which(deg2 > 0)
        
        # L.sonn.good <- as(A.sonn[good.node, good.node], 'dMatrix')
        L.sonn.good <- A.sonn[good.node, good.node] + 1 - 1
        d.sonn.good <- sparseMatrix( i = 1:length(good.node), j = 1:length(good.node),
                                    x = 1/sqrt(deg2[good.node]))
        
        if(laplace){
          L.sonn.good <- tcrossprod(crossprod(d.sonn.good, L.sonn.good), 
                                    d.sonn.good)
          L.sonn.good[is.na(L.sonn.good)] <- 0
        }
        
        # eig.max.good <- irlba::partial_eigen(x = L.sonn.good, n = K,
        #                                      symmetric = T)$vectors
                
        eig.max.good <- RSpectra::eigs_sym(L.sonn.good, K, which = "LM",
                          opts = list(v0 = rep(1, length(good.node))))$vectors
        
        rn.eig.good <- eig.max.good
        
        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max.good^2))
          rownorm[rownorm == 0] <- 1
          
          rn.eig.good <- eig.max.good/rownorm
          
          out.BM[[tt]][good.node] <- as.integer(pam(rn.eig.good, K,
                                                    metric = "euclidean",
                                                    do.swap = F, cluster.only = T,
                                                    pamonce = 6))
        }else{
          out.BM[[tt]][good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                                      iter.max = 10^7)$cluster)
        }
        
        
      }
    }
    
    return(list('BM' = out.BM))
    # 'psi' = psi.hat))
  },
  mc.cores = ncore)
  
  tau.size <- length(tau.cand)
  
  est.out <- mclapply(1:(tau.size*nrow(raw.ind)), function(ii){
    
    tt <- ii %% tau.size
    tt <- ifelse(tt == 0, tau.size, tt)
    
    rot <- ceiling(ii / tau.size)
    
    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]
    
    #message(paste0("Est. at s=",q, " started"))
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    
    out.BM.std <- raw.out[raw.ind[,1] == r][[1]]$BM[[tt]]
    
    out.BM <- raw.out[raw.ind[,1] == r][[q]]$BM[[tt]]
    
    work.tau <- tau.cand[tt]
    
    if(K == 1){
      mat.BM <- rep(1, m)
      
      return(list('gBM' = mat.BM))
    }
    
    E.BM.kc <- best.perm.label.match(out.BM[1:o], 
                                     out.BM.std[1:o],
                                     o, K)
    
    # tmp.BM <- sparseMatrix(i = 1:m, j = out.BM[(o+1):(o+m)], 
    #                        dims = c(m, K))
    tmp.BM <- sparseMatrix(i = 1:(o+m), j = out.BM, dims = c((o+m), K))
    
    mat.BM <- as.vector(tcrossprod(tcrossprod(tmp.BM, E.BM.kc),
                                   rbind(1:K)))
    
    message(paste0("Est. at s=",q, " finished"))
    
    return(list('gBM' = mat.BM[(o+1):(o+m)], 'ggBM' = mat.BM))
  },
  mc.cores = ncore)
  
  g.BM <- list()
  
  gg.BM <- list()
  
  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = tau.size), ], 
                   rep(1:tau.size, nrow(raw.ind)))
  
  for(r in 1:R){
    g.BM[[r]] <- list()
    gg.BM[[r]] <- list()
    
    for(tt in seq_along(tau.cand)){
      tmp.est <- est.out[which(raw.mat[,3] == tt & raw.mat[,1] == r)]
      g.BM[[r]][[tt]] <- list()
      
      gg.BM[[r]][[tt]] <- vector(length = n)
      
      gg.BM[[r]][[tt]][over[[r]]] <- tmp.est[[1]]$ggBM[1:o]
      
      for(q in 1:s){
        g.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$gBM
        
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
    
    # deg <- rowSums(A)
    # avg.deg <- mean(deg)
      if(tau.cand[tt] > 0){
      A.tau <- A + tau.cand[tt]*avg.deg/n
      
      d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                                  x = 1/sqrt(deg + tau.cand[tt]*avg.deg))
      
      L.tau <- A.tau + 1 - 1
      
      if(laplace){
        L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
        L.tau[is.na(L.tau)] <- 0
      }
      
      # eig.max <- irlba::partial_eigen(x = L.tau, n = K,
      #                                 symmetric = T)$vectors
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
        
        out.BM <- as.integer(pam(rn.eig, K,
                                 metric = "euclidean",
                                 do.swap = F, cluster.only = T,
                                 pamonce = 6))
        return(out.BM)
      }else{
        out.BM <- as.integer(kmeans(eig.max, K, nstart = 100,
                                    iter.max = 10^7)$cluster)
        return(out.BM)
      }
      
      
      
      # out.BM <- as.integer(kmeans(rn.eig, K, nstart = 100,
      #                                   iter.max = 10^7)$cluster)
    
    # return(list('BM' = out.BM))
      }
      
      if(tau.cand[tt] == 0){
        out.BM <- vector(length = n)
        
        bad.node <- which(deg == 0)
        out.BM[bad.node] <- sample(1:K, length(bad.node), replace = T)
        
        good.node <- which(deg > 0)
        
        L.good <- A[good.node, good.node] + 1 - 1
        d.good <- sparseMatrix( i = 1:length(good.node), j = 1:length(good.node),
                                     x = 1/sqrt(deg[good.node]))
        
        if(laplace){
          L.good <- tcrossprod(crossprod(d.good, L.good), d.good)
          L.good[is.na(L.good)] <- 0
        }
        
        # eig.max.good <- irlba::partial_eigen(x = L.good, n = K,
        #                                      symmetric = T)$vectors
        
        eig.max.good <- RSpectra::eigs_sym(L.good, K, which = "LM",
                          opts = list(v0 = rep(1, length(good.node))))$vectors
        
        rn.eig.good <- eig.max.good
        
        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max.good^2))
          rownorm[rownorm == 0] <- 1
          
          rn.eig.good <- eig.max.good/rownorm
          
          out.BM[good.node] <- as.integer(pam(rn.eig.good, K,
                                                    metric = "euclidean",
                                                    do.swap = F, cluster.only = T,
                                                    pamonce = 6))
        }else{
          out.BM[good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                                       iter.max = 10^7)$cluster)
        }
        
        return(out.BM)
      }
  },
  mc.cores = ncore)
  
  
  
  L.all <- mclapply(1:nrow(non.mat), function(ii){
    r <- non.mat[ii, 1]
    tt <- non.mat[ii, 2]
    p <- non.mat[ii, 3]
    q <- non.mat[ii, 4]
    
    p.non <- non.over[[r]][((p-1)*m+1):(p*m)]
    q.non <- non.over[[r]][((q-1)*m+1):(q*m)]
    
    L.temp <- matrix(0, nrow = length(loss), ncol = 1)
    row.names(L.temp) <- loss
    colnames(L.temp) <- as.character(tau.cand[tt])
    
    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        do.call(loss[lq], list(g.BM[[r]][[tt]][[p]], 
                               g.BM[[r]][[tt]][[q]], 
                               all.out[[tt]][p.non], 
                               all.out[[tt]][q.non]))
    }
    
    return(L.temp)},
    mc.cores = ncore)
  
  L <- list()
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:tau.size, function(tt){
      Reduce('+', L.all[which(non.mat[,2] == tt & non.mat[,1] == r)])
    }))
    
    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(tau.cand)
  }
  
  obj <- data.table(`Candidate Tau` = tau.cand)
  
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
      
      out.BM <- as.integer(pam(rn.eig, K,
                                     metric = "euclidean",
                                     do.swap = F, cluster.only = T,
                                     pamonce = 6))
      }
      # out.BM <- as.integer(kmeans(rn.eig, K, nstart = 100,
      #                             iter.max = 10^7)$cluster)
      
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

l2.tune.regsp <- function(A, K, tau.cand,
                                 DCBM = F,
                                 s, o, R,
                                 laplace = F,
                                 dc.est = 2,
                                 loss = c("l2", "bin.dev", "AUC"),
                                 ncore = 1){
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
    
    deg <- rowSums(A.sonn)
    avg.deg <- mean(deg)
    
    out.BM <- list()
    
    for(tt in seq_along(tau.cand)){
      if(K == 1){
        out.BM[[tt]] <- rep(1, o+m)
        next
      }
      
      if(tau.cand[tt] > 0){
        A.sonn.tau <- A.sonn + tau.cand[tt]*avg.deg/(o+m)
        
        d.sonn.tau <- sparseMatrix( i = 1:(o+m), j = 1:(o+m),
                                    x = 1/sqrt(deg + tau.cand[tt]*avg.deg))
        
        L.sonn <- A.sonn.tau + 1 - 1
        
        if(laplace){
          L.sonn <- tcrossprod(crossprod(d.sonn.tau, A.sonn.tau), d.sonn.tau)
          L.sonn[is.na(L.sonn)] <- 0
        }
        
        # eig.max <- irlba::partial_eigen(x = L.sonn, n = K,
        #                                 symmetric = T)$vectors
        eig.max <- RSpectra::eigs_sym(L.sonn, K, which = "LM",
                                      opts = list(v0 = rep(1, o+m)))$vectors
        
        rn.eig <- eig.max
        
        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max^2))
          rownorm[rownorm == 0] <- 1
          
          rn.eig <- eig.max/rownorm
          
          out.BM[[tt]] <- as.integer(pam(rn.eig, K,
                                         metric = "euclidean",
                                         do.swap = F, cluster.only = T,
                                         pamonce = 6))
        }else{
          out.BM[[tt]] <- as.integer(kmeans(eig.max, K, nstart = 100,
                                            iter.max = 10^7)$cluster)
        }
        
        
        
        # out.BM[[tt]] <- as.integer(kmeans(rn.eig, K, nstart = 100,
        #                                   iter.max = 10^7)$cluster)
      }
      if(tau.cand[tt] == 0){
        out.BM[[tt]] <- vector(length = o+m)
        
        bad.node <- which(deg == 0)
        out.BM[[tt]][bad.node] <- sample(1:K, length(bad.node), replace = T)
        
        good.node <- which(deg > 0)
        
        L.sonn.good <- A.sonn[good.node, good.node] + 1 - 1
        d.sonn.good <- sparseMatrix( i = 1:length(good.node), j = 1:length(good.node),
                                     x = 1/sqrt(deg[good.node]))
        
        if(laplace){
          L.sonn.good <- tcrossprod(crossprod(d.sonn.good, L.sonn.good), 
                                    d.sonn.good)
          L.sonn.good[is.na(L.sonn.good)] <- 0
        }
        
        # eig.max.good <- irlba::partial_eigen(x = L.sonn.good, n = K,
        #                                      symmetric = T)$vectors
        
        eig.max.good <- RSpectra::eigs_sym(L.sonn.good, K, which = "LM",
                          opts = list(v0 = rep(1, length(good.node))))$vectors
        
        rn.eig.good <- eig.max.good
        
        if(DCBM){
          rownorm <- sqrt(rowSums(eig.max.good^2))
          rownorm[rownorm == 0] <- 1
          
          rn.eig.good <- eig.max.good/rownorm
          
          out.BM[[tt]][good.node] <- as.integer(pam(rn.eig.good, K,
                                                    metric = "euclidean",
                                                    do.swap = F, cluster.only = T,
                                                    pamonce = 6))
        }else{
          out.BM[[tt]][good.node] <- as.integer(kmeans(eig.max.good, K, nstart = 100,
                                                       iter.max = 10^7)$cluster)
        }
        
        
      }
    }
    
    return(list('BM' = out.BM))
    # 'psi' = psi.hat))
  },
  mc.cores = ncore)
  
  tau.size <- length(tau.cand)
  
  est.out <- mclapply(1:(tau.size*nrow(raw.ind)), function(ii){
    
    tt <- ii %% tau.size
    tt <- ifelse(tt == 0, tau.size, tt)
    
    rot <- ceiling(ii / tau.size)
    
    q <- raw.ind[rot, 2]
    r <- raw.ind[rot, 1]
    
    #message(paste0("Est. at s=",q, " started"))
    sonn <- c(over[[r]], non.over[[r]][((q-1)*m+1):(q*m)])
    A.sonn <- A[sonn, sonn]
    
    out.BM.std <- raw.out[raw.ind[,1] == r][[1]]$BM[[tt]]
    
    out.BM <- raw.out[raw.ind[,1] == r][[q]]$BM[[tt]]
    
    work.tau <- tau.cand[tt]
    
    if(K == 1){
      mat.BM <- rep(1, m)
      
      if(!DCBM){
        B.BM <- fast.SBM.est(A.sonn, rep(1,o+m), o+m, 1)
        mat.BM <- rep(1, m)
        psi.BM <- rep(1, m)
        
        return(list('gBM' = mat.BM, 'BBM' = B.BM,
                    'psiBM' = psi.BM))
      }
      
      if(DCBM){
        if(dc.est == 2){
          tmp <- fast.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                               p.sample = 1)
        }else{
          #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
          tmp <- eigen.DCBM.est(A.sonn, rep(1,o+m), o+m, 1, o,
                                p.sample = 1)
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
    
    tmp.BM <- sparseMatrix(i = 1:(o+m), j = out.BM, 
                           dims = c((o+m), K))
    
    mat.BM <- as.vector(tcrossprod(tcrossprod(tmp.BM, E.BM.kc),
                                   rbind(1:K)))
    if(!DCBM){
      B.BM <- fast.SBM.est(A.sonn, mat.BM, o+m, K)
      mat.BM <- mat.BM[-(1:o)]
      psi.BM <- rep(1, m)
      
      return(list('gBM' = mat.BM, 'BBM' = B.BM,
                  'psiBM' = psi.BM))
    }
    
    if(dc.est == 2){
      tmp <- fast.DCBM.est(A.sonn, mat.BM, o+m, K, o, 
                           p.sample = 1)
    }else{
      #psi.hat <- raw.out[raw.ind[,1] == r][[q]]$psi[[k.cand]]
      tmp <- eigen.DCBM.est(A.sonn, mat.BM, o+m, K, o, 
                            p.sample = 1)
    }
    
    B.BM <- tmp$Bsum
    psi.BM <- tmp$psi
    mat.BM <- mat.BM[-(1:o)]
    
    message(paste0("Est. at s=",q, " finished"))
    
    return(list('gBM' = mat.BM,
                'BBM' = B.BM,
                'psiBM' = psi.BM))
  },
  mc.cores = ncore)
  
  g.BM <- list()
  B.BM <- list()
  psi.BM <- list()
  
  raw.mat <- cbind(raw.ind[rep(1:nrow(raw.ind), each = tau.size), ], 
                   rep(1:tau.size, nrow(raw.ind)))
  
  for(r in 1:R){
    g.BM[[r]] <- list()
    B.BM[[r]] <- list()
    psi.BM[[r]] <- list()
    
    for(tt in seq_along(tau.cand)){
      tmp.est <- est.out[which(raw.mat[,3] == tt & raw.mat[,1] == r)]
      B.BM[[r]][[tt]] <- 0
      g.BM[[r]][[tt]] <- list()
      psi.BM[[r]][[tt]] <- list()
      
      for(q in 1:s){
        B.BM[[r]][[tt]] <- B.BM[[r]][[tt]] + 
          tmp.est[[q]]$BBM/s
        
        g.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$gBM
        
        psi.BM[[r]][[tt]][[q]] <- tmp.est[[q]]$psiBM
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
    
    P.BM <- B.BM[[r]][[tt]][g.BM[[r]][[tt]][[p]],
                            g.BM[[r]][[tt]][[q]] ] *
      tcrossprod(psi.BM[[r]][[tt]][[p]],
                 psi.BM[[r]][[tt]][[q]])
    
    P.BM[P.BM < 1e-6] <- 1e-6
    P.BM[P.BM > 1- 1e-6] <- 1 - 1e-6
    
    for(lq in seq_along(loss)){
      tmp.nm <- loss[lq]
      L.temp[tmp.nm, 1] <-  L.temp[tmp.nm, 1] +
        do.call(loss[lq], list(as.numeric(A.non), P.BM))/(s*(s-1)*0.5)
    }
    
    return(L.temp)},
    mc.cores = ncore)
  
  L <- list()
  
  for(r in 1:R){
    L[[r]] <- do.call('cbind', lapply(1:tau.size, function(tt){
      Reduce('+', L.all[which(non.mat[,2] == tt & non.mat[,1] == r)])
    }))
    
    row.names(L[[r]]) <- loss
    colnames(L[[r]]) <- as.character(tau.cand)
  }
  
  obj <- data.table(`Candidate Tau` = tau.cand)
  
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
  
  return(c(list('loss' = obj), 
           obj2))
}

DKest <- function(A, K, true.g = NULL, tau.cand, DCBM = T, DC.est = 2,
                  ncore = 1){
  n <- nrow(A)
  
  deg <- rowSums(A)
  
  DK.out <- mclapply(seq_along(tau.cand), function(tt){
    A.tau <- A + tau.cand[tt]*mean(deg)/n
    
    d.tau <- sparseMatrix( i = 1:n, j = 1:n,
                          x = 1/sqrt(deg + tau.cand[tt]*mean(deg)/n))
    
    L.tau <- tcrossprod(crossprod(d.tau, A.tau), d.tau)
    
    eig.max <- irlba::partial_eigen(x = L.tau, n = K,
                                    symmetric = T)$vectors
    
    rn.eig <- eig.max
    
    if(DCBM){
      rownorm <- sqrt(rowSums(eig.max^2))
      rownorm[rownorm == 0] <- 1
      
      rn.eig <- eig.max/rownorm
    }
    
    out.comm <- as.integer(kmeans(rn.eig, K, nstart = 100,
                                  iter.max = 10^7)$cluster)
    
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
    
    P.hat[P.hat > 1] <- 1
    P.hat[P.hat < 0] <- 0
    
    deg.hat <- rowSums(P.hat)
    P.hat.tau <- P.hat + tau.cand[tt]*mean(deg)/n
    
    d.tau.hat <- sparseMatrix( i = 1:n, j = 1:n,
                              x = 1/sqrt(deg.hat + tau.cand[tt]*mean(deg)/n))
    
    L.hat.tau <- tcrossprod(crossprod(d.tau.hat, P.hat.tau), d.tau.hat)
    
    numerator <- tryCatch({irlba::irlba(A = abs(L.tau - L.hat.tau), 
                                            nv = 1)$d},
                          error = function(e){
                            return(1)
                          }
                          )
    
    denominator <- tryCatch({irlba::partial_eigen(x = L.hat.tau, n = K,
                                        symmetric = T)$values[K]},
                            error = function(e){
                              return(0)
                            }
                            )
    
    return(c(tau.cand = tau.cand[tt], DK.stat = abs(numerator/denominator),
                accu = mean(matched.lab(out.comm, g) == g)
                ))
    
  }, mc.cores = ncore)
  
  DK.out <- do.call('rbind', DK.out)
  
}


net <- randnet::BlockModel.Gen(lambda = 5, n = 600, beta = 0.2,
                               K = 3, rho = 1, simple = T, power = T,
                               alpha = 5)

DK.out <- DKest(net$A, 3, net$g, tau.cand = 0.1*(0:20), DCBM = T,
                DC.est = 1, ncore = 8)

#### add conditional on dcbm to choose from kmeans or pam


##################################################
# psi.tmp <- poweRlaw::rplcon(n = 300, xmin = 1, alpha = 5)
psi.tmp <- rbeta(300, 1, 4)
psi2 <- sample(psi.tmp, 1000, replace = T)
  
net <- blockmodel.gen.fast(n = 1000, K = 5, 
                           B = matrix(0.1, 5, 5) + diag(0.2, 5),
                           psi = psi2, ncore = 6)
mean(rowSums(net$A))

time.out <- system.time(
  tune.out <- netcrop.tune.regsp(A = net$A, K = 5, 
                                    tau.cand = 0.1*(0:20),
                                    DCBM = T,
                                    s = 5, o = 500, R = 5,
                                    dc.est = 2,
                                    laplace = T,
                                    loss = c("pair.hemming.loss", 
                                             "pair.NMI.loss", "l2", "bin.dev",
                                             "AUC"),
                                    true.g = net$member,
                                    ncore = 6 ))[3]


system.time(
  dk.out <- DKest(A = net$A, K = 5,
                  true.g = net$member, tau.cand = 0.1*(0:20), laplace = T,
                  DCBM = T,
                  DC.est = 2, ncore = 6)
)

dk.out
which.min(dk.out[,2])

final.out <- list()
final.time <- vector()

psi.tmp <- rbeta(300, 1, 4)
psi2 <- sample(psi.tmp, 10000, replace = T)
deg.out <- vector()
n <- 10000; K <- 5
B <- matrix(0.1, K, K) + diag(0.2, K)
ncore <- detectCores() - 1

for(enum in 1:100){
  net <- blockmodel.gen.fast(n = 10000, K = 5, B = B,
                             psi = psi2, ncore = ncore)
  
  deg.out[enum] <- mean(rowSums(net$A))
  final.time[enum] <- system.time(
    final.out[[enum]] <- croissant.tune.regsp(A = net$A, K = K, 
                                      tau.cand = 0.1*(0:20),
                                      DCBM = T,
                                      s = 2, o = 4000, R = 5,
                                      laplace = T,
                                      loss = c("pair.hemming.loss", 
                                               "pair.NMI.loss"),
                                      true.g = net$member,
                                      ncore = ncore ))[3]
  
  cat(enum, " done\n")
  
  save.image(file = "~/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT/temp_reg_out.RData")  

}

################################################################################
## Plotting

nsim <- length(final.out)

clear.out <- list()
mat.out <- list()

for(ii in 1:nsim){
  clear.out[[ii]] <- data.table(
      tau = c(as.character(final.out[[ii]]$obj2$`Candidate Tau`),
              "\u03c4",
              paste0("\u03c4", "(mean)"),
              paste0("\u03c4", "(mode)")),
      NMI.accuracy = 100*c(final.out[[ii]]$all.accu,
                   final.out[[ii]]$croissant.all.accu["pair.NMI.loss"],
                   final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mean"],
                   final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mode"]),
      Hemming.accuracy = 100*c(final.out[[ii]]$all.accu,
                   final.out[[ii]]$croissant.all.accu["pair.hemming.loss"],
                   final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mean"],
                   final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mode"])
    )
  
  mat.out[[ii]] <- bind_cols(
    NMI = as.numeric(100*c(final.out[[ii]]$all.accu,
          final.out[[ii]]$croissant.all.accu["pair.NMI.loss"],
          final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mean"],
          final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mode"])),
    Hemming = as.numeric(100*c(final.out[[ii]]$all.accu,
          final.out[[ii]]$croissant.all.accu["pair.hemming.loss"],
          final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mean"],
          final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mode"]))
  )
                               
}

mat.mean <- Reduce('+', mat.out)/nsim
mat.sd <- sqrt((Reduce('+', lapply(mat.out, function(x) x^2)) - 
             nsim*mat.mean^2)/(nsim-1))

plot.out <- data.table(
  tau = rep(clear.out[[1]]$tau, 2),
  accuracy = c(mat.mean[,"NMI"], mat.mean[,"Hemming"]),
  sd = c(mat.sd[,"NMI"], mat.sd[,"Hemming"]),
  loss = rep(c("Negative Pairwise NMI", "Pairwise Hemming Loss"),
             each = length(clear.out[[1]]$tau)),
  method = rep(c(
    rep("Pre-fixed", length(clear.out[[1]]$tau) - 3),
    "NETCROP \u03c4",
    "Mean of NETCROP \u03c4",
    "Mode of NETCROP \u03c4"
  ), 2)
)


# error bar plot
library(ggplot2)
library(ggpubr)

plot.out %>%
  filter(loss == "Pairwise Hemming Loss") %>%
  ggplot(aes(x = tau, y = accuracy)) +
  geom_point(aes(color = method), size = 1.5) +
  geom_errorbar(aes(ymin = pmax(accuracy - sd, 0),
                    ymax = pmin(accuracy + sd, 100),
                    color = method)
                ) +
  theme_pubclean()

################################################################################
## new plotting
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggtext)
library(extrafont)
library(RColorBrewer)

font_import()  # This may take a little time
# loadfonts(device = "win")  # On Windows
loadfonts(device = "pdf")  # On Mac or Linux

load('/Users/sayanchakrabarty/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT/partune_clus_out/partune/l2/l2_final_par_tune_out.RData')
load('/Users/sayanchakrabarty/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT/partune_clus_out/partune/DK/final_DK_par_tune_out.RData')

nsim <- length(final.out)

mat.out <- list()

for(ii in 1:nsim){

  
    
  mat.out[[ii]] <- 100*c(
    final.out[[ii]]$all.accu[1],
    max(final.out[[ii]]$all.accu),
    dk.out[[ii]][,'accu'][which.min(dk.out[[ii]][,'DK.stat'])],
    final.out[[ii]]$croissant.all.accu["l2"],
    final.out[[ii]]$croissant.all.accu["l2.mean"],
    final.out[[ii]]$croissant.all.accu["l2.mode"],
    final.out[[ii]]$croissant.all.accu["bin.dev"],
    final.out[[ii]]$croissant.all.accu["bin.dev.mean"],
    final.out[[ii]]$croissant.all.accu["bin.dev.mode"],
    final.out[[ii]]$croissant.all.accu["AUC"],
    final.out[[ii]]$croissant.all.accu["AUC.mean"],
    final.out[[ii]]$croissant.all.accu["AUC.mode"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mean"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mode"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mean"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mode"]
  )
  
  names(mat.out[[ii]]) <- c("0", "Oracle", 
            "Davis-Kahan Estimator",
            "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode",
            "NETCROP(bd)", "NETCROP(bd)-Mean", "NETCROP(bd)-Mode",
            "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode",
            "NETCROP(NMI)", "NETCROP(NMI)-Mean", "NETCROP(NMI)-Mode",
            "NETCROP(Hemming)", "NETCROP(Hemming)-Mean", "NETCROP(Hemming)-Mode"
            )
  
}


mat.all <- do.call(cbind, mat.out)

mat.median <- apply(mat.all, 1, mean)
mat.med.sd <- apply(mat.all, 1, sd)

plot.out <- data.table(
  tau = factor(rownames(mat.all),
               levels = c("0", "Oracle", 
                          "Davis-Kahan Estimator",
                          "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode",
                          "NETCROP(bd)", "NETCROP(bd)-Mean", 
                          "NETCROP(bd)-Mode",
                          "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode",
                          "NETCROP(NMI)", "NETCROP(NMI)-Mean",
                          "NETCROP(NMI)-Mode", "NETCROP(Hemming)",
                          "NETCROP(Hemming)-Mean", "NETCROP(Hemming)-Mode")),
  accuracy = mat.median,
  sd = mat.med.sd,
  xxmin = mat.median - mat.med.sd,
  xxmax = mat.median + mat.med.sd
)

setwd('/Users/sayanchakrabarty/Library/CloudStorage/OneDrive-Umich/From UIUC/Research/Model Selection Codes/CROISSANT/partune_clus_out')
#make a directory named plots if it already doesn't exist
if(!dir.exists("plots")) dir.create("plots")

## black and white
bw.all <- plot.out %>%
  ggplot(aes(x = accuracy, y = tau)) +
  # geom_point(aes(color = tau), size = 1.5) +
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(limits = rev(levels(plot.out$tau)),
                   labels = c(
                     "0" = "0",
                     "Oracle" = "Oracle",
                     "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
                     "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
                     "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
                     "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode")),
                     "NETCROP(bd)" = expression(paste("NETCROP(", italic(l[bin]), ")")),
                     "NETCROP(bd)-Mean" = expression(paste("NETCROP(", italic(l[bin]), ")-Mean")),
                     "NETCROP(bd)-Mode" = expression(paste("NETCROP(", italic(l[bin]), ")-Mode")),
                     "NETCROP(AUC)" = "NETCROP(AUC)",
                     "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
                     "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode",
                     "NETCROP(NMI)" = "NETCROP(NMI)",
                     "NETCROP(NMI)-Mean" = "NETCROP(NMI)-Mean",
                     "NETCROP(NMI)-Mode" = "NETCROP(NMI)-Mode",
                     "NETCROP(Hemming)" = "NETCROP(Hamming)",
                     "NETCROP(Hemming)-Mean" = "NETCROP(Hamming)-Mean",
                     "NETCROP(Hemming)-Mode" = "NETCROP(Hamming)-Mode"
                   )) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
       y = expression(tau)) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold"),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

bw.all

ggsave("plots/all_bw.png", bw.all, width = 8, height = 4.5, dpi = 600, 
       device = "png")

## color version
# Creating a new variable for the color groups
plot.out$group <- factor(
  ifelse(plot.out$tau %in% c("0", "Oracle", "Davis-Kahan Estimator"), "Group 1",  # Group for 0, Oracle, Davis-Kahan
         ifelse(plot.out$tau %in% c("NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode"), "Group 2",  # Group for l2
                ifelse(plot.out$tau %in% c("NETCROP(bd)", "NETCROP(bd)-Mean", "NETCROP(bd)-Mode"), "Group 3",  # Group for bd
                       ifelse(plot.out$tau %in% c("NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode"), "Group 4",  # Group for AUC
                              ifelse(plot.out$tau %in% c("NETCROP(NMI)", "NETCROP(NMI)-Mean", "NETCROP(NMI)-Mode"), "Group 5",  # Group for NMI
                                     "Group 6"))))))  # Group for Hemming

color.all <- plot.out %>%
  ggplot(aes(x = accuracy, y = tau, color = group)) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.out$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode")),
      "NETCROP(bd)" = expression(paste("NETCROP(", italic(l[bin]), ")")),
      "NETCROP(bd)-Mean" = expression(paste("NETCROP(", italic(l[bin]), ")-Mean")),
      "NETCROP(bd)-Mode" = expression(paste("NETCROP(", italic(l[bin]), ")-Mode")),
      "NETCROP(AUC)" = "NETCROP(AUC)",
      "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
      "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode",
      "NETCROP(NMI)" = "NETCROP(NMI)",
      "NETCROP(NMI)-Mean" = "NETCROP(NMI)-Mean",
      "NETCROP(NMI)-Mode" = "NETCROP(NMI)-Mode",
      "NETCROP(Hemming)" = "NETCROP(Hamming)",
      "NETCROP(Hemming)-Mean" = "NETCROP(Hamming)-Mean",
      "NETCROP(Hemming)-Mode" = "NETCROP(Hamming)-Mode"
    )
  ) +
  scale_color_manual(
    values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
               "Group 2" = "red",   # Color for Group 2 (l2)
               "Group 3" = "green", # Color for Group 3 (bd)
               "Group 4" = "purple",# Color for Group 4 (AUC)
               "Group 5" = "orange",# Color for Group 5 (NMI)
               "Group 6" = "brown") # Color for Group 6 (Hemming)
  ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  ) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               color = c(
                                 rep("brown", 3),
                                 rep("orange", 3),
                                 rep("purple", 3),   
                                 rep("darkgreen", 3),
                                 rep("darkred", 3),
                                 rep("darkblue", 3)
                               )),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )
  
color.all

ggplot2::ggsave("plots/all_color.png", color.all, width = 8, height = 4.5, 
                dpi = 600, device = "png")


################################################################################
## only 0, oracle, l2, bd, auc
plot.simple <- plot.out %>% filter(
  tau != tau[grepl("NMI|Hemming", tau)]
) %>%
  mutate(
    tau = factor(tau,
                 levels = c("0", "Oracle", 
                            "Davis-Kahan Estimator",
                            "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode",
                            "NETCROP(bd)", "NETCROP(bd)-Mean", 
                            "NETCROP(bd)-Mode",
                            "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode"
                            )
  ))

plot.simple$group <- factor(
  ifelse(plot.simple$tau %in% c("0", "Oracle", "Davis-Kahan Estimator"), "Group 1",  # Group for 0, Oracle, Davis-Kahan
         ifelse(plot.simple$tau %in% c("NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode"), "Group 2",  # Group for l2
                ifelse(plot.simple$tau %in% c("NETCROP(bd)", "NETCROP(bd)-Mean", "NETCROP(bd)-Mode"), "Group 3",  # Group for bd
                       "Group 4")))  # Group for AUC
)

bw.simple <- plot.simple %>%
  ggplot(aes(x = accuracy, y = tau, 
             # color = group
             )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.simple$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode"),
      "NETCROP(bd)" = expression(paste("NETCROP(", italic(l[bin]), ")")),
      "NETCROP(bd)-Mean" = expression(paste("NETCROP(", italic(l[bin]), ")-Mean")),
      "NETCROP(bd)-Mode" = expression(paste("NETCROP(", italic(l[bin]), ")-Mode")),
      "NETCROP(AUC)" = "NETCROP(AUC)",
      "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
      "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
    )
  )) +
  # scale_color_manual(
  #   values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
  #              "Group 2" = "red",   # Color for Group 2 (l2)
  #              "Group 3" = "green", # Color for Group 3 (bd)
  #              "Group 4" = "purple") # Color for Group 4 (AUC)
  # ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  ) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               # color = c(
                               #   rep("brown", 3),
                               #   rep("orange", 3),
                               #   rep("purple", 3),   
                               #   rep("darkgreen", 3),
                               #   rep("darkred", 3),
                               #   rep("darkblue", 3)
                               # )
                               ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

bw.simple  

ggsave("plots/simple_bw.png", bw.simple, width = 8, height = 4.5, dpi = 600, 
       device = "png")

color.simple <- plot.simple %>%
  ggplot(aes(x = accuracy, y = tau, 
             color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.simple$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode"),
                                      "NETCROP(bd)" = expression(paste("NETCROP(", italic(l[bin]), ")")),
                                      "NETCROP(bd)-Mean" = expression(paste("NETCROP(", italic(l[bin]), ")-Mean")),
                                      "NETCROP(bd)-Mode" = expression(paste("NETCROP(", italic(l[bin]), ")-Mode")),
                                      "NETCROP(AUC)" = "NETCROP(AUC)",
                                      "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
                                      "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
      )
    )) +
  scale_color_manual(
    values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
               "Group 2" = "red",   # Color for Group 2 (l2)
               "Group 3" = "green", # Color for Group 3 (bd)
               "Group 4" = "purple") # Color for Group 4 (AUC)
  ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  ) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               color = c(
                                 rep("purple", 3),   
                                 rep("darkgreen", 3),
                                 rep("darkred", 3),
                                 rep("darkblue", 3)
                               )),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

color.simple

ggsave("plots/simple_color.png", color.simple, width = 8, height = 4.5, dpi = 600, 
       device = "png")


################################################################################
## only 0, oracle, l2
plot.l2 <- plot.out %>% filter(
  tau %notin% tau[grepl("NMI|Hemming|bd|AUC", tau)]
) %>%
  mutate(
    tau = factor(tau,
                 levels = c("0", "Oracle", 
                            "Davis-Kahan Estimator",
                            "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode"
                            )
  ))

plot.l2$group <- factor(
  if_else(
    plot.l2$tau %in% c("0", "Oracle", "Davis-Kahan Estimator"), 
    "Group 1", "Group 2"
    )
)

bw.l2 <- plot.l2 %>%
  ggplot(aes(x = accuracy, y = tau, 
             # color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.l2$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode"))
    )) +
  # scale_color_manual(
  #   values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
  #              "Group 2" = "red",   # Color for Group 2 (l2)
  #              "Group 3" = "green", # Color for Group 3 (bd)
  #              "Group 4" = "purple") # Color for Group 4 (AUC)
  # ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  ) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               # color = c(
                               #   rep("brown", 3),
                               #   rep("orange", 3),
                               #   rep("purple", 3),   
                               #   rep("darkgreen", 3),
                               #   rep("darkred", 3),
                               #   rep("darkblue", 3)
                               # )
                               ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

bw.l2 

ggsave("plots/l2_bw.png", bw.l2, width = 8, height = 4.5, dpi = 600, 
       device = "png")

color.l2 <- plot.l2 %>%
  ggplot(aes(x = accuracy, y = tau, 
             color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.l2$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode"))
                                      # "NETCROP(bd)" = "NETCROP(bd)",
                                      # "NETCROP(bd)-Mean" = "NETCROP(bd)-Mean",
                                      # "NETCROP(bd)-Mode" = "NETCROP(bd)-Mode",
                                      # "NETCROP(AUC)" = "NETCROP(AUC)",
                                      # "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
                                      # "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
    )) +
  scale_color_manual(
    values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
               "Group 2" = "red"  # Color for Group 2 (l2)
  )) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  ) +
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               color = c(
                                 rep("darkred", 3),
                                 rep("darkblue", 3)
                                 )
                               ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

color.l2

ggsave("plots/l2_color.png", color.l2, width = 8, height = 4.5, dpi = 600, 
       device = "png")

################################################################################
## only 0, oracle, auc
plot.auc <- plot.out %>% filter(
  tau %notin% tau[grepl("NMI|Hemming|bd|l2", tau)]
) %>%
  mutate(
    tau = factor(tau,
                 levels = c("0", "Oracle", 
                            "Davis-Kahan Estimator",
                            "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode"
                            )
  ))

plot.auc$group <- factor(
  if_else(
    plot.auc$tau %in% c("0", "Oracle", "Davis-Kahan Estimator"), 
    "Group 1", "Group 2"
    )
)

bw.auc <- plot.auc %>%
  ggplot(aes(x = accuracy, y = tau, 
             # color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.auc$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(AUC)" = "NETCROP(AUC)",
      "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
      "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
      )
    ) +
  # scale_color_manual(
  #   values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
  #              "Group 2" = "red",   # Color for Group 2 (l2)
  #              "Group 3" = "green", # Color for Group 3 (bd)
  #              "Group 4" = "purple") # Color for Group 4 (AUC)
  # ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
                   y = expression(tau)
    )+
                   theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               # color = c(
                               #   rep("brown", 3),
                               #   rep("orange", 3),
                               #   rep("purple", 3),   
                               #   rep("darkgreen", 3),
                               #   rep("darkred", 3),
                               #   rep("darkblue", 3)
                               # )
    ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

bw.auc                   

ggsave("plots/auc_bw.png", bw.auc, width = 8, height = 4.5, dpi = 600, 
       device = "png")



color.auc <- plot.auc %>%
  ggplot(aes(x = accuracy, y = tau, 
             color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbarh(aes(xmin = xxmin, xmax = xxmax)) +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.auc$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(AUC)" = "NETCROP(AUC)",
      "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
      "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
    )
  ) +
  scale_color_manual(
    values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
               "Group 2" = "purple") # Color for Group 4 (AUC)
  ) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
  )+
  theme_pubclean() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "Times New Roman", size = 13,
                               face = "bold",
                               color = c(
                                 rep("purple", 3),
                                 rep("darkblue", 3)
                               )
    ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

color.auc                   

ggsave("plots/auc_color.png", color.auc, width = 8, height = 4.5, dpi = 600, 
       device = "png")




f <- function(n){
  n/(456-n)
}

curve(f, from = 1, to = 456)

g <- function(x){
  ceiling(456*x/(x+1))
}

g(2)
f(201)




















