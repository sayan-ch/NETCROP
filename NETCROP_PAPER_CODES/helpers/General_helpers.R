library(Matrix)
library(parallel)
################################################################################
## General helpers: usually used in all functions=
modal <- function(x) { # computes modal value in a vector
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

l2 <- function(x, y) { #squared l2 distance 
  sum((x - y)^2)
}

bin.dev <- function(x, y) { # binomial deviance loss: y: prob, x: binary response
  y[y < 1e-5] <- 1e-5
  y[y > 1 - 1e-5] <- 1 - 1e-5
  
  tmp <- -x * log(y) - (1 - x) * log(1 - y)
  
  return(sum(tmp, rm.na = T))
}

AUC <- function(AA, PP) { # returns -ve so that it can be used as a loss function to be minimized
  -Rfast::auc(as(AA, "numeric"), as(PP, "numeric"))
}

softplus <- function(x) pmax(x,0) + log1p(exp(-abs(x)))

sigmoid  <- function(x) 1/(1+exp(-x))
################################################################################
# NETCROP parameter selection
netcrop_param <- function(p.test = 0.02, n = NULL,
                          o.range = c(0, 0.8)){
  over.upper <- (1 - sqrt(p.test))
  over.lower <- 1 - sqrt(2*p.test)
  
  o.range[o.range == 1] <- 0.8
  
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

