interpret.ptheta <- function(th, k) {
    ii <- 1:k
    ij <- k + 1:(k * (k-1)/2)
    L.prec <- diag(exp(th[ii]))
    L.prec[lower.tri(L.prec)] <- th[ij]
    ret <- list(mprec = tcrossprod(L.prec))
    ret$covar <- chol2inv(t(L.prec))
    ret$sd <- sqrt(diag(ret$covar))
    ret$correl <- t(ret$covar/ret$sd)/ret$sd
    return(ret)
}

interpret.vtheta <- function(th, k) {
    ii <- 1:k
    ij <- k + 1:(k * (k-1)/2)
    Lcov <- diag(exp(-th[ii]))
    Lcov[lower.tri(Lcov)] <- th[ij]
    ret <- list(mprec = chol2inv(t(Lcov)))
    ret$covar <- tcrossprod(Lcov)
    ret$sd <- sqrt(diag(ret$covar))
    ret$correl <- t(ret$covar/ret$sd)/ret$sd
    return(ret)
}

rgeneric.MBesag <- 
    function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                     "log.prior", "quit"), theta = NULL)
{
    envir <- parent.env(environment())
    if (exists("Args", envir=envir)) {
        Args <- get("Args", envir = envir)
    } else {
        assign("Args", Args, envir = envir)
    }
    K <- Args$NComps
    
    interpret.theta <- function(th, k) {
        ii <- 1:k
        ij <- k + 1:(k * (k-1)/2)
        L.prec <- diag(exp(th[ii]))
        L.prec[lower.tri(L.prec)] <- th[ij]
        ret <- list(mprec = tcrossprod(L.prec))
        ret$covar <- chol2inv(t(L.prec))
        ret$sd <- sqrt(diag(ret$covar))
        ret$correl <- t(ret$covar/ret$sd)/ret$sd
        return(ret)
    }

    graph <- function() {
        return(Q()) 
    }
    
    Q <- function() {        
        mprec <- interpret.theta(theta, K)$mprec 
        if(!is.null(Args$Structure))
            mprec <- kronecker(mprec, Args$Structure)
        return(mprec)
    }
    
    mu <- function() {
        return(numeric(0))
    }
    
    log.norm.const <- function() {
        return(numeric(0))
    }
    
    log.prior <- function() {
        params <- interpret.theta(theta, K)
        stopifnot(Args$LKJeta>0)
        eta <- Args$LKJeta
        k <- 1:(K-1)
        a <- sum((2*(eta -1) + K-k)*(K-k)) * log(2)
        b <- sum(lbeta(a = eta + (K-k-1)/2,
                       b = eta + (K-k-1)/2) * (K-k))
        cc <- (eta-1)*sum(log(diag(chol(params$correl)))*2)

        stopifnot(all(Args$pcprec.lambda>0))
        lambs <- Args$pcprec.lambda 
        precs <- diag(params$mprec)

        val <- a + b + cc +
          sum(dgamma(log(precs), shape = 1, rate = 1, log = TRUE))
        
        return(val)
    }
    
    initial <- function() {
        ret <- c(rep(K/2, K),
                 rep(1, K * (K -1)/2))
        return(ret)
    }
    
    quit <- function() {
        return(invisible())
    }
    
    if (!length(theta)) {
        theta <- initial()
    }
    
    val <- do.call(match.arg(cmd), args = list())
    return(val)
}

rgeneric.CorGraphsBesag <- 
    function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                     "log.prior", "quit"), theta = NULL)
{
  envir <- parent.env(environment())
  library(Rgraphviz)
  library(numDeriv)
  if (exists("GS", envir=envir)) {
    GS <- get("GS", envir=envir) # need to get or already in environment?
  } else {
    GS <- Argm$SP # graph structure
    assign("GS", GS, envir=envir)
  }
  

  graph <- function() {
    return(Q()) # does not depend on theta
  }
  
  Q <- function() {
      jprec <- 1:GS$NC
      q <- exp(theta[-jprec])
      ## numerical precision matrix
      Q <- hessian(GS$JD[[1]],
                   rep(0, length(GS$STR[[1]])),
                   SDev = q) * GS$Pmat 
      diag(Q) <- diag(Q) + 1e-12 ## safe
      VC.all <- solve(Q) ## numerical variance-covariance
      VC <- VC.all[jprec, jprec]
      num.isd <- 1/sqrt(diag(VC))
      mcor <- t(VC * num.isd) * num.isd ## implied correlation
      sd.m <- exp(-0.5 * theta[jprec])
      COV <- t(mcor * sd.m) * sd.m ## model variance-covariance
      Q.k <- solve(COV) ## model precision
      if(!is.null(GS$Structure)) {
          Q.k <- kronecker(Q.k, GS$Structure)
      }
      return(Q.k)
  }
  
  mu <- function() {
    return(numeric(0))
  }

  log.norm.const <- function() {
    return(numeric(0))
  }

  log.prior <- function() {
    q <- theta[-(1:GS$NC)]
    val <- (sum(dgamma(exp(theta[1:GS$NC]), shape = 1, rate = 1, log = TRUE)) +
            sum(theta[1:GS$NC]) +
            sum(unlist(Argm$GraphPrior(Argm$S, q, Argm$lambda, GS, Argm$Tdist)))+
            sum(theta[-(1:GS$NC)]))
    return(val)
  }

    initial <- function() {
        ret <- c(rep(4, GS$NC), rep(Argm$init, GS$NP))
        return(ret)
  }

  quit <- function() {
    return(invisible())
  }

  if (!length(theta)) {
    theta <- initial()
  }

  val <- do.call(match.arg(cmd), args = list())
  return(val)
}

