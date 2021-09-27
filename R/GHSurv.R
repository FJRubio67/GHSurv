#--------------------------------------------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) hazard function.
#' http://rpubs.com/FJRubio/PGW
#--------------------------------------------------------------------------------------------------------------------------
#' @param eta     : scale parameter
#' @param nu      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the PGW hazard function
#' @export

hpgw <- function(t, eta, nu, delta, log = FALSE){
  val <- log(nu) - log(delta) - nu*log(eta) + (nu-1)*log(t) +
    (1/delta - 1)*log( 1 + (t/eta)^nu )
  if(log) return(val) else return(exp(val))
}

#--------------------------------------------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) cumulative hazard function.
#' http://rpubs.com/FJRubio/PGW
#--------------------------------------------------------------------------------------------------------------------------
#' @param eta     : scale parameter
#' @param nu      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @return the value of the PGW cumulative hazard function
#' @export

chpgw <- function(t, eta, nu, delta){
  val <- -1 + ( 1 + (t/eta)^nu )^(1/delta)
  return(val)
}

#--------------------------------------------------------------------------------------------------------------------------
#' Exponentiated Weibull (EW) hazard function.
#' http://rpubs.com/FJRubio/EW
#--------------------------------------------------------------------------------------------------------------------------
#' @param lambda     : scale parameter
#' @param kappa      : shape parameter
#' @param alpha   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the EW hazard function
#' @export

hexpweibull <- function(t,lambda,kappa,alpha,log=FALSE){
  log.pdf <-  log(alpha) + (alpha-1)*pweibull(t,scale=lambda,shape=kappa,log.p=TRUE) +
    dweibull(t,scale=lambda,shape=kappa,log=TRUE)
  cdf <- exp(alpha*pweibull(t,scale=lambda,shape=kappa,log.p=TRUE) )
  log.h <- log.pdf - log(1-cdf)
  ifelse(log, return(log.h), return(exp(log.h)))
}

#--------------------------------------------------------------------------------------------------------------------------
#' Exponentiated Weibull (EW) cumulative hazard function.
#' http://rpubs.com/FJRubio/EW
#--------------------------------------------------------------------------------------------------------------------------
#' @param lambda     : scale parameter
#' @param kappa      : shape parameter
#' @param alpha   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the EW cumulative hazard function
#' @export

chexpweibull <- function(t,lambda,kappa,alpha,log.p=FALSE){
  cdf <- exp(alpha*pweibull(t,scale=lambda,shape=kappa,log.p=TRUE) )
  return(-log(1-cdf))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) Probability Density Function.
#' http://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param theta   : scale parameter
#' @param kappa      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the GG probability density function
#' @export

dggamma <- function(t, theta, kappa, delta, log = FALSE){
  val <- log(delta) - kappa*log(theta) - lgamma(kappa/delta) + (kappa - 1)*log(t) -
    (t/theta)^delta
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) Cumulative Distribution Function.
#' http://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param theta   : scale parameter
#' @param kappa      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' #' @param log.p: log scale (TRUE or FALSE)
#' @return the value of the GG cumulative distribution function
#' @export

pggamma <- function(t, theta, kappa, delta, log.p = FALSE){
  val <- pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE)
  if(log.p) return(val) else return(exp(val))
}


#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) Survival Function.
#' http://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param theta   : scale parameter
#' @param kappa      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the GG survival function
#' @export
sggamma <- function(t, theta, kappa, delta, log.p = FALSE){
  val <- pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE, lower.tail =  FALSE)
  if(log.p) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) Hazard Function.
#' http://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param theta   : scale parameter
#' @param kappa      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the GG hazard function
#' @export
hggamma <- function(t, theta, kappa, delta, log = FALSE){
  val <- dggamma(t, theta, kappa, delta, log = TRUE) - sggamma(t, theta, kappa, delta, log.p = TRUE)
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Generalised Gamma (GG) Cumulative Hazard Function.
#' http://rpubs.com/FJRubio/GG
#----------------------------------------------------------------------------------------
#' @param theta   : scale parameter
#' @param kappa      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @return the value of the GG cumulative hazard function
#' @export
chggamma <- function(t, theta, kappa, delta){
  val <- -pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE, lower.tail =  FALSE)
  return(val)
}


#----------------------------------------------------------------------------------------
#' Lognormal (LN) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LN hazard function
#' @export

hlnorm <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlnorm(t,mu,sigma, log = T)
  ls0 <- plnorm(t,mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Lognormal (LN) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LN cumulative hazard function
#' @export
chlnorm <- function(t,mu,sigma){
  H0 <- -plnorm(t,mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Log-logistic (LL) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LL hazard function
#' @export
hllogis <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlogis(log(t),mu,sigma, log = T) - log(t)
  ls0 <- plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Lognormal (LL) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LL cumulative hazard function
#' @export
chllogis <- function(t,mu,sigma){
  H0 <- -plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Gamma (G) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param shape   : shape parameter
#' @param scale      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the G hazard function
#' @export
hgamma <- function(t, shape, scale, log = FALSE){
  lpdf0 <-  dgamma(t, shape = shape, scale = scale, log = T)
  ls0 <- pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Gamma (G) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param shape   : shape parameter
#' @param scale      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the G cumulative hazard function
#' @export
chgamma <- function(t, shape, scale){
  H0 <- -pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

###############################################################################################
###############################################################################################
###############################################################################################
#' Overall Survival GH model.
###############################################################################################
###############################################################################################
###############################################################################################

########################################################################################################
#' Log likelihood and MLE for the GH hazards model.
#' Baseline hazards: Lognormal, Log-logistic, Gamma, PGW, EW, GG
########################################################################################################
#' @param init  : initial points for optimisation
#' @param des_t : design matrix for time-dependent effects (q x n), q >= 1
#' @param des  : design matrix for hazard-level effects (p x n), p >= 1
#' @param status : vital status (1 - dead, 0 - alive)
#' @param times  : survival times
#' @param hstr : hazard structure including baseline (LNGH,LLGH,GGH,PGWGH,EWGH,GGGH)
#' @param method: "nlminb" or a method from "optim"
#' @param maxit: The maximum number of iterations. Defaults to 1000
#' @return a list containing the output of the optimisation (OPT) and the log-likelihood function (loglik)
#' @export

GHMLE <- function(init, times, status, hstr, des_t, des, method = "Nelder-Mead", maxit = 1000){
# Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  times.obs <- times[status]
  if(!is.null(des_t))  des.obs <- des[status,]
  if(!is.null(des_t))  des_t.obs <- des_t[status,]
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]

# PGW - GH Model
  if(hstr == "PGWGH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

# EW - GH Model
  if(hstr == "EWGH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hexpweibull(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chexpweibull(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

# GG - GH Model
  if(hstr == "GGGH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hggamma(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chggamma(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

# Gamma - GH Model
  if(hstr == "GGH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

# Lognormal - GH Model
  if(hstr == "LNGH"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

# Loglogistic - GH Model
  if(hstr == "LLGH"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs
      val <- - sum(lhaz0) + sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, loglik = log.lik)
  return(OUT)
}


###############################################################################################
###############################################################################################
###############################################################################################
#' Relative Survival GH model.
###############################################################################################
###############################################################################################
###############################################################################################

########################################################################################################
#' Log likelihood and MLE for the GH excess hazards model.
#' Baseline hazards: Lognormal, Log-logistic, Gamma, PGW, EW, GG
########################################################################################################
#' @param init  : initial points for optimisation
#' @param des_t  : design matrix for time-dependent effects (q x n), q >= 1
#' @param des  : design matrix for hazard-level effects (p x n), p >= 1
#' @param status : vital status (1 - dead, 0 - alive)
#' @param times  : survival times
#' @param hstr : hazard structure including baseline (LNGEH,LLGEH,GGEH,PGWGEH,EWGEH,GGGEH)
#' @param hp.obs  : population hazards (for uncensored individuals)
#' @param method: "nlminb" or a method from "optim"
#' @param maxit: The maximum number of iterations. Defaults to 1000
#' @return a list containing the output of the optimisation (OPT) and the log-likelihood function (loglik)
#' @export

GEHMLE <- function(init, times, status, hstr, des_t, des, hp.obs, method = "Nelder-Mead", maxit = 1000){
  # Required variables
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  times.obs <- times[status]
  if(!is.null(des_t))  des.obs <- des[status,]
  if(!is.null(des_t))  des_t.obs <- des_t[status,]
  hp.obs <- as.vector(hp.obs)
  p0 <- dim(des_t)[2]
  p1 <- dim(des)[2]

  # PGW - GEH Model
  if(hstr == "PGWGEH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

  # EW - GEH Model
  if(hstr == "EWGEH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hexpweibull(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chexpweibull(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

  # GG - GEH Model
  if(hstr == "GGGEH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);  ce0 <- exp(par[3]); alpha <- par[4:(3+p0)]; beta <- par[(4+p0):(3+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hggamma(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chggamma(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

  # Gamma - GEH Model
  if(hstr == "GGEH"){
    log.lik <- function(par){
      ae0 <- exp(par[1]); be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

  # Lognormal - GEH Model
  if(hstr == "LNGEH"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }

  # Log-logistic - GEH Model
  if(hstr == "LLGEH"){
    log.lik <- function(par){
      ae0 <- par[1]; be0 <- exp(par[2]);   alpha <- par[3:(2+p0)]; beta <- par[(3+p0):(2+p0+p1)]
      x.alpha <- as.vector(des_t%*%alpha)
      x.beta <- as.vector(des%*%beta)
      exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
      exp.x.alpha <- exp(x.alpha)
      exp.x.alpha.obs <- exp(x.alpha[status])
      x.beta.obs <- x.beta[status]
      lhaz0 <- log( hp.obs + hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = FALSE)*exp(x.beta.obs) )
      val <- - sum(lhaz0) + sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
      return(sum(val))
    }
  }
  if(method != "nlminb") OPT <- optim(init,log.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,log.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, loglik = log.lik)
  return(OUT)
}


###############################################################################################
###############################################################################################
###############################################################################################
#' Confidence intervals.
###############################################################################################
###############################################################################################
###############################################################################################

###########################################################################################
#' Function to calculate the normal confidence intervals.
#' The parameters indicated with "index" are transformed to the real line using log().
###########################################################################################
#' @param FUN   : minus log-likelihood function to be used to calculate the confidence intervals
#' @param MLE   : maximum likelihood estimator of the parameters of interest
#' @param level : confidence level
#' @param index : position of the positive parameters under the original parameterisation
#' @return a list containing the upper and lower conf.int limits, the transformed MLE, and std errors
#' @export

Conf.Int <- function(FUN,MLE,level=0.95,index=NULL){
  sd.int <- abs(qnorm(0.5*(1-level)))
  tempf <- function(par){
    par[index] = exp(par[index])
    return(FUN( par ))
  }
  r.MLE <- MLE
  r.MLE[index] <- log(MLE[index])
  HESS <- hessian(tempf,x=r.MLE)
  Fisher.Info <- solve(HESS)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- r.MLE + sd.int*Sigma
  L<- r.MLE - sd.int*Sigma
  C.I <- cbind(L,U,r.MLE, Sigma)
  names.row <- paste0("par", seq_along(1:length(MLE)))
  names.row[index] <- paste0("log.par", seq_along(index))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

