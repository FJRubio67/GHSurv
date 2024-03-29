% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHSurv.R
\name{GEHMLE}
\alias{GEHMLE}
\title{Relative Survival GH model.
Log likelihood and MLE for the GH excess hazards model.
Baseline hazards: Lognormal, Log-logistic, Gamma, PGW, EW, GG}
\usage{
GEHMLE(
  init,
  times,
  status,
  hstr,
  des_t,
  des,
  hp.obs,
  method = "Nelder-Mead",
  maxit = 1000
)
}
\arguments{
\item{init}{: initial points for optimisation}

\item{times}{: survival times}

\item{status}{: vital status (1 - dead, 0 - alive)}

\item{hstr}{: hazard structure including baseline (LNGEH,LLGEH,GGEH,PGWGEH,EWGEH,GGGEH)}

\item{des_t}{: design matrix for time-dependent effects (q x n), q >= 1}

\item{des}{: design matrix for hazard-level effects (p x n), p >= 1}

\item{hp.obs}{: population hazards (for uncensored individuals)}

\item{method:}{"nlminb" or a method from "optim"}

\item{maxit:}{The maximum number of iterations. Defaults to 1000}
}
\value{
a list containing the output of the optimisation (OPT) and the log-likelihood function (loglik)
}
\description{
Relative Survival GH model.
Log likelihood and MLE for the GH excess hazards model.
Baseline hazards: Lognormal, Log-logistic, Gamma, PGW, EW, GG
}
