% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GHSurv.R
\name{Conf.Int}
\alias{Conf.Int}
\title{Confidence intervals.
Function to calculate the normal confidence intervals.
The parameters indicated with "index" are transformed to the real line using log().}
\usage{
Conf.Int(FUN, MLE, level = 0.95, index = NULL)
}
\arguments{
\item{FUN}{: minus log-likelihood function to be used to calculate the confidence intervals}

\item{MLE}{: maximum likelihood estimator of the parameters of interest}

\item{level}{: confidence level}

\item{index}{: position of the positive parameters under the original parameterisation}
}
\value{
a list containing the upper and lower conf.int limits, the transformed MLE, and std errors
}
\description{
Confidence intervals.
Function to calculate the normal confidence intervals.
The parameters indicated with "index" are transformed to the real line using log().
}
