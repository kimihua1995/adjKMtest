\name{km.est.pi}
\alias{km.est.pi}
\title{Inverse Probabilty of Treatment Weighting}
\description{Estimate survival probabilities by IPTW}
\usage{
km.est.pi(y,delta,treat,x,treat.select,t,standardize=FALSE)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{treat.select}{a value of 1 (treatment) or 0 (control)}
\item{t}{a numeric vector of time points for estimation}
\item{standardize}{a logical object of whether the covariance matrix \code{x} would be centered and scaled, default by \code{FALSE}}
}

\value{
A list containing the following components:
\item{St}{point estimates of survival probability (S(t)) at \code{t}}
\item{sd}{standard deviation of survival probability (S(t)) at \code{t}}
}
\author{}
\examples{
library(survival)
y <- as.matrix(gbsg$rfstime/(365.25/12))
treat <- data.matrix(gbsg$hormon)
delta <- as.matrix(gbsg$status)
x <- data.matrix(gbsg[,c("age","meno","size","grade","nodes","pgr","er")])
t <- seq(0,80,2)

km.est.pi(y=y,delta=delta,treat=treat,x=x,treat.select=1,t=t,standardize=FALSE)
}
