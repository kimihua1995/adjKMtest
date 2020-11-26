\name{el.est}
\alias{el.est}
\title{Empirical Likelihood-Based Estimates}
\description{Estimate survival probabilities by empirical likelihood-based (EL) method}
\usage{
el.est(y,delta,treat,x,psix_moment = c("first","second"),
treat.select,t,get.sd=TRUE,Nboot=500,standardize=FALSE)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix_moment}{\code{first} or \code{second}, which correspond to adjusting first or second moments of the covariates matrix \code{x}}
\item{treat.select}{a value of 1 (treatment) or 0 (control)}
\item{t}{a numeric vector of time points for estimation}
\item{get.sd}{a logical object of whether return bootstrap standard deviation, default by \code{TRUE}}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}}
\item{standardize}{a logical object of whether the covariance matrix \code{x} would be centered and scaled, default by \code{FALSE}}
}

\value{
A numeric vector with the same length of \code{treat}, which contains the propensity scores for each observations.
}
\author{}
\examples{
library(survival)
y <- as.matrix(gbsg$rfstime/(365.25/12))
treat <- data.matrix(gbsg$hormon)
delta <- as.matrix(gbsg$status)
x <- data.matrix(gbsg[,c("age","meno","size","grade","nodes","pgr","er")])
t <- seq(0,80,2)

el.est(y=y,delta=delta,treat=treat,x=x,psix_moment="first",treat.select=1,t=t,
                   get.sd=TRUE,Nboot=500,standardize=FALSE)
}
\keyword{}
