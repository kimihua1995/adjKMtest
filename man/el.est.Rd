\name{el.est}
\alias{el.est}
\title{Empirical Likelihood-Based Estimates}
\description{Estimate survival probabilities by empirical likelihood-based (EL) method}
\usage{
el.est(y,delta,treat,x,psix,treat.select,t,get.sd=TRUE,Nboot=500)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix}{a function of covariance matrix \code{x}}
\item{treat.select}{a value of 1 (treatment) or 0 (control)}
\item{t}{a numeric vector of time points for estimation}
\item{get.sd}{a logical object of whether return bootstrap standard deviation, default by \code{TRUE}}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}}
}

\value{
A numeric vector with the same length of \code{treat}, which contains the propensity scores for each observations.
}
\author{}
\examples{
}
\keyword{}
