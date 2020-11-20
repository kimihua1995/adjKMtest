\name{el.est.table}
\alias{el.est.table}
\title{Table for Empirical Likelihood-Based (EL) Method}
\description{Generate table for estimated survival probabilities by EL for treatment, control and their difference. Require \code{formattable} package to run.}
\usage{
el.est.table(y,delta,treat,x,psix,t,Nboot=500)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix}{a function of covariance matrix \code{x}}
\item{t}{a numeric vector of time points for estimation}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}}
}

\value{
A table (\code{formattable} object) containing the estimated S(t) with standard deviation at each time point.
}
\author{}
\examples{
}
\keyword{}
