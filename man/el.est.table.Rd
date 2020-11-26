\name{el.est.table}
\alias{el.est.table}
\title{Table for Empirical Likelihood-Based (EL) Method}
\description{Generate table for estimated survival probabilities by EL for treatment, control and their difference. Require \code{formattable} package to run.}
\usage{
el.est.table(y,delta,treat,x,psix_moment = c("first","second"),
t,Nboot=500,standardize=FALSE)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix_moment}{\code{first} or \code{second}, which correspond to adjusting first or second moments of the covariates matrix \code{x}}
\item{t}{a numeric vector of time points for estimation}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}}
\item{standardize}{a logical object of whether the covariance matrix \code{x} would be centered and scaled, default by \code{FALSE}}
}

\value{
A table (\code{formattable} object) containing the estimated S(t) with standard deviation at each time point.
}
\author{}
\examples{
library(survival)
y <- as.matrix(gbsg$rfstime/(365.25/12))
treat <- data.matrix(gbsg$hormon)
delta <- as.matrix(gbsg$status)
x <- data.matrix(gbsg[,c("age","meno","size","grade","nodes","pgr","er")])
t <- seq(0,80,2)

el.est.table(y=y,delta=delta,treat=treat,x=x,psix_moment="first",t=t,
             Nboot=500,standardize=FALSE)
}
\keyword{}
