\name{km.est.pi.table}
\alias{km.est.pi.table}
\title{Table for Inverse Probabilty of Treatment Weighting}
\description{Generate table for estimated survival probabilities by IPTW for treatment, control and their difference. Require \code{formattable} package to run.}
\usage{
km.est.pi.table(y,delta,treat,x,t,standardize=FALSE)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{t}{a numeric vector of time points for estimation}
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

km.est.pi.table(y=y,delta=delta,treat=treat,x=x,t=t,standardize=FALSE)
}
