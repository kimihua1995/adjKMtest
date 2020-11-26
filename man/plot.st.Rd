\name{plot.st}
\alias{plot.st}
\title{Plots for Survival Curves}
\description{Plots of survival curves for treatment, control and their difference by 3 methods (KM, IPTW, EL).}
\usage{
plot.st(y,delta,treat,x,psix_moment=c("first","second"),t,
method = c("KM","IPTW","EL"),Nboot=500,
standardize=FALSE,alpha=0.05)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix_moment}{\code{first} or \code{second}, which correspond to adjusting first or second moments of the covariates matrix \code{x}. Only applied when \code{method="EL"}}
\item{t}{a numeric vector of time points for estimation}
\item{method}{the estimation method selected from KM, IPTW or EL}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}. Only applied when \code{method="EL"}}
\item{standardize}{a logical object of whether the covariance matrix \code{x} would be centered and scaled, default by \code{FALSE}}
\item{alpha}{alpha level for the conficendence interval, default by \code{0.05}}
}

\value{
A plot object with survival curves for treatment, control and their difference with pointwise 95\% confidence intervals. The CIs for KM and IPTW are asymptotic estimates, for EL is estimated by the quantiles of Bootstrap
}
\author{}
\examples{
library(survival)
y <- as.matrix(gbsg$rfstime/(365.25/12))
treat <- data.matrix(gbsg$hormon)
delta <- as.matrix(gbsg$status)
x <- data.matrix(gbsg[,c("age","meno","size","grade","nodes","pgr","er")])
t <- seq(0,80,2)

plot.st(y=y,delta=delta,treat=treat,x=x,t=t,method="KM",alpha=0.05)
plot.st(y=y,delta=delta,treat=treat,x=x,t=t,method="IPTW",alpha=0.05)
plot.st(y=y,delta=delta,treat=treat,x=x,psix_moment="second",t=t,
method="EL",Nboot=500,standardize=TRUE,alpha=0.05)

}
\keyword{}
