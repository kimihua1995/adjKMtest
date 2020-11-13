\name{plot.st}
\alias{plot.st}
\title{Plots for Survival Curves}
\description{Plots of survival curves for treatment, control and their difference by 3 methods (KM, IPTW, EL).}
\usage{
plot.st(y,delta,treat,x,psix=x,t,method = c("KM","IPTW","EL"),Nboot=500)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariance matrix containing the variables in the model}
\item{psix}{a function of covariance matrix \code{x}}
\item{t}{a numeric vector of time points for estimation}
\item{method}{the estimation method selected from KM, IPTW or EL}
\item{Nboot}{number of iterations in Bootstrap, default by \code{500}}
}

\value{
A plot object with survival curves for treatment, control and their difference with pointwise 95\% confidence intervals. The CIs for KM and IPTW are asymptotic estimates, for EL is estimated by the quantiles of Bootstrap
}
\author{}
\examples{
}
\keyword{}
