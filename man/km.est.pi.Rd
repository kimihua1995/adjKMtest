\name{km.est.pi}
\alias{km.est.pi}
\title{Inverse Probabilty of Treatment Weighting}
\description{Estimate survival probability by IPTW}
\usage{
km.est.pi(y,delta,treat,ex,treat.select,t)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{ex}{a numeric vector of propensity score}
\item{treat.select}{a value of 1 (treatment) or 0 (control)}
\item{t}{a numeric vector of time points for estimation}
}

\value{
A list containing the following components:
\item{St}{point estimates of survival probability (S(t)) at \code{t}}
\item{sd}{standard deviation of survival probability (S(t)) at \code{t}}
}
\author{}
\examples{
}
