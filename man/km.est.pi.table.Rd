\name{km.est.pi.table}
\alias{km.est.pi.table}
\title{Table for Inverse Probabilty of Treatment Weighting}
\description{Generate table for estimated survival probability by IPTW for treatment, control and their difference. Require \code{formattable} package to run.}
\usage{
km.est.pi.table(y,delta,treat,ex,t)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of the treatment indicator}
\item{ex}{a numeric vector of propensity score}
\item{t}{a numeric vector of time points for estimation}
}

\value{
A table (\code{formattable} object) containing the estimated S(t) with standard deviation at each time point.
}
\author{}
\examples{
}
