\name{km.est}
\alias{km.est}
\title{Kaplan-Meier Estimates}
\description{Estimate survival probabilities by standard KM method}
\usage{
km.est(y,delta,treat,treat.select,t)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of indicator for treatment}
\item{treat.select}{a numeric value of treatment selected for estimation, 1 for treatment group, 0 for control group}
\item{t}{a numeric vector of time points for estimation}
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
t <- seq(0,80,2)

km.est(y=y,delta=delta,treat=treat,treat.select=1,t=t)
}
\keyword{}
