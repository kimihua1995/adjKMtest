\name{km.est.table}
\alias{km.est.table}
\title{Table for Kaplan-Meier Estimates}
\description{Generate table from estimated survival probabilities by standard KM method for treatment, control and their difference. Require \code{formattable} package to run.}
\usage{
km.est.table(y,delta,treat,t)
}
\arguments{
\item{y}{a numeric vector of observed time}
\item{delta}{a factor vector of indicator for events}
\item{treat}{a factor vector of indicator for treatment, 1 for treatment group, 0 for control group}
\item{t}{a numeric vector of time points for estimation}
}

\value{
A table (\code{formattable} object) containing the estimated S(t) with standard deviation at each time point.
}
\author{}
\examples{
}
\keyword{}
