\name{ps}
\alias{ps}
\title{Propensity Scores}
\description{Calculate the propensity scores.}
\usage{
ps(treat,x)
}
\arguments{
\item{treat}{a factor vector of the treatment indicator}
\item{x}{an optional covariace matrix containing the variables in the model}
}

\value{
A numeric vector with the same length of \code{treat}, which contains the propensity scores for each observations.
}
\author{}
\examples{
}
\keyword{}
