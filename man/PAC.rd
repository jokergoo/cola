\name{PAC}
\alias{PAC}
\title{
PAC score
}
\description{
PAC score
}
\usage{
PAC(consensus_mat, original = FALSE)
}
\arguments{

  \item{consensus_mat}{a consensus matrix}

}
\details{
This a variant of the orignial PAC (proportion of ambiguous clustering) method.

Assume x_1 in [0, 0.3] and x_2 in [0.7, 1], we calculate s = (\\int_{x_1}^{x_2} F(x) - F(x1)*(x_2 - x_1))/(x_2 - x_1) where F(x) is the CDF of the consensus matrix
and the mean value is taken as the final value.
}
\value{
A single PAC score.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
