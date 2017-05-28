\name{PAC}
\alias{PAC}
\title{
PAC scores
}
\description{
PAC scores
}
\usage{
PAC(consensus_mat)
}
\arguments{

  \item{consensus_mat}{consensus matrix}

}
\details{
This a variant of the orignial PAC (proportion of ambiguous clustering) method.

Assume x_1 in [0, 0.3] and x_2 in [0.7, 1], we calculate s = (\\int_{x_1}^{x_2} F(x) - F(x1)*(x_2 - x_1))/(x_2 - x_1),
and the mean value is taken as the final value.
}
\value{
A single PAC scores.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
