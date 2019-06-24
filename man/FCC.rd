\name{FCC}
\alias{FCC}
\title{
Flatness of the CDF curve
}
\description{
Flatness of the CDF curve
}
\usage{
FCC(consensus_mat, diff = 0.1)
}
\arguments{

  \item{consensus_mat}{a consensus matrix.}
  \item{diff}{Difference of F(b) - F(a)}

}
\details{
For a in [0, 0.5] and for b in [0.5, 1], the flatness measures
the flatness of the CDF curve of the consensus matrix, it is 
calculated as the maximum width that fits F(b) - F(a) <= diff

A flatness larger than 0.9 is treated as stable partitions.
}
\examples{
# There is no example
NULL

}
