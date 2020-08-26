\name{aPAC}
\alias{aPAC}
\title{
Adapted PAC scores
}
\description{
Adapted PAC scores
}
\usage{
aPAC(consensus_mat)
}
\arguments{

  \item{consensus_mat}{A consensus matrix.}

}
\details{
For the consensus values x, it is transformed to 1 - x if x < 0.5.
After the transformation, for any pair of samples in the consensus matrix,
If they are always in a same group or always in different groups, the 
value x is both to 1. Thus, if the consensus matrix shows stable partitions,
values x will be all close to 1. Reflected in the CDF of x, the curve is 
shifted to the right and the area under CDF curve should be very small.

An aPAC value less than 0.05 is considered as the stable partition, which can
be thought the proportion of abmiguous partitioning is less than 0.05.
}
\value{
A numeric value.
}
\examples{
data(golub_cola)
aPAC(get_consensus(golub_cola[1, 1], k = 2))
aPAC(get_consensus(golub_cola[1, 1], k = 3))
aPAC(get_consensus(golub_cola[1, 1], k = 4))
aPAC(get_consensus(golub_cola[1, 1], k = 5))
aPAC(get_consensus(golub_cola[1, 1], k = 6))
}
