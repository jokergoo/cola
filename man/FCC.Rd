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

  \item{consensus_mat}{A consensus matrix.}
  \item{diff}{Difference of F(b) - F(a).}

}
\details{
For a in [0, 0.5] and for b in [0.5, 1], the flatness measures
the flatness of the CDF curve of the consensus matrix. It is 
calculated as the maximum width that fits F(b) - F(a) <= diff
}
\value{
A numeric value.
}
\examples{
data(golub_cola)
FCC(get_consensus(golub_cola[1, 1], k = 2))
FCC(get_consensus(golub_cola[1, 1], k = 3))
FCC(get_consensus(golub_cola[1, 1], k = 4))
FCC(get_consensus(golub_cola[1, 1], k = 5))
FCC(get_consensus(golub_cola[1, 1], k = 6))
}
