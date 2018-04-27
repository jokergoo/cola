\name{concordance}
\alias{concordance}
\title{
Concordance of partitions to the consensus partition
}
\description{
Concordance of partitions to the consensus partition
}
\usage{
concordance(membership_each, class)
}
\arguments{

  \item{membership_each}{all repetetive parititions.}
  \item{class}{consensus class ids.}

}
\details{
Class ids in \code{membership_meach} have already be adjusted to the consensus class ids
to let \code{sum(x1 == x_consensus)} to get maximum.

The concordance score is the mean probability of fiting the consensus class ids in all
partitions.

This function is used internally.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
