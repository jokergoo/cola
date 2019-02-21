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

  \item{membership_each}{a matrix which contains partitions in every single runs.}
  \item{class}{consensus class IDs.}

}
\details{
Class IDs in \code{membership_meach} have already be adjusted to the consensus class IDs
to let \code{sum(x_single == x_consensus)} reach maximum.

The concordance score is the mean probability of fitting the consensus class IDs in all
partitions.

This function is used internally.
}
\value{
A numeric value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
membership_each = get_membership(cola_rl["sd", "kmeans"], each = TRUE, k = 3)
consensus_classes = get_classes(cola_rl["sd", "kmeans"], k = 3)$class
concordance(membership_each, consensus_classes)
}
