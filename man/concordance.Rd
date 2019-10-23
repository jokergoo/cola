\name{concordance}
\alias{concordance}
\title{
Concordance to the consensus partition
}
\description{
Concordance to the consensus partition
}
\usage{
concordance(membership_each, class)
}
\arguments{

  \item{membership_each}{A matrix which contains partitions in every single runs where columns correspond to runs.}
  \item{class}{Consensus class IDs.}

}
\details{
Note class IDs in \code{membership_each} should already be adjusted to the consensus class IDs
to let \code{sum(x_single == x_consensus)} reach maximum.

The concordance score is the mean proportion of samples having the same class ID as the consensus class ID
among runs.

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
