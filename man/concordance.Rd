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

  \item{membership_each}{A matrix which contains partitions in every single runs where columns correspond to runs. The object can be get from \code{get_membership(..., each = TRUE)}.}
  \item{class}{Consensus subgroup labels.}

}
\details{
Note subgroup labels in \code{membership_each} should already be adjusted to the consensus labels, i.e. by \code{\link{relabel_class}}.

The concordance score is the mean proportion of samples having the same subgroup labels as the consensus labels
among individual partition runs.
}
\value{
A numeric value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
membership_each = get_membership(golub_cola["SD", "kmeans"], each = TRUE, k = 3)
consensus_classes = get_classes(golub_cola["SD", "kmeans"], k = 3)$class
concordance(membership_each, consensus_classes)
}
