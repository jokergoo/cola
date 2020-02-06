\name{get_classes-ConsensusPartitionList-method}
\alias{get_classes,ConsensusPartitionList-method}
\title{
Get subgroup labels
}
\description{
Get subgroup labels
}
\usage{
\S4method{get_classes}{ConsensusPartitionList}(object, k)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{Number of subgroups.}

}
\details{
The subgroup labels are inferred by merging partitions from all methods
by weighting the mean silhouette scores in each method.
}
\value{
A data frame with subgroup labels and other columns which are entropy of the percent membership matrix
and the silhouette scores which measure the stability of a sample to stay in its group.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
get_classes(cola_rl, k = 2)
}
