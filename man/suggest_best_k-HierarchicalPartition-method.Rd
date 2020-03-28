\name{suggest_best_k-HierarchicalPartition-method}
\alias{suggest_best_k,HierarchicalPartition-method}
\title{
Guess the best number of partitions
}
\description{
Guess the best number of partitions
}
\usage{
\S4method{suggest_best_k}{HierarchicalPartition}(object, jaccard_index_cutoff = 0.95)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{jaccard_index_cutoff}{The cutoff for Jaccard index for comparing to previous k.}

}
\details{
It basically gives the best k at each node.
}
\value{
A data frame with the best k and other statistics for each node.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
suggest_best_k(cola_rh)
}
