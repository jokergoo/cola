\name{guess_best_k-ConsensusPartition-method}
\alias{guess_best_k,ConsensusPartition-method}
\title{
Guess the best number of partitions
}
\description{
Guess the best number of partitions
}
\usage{
\S4method{guess_best_k}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}

}
\details{
It looks for the best k with highest cophenetic correlation coefficient
or lowest PAC score or highest mean silhouette value.
}
\value{
The best k
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
rl = readRDS(system.file("extdata/example.rds", package = "cola"))
obj = rl["sd", "kmeans"]
guess_best_k(obj)
}
