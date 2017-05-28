\name{get_param-ConsensusPartition-method}
\alias{get_param,ConsensusPartition-method}
\alias{get_param}
\title{
Get parameters
}
\description{
Get parameters
}
\usage{
\S4method{get_param}{ConsensusPartition}(object, k, unique = TRUE)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}
  \item{unique}{whether apply \code{\link[base]{unique}} to rows}

}
\value{
A data frame of parameters corresponding to the current k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
