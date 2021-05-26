\name{get_signatures-DownSamplingConsensusPartition-method}
\alias{get_signatures,DownSamplingConsensusPartition-method}
\title{
Get signature rows
}
\description{
Get signature rows
}
\usage{
\S4method{get_signatures}{DownSamplingConsensusPartition}(object, k,
    p_cutoff = 1, ...)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{p_cutoff}{Cutoff for p-values of class label prediction. Samples with values  higher than it are not used for finding signature rows.}
  \item{...}{Other arguments passed to \code{\link{get_signatures,ConsensusPartition-method}}.}

}
\details{
This function is very similar as \code{\link{get_signatures,ConsensusPartition-method}}.
}
\examples{
\donttest{
data(golub_cola_ds)
get_signatures(golub_cola_ds, k = 2)
get_signatures(golub_cola_ds, k = 3)
}
}
