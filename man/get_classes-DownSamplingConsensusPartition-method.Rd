\name{get_classes-DownSamplingConsensusPartition-method}
\alias{get_classes,DownSamplingConsensusPartition-method}
\title{
Get subgroup labels
}
\description{
Get subgroup labels
}
\usage{
\S4method{get_classes}{DownSamplingConsensusPartition}(object, k = object@k, p_cutoff = 0.05)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{p_cutoff}{Cutoff of p-values of class label prediction. It is only used when \code{k} is a vector.}

}
\value{
If \code{k} is a scalar, it returns a data frame with two columns:

\itemize{
  \item the class labels
  \item the p-value for the prediction of class labels.
}

If \code{k} is a vector, it returns a data frame of class labels for each k. The class
label with prediction p-value > \code{p_cutoff} is set to \code{NA}.
}
\examples{
data(golub_cola_ds)
get_classes(golub_cola_ds, k = 3)
get_classes(golub_cola_ds)
}
