\name{register_partition_fun}
\alias{register_partition_fun}
\title{
Register user-defined partition functions
}
\description{
Register user-defined partition functions
}
\usage{
register_partition_fun(..., scale_method = c("standardization", "rescale", "none"))
}
\arguments{

  \item{...}{a named list of functions}
  \item{scale_method}{normally, data matrix are scaled by rows before sent to the partition function. The default scaling is apply by \code{\link[base]{scale}}. However, some partition function may not accept negative values which  are produced by \code{\link[base]{scale}}. Here \code{scale_method} can be set to \code{rescale} which scales rows by \code{(x - min)/(max - min)}. Note here \code{scale_method} only means how to scale rows when \code{scale_rows} is set to \code{TRUE} in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}. The value for \code{scale_method} can be a vector.}

}
\details{
The user-defined function should only accept three arguments. The first two arguments are the data
matrix and the number of partitions. The third argument should always be \code{\link{...}} so that parameters
for the partition function can be passed by \code{partition_param} from \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.

The function should return a vector of partitions (or class labels).

The partition function is applied on rows.

The registered partition method will be used in \code{\link{run_all_consensus_partition_methods}}.

To remove a partition method, use \code{\link{remove_partition_method}}.
}
\value{
No value is returned.
}
\seealso{
\code{\link{all_partition_methods}} lists all registered partition methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
all_partition_methods()
register_partition_fun(random = function(mat, k) sample(k, nrow(mat), replace = TRUE))
all_partition_methods()
remove_partition_method("random")
}
