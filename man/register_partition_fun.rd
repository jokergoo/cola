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
  \item{scale_method}{normally, data matrix are scaled by rows before sent to the partition function. The default scaling is apply by \code{\link[base]{scale}}. However, some partition function may not accept negative values which may be produced by \code{\link[base]{scale}}. Here \code{scale_method} can be set to \code{rescale} which scale rows by \code{(x - min)/(max - min)}.}

}
\details{
The user-defined function should only accept two arguments which are the data
matrix and the number of partitions. The function should return a vector of
partitions (or group classes).

The partition function is applied on rows.

To remove a partition method, use \code{\link{remove_partition_method}}.
}
\examples{
ALL_PARTITION_METHOD()
register_partition_fun(random = function(mat, k) sample(k, nrow(mat), replace = TRUE))
ALL_PARTITION_METHOD()
remove_partition_method("random")
}
