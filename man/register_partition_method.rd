\name{register_partition_method}
\alias{register_partition_method}
\title{
Register user-defined partition functions
}
\description{
Register user-defined partition functions
}
\usage{
register_partition_method(..., scale_method = c("standardization", "rescale", "none"))
}
\arguments{

  \item{...}{a named list of functions.}
  \item{scale_method}{normally, data matrix are scaled by rows before sent to the partition function. The default scaling is applied by \code{\link[base]{scale}}. However, some partition functions may not accept negative values which  are produced by \code{\link[base]{scale}}. Here \code{scale_method} can be set to \code{rescale} which scales rows by \code{(x - min)/(max - min)}. Note here \code{scale_method} only means the method to scale rows. When \code{scale_rows} is set to \code{FALSE} in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}, there wil be no row scaling when doing partition. The value for \code{scale_method} can be a vector if user specifies more than one partition function.}

}
\details{
The user-defined function should accept at least two arguments. The first two arguments are the data
matrix and the number of partitions. The third optional argument should always be \code{...} so that parameters
for the partition function can be passed by \code{partition_param} from \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.
If users forget to add \code{...} in the end, it is added internally.

The function should return a vector of partitions (or class labels) or an object which can be recognized by \code{\link[clue]{cl_membership}}.

The partition function should be applied on columns (Users should be careful with this because some of the R functions apply on rows and
some of the R functions apply on columns). E.g. following is how we register kmeans partition method:

  \preformatted{
  register_partition_method(
      kmeans = function(mat, k, ...) \{
          kmeans(t(mat), centers = k, ...)$centers
      \}
  )  }

The registered partition methods will be used as defaults in \code{\link{run_all_consensus_partition_methods}}.

To remove a partition method, use \code{\link{remove_partition_method}}.

There are following default partition methods:

\describe{
  \item{"hclust"}{hierarchcial clustering with Euclidean distance, later columns are partitioned by \code{\link[stats]{cutree}}. If users want to use another distance metric, consider to register a new partition method. E.g. \code{register_partition_method(hclust_cor = function(mat, k) hc = cutree(hclust(as.dist(cor(mat)))))}.}
  \item{"kmeans"}{by \code{\link[stats]{kmeans}}.}
  \item{"skmeans"}{by \code{\link[skmeans]{skmeans}}.}
  \item{"pam"}{by \code{\link[cluster]{pam}}.}
  \item{"mclust"}{by \code{\link[mclust]{Mclust}}. mclust is applied to the first three principle components from rows.}
  \item{"som"}{by \code{\link[kohonen]{som}}. The SOM map is organized as \code{kr x kr} grids where \code{kr = floor(sqrt(ncol(mat)))}.}
}
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
register_partition_method(
    random = function(mat, k) sample(k, ncol(mat), replace = TRUE)
)
all_partition_methods()
remove_partition_method("random")
}
