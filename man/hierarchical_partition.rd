\name{hierarchical_partition}
\alias{hierarchical_partition}
\title{
Hierarchical partition
}
\description{
Hierarchical partition
}
\usage{
hierarchical_partition(data, top_value_method = "MAD", partition_method = "kmeans",
    PAC_cutoff = 0.1, silhouette_cutoff = 0.5,
    min_samples = 6, min_signatures = c(50, 0.05), max_k = 4, verbose = TRUE,
    mc.cores = 1, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{a single top-value method. Available methods are in \code{\link{all_top_value_methods}}.}
  \item{partition_method}{a single partition method. Available methods are in \code{\link{all_partition_methods}}.}
  \item{PAC_cutoff}{the cutoff of PAC scores to determine whether to continue looking for subgroups.}
  \item{silhouette_cutoff}{cutoff for silhouette scores.}
  \item{min_samples}{the cutoff of number of samples to determine whether to continue looking for subgroups.}
  \item{min_signatures}{minimal number of signatures to determine whether to continue looking for subgroups.}
  \item{max_k}{maximal number of partitions to try. The function will try \code{2:max_k} partitions. Note this is the number of partitions that will be tried out on each node of the hierarchical partition. Since more subgroups will be found in the whole partition hierarchy, on each node, \code{max_k} should not be set to a large value.}
  \item{verbose}{whether print message.}
  \item{mc.cores}{multiple cores to use. }
  \item{...}{pass to \code{\link{consensus_partition}}}

}
\details{
The function looks for subgroups in a hierarchical way.

There is a special way to encode the node in the hierarchy. The length of the node name
is the depth of the node in the hierarchy and the substring excluding the last digit is the name
node of the parent node. E.g. for the node \code{0011}, the depth is 4 and the parent node is \code{001}.
}
\value{
A \code{\link{HierarchicalPartition-class}} object. Simply type object in the interactive R session
to see which functions can be applied on it.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
set.seed(123)
m = cbind(rbind(matrix(rnorm(20*20, mean = 2, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
          rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
          rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
                matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20))
         ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
cola_rh = hierarchical_partition(m, top_n = c(20, 30, 40), PAC_cutoff = 0.3)
}
data(cola_rh)
cola_rh
}
