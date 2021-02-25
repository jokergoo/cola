\name{hierarchical_partition}
\alias{hierarchical_partition}
\title{
Hierarchical partition
}
\description{
Hierarchical partition
}
\usage{
hierarchical_partition(data,
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(3000, round(nrow(data)*0.3)),
    length.out = 3),
    top_value_method = c("SD", "ATC"),
    partition_method = c("kmeans", "skmeans"),
    combination_method =  expand.grid(top_value_method, partition_method),
    anno = NULL, anno_col = NULL,
    PAC_cutoff = 0.2, min_samples = max(6, round(ncol(data)*0.01)), subset = Inf,
    min_n_signatures = round(nrow(data)*min_p_signatures),
    min_p_signatures = 0.01,
    max_k = 4, verbose = TRUE, mc.cores = 1, cores = mc.cores, help = TRUE, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_n}{Number of rows with top values. Since the number of rows of sub-matrieces might be different, here \code{top_n} can be set as a vector of values less than 1 which are treated as the fraction of the rows from the input matrix.}
  \item{top_value_method}{a single top-value method. Available methods are in \code{\link{all_top_value_methods}}.}
  \item{partition_method}{a single partition method. Available methods are in \code{\link{all_partition_methods}}.}
  \item{combination_method}{A list of combinations of top-value methods and partitioning methods. The value can be a two-column data frame where the first column is the top-value methods and the second column is the partitioning methods. Or it can be a vector of combination names in a form of "top_value_method:partitioning_method".}
  \item{anno}{A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation to predicted subgroups will be tested.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{PAC_cutoff}{the cutoff of PAC scores to determine whether to continue looking for subgroups.}
  \item{min_samples}{the cutoff of number of samples to determine whether to continue looking for subgroups.}
  \item{subset}{Number of columns to randomly sample.}
  \item{min_n_signatures}{Minimal number of signatures under the best classification.}
  \item{min_p_signatures}{Minimal proportion of signatures under the best classification. If the corresponding values are smaller than both \code{min_n_signatures} and \code{min_p_signatures}, the hierarchical partitioning stops on that node.}
  \item{max_k}{maximal number of partitions to try. The function will try \code{2:max_k} partitions. Note this is the number of partitions that will be tried out on each node of the hierarchical partition. Since more subgroups will be found in the whole partition hierarchy, on each node, \code{max_k} should not be set to a large value.}
  \item{verbose}{whether print message.}
  \item{mc.cores}{multiple cores to use. This argument will be removed in future versions.}
  \item{cores}{multiple cores to use. }
  \item{help}{Whether to show the help message.}
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
rh = hierarchical_partition(m, top_value_method = "SD", partition_method = "kmeans")
}
}
