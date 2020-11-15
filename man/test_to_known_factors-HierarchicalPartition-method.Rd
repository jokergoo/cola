\name{test_to_known_factors-HierarchicalPartition-method}
\alias{test_to_known_factors,HierarchicalPartition-method}
\title{
Test correspondance between predicted classes and known factors
}
\description{
Test correspondance between predicted classes and known factors
}
\usage{
\S4method{test_to_known_factors}{HierarchicalPartition}(object, known = get_anno(object[1]),
    filter_node = filter_node_param(), verbose = FALSE)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{filter_node}{Parameters to merge sub-dendrograms, see \code{\link{filter_node_param}}.}
  \item{known}{A vector or a data frame with known factors. By default it is the annotation table set in \code{\link{hierarchical_partition}}.}
  \item{verbose}{Whether to print messages.}

}
\value{
A data frame with columns:

\itemize{
  \item number of samples
  \item p-values from the tests
  \item number of classes
}

The classifications are extracted for each depth.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
# golub_cola_rh already has known annotations, so test_to_known_factors()
# can be directly applied
test_to_known_factors(golub_cola_rh)
}
