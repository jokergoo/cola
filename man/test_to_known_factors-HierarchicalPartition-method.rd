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
    depth = 2:max_depth(object), verbose = FALSE)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}
  \item{known}{a vector or a data frame with known factors. By default it is the annotation table set in \code{\link{hierarchical_partition}}.}
  \item{verbose}{whether to print messages.}

}
\details{
The function test correlation between classes and known annotations at each node in the hierarchy.

}
\value{
A data frame with columns:

\itemize{
  \item number of samples
  \item p-values from the tests
  \item number of classes
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
test_to_known_factors(cola_rh, known = 1:60)
}
