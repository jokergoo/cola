\name{max_depth-HierarchicalPartition-method}
\alias{max_depth,HierarchicalPartition-method}
\alias{max_depth}
\title{
Max depth of the hierarchy
}
\description{
Max depth of the hierarchy
}
\usage{
\S4method{max_depth}{HierarchicalPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}

}
\value{
A numeric value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
max_depth(cola_rh)
}
