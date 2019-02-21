\name{collect_classes-HierarchicalPartition-method}
\alias{collect_classes,HierarchicalPartition-method}
\title{
Collect classes from HierarchicalPartition object
}
\description{
Collect classes from HierarchicalPartition object
}
\usage{
\S4method{collect_classes}{HierarchicalPartition}(object, depth = max_depth(object),
    anno = get_anno(object[1]), anno_col = get_anno_col(object[1]))
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}
  \item{anno}{a data frame with known annotation of the mtarix columns.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}

}
\details{
The function plots the hierarchy of the subgroups.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
collect_classes(cola_rh)
collect_classes(cola_rh, depth = 2)
}
