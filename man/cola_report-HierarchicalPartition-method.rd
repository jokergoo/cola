\name{cola_report-HierarchicalPartition-method}
\alias{cola_report,HierarchicalPartition-method}
\title{
Make report for the HierarchicalPartition object
}
\description{
Make report for the HierarchicalPartition object
}
\usage{
\S4method{cola_report}{HierarchicalPartition}(object, output_dir, env = parent.frame())
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{output_dir}{the output directory where the report is put.}
  \item{env}{where the objects in the report are found, internally used.}

}
\details{
This function generates a HTML report which contains all plots for all nodes
in the partition hierarchy.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rh)
cola_report(cola_rh, output_dir = "~/test2")
}
}
