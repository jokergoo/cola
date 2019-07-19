\name{cola_report-HierarchicalPartition-method}
\alias{cola_report,HierarchicalPartition-method}
\title{
Make HTML report from the HierarchicalPartition object
}
\description{
Make HTML report from the HierarchicalPartition object
}
\usage{
\S4method{cola_report}{HierarchicalPartition}(object, output_dir, mc.cores = 1, env = parent.frame())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{output_dir}{The output directory where the report is put.}
  \item{mc.cores}{Multiple cores to use.}
  \item{env}{where The objects in the report are found, internally used.}

}
\details{
This function generates a HTML report which contains all plots for all nodes
in the partition hierarchy.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rh)
cola_report(cola_rh, output_dir = "~/test_cola_rh_report")
}
}
