\name{cola_report-HierarchicalPartition-method}
\alias{cola_report,HierarchicalPartition-method}
\title{
Make HTML report from the HierarchicalPartition object
}
\description{
Make HTML report from the HierarchicalPartition object
}
\usage{
\S4method{cola_report}{HierarchicalPartition}(object, output_dir = getwd(),
    title = qq("cola Report for Hierarchical Partitioning (@{object[1]@top_value_method}:@{object[1]@partition_method})"),
    env = parent.frame())
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{output_dir}{the output directory where the report is put.}
  \item{mc.cores}{multiple cores to use.}
  \item{title}{Title of the report.}
  \item{env}{where the objects in the report are found, internally used.}

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
# There is no example
NULL

}
