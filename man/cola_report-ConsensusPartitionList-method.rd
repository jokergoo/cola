\name{cola_report-ConsensusPartitionList-method}
\alias{cola_report,ConsensusPartitionList-method}
\title{
Make HTML report from the ConsensusPartitionList object
}
\description{
Make HTML report from the ConsensusPartitionList object
}
\usage{
\S4method{cola_report}{ConsensusPartitionList}(object, output_dir = getwd(), mc.cores = 1, env = parent.frame())
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{output_dir}{the output directory where the report is put.}
  \item{mc.cores}{number of cores. On OSX it is enforced to be 1.}
  \item{env}{where the objects in the report are found, internally used.}

}
\details{
The \code{\link{ConsensusPartitionList-class}} object contains results for all top-value methods and all partition methods.
This function generates a HTML report which contains all plots and tables for every combination
of top-value method and partition method.

The report generation may take a while because it generates A LOT of heatmaps.

Icon (\url{https://www.flaticon.com/free-icon/can_1366373} ) of the HTML page is made by photo3idea_studio (\url{https://www.flaticon.com/authors/photo3idea-studio} ) from \url{http://www.flaticon.com} licensed by CC 3.0 BY.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rl)
cola_report(cola_rl[c("sd", "MAD"), c("hclust", "skmeans")], output_dir = "~/test_cola_cl_report")
}
}
