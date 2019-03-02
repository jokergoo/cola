\name{cola_report-ConsensusPartitionList-method}
\alias{cola_report,ConsensusPartitionList-method}
\title{
Make report from the ConsensusPartitionList object
}
\description{
Make report from the ConsensusPartitionList object
}
\usage{
\S4method{cola_report}{ConsensusPartitionList}(object, output_dir = getwd(), env = parent.frame())
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{output_dir}{the output directory where the report is put.}
  \item{env}{where the objects in the report are found, internally used.}

}
\details{
The \code{\link{ConsensusPartitionList-class}} object contains results for all top-value methods and all partition methods.
This function generates a HTML report which contains all plots and tables for every combination
of top-value method and partition method.

The report generation may take a while.

Icon (\url{https://www.flaticon.com/free-icon/can_1366373} ) of the HTML page is made by photo3idea_studio (\url{https://www.flaticon.com/authors/photo3idea-studio} ) from www.flaticon.com is licensed by CC 3.0 BY.
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
cola_report(cola_rl[c("sd", "MAD"), c("hclust", "skmeans")], output_dir = "~/test")
}
}
