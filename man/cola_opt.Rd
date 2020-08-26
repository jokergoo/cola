\name{cola_opt}
\alias{cola_opt}
\title{
Global parameters
}
\description{
Global parameters
}
\usage{
cola_opt(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE)
}
\arguments{

  \item{...}{Arguments for the parameters, see "details" section.}
  \item{RESET}{Whether to reset to default values.}
  \item{READ.ONLY}{Please ignore.}
  \item{LOCAL}{Please ignore.}
  \item{ADD}{Please ignore.}

}
\details{
There are following global parameters:

\describe{
  \item{\code{group_diff}}{Used in \code{\link{get_signatures,ConsensusPartition-method}} to globally control the minimal difference between subgroups.}
  \item{\code{fdr_cutoff}}{Used in \code{\link{get_signatures,ConsensusPartition-method}} to globally control the cutoff of FDR for the differential signature tests.}
  \item{\code{color_set_2}}{Colors for the predicted subgroups.}
}
}
\examples{
cola_opt
cola_opt$group_diff = 0.2  # e.g. for methylation datasets
cola_opt$fdr_cutoff = 0.1  # e.g. for methylation datasets
cola_opt
cola_opt(RESET = TRUE)
}
