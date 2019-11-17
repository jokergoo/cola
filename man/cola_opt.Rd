\name{cola_opt}
\alias{cola_opt}
\title{
Global Parameters
}
\description{
Global Parameters
}
\usage{
cola_opt(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE)
}
\arguments{

  \item{...}{Arguments for the parameters, see "details" section}
  \item{RESET}{reset to default values}
  \item{READ.ONLY}{please ignore}
  \item{LOCAL}{please ignore}
  \item{ADD}{please ignore}

}
\details{
There are following global parameters:

\describe{
  \item{\code{group_diff}}{Used in \code{\link{get_signatures,ConsensusPartition-method}}.}
  \item{\code{fdr_cutoff}}{Used in \code{\link{get_signatures,ConsensusPartition-method}}.}
  \item{\code{color_set_2}}{Colors for the predicted classes.}
}
}
\examples{
cola_opt
cola_opt$group_diff = 0.2  # e.g. for methylation datasets
cola_opt$fdr_cutoff = 0.1  # e.g. for methylation datasets
cola_opt
cola_opt(RESET = TRUE)
}
