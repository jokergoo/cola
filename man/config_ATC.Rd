\name{config_ATC}
\alias{config_ATC}
\title{
Adjust parameters for default ATC method
}
\description{
Adjust parameters for default ATC method
}
\usage{
config_ATC(cor_fun = stats::cor, min_cor = 0, power = 1, k_neighbours = -1, group = NULL, cores = 1, ...)
}
\arguments{

  \item{cor_fun}{A function that calculates correlations from a matrix (on matrix rows).}
  \item{min_cor}{Cutoff for the minimal absolute correlation.}
  \item{power}{Power on the correlation values.}
  \item{k_neighbours}{Number of the closest neighbours to use.}
  \item{group}{A categorical variable.}
  \item{cores}{Number of cores.}
  \item{...}{Other arguments passed to \code{\link{ATC}}.}

}
\details{
This function changes the default parameters for ATC method. All the arguments in this function all pass to \code{\link{ATC}}.
}
\examples{
# use Spearman correlation
config_ATC(cor_fun = function(m) stats::cor(m, method = "spearman"))
# use knn
config_ATC(k_neighbours = 100)
}
