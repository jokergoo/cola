\name{ATC_approx}
\alias{ATC_approx}
\title{
Ability to correlate to other rows - an approximated method
}
\description{
Ability to correlate to other rows - an approximated method
}
\usage{
ATC_approx(mat, cor_fun = stats::cor, min_cor = 0, power = 1, k_neighbours = -1,
    mc.cores = 1, cores = mc.cores, n_sampling = c(1000, 500),
    group = NULL, ...)
}
\arguments{

  \item{mat}{A numeric matrix. ATC score is calculated by rows.}
  \item{cor_fun}{A function which calculates correlations on matrix rows.}
  \item{min_cor}{Cutoff for the minimal absolute correlation.}
  \item{power}{Power on the correlation values.}
  \item{k_neighbours}{Nearest k neighbours. Note when this argument is set, there won't be subset sampling for calculating correlations, whihc means, it will calculate correlation to all other rows.}
  \item{mc.cores}{Number of cores. This argument will be removed in future versions.}
  \item{cores}{Number of cores.}
  \item{n_sampling}{When there are too many rows in the matrix, to get the curmulative distribution of how one row correlates other rows, actually we don't need to use all the rows in the matrix, e.g. 1000 rows can already give a very nice estimation.}
  \item{group}{A categorical variable. If it is specified, the correlation is only calculated for the rows in the same group as current row.}
  \item{...}{Pass to \code{cor_fun}.}

}
\details{
For a matrix with huge number of rows. It is not possible to calculate correlation to all other rows, thus the correlation is only
calculated for a randomly sampled subset of othe rows.

With small numbers of rows of the matrix, \code{\link{ATC}} should be used which calculates the "exact" ATC value, but the value of \code{\link{ATC}} and \code{\link{ATC_approx}}
should be very similar.
}
\examples{
# There is no example
NULL

}
