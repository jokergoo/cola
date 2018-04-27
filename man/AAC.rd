\name{AAC}
\alias{AAC}
\title{
AAC score
}
\description{
AAC score
}
\usage{
AAC(mat, cor_method = "pearson", min_cor = 0, max_cor = 1,
    mc.cores = 1, n_sampling = 1000, q_sd = 0)
}
\arguments{

  \item{mat}{a numeric matrix. AAC score is calculated by rows.}
  \item{cor_method}{pass to \code{\link[stats]{cor}}.}
  \item{min_cor}{minimal absolute correlation.}
  \item{max_cor}{maximal absolute correlation.}
  \item{mc.cores}{number of cores.}
  \item{n_sampling}{when the number of columns are too big, to get the curmulative distribution, actually we don't need to use all the rows, e.g. 1000 rows can already give a farely nice estimation for the distribution.}
  \item{q_sd}{percential of the sd for the rows to ignore.}

}
\details{
For a given row in a matrix, the AAC score is the area above the curve of the curmulative density
distribution of the absolute correlation to all other rows. Formally, if \code{F_i(x)} is the 
CDF of the absolute correlation for row i, \code{AAC_i = 1 - \\int_{min_cor}^{max_cor} F_i(x)}.
}
\value{
A vector of numeric values.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(12345)
nr1 = 100
mat1 = matrix(rnorm(100*nr1), nrow = nr1)

nr2 = 10
require(mvtnorm)
sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))

nr3 = 50
sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))

mat = rbind(mat1, mat2, mat3)
AAC_score = AAC(mat)
plot(AAC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
}
