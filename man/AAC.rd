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
    mc.cores = 1, n_sampling = 1000)
}
\arguments{

  \item{mat}{a numeric matrix. AAC score is calculated by columns.}
  \item{cor_method}{pass to \code{\link[stats]{cor}}.}
  \item{min_cor}{minimal absolute correlation.}
  \item{max_cor}{maximal absolute correlation.}
  \item{mc.cores}{number of cores.}
  \item{n_sampling}{when the number of columns are too high, to get the curmulative distribution, actually we don't need to use all the columns, e.g. 1000 columns can already give a farely nice estimation for the distribution.}

}
\details{
AAC score for a given item is the area above the curve of the curmulative 
distribution of the absolute correlation to all other items with \code{x >= min_cor}.
}
\value{
A vector of AAC scores.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(12345)
require(matrixStats)
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
AAC_score = AAC(t(mat))
plot(AAC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
}
