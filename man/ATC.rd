\name{ATC}
\alias{ATC}
\title{
Ability to correlate other rows in the matrix (ATC score)
}
\description{
Ability to correlate other rows in the matrix (ATC score)
}
\usage{
ATC(mat, cor_fun = stat::cor, min_cor = 0, power = 1,
    mc.cores = 1, n_sampling = 1000, q_sd = 0, ...)
}
\arguments{

  \item{mat}{a numeric matrix. ATC score is calculated by rows.}
  \item{cor_fun}{a function which calculates correlation.}
  \item{min_cor}{minimal absolute correlation.}
  \item{power}{power on the correlation values.}
  \item{mc.cores}{number of cores.}
  \item{n_sampling}{when there are too many rows in the matrix, to get the curmulative distribution of how one row correlates other rows, actually we don't need to use all the rows in the matrix, e.g. 1000 rows can already give a farely nice estimation.}
  \item{q_sd}{percentile of the standard deviation for the rows. Rows with values less than it are ignored.}
  \item{...}{pass to \code{cor_fun}, e.g. \code{method = 'spearman'} can be passed to \code{cor_fun} if the correlation function is \code{\link[stats]{cor}}.}

}
\details{
For a given row in a matrix, the ATC score is the area above the curve of the curmulative density
distribution of the absolute correlation to all other rows. Formally, if \code{F_i(X)} is the 
cumulative distribution function of \code{X} where \code{X} is the absolute correlation for row i with power \code{power} (i.e. \code{X^power}),
\code{ATC_i = 1 - \\int_{min_cor}^{1} F_i(X)}.

By default the ATC scores are calculated by Pearson correlation, to use Spearman correlation, you can register
the top-value method by:

  \preformatted{
    register_top_value_methods("ATC_spearman" = function(m) ATC(m, method = "spearman"))  }

Similarly, to use a robust correlation method, e.g. \code{\link[WGCNA]{bicor}} function, you can do like:

  \preformatted{
    register_top_value_methods("ATC_bicor" = function(m) ATC(m, cor_fun = WGCNA::bicor))  }
}
\value{
A vector of numeric values with the same order as rows in the input matrix.
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
ATC_score = ATC(mat)
plot(ATC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
}
