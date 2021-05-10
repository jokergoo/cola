\name{adjust_matrix}
\alias{adjust_matrix}
\title{
Remove rows with low variance and impute missing values
}
\description{
Remove rows with low variance and impute missing values
}
\usage{
adjust_matrix(m, sd_quantile = 0.05, max_na = 0.25, verbose = TRUE)
}
\arguments{

  \item{m}{A numeric matrix.}
  \item{sd_quantile}{Cutoff of the quantile of standard deviation. Rows with standard deviation less than it are removed.}
  \item{max_na}{Maximum NA fraction in each row. Rows with NA fraction larger than it are removed.}
  \item{verbose}{Whether to print messages.}

}
\details{
The function uses \code{\link[impute]{impute.knn}} to impute missing values, then
uses \code{\link{adjust_outlier}} to adjust outliers and 
removes rows with low standard deviations.
}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
m = matrix(rnorm(100), nrow = 10)
m[sample(length(m), 5)] = NA
m[1, ] = 0
m
m2 = adjust_matrix(m)
m2
}
