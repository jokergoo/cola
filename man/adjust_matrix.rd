\name{adjust_matrix}
\alias{adjust_matrix}
\title{
Remove rows with low variance and impute missing data
}
\description{
Remove rows with low variance and impute missing data
}
\usage{
adjust_matrix(m, sd_quantile = 0.05, max_na = 0.25)
}
\arguments{

  \item{m}{a numeric matrix}
  \item{sd_quantile}{cutoff the quantile of standard variation}
  \item{max_na}{maximum NA rate for rows}

}
\details{
The function uses \code{\link[impute]{impute.knn}} to impute missing data, then
uses \code{\link{adjust_outlier}} to adjust outliers and 
removes rows with low standard variation.
}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
m = matrix(rnorm(200), 10)
m[1, 1] = 1000
range(m)
m2 = adjust_matrix(m)
range(m2)
}
