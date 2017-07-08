\name{adjust_matrix}
\alias{adjust_matrix}
\title{
Remove rows with low variance
}
\description{
Remove rows with low variance
}
\usage{
adjust_matrix(m, sd_quantile = 0.05)
}
\arguments{

  \item{m}{a numeric matrix}
  \item{sd_quantile}{cutoff the quantile of standard variation}

}
\details{
The function first uses \code{\link{adjust_outlier}} to adjust outliers and 
then removes rows with low standard variation.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
