\name{adjust_outlier}
\alias{adjust_outlier}
\title{
Adjust outliers
}
\description{
Adjust outliers
}
\usage{
adjust_outlier(x, q = 0.05)
}
\arguments{

  \item{x}{A numeric vector.}
  \item{q}{Percentile to adjust.}

}
\details{
Vaules larger than percentile \code{1 - q} are adjusted to the \code{1 - q} percentile and 
values smaller than percentile \code{q} are adjusted to the \code{q} percentile
}
\value{
A numeric vector with same length as the original one.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
x = rnorm(40)
x[1] = 100
adjust_outlier(x)
}
