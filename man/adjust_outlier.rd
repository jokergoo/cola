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

  \item{x}{a numeric vector.}
  \item{q}{quantile to adjust.}

}
\details{
Vaules larger than quantile \code{1 - q} are adjusted to the \code{1 - q} quantile and 
values smaller than quantile \code{q} are adjusted to the \code{q} quantile.
}
\value{
A numeric vector with same length as the original one.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
x = rnorm(10)
x[1] = 100
adjust_outlier(x)
}
