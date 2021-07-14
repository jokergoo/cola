\name{knee_finder2}
\alias{knee_finder2}
\title{
Find the knee/elbow of a list of sorted points
}
\description{
Find the knee/elbow of a list of sorted points
}
\usage{
knee_finder2(x, plot = FALSE)
}
\arguments{

  \item{x}{A numeric vector.}
  \item{plot}{Whether to make the plot.}

}
\value{
A vector of two numeric values. One for the left knee and the second for the right knee.
}
\examples{
x = rnorm(1000)
knee_finder2(x, plot = TRUE)
}
