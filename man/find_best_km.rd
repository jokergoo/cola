\name{find_best_km}
\alias{find_best_km}
\title{
Find a best k for the k-means clustering
}
\description{
Find a best k for the k-means clustering
}
\usage{
find_best_km(mat, max_km = 15)
}
\arguments{

  \item{mat}{A matrix where k-means clustering is executed by rows.}
  \item{max_km}{Maximal k to try.}

}
\details{
The best k is determined by looking for the knee/elbow of the WSS curve (within-cluster sum of square).

Note this function is only for a rough and quick determination of the best k.
}
\examples{
# There is no example
NULL

}
