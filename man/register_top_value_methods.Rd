\name{register_top_value_methods}
\alias{register_top_value_methods}
\title{
Register user-defined top-value methods
}
\description{
Register user-defined top-value methods
}
\usage{
register_top_value_methods(..., validate = TRUE)
}
\arguments{

  \item{...}{A named list of functions.}
  \item{validate}{Whether validate the functions.}

}
\details{
The user-defined function should accept one argument which is the data
matrix where the scores are calculated by rows. Rows with top scores are treated
as "top rows" in cola analysis. Following is how we register "SD" (standard deviation) top-value method:

  \preformatted{
  register_top_value_methods(SD = function(mat) apply(mat, 1, sd))  }

Of course, you can use \code{\link[matrixStats]{rowSds}} to give a faster calculation of row SD:

  \preformatted{
  register_top_value_methods(SD = rowSds)  }

The registered top-value method will be used as defaults in \code{\link{run_all_consensus_partition_methods}}.

To remove a top-value method, use \code{\link{remove_top_value_methods}}.

There are four default top-value methods:

\describe{
  \item{"SD"}{standard deviation, by \code{\link[matrixStats]{rowSds}}.}
  \item{"CV"}{coefficient variance, calculated as \code{sd/(mean+s)} where \code{s} is the 10^th percentile of all row means.}
  \item{"MAD"}{median absolute deviation, by \code{\link[matrixStats:rowSds]{rowMads}}.}
  \item{"ATC"}{the \code{\link{ATC}} method.}
}
}
\value{
No value is returned.
}
\seealso{
\code{\link{all_top_value_methods}} lists all registered top-value methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
all_top_value_methods()
register_top_value_methods(
    ATC_spearman = function(mat) ATC(mat, method = "spearman")
)
all_top_value_methods()
remove_top_value_methods("ATC_spearman")
}
