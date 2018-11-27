\name{register_top_value_method}
\alias{register_top_value_method}
\title{
Register user-defined top-value methods
}
\description{
Register user-defined top-value methods
}
\usage{
register_top_value_method(...)
}
\arguments{

  \item{...}{a named list of functions.}

}
\details{
The user-defined function should accept one argument which is the data
matrix and the scores are calculated by rows. Rows with top scores are treated
as "top rows". Follow is how we register "sd" top-value method:

  \preformatted{
  register_top_value_method(sd = function(mat), apply(mat, 1, sd))  }

Of course, you can use \code{\link[matrixStats]{rowSds}} to give a faster calculation of row sd:

  \preformatted{
  register_top_value_method(sd = rowSds)  }

The registered top-value method will be used as defaults in \code{\link{run_all_consensus_partition_methods}}.

To remove a top-value method, use \code{\link{remove_top_value_method}}.

There are four default top-value methods:

\describe{
  \item{"sd"}{standard deviation, by \code{\link[matrixStats]{rowSds}}.}
  \item{"cv"}{coefficient variance, calculated as \code{sd/(mean+s)} where \code{s} is the 10th quantile of all row means.}
  \item{"MAD"}{median absolute deviation, by \code{\link[matrixStats:rowSds]{rowMads}}.}
  \item{"AAC"}{the \code{\link{AAC}} method.}
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
register_top_value_method(
    AAC_spearman = function(mat) AAC(mat, cor_method = "spearman"),
    AAC_multicore = function(mat) AAC(mat, mc.cores = 2)
)
all_top_value_methods()
remove_top_value_method(c("AAC_spearman", "AAC_multicore"))
}
