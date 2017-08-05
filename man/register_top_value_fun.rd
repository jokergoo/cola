\name{register_top_value_fun}
\alias{register_top_value_fun}
\title{
Register user-defined top value functions
}
\description{
Register user-defined top value functions
}
\usage{
register_top_value_fun(...)
}
\arguments{

  \item{...}{a named list of functions.}

}
\details{
The user-defined function should only accept one argument which is the data
matrix and the scores are calculated by rows.

The registered top method will be used in \code{\link{run_all_consensus_partition_methods}}.

To remove a top method, use \code{\link{remove_top_value_method}}.
}
\value{
No value is returned.
}
\seealso{
\code{\link{all_top_value_methods}} lists all registered top methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
all_top_value_methods()
register_top_value_fun(AAC_spearman = function(mat) AAC(t(mat), cor_method = "spearman"),
                       AAC_multicore = function(mat) AAC(t(mat), mc.cores = 2))
all_top_value_methods()
remove_top_value_method(c("AAC_spearman", "AAC_multicore"))
}
