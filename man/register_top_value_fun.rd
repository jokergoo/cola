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

To remove a top method, use \code{\link{remove_top_value_method}}.
}
\examples{
ALL_TOP_VALUE_METHOD()
register_top_value_fun(mean = function(mat) rowMeans(mat),
                       median = function(mat) rowMedians(mat))
ALL_TOP_VALUE_METHOD()
remove_top_value_method(c("mean", "median"))
}
