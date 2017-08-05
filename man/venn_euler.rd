\name{venn_euler}
\alias{venn_euler}
\title{
Make Venn Euler diagram from a list
}
\description{
Make Venn Euler diagram from a list
}
\usage{
venn_euler(lt, ...)
}
\arguments{

  \item{lt}{a list of items}
  \item{...}{other arguments passed to \code{\link[graphics]{plot.default}}}

}
\details{
The function calls \code{\link[gplots]{venn}} to reformat the data and
then calls \code{\link[venneuler]{venneuler}} to make the plot.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
lt = list(foo = sample(100, 50), bar = sample(100, 50))
venn_euler(lt)
}
