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

  \item{lt}{a list of vectors.}
  \item{...}{other arguments passed to \code{\link[graphics]{plot.default}}.}

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
lt = list(a = sample(letters, 13),
          b = sample(letters, 13),
          c = sample(letters, 13))
if(requireNamespace("venneuler")) {
	venn_euler(lt)
}
}
