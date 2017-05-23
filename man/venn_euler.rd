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
  \item{...}{other arguments}

}
\details{
The function calls \code{\link[venneuler]{venneuler}} to make the plot
}
\examples{
lt = list(foo = sample(100, 50), bar = sample(100, 50))
venn_euler(lt)
}
