\name{knitr_add_tab_item}
\alias{knitr_add_tab_item}
\title{
Add one JavaScript tab in the report
}
\description{
Add one JavaScript tab in the report
}
\usage{
knitr_add_tab_item(code, header, desc = "", opt = NULL, message = NULL)
}
\arguments{

  \item{code}{R code to execute.}
  \item{header}{header or the title for the tab.}
  \item{desc}{decription in the tab.}
  \item{opt}{options for the knitr chunk.}
  \item{message}{message to print.}

}
\details{
Each tab contains the R source code and results generated from it (figure, tables, text, ...).

This function in only for internal use.
}
\seealso{
\code{\link{knitr_insert_tabs}} produces a complete HTML fragment.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
