\name{knitr_add_tab_item}
\alias{knitr_add_tab_item}
\title{
Add one JavaScript tab in the report
}
\description{
Add one JavaScript tab in the report
}
\usage{
knitr_add_tab_item(code, header, prefix, desc = "", opt = NULL,
    message = NULL, hide_and_show = FALSE)
}
\arguments{

  \item{code}{R code to execute.}
  \item{header}{Header or the title for the tab.}
  \item{prefix}{Prefix of the chunk label.}
  \item{desc}{Decription in the tab.}
  \item{opt}{Options for the knitr chunk.}
  \item{message}{Message to print.}
  \item{hide_and_show}{Whether to hide the code output.}

}
\details{
Each tab contains the R source code and results generated from it (figure, tables, text, ...).

This function is only for internal use.
}
\value{
No value is returned.
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
