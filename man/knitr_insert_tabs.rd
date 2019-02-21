\name{knitr_insert_tabs}
\alias{knitr_insert_tabs}
\title{
Generate the HTML fragment for the JavaScript tabs.
}
\description{
Generate the HTML fragment for the JavaScript tabs.
}
\usage{
knitr_insert_tabs()
}
\details{
The jQuery UI is used to generate html tabs (\url{https://jqueryui.com/tabs/} ).

\code{knitr_insert_tabs} should be used after several callings of \code{\link{knitr_add_tab_item}}
to generate a complete HTML fragment for all tabs with all necessary Javascript and css code.

This function is only for internal use.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
