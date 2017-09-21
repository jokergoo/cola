

KNITR_TAB_ENV = environment()
KNITR_TAB_ENV$current_tab_index = 0
KNITR_TAB_ENV$current_div_index = 0
KNITR_TAB_ENV$header = NULL
KNITR_TAB_ENV$current_html = ""
KNITR_TAB_ENV$random_str = round(runif(1, min = 1, max = 1e8))
KNITR_TAB_ENV$css_added = FALSE

knitr_add_tab_item = function(code, header, desc = "") {
	KNITR_TAB_ENV$current_tab_index = KNITR_TAB_ENV$current_tab_index + 1
	tab = qq("tab-@{KNITR_TAB_ENV$random_str}-@{KNITR_TAB_ENV$current_tab_index}")
	knitr_text = qq(
"@{strrep('`', 3)}{r @{tab}}
@{code}
@{strrep('`', 3)}

@{desc}
")	
	md = knit(text = knitr_text, quiet = TRUE, envir = parent.frame())
	html = markdownToHTML(text = md, fragment.only = TRUE, options = setdiff(options('markdown.HTML.options')[[1]], "toc"))
	html = qq("<div id='@{tab}'>\n@{html}\n</div>\n")
	
	KNITR_TAB_ENV$header = c(KNITR_TAB_ENV$header, header)
	KNITR_TAB_ENV$current_html = paste0(KNITR_TAB_ENV$current_html, html)
	return(invisible(NULL))
}

knitr_insert_tabs = function() {
	KNITR_TAB_ENV$current_div_index = KNITR_TAB_ENV$current_div_index + 1

	if(!KNITR_TAB_ENV$css_added) {
		css = readLines("https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css")
		cat("<style type='text/css'>\n")
		cat(css, sep = "\n")
		cat("</style>\n")
		cat("<script src='https://code.jquery.com/jquery-1.12.4.js'></script>\n")
		cat("<script src='https://code.jquery.com/ui/1.12.1/jquery-ui.js'></script>\n")
	}
	qqcat("
<script>
$( function() {
	$( '#tabs@{KNITR_TAB_ENV$current_div_index}' ).tabs();
} );
</script>
")
	qqcat("<div id='tabs@{KNITR_TAB_ENV$current_div_index}'>\n")
	cat("<ul>\n")
	qqcat("<li><a href='#tab-@{KNITR_TAB_ENV$random_str}-@{seq_len(KNITR_TAB_ENV$current_tab_index)}'>@{KNITR_TAB_ENV$header}</a></li>\n")
	cat("</ul>\n")
	cat(KNITR_TAB_ENV$current_html)
	cat("</div>\n")

	KNITR_TAB_ENV$current_tab_index = 0
	KNITR_TAB_ENV$header = NULL
	KNITR_TAB_ENV$current_html = ""
	KNITR_TAB_ENV$random_str = round(runif(1, min = 1, max = 1e8))
	KNITR_TAB_ENV$css_added = TRUE
	
	return(invisible(NULL))
}

cola_report = function(object, output) {

	txt = readLines("~/project/development/cola/inst/extdata/cola_report_template.Rmd")
	txt = paste(txt, collapse = "\n")
	txt = qq(txt, code.pattern = "\\[%CODE%\\]")

	tempfile = tempfile(tmpdir = getwd(), fileext = ".Rmd")
	writeLines(txt, tempfile)
	knit2html(tempfile, output)
	#file.remove(tempfile)
}
