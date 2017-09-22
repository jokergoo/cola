

KNITR_TAB_ENV = environment()
KNITR_TAB_ENV$current_tab_index = 0
KNITR_TAB_ENV$current_div_index = 0
KNITR_TAB_ENV$header = NULL
KNITR_TAB_ENV$current_html = ""
KNITR_TAB_ENV$random_str = round(runif(1, min = 1, max = 1e8))
KNITR_TAB_ENV$css_added = FALSE

# == title
# Add one tab item
#
# == param
# -code R code
# -header header for the tab
# -desc decription
# -opt options for knitr
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
knitr_add_tab_item = function(code, header, desc = "", opt = NULL) {
	KNITR_TAB_ENV$current_tab_index = KNITR_TAB_ENV$current_tab_index + 1
	tab = qq("tab-@{KNITR_TAB_ENV$random_str}-@{KNITR_TAB_ENV$current_tab_index}")
	knitr_text = qq(
"@{strrep('`', 3)}{r @{tab}@{ifelse(is.null(opt), '', paste(', ', opt))}}
@{code}
@{strrep('`', 3)}

@{desc}
")	
	op = getOption("markdown.HTML.options")
	options(markdown.HTML.options = setdiff(op, c("base64_images", "toc")))
	md = knit(text = knitr_text, quiet = TRUE, envir = parent.frame())
	html = markdownToHTML(text = md, fragment.only = TRUE)
	html = qq("<div id='@{tab}'>\n@{html}\n</div>\n")
	options(markdown.HTML.options = op)
	KNITR_TAB_ENV$header = c(KNITR_TAB_ENV$header, header)
	KNITR_TAB_ENV$current_html = paste0(KNITR_TAB_ENV$current_html, html)
	return(invisible(NULL))
}

# == title
# Indert the HTML tabs
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
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

# == title
# Make report for the ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object
# -output_dir the path for the output html file
# -env where object is found, internally used
#
# == details
# A html report which contains all plots
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "cola_report",
	signature = "ConsensusPartitionList",
	definition = function(object, output_dir, env = parent.frame()) {

	var_name = deparse(substitute(object, env = env))

	txt = readLines("~/project/development/cola/inst/extdata/cola_report_template.Rmd")
	txt = paste(txt, collapse = "\n")
	txt = qq(txt, code.pattern = "\\[%CODE%\\]")

	tempfile = tempfile(tmpdir = output_dir, fileext = ".Rmd")
	writeLines(txt, tempfile)
	op = getOption("markdown.HTML.options")
	options(markdown.HTML.options = setdiff(op, "base64_images"))
	md_file = gsub("Rmd$", "md", tempfile)
	knit(tempfile, md_file)
	markdownToHTML(md_file, paste0(output_dir, "/", "cola_report.html"))
	file.remove(c(tempfile, md_file))
	options(markdown.HTML.options = op)
	qqcat("report is at @{output_dir}/cola_report.html\n")
	return(invisible(NULL))
})


# == title
# Make report for the ConsensusPartition object
#
# == param
# -object a `ConsensusPartition-class` object
# -output_dir the path for the output html file
#
# == details
# A html report which contains all plots
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "cola_report",
	signature = "ConsensusPartition",
	definition = function(object, output_dir) {

	qqcat("Please call `cola_report()` on `ConsensusPartitionList` object directly.\n")
	return(invisible(NULL))
})
