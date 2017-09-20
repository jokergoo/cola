cola_report = function(object, output) {
	txt = readLines("~/project/cola/inst/extdata/cola_report_template.Rmd")
	txt = paste(txt, collapse = "\n")
	txt = qq(txt, code.pattern = "\\[%CODE%\\]")

	tempfile = tempfile(tmpdir = getwd(), fileext = ".Rmd")
	writeLines(txt, tempfile)
	knit2html(tempfile, ouptput)
	file.remove(tempfile)
}
