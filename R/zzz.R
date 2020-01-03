.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/cola/
Github page: https://github.com/jokergoo/cola
Documentation: https://jokergoo.github.io/cola/
Examples: https://jokergoo.github.io/cola_collection/

This message can be suppressed by:
  suppressPackageStartupMessages(library(cola))
========================================
")	

    packageStartupMessage(msg)
}
