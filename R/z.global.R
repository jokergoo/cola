

brewer_pal_set1_col = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
names(brewer_pal_set1_col) = 1:17
brewer_pal_set2_col = c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
names(brewer_pal_set2_col) = 1:16


# == title
# Global parameters
#
# == param
# -... Arguments for the parameters, see "details" section.
# -RESET Whether to reset to default values.
# -READ.ONLY Please ignore.
# -LOCAL Please ignore.
# -ADD Please ignore.
# 
# == details
# There are following global parameters:
#
# -``group_diff`` Used in `get_signatures,ConsensusPartition-method` to globally control the minimal difference between subgroups.
# -``fdr_cutoff`` Used in `get_signatures,ConsensusPartition-method` to globally control the cutoff of FDR for the differential signature tests.
# -``color_set_2`` Colors for the predicted subgroups.
#
# == example
# cola_opt
# cola_opt$group_diff = 0.2  # e.g. for methylation datasets
# cola_opt$fdr_cutoff = 0.1  # e.g. for methylation datasets
# cola_opt
# cola_opt(RESET = TRUE)
cola_opt = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE) {}
cola_opt = setGlobalOptions(
	group_diff = 0,
	fdr_cutoff = 0.05,
	color_set_1 = brewer_pal_set1_col,
	color_set_2 = brewer_pal_set2_col,
	message = TRUE
)

TEMPLATE_DIR = NULL

# Only for testing
if(Sys.info()["user"] == "jokergoo") {
	TEMPLATE_DIR = "~/project/cola/inst/extdata"
} else {
	TEMPLATE_DIR = "/desktop-home/guz/project/development/cola/inst/extdata"
}


.onLoad = function(...) {
	TEMPLATE_DIR <<- system.file("extdata", package = "cola")
}


STAT_USED = c("1-PAC", "mean_silhouette", "concordance")
STAT_ALL = c("1-PAC", "mean_silhouette", "concordance", "cophcor", "aPAC", "FCC")
