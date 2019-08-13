

brewer_pal_set1_col = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
names(brewer_pal_set1_col) = 1:17
brewer_pal_set2_col = c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
names(brewer_pal_set2_col) = 1:16


# == title
# Global Parameters
#
# == param
# -... Arguments for the parameters, see "details" section
# -RESET reset to default values
# -READ.ONLY please ignore
# -LOCAL please ignore
# -ADD please ignore
# 
# == details
# There are following global parameters:
#
# -``group_diff`` Used in `get_signatures,ConsensusPartition,method`.
#
cola_opt = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE) {}
cola_opt = setGlobalOptions(
	raster_resize = list(
		.value = FALSE,
		.visible = FALSE
	),
	group_diff = 0,
	fdr_cutoff = 0.05
)

TEMPLATE_DIR = NULL
.onLoad = function(...) {
	TEMPLATE_DIR <<- system.file("extdata", package = "cola")
}


STAT_USED = c("1-PAC", "mean_silhouette", "concordance")
STAT_ALL = c("1-PAC", "mean_silhouette", "concordance", "cophcor", "aPAC", "FCC")
