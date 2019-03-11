

brewer_pal_set1_col = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
names(brewer_pal_set1_col) = 1:17
brewer_pal_set2_col = c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
names(brewer_pal_set2_col) = 1:16


cola_opt = setGlobalOptions(
	raster_resize = FALSE
)

TEMPLATE_DIR = NULL
.onLoad = function(...) {
	TEMPLATE_DIR <<- system.file("extdata", package = "cola")
}
