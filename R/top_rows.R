
# == title
# Overlap of top rows from different top-value methods
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -top_n Number of top rows.
# -method ``euler``: plot Euler diagram by `eulerr::euler`; ``venn``: plot Venn diagram by `gplots::venn`; 
#         ``correspondance``: use `correspond_between_rankings`.
# -... Additional arguments passed to `eulerr::plot.euler` or `correspond_between_rankings`.
#
# == value
# No value is returned.
#
# == seealso
# `top_elements_overlap`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# top_rows_overlap(cola_rl, method = "venn")
# top_rows_overlap(cola_rl, method = "correspondance")
setMethod(f = "top_rows_overlap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = min(object@list[[1]]@top_n), 
		method = c("euler", "venn", "correspondance"), ...) {

	all_top_value_list = object@.env$all_top_value_list[object@top_value_method]

	top_elements_overlap(all_top_value_list, top_n = top_n, method = method, ...)
})

# == title
# Overlap of top rows from different top-value methods
#
# == param
# -object A numeric matrix.
# -top_value_method Methods defined in `all_top_value_methods`.
# -top_n Number of top rows.
# -method ``euler``: plot Euler diagram by `eulerr::euler`; ``venn``: plot Venn diagram by `gplots::venn`; 
#         ``correspondance``: use `correspond_between_rankings`.
# -... Additional arguments passed to `eulerr::plot.euler` or `correspond_between_rankings`.
#
# == details
# It first calculates scores for every top-value method and make plot by `top_elements_overlap`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `top_elements_overlap`
#
# == example
# set.seed(123)
# mat = matrix(rnorm(1000), nrow = 100)
# top_rows_overlap(mat, top_n = 25)
setMethod(f = "top_rows_overlap",
	signature = "matrix",
	definition = function(object, top_value_method = all_top_value_methods(), 
		top_n = round(0.25*nrow(object)), 
		method = c("euler", "venn", "correspondance"), ...) {

	all_top_value_list = lapply(top_value_method, function(x) {
		get_top_value_fun = get_top_value_method(x)
		all_top_value = get_top_value_fun(object)
		all_top_value[is.na(all_top_value)] = -Inf
		all_top_value
	})
	names(all_top_value_list) = top_value_method

	top_elements_overlap(all_top_value_list, top_n = top_n, method = method, ...)
})


# == title
# Overlap of top elements from different metrics
#
# == param
# -object A list which contains values from different metrics.
# -top_n Number of top rows.
# -method ``euler``: plot Euler diagram by `eulerr::euler`; ``venn``: plot Venn diagram by `gplots::venn`; 
#         ``correspondance``: use `correspond_between_rankings`.
# -... Additional arguments passed to `eulerr::plot.euler` or `correspond_between_rankings`.
#
# == details
# The i^th value in every vectors in ``object`` should correspond to the same element from the original data.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(matrixStats)
# set.seed(123)
# mat = matrix(rnorm(1000), nrow = 100)
# lt = list(sd = rowSds(mat), mad = rowMads(mat))
# top_elements_overlap(lt, top_n = 25, method = "venn")
# top_elements_overlap(lt, top_n = 25, method = "correspondance")
top_elements_overlap = function(object, top_n = round(0.25*length(object[[1]])), 
		method = c("euler", "venn", "correspondance"), ...) {

	if(length(unique(sapply(object, length))) > 1) {
		stop_wrap("Length of all vectors in the input list should be the same.")
	}

	lt = lapply(object, function(x) order(x, decreasing = TRUE)[1:top_n])

	if(length(lt) == 1) {
		stop_wrap("Expect at least two lists.")
	}
    
    method = match.arg(method)

    if(method == "venn") {
		gplots::venn(lt, ...)
		title(qq("top @{top_n} rows"))
	} else if(method == "euler") {
		print(plot(eulerr::euler(lt), main = qq("top @{top_n} rows"), legend = TRUE, ...))
	} else if(method == "correspondance") {
		correspond_between_rankings(object, top_n = top_n, ...)
	}
	return(invisible(NULL))
}

# == title
# Heatmap of top rows from different top-value methods
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -top_n Number of top rows.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -scale_rows Wether scale rows. 
# -... Pass to `top_rows_heatmap,matrix-method`
#
# == value
# No value is returned.
#
# == seealso
# `top_rows_heatmap,matrix-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "top_rows_heatmap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = min(object@list[[1]]@top_n), 
	anno = get_anno(object), anno_col = get_anno_col(object),
	scale_rows = object@list[[1]]@scale_rows, ...) {

	all_top_value_list = object@.env$all_top_value_list[object@top_value_method]
    
    mat = object@.env$data

    if(is.null(anno)) {
		bottom_anno = NULL
	} else {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		}

		if(is.null(anno_col)) {
			bottom_anno = HeatmapAnnotation(df = anno,
				show_annotation_name = TRUE, annotation_name_side = "right")
		} else {
			bottom_anno = HeatmapAnnotation(df = anno, col = anno_col,
				show_annotation_name = TRUE, annotation_name_side = "right")
		}
	}

    top_rows_heatmap(mat, all_top_value_list = all_top_value_list, top_n = top_n, 
    	scale_rows = scale_rows, bottom_annotation = bottom_anno, ...)
})

# == title
# Heatmap of top rows from different top-value methods
#
# == param
# -object A numeric matrix.
# -all_top_value_list Top-values that have already been calculated from the matrix. If it is ``NULL``
#              the values are calculated by methods in ``top_value_method`` argument.
# -top_value_method Methods defined in `all_top_value_methods`.
# -bottom_annotation A `ComplexHeatmap::HeatmapAnnotation-class` object.
# -top_n Number of top rows to show in the heatmap.
# -scale_rows Whether scale rows.
#
# == details
# The function makes heatmaps where the rows are scaled (or not scaled) for the top n rows
# from different top-value methods.
#
# The top n rows are used for subgroup classification in cola analysis, so the heatmaps show which 
# top-value method gives better candidate rows for the classification.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == example
# set.seed(123)
# mat = matrix(rnorm(1000), nrow = 100)
# top_rows_heatmap(mat, top_n = 25)
setMethod(f = "top_rows_heatmap",
	signature = "matrix",
	definition = function(object, all_top_value_list = NULL, 
	top_value_method = all_top_value_methods(), 
	bottom_annotation = NULL,
	top_n = round(0.25*nrow(object)), scale_rows = TRUE) {

	if(is.null(all_top_value_list)) {
		all_top_value_list = lapply(top_value_method, function(x) {
			get_top_value_fun = get_top_value_method(x)
			all_top_value = get_top_value_fun(object)
			all_top_value[is.na(all_top_value)] = -Inf
			all_top_value
		})
		names(all_top_value_list) = top_value_method
	} else {
		top_value_method = names(all_top_value_list)
	}

	lt = lapply(all_top_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = length(lt),
    	heights = unit.c(2*max_text_height("foo"), unit(1, "null")))))
    for(i in seq_along(lt)) {
    	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
    	grid.text(qq("top @{top_n} rows of @{top_value_method[i]}"))
    	popViewport()
	}
    image_width = 500
	image_height = 500
    for(i in seq_along(lt)) {
    	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
		file_name = tempfile()
        png(file_name, width = image_width, height = image_height)
		
		mat = object[lt[[i]], ]
		if(nrow(mat) > 5000) {
			mat = mat[sample(nrow(mat), 5000), ]
		}
		if(scale_rows) {
			mat = t(scale(t(mat)))
		}
		if(scale_rows) {
			mat_range = quantile(abs(mat), 0.95, na.rm = TRUE)
			col_fun = colorRamp2(c(-mat_range, 0, mat_range), c("green", "white", "red"))
			heatmap_name = "Z-score"
		} else {
			mat_range = quantile(mat, c(0.05, 0.95))
			col_fun = colorRamp2(c(mat_range[1], mean(mat_range), mat_range[2]), c("blue", "white", "red"))
			heatmap_name = "expr"
		}
		oe = try(
			draw(Heatmap(mat, name = heatmap_name, col = col_fun, show_row_names = FALSE,
				column_title = NULL,
				show_row_dend = FALSE, show_column_names = FALSE, 
				bottom_annotation = bottom_annotation,
				use_raster = TRUE, raster_quality = 2),
				merge_legend = TRUE)
		)
		dev.off2()
	    if(!inherits(oe, "try-error")) {
	    	grid.raster(readPNG(file_name))
	    }
	    grid.rect(gp = gpar(fill = "transparent"))
	    upViewport()
	    if(file.exists(file_name)) file.remove(file_name)
	}
	upViewport()
})

