
# == title
# Overlap of top rows from different top value methods
#
# == param
# -object a `ConsensusPartitionList-class` object.
# -top_n number of top rows.
# -method ``venn``: use Venn diagram; ``venneuler``: use Venn Euler diagram; ``correspondance``: use `correspond_between_rankings`.
# -... additional arguments passed to `venn_euler` or `correspond_between_rankings`.
#
# == value
# No value is returned.
#
# == seealso
# `top_rows_overlap,list-method`
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
		method = c("venn", "venneuler", "correspondance"), ...) {

	all_top_value_list = object@.env$all_top_value_list

	top_rows_overlap(all_top_value_list, top_n = top_n, method = method, ...)
})

# == title
# Overlap of top rows from different top value methods
#
# == param
# -object a numeric matrix.
# -top_value_method methods defined in `all_top_value_methods`.
# -top_n number of top rows.
# -method ``venn``: use Venn diagram; ``venneuler``: use Venn Euler diagram; ``correspondance``: use `correspond_between_rankings`.
# -... additional arguments passed to `venn_euler` or `correspond_between_rankings`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == seealso
# `top_rows_overlap,list-method`
#
# == example
# set.seed(123)
# mat = matrix(rnorm(1000), nrow = 100)
# top_rows_overlap(mat, top_n = 25)
setMethod(f = "top_rows_overlap",
	signature = "matrix",
	definition = function(object, top_value_method = all_top_value_methods(), 
		top_n = round(0.25*nrow(object)), method = c("venn", "venneuler", "correspondance"), ...) {

	all_top_value_list = lapply(top_value_method, function(x) {
		get_top_value_fun = get_top_value_method(x)
		all_top_value = get_top_value_fun(object)
		all_top_value[is.na(all_top_value)] = -Inf
		all_top_value
	})
	names(all_top_value_list) = top_value_method

	top_rows_overlap(all_top_value_list, top_n = top_n, method = method, ...)
})


# == title
# Overlap of top rows from different top value methods
#
# == param
# -object a list which contains rankings from different metrics.
# -top_n number of top rows.
# -method ``venn``: use Venn diagram; ``venneuler``: use Venn Euler diagram; ``correspondance``: use `correspond_between_rankings`.
# -... additional arguments passed to `venn_euler` or `correspond_between_rankings`.
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
# top_rows_overlap(lt, top_n = 25)
setMethod(f = "top_rows_overlap",
	signature = "list",
	definition = function(object, top_n = round(0.25*length(object[[1]])), 
		method = c("venn", "venneuler", "correspondance"), ...) {

	lt = lapply(object, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    method = match.arg(method)
    if(method == "venn") {
		venn(lt, ...)
		title(qq("top @{top_n} rows"))
	} else if(method == "venneuler") {
		venn_euler(lt, main = qq("top @{top_n} rows"), ...)
	} else if(method == "correspondance") {
		correspond_between_rankings(object, top_n = top_n, ...)
	}
})

# == title
# Heatmap of top rows from different top value methods
#
# == param
# -object a `ConsensusPartitionList-class` object.
# -top_n number of top rows.
# -scale_rows wether scale rows. 
# -... pass to `top_rows_heatmap,matrix-method`
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
		scale_rows = object@list[[1]]@scale_rows, ...) {

	all_top_value_list = object@.env$all_top_value_list[object@top_value_method]
    
    mat = object@.env$data

    top_rows_heatmap(mat, all_top_value_list = all_top_value_list, top_n = top_n, 
    	scale_rows = scale_rows, ...)
})

# == title
# Heatmap of top rows from different top methods
#
# == param
# -object a numeric matrix.
# -all_top_value_list top values that have already been calculated from the matrix. If it is ``NULL``
#              the values are calculated by methods in ``top_value_method`` argument.
# -top_value_method methods defined in `all_top_value_methods`.
# -top_n number of top rows to show.
# -layout_nr number of rows in the layout.
# -scale_rows whether scale rows.
#
# == details
# The function makes heatmaps where the rows are scaled (or not scaled) for the top n rows
# from different top value methods.
#
# The top n rows are used for subgroup classification, so the heatmaps show which top
# value method gives better candidate rows for the classification.
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
		top_n = round(0.25*nrow(object)), layout_nr = 1, scale_rows = TRUE) {

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
    pushViewport(viewport(layout = grid.layout(nrow = layout_nr, ncol = ceiling(length(lt)/layout_nr))))
    image_width = 500
	image_height = 500
    for(i in seq_along(lt)) {
    	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
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
				column_title = qq("top @{top_n} rows of @{top_value_method[i]}"),
				show_row_dend = FALSE, show_column_names = FALSE, 
				use_raster = TRUE, raster_quality = 2))
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


# == title
# Make Venn Euler diagram from a list
#
# == param
# -lt a list of vectors.
# -... other arguments passed to `graphics::plot.default`.
#
# == details
# The function calls `gplots::venn` to reformat the data and
# then calls `venneuler::venneuler` to make the plot.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# lt = list(foo = sample(100, 50), bar = sample(100, 50))
# venn_euler(lt)
venn_euler = function(lt, ...) {

	oe = try(has_venneuler <- requireNamespace("venneuler"))
	if(inherits(oe, "try-error")) {
		message("You need to install venneuler package.")
		return(invisible(NULL))
	}
    df = venn(lt, show.plot = FALSE)
    df = df[-1, ]
    set = df[, "num"]
    category = df[, -1]
    names(set) = apply(category, 1, function(x) {
        paste(colnames(category)[as.logical(x)], collapse = "&")
    })
    oe = try(fun <- getFromNamespace("venneuler", "venneuler"))
    if(inherits(oe, "try-error")) {
		return(invisible(NULL))
	}
    plot(fun(set), ...)
}
