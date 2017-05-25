
# == title
# Overlap of top rows from different top methods
#
# == param
# -object a `ConsensusPartitionList-class` object
# -top_n number of top rows
# -type ``venn``: use venn euler plots; ``correspondance``: use `correspond_between_rankings`.
#
setMethod(f = "top_rows_overlap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = round(0.25*length(all_value_list[[1]])), 
		type = c("venn", "correspondance")) {

	all_value_list = object@.env$all_value_list

	top_rows_overlap(all_value_list, type = type)
})

# == title
# Overlap of top rows from different top methods
#
# == param
# -object a numeric matrix
# -top_method methods defined in `ALL_TOP_VALUE_METHOD`.
# -top_n number of top rows
# -type ``venn``: use venn euler plots; ``correspondance``: use `correspond_between_rankings`.
#
# == example
# mat = matrix(rnorm(1000), nrow = 100)
# top_rows_overlap(mat, top_n = 25)
setMethod(f = "top_rows_overlap",
	signature = "matrix",
	definition = function(object, top_method = ALL_TOP_VALUE_METHOD(), top_n = round(0.25*nrow(object)),
		type = c("venn", "correspondance")) {

	all_value_list = lapply(top_method, function(x) {
		get_value_fun = get_top_value_fun(x)
		all_value = get_value_fun(object)
		all_value[is.na(all_value)] = -Inf
		all_value
	})
	names(all_value_list) = top_method

	top_rows_overlap(all_value_list, top_n = top_n, type = type)
})


# == title
# Overlap of top rows from different top methods
#
# == param
# -object a list which contains rankings from different metrics.
# -top_n number of top rows
# -type ``venn``: use venn euler plots; ``correspondance``: use `correspond_between_rankings`.
#
# == example
# require(matrixStats)
# mat = matrix(rnorm(1000), nrow = 100)
# lt = list(sd = rowSds(mat), mad = rowMads(mat))
# top_rows_overlap(lt, top_n = 25)
setMethod(f = "top_rows_overlap",
	signature = "list",
	definition = function(object, top_n = round(0.25*length(object[[1]])), 
		type = c("venn", "correspondance")) {

	lt = lapply(object, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    type = match.arg(type)
    if(type == "venn") {
    	oe = try(has_venneuler <- requireNamespace("venneuler"))
    	if(inherits(oe, "try-error")) {
    		venn(lt)
    		title(qq("top @{top_n} rows"))
    	} else if(!has_venneuler) {
    		venn(lt)
    		title(qq("top @{top_n} rows"))
    	} else {
   			venn_euler(lt, main = qq("top @{top_n} rows"))
   		}
	} else if(type == "correspondance") {
		correspond_between_rankings(object, top_n = top_n)
	}
})

# == title
# Heatmap of top rows from different top methods
#
# == param
# -object a `ConsensusPartitionList-class` object
# -top_n number of top rows
#
setMethod(f = "top_rows_heatmap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = round(0.25*length(all_value_list[[1]]))) {

	all_value_list = object@.env$all_value_list
    
    mat = object@.env$data

    top_rows_heatmap(mat, all_value_list = all_value_list, top_n = top_n)
})

# == title
# Heatmap of top rows from different top methods
#
# == param
# -object a numeric matrix
# -all_value_list scores that have already been calculated from the matrix
# -top_method methods defined in `ALL_TOP_VALUE_METHOD`.
# -top_n number of top rows
#
# == example
# mat = matrix(rnorm(1000), nrow = 100)
# top_rows_heatmap(mat, top_n = 25)
setMethod(f = "top_rows_heatmap",
	signature = "matrix",
	definition = function(object, all_value_list = NULL, top_method = ALL_TOP_VALUE_METHOD(), 
		top_n = round(0.25*nrow(object))) {

	if(is.null(all_value_list)) {
		all_value_list = lapply(top_method, function(x) {
			get_value_fun = get_top_value_fun(x)
			all_value = get_value_fun(object)
			all_value[is.na(all_value)] = -Inf
			all_value
		})
		names(all_value_list) = top_method
	} else {
		top_method = names(all_value_list)
	}

	lt = lapply(all_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    for(i in seq_along(lt)) {
		if(dev.interactive() && interactive()) {
			readline(prompt = "press enter to load next plot: ")
		}
		mat = object[lt[[i]], ]
		draw(Heatmap(t(scale(t(mat))), name = "scaled_expr", show_row_names = FALSE, 
			column_title = qq("top @{top_n} rows of @{top_method[i]}"),
			show_row_dend = FALSE, show_column_names = FALSE))
	}
})


# == title
# Make Venn Euler diagram from a list
#
# == param
# -lt a list of items
# -... other arguments
#
# == details
# The function calls `venneuler::venneuler` to make the plot
#
# == example
# lt = list(foo = sample(100, 50), bar = sample(100, 50))
# venn_euler(lt)
venn_euler = function(lt, ...) {

	oe = try(has_venneuler <- requireNamespace("venneuler"))
	if(inherits(oe, "try-error")) {
		message("You need to install venneuler package.")
		return(invisible(NULL))
	}
    foo = venn(lt, show.plot = FALSE)
    foo = foo[-1, ]
    set = foo[, "num"]
    category = foo[, -1]
    names(set) = apply(category, 1, function(x) {
        paste(colnames(category)[as.logical(x)], collapse = "&")
    })
    oe = try(fun <- getFromNamespace("venneuler", "venneuler"))
    if(inherits(oe, "try-error")) {
		return(invisible(NULL))
	}
    plot(fun(set), ...)
}
