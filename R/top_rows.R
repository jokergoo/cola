
#' Overlap of top rows from different top methods
#'
#' @param res_list a `consensus_partition_all_methods` object
#' @param top_n number of top rows
#' @param type venn: use venn euler plots; correspondance: use [correspond_between_rankings()].
#'
#' @export

setMethod("top_rows_overlap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = 2000, type = c("venn", "correspondance")) {

	all_value_list = object@.env$all_value_list

	top_rows_overlap(all_value_list)
}

setMethod("top_rows_overlap",
	signature = "matrix",
	definition = function(object, top_method = ALL_TOP_VALUE_METHOD(), top_n = 2000,
		type = c("venn", "correspondance")) {

	all_value_list = lapply(top_method, function(x) {
		get_value_fun = get_top_value_fun(x)
		all_value = get_value_fun(object)
		all_value[is.na(all_value)] = -Inf
	})
	names(all_value_list) = top_method

	top_rows_overlap(all_value_list, top_n = top_n, type = type)
})

setMethod("top_rows_overlap",
	signature = "list",
	definition = function(object, top_n = 2000, type = c("venn", "correspondance")) {

	lt = lapply(object, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    if(type == "venn") {
    	if(!requireNamespace("venneuler")) {
    		venn(lt)
    		title(qq("top @{top_n} rows"))
    	} else {
   			venn_euler(lt, main = qq("top @{top_n} rows"))
   		}
	} else if(type == "correspondance") {
		correspond_between_rankings(all_value_list, top_n = top_n)
	}
})


setMethod("top_rows_heatmap",
	signature = "ConsensusPartitionList",
	definition = function(object, top_n = 2000) {

	all_value_list = object@.env$all_value_list
    
    mat = object@.env$data

    top_rows_heatmap(mat, all_value_list = all_value_list, top_n = top_n)
})

setMethod("top_rows_heatmap",
	signature = "matrix",
	definition = function(object, all_value_list = NULL, top_method = ALL_TOP_VALUE_METHOD(), top_n = 2000) {

	if(is.null(all_value_list)) {
		all_value_list = lapply(top_method, function(x) {
			get_value_fun = get_top_value_fun(x)
			all_value = get_value_fun(object)
			all_value[is.na(all_value)] = -Inf
		})
		names(all_value_list) = top_method
	} else {
		top_method = names(all_value_list)
	}

	lt = lapply(all_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    for(i in seq_along(lt)) {
		if(dev.interactive()) {
			readline(prompt = "press enter to load next plot: ")
		}
		mat = object[lt[[i]], ]
		draw(Heatmap(t(scale(t(mat))), name = "scaled_expr", show_row_names = FALSE, 
			column_title = qq("top @{top_n} rows of @{top_method[i]}"),
			show_row_dend = FALSE, show_column_names = FALSE))
	}
})


#' Make Venn Euler diagram from a list
#'
#' @param lt a list of items
#' @param ... other arguments
#'
#' @export
#' @importFrom gplots venn
venn_euler = function(lt, ...) {

	if(!requireNamespace("venneuler")) {
		stop("You need to install venneuler package.")
	}
    foo = venn(lt, show.plot = FALSE)
    foo = foo[-1, ]
    set = foo[, "num"]
    category = foo[, -1]
    names(set) = apply(category, 1, function(x) {
        paste(colnames(category)[as.logical(x)], collapse = "&")
    })
    plot(getFromNamespace("venneuler", "venneuler")(set), ...)
}
