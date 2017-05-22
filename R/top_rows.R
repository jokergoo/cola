
#' Overlap of top rows from different top methods
#'
#' @param res_list a `consensus_partition_all_methods` object
#' @param top_n number of top rows
#' @param type venn: use venn euler plots; correspondance: use [correspond_between_rankings()].
#'
#' @export
top_rows_overlap = function(res_list, top_n = 2000, type = c("venn", "correspondance")) {

	all_value_list = res_list$list[[1]]$.env$all_value_list

	type = match.arg(type)

    lt = lapply(all_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    if(type == "venn") {
    	if(!requireNamespace("venneuler")) {
    		venn(lt)
    		title(qq("top @{top_n} rows"))
    	} else 
   			venn_euler(lt, main = qq("top @{top_n} rows"))
   		}
	} else if(type == "correspondance") {
		correspond_between_rankings(all_value_list, top_n = top_n)
	}
}

#' Heatmap for the top rows
#'
#' @param res_list a `consensus_partition_all_methods` object
#' @param top_n number of top rows
#'
#' @export
#' @import ComplexHeatmap
top_rows_heatmap = function(res_list, top_n = 2000) {
	
	all_value_list = res_list$list[[1]]$.env$all_value_list
    lt = lapply(all_value_list, function(x) order(x, decreasing = TRUE)[1:top_n])
    
    for(i in seq_along(lt)) {
		if(dev.interactive()) {
			readline(prompt = "press enter to load next plot: ")
		}
		mat = res_list$.env$data[lt[[i]], ]
		draw(Heatmap(t(scale(t(mat))), name = "scaled_expr", show_row_names = FALSE, 
			column_title = qq("top @{length(lt[[i]])} rows of @{res_list$top_method[i]}"),
			show_row_dend = FALSE, show_column_names = FALSE))
	}
}

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
