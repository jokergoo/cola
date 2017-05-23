

#' Collect plots from consensus_partition_all_methods object
#'
#' @param x a `consensus_partition_all_methods` object from `run_all()`.
#' @param k number of partitions.
#' @param fun function used to generate plots. Valid functions are [consensus_heatmap()],
#'        [plot_ecdf()], [membership_heatmap()] and [get_signatures()].
#' @param top_method a vector of top methods.
#' @param partition_method a vector of partition methods.
#' @param ... other arguments passed to corresponding `fun`.
#'
#' @export
#' @import grid
#' @import png
#' @import grDevices
#' @import utils
setMethod("collect_plots",
	signature = "ConsensusClusterList",
	definition = function(object, k = 2, fun = consensus_heatmap,
	top_method = object@top_method, partition_method = object@partition_method, ...) {

	reference_class = NULL
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = length(top_method)+1, 
	    ncol = length(partition_method)+1,
	    widths = unit.c(2*grobHeight(textGrob("foo")), unit(rep(1, length(partition_method)), "null")),
	    heights = unit.c(2*grobHeight(textGrob("foo")), unit(rep(1, length(top_method)), "null")))))
	for(i in seq_along(top_method)) {
	    pushViewport(viewport(layout.pos.row = i+1, layout.pos.col = 1))
	    grid.text(top_method[i], rot = 90)
	    upViewport()
	}
	for(j in seq_along(partition_method)) {
	    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j+1))
	    grid.text(partition_method[j])
	    upViewport()
	}
	
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(object, top_method = top_method[i], partition_method = partition_method[j])

	        file_name = tempfile()
	        png(file_name)
	        oe = try(fun(res, k = k, show_legend = FALSE, ...))
	        dev.off()
	        if(!inherits(oe, "try-error")) {
		        pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
		        grid.raster(readPNG(file_name))
		        grid.rect(gp = gpar(fill = "transparent"))
		        upViewport()
		    }
		    file.remove(file_name)
	    }
	}

	upViewport()
}


#' Collect plots from consensus_partition object
#'
#' @param x a `consensus_partition` object
#' @param ... other arguments
#' 
#' @details
#' Plots by [select_k()], [collect_classes.consensus_partition], [consensus_heatmap()], [membership_heatmap()] and [get_signatures()]
#' are arranged in one single page.
#'
#' @export
#' @import grid
#' @import png
#'
setMethod("collect_plots",
	signature = "ConsensusCluster",
	definition = function(object, ...) {

	all_k = object$k
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = max(c(2, length(all_k))))))
	
	# ecdf
	file_name = tempfile()
    png(file_name)
    oe = try(plot_ecdf(object))
    dev.off()
    if(!inherits(oe, "try-error")) {
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
        grid.raster(readPNG(file_name))
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport()
    }
    file.remove(file_name)
	
	file_name = tempfile()
    png(file_name)
    oe = try(collect_classes(object, show_legend = FALSE))
    dev.off()
    if(!inherits(oe, "try-error")) {
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
        grid.raster(readPNG(file_name))
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport()
    }
    file.remove(file_name)

	for(i in seq_along(all_k)) {
		file_name = tempfile()
        png(file_name)
        oe = try(consensus_heatmap(object, k = all_k[i], show_legend = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
	        grid.raster(readPNG(file_name))
	        grid.rect(gp = gpar(fill = "transparent"))
	        upViewport()
	    }
	    file.remove(file_name)

	    file_name = tempfile()
        png(file_name)
        oe = try(membership_heatmap(object, k = all_k[i], show_legend = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
	        grid.raster(readPNG(file_name))
	        grid.rect(gp = gpar(fill = "transparent"))
	        upViewport()
	    }
	    file.remove(file_name)

	    file_name = tempfile()
        png(file_name)
        oe = try(get_signatures(object, k = all_k[i], show_legend = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        pushViewport(viewport(layout.pos.row = 4, layout.pos.col = i))
	        grid.raster(readPNG(file_name))
	        grid.rect(gp = gpar(fill = "transparent"))
	        upViewport()
	    }
	    file.remove(file_name)
	}
	upViewport()

}


#' Collect classes from consensus_partition_all_methods object
#'
#' @param x a `consensus_partition_all_methods` object returned by [run_all()].
#' @param k number of partitions
#' @param top_method a vector of top methods
#' @param partition_method a vector of partition methods
#' @param ... other arguments.
#'
#' @export
#' @import ComplexHeatmap
setMethod("collect_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k, 
	top_method = object@top_method, partition_method = object@partition_method, ...) {

	res_list = x

	top_method_vec = NULL
	partition_method_vec = NULL
	class_df = NULL
	reference_class = NULL
	time_used = NULL
	mem_used = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(object, top_method = top_method[i], partition_method = partition_method[j])

	        # relabel the class according to the class in the first object
	        ik = which(res@k == k)
	        if(is.null(reference_class)) {
	        	reference_class = res$object_list[[ik]]$classification$class
	        } else {
	        	map = relabel_class(res$object_list[[ik]]$classification$class, reference_class, 1:k)
	        	res$object_list[[ik]]$classification$class = as.numeric(map[as.character(res$object_list[[ik]]$classification$class)])
	        }

	        top_method_vec = c(top_method_vec, top_method[i])
	        partition_method_vec = c(partition_method_vec, partition_method[j])
	        class_df = cbind(class_df, res$object_list[[ik]]$classification$class)
	        time_used = c(time_used, attr(res, "system.time")[3])
	        mem_used = c(mem_used, attr(res, "mem_used"))
	    }
	}

	class_df = as.matrix(class_df)
	ht = Heatmap(class_df, name = "class", col = res$object_list[[ik]]$class_color, 
		top_annotation = HeatmapAnnotation(top_method = top_method_vec, partition_method = partition_method_vec,
			time_used = anno_barplot(time_used, axis = TRUE, axis_side = "right", baseline = 0),
			mem_used = anno_barplot(round(mem_used/1024/1024), axis = TRUE, axis_side = "right", baseline = 0),
			show_annotation_name = c(TRUE, TRUE, FALSE, FALSE), annotation_name_side = "right",
			col = list(top_method = structure(names = top_method, brewer.pal(length(top_method), "Set1")),
			           partition_method = structure(names = partition_method, brewer.pal(length(partition_method), "Set2"))),
			annotation_height = unit(c(5, 5, 20, 20), "mm")),
		show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE, cluster_columns = TRUE,
		show_column_dend = FALSE,
		column_order = order(partition_method_vec, top_method_vec))

	class_df2 = get_class(res_list, k = k)

	ht = ht + Heatmap(class_df2[, -ncol(class_df2)], name = "membership", col = colorRamp2(c(0, 1), c("white", "red")),
		show_row_names = FALSE, cluster_columns = FALSE)

	map = relabel_class(class_df2[, ncol(class_df2)], class_df[, 1])
	consensus_class = map[as.character(class_df2[, ncol(class_df2)])]
	ht = ht + Heatmap(consensus_class, col =  res$object_list[[ik]]$class_color, name = "consensus_class", show_row_names = FALSE)

	if(!is.null(res_list$list[[1]]$known)) {
		ht = ht + Heatmap(res_list$list[[1]]$known, name = "known", col = res_list$list[[1]]$known_color,
			show_row_names = FALSE)
	}

	draw(ht, column_title = "classification from all methods", row_title = "samples", row_order = order(consensus_class), cluster_rows = FALSE)
	decorate_annotation("time_used", {
		grid.text("time_used", x = unit(1, "npc") + unit(1.5, "cm"), just = "left")
	})
	decorate_annotation("mem_used", {
		grid.text("mem_used (MB)", x = unit(1, "npc") + unit(1.5, "cm"), just = "left")
	})
}

#' Collect classes from consensus_partition object
#'
#' @param x a `consensus_partition` object
#' @param show_legend whether show legend.
#' @param ... other arguments.
#'
#' @export
#' @import ComplexHeatmap
collect_classes.consensus_partition = function(x, show_legend = TRUE, ...) {
	res = x
	all_k = res$k

	ht_list = NULL
	gap = NULL
	class_mat = NULL
	for(i in seq_along(all_k)) {
		membership = res$object_list[[i]]$membership
		class = res$object_list[[i]]$classification$class
		if(i > 1) {
			map = relabel_class(class, previous_class)
			class = map[as.character(class)]
		}
		previous_class = class
		ht_list = ht_list + Heatmap(membership, col = colorRamp2(c(0, 1), c("white", "red")),
			show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, show_heatmap_legend = i == 1) + 
			Heatmap(class, col = res$object_list[[i]]$class_color, 
				show_row_names = FALSE, show_heatmap_legend = i == length(all_k), name = paste(all_k[i], "_classes"))
		gap = c(gap, c(0, 4))
		class_mat = cbind(class_mat, as.numeric(res$object_list[[i]]$classification$class))
	}
	
    draw(ht_list,gap = unit(gap, "mm"), row_order = do.call("order", as.data.frame(class_mat)),
    	column_title = qq("classes from k = '@{paste(all_k, collapse = ', ')}'"),
    	show_heatmap_legend = show_legend, show_annotation_legend = show_legend)
}


#' Cophenetic correlation coefficient in all methods
#'
#' @param res_list a `consensus_partition_all_methods` object
#' @param top_method a vector of top methods
#' @param partition_method a vector of partition methods
#' 
#' @details
#' A data frame with following columns:
#' 
#' -`cophcor`: the cophenetic correlation coefficient
#' -`mean_silhouette`: mean silhouette
#' -`concordance_to_known`: correlation to known factors if it is provided.
#'
#' @export
#' @importFrom NMF cophcor
collect_stat = function(res_list, top_method = res_list$top_method, 
	partition_method = res_list$partition_method) {
	
	df = data.frame(top_method = character(0), partition_method = character(0), k = numeric(0), cophcor = numeric(0),
		mean_silhouette = numeric(0), concordance_to_known = numeric(0))
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {    
	    	res = get_single_run(res_list, top_method = top_method[i], partition_method = partition_method[j])

	        for(ik in seq_along(res$k)) {
	        	k = res$k[ik]
	        	df = rbind(df, data.frame(top_method = top_method[i], partition_method = partition_method[j], k = k, 
	        		cophcor = cophcor(res$object_list[[ik]]$consensus),
	        		mean_silhouette = mean(res$object_list[[ik]]$classification$silhouette),
	        		concordance_to_known = res$object_list[[ik]]$concordance_to_known))
	        }
	    }
	}
	return(df)
}

