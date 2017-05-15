
#' General method for collect_plots
#'
#' @param x x
#' @param ... other arguments 
#'
#' @export
collect_plots = function(x, ...) {
	UseMethod("collect_plots", x)
}


#' Collect plots from run_all object
#'
#' @param x a `run_all` object from `run_all()`.
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
collect_plots.run_all = function(x, k = 2, fun = consensus_heatmap,
	top_method = x$top_method, partition_method = x$partition_method, ...) {

	res_list = x

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
	max_mean_score = 0
	max_i = NULL
	max_j = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(res_list, top_method = top_method[i], partition_method = partition_method[j])

	        # relabel the class according to the class in the first object
	        ik = which(res$k == k)
	        if(is.null(reference_class)) {
	        	reference_class = res$object_list[[ik]]$classification$class
	        } else {
	        	# following elements need to be relabeled
	        	# - res$object_list[[ik]]$classification$class
	        	# - column order of res$object_list[[ik]]$membership
	        	# - res$object_list[[ik]]$membership_each
	        	map = relabel_class(res$object_list[[ik]]$classification$class, reference_class, 1:k)
	        	map2 = structure(names(map), names = map)
	        	res$object_list[[ik]]$classification$class = as.numeric(map[as.character(res$object_list[[ik]]$classification$class)])
	        	res$object_list[[ik]]$membership = res$object_list[[ik]]$membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(res$object_list[[ik]]$membership) = paste0("p", 1:k)
				odim = dim(res$object_list[[ik]]$membership_each)
				res$object_list[[ik]]$membership_each = as.numeric(map[as.character(res$object_list[[ik]]$membership_each)])
				dim(res$object_list[[ik]]$membership_each) = odim
	        }

	        consensus_mat = res$object_list[[ik]]$consensus
			mean_consensus = tapply(seq_along(res$object_list[[ik]]$classification$class), 
				res$object_list[[ik]]$classification$class, function(ind) {
				m = consensus_mat[ind, ind, drop = FALSE]
				if(nrow(m) > 1) {
					mean(m[lower.tri(m)])
				} else {
					NA
				}
			})
			tb = table(res$object_list[[ik]]$classification$class)
			l = is.na(mean_consensus)
			mean_consensus = mean_consensus[!l]
			tb = tb[!l]
			mean_score = sum(mean_consensus*tb)/sum(tb)
			if(mean_score > max_mean_score) {
				max_mean_score = mean_score
				max_i = i
				max_j = j
			}

	        file_name = tempfile()
	        png(file_name)
	        oe = try(fun(res, k = k, show_legend = FALSE, ...))
	        dev.off()
	        if(!inherits(oe, "try-error")) {
		        pushViewport(viewport(layout.pos.row = i+1, layout.pos.col = j+1))
		        grid.raster(readPNG(file_name))
		        grid.rect(gp = gpar(fill = "transparent"))
		        upViewport()
		    }
		    file.remove(file_name)
	    }
	}

	if(!is.null(max_i)) {
		pushViewport(viewport(layout.pos.row = max_i+1, layout.pos.col = max_j+1))
        grid.rect(gp = gpar(fill = "transparent", col = "red", lwd = 2))
        upViewport()
	}
	upViewport()
}


#' Collect plots from consensus_partition object
#'
#' @param x a `consensus_partition` object
#' @param ... other arguments
#'
#' @export
#' @import grid
#' @import png
#'
collect_plots.consensus_partition = function(x, ...) {

	res = x
	all_k = res$k
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = max(c(5, length(all_k))))))

	for(i in 1:4) {
		file_name = tempfile()
        png(file_name)
        oe = try(select_k(res, plot = i))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
	        grid.raster(readPNG(file_name))
	        grid.rect(gp = gpar(fill = "transparent"))
	        upViewport()
	    }
	    file.remove(file_name)
	}

	file_name = tempfile()
    png(file_name)
    oe = try(collect_classes(res))
    dev.off()
    if(!inherits(oe, "try-error")) {
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
        grid.raster(readPNG(file_name))
        grid.rect(gp = gpar(fill = "transparent"))
        upViewport()
    }
    file.remove(file_name)

	for(i in seq_along(all_k)) {
		file_name = tempfile()
        png(file_name)
        oe = try(consensus_heatmap(res, k = all_k[i], show_legend = FALSE, ...))
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
        oe = try(membership_heatmap(res, k = all_k[i], show_legend = FALSE, ...))
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
        oe = try(get_signatures(res, k = all_k[i], show_legend = FALSE, ...))
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


#' General method for collect_classes
#'
#' @param x x
#' @param ... other arguments
#'
#' @export
collect_classes = function(x, ...) {
	UseMethod("collect_classes", x)
}


#' Collect classes from run_all object
#'
#' @param x a `run_all` object returned by [run_all()].
#' @param k number of partitions
#' @param top_method a vector of top methods
#' @param partition_method a vector of partition methods
#' @param ... other arguments.
#'
#' @export
#' @import ComplexHeatmap
collect_classes.run_all = function(x, k, 
	top_method = x$top_method, partition_method = x$partition_method, ...) {

	res_list = x

	top_method_vec = NULL
	partition_method_vec = NULL
	class_df = NULL
	reference_class = NULL
	time_used = NULL
	mem_used = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(res_list, top_method = top_method[i], partition_method = partition_method[j])

	        # relabel the class according to the class in the first object
	        ik = which(res$k == k)
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
#' @param ... other arguments.
#'
#' @export
#' @import ComplexHeatmap
collect_classes.consensus_partition = function(x, ...) {
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
			show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, show_heatmap_legend = FALSE) + 
			Heatmap(class, col = res$object_list[[i]]$class_color, 
				show_row_names = FALSE, show_heatmap_legend = FALSE, name = paste(all_k[i], "_classes"))
		gap = c(gap, c(0, 4))
		class_mat = cbind(class_mat, as.numeric(res$object_list[[i]]$classification$class))
	}
	
    draw(ht_list,gap = unit(gap, "mm"), row_order = do.call("order", as.data.frame(class_mat)))
}
