
# == title
# Collect plots from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object from `run_all_consensus_partition_methods`.
# -k number of partitions.
# -fun function used to generate plots. Valid functions are `consensus_heatmap`,
#        `plot_ecdf`, `membership_heatmap`, `get_signatures` and `dimension_reduction`.
# -top_value_method a vector of top value methods.
# -partition_method a vector of partition methods.
# -... other arguments passed to corresponding ``fun``.
#
# == details
# Plots for all combinations of top value methods and parittion methods are arranged in one single page.
#
# This function makes it easy to directly compare results from multiple methods. 
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# collect_plots(cola_rl, k = 3)
# \dontrun{
# collect_plots(cola_rl, k = 3, fun = membership_heatmap)
# collect_plots(cola_rl, k = 3, fun = get_signatures)
# }
setMethod(f = "collect_plots",
	signature = "ConsensusPartitionList",
	definition = function(object, k = 2, fun = consensus_heatmap,
	top_value_method = object@top_value_method, 
	partition_method = object@partition_method, ...) {

	nv = length(dev.list())
	op = cola_opt$raster_resize
	cola_opt$raster_resize = TRUE
	on.exit({
		nv2 = length(dev.list())
		while(nv2 > nv & nv2 > 1) {
			dev.off2()
			nv2 = length(dev.list())
		}
		cola_opt$raster_resize = op
	})

	fun_name = deparse(substitute(fun))
	grid.newpage()
	pushViewport(viewport(width = unit(1, "npc") - unit(2, "mm"), height = unit(1, "npc") - unit(2, "mm")))
	pushViewport(viewport(layout = grid.layout(nrow = length(top_value_method)+1, 
	    ncol = length(partition_method)+1,
	    widths = unit.c(2*grobHeight(textGrob("foo")), unit(rep(1, length(partition_method)), "null")),
	    heights = unit.c(2*grobHeight(textGrob("foo")), unit(rep(1, length(top_value_method)), "null")))))
	for(i in seq_along(top_value_method)) {
	    pushViewport(viewport(layout.pos.row = i+1, layout.pos.col = 1))
	    grid.text(top_value_method[i], rot = 90)
	    upViewport()
	}
	for(j in seq_along(partition_method)) {
	    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j+1))
	    grid.text(qq("@{partition_method[j]} (k = @{k})"))
	    upViewport()
	}
	
	highlight_row = NULL
	highlight_col = NULL
	for(i in seq_along(top_value_method)) {
	    for(j in seq_along(partition_method)) {
	    	qqcat("* applying @{fun_name} for @{top_value_method[i]}:@{partition_method[j]}\n")
	    	res = object[top_value_method[i], partition_method[j]]
	    	# if(!missing(k)) {
		    # 	if(get_stat(res, k = k)[, "PAC"] < 0.05) {
		    # 		highlight_row = c(highlight_row, i + 1)
		    # 		highlight_col = c(highlight_col, j + 1)
		    # 	}
		    # }
	    	pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
	    	# image_width = convertWidth(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		# image_height = convertHeight(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		image_width = 800
    		image_height = 800
			if(is.null(.ENV$TEMP_DIR)) {
				file_name = tempfile(fileext = ".png", tmpdir = ".")
		        png(file_name, width = image_width, height = image_height)
		        oe = try(fun(res, k = k, internal = TRUE, use_raster = TRUE, ...))
		        dev.off2()
		        if(!inherits(oe, "try-error")) {
					grid.raster(readPNG(file_name))
			    } else {
			    	qqcat("* Caught an error for @{top_value_method[i]}:@{partition_method[j]}:\n@{oe}\n")
			    }
			    if(file.exists(file_name)) file.remove(file_name)
			} else {
				file_name = paste0(.ENV$TEMP_DIR, qq("/@{top_value_method[i]}_@{partition_method[j]}_@{fun_name}_@{k}.png"))
				if(file.exists(file_name)) {
					if(verbose) qqcat("  - use cache png: @{top_value_method[i]}_@{partition_method[j]}_@{fun_name}_@{k}.png\n")
					grid.raster(readPNG(file_name))
				} else {
					png(file_name, width = image_width, height = image_height)
			        oe = try(fun(res, k = k, internal = TRUE, use_raster = TRUE, ...))
			        dev.off2()
			        if(!inherits(oe, "try-error")) {
						grid.raster(readPNG(file_name))
				    } else {
				    	qqcat("* Caught an error for @{top_value_method[i]}:@{partition_method[j]}:\n@{oe}\n")
				    }
				}
			}
			
		    grid.rect(gp = gpar(fill = "transparent", col = "black"))
		    upViewport()
	    }
	}
	# if(!missing(k)) {
	# 	for(i in seq_along(highlight_row)) {
	# 		pushViewport(viewport(layout.pos.row = highlight_row[i], layout.pos.col = highlight_col[i]))
	# 		grid.rect(gp = gpar(fill = "transparent", col = "red", lwd = 2))
	# 		upViewport()
	# 	}
	# }

	upViewport()
	upViewport()
})

# == title
# Collect plots from ConsensusPartition object
#
# == param
# -object a `ConsensusPartition-class` object.
# -verbose whether print messages.
# 
# == details
# Plots by `plot_ecdf`, `collect_classes,ConsensusPartition-method`, `consensus_heatmap`, `membership_heatmap` 
# and `get_signatures` are arranged in one single page, for all avaialble k.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# data(cola_rl)
# collect_plots(cola_rl["sd", "kmeans"])
# }
setMethod(f = "collect_plots",
	signature = "ConsensusPartition",
	definition = function(object, verbose = TRUE) {

	op = cola_opt$raster_resize
	cola_opt$raster_resize = TRUE
	nv = length(dev.list())
	on.exit({
		nv2 = length(dev.list())
		while(nv2 > nv && nv2 > 1) {
			dev.off2()
			nv2 = length(dev.list())
		}
		cola_opt$raster_resize = op
	})

	all_k = object@k
	grid.newpage()
	text_height = grobHeight(textGrob("foo"))
	layout_ncol = 1+max(c(2, length(all_k)))
	pushViewport(viewport(width = unit(1, "npc") - unit(2, "mm"), height = unit(1, "npc") - unit(2, "mm")))
	pushViewport(viewport(layout = grid.layout(nrow = 4+2, ncol = layout_ncol,
		widths = unit.c(2*text_height, unit(rep(1, layout_ncol - 1), "null")),
		heights = unit.c(2*text_height, unit(1, "null"), 2*text_height, unit(rep(1, 3), "null")))))
	
	# first row are two names
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.text("ecdf")
	upViewport()
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
	grid.text("classes")
	upViewport()

	# ecdf
	if(verbose) cat("* plotting empirical cumulative distribution curves of the consensus matrix\n")
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
	file_name = tempfile()
	# image_width = convertWidth(unit(1, "npc"), "bigpts", valueOnly = TRUE)
 #    image_height = convertHeight(unit(1, "npc"), "bigpts", valueOnly = TRUE)
	image_width = 800
	image_height = 800
    png(file_name, width = image_width, height = image_height)
    oe = try(plot_ecdf(object))
    dev.off2()
    if(!inherits(oe, "try-error")) {
    	grid.raster(readPNG(file_name))
    }
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
    if(file.exists(file_name)) file.remove(file_name)
	
	if(verbose) cat("* plotting classes for all k\n")
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
	file_name = tempfile()
    png(file_name, width = image_width, height = image_height)
    oe = try(collect_classes(object, internal = TRUE))
    dev.off2()
    if(!inherits(oe, "try-error")) {
        grid.raster(readPNG(file_name))
    }
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
    if(file.exists(file_name)) file.remove(file_name)

    # pac = get_stat(object, k = all_k)[, "PAC"]
    # border_color = ifelse(pac < 0.1, "red", "black")
    border_color = rep("black", length(all_k))
	
	pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
	grid.text("consensus heatmap", rot = 90)
	upViewport()
	pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 1))
	grid.text("member heatmap", rot = 90)
	upViewport()
	pushViewport(viewport(layout.pos.row = 6, layout.pos.col = 1))
	grid.text("signature heatmap", rot = 90)
	upViewport()

	top_value_method = object@top_value_method
	partition_method = object@partition_method
	for(i in seq_along(all_k)) {
		pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i + 1))
		grid.text(qq("k = @{all_k[i]}"))
		upViewport()

		qqcat("* making consensus heatmap for k = @{all_k[i]}\n")
		pushViewport(viewport(layout.pos.row = 4, layout.pos.col = i + 1))

		if(is.null(.ENV$TEMP_DIR)) {
			file_name = tempfile(fileext = ".png", tmpdir = ".")
	        png(file_name, width = image_width, height = image_height)
	        oe = try(consensus_heatmap(object, k = all_k[i], internal = TRUE))
	        dev.off2()
	        if(!inherits(oe, "try-error")) {
				grid.raster(readPNG(file_name))
		    } else {
		    	qqcat("* Caught an error for consensus_heatmap:@{top_value_method}:@{partition_method}:\n@{oe}\n")
		    }
		    if(file.exists(file_name)) file.remove(file_name)
		} else {
			file_name = paste0(.ENV$TEMP_DIR, qq("/@{top_value_method}_@{partition_method}_consensus_heatmap_@{all_k[i]}.png"))
			if(file.exists(file_name)) {
				if(verbose) qqcat("  - use cache png: @{top_value_method}_@{partition_method}_consensus_heatmap_@{all_k[i]}.png\n")
				grid.raster(readPNG(file_name))
			} else {
				png(file_name, width = image_width, height = image_height)
		        oe = try(consensus_heatmap(object, k = all_k[i], internal = TRUE))
		        dev.off2()
		        if(!inherits(oe, "try-error")) {
					grid.raster(readPNG(file_name))
			    } else {
			    	qqcat("* Caught an error for consensus_heatmap:@{top_value_method}:@{partition_method}:\n@{oe}\n")
			    }
			}
		}

		grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()

	    qqcat("* making membership heatmap for k = @{all_k[i]}\n")
	    pushViewport(viewport(layout.pos.row = 5, layout.pos.col = i + 1))
	    
	    if(is.null(.ENV$TEMP_DIR)) {
			file_name = tempfile(fileext = ".png", tmpdir = ".")
	        png(file_name, width = image_width, height = image_height)
	        oe = try(membership_heatmap(object, k = all_k[i], internal = TRUE))
	        dev.off2()
	        if(!inherits(oe, "try-error")) {
				grid.raster(readPNG(file_name))
		    } else {
		    	qqcat("* Caught an error for membership_heatmap:@{top_value_method}:@{partition_method}:\n@{oe}\n")
		    }
		    if(file.exists(file_name)) file.remove(file_name)
		} else {
			file_name = paste0(.ENV$TEMP_DIR, qq("/@{top_value_method}_@{partition_method}_membership_heatmap_@{all_k[i]}.png"))
			if(file.exists(file_name)) {
				if(verbose) qqcat("  - use cache png: @{top_value_method}_@{partition_method}_membership_heatmap_@{all_k[i]}.png\n")
				grid.raster(readPNG(file_name))
			} else {
				png(file_name, width = image_width, height = image_height)
		        oe = try(membership_heatmap(object, k = all_k[i], internal = TRUE))
		        dev.off2()
		        if(!inherits(oe, "try-error")) {
					grid.raster(readPNG(file_name))
			    } else {
			    	qqcat("* Caught an error for membership_heatmap:@{top_value_method}:@{partition_method}:\n@{oe}\n")
			    }
			}
		}

		grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()

	    qqcat("* making signature heatmap for k = @{all_k[i]}\n")
	    pushViewport(viewport(layout.pos.row = 6, layout.pos.col = i + 1))
	    
	    if(is.null(.ENV$TEMP_DIR)) {
			file_name = tempfile(fileext = ".png", tmpdir = ".")
	        png(file_name, width = image_width, height = image_height)
	        oe = try(get_signatures(object, k = all_k[i], internal = TRUE, use_raster = TRUE))
	        dev.off2()
	        if(!inherits(oe, "try-error")) {
				grid.raster(readPNG(file_name))
		    } else {
		    	qqcat("* Caught an error for get_signatures:@{top_value_method}:@{partition_method}:\n@{oe}\n")
		    }
		    if(file.exists(file_name)) file.remove(file_name)
		} else {
			file_name = paste0(.ENV$TEMP_DIR, qq("/@{top_value_method}_@{partition_method}_get_signatures_@{all_k[i]}.png"))
			if(file.exists(file_name)) {
				if(verbose) qqcat("  - use cache png: @{top_value_method}_@{partition_method}_get_signatures_@{all_k[i]}.png\n")
				grid.raster(readPNG(file_name))
			} else {
				png(file_name, width = image_width, height = image_height)
		        oe = try(get_signatures(object, k = all_k[i], internal = TRUE, use_raster = TRUE))
		        dev.off2()
		        if(!inherits(oe, "try-error")) {
					grid.raster(readPNG(file_name))
			    } else {
			    	qqcat("* Caught an error for get_signatures:@{top_value_method}:@{partition_method}:\n@{oe}\n")
			    }
			}
		}

		grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()
	}
	upViewport()
	upViewport()
})

# == title
# Collect classes from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object returned by `run_all_consensus_partition_methods`.
# -k number of partitions
# -top_value_method a vector of top value methods
# -partition_method a vector of partition methods
#
# == details
# There are following panels in the plot:
#
# - a heatmap shows partitions predicted from all methods where the top annotation
#   is the consensus partition summarized from partitions from all methods, weighted
#   by mean silhouette scores.
# - a row barplot annotation showing the mean silhouette scores for different methods.
# - a heatmap shows the similarities of the partitions of pairwise methods, calculated
#    by `clue::cl_dissimilarity` with the ``comembership`` method.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# data(cola_rl)
# collect_classes(cola_rl, k = 3)
# }
setMethod(f = "collect_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k, 
	top_value_method = object@top_value_method, partition_method = object@partition_method) {

	top_value_method_vec = NULL
	partition_method_vec = NULL
	class_mat = NULL
	silhouette_mat = NULL
	for(i in seq_along(top_value_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = object[top_value_method[i], partition_method[j]]

	        top_value_method_vec = c(top_value_method_vec, top_value_method[i])
	        partition_method_vec = c(partition_method_vec, partition_method[j])
	        class_df = get_classes(res, k)
	        class_mat = cbind(class_mat, class_df[, "class"])
	        silhouette_mat = cbind(silhouette_mat, class_df[, "silhouette"])
	    }
	}

	class_mat = as.matrix(class_mat)
	colnames(class_mat) = paste(top_value_method_vec, partition_method_vec, sep = ":")
	ik = which(res@k == k)
	
	silhouette_mat = as.matrix(silhouette_mat)
	silhouette_mat[silhouette_mat < 0] = 0

	adjust_by_transparency = function(col, transparency) {
		rgb( 1 - (1 - t(col2rgb(col)/255)) * (1 - transparency))
	}

	pac = get_stat(object, k)[, "PAC"][colnames(class_mat)]
	consensus_class = get_classes(object, k = k)$class
	m = t(class_mat)
	column_order = column_order_by_group(consensus_class, m)
	if(is.null(object@list[[1]]@anno)) {
		bottom_annotation = NULL
	} else {
		bottom_annotation = HeatmapAnnotation(df = object@list[[1]]@anno, 
			col = object@list[[1]]@anno_col, show_annotation_name = TRUE,
			annotation_name_side = "left")
	}
	pl = lapply(object@list[paste(top_value_method_vec, partition_method_vec, sep = ":")], function(x) as.cl_partition(get_membership(x, k = k)))
	clen = cl_ensemble(list = pl)
	m_diss = cl_dissimilarity(clen, method = "comembership")

	ht = Heatmap(m, name = "Class", col = brewer_pal_set2_col, column_order = column_order,
		column_title = qq("classification from all methods, k = @{k}"),
		row_names_side = "left", cluster_rows = hclust(m_diss), cluster_columns = FALSE,
		show_column_dend = FALSE, rect_gp = gpar(type = "none"),
		cell_fun = function(j, i, x, y, w, h, fill) {
			col = adjust_by_transparency(fill, 1 - silhouette_mat[j, i])
			grid.rect(x, y, w, h, gp = gpar(fill = col, col = col))
		}, top_annotation = HeatmapAnnotation(consensus_class = consensus_class, 
			col = list(consensus_class = brewer_pal_set2_col),
			show_annotation_name = TRUE, annotation_name_side = "left", show_legend = FALSE),
		bottom_annotation = bottom_annotation, width = 1) + 
		rowAnnotation(mean_silhouette = row_anno_barplot(get_stat(object, k = k)[colnames(class_mat), "mean_silhouette"], baseline = 0, axis = TRUE),
			width = unit(2, "cm")) +
		rowAnnotation("Top value method" = top_value_method_vec, 
			"Partition method" = partition_method_vec,
			show_annotation_name = FALSE, annotation_name_side = "bottom",
			col = list("Top value method" = structure(names = top_value_method, brewer_pal_set1_col[seq_along(top_value_method)]),
			           "Partition method" = structure(names = partition_method, brewer_pal_set2_col[seq_along(partition_method)])),
			width = unit(10, "mm"))
	
	draw(ht, padding = unit(c(15, 1, 1, 1), "mm"))
	decorate_annotation("mean_silhouette", {
		grid.text("mean_silhouette", y = unit(-1, "cm"))
	})
})


# == title
# Collect classes from ConsensusPartition object
#
# == param
# -object a `ConsensusPartition-class` object.
# -internal used internally.
# -show_row_names whether show row names
#
# == details
# Membership matrix and the classes with each k are plotted in the heatmap.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# collect_classes(cola_rl["sd", "kmeans"])
setMethod(f = "collect_classes",
	signature = "ConsensusPartition",
	definition = function(object, internal = FALSE, show_row_names = FALSE) {

	all_k = object@k

	ht_list = NULL
	gap = NULL
	class_mat = NULL
	for(i in seq_along(all_k)) {
		membership = object@object_list[[i]]$membership
		class = object@object_list[[i]]$class_df$class

		ht_list = ht_list + Heatmap(membership, col = colorRamp2(c(0, 1), c("white", "red")),
			show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, show_heatmap_legend = i == 1,
			heatmap_legend_param = list(title = "Prob"),
			column_title = qq("k = @{all_k[i]}"),
			name = paste0("membership_", all_k[i])) + 
			Heatmap(class, col = brewer_pal_set2_col, 
				show_row_names = FALSE, show_heatmap_legend = i == length(all_k), 
				heatmap_legend_param = list(title = "Class"),
				name = paste(all_k[i], "_classes"))
		gap = c(gap, c(0, 4))
		class_mat = cbind(class_mat, as.numeric(class))
	}

	if(!internal & show_row_names) {
		rn = rownames(membership)
		ht_list = ht_list + rowAnnotation(nm = row_anno_text(rn, location = unit(0, "npc"), just = "left"), width = max_text_width(rn))
	}

    draw(ht_list, gap = unit(gap, "mm"), row_order = do.call("order", as.data.frame(class_mat)),
    	# column_title = qq("classes from k = '@{paste(all_k, collapse = ', ')}'"),
    	show_heatmap_legend = !internal, show_annotation_legend = !internal)

    for(k in all_k) {
    	ik = which(all_k == k )
    	# border_color = ifelse(object@object_list[[ik]]$stat$PAC < 0.1, "red", "black")
    	border_color = rep("black", length(all_k))
    	decorate_heatmap_body(paste0("membership_", k), {
    		grid.rect(0, width = unit(1+1/k, "npc"), just = "left", gp = gpar(col = border_color, fill = "transparent"))
    	})
    }
})


