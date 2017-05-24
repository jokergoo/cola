
# == title
# Collect plots from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object from `run_all_consensus_partition_methods`.
# -k number of partitions.
# -fun function used to generate plots. Valid functions are `consensus_heatmap,ConsensusPartition-method`,
#        `plot_ecdf,ConsensusPartition-method`, `membership_heatmap,ConsensusPartition-method` and 
#        `get_signatures,ConsensusPartition-method`.
# -top_method a vector of top methods.
# -partition_method a vector of partition methods.
# -... other arguments passed to corresponding ``fun``.
#
setMethod(f = "collect_plots",
	signature = "ConsensusPartitionList",
	definition = function(object, k, fun = consensus_heatmap,
	top_method = object@top_method, partition_method = object@partition_method, ...) {

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

	    	pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
	    	# image_width = convertWidth(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		# image_height = convertHeight(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		image_width = 800
    		image_height = 800
	        file_name = tempfile()
	        png(file_name, width = image_width, height = image_height)
	        oe = try(fun(res, k = k, show_legend = FALSE, show_column_names = FALSE, use_raster = FALSE, ...))
	        dev.off()
	        if(!inherits(oe, "try-error")) {
		        grid.raster(readPNG(file_name))
		    }
		    grid.rect(gp = gpar(fill = "transparent"))
		    upViewport()
		    file.remove(file_name)
	    }
	}

	upViewport()
})

# == title
# Collect plots from ConsensusPartition object
#
# == param
# -object a `ConsensusPartition-class` object
# -... other arguments
# 
# == details
# Plots by `plot_ecdf,ConsensusPartition-method`, `collect_classes,ConsensusPartition-method`, `consensus_heatmap,ConsensusPartition-method`, `membership_heatmap,ConsensusPartition-method` 
# and `get_signatures,ConsensusPartition-method` are arranged in one single page.
#
setMethod(f = "collect_plots",
	signature = "ConsensusPartition",
	definition = function(object, ...) {

	all_k = object@k
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = max(c(2, length(all_k))))))
	
	# ecdf
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
	file_name = tempfile()
	# image_width = convertWidth(unit(1, "npc"), "bigpts", valueOnly = TRUE)
 #    image_height = convertHeight(unit(1, "npc"), "bigpts", valueOnly = TRUE)
	image_width = 800
	image_height = 800
    png(file_name, width = image_width*2, height = image_height*2)
    oe = try(plot_ecdf(object))
    dev.off()
    if(!inherits(oe, "try-error")) {
    	grid.raster(readPNG(file_name))
    }
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
    file.remove(file_name)
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
	file_name = tempfile()
    png(file_name, width = image_width*2, height = image_height*2)
    oe = try(collect_classes(object, show_legend = FALSE))
    dev.off()
    if(!inherits(oe, "try-error")) {
        grid.raster(readPNG(file_name))
    }
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
    file.remove(file_name)

	for(i in seq_along(all_k)) {
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
		file_name = tempfile()
        png(file_name, width = image_width*2, height = image_height*2)
        oe = try(consensus_heatmap(object, k = all_k[i], show_legend = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))  
	    }
	    grid.rect(gp = gpar(fill = "transparent"))
	    upViewport()
	    file.remove(file_name)

	    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
	    file_name = tempfile()
        png(file_name, width = image_width*2, height=  image_height*2)
        oe = try(membership_heatmap(object, k = all_k[i], show_legend = FALSE, show_column_names = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))
	    }
	    grid.rect(gp = gpar(fill = "transparent"))
	    upViewport()
	    file.remove(file_name)

	    pushViewport(viewport(layout.pos.row = 4, layout.pos.col = i))
	    file_name = tempfile()
        png(file_name, width = image_width*2, height=  image_height*2)
        oe = try(get_signatures(object, k = all_k[i], show_legend = FALSE, show_column_names = FALSE, use_raster = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))  
	    }
	    grid.rect(gp = gpar(fill = "transparent"))
	    upViewport()
	    file.remove(file_name)
	}
	upViewport()

})

# == title
# Collect classes from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object returned by `run_all_consensus_partition_methods`.
# -k number of partitions
# -top_method a vector of top methods
# -partition_method a vector of partition methods
# -... other arguments.
#
setMethod(f = "collect_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k, 
	top_method = object@top_method, partition_method = object@partition_method, ...) {

	top_method_vec = NULL
	partition_method_vec = NULL
	class_df = NULL
	time_used = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(object, top_method = top_method[i], partition_method = partition_method[j])

	        top_method_vec = c(top_method_vec, top_method[i])
	        partition_method_vec = c(partition_method_vec, partition_method[j])
	        class_df = cbind(class_df, get_class(res, k)[, "class"])
	        time_used = c(time_used, attr(res, "system.time")[3])
	    }
	}

	class_df = as.matrix(class_df)
	colnames(class_df) = paste(top_method_vec, partition_method_vec, sep = "/")
	class_col = structure(unique(res@object_list[[1]]$class_df$class_col), names = unique(res@object_list[[1]]$class_df$class))
	class_col = class_col[order(names(class_col))]
	ht = Heatmap(class_df, name = "class", col = class_col, 
		top_annotation = HeatmapAnnotation(top_method = top_method_vec, partition_method = partition_method_vec,
			time_used = anno_barplot(time_used, axis = TRUE, axis_side = "right", baseline = 0),
			show_annotation_name = c(TRUE, TRUE, FALSE, FALSE), annotation_name_side = "right",
			col = list(top_method = structure(names = top_method, brewer_pal_set1_col[seq_along(top_method)]),
			           partition_method = structure(names = partition_method, brewer_pal_set2_col[seq_along(partition_method)])),
			annotation_height = unit(c(5, 5, 20), "mm")),
		show_row_names = FALSE, show_row_dend = FALSE, cluster_columns = TRUE,
		show_column_dend = FALSE,
		column_order = order(partition_method_vec, top_method_vec))

	class_df2 = get_class(object, k = k)

	ht = ht + Heatmap(class_df2[, -ncol(class_df2)], name = "membership", col = colorRamp2(c(0, 1), c("white", "red")),
		show_row_names = FALSE, cluster_columns = FALSE)

	ht = ht + Heatmap(class_df2[, "class"], col =  class_col, name = "consensus_class", show_row_names = FALSE)

	draw(ht, column_title = "classification from all methods", row_title = qq("@{nrow(class_df)} samples"), row_order = order(class_df2[, "class"]), 
		cluster_rows = FALSE)
	decorate_annotation("time_used", {
		grid.text("time_used (s)", x = unit(1, "npc") + unit(1, "cm"), just = "left")
	})
})

# == title
# Collect classes from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object
# -show_legend whether show legend.
# -... other arguments.
#
setMethod(f = "collect_classes",
	signature = "ConsensusPartition",
	definition = function(object, show_legend = TRUE, ...) {

	all_k = object@k

	ht_list = NULL
	gap = NULL
	class_mat = NULL
	for(i in seq_along(all_k)) {
		membership = object@object_list[[i]]$membership
		class = object@object_list[[i]]$class_df$class
		class_col = object@object_list[[i]]$class_df$class_col
		class_col = structure(unique(class_col), names = unique(class))
		class_col = class_col[order(names(class_col))]

		ht_list = ht_list + Heatmap(membership, col = colorRamp2(c(0, 1), c("white", "red")),
			show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, show_heatmap_legend = i == 1,
			name = paste0("membership_", all_k[i])) + 
			Heatmap(class, col = class_col, 
				show_row_names = FALSE, show_heatmap_legend = i == length(all_k), name = paste(all_k[i], "_classes"))
		gap = c(gap, c(0, 4))
		class_mat = cbind(class_mat, as.numeric(class))
	}

    draw(ht_list,gap = unit(gap, "mm"), row_order = do.call("order", as.data.frame(class_mat)),
    	column_title = qq("classes from k = '@{paste(all_k, collapse = ', ')}'"),
    	show_heatmap_legend = show_legend, show_annotation_legend = show_legend)

    for(k in all_k) {
    	decorate_heatmap_body(paste0("membership_", k), {
    		grid.rect(0, width = unit(1+1/k, "npc"), just = "left", gp = gpar(col = "black", fill = "transparent"))
    	})
    }
})


