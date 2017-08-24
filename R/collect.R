
# == title
# Collect plots from ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object from `run_all_consensus_partition_methods`.
# -k number of partitions.
# -fun function used to generate plots. Valid functions are `consensus_heatmap,ConsensusPartition-method`,
#        `plot_ecdf,ConsensusPartition-method`, `membership_heatmap,ConsensusPartition-method`,
#        `get_signatures,ConsensusPartition-method` and `dimension_reduction,ConsensusPartition-method`.
# -top_method a vector of top methods.
# -partition_method a vector of partition methods.
# -... other arguments passed to corresponding ``fun``.
#
# == details
# Plots for all combinations of top methods and parittion methods are arranged in one page.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
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
	
	highlight_row = NULL
	highlight_col = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(object, top_method = top_method[i], partition_method = partition_method[j])
	    	if(get_stat(res, k = k)[, "PAC"] < 0.1) {
	    		highlight_row = c(highlight_row, i + 1)
	    		highlight_col = c(highlight_col, j + 1)
	    	}
	    	pushViewport(viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
	    	# image_width = convertWidth(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		# image_height = convertHeight(unit(1, "npc"), "bigpts", valueOnly = TRUE)
    		image_width = 800
    		image_height = 800
	        file_name = tempfile(fileext = ".png")
	        png(file_name, width = image_width, height = image_height)
	        oe = try(fun(res, k = k, show_legend = FALSE, show_column_names = FALSE, use_raster = FALSE, ...))
	        dev.off()
	        if(!inherits(oe, "try-error")) {
		        grid.raster(readPNG(file_name))
		    } else {
		    	qqcat("Caught an error for @{top_method[i]}:@{partition_method[j]}\n")
		    }
		    grid.rect(gp = gpar(fill = "transparent", col = "black"))
		    upViewport()
			if(file.exists(file_name)) file.remove(file_name)
	    }
	}
	for(i in seq_along(highlight_row)) {
		pushViewport(viewport(layout.pos.row = highlight_row[i], layout.pos.col = highlight_col[i]))
		grid.rect(gp = gpar(fill = "transparent", col = "red", lwd = 2))
		upViewport()
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
# and `get_signatures,ConsensusPartition-method` are arranged in one single page, for all avaialble k.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
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
    if(file.exists(file_name)) file.remove(file_name)
	
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
    if(file.exists(file_name)) file.remove(file_name)

    pac = get_stat(object, k = all_k)[, "PAC"]
    border_color = ifelse(pac < 0.1, "red", "black")

	for(i in seq_along(all_k)) {
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
		file_name = tempfile()
        png(file_name, width = image_width*2, height = image_height*2)
        oe = try(consensus_heatmap(object, k = all_k[i], show_legend = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))  
	    }
	    grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()
	    if(file.exists(file_name)) file.remove(file_name)

	    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
	    file_name = tempfile()
        png(file_name, width = image_width*2, height=  image_height*2)
        oe = try(membership_heatmap(object, k = all_k[i], show_legend = FALSE, show_column_names = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))
	    }
	    grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()
	    if(file.exists(file_name)) file.remove(file_name)

	    pushViewport(viewport(layout.pos.row = 4, layout.pos.col = i))
	    file_name = tempfile()
        png(file_name, width = image_width*2, height=  image_height*2)
        oe = try(get_signatures(object, k = all_k[i], show_legend = FALSE, show_column_names = FALSE, use_raster = FALSE, ...))
        dev.off()
        if(!inherits(oe, "try-error")) {
	        grid.raster(readPNG(file_name))  
	    }
	    grid.rect(gp = gpar(fill = "transparent", col = border_color[i]))
	    upViewport()
	    if(file.exists(file_name)) file.remove(file_name)
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
# == details
# The function plots class IDs for all methods and the similarity between methods.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "collect_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k, 
	top_method = object@top_method, partition_method = object@partition_method, ...) {

	top_method_vec = NULL
	partition_method_vec = NULL
	class_mat = NULL
	silhouette_mat = NULL
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {  
	    	res = get_single_run(object, top_method = top_method[i], partition_method = partition_method[j])

	        top_method_vec = c(top_method_vec, top_method[i])
	        partition_method_vec = c(partition_method_vec, partition_method[j])
	        class_df = get_class(res, k)
	        class_mat = cbind(class_mat, class_df[, "class"])
	        silhouette_mat = cbind(silhouette_mat, class_df[, "silhouette"])
	    }
	}

	class_mat = as.matrix(class_mat)
	colnames(class_mat) = paste(top_method_vec, partition_method_vec, sep = ":")
	ik = which(res@k == k)
	
	silhouette_mat = as.matrix(silhouette_mat)
	silhouette_mat[silhouette_mat < 0] = 0

	adjust_by_transparency = function(col, transparency) {
		rgb( 1 - (1 - t(col2rgb(col)/255)) * (1 - transparency))
	}

	pac = get_stat(object, k)[, "PAC"][colnames(class_mat)]
	consensus_class = get_class(object, k = k)$class
	m = t(class_mat)
	column_order = column_order_by_group(consensus_class, m)
	if(is.null(object@list[[1]]@known_anno)) {
		bottom_annotation = NULL
	} else {
		bottom_annotation = HeatmapAnnotation(df = object@list[[1]]@known_anno, 
			col = object@list[[1]]@known_col, show_annotation_name = TRUE,
			annotation_name_side = "left")
	}
	ht = Heatmap(m, name = "class", col = brewer_pal_set2_col, column_order = column_order,
		row_names_side = "left", show_row_dend = FALSE, cluster_columns = FALSE,
		show_column_dend = FALSE, rect_gp = gpar(type = "none"),
		cell_fun = function(j, i, x, y, w, h, fill) {
			col = adjust_by_transparency(fill, 1 - silhouette_mat[j, i])
			grid.rect(x, y, w, h, gp = gpar(fill = col, col = col))
		}, top_annotation = HeatmapAnnotation(consensus_class = consensus_class, 
			col = list(consensus_class = brewer_pal_set2_col),
			show_annotation_name = TRUE, annotation_name_side = "left", show_legend = FALSE),
		bottom_annotation = bottom_annotation) + 
		rowAnnotation(mean_silhouette = row_anno_barplot(get_stat(object, k = k)[, "mean_silhouette"], baseline = 0, axis = TRUE),
			width = unit(2, "cm")) +
		Heatmap((pac < 0.1) + 0, name = "PAC < 0.1", col = c("1" = "orange", "0" = "white"),
			show_row_names = FALSE, width = unit(2, "mm")) +
		rowAnnotation(top_method = top_method_vec, 
			partition_method = partition_method_vec,
			show_annotation_name = TRUE, annotation_name_side = "bottom",
			col = list(top_method = structure(names = top_method, brewer_pal_set1_col[seq_along(top_method)]),
			           partition_method = structure(names = partition_method, brewer_pal_set2_col[seq_along(partition_method)])),
			width = unit(10, "mm"))
	
	pl = lapply(object@list[paste(top_method_vec, partition_method_vec, sep = ":")], function(x) as.cl_partition(get_membership(x, k = k)))
	clen = cl_ensemble(list = pl)
	m_diss = cl_dissimilarity(clen, method = "comembership")
	m_diss = as.matrix(m_diss)

	ht = ht + Heatmap(m_diss, name = "dissimilarity", show_row_names = FALSE, show_column_names = FALSE,
		show_row_dend = FALSE, show_column_dend = FALSE, width = unit(6, "cm"),
		col = colorRamp2(quantile(m_diss, c(0, 0.5, 1)), c("red", "white", "blue")))

	draw(ht, main_heatmap = "dissimilarity", column_title = qq("classification from all methods, k = @{k}"))
	decorate_annotation("mean_silhouette", {
		grid.text("mean_silhouette", y = unit(-1, "cm"))
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
# == details
# Membership matrix and the classes with each k are plotted in the heatmap.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
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

		ht_list = ht_list + Heatmap(membership, col = colorRamp2(c(0, 1), c("white", "red")),
			show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, show_heatmap_legend = i == 1,
			name = paste0("membership_", all_k[i])) + 
			Heatmap(class, col = brewer_pal_set2_col, 
				show_row_names = FALSE, show_heatmap_legend = i == length(all_k), name = paste(all_k[i], "_classes"))
		gap = c(gap, c(0, 4))
		class_mat = cbind(class_mat, as.numeric(class))
	}

    draw(ht_list,gap = unit(gap, "mm"), row_order = do.call("order", as.data.frame(class_mat)),
    	column_title = qq("classes from k = '@{paste(all_k, collapse = ', ')}'"),
    	show_heatmap_legend = show_legend, show_annotation_legend = show_legend)

    for(k in all_k) {
    	ik = which(all_k == k )
    	border_color = ifelse(object@object_list[[ik]]$stat$PAC < 0.1, "red", "black")
    	decorate_heatmap_body(paste0("membership_", k), {
    		grid.rect(0, width = unit(1+1/k, "npc"), just = "left", gp = gpar(col = border_color, fill = "transparent"))
    	})
    }
})


