
# == title
# Consensus partition
#
# == param
# -data a numeric matrix where subgroups are found by columns.
# -top_method a single top method. Avaialble methods are in `all_top_value_methods`.
# -top_n number of rows with top values.
# -partition_method a single partition method. Avaialble methods are in `all_partition_methods`.
# -k number of partitions. The value is a vector.
# -p_sampling proportion of the top n rows to sample.
# -partition_repeat number of repeats for the random sampling.
# -partition_param parameters for the partition method.
# -known_anno a data frame with known annotation of columns.
# -known_col a list of colors for the annotations in ``known_anno``.
# -column_index a subst of columns in the matrix, only for internal use.
# -.env an environment, internally used.
#
# == return
# A `ConsensusPartition-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
consensus_partition = function(data,
	top_method = all_top_value_methods()[1],
	top_n = seq(min(2000, round(nrow(data)*0.2)), min(c(6000, round(nrow(data)*0.6))), length.out = 5),
	partition_method = all_partition_methods()[1],
	k = 2:6, p_sampling = 0.8,
	partition_repeat = 50,
	partition_param = list(),
	known_anno = NULL,
	known_col = NULL,
	column_index = seq_len(ncol(data)),
	.env) {

	if(missing(.env)) {
		if(is.data.frame(data)) data = as.matrix(data)
		if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		.env = new.env()
		.env$data = data
		.env$column_index = column_index
	} else if(is.null(.env$data)) {
		if(is.data.frame(data)) data = as.matrix(data)
		if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		.env$data = data
		.env$column_index = column_index
	} else {
		data = .env$data
	}

	data = data[, .env$column_index, drop = FALSE]

	k = sort(k)
	k = k[k <= ncol(data)]
	
	top_n = round(top_n)
	top_n = top_n[top_n < nrow(data)]

	partition_fun = get_partition_fun(partition_method, partition_param)

	get_value_fun = get_top_value_fun(top_method)

	param = data.frame(top_n = numeric(0), k = numeric(0), n_row = numeric(0))
	partition_list = list()
	
	# .env is already defined
	if(is.null(.env$all_value_list)) {
		all_value = get_value_fun(data)
		all_value[is.na(all_value)] = -Inf
		.env$all_value_list = list()
		.env$all_value_list[[top_method]] = all_value
	} else if(is.null(.env$all_value_list[[top_method]])) {
		all_value = get_value_fun(data)
		all_value[is.na(all_value)] = -Inf
		.env$all_value_list[[top_method]] = all_value
	} else {
		all_value = .env$all_value_list[[top_method]]
	}

	scale_method = attr(partition_fun, "scale_method")
	if("standardization" %in% scale_method) {
		cat("rows are scaled before sent to partition: standardization\n")
		data = t(scale(t(data)))
	} else if("rescale" %in% scale_method) {
		cat("rows are scaled before sent to partition: rescale\n")
		row_min = rowMeans(data)
		row_max = rowMaxs(data)
		row_range = row_max - row_min
		data = apply(data, 2, function(x) (x - row_min)/row_range)
	}

	for(i in seq_along(top_n)) {
		qqcat("get top @{top_n[i]} rows by @{top_method} method\n")
		ind = order(all_value, decreasing = TRUE)[1:top_n[i]]

		if(length(ind) > 5000) {
			ind = sample(ind, 5000)
		}

		for(j in 1:partition_repeat) {
			ind_sub = sample(ind, round(p_sampling*length(ind)))
			mat = data[ind_sub, , drop = FALSE]
			for(y in k) {
				if(interactive()) cat(strrep("\b", 100))
				if(interactive()) qqcat("  [k = @{y}] @{partition_method} repeated for @{j}th sampling from top @{top_n[i]} rows.")
				partition_list = c(partition_list, list(list(partition_fun(mat, y))))
				param = rbind(param, data.frame(top_n = top_n[i], k = y, n_row = nrow(mat)))
			}
		}
		if(interactive()) cat("\n")
	}

	construct_consensus_object = function(param, partition_list, k, prefix = "  ") {

		partition_list = do.call("c", partition_list)
		partition_list = cl_ensemble(list = partition_list)

		qqcat("@{prefix}merging @{length(partition_list)} partitions into a single ensemble object.\n")
		partition_consensus = cl_consensus(partition_list)

		class_ids = as.vector(cl_class_ids(partition_consensus))

		membership_mat = cl_membership(partition_consensus)

		class(membership_mat) = "matrix"
		colnames(membership_mat) = paste0("p", 1:ncol(membership_mat))
		attr(membership_mat, "n_of_classes") = NULL
		attr(membership_mat, "is_cl_hard_partition") = NULL

		qqcat("@{prefix}calculate consensus matrix for samples clustered in a same group.\n")
		membership_each = do.call("cbind", lapply(partition_list, function(x) {
			class = as.vector(cl_class_ids(x))
			map = relabel_class(class, class_ids)
			class = as.numeric(map[as.character(class)])
			class
		}))
		rownames(membership_each) = rownames(membership_mat)

		consensus_mat = matrix(1, nrow = nrow(membership_mat), ncol = nrow(membership_mat))
		for(i in 1:(nrow(membership_each)-1)) {
			for(j in (i+1):nrow(membership_each)) {
				consensus_mat[i, j] = sum(membership_each[i, ] == membership_each[j, ])/ncol(membership_each)
				consensus_mat[j, i] = consensus_mat[i, j]
			}
	 	}
	 	rownames(consensus_mat) = rownames(membership_mat)
	 	colnames(consensus_mat) = rownames(membership_mat)

	 	class_df = data.frame(
	 		class = class_ids,
	 		entropy = apply(membership_mat, 1, entropy),
	 		stringsAsFactors = FALSE
	 	)
	 	rownames(class_df) = colnames(data)

	 	if(length(unique(class_ids)) == 1) {
	 		class_df$silhouette = rep(0, length(class_ids))
	 	} else {
			class_df$silhouette = silhouette(class_ids, dist(t(consensus_mat)))[, "sil_width"]
		}

		ind = order(all_value, decreasing = TRUE)[1:max(top_n)]
		stat = list(
			ecdf = ecdf(consensus_mat[lower.tri(consensus_mat)]),
			cophcor =  cophcor(consensus_mat),
			PAC = PAC(consensus_mat),
			mean_silhouette = mean(class_df$silhouette),
			tot_withinss = tot_withinss(class_df$class, data[ind, , drop = FALSE])
		)
		
		return(list(
			class_df = class_df, 
			membership = membership_mat, 
			consensus = consensus_mat, 
			param = param, 
			membership_each = membership_each,
			stat = stat
		))
	}

	object_list = lapply(k, function(y) {
		l = param$k == y
		qqcat("wrapping results for k = @{y}\n")
		construct_consensus_object(param[l, ], partition_list[l], y)
	})
	names(object_list) = as.character(k)

	rm(partition_list)
	gc(verbose = FALSE)

	## adjust class labels
	reference_class = object_list[[1]]$class_df$class
	for(i in seq_along(k)[-1]) {
		class_df = object_list[[i]]$class_df
    	class = class_df[, "class"]

    	map = relabel_class(class, reference_class)
    	map2 = structure(names(map), names = map)
    	object_list[[i]]$class_df$class = as.numeric(map[as.character(class)])
    	
    	object_list[[i]]$membership = object_list[[i]]$membership[, as.numeric(map2[as.character(1:k[i])]) ]
		colnames(object_list[[i]]$membership) = paste0("p", 1:k[i])
		
		odim = dim(object_list[[i]]$membership_each)
		object_list[[i]]$membership_each = as.numeric(map[as.character(object_list[[i]]$membership_each)])
		dim(object_list[[i]]$membership_each) = odim

		reference_class = object_list[[i]]$class_df$class
	}

	ak = sapply(object_list, function(obj) {
		f = obj$stat$ecdf
		x = seq(0, 1, length = 100)
		n = length(x)
		sum((x[2:n] - x[1:(n-1)])*f(x[2:n]))
	})
	delta_k = ak
	for(i in seq_along(k)[-1]) {
		delta_k[i] = (ak[i] - ak[i-1])/ak[i-1]
	}
	for(i in seq_along(object_list)) {
		object_list[[i]]$stat$area_increased = delta_k[i]
	}

	if(!is.null(known_anno)) {
		if(is.atomic(known_anno)) {
			known_nm = deparse(substitute(known_anno))
			known_anno = data.frame(known_anno)
			colnames(known_anno) = known_nm
			if(!is.null(known_col)) {
				known_col = list(known_col)
				names(known_col) = known_nm
			}
		}
	}

	res = ConsensusPartition(object_list = object_list, k = k, n_partition = partition_repeat * length(top_n),  
		partition_method = partition_method, top_method = top_method, top_n = top_n,
		known_anno = known_anno, known_col = known_col, .env = .env)
	
	return(res)
}

# == title
# Print the consensus_partition object
#
# == param
# -object a `ConsensusPartition-class` object
# -... other arguments
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "show",
	signature = "ConsensusPartition",
	definition = function(object) {
	qqcat("A 'ConsensusPartition' object with k = @{paste(object@k, collapse = ', ')}.\n")
	qqcat("  top rows (@{paste(object@top_n, collapse = ', ')}) are extracted by '@{object@top_method}' method.\n")
	qqcat("  subgroups are detected by '@{object@partition_method}' method.\n")
})

# == title
# Plot the ecdf of the consensus matrix
#
# == param
# -object a `ConsensusPartition-class` object.
# -... other arguments.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "plot_ecdf",
	signature = "ConsensusPartition",
	definition = function(object, ...) {
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "P(X >= x)")
	for(i in seq_along(object@k)) {
		consensus_mat = get_consensus(object, k = object@k[i])
		f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		x = seq(0, 1, length = 100)
		lines(x, f(x), col = i, lwd = 4)
	}
	legend("bottomright", pch = 15, legend = paste0("k = ", object@k), col = seq_along(object@k))
})

# == title
# Several plots for determine the optimized number of partitions
#
# == param
# -object a `ConsensusPartition-class` object
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "select_partition_number",
	signature = "ConsensusPartition",
	definition = function(object) {
	op = par(no.readonly = TRUE)

	m = get_stat(object)
	nm = colnames(m)

	par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))

	plot_ecdf(object)

	for(i in seq_len(ncol(m))) {
		plot(object@k, m[, i], type = "b", xlab = "k", ylab = nm[i])
	}

	par(op)
})

# == title
# Get the best number of partitions
#
# == param
# -obje# -object a `ConsensusPartition-class` object
#
# == value
# The best k
#
# == author
# Zuguang Gu
setMethod(f = "get_best_k",
	signature = "ConsensusPartition",
	definition = function(object) {

	stat = get_stat(object)
	dec = c(which.max(stat[, "cophcor"]), 
		    which.min(stat[, "PAC"]),
		    which.max(stat[, "mean_silhouette"]))
	if(length(unique(dec)) == 3) {
		x = rownames(stat)[dec[1]]
	} else {
		tb = table(dec)
		x = rownames(stat)[as.numeric(names(tb[which.max(tb)]))]
	}
	as.numeric(x)
})

# == title
# Heatmap for the consensus matrix
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
# -show_legend whether show heatmap and annotation legends.
# -anno a data frame with column annotations
# -anno_col colors for the annotations
# -show_row_names whether plot row names on the consensus heatmap
# -... other arguments
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "consensus_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, show_legend = TRUE,
	anno = object@known_anno, 
	anno_col = if(missing(anno)) object@known_col else NULL, 
	show_row_names = FALSE, ...) {

	class_df = get_class(object, k)
	class_ids = class_df$class

	consensus_mat = get_consensus(object, k)

	mat_col_od = column_order_by_group(class_ids, consensus_mat)

	membership_mat = get_membership(object, k)

	ht_list = Heatmap(membership_mat, name = "membership", cluster_columns = FALSE, show_row_names = FALSE, 
		width = unit(5, "mm")*k, col = colorRamp2(c(0, 1), c("white", "red"))) + 
	Heatmap(class_df$silhouette, name = "silhouette", width = unit(5, "mm"), 
		show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "purple"))) +
	Heatmap(class_ids, name = "class", col = brewer_pal_set2_col, 
		show_row_names = FALSE, width = unit(5, "mm"))
	
	ht_list = ht_list +	Heatmap(consensus_mat, name = "consensus", show_row_names = show_row_names, show_row_dend = FALSE,
		col = colorRamp2(c(0, 1), c("white", "blue")), row_order = mat_col_od, column_order = mat_col_od,
		cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE)

	if(!is.null(anno)) {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		}
		if(is.null(anno_col))
			ht_list = ht_list + rowAnnotation(df = anno, show_annotation_name = TRUE,
				annotation_name_side = "bottom", width = unit(ncol(anno)*5, "mm"))
		else {
			ht_list = ht_list + rowAnnotation(df = anno, col = anno_col, show_annotation_name = TRUE,
				annotation_name_side = "bottom", width = unit(ncol(anno)*5, "mm"))
		}
	}
	draw(ht_list, main_heatmap = "consensus", column_title = qq("consensus @{object@partition_method} with @{k} groups from @{object@n_partition} partitions"),
		show_heatmap_legend = show_legend, show_annotation_legend = show_legend)
})

# == title
# Heatmap of membership of columns in each random sampling
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions.
# -show_legend whether show heatmap and annotation legends.
# -anno a data frame with column annotations
# -anno_col colors for the annotations
# -show_column_names whether show column names in the heatmap
# -... other arguments
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#`
setMethod(f = "membership_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, show_legend = TRUE, 
	anno = object@known_anno, 
	anno_col = if(missing(anno)) object@known_col else NULL,
	show_column_names = TRUE, ...) {

	class_df = get_class(object, k)
	class_ids = class_df$class

	membership_mat = get_membership(object, k)
	col_list = rep(list(colorRamp2(c(0, 1), c("white", "red"))), k)
	names(col_list) = colnames(membership_mat)

	membership_each = get_membership(object, k, each = TRUE)
	membership_each = t(membership_each)
	mat_col_od = column_order_by_group(class_ids, membership_each)

	col = brewer_pal_set1_col[1:k]

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

	param = get_param(object, k, unique = FALSE)

	n_row_level = unique(param$n_row)
	n_row_col = structure(brewer_pal_set2_col[seq_along(n_row_level)], names = n_row_level)
	
	ht = Heatmap(membership_each, name = "cluster", show_row_dend = FALSE, show_column_dend = FALSE, col = col,
		column_title = qq("membership heatmap, k = @{k}"), column_order = mat_col_od, cluster_columns = FALSE,
		split = param$n_row,
		top_annotation = HeatmapAnnotation(df = as.data.frame(membership_mat),
			class = class_ids, col = c(list(class = brewer_pal_set2_col), col_list),
			show_annotation_name = TRUE, annotation_name_side = "right",
			show_legend = c(TRUE, rep(FALSE, k - 1), TRUE)),
		bottom_annotation = bottom_anno,
		combined_name_fun = function(x) paste0(x, " rows"),
		show_column_names = show_column_names
		) + 
	Heatmap(as.character(param$n_row), name = "n_row", col = n_row_col,
		width = unit(5, "mm"), show_row_names = FALSE)

	draw(ht, row_title = qq("@{round(object@n_partition/length(n_row_level))} x @{length(n_row_level)} random samplings"),
		show_heatmap_legend = show_legend, show_annotation_legend = show_legend)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
# -top_n top n rows to use.
# -method which method to reduce the dimension of the data. ``mds`` uses `stats::cmdscale`,
#         ``pca`` uses ``stats::prcomp`` and ``tsne`` uses `Rtsne::Rtsne`.
# -silhouette_cutoff cutoff of silhouette. Data points with values less
#        than it will be mapped to small points.
# -remove whether to remove columns which have less silhouette values than
#        the cutoff.
# -tsne_param parameters pass to `Rtsne::Rtsne`
# -... other arguments
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "dimension_reduction",
	signature = "ConsensusPartition",
	definition = function(object, k, top_n = NULL,
	method = c("mds", "pca", "tsne"),
	silhouette_cutoff = 0.5, remove = FALSE,
	tsne_param = list(), ...) {

	method = match.arg(method)
	data = object@.env$data

	if(!is.null(top_n)) {
		top_n = min(c(top_n, nrow(data)))
		all_value = object@.env$all_value_list[[object@top_method]]
		ind = order(all_value)[1:top_n]
		data = data[ind, , drop = FALSE]
	}

	data = t(scale(t(data)))

	ik = which(object@k == k)

	class_df = get_class(object, k)
	class_ids = class_df$class

	l = class_df$silhouette >= silhouette_cutoff
	if(method == "mds") {
		loc = cmdscale(dist(t(data)))
	} else if(method == "pca") {
		loc = prcomp(t(data))$x[, 1:2]
	} else if(method == "tsne") {
		loc = do.call("Rtsne", c(list(X = as.matrix(t(data))), tsne_param))$Y
	}

	colnames(loc) = c("P1", "P2")
	loc = as.data.frame(loc)

	if(remove) {
		plot(loc[l, ], pch = 16, col = brewer_pal_set2_col[as.character(class_df$class[l])],
			cex = 1)
	} else {
		plot(loc, pch = ifelse(l, 16, 4), col = brewer_pal_set2_col[as.character(class_df$class)],
			cex = ifelse(l, 1, 0.7))
	}
	title(qq("Dimension reduction by @{method}, @{sum(l)}/@{length(l)} samples"))
})
