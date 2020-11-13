
# == title
# Get class IDs from the HierarchicalPartition object
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
#
# == return
# A data frame of classes IDs. The class IDs are the node IDs where the subgroup sits in the hierarchy.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# get_classes(golub_cola_rh)
setMethod(f = "get_classes",
	signature = "HierarchicalPartition",
	definition = function(object, filter_node = filter_node_param()) {

	subgroup = object@subgroup
	all_leaves = all_leaves(object, filter_node = filter_node)
	# all_leaves should be parent node of subgroup
	map = rep(NA, length(unique(subgroup)))
	names(map) = unique(subgroup)
	for(nm in names(map)) {
		for(leaf in all_leaves) {
			if(grepl(qq("^@{leaf}"), nm)) {
				map[nm] = leaf
			}
		}
	}

	subgroup = map[subgroup]
	names(subgroup) = colnames(object)
	return(subgroup)
})


# == title
# Get signatures rows
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
# -group_diff Cutoff for the maximal difference between group means.
# -row_km Number of groups for performing k-means clustering on rows. By default it is automatically selected.
# -scale_rows whether apply row scaling when making the heatmap.
# -anno a data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `hierarchical_partition`.
# -anno_col a list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_column_names whether show column names in the heatmap.
# -column_names_gp Graphic parameters for column names.
# -verbose whether to print messages.
# -plot whether to make the plot.
# -seed Random seed.
# -... other arguments pass to `get_signatures,ConsensusPartition-method`.
# 
# == details
# The function calls `get_signatures,ConsensusPartition-method` to find signatures at
# each node of the partition hierarchy.
#
# == return 
# A data frame with more than two columns:
#
# -``which_row``: row index corresponding to the original matrix.
# -``km``: the k-means groups if ``row_km`` is set.
# -other_columns: the mean value (depending rows are scaled or not) in each subgroup.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# \donttest{
# data(golub_cola_rh)
# tb = get_signatures(golub_cola_rh)
# head(tb)
# }
setMethod(f = "get_signatures",
	signature = "HierarchicalPartition",
	definition = function(object, filter_node = filter_node_param(),
	group_diff = cola_opt$group_diff,
	row_km = NULL,
	scale_rows = object[1]@scale_rows, 
	anno = get_anno(object), 
	anno_col = get_anno_col(object),
	show_column_names = FALSE, column_names_gp = gpar(fontsize = 8),
	verbose = TRUE, plot = TRUE, seed = 888,
	...) {

	if(!has_hierarchy(object)) {
		cat("No hierarchy found.")
		return(invisible(NULL))
	}

	alf = all_leaves(object, filter_node)
	ap = setdiff(all_nodes(object, filter_node), alf)

	sig_lt = list()
	.env = object@list[[1]]@.env
	for(p in ap) {
		best_k = suggest_best_k(object[[p]])
		if(verbose) qqcat("* get signatures at node @{p} with @{best_k} subgroups.\n")
		sig_tb = get_signatures(object[[p]], k = best_k, prefix = "  ", verbose = TRUE, plot = FALSE, simplify = TRUE, seed = seed, ...)
		if(is.null(.env$signature_hash)) {
    		.env$signature_hash = list()
    	}
    	.env$signature_hash[[p]] = attr(sig_tb, "hash")
    	
		sig_lt[[p]] = sig_tb
		# if(verbose) qqcat("  * find @{nrow(sig_tb)} signatures at node @{p}\n")
	}

	all_index = sort(unique(unlist(lapply(sig_lt, function(x) x[, 1]))))

	returned_df = data.frame(which_row = all_index)

	# filter by group_diff
	mat = object@.env$data[all_index, , drop = FALSE]
	class = get_classes(object, filter_node)

	mat1 = mat
	if(nrow(mat) == 1) {
		group_mean = rbind(tapply(mat1, class, mean))
	} else {
		group_mean = do.call("cbind", tapply(seq_len(ncol(mat1)), class, function(ind) {
			rowMeans(mat1[, ind, drop = FALSE])
		}))
	}
	colnames(group_mean) = paste0("mean_", colnames(group_mean))
	returned_df = cbind(returned_df, group_mean)
	returned_df$group_diff = apply(group_mean, 1, function(x) max(x) - min(x))

	if(group_diff > 0) {
		l_diff = returned_df$group_diff >= group_diff
		mat = mat[l_diff, , drop = FALSE]
		mat1 = mat1[l_diff, , drop = FALSE]
		returned_df = returned_df[l_diff, , drop = FALSE]
	}

	if(scale_rows) {
		mat1_scaled = t(scale(t(mat)))
		if(nrow(mat) == 1) {
			group_mean_scaled = rbind(tapply(mat1_scaled, class, mean))
		} else {
			group_mean_scaled = do.call("cbind", tapply(seq_len(ncol(mat1_scaled)), class, function(ind) {
				rowMeans(mat1_scaled[, ind, drop = FALSE])
			}))
		}
		colnames(group_mean_scaled) = paste0("scaled_mean_", colnames(group_mean_scaled))
		returned_df = cbind(returned_df, group_mean_scaled)
	}

	returned_obj = returned_df
	rownames(returned_obj) = NULL

	## add k-means
	row_km_fit = NULL
	if(nrow(mat1) > 10) {
		if(scale_rows) {
			mat_for_km = t(scale(t(mat1)))
		} else {
			mat_for_km = mat1
		}

		if(nrow(mat_for_km) > 5000) {
			set.seed(seed)
			mat_for_km2 = mat_for_km[sample(nrow(mat_for_km), 5000), , drop = FALSE]
		} else {
			mat_for_km2 = mat_for_km
		}

		set.seed(seed)
		if(is.null(row_km)) {
			row_km = guess_best_km(mat_for_km2)
			if(length(unique(class)) == 1) row_km = 1
			if(length(unique(class)) == 2) row_km = min(row_km, 2)
		}
		if(row_km > 1) {
			row_km_fit = kmeans(mat_for_km2, centers = row_km)
			returned_obj$km = apply(pdist(row_km_fit$centers, mat_for_km, as.integer(1)), 2, which.min)
		}
		if(verbose) qqcat("* split rows into @{row_km} groups by k-means clustering.\n")
	}

	if(verbose) {
		qqcat("* found @{nrow(mat)} signatures (@{sprintf('%.1f',nrow(mat)/nrow(object)*100)}%).\n")
	}

	if(nrow(mat) == 0) {
		if(plot) {
			grid.newpage()
			fontsize = convertUnit(unit(0.1, "npc"), "char", valueOnly = TRUE)*get.gpar("fontsize")$fontsize
			grid.text("no sigatures", gp = gpar(fontsize = fontsize))
		}
		return(invisible(data.frame(which_row = integer(0))))
	}
	
	if(!plot) {
		return(invisible(returned_obj))
	}

	if(nrow(mat) > 2000) {
		set.seed(seed)
		if(verbose) qqcat("* randomly sample @{2000} rows from @{nrow(mat)} total rows.\n")
		row_index = sample(nrow(mat), 2000)
	} else {
		row_index = seq_len(nrow(mat))
	}
	mat1 = mat[row_index, , drop = FALSE]

	base_mean = rowMeans(mat1)
	if(nrow(mat) == 1) {
		group_mean = matrix(tapply(mat1, class, mean), nrow = 1)
	} else {
		group_mean = do.call("cbind", tapply(seq_len(ncol(mat1)), class, function(ind) {
			rowMeans(mat1[, ind, drop = FALSE])
		}))
	}
	rel_diff = (rowMaxs(group_mean) - rowMins(group_mean))/base_mean/2

	if(is.null(anno)) {
		bottom_anno1 = NULL
	} else {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
			}
		}

		if(is.null(anno_col)) {
			bottom_anno1 = HeatmapAnnotation(df = anno,
				show_annotation_name = TRUE, annotation_name_side = "right")
		} else {
			bottom_anno1 = HeatmapAnnotation(df = anno, col = anno_col,
				show_annotation_name = TRUE, annotation_name_side = "right")
		}
	}

	if(scale_rows) {
		scaled_mean = base_mean
		scaled_sd = rowSds(mat1)
		scaled_mat1 = t(scale(t(mat1)))

		use_mat1 = scaled_mat1
		mat_range = quantile(abs(scaled_mat1), 0.95, na.rm = TRUE)
		col_fun = colorRamp2(c(-mat_range, 0, mat_range), c("green", "white", "red"))
		heatmap_name = "z-score"
	} else {
		use_mat1 = mat1
		mat_range = quantile(mat1, c(0.05, 0.95))
		col_fun = colorRamp2(c(mat_range[1], mean(mat_range), mat_range[2]), c("blue", "white", "red"))
		heatmap_name = "expr"
	}

	row_split = factor(returned_obj$km[row_index], levels = sort(unique(returned_obj$km[row_index])))

	if(verbose) qqcat("* making heatmaps for signatures\n")

	ha1 = HeatmapAnnotation(
		Class = class,
		col = list(Class = object@subgroup_col))

	# dend = cluster_within_group(use_mat1, class)
	dend = calc_dend(object, filter_node, use_mat1)

	ht_list = Heatmap(use_mat1, top_annotation = ha1,
		name = heatmap_name, show_row_names = FALSE, 
		show_column_names = show_column_names, column_names_gp = column_names_gp,
		col = col_fun,
		use_raster = TRUE, row_split = row_split,
		show_row_dend = FALSE, cluster_columns = dend, column_split = length(unique(class)),
		column_title = qq("@{length(unique(class))} groups, @{length(all_index)} signatures"),
		bottom_annotation = bottom_anno1,
		row_title = {if(length(unique(row_split)) <= 1) NULL else qq("k-means with @{length(unique(row_split))} groups")})

	all_value_positive = !any(mat1 < 0)
 	if(scale_rows && all_value_positive) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm")) +
			Heatmap(rel_diff, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
				show_row_names = FALSE, name = "rel_diff", width = unit(5, "mm"))
	}

	draw(ht_list)
	return(invisible(returned_obj))
})


# == title
# Compare Signatures from Different Nodes
#
# == param
# -object A `HierarchicalPartition-class` object. 
# -depth Depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
# -method Method to visualize.
# -upset_max_comb_sets Maximal number of combination sets to show.
# -verbose Whether to print message.
# -... Other arguments passed to `get_signatures,HierarchicalPartition-method`.
#
# == details
# It plots an Euler diagram or a UpSet plot showing the overlap of signatures from different nodes.
# On each node, the number of subgroups is inferred by `suggest_best_k,ConsensusPartition-method`.
#
# == example
# data(golub_cola_rh)
# compare_signatures(golub_cola_rh)
setMethod(f = "compare_signatures",
	signature = "HierarchicalPartition",
	definition = function(object, filter_node = filter_node_param(),
	method = c("euler", "upset"), upset_max_comb_sets = 20,
	verbose = interactive(), ...) {

	if(!has_hierarchy(object)) {
		cat("No hierarchy found.")
		return(invisible(NULL))
	}

	nodes = setdiff(all_nodes(object, filter_node), all_leaves(object, filter_node))
	lt = object@list[nodes]

	sig_list = lapply(lt, function(x) {
		tb = get_signatures(x, k = suggest_best_k(x), verbose = verbose, ..., plot = FALSE)
		if(is.null(tb)) {
			return(integer(0))
		} else {
			return(tb$which_row)
		}
	})

	l = sapply(sig_list, length) > 0
	if(any(!l) && verbose) {
		qqcat("Following nodes have no signature found: \"@{paste(names(sig_list)[!l], collapse=', ')}\"\n")
	}
	sig_list = sig_list[l]

	if(missing(method)) {
		if(length(sig_list) <= 6) {
			method = "euler"
		} else {
			method = "upset"
		}
	} else {
		method = match.arg(method)[1]
	}

	if(method == "euler") {
		plot(eulerr::euler(sig_list), legend = TRUE, quantities = TRUE, main = "Signatures from different nodes")
	} else {
		m = make_comb_mat(sig_list)
		if(length(comb_size(m)) > upset_max_comb_sets) {
			m = m[order(comb_size(m), decreasing = TRUE)[1:upset_max_comb_sets]]
		} else {
			m = m[order(comb_size(m), decreasing = TRUE)]
		}
		draw(UpSet(m, column_title = "Signatures from different nodes"))
	}

})

# == title
# Collect classes from HierarchicalPartition object
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
# -show_row_names Whether to show the row names.
# -row_names_gp Graphic parameters for row names.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `hierarchical_partition`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
#
# == details
# The function plots the hierarchy of the classes.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# collect_classes(golub_cola_rh)
# collect_classes(golub_cola_rh, depth = 2)
setMethod(f = "collect_classes",
	signature = "HierarchicalPartition",
	definition = function(object, filter_node = filter_node_param(),
	show_row_names = FALSE, row_names_gp = gpar(fontsize = 8),
	anno = get_anno(object[1]), anno_col = get_anno_col(object[1])) {

	if(!has_hierarchy(object)) {
		cat("No hierarchy found.")
		return(invisible(NULL))
	}

	cl = get_classes(object, filter_node)
	dend = calc_dend(object, filter_node)

	ht_list = Heatmap(cl, name = "Class", col = object@subgroup_col, width = unit(5, "mm"),
		row_title_rot = 0, cluster_rows = dend, row_dend_width = unit(2, "cm"),
		row_split = length(unique(cl)), show_row_names = FALSE, row_title = NULL)
	if(!is.null(anno)) {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
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
	if(show_row_names) {
		ht_list = ht_list + rowAnnotation(rn = anno_text(colnames(object), gp = row_names_gp))
	}

	draw(ht_list)
})


# == title
# Test correspondance between predicted classes and known factors
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
# -known A vector or a data frame with known factors. By default it is the annotation table set in `hierarchical_partition`.
# -verbose Whether to print messages.
#
# == value
# A data frame with columns:
#
# - number of samples
# - p-values from the tests
# - number of classes
#
# The classifications are extracted for each depth.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# # golub_cola_rh already has known annotations, so test_to_known_factors()
# # can be directly applied
# test_to_known_factors(golub_cola_rh)
setMethod(f = "test_to_known_factors",
	signature = "HierarchicalPartition",
	definition = function(object, known = get_anno(object[1]),
	filter_node = filter_node_param(), verbose = FALSE) {

	if(!has_hierarchy(object)) {
		cat("No hierarchy found.")
		return(invisible(NULL))
	}

	if(!is.null(known)) {
		if(is.atomic(known)) {
			df = data.frame(known)
			colnames(df) = deparse(substitute(known))
			known = df
		}
	} else {
		stop_wrap("Known factors should be provided.")
	}

	class = get_classes(object, filter_node)
	m = test_between_factors(class, known, verbose = verbose)
	return(m)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
# -min_n_signatures Minimal number of signatures on the node.
# -top_n Top n rows to use. By default it uses all rows in the original matrix.
# -parent_node Parent node. If it is set, the function call is identical to ``dimension_reduction(object[parent_node])``
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -scale_rows Whether to perform scaling on matrix rows.
# -verbose Whether print messages.
# -... Other arguments passed to `dimension_reduction,ConsensusPartition-method`.
#
# == details
# The class IDs are extract at ``depth``.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# dimension_reduction(golub_cola_rh)
setMethod(f = "dimension_reduction",
	signature = "HierarchicalPartition",
	definition = function(object, filter_node = filter_node_param(),
	parent_node, top_n = NULL, top_value_method = object@list[[1]]@top_value_method,
	method = c("PCA", "MDS", "t-SNE", "UMAP"),
	scale_rows = TRUE, verbose = TRUE, ...) {

	if(!has_hierarchy(object)) {
		cat("No hierarchy found.")
		return(invisible(NULL))
	}

	cl = as.list(match.call())
	# default value
	if(! "method" %in% names(cl)) {
		method = NULL
		if(is.null(method)) {
			oe = try(loadNamespace("umap"), silent = TRUE)
			if(inherits(oe, "try-error")) {
				if(verbose) cat("umap package is not installed.\n")
			} else {
				if(verbose) cat("use UMAP\n")
				method = "UMAP"
			}
		}
		if(is.null(method)) {
			oe = try(loadNamespace("Rtsne"), silent = TRUE)
			if(inherits(oe, "try-error")) {
				if(verbose) cat("Rtsne package is not installed.\n")
			} else {
				if(verbose) cat("use t-SNE\n")
				method = "t-SNE"
			}
		}
		if(is.null(method)) {
			if(verbose) cat("use PCA\n")
			method = "PCA"
		}
	}

	method = match.arg(method)
	data = object@list[[1]]@.env$data

	op = par(c("mar", "xpd"))
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)

	if(missing(parent_node)) {
		if(!is.null(top_n)) {
			top_n = min(c(top_n, nrow(data)))
			if(top_value_method %in% names(object@.env$node_0_top_value_list)) {
				all_top_value = object@.env$node_0_top_value_list[top_value_method]
			} else {
				if(verbose) qqcat("calculating @{top_value_method} values.\n")
				all_top_value = get_top_value_method(top_value_method)(data)
			}
			ind = order(all_top_value, decreasing = TRUE)[1:top_n]
			data = data[ind, , drop = FALSE]
		} else {
			top_n = nrow(data)
		}
		class = get_classes(object, filter_node)
		n_class = length(unique(class))
		dimension_reduction(data, pch = 16, col = object@subgroup_col[class],
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{ncol(data)} samples with @{n_class} classes"),
			method = method, scale_rows = scale_rows, ...)
		class_level = sort(unique(class))
		legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(class_level, "ambiguous"), 
			pch = c(rep(16, n_class), 0),
			col = c(object@subgroup_col[class_level], "white"), xjust = 0, yjust = 0.5,
			title = "Class", title.adj = 0.1, bty = "n",
			text.col = c(rep("black", n_class), "white"))
	} else {
		if(!parent_node %in% setdiff(all_nodes(object), all_leaves(object))) {
			stop_wrap(qq("@{parent_node} has no children nodes."))
		}
		obj = object[parent_node]
		dimension_reduction(obj, k = suggest_best_k(obj), top_n = top_n, method = method,
			scale_rows = scale_rows, ...)
		legend(x = par("usr")[2], y = par("usr")[4], legend = qq("node @{parent_node}"))
	}

	par(op)
})



# == title
# Guess the best number of partitions
#
# == param
# -object A `HierarchicalPartition-class` object.
# -jaccard_index_cutoff The cutoff for Jaccard index for comparing to previous k.
#
# == details
# It basically gives the best k at each node.
#
# == value
# A data frame with the best k and other statistics for each node.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# suggest_best_k(golub_cola_rh)
setMethod(f = "suggest_best_k",
	signature = "HierarchicalPartition",
	definition = function(object, jaccard_index_cutoff = 0.95) {

	best_k = NULL
	stability = NULL
	mean_silhouette = NULL
	concordance = NULL
	n_sample = NULL
	method = NULL
	for(i in seq_along(object@list)) {
		obj = object@list[[i]]
		if(inherits(obj, "ConsensusPartition")) {
			k = suggest_best_k(obj, jaccard_index_cutoff)
			best_k[i] = k

			if(is.na(best_k[i])) {
				stability[i] = NA
				mean_silhouette[i] = NA
				concordance[i] = NA
			} else {
				stat = get_stats(obj, k = best_k[i])
				stability[i] = stat[1, "1-PAC"]
				mean_silhouette[i] = stat[1, "mean_silhouette"]
				concordance[i] = stat[1, "concordance"]
			}
			n_sample[i] = ncol(obj)
			method[i] = paste0(obj@top_value_method, ":", obj@partition_method)
		} else {
			best_k[i] = NA
			stability[i] = NA
			mean_silhouette[i] = NA
			concordance[i] = NA
			n_sample[i] = length(attr(obj, "column_index"))
			method[i] = "not applied"
		}
	}

	tb = data.frame(
		node = names(object@list),
		best_method = method,
		is_leaf = names(object@list) %in% all_leaves(object),
		best_k = best_k,
		"1-PAC" = stability,
		mean_silhouette = mean_silhouette,
		concordance = concordance,
		n_sample = n_sample,
		check.names = FALSE)

	rntb = rownames(tb)
	l = tb$`1-PAC` >= 0.9 & !is.na(tb$best_k)

	tb = cbind(tb, ifelse(l, ifelse(tb$`1-PAC` <= 0.95, "*", "**"), ""), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""
	
	stop_reason = lapply(object@list, function(obj) {
		attr(obj, "stop_reason")
	})
	attr(tb, "stop_reason") = stop_reason

	class(tb) = c("hc_table_suggest_best_k", class(tb))
	tb
})

# == title
# Print the hc_table_suggest_best_k object
#
# == param
# -x A ``hc_table_suggest_best_k`` object from `suggest_best_k,HierarchicalPartition-method`.
# -... Other arguments.
#
print.hc_table_suggest_best_k = function(x, ...) {
	stop_reason = attr(x, "stop_reason")
	stop_reason = sapply(stop_reason, function(x) {
		if(is.null(x)) {
			return(NA)
		} else {
			return(STOP_REASON_INDEX[x])
		}
	})

	x$is_leaf = ifelse(x$is_leaf, "\u2713", "")
	x$is_leaf = ifelse(is.na(stop_reason), x$is_leaf, paste0(x$is_leaf, "(", stop_reason, ")"))
	x$`1-PAC` = round(x$`1-PAC`, 3)
	x$mean_silhouette = round(x$mean_silhouette, 3)
	x$concordance = round(x$concordance, 3)
	print.data.frame(x, digits = 3, row.names = FALSE)
	cat(strrep("-", sum(sapply(colnames(x), nchar))+ ncol(x) - 4 + max(sapply(x$node, nchar))), "\n")

	if(any(!is.na(stop_reason))) {
		cat("Stop reason:\n")
	}
	for(a in sort(unique(stop_reason[!is.na(stop_reason)]))) {
		cat("  ", a, ") ", names(which(STOP_REASON_INDEX == a)), "\n", sep = "")
	}
}


# == title
# Make HTML report from the HierarchicalPartition object
#
# == param
# -object A `HierarchicalPartition-class` object.
# -output_dir The output directory where the report is put.
# -mc.cores Multiple cores to use.
# -title Title of the report.
# -env Where the objects in the report are found, internally used.
#
# == details
# This function generates a HTML report which contains all plots for all nodes
# in the partition hierarchy.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# if(FALSE) {
# # the following code is runnable
# data(golub_cola_rh)
# cola_report(golub_cola_rh, output_dir = "~/test_cola_rh_report")
# }
setMethod(f = "cola_report",
	signature = "HierarchicalPartition",
	definition = function(object, output_dir = getwd(), mc.cores = 1,
	title = qq("cola Report for Hierarchical Partitioning"), 
	env = parent.frame()) {

	if(max_depth(object) <= 1) {
		cat("No hierarchy is detected, no report is generated.\n")

		dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
		output_dir = normalizePath(output_dir, mustWork = FALSE)

		qqcat("<html><head><title>@{title}</title></head><body><p>No hierarchy is detected, no report is generated.</p></body></html>", file = qq("@{output_dir}/cola_hc.html"))
	
		return(invisible(NULL))
	}

	check_pkg("genefilter", bioc = TRUE)

	var_name = deparse(substitute(object, env = env))
	make_report(var_name, object, output_dir, mc.cores = mc.cores, title = title, class = "HierarchicalPartition")
})


# == title
# Heatmap of top rows from different top-value methods
#
# == param
# -object A `HierarchicalPartition-class` object.
# -top_n Number of top rows.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `hierarchical_partition`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -scale_rows Wether to scale rows. 
# -... Pass to `top_rows_heatmap,matrix-method`
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
	signature = "HierarchicalPartition",
	definition = function(object, top_n = min(object@list[[1]]@top_n), 
	anno = get_anno(object), anno_col = get_anno_col(object),
	scale_rows = object@list[[1]]@scale_rows, ...) {

	all_top_value_list = object@.env$node_0_top_value_list
    
    mat = object@.env$data

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

    top_rows_heatmap(mat, all_top_value_list = all_top_value_list, top_n = top_n, 
    	scale_rows = scale_rows, bottom_annotation = bottom_anno, ...)
})

