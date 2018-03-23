
# == title
# The HierarchicalPartition class
#
# == alias
# HierarchicalPartition
#
# == methods
# The `HierarchicalPartition-class` has following methods:
#
# -`hierarchical_partition`: constructor method.
# -`collect_classes,HierarchicalPartition-method`: plot the hierarchy of subgroups predicted.
# -`get_class,HierarchicalPartition-method`: get the class IDs of subgroups.
# -`get_signatures,HierarchicalPartition-method`: get the signatures for each subgroup.
# -`get_single_run,HierarchicalPartition-method`: get a `ConsensusPartition-class` object at a specified hierarchy level.
# -`test_to_known_factors,HierarchicalPartition-method`: test correlation between predicted subgrouping and known annotations, if available.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
HierarchicalPartition = setClass("HierarchicalPartition",
    slots = list(
        list = "list",
        hierarchy = "matrix",
        subgroup = "character",
        subgroup_col = "character",
        calling = "ANY"
    )
)

# == title
# Hierarchical detection of subgroups
#
# == param
# -data a numeric matrix where subgroups are found by samples.
# -top_method a single top method. Available methods are in `all_top_value_methods`.
# -partition_method a single partition method. Available ialble methods are in `all_partition_methods`.
# -PAC_cutoff the cutoff of PAC scores to determine whether to continuou looking to subgroups.
# -min_samples the cutoff of number of samples to determine whether to continuou looking to subgroups.
# -k a list number of partitions.
# -... pass to `consensus_partition`
#
# == details
# The function looks for subgroups in a hierarchical way.
#
# == return
# A `HierarchicalPartition-class` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
hierarchical_partition = function(data, top_method = "MAD", partition_method = "kmeans",
	PAC_cutoff = 0.2, min_samples = 6, k = 2:4, ...) {

	cl = match.call()
	
	.h_obj = new.env()
	.h_obj$hierarchy = matrix(nrow = 0, ncol = 2)
	.h_obj$list = list()
	.h_obj$subgroup = rep("", ncol(data))

	.hierarchical_partition = function(.env, column_index, PAC_cutoff = 0.2, node_id = '0', parent = NULL, min_samples = 6, k = 2:4, ...) {

		if(node_id != "0") {
			cat("\n")
			.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(parent, node_id))
		}
		qqcat("submatrix with @{length(column_index)} columns and @{nrow(data)} rows, node_id: @{node_id}.\n")
	   	.env$all_value_list = NULL
	   	.env$column_index = column_index
		part = consensus_partition(verbose = FALSE, .env = .env, k = k, 
			top_method = top_method, partition_method = partition_method, ...)

	    .h_obj$list[[node_id]] = part
	    .h_obj$subgroup[column_index] = node_id

	    best_k = get_best_k(part)
	    cl = get_class(part, k = best_k)

		oe = try(sig <- get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE))
		if(inherits(oe, "try-error")) {
			cat("caught an error, stop.\n")
			return(NULL)
		}
		if(length(sig$group) < 100) {
			qqcat("Number of signatures are too small (< 100), stop.\n")
			return(NULL)
		}

	    mat = .env$data[, column_index, drop = FALSE]
	    mat = t(scale(t(mat)))
	    mat = mat[order(.env$all_value_list[[1]], decreasing = TRUE)[1:max(part@top_n)], , drop = FALSE]
	    s = most_tight_subgroup(mat, cl$class)
	    set1 = which(cl$class == s)
	    set2 = which(cl$class != s)

	    qqcat("best k = @{best_k}, most tight subgroup: @{s}\n")
		PAC_score = get_stat(part, k = best_k)[, "PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	qqcat("PAC score too big @{PAC_score}, stop.\n")
	    	return(NULL)
	    }

	    # dist_decrease = mean_dist_decrease(mat, set1, set2)
	    # if(dist_decrease < reduce) {
	    # 	qqcat("mean distance does not decrease too much @{dist_decrease}, stop.\n")
	    # 	return(NULL)
	    # } else {

	    	if(length(set1) <= min_samples || length(set2) <= min_samples) {
	    		cat("subgroups have too few columns, stop.\n")
	    		return(NULL)
	    	}

	    	qqcat("partitioned into two subgroups with @{length(set1)} and @{length(set2)} columns.\n")
	    	# insert the two subgroups into the hierarchy
	    	sub_node_1 = paste0(node_id, s)
	    	sub_node_2 = paste0(node_id, "0")

	    	if(length(set1) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set1], node_id = sub_node_1, parent = node_id,
	    			PAC_cutoff = PAC_cutoff, min_samples = min_samples, k = k, ...)
	    	}

	    	if(length(set2) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set2], node_id = sub_node_2, parent = node_id,
	    			PAC_cutoff = PAC_cutoff, min_samples = min_samples, k = k, ...)
	    	}
	    	
	    # }

	    return(NULL)
	}

	.env = new.env()
	.env$data = data
	.hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), PAC_cutoff = PAC_cutoff, min_samples = min_samples, node_id = "0", parent = NULL, k = k, ...)

	hp = new("HierarchicalPartition")
	hp@hierarchy = .h_obj$hierarchy
	hp@list = .h_obj$list
	hp@subgroup = .h_obj$subgroup
	names(hp@subgroup) = colnames(data)

	le = unique(as.vector(hp@hierarchy))
	hp@subgroup_col = structure(rand_color(length(le), luminosity = "bright"), names = le)
	hp@calling = cl

	return(hp)
}

calc_dend = function(object, hierarchy = object@hierarchy) {

	lt = list()
	lt[["0"]] = Node$new("all samples")
	cn = colnames(object@list[["0"]]@.env$data)
	max_depth = max(nchar(hierarchy))
	lt[["0"]]$node_height = max_depth + 1
	for(i in seq_len(nrow(hierarchy))) {
		lt[[ hierarchy[i, 2] ]] = lt[[ hierarchy[i, 1] ]]$AddChildNode({
			node = Node$new(hierarchy[i, 2])
			node$node_height = max_depth - nchar(hierarchy[i, 2]) + 1
			node
		})
		l = hierarchy[, 1] == hierarchy[i, 2]
		if(sum(l) == 0) {
			cn2 = cn[object@list[[ hierarchy[i, 2] ]]@column_index]

			sample_node = lt[[ hierarchy[i, 2] ]]
			sample_node$node_height = 0
			
			sample_node$AddChildNode({
				node = Node$new(cn2[1])
				node$node_height = 0
				node
			})
			if(length(cn2) == 2) {
				sample_node$AddSiblingNode({
					node = Node$new(cn2[2])
					node$node_height = 0
					node
				})
			} else {

				foo = sample_node$AddChildNode({
					node = Node$new(paste0(cn2[1], "_node"))
					node$node_height = 0
					node
				})
				for(k in seq_along(cn2)[-1]) {
					nm = cn2[k]
					if(k < length(cn2) - 1) {
						foo$AddChildNode({
							node = Node$new(nm)
							node$node_height = 0
							node
						})
						foo = foo$AddChildNode({
							node = Node$new(paste0(nm, "_node"))
							node$node_height = 0
							node
						})
						
					} else if(k == length(cn2) - 1) {
						foo$AddChildNode({
							node = Node$new(cn2[k])
							node$node_height = 0
							node
						})$AddSiblingNode({
							node = Node$new(cn2[k+1])
							node$node_height = 0
							node
						})
					}
				}
			}

		}
	}
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height")
	od = structure(seq_along(cn), names = cn)
	order.dendrogram(dend) = od[labels(dend)]

	od = structure(1:nleaves(dend), names = labels(dend))
	dend_env = new.env()
	dend_env$dend = dend

	update_midpoint = function(index = NULL) {
		if(is.null(index)) {
			if(is.leaf(dend_env$dend)) {
				pos = od[attr(dend_env$dend, "label")]
				midpoint = 0
			} else {
				x = NULL
				for(i in seq_len(length(dend_env$dend))) {
					if(is.null(attr(dend_env$dend[[i]], "x"))) {
						update_midpoint(i)
					}
					x[i] = attr(dend_env$dend[[i]], "x")
				}
				pos = (max(x) + min(x))/2
				midpoint = (max(x) - min(x))/2
			}
		} else {
			if(is.leaf(dend_env$dend[[index]])) {
				pos = od[attr(dend_env$dend[[index]], "label")]
				midpoint = 0
			} else {
				x = NULL
				for(i in seq_len(length(dend_env$dend[[index]]))) {
					if(is.null(attr(dend_env$dend[[c(index, i)]], "x"))) {
						update_midpoint(c(index, i))
					}
					x[i] = attr(dend_env$dend[[c(index, i)]], "x")
				}
				pos = (max(x) + min(x))/2
				midpoint = (max(x) - min(x))/2
			}
		}
		if(is.null(index)) {
			attr(dend_env$dend, "x") = pos
		} else {
			attr(dend_env$dend[[index]], "x") = pos
			attr(dend_env$dend[[index]], "midpoint") = midpoint
		}
	}
	update_midpoint()

	dend_env$dend
}

most_tight_subgroup = function(mat, subgroup) {
	
	mean_dist = tapply(seq_len(ncol(mat)), subgroup, function(ind) {
		n = length(ind)
		if(n == 1) {
			return(Inf)
		}
		sum(dist(t(mat[, ind, drop = FALSE]))^2)/(n*(n-1)/2)
	})
	as.numeric(names(mean_dist[which.min(mean_dist)]))
}


# mat = matrix(rnorm(100), 10)
# mat2 = cbind(matrix(rnorm(50, mean = -1), nrow = 10), matrix(rnorm(50, mean = 1), nrow = 10))
# subset = sample(10, 10)
# mean_dist_decrease(mat, subset[1:5], subset[6:10])
mean_dist_decrease = function(mat, subset1, subset2) {
	n = ncol(mat)
	global_mean_dist = sum(dist(t(mat))^2) / (n*(n-1)/2)

	n1 = length(subset1)
	n2 = length(subset2)
	mean_dist = (sum(dist(t(mat[, subset1, drop = FALSE]))^2) +  sum(dist(t(mat[, subset2, drop = FALSE]))^2)) / ( n1*(n1-1)/2 + n2*(n2-1)/2)

	(global_mean_dist - mean_dist)/global_mean_dist
}

# == title
# Get class from the HierarchicalPartition object
#
# == param
# -object a `HierarchicalPartition-class` object
# -depth minimal depth of the hierarchy
#
# == return
# A vector of predicted classes.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_class",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL) {

	subgroup = object@subgroup
	if(!is.null(depth)) {
		l = nchar(subgroup) > depth
		subgroup[l] = substr(subgroup[l], 1, depth)
	}
	subgroup
})

# == title
# Print the HierarchicalPartition object
#
# == param
# -object a `HierarchicalPartition-class` object
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "show",
	signature = "HierarchicalPartition",
	definition = function(object) {

	qqcat("A 'HierarchicalPartition' object with '@{object@list[[1]]@top_method}:@{object@list[[1]]@partition_method}' method.\n")
	cat("\n")

	hierarchy = object@hierarchy
	nodes = hierarchy[, 2]
	nc = nchar(nodes)
	names(nc) = nodes
	n = length(nc)

	parent = structure(hierarchy[, 1], names = hierarchy[, 2])

	lines = character(n)
	for(i in seq_len(n)) {
		lines[i] = paste0(strrep("    ", nc[i] - 2), ifelse(grepl("0$", nodes[i]), "`-", "|-") ,"- ", nodes[i], qq(", @{length(object@list[[nodes[i]]]@column_index)} cols"))
		p = nodes[i]
		while(p != "0") {
			p = parent[p]
			if(!grepl("0$", p)) {
				substr(lines[i], (nc[p] - 2)*4+1, (nc[p] - 2)*4+1) = "|"
			}
		}
	}
	# substr(lines[1], 1, 1) = "+"
	lines = c(qq("0, @{length(object@list[['0']]@column_index)} cols"), lines)
	cat(lines, sep = "\n")
	qqcat("\n")
	qqcat("Following methods can be applied to this 'HierarchicalPartition' object:\n")
	txt = showMethods(classes = "HierarchicalPartition", where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(fname)
})

# == title
# Get signatures rows
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth minimal depth of the hierarchy
# 
# == details
# The function called `get_signatures,ConsensusPartition-method` to find signatures at
# each level of the partition hierarchy.
#
# == value
# A list of signature names (row names of the original data matrix)
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_signatures",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL, force_specific = FALSE, ...) {

	all_subgroup_labels = unique(object@subgroup)
	if(!is.null(depth)) {
		all_subgroup_labels = all_subgroup_labels[nchar(all_subgroup_labels) <= depth]
	}

	sig_lt = list()
	for(nm in all_subgroup_labels) {
		nc = nchar(nm)
		which_group = as.numeric(substr(nm, nc, nc))
		parent = substr(nm, 1, nc - 1)
		obj = object@list[[parent]]
		best_k = get_best_k(obj)
		qqcat("get signatures at node @{parent} for @{nm} with @{best_k} subgroups.\n")
		sig = get_signatures(obj, k = best_k, verbose = FALSE, plot = FALSE, ...)$group
		if(which_group != 0) {
			sig_lt[[nm]] = names(sig[sig == which_group])
		} else if(nm == strrep("0", nchar(parent) + 1)) {
			l = grepl(qq("^@{parent}\\d$"), all_subgroup_labels)
			same_level = all_subgroup_labels[l]
			other_group = setdiff(seq_len(best_k), substr(same_level, nc, nc))
			qqcat("  use subgroup @{paste(other_group, collapse = ', ')}\n")
			sig_lt[[nm]] = names(sig[sig %in% other_group])
		} else {
			l = grepl(qq("^@{parent}[1-9]$"), all_subgroup_labels)
			if(any(l)) {
				px = as.numeric(substr(all_subgroup_labels[l], nc, nc))
				sig_lt[[nm]] = names(sig[! sig %in% px])
			}
		}
	}
	all_sig = unique(unlist(sig_lt))

	data = object@list[[1]]@.env$data
	data = t(scale(t(data)))
	m = data[all_sig, , drop = FALSE]
	
	m_sig = matrix(0, nrow = length(all_sig), ncol = length(sig_lt))
	rownames(m_sig) = all_sig
	colnames(m_sig) = names(sig_lt)
	for(i in seq_along(sig_lt)) {
		m_sig[sig_lt[[i]], i] = 1
	}
	col2 = lapply(names(sig_lt), function(x) {
		foo = c("1" = object@subgroup_col[x], "0" = "white")
		names(foo) = c("1", "0")
		foo
	})
	names(col2) = colnames(m_sig)

	if(force_specific) {
		l = rowSums(m_sig) == 1
		m = m[l, ]
		m_sig = m_sig[l, ]
	}

	if(nrow(m) > 5000) {
		qqcat("randomly sample 5000 rows from @{nrow(m)} total rows.\n")
		ind = sample(nrow(m), 5000)
	} else {
		ind = seq_len(nrow(m))
	}
	subgroup = get_class(object, depth = depth)
	col_order = column_order_by_group(subgroup, m)
	mat_range = quantile(abs(m[ind, ]), 0.95)
	col_fun = colorRamp2(c(-mat_range, 0, mat_range), c("green", "white", "red"))

	ht_list = Heatmap(m[ind, ], top_annotation = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = object@subgroup_col), show_annotation_name = TRUE),
		name = "scaled_expr", show_row_names = FALSE, cluster_columns = FALSE, column_order = col_order, col = col_fun,
		show_row_dend = FALSE)
	
	ht_list = ht_list + rowAnnotation(df = m_sig[ind, , drop = FALSE], col = col2, width = unit(0.5*length(sig_lt), "cm"), show_legend = FALSE)
	
	if(force_specific) {
		group = character(nrow(m_sig))
		for(i in seq_len(ncol(m_sig))) {
			group[as.logical(m_sig[, i])] = colnames(m_sig)[i]
		}
		draw(ht_list, split = group[ind])
		return(invisible(list(mat = m, group = group)))
	} else {
		draw(ht_list)
		return(invisible(sig_lt))
	}
})

# == title
# Collect classes from hierarchical_partition object
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth minimal depth of the hierarchy
# -anno a data frame with known annotation of samples.
# -anno_col a list of colors for the annotations in ``anno``.
# -... other arguments.
#
# == details
# The function plots the hierarchy of the subgroups.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "collect_classes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL, anno = object@list[[1]]@known_anno, 
	anno_col = if(missing(anno)) object@list[[1]]@known_col else NULL,
	...) {

	data = object@list[[1]]@.env$data
	if(nrow(data) > 5000) {
		data = data[order(rowSds(data), decreasing = TRUE)[1:5000], ]
	}
	cl = get_class(object, depth = depth)

	hierarchy = object@hierarchy
	if(!is.null(depth)) {
		hierarchy = hierarchy[nchar(hierarchy[, 2]) <= depth , , drop = FALSE]
	}
	dend = calc_dend(object, hierarchy)

	ht_list = Heatmap(cl, name = "subgroups", width = unit(1, "cm"), col = object@subgroup_col,
		row_title_rot = 0, cluster_rows = dend, row_dend_width = unit(2, "cm"))
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
	draw(ht_list)
})

# == title
# Get result for a specified level in the partition hierarchy
#
# == param
# -object a `HierarchicalPartition-class` object
# -node node labal, see `hierarchical_partition` for explaination.
#
# == return 
# A `ConsensusPartition-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_single_run",
	signature = "HierarchicalPartition",
	definition = function(object, node = "0") {
	object@list[[node]]
})

# == title
# Test correspondance between predicted and known classes
#
# == param
# -object a `HierarchicalPartition-class` object
# -depth minimal depth of the hierarchy
# -known a vector or a data frame with known factors
#
# == details
# The function test correlation between classes and known annotations at each node in the hierarchy.
#
# == value
# A matrix of p-values.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "test_to_known_factors",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL, known = object@list[[1]]@known_anno) {

	if(!is.null(known)) {
		if(is.atomic(known)) {
			df = data.frame(known)
			colnames(df) = deparse(substitute(known))
			knwon = df
		}
	}

	all_nodes = sort(unique(as.vector(object@hierarchy)))
	all_nodes = all_nodes[all_nodes != "0"]
	if(!is.null(depth)) {
		all_nodes = all_nodes[nchar(all_nodes) <= depth]
	}

	p = NULL
	for(nm in all_nodes) {
		nc = nchar(nm)
		which_group = as.numeric(substr(nm, nc, nc))
		parent = substr(nm, 1, nc - 1)
		obj = object@list[[parent]]
		best_k = get_best_k(obj)
		
		cl = get_class(obj, k = best_k)$class
		cl[cl != which_group] = 0

		column_index = obj@column_index

		p = rbind(p, test_between_factors(factor(cl), known[column_index, , drop = FALSE]))
	}
	rownames(p) = all_nodes

	class = get_class(object)
	m = test_between_factors(class, known, verbose = FALSE)
	p = rbind(p, overall = m)
	return(p)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object a numeric matrix
# -merge whether merge all samples into one single plot or make plots for every hierarchy
# -depth depth of the hierarchy
# -top_n top n genes to use
# -method which method to reduce the dimension of the data. ``mds`` uses `stats::cmdscale`,
#         ``pca`` uses `stats::prcomp` and ``tsne`` uses `Rtsne::Rtsne`.
# -silhouette_cutoff silhouette cutoff
# -tsne_param parameters pass to `Rtsne::Rtsne`
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "dimension_reduction",
	signature = "HierarchicalPartition",
	definition = function(object, merge = FALSE, depth = NULL,
	top_n = NULL, method = c("mds", "pca", "tsne"),
	silhouette_cutoff = 0.5, tsne_param = list()) {

	method = match.arg(method)
	data = data = object@list[[1]]@.env$data

	if(!is.null(top_n)) {
		top_n = min(c(top_n, nrow(data)))
		all_value = object@list[["0"]]@.env$all_value_list[[object@list[["0"]]@top_method]]
		ind = order(all_value)[1:top_n]
		data = data[ind, , drop = FALSE]
	}

	class = get_class(object, depth = depth)
	if(merge) {
		dimension_reduction(data, pch = 16, col = object@subgroup_col[class],
			cex = 1, main = qq("Dimension reduction by @{method}, @{ncol(data)} samples"),
			method = method, tsne_param = tsne_param)
	} else {
		all_nodes = sort(unique(as.vector(object@hierarchy)))
		all_nodes = all_nodes[all_nodes != "0"]
		if(!is.null(depth)) {
			all_nodes = all_nodes[nchar(all_nodes) <= depth]
		}
		all_parents = all_nodes[nchar(all_nodes) < max(nchar(all_nodes))]
		for(p in all_parents) {
			obj = get_single_run(object, node = p)
			dimension_reduction(obj, k = get_best_k(obj), top_n = top_n, method = method,
				silhouette_cutoff = silhouette_cutoff, tsne_param = tsne_param)
			legend("topleft", legend = qq("node @{p}"))
		}
	}
})

# == title
# max depth of the hierarchy
#
# == param
# -object the object
#
setMethod(f = "max_depth",
	signature = "HierarchicalPartition",
	definition = function(object) {

	max(nchar(object@hierarchy[, 2]))
})

# == title
# all nodes in the hierarchy
#
# == param
# -object the object
# -depth depth
#
setMethod(f = "all_nodes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL) {

	all_nodes = unique(as.vector(object@hierarchy))
	if(!is.null(depth)) {
		all_nodes = all_nodes[nchar(all_nodes) <= depth]
	}
	return(all_nodes)
})

# == title
# all leaves in the hierarchy
#
# == param
# -object the object
# -depth depth
#
setMethod(f = "all_leaves",
	signature = "HierarchicalPartition",
	definition = function(object, depth = NULL) {

	hierarchy = unique(object@hierarchy)
	if(!is.null(depth)) {
		hierarchy = hierarchy[nchar(hierarchy[, 2]) <= depth, , drop = FALSE]
	}
	
	tb = table(hierarchy)
	tb = tb[tb <= 1]
	names(tb)	
})
