
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
# -`get_classes,HierarchicalPartition-method`: get the class IDs of subgroups.
# -`get_signatures,HierarchicalPartition-method`: get the signatures for each subgroup.
# -`test_to_known_factors,HierarchicalPartition-method`: test correlation between predicted subgrouping and known annotations, if available.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
HierarchicalPartition = setClass("HierarchicalPartition",
    slots = list(
        list = "list",
        best_k = "numeric",
        hierarchy = "matrix",
        subgroup = "character",
        subgroup_col = "character",
        calling = "ANY",
        .env = "environment"
    )
)

# == title
# Hierarchical detection of subgroups
#
# == param
# -data a numeric matrix where subgroups are found by samples.
# -top_value_method a single top method. Available methods are in `all_top_value_methods`.
# -partition_method a single partition method. Available ialble methods are in `all_partition_methods`.
# -concordance_cutoff the cutoff of PAC scores to determine whether to continuou looking to subgroups.
# -PAC_cutoff PAC_cutoff
# -min_samples the cutoff of number of samples to determine whether to continuou looking to subgroups.
# -max_k a list number of partitions.
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
hierarchical_partition = function(data, top_value_method = "MAD", partition_method = "kmeans",
	concordance_cutoff = 0.9, PAC_cutoff = 0.15, min_samples = 6, max_k = 4, ...) {

	cl = match.call()
	
	.h_obj = new.env()
	.h_obj$hierarchy = matrix(nrow = 0, ncol = 2)
	.h_obj$list = list()
	.h_obj$subgroup = rep("", ncol(data))

	.hierarchical_partition = function(.env, column_index, concordance_cutoff = 0.9, node_id = '0', 
		parent = NULL, min_samples = 6, max_k = 4, ...) {

		if(node_id != "0") {
			cat("\n")
			.h_obj$hierarchy = rbind(.h_obj$hierarchy, c(parent, node_id))
		}
		qqcat("* submatrix with @{length(column_index)} columns and @{nrow(data)} rows, node_id: @{node_id}.\n")
	   	.env$all_top_value_list = NULL
	   	.env$column_index = column_index
		part = consensus_partition(verbose = FALSE, .env = .env, max_k = max_k, 
			top_value_method = top_value_method, partition_method = partition_method, ...)

	    .h_obj$list[[node_id]] = part
	    .h_obj$subgroup[column_index] = node_id

	    best_k = guess_best_k(part)
	    cl = get_classes(part, k = best_k)
	    .h_obj$best_k[node_id] = best_k

		# oe = try(sig <- get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE))
		# if(inherits(oe, "try-error")) {
		# 	cat("caught an error, stop.\n")
		# 	return(NULL)
		# }
		# if(length(sig$group) < 100) {
		# 	qqcat("Number of signatures are too small (< 100), stop.\n")
		# 	return(NULL)
		# }

	    mat = .env$data[, column_index, drop = FALSE]
	    if(part@scale_rows) {
	    	mat = t(scale(t(mat)))
	    }
	    mat = mat[order(.env$all_top_value_list[[1]], decreasing = TRUE), , drop = FALSE]
	    s = farthest_subgroup(mat, cl$class, part@top_n)
	    set1 = which(cl$class == s)
	    set2 = which(cl$class != s)
	    .env$all_top_value_list = NULL

	    qqcat("* best k = @{best_k}, the farthest subgroup: @{s}\n")
		PAC_score = get_stat(part, k = best_k)[, "PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	qqcat("* PAC score is too big @{sprintf('%.2f', PAC_score)}, stop.\n")
	    	return(NULL)
	    }

	    # dist_decrease = mean_dist_decrease(mat, set1, set2)
	    # if(dist_decrease < reduce) {
	    # 	qqcat("mean distance does not decrease too much @{dist_decrease}, stop.\n")
	    # 	return(NULL)
	    # } else {

	    	if(length(set1) <= min_samples || length(set2) <= min_samples) {
	    		cat("* subgroups have too few columns, stop.\n")
	    		return(NULL)
	    	}

	    	qqcat("* partitioned into two subgroups with @{length(set1)} and @{length(set2)} columns.\n")
	    	# insert the two subgroups into the hierarchy
	    	sub_node_1 = paste0(node_id, s)
	    	sub_node_2 = paste0(node_id, "0")

	    	if(length(set1) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set1], node_id = sub_node_1, parent = node_id,
	    			concordance_cutoff = concordance_cutoff, min_samples = min_samples, max_k = max_k, ...)
	    	}

	    	if(length(set2) > min_samples) {
	    		.hierarchical_partition(.env, column_index = column_index[set2], node_id = sub_node_2, parent = node_id,
	    			concordance_cutoff = concordance_cutoff, min_samples = min_samples, max_k = max_k, ...)
	    	}
	    	
	    # }

	    return(NULL)
	}

	.env = new.env()
	.env$data = data
	.hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), concordance_cutoff = concordance_cutoff, min_samples = min_samples, 
		node_id = "0", parent = NULL, max_k = max_k, ...)

	hp = new("HierarchicalPartition")
	hp@hierarchy = .h_obj$hierarchy
	hp@list = .h_obj$list
	hp@best_k = .h_obj$best_k
	hp@subgroup = .h_obj$subgroup
	names(hp@subgroup) = colnames(data)

	le = unique(as.vector(hp@hierarchy))
	if(length(le) <= 16) {
		hp@subgroup_col = structure(brewer_pal_set2_col[seq_along(le)], names = le)
	} else {
		hp@subgroup_col = structure(rand_color(length(le), luminosity = "bright"), names = le)
	}
	hp@calling = cl
	hp@.env = .h_obj$list[[1]]@.env

	return(hp)
}

subgroup_dend = function(object, hierarchy = object@hierarchy) {

	lt = list()
	lt[["0"]] = Node$new("all samples")
	cn = colnames(object@list[["0"]]@.env$data)
	max_depth = max(nchar(hierarchy))
	lt[["0"]]$node_height = max_depth + 1
	for(i in seq_len(nrow(hierarchy))) {
		lt[[ hierarchy[i, 2] ]] = lt[[ hierarchy[i, 1] ]]$AddChildNode({
			node = Node$new(hierarchy[i, 2])
			node$node_height = max_depth - nchar(hierarchy[i, 2])
			node
		})
		l = hierarchy[, 1] == hierarchy[i, 2]
	}
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height")

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

	# reorder(dend_env$dend, wts = order(labels(dend_env$dend)))
	dend = dend_env$dend
	dendrapply(dend, function(d) {
		if(is.leaf(d)) {
			attr(d, "height") = 0
		}
		d
	})
}

get_hierarchy = function(object, depth = max_depth(object)) {

	hierarchy = object@hierarchy
	if(!is.null(depth)) {
		hierarchy = hierarchy[nchar(hierarchy[, 2]) <= depth , , drop = FALSE]
	}

	dend = subgroup_dend(object, hierarchy)
	order.dendrogram(dend) = 1:nobs(dend)
	dend
}

random_dend = function(n) {
	x = rnorm(n)
	dend = as.dendrogram(hclust(dist(1:n)))
	# set height to zero

	dendrapply(dend, function(x) {attr(x, "height") = 0; x})
}


calc_dend = function(object, depth = max_depth(object)) {

	pd = get_hierarchy(object, depth = depth)
	classes = get_classes(object, depth = depth)
	tb = table(classes)
	cd_list = lapply(tapply(names(classes), classes, function(x) x), function(x) {
		d = random_dend(length(x))
		labels(d) = x
		d
	})
	cd_list = cd_list[labels(pd)]


	dend = merge_dendrogram(pd, cd_list)
	dend = adjust_dend_by_x(dend)
	order.dendrogram(dend) = structure(1:length(classes), names = names(classes))[labels(dend)]
	dend
}

tightest_subgroup = function(mat, subgroup, top_n) {
	x = NULL
	for(n in top_n) {
		m2 = mat[1:n, , drop = FALSE]
		mean_dist = tapply(seq_len(ncol(m2)), subgroup, function(ind) {
			n = length(ind)
			if(n == 1) {
				return(Inf)
			}
			sum(dist(t(m2[, ind, drop = FALSE]))^2)/(n*(n-1)/2)
		})
		x = c(x, as.numeric(names(mean_dist[which.min(mean_dist)])))
	}
	tb = table(x)
	as.numeric(names(which.max(tb)))
}


farthest_subgroup = function(mat, subgroup, top_n) {
	x = NULL
	for(n in top_n) {
		m2 = mat[1:n, , drop = FALSE]
		group_mean = do.call("cbind", tapply(seq_len(ncol(m2)), subgroup, function(ind) {
			rowMeans(m2[, ind, drop = FALSE])
		}))
		d = as.matrix(dist(t(group_mean)))
		x = c(x, as.numeric(names(which.max(rowSums(d)))))
	}
	tb = table(x)
	as.numeric(names(which.max(tb)))
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
# == example
# data(cola_rh)
# get_classes(cola_rh)
# get_classes(cola_rh, depth = 2)
setMethod(f = "get_classes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object)) {

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
# == example
# data(cola_rh)
# cola_rh
setMethod(f = "show",
	signature = "HierarchicalPartition",
	definition = function(object) {

	qqcat("A 'HierarchicalPartition' object with '@{object@list[[1]]@top_value_method}:@{object@list[[1]]@partition_method}' method.\n")
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{ncol(object@.env$data)} columns.\n")
	qqcat("  Performed in total @{object@list[[1]]@n_partition*length(object@list)} partitions.\n")
	qqcat("  There are @{length(all_leaves(object))} groups.\n")
	cat("\n")

	cat("Hierarchy of the partition:\n")
	hierarchy = object@hierarchy
	nodes = hierarchy[, 2]
	nc = nchar(nodes)
	names(nc) = nodes
	n = length(nc)

	parent = structure(hierarchy[, 1], names = hierarchy[, 2])

	lines = character(n)
	for(i in seq_len(n)) {
		lines[i] = paste0("  ", strrep("    ", nc[i] - 2), ifelse(grepl("0$", nodes[i]), "`-", "|-") ,"- ", nodes[i], qq(", @{length(object@list[[nodes[i]]]@column_index)} cols"))
		p = nodes[i]
		while(p != "0") {
			p = parent[p]
			if(!grepl("0$", p)) {
				substr(lines[i], (nc[p] - 2)*4+3, (nc[p] - 2)*4+3) = "|"
			}
		}
	}
	# substr(lines[1], 1, 1) = "+"
	lines = c(qq("  0, @{length(object@list[['0']]@column_index)} cols"), lines)
	cat(lines, sep = "\n")
	qqcat("\n")
	qqcat("Following methods can be applied to this 'HierarchicalPartition' object:\n")
	txt = showMethods(classes = "HierarchicalPartition", where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(fname)
	cat("\n")
	cat("You can get result for a single node by e.g. object[\"01\"]\n")
})

# == title
# Get signatures rows
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth minimal depth of the hierarchy
# -scale_rows whether scale rows
# -anno annotation
# -anno_col annotation color
# -show_column_names whether show column names
# -verbose whether print messages
# -plot whether make the heatmap
# -silhouette_cutoff silhouette_cutoff
# -... other arguments
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
# == example
# \dontrun{
# data(cola_rh)
# get_signatures(cola_rh)
# }
setMethod(f = "get_signatures",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object), 
	scale_rows = object[1]@scale_rows, 
	anno = get_anno(object), 
	anno_col = get_anno_col(object),
	show_column_names = FALSE, 
	verbose = TRUE, plot = TRUE,
	silhouette_cutoff = -Inf, 
	...) {

	if(depth <= 1) {
		stop("depth should be at least larger than 1.")
	}

	alf = all_leaves(object, depth = depth)
	ap = setdiff(all_nodes(object, depth = depth), all_leaves(object, depth = depth))

	sig_lt = list()
	for(p in ap) {
		best_k = guess_best_k(object[[p]])
		if(verbose) qqcat("* get signatures at node @{p} with @{best_k} subgroups.\n")
		sig_tb = get_signatures(object[[p]], k = best_k, verbose = FALSE, plot = FALSE, silhouette_cutoff = silhouette_cutoff, ...)
		
		sig_lt[[p]] = sig_tb
		if(verbose) qqcat("  - find @{nrow(sig_tb)} signatures at node @{p}\n")
	}

	all_index = sort(unique(unlist(lapply(sig_lt, function(x) x[, 1]))))

	sig = data.frame(which_row = all_index)
	
	if(!plot) {
		return(invisible(sig))
	}

	all_sig = all_index
	if(verbose) qqcat("* in total @{length(all_sig)} signatures in all classes found\n")

	m = object@.env$data[all_sig, , drop = FALSE]

	class = get_classes(object, depth = depth)

	if(nrow(m) > 2000) {
		if(verbose) qqcat("* randomly sample @{2000} rows from @{nrow(m)} total rows.\n")
		ind = sample(nrow(m), 2000)
	} else {
		ind = seq_len(nrow(m))
	}
	mat1 = m[ind, , drop = FALSE]

	base_mean = rowMeans(mat1)
	if(nrow(mat1) == 1) {
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

	if(verbose) qqcat("* making heatmaps for signatures\n")

	ha1 = HeatmapAnnotation(
		class = class,
		col = list(class = object@subgroup_col))

	ht_list = Heatmap(use_mat1, top_annotation = ha1,
		name = heatmap_name, show_row_names = FALSE, 
		show_column_names = show_column_names,
		column_split = class, col = col_fun,
		use_raster = TRUE,
		show_row_dend = FALSE, show_column_dend = FALSE,
		column_title = qq("@{length(unique(class))} groups, @{length(all_sig)} signatures"),
		bottom_annotation = bottom_anno1)

	all_value_positive = !any(m < 0)
 	if(scale_rows && all_value_positive) {
		ht_list = ht_list + Heatmap(base_mean, show_row_names = FALSE, name = "base_mean", width = unit(5, "mm")) +
			Heatmap(rel_diff, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
				show_row_names = FALSE, name = "rel_diff", width = unit(5, "mm"))
	}

	draw(ht_list)
	return(invisible(sig))
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
# == example
# data(cola_rh)
# collect_classes(cola_rh)
# collect_classes(cola_rh, depth = 2)
setMethod(f = "collect_classes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object), 
	anno = get_anno(object[1]), anno_col = get_anno_col(object[1]), ...) {

	data = object@list[[1]]@.env$data
	if(nrow(data) > 5000) {
		data = data[order(rowSds(data), decreasing = TRUE)[1:5000], ]
	}
	cl = get_classes(object, depth = depth)

	dend = calc_dend(object, depth = depth)

	ht_list = Heatmap(cl, name = "Groups", width = ht_opt$anno_simple_size, col = object@subgroup_col,
		row_title_rot = 0, cluster_rows = dend, row_dend_width = unit(2, "cm"),
		show_row_names = FALSE)
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
# Subset the HierarchicalPartition object
#
# == param
# -x a `HierarchicalPartition-class` object
# -i index
#
# == details
# Note you cannot get a sub-hierarchy of the partition.
#
# == value
# A `ConsensusPartition-class` object.
#
# == example
# data(cola_rh)
# cola_rh["01"]
# cola_rh["2"]
"[.HierarchicalPartition" = function(x, i) {
	if(length(i) > 1) {
		stop("length of the index should only be one.")
	}
	if(is.numeric(i)) {
		i = names(x@list)[i]
	}
	x@list[[i]]
}

# == title
# Subset the HierarchicalPartition object
#
# == param
# -x a `HierarchicalPartition-class` object
# -i index
#
# == details
# Note you cannot get a sub-hierarchy of the partition.
#
# == value
# A `ConsensusPartition-class` object.
#
# == example
# data(cola_rh)
# cola_rh[["01"]]
# cola_rh[["2"]]
"[[.HierarchicalPartition" = function(x, i) {
	x[i]
}


# == title
# Test correspondance between predicted and known classes
#
# == param
# -object a `HierarchicalPartition-class` object.
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
	definition = function(object, depth = max_depth(object), 
	known = get_anno(object[1])) {

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
		best_k = guess_best_k(obj)
		
		cl = get_classes(obj, k = best_k)$class
		cl[cl != which_group] = 0

		column_index = obj@column_index

		p = rbind(p, test_between_factors(factor(cl), known[column_index, , drop = FALSE]))
	}
	rownames(p) = all_nodes

	class = get_classes(object)
	m = test_between_factors(class, known, verbose = FALSE)
	p = rbind(p, overall = m)
	return(p)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
# -top_n top n genes to use.
# -parent_node parent node. If it is set, the function calls is identical to ``dimension_reduction(object[parent_node])``
# -method which method to reduce the dimension of the data. ``mds`` uses `stats::cmdscale`,
#         ``pca`` uses `stats::prcomp`.
# -silhouette_cutoff silhouette cutoff
#
# == details
# The classes is extract at depth ``depth``.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# dimension_reduction(cola_rh)
# dimension_reduction(cola_rh, parent_node = "00")
setMethod(f = "dimension_reduction",
	signature = "HierarchicalPartition",
	definition = function(object,
	depth = max_depth(object), parent_node,
	top_n = NULL, method = c("PCA", "MDS"),
	silhouette_cutoff = 0.5) {

	method = match.arg(method)
	data = object@list[[1]]@.env$data

	op = par(c("mar", "xpd"))
	par(mar = c(4.1, 4.1, 4.1, 5.1), xpd = NA)

	if(missing(parent_node)) {
		if(!is.null(top_n)) {
			top_n = min(c(top_n, nrow(data)))
			all_top_value = object@list[["0"]]@.env$all_top_value_list[[object@list[["0"]]@top_value_method]]
			ind = order(all_top_value)[1:top_n]
			data = data[ind, , drop = FALSE]
		}
		class = get_classes(object, depth = depth)
		dimension_reduction(data, pch = 16, col = object@subgroup_col[class],
			cex = 1, main = qq("Dimension reduction by @{method}, @{ncol(data)} samples"),
			method = method)
		class_level = sort(unique(class))
		legend(x = par("usr")[2], y = par("usr")[4], legend = paste0("group", class_level), pch = 16,
				col = object@subgroup_col[class_level], adj = c(0, 1))
	} else {
		if(!parent_node %in% setdiff(all_nodes(object), all_leaves(object))) {
			stop(qq("@{parent_node} has no children nodes.\n"))
		}
		obj = object[parent_node]
		dimension_reduction(obj, k = guess_best_k(obj), top_n = top_n, method = method,
			silhouette_cutoff = silhouette_cutoff)
		legend("topleft", legend = qq("node @{parent_node}"))
	}

	par(op)
})

# == title
# Max depth of the hierarchy
#
# == param
# -object a `HierarchicalPartition-class` object.
#
# == value
# A numeric value.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# max_depth(cola_rh)
setMethod(f = "max_depth",
	signature = "HierarchicalPartition",
	definition = function(object) {

	max(nchar(object@hierarchy[, 2]))
})

# == title
# All nodes in the hierarchy
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth in the hierarchy.
#
# == value
# A vector of node ids.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# all_nodes(cola_rh)
setMethod(f = "all_nodes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object)) {

	all_nodes = unique(as.vector(object@hierarchy))
	if(!is.null(depth)) {
		all_nodes = all_nodes[nchar(all_nodes) <= depth]
	}
	return(all_nodes)
})

# == title
# All leaves in the hierarchy
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth in the hierarchy.
#
# == value
# A vector of node ids.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# all_leaves(cola_rh)
setMethod(f = "all_leaves",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object)) {

	hierarchy = unique(object@hierarchy)
	if(!is.null(depth)) {
		hierarchy = hierarchy[nchar(hierarchy[, 2]) <= depth, , drop = FALSE]
	}
	
	tb = table(hierarchy)
	tb = tb[tb <= 1]
	names(tb)	
})

get_children = function(object, node = "0") {
	hierarchy = unique(object@hierarchy)
	hierarchy[hierarchy[, 1] == node, 2]
}

# == title
# Get the original matrix
#
# == param
# -object a `HierarchicalPartition-class` object.
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_matrix",
	signature = "HierarchicalPartition",
	definition = function(object) {
	object@.env$data
})


# == title
# Get annotations
#
# == param
# -object a `HierarchicalPartition-class` object
#
# == value
# A data frame if ``anno`` was specified in `hierarchical_partition`, or ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno",
	signature = "HierarchicalPartition",
	definition = function(object) {
	get_anno(object[1])
})



# == title
# Get annotation colors
#
# == param
# -object a `HierarchicalPartition-class` object
#
# == value
# A list of color vectors or ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno_col",
	signature = "HierarchicalPartition",
	definition = function(object) {
	get_anno_col(object[1])
})
