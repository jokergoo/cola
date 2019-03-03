
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
# -`guess_best_k,HierarchicalPartition-method`: guess the best number of partitions for each node.
# -`get_matrix,HierarchicalPartition-method`: get the original matrix.
# -`get_signatures,HierarchicalPartition-method`: get the signatures for each subgroup.
# -`dimension_reduction,HierarchicalPartition-method`: make dimension reduction plots.
# -`test_to_known_factors,HierarchicalPartition-method`: test correlation between predicted subgrouping and known annotations, if available.
# -`cola_report,HierarchicalPartition-method`: generate a HTML report for the whole analysis.
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
# Hierarchical partition
#
# == param
# -data a numeric matrix where subgroups are found by columns.
# -top_value_method a single top-value method. Available methods are in `all_top_value_methods`.
# -partition_method a single partition method. Available methods are in `all_partition_methods`.
# -concordance_cutoff the cutoff of concordance scores to determine whether to continue looking for subgroups. Currently it is not used.
# -PAC_cutoff the cutoff of PAC scores to determine whether to continue looking for subgroups.
# -min_samples the cutoff of number of samples to determine whether to continue looking for subgroups.
# -max_k maximal number of partitions to try. The function will try ``2:max_k`` partitions. Note this is the number of
#        partitions that will be tried out on each node of the hierarchical partition. Since more subgroups will be found
#        in the whole partition hierarchy, on each node, ``max_k`` should not be set to a large value.
# -... pass to `consensus_partition`
#
# == details
# The function looks for subgroups in a hierarchical way.
#
# There is a special way to encode the node in the hierarchy. The length of the node name
# is the depth of the node in the hierarchy and the substring excluding the last digit is the name
# node of the parent node. E.g. for the node ``0011``, the depth is 4 and the parent node is ``001``.
#
# == return
# A `HierarchicalPartition-class` object. Simply type object in the interactive R session
# to see which functions can be applied on it.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# set.seed(123)
# m = cbind(rbind(matrix(rnorm(20*20, mean = 2, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
#           rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20)),
#           rbind(matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 0, sd = 0.3), nr = 20),
#                 matrix(rnorm(20*20, mean = 1, sd = 0.3), nr = 20))
#          ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
# cola_rh = hierarchical_partition(m, top_n = c(20, 30, 40), PAC_cutoff = 0.3)
# }
# data(cola_rh)
# cola_rh
hierarchical_partition = function(data, top_value_method = "MAD", partition_method = "kmeans",
	concordance_cutoff = 0.9, PAC_cutoff = 0.2, min_samples = 6, max_k = 4, ...) {

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
		PAC_score = get_stats(part, k = best_k)[, "PAC"]
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

	if(!requireNamespace("data.tree")) {
		stop_wrap("You need to install data.tree package.")
	}
	lt = list()
	lt[["0"]] = data.tree::Node$new("all samples")
	cn = colnames(object@list[["0"]]@.env$data)
	max_depth = max(nchar(hierarchy))
	lt[["0"]]$node_height = max_depth - 1
	for(i in seq_len(nrow(hierarchy))) {
		lt[[ hierarchy[i, 2] ]] = lt[[ hierarchy[i, 1] ]]$AddChildNode({
			node = data.tree::Node$new(hierarchy[i, 2])
			node$node_height = max_depth - nchar(hierarchy[i, 2])
			node
		})
		l = hierarchy[, 1] == hierarchy[i, 2]
	}
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height")

	od = structure(1:length(order.dendrogram(dend)), names = labels(dend))
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
	dend = dendextend::`order.dendrogram<-`(dend, value = 1:nobs(dend))
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
	classes = get_classes(object, depth = depth)[, 1]
	if(is.null(names(classes))) names(classes) = seq_along(classes)
	cd_list = lapply(tapply(names(classes), classes, function(x) x), function(x) {
		d = random_dend(length(x))
		d = dendextend::`labels<-`(d, value = x)
		d
	})
	cd_list = cd_list[labels(pd)]


	dend = merge_dendrogram(pd, cd_list)
	dend = adjust_dend_by_x(dend)
	dend = dendextend::`order.dendrogram<-`(dend, value = structure(1:length(classes), names = names(classes))[labels(dend)])
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
# Get class IDs from the HierarchicalPartition object
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
#
# == return
# A data frame of classes IDs. The class IDs are the node IDs where the subgroup sits in the hierarchy.
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

	if(length(depth) > 1) {
		df = do.call(cbind, lapply(depth, function(d) get_classes(object, d)))
		colnames(df) = paste0("depth=", depth)
		return(df)
	}
	subgroup = object@subgroup
	if(!is.null(depth)) {
		l = nchar(subgroup) > depth
		subgroup[l] = substr(subgroup[l], 1, depth)
	}
	data.frame(class = subgroup, stringsAsFactors = FALSE)
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
	
	# ne = nodes[which.max(nc)[1]]
	# qqcat("e.g. a node '@{ne}' which a parent node '@{gsub(\'.$\', \'\', ne)}'\n")

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
# -depth depth of the hierarchy.
# -scale_rows whether apply row scaling when making the heatmap.
# -anno a data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `hierarchical_partition`.
# -anno_col a list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_column_names whether show column names in the heatmap.
# -verbose whether to print messages.
# -plot whether to make the plot.
# -silhouette_cutoff cutoff for silhouette scores. Samples with values 
#        less than it are not used for finding signature rows. For selecting a 
#        proper silhouette cutoff, please refer to https://www.stat.berkeley.edu/~s133/Cluster2a.html#tth_tAb1.
# -... other arguments
# 
# == details
# The function calls `get_signatures,ConsensusPartition-method` to find signatures at
# each node of the partition hierarchy.
#
# == value
# A list of row indices where rows are significantly different between subgroups in at least one node.
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
	silhouette_cutoff = 0.5, 
	...) {

	if(depth <= 1) {
		stop_wrap("depth should be at least larger than 1.")
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

	class = get_classes(object, depth = depth)[, 1]

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

	row_split = NULL
	if(nrow(use_mat1) > 10) {
		wss <- (nrow(use_mat1)-1)*sum(apply(use_mat1,2,var))
		for (i in 2:15) wss[i] <- sum(kmeans(use_mat1, centers = i, iter.max = 50)$withinss)
		row_km = min(elbow_finder(1:15, wss)[1], knee_finder(1:15, wss)[1])
		if(length(unique(class)) == 1) row_km = 1
		if(length(unique(class)) == 2) row_km = min(row_km, 2)
		if(row_km > 1) {
			row_split = kmeans(use_mat1, centers = row_km)$cluster
		}
		if(verbose) qqcat("  - split rows into @{row_km} groups by k-means clustering.\n")
	}
			

	if(verbose) qqcat("* making heatmaps for signatures\n")

	ha1 = HeatmapAnnotation(
		class = class,
		col = list(class = object@subgroup_col))

	dend = cluster_within_group(use_mat1, class)
	
	ht_list = Heatmap(use_mat1, top_annotation = ha1,
		name = heatmap_name, show_row_names = FALSE, 
		show_column_names = show_column_names, col = col_fun,
		use_raster = TRUE, row_split = row_split,
		show_row_dend = FALSE, cluster_columns = dend,
		column_title = qq("@{length(unique(class))} groups, @{length(all_sig)} signatures"),
		bottom_annotation = bottom_anno1,
		row_title = {if(length(unique(row_split)) <= 1) NULL else qq("k-means with @{length(unique(row_split))} groups")})

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
# Collect classes from HierarchicalPartition object
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
# -anno a data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `hierarchical_partition`.
# -anno_col a list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
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
# data(cola_rh)
# collect_classes(cola_rh)
# collect_classes(cola_rh, depth = 2)
setMethod(f = "collect_classes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object), 
	anno = get_anno(object[1]), anno_col = get_anno_col(object[1])) {

	cl = get_classes(object, depth = depth)[, 1]

	dend = calc_dend(object, depth = depth)

	ht_list = Heatmap(cl, name = "Class", col = object@subgroup_col, width = unit(5, "mm"),
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
	draw(ht_list)
})

# == title
# Subset the HierarchicalPartition object
#
# == param
# -x a `HierarchicalPartition-class` object.
# -i index. The value should be numeric or a node ID.
#
# == details
# On each node, there is a `ConsensusPartition-class` object.
#
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
		stop_wrap("length of the index should only be one.")
	}
	if(is.numeric(i)) {
		i = names(x@list)[i]
	}
	if(is.null(x@list[[i]])) {
		stop_wrap("index for the HierarchicalPartition object cannot be found.")
	}
	return(x@list[[i]])
}

# == title
# Subset the HierarchicalPartition object
#
# == param
# -x a `HierarchicalPartition-class` object
# -i index. The value should be numeric or a node ID.
#
# == details
# On each node, there is a `ConsensusPartition-class` object.
#
# Note you cannot get a sub-hierarchy of the partition.
#
# == value
# A `ConsensusPartition-class` object.
#
# == example
# data(cola_rh)
# cola_rh[["01"]]
# cola_rh[[2]]
"[[.HierarchicalPartition" = function(x, i) {
	x[i]
}


# == title
# Test correspondance between predicted classes and known factors
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
# -known a vector or a data frame with known factors. By default it is the annotation table set in `hierarchical_partition`.
# -verbose whether to print messages.
#
# == value
# A data frame with columns:
#
# - number of samples
# - p-values from the tests
# - number of classes
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rh)
# test_to_known_factors(cola_rh, known = 1:60)
setMethod(f = "test_to_known_factors",
	signature = "HierarchicalPartition",
	definition = function(object, known = get_anno(object[1]),
	depth = 2:max_depth(object), verbose = FALSE) {

	if(!is.null(known)) {
		if(is.atomic(known)) {
			df = data.frame(known)
			colnames(df) = deparse(substitute(known))
			known = df
		}
	} else {
		stop_wrap("Known factors should be provided.")
	}

	class = get_classes(object, depth)
	m = test_between_factors(class, known, verbose = verbose)
	colnames(m) = paste0(colnames(m), "(p-value)")
	df = cbind(n = nrow(class), m, n_class = apply(class, 2, function(x) length(unique(x))))
	return(df)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object a `HierarchicalPartition-class` object.
# -depth depth of the hierarchy.
# -top_n top n rows to use. By default it uses all rows in the original matrix.
# -parent_node parent node. If it is set, the function call is identical to ``dimension_reduction(object[parent_node])``
# -method which method to reduce the dimension of the data. ``mds`` uses `stats::cmdscale`,
#         ``pca`` uses `stats::prcomp`.
# -silhouette_cutoff cutoff of silhouette score. Data points with values less
#        than it will be mapped to small points.
# -scale_rows whether perform scaling on matrix rows.
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
# data(cola_rh)
# dimension_reduction(cola_rh)
# dimension_reduction(cola_rh, parent_node = "00")
setMethod(f = "dimension_reduction",
	signature = "HierarchicalPartition",
	definition = function(object,
	depth = max_depth(object), parent_node,
	top_n = NULL, method = c("PCA", "MDS"),
	silhouette_cutoff = 0.5, scale_rows = TRUE) {

	method = match.arg(method)
	data = object@list[[1]]@.env$data

	op = par(c("mar", "xpd"))
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)

	if(missing(parent_node)) {
		if(!is.null(top_n)) {
			top_n = min(c(top_n, nrow(data)))
			all_top_value = object@list[["0"]]@.env$all_top_value_list[[object@list[["0"]]@top_value_method]]
			ind = order(all_top_value)[1:top_n]
			data = data[ind, , drop = FALSE]
		} else {
			top_n = nrow(data)
		}
		class = get_classes(object, depth = depth)[, 1]
		n_class = length(unique(class))
		dimension_reduction(data, pch = 16, col = object@subgroup_col[class],
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{object@list[[1]]@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{ncol(data)} samples at depth @{depth} with @{n_class} classes"),
			method = method, scale_rows = scale_rows)
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
		dimension_reduction(obj, k = guess_best_k(obj), top_n = top_n, method = method,
			silhouette_cutoff = silhouette_cutoff, scale_rows = scale_rows)
		legend(x = par("usr")[2], y = par("usr")[4], legend = qq("node @{parent_node}"))
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
# A vector of node ID.
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
# A vector of node ID.
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
# -object a `HierarchicalPartition-class` object.
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
# -object a `HierarchicalPartition-class` object.
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

# == title
# Guess the best number of partitions
#
# == param
# -object a `HierarchicalPartition-class` object.
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
# data(cola_rh)
# guess_best_k(cola_rh)
setMethod(f = "guess_best_k",
	signature = "HierarchicalPartition",
	definition = function(object) {

	best_k = NULL
	cophcor = NULL
	PAC = NULL
	mean_silhouette = NULL
	concordance = NULL
	for(nm in names(object@list)) {
		obj = object@list[[nm]]
		best_k[nm] = guess_best_k(obj)
		stat = get_stats(obj, k = best_k[nm])
		cophcor[nm] = stat[1, "cophcor"]
		PAC[nm] = stat[1, "PAC"]
		mean_silhouette[nm] = stat[1, "mean_silhouette"]
		concordance[nm] = stat[1, "concordance"]
	}
	tb = data.frame(best_k = best_k,
		cophcor = cophcor,
		PAC = PAC,
		mean_silhouette = mean_silhouette,
		concordance = concordance)

	l = rownames(tb) %in% all_leaves(object)
	tb = cbind(tb, ifelse(l, "leaf", "node"), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""

	return(tb)
})

