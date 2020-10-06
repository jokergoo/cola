
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
# -`suggest_best_k,HierarchicalPartition-method`: guess the best number of partitions for each node.
# -`get_matrix,HierarchicalPartition-method`: get the original matrix.
# -`get_signatures,HierarchicalPartition-method`: get the signatures for each subgroup.
# -`compare_signatures,HierarchicalPartition-method`: compare signatures from different nodes.
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
        n_columns = "numeric",
        n_signatures = "numeric",
        call = "ANY",
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
# -PAC_cutoff the cutoff of PAC scores to determine whether to continue looking for subgroups.
# -min_samples the cutoff of number of samples to determine whether to continue looking for subgroups.
# -subset Number of columns to randomly sample.
# -min_n_signatures Minimal number of signatures under the best classification.
# -min_p_signatures Minimal proportion of signatures under the best classification. If the corresponding values
#     are smaller than both ``min_n_signatures`` and ``min_p_signatures``, the hierarchical partitioning stops on that node.
# -max_k maximal number of partitions to try. The function will try ``2:max_k`` partitions. Note this is the number of
#        partitions that will be tried out on each node of the hierarchical partition. Since more subgroups will be found
#        in the whole partition hierarchy, on each node, ``max_k`` should not be set to a large value.
# -verbose whether print message.
# -mc.cores multiple cores to use. 
# -help Whether to show the help message.
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
# == examples
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
# rh = hierarchical_partition(m, top_value_method = "SD", partition_method = "kmeans")
# }
hierarchical_partition = function(data, 
	top_value_method = "ATC", 
	partition_method = "skmeans",
	combination_method =  expand.grid(top_value_method, partition_method),
	PAC_cutoff = 0.2, min_samples = 6, subset = Inf,
	min_n_signatures = round(nrow(data)*min_p_signatures), 
	min_p_signatures = 0.01,
	max_k = 4, verbose = TRUE, mc.cores = 1, help = TRUE, ...) {

	if(help) {
		message_wrap("We suggest to try both 'ATC/skmeans' and 'SD/kmeans' for 'top_value_method' and 'partition_method' parameters. These two combinations of methods are correlation-based and Euclidean distance-based respectively and they generate different results that are all worth to look at. Set the argument 'help = FALSE' to turn off this message.")
	}

	cl = match.call()

	if(min_samples < 3) {
		if(verbose) qqcat("! 'min_samples' was reset to 3.\n")
		min_samples = 3
	}

	check_pkg("Polychrome", bioc = FALSE)

	if(verbose) {
		qqcat("* on a @{nrow(data)}x@{ncol(data)} matrix.\n")
		qqcat("* hierarchical partition by @{top_value_method}:@{partition_method}.\n")
	}

	if(is.data.frame(combination_method)) combination_method = as.matrix(combination_method)
	# convert combination_method to a list
	if(is.matrix(combination_method)) {
		combination_method = lapply(seq_len(nrow(combination_method)), function(i) combination_method[i, ])
	} else if(is.atomic(combination_method)) {
		if(!all(grepl(":+", combination_method))) {
			stop_wrap("If `combination_method` is specified as a character vector, the individual values in it should be in the form: 'top_value_method:partition_method', e.g. 'SD:kmeans'.")
		}
		combination_method = strsplit(combination_method, ":+")
	} else {
		stop_wrap("Wrong format of `combination_method`.")
	}
	
	.hierarchical_partition = function(.env, column_index, node_id = '0', 
		min_samples = 6, max_k = 4, verbose = TRUE, mc.cores = 1, ...) {

		prefix = ""
		if(node_id != "0") {
			prefix = paste(rep("  ", nchar(node_id) - 1), collapse = "")
		}

		if(verbose) qqcat("@{prefix}=========================================================\n")
		if(verbose) qqcat("@{prefix}* submatrix with @{length(column_index)} columns, node_id: @{node_id}.\n")

		if(length(column_index) < 2*min_samples) {
			if(verbose) qqcat("@{prefix}* number of samples is not enough to perform partitioning.\n")
			part = "Not enough samples"
			# we need the following two values for other functions
			attr(part, "node_id") = node_id
			attr(part, "column_index") = column_index
			return(list(obj = part))
		}

		## all_top_value_list is only used in run_all_consensus_partition_methods(), we remove it here
	   	.env$all_top_value_list = NULL
	   	.env$column_index = column_index  #note .env$column_index is only for passing to `consensus_partition()` function
		
		if(node_id != "0") {
			row_sd = rowSds(data[, column_index, drop = FALSE])
			qa = quantile(unique(row_sd[row_sd > 1e-10]), 0.05, na.rm = TRUE)
			l = row_sd > qa
			.env$row_index = which(l)
			if(verbose) qqcat("@{prefix}* @{sum(!l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")
		}
		if(mc.cores > 1 && verbose) {
			qqcat("@{prefix}* running consensus partitioning with @{mc.cores} cores.\n")
		}

		part_list = list(length(combination_method))
		if(length(column_index) <= subset) {
			for(i in seq_along(combination_method)) {
				if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
				part_list[[i]] = consensus_partition(verbose = TRUE, .env = .env, max_k = max_k, prefix = prefix,
					top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					mc.cores = mc.cores, ...)
			}
		} else {
			for(i in seq_along(combination_method)) {
				if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
				part_list[[i]] = consensus_partition_by_down_sampling(subset = subset, verbose = TRUE, .env = .env, max_k = max_k, prefix = prefix,
					top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					mc.cores = mc.cores, ...)
			}
		}
		if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")

		ind = which.max(sapply(part_list, function(part) {
			stat_df = get_stats(part)
			best_k = suggest_best_k(part)
			if(is.na(best_k)) {
				0
			} else {
				stat_df[stat_df[, "k"] == best_k, "1-PAC"]
			}
		}))
		part = part_list[[ind]]

		if(verbose) qqcat("@{prefix}* select @{part@top_value_method}:@{part@partition_method} with the highest 1-PAC among methods.\n")

		attr(part, "node_id") = node_id

		lt = list(obj = part)

	    best_k = suggest_best_k(part)
	    if(is.na(best_k)) {
	    	attr(lt$obj, "stop_reason") = "Rand indices for all k were too high."
	    	if(verbose) qqcat("@{prefix}* Rand index is too high, no meaningful subgroups, stop.\n")
	    	return(lt)
	    }
	    cl = get_classes(part, k = best_k)

	    # check PAC score
	    PAC_score = 1 - get_stats(part, k = best_k)[, "1-PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	if(verbose) qqcat("@{prefix}* PAC score is too big (@{sprintf('%.2f', PAC_score)}), stop.\n")
	    	attr(lt$obj, "stop_reason") = "PAC score was too big."
	    	return(lt)
	    }

	    .env$all_top_value_list = NULL

	    ## test the number of samples
	    sample_too_small = FALSE
	    if(best_k == 2) {
	    	kg1 = 1
	    	kg2 = 2
	    	set1 = which(cl$class == kg1)
	    	set2 = which(cl$class == kg2)

	    	if(length(set1) < min_samples || length(set2) < min_samples) {
	    		sample_too_small = TRUE
	    	}
	    } else {
	    	tb = table(cl$class)

	    	# method1: take the largest group as group 1
	    	# kg1 = as.numeric(names(tb[which.max(tb)[1]]))
	    	# kg2 = sort(setdiff(cl$class, kg1))
	    	# set1 = which(cl$class %in% kg1)
	    	# set2 = which(cl$class %in% kg2)

	    	# method2: take the group with minimal within-group mean distance as group 1
	    	#   the class labels were already sorted by that
	    	sample_too_small = TRUE
	    	for(ik in sort(as.numeric(names(tb)))) {
	    		kg1 = ik
	    		kg2 = sort(setdiff(cl$class, kg1))

	    		set1 = which(cl$class %in% kg1)
	    		set2 = which(cl$class %in% kg2)

	    		if(length(set1) >= min_samples && length(set2) >= min_samples) {
		    		sample_too_small = FALSE
		    		break
		    	}
	    	}
	    }

    	if(sample_too_small) {
    		if(verbose) qqcat("@{prefix}* some of the subgroups have too few columns (< @{min_samples}), won't split.\n")
    		attr(lt$obj, "stop_reason") = "Subgroup had too few columns."
    		return(lt)
    	}

    	# check the numbers of signatures
    	if(verbose) qqcat("@{prefix}* checking number of signatures in the best classification.\n")
    	sig_df = get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE, simplify = TRUE)
    	n_sig = nrow(sig_df)
    	p_sig = n_sig/nrow(part)
    	if(n_sig <= min_n_signatures && p_sig <= min_p_signatures) {
    		if(verbose) qqcat("@{prefix}* too few signatures (n = @{n_sig}, p = @{sprintf('%.3f', p_sig)}) under the best classification, stop.\n")
    		attr(lt$obj, "stop_reason") = "There were too few signatures."
    		return(lt)
    	}

    	if(verbose) qqcat("@{prefix}* best k = @{best_k}, split into two groups with class IDs (@{paste(kg1, collapse=',')}) and (@{paste(kg2, collapse=',')})\n")
		if(verbose) qqcat("@{prefix}* partition into two subgroups with @{length(set1)} and @{length(set2)} columns.\n")
    	# insert the two subgroups into the hierarchy
    	sub_node_1 = paste0(node_id, 1)
    	sub_node_2 = paste0(node_id, 0)

    	lt2 = lapply(1:2, function(ind) {
	    	if(ind == 1) {
	    		return(.hierarchical_partition(.env, column_index = column_index[set1], node_id = sub_node_1,
	    			min_samples = min_samples, max_k = min(max_k, length(set1)-1), mc.cores = mc.cores, verbose = verbose, ...))
	    	}

	    	if(ind == 2) {
	    		return(.hierarchical_partition(.env, column_index = column_index[set2], node_id = sub_node_2,
	    			min_samples = min_samples, max_k = min(max_k, length(set2)-1), mc.cores = mc.cores, verbose = verbose, ...))
	    	}

	    	return(NULL)
	    })

	    if(!is.null(lt2[[1]])) lt$child1 = lt2[[1]]
	    if(!is.null(lt2[[2]])) lt$child2 = lt2[[2]]

	    return(lt)
	}

	.env = new.env(parent = emptyenv())
	.env$data = data
	lt = .hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), min_samples = min_samples, 
		node_id = "0", max_k = min(max_k, ncol(data)-1), verbose = verbose, mc.cores = mc.cores, ...)

	# reformat lt
	.e = new.env(parent = emptyenv())
	.e$hierarchy = matrix(nrow = 0, ncol = 2)
	.e$lt = list()
	reformat_lt = function(lt, .e) {
		nm = names(lt)
		parent_id = attr(lt$obj, "node_id")
		.e$lt[[parent_id]] = lt$obj
		if("child1" %in% nm) {
			child_id = attr(lt$child1$obj, "node_id")
			.e$hierarchy = rbind(.e$hierarchy, c(parent_id, child_id))
			reformat_lt(lt$child1, .e)
		}
		if("child2" %in% nm) {
			child_id = attr(lt$child2$obj, "node_id")
			.e$hierarchy = rbind(.e$hierarchy, c(parent_id, child_id))
			reformat_lt(lt$child2, .e)
		}
	}
	reformat_lt(lt, .e)

	hp = new("HierarchicalPartition")
	hp@hierarchy = .e$hierarchy
	hp@list = .e$lt
	hp@best_k = sapply(.e$lt, function(x) {
		if(inherits(x, "ConsensusPartition")) {
			suggest_best_k(x)
		} else {
			NA
		}
	})
	leaves = all_leaves(hp)
	subgroup = rep("0", ncol(data))
	for(le in leaves) {
		if(inherits(.e$lt[[le]], "ConsensusPartition")) {
			subgroup[ .e$lt[[le]]@column_index ] = le
		} else {
			subgroup[ attr(.e$lt[[le]], "column_index") ] = le
		}
	}
	hp@subgroup = subgroup
	names(hp@subgroup) = colnames(data)

	n = length(hp@list)
	n_columns = numeric(n); names(n_columns) = names(hp@list)
	n_signatures = rep(NA, n); names(n_signatures) = names(hp@list)
	for(i in seq_len(n)) {
		if(inherits(hp@list[[i]], "ConsensusPartition")) {
			n_columns[i] = length(hp@list[[i]]@column_index)
		} else {
			n_columns[i] = length(attr(hp@list[[i]], "column_index"))
		}
		if(nodes[i] %in% leaves) {
			n_signatures[i] = NA
		} else {
			sig_tb = get_signatures(hp@list[[i]], k = suggest_best_k(hp@list[[i]]), verbose = FALSE, plot = FALSE, simplify = TRUE)
			n_signatures[i] = nrow(sig_tb)
		}
	}

	hp@n_columns = n_columns
	hp@n_signatures = n_signatures

	le = unique(as.vector(hp@hierarchy))
	col_pal = Polychrome::kelly.colors(22)
	col_pal = col_pal[!(names(col_pal) %in% c("white", "black"))]
	if(length(le) <= 20) {
		hp@subgroup_col = structure(col_pal[seq_along(le)], names = le)
	} else {
		hp@subgroup_col = structure(c(col_pal, rand_color(length(le) - length(col_pal))), names = le)
	}
	hp@call = cl
	hp@.env = hp@list[[1]]@.env

	return(hp)
}

has_hierarchy = function(object) {
	nrow(object@hierarchy) > 0
}

subgroup_dend = function(object, hierarchy = object@hierarchy) {

	check_pkg("data.tree", bioc = FALSE)

	lt = list()
	lt[["0"]] = data.tree::Node$new("0")
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
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height", edgetext = TRUE)

	dend = dendextend::`order.dendrogram<-`(dend, value = 1:nobs(dend))

	dend = edit_node(dend, function(d, index) {
		if(is.leaf(d)) {
			attr(d, "node_id") = attr(d, "label")
		} else {
			attr(d, "node_id") = attr(d, "edgetext")
		}
		d
	})

	# make sure all nodes have a node_id attribute
	# depth first
	.get_node_id = function(d) {
		node_id = attr(d, "node_id")
		if(is.null(node_id)) {
			child_node_id = .get_node_id(d[[1]])
			if(is.null(child_node_id)) {
				child_node_id = .get_node_id(d[[2]])
			}
			return(gsub("\\d$", "", child_node_id))
		}
		return(node_id)
	}

	dend = edit_node(dend, function(d, index) {
		node_id = attr(d, "node_id")
		if(is.null(node_id)) {
			node_id = .get_node_id(d)
			attr(d, "node_id") = node_id
		}
		d
	})

	# od = structure(1:length(order.dendrogram(dend)), names = labels(dend))
	# dend_env = new.env(parent = emptyenv())
	# dend_env$dend = dend

	# update_midpoint = function(index = NULL) {
	# 	if(is.null(index)) {
	# 		if(is.leaf(dend_env$dend)) {
	# 			pos = od[attr(dend_env$dend, "label")]
	# 			midpoint = 0
	# 		} else {
	# 			x = NULL
	# 			for(i in seq_len(length(dend_env$dend))) {
	# 				if(is.null(attr(dend_env$dend[[i]], "x"))) {
	# 					update_midpoint(i)
	# 				}
	# 				x[i] = attr(dend_env$dend[[i]], "x")
	# 			}
	# 			pos = (max(x) + min(x))/2
	# 			midpoint = (max(x) - min(x))/2
	# 		}
	# 	} else {
	# 		if(is.leaf(dend_env$dend[[index]])) {
	# 			pos = od[attr(dend_env$dend[[index]], "label")]
	# 			midpoint = 0
	# 		} else {
	# 			x = NULL
	# 			for(i in seq_len(length(dend_env$dend[[index]]))) {
	# 				if(is.null(attr(dend_env$dend[[c(index, i)]], "x"))) {
	# 					update_midpoint(c(index, i))
	# 				}
	# 				x[i] = attr(dend_env$dend[[c(index, i)]], "x")
	# 			}
	# 			pos = (max(x) + min(x))/2
	# 			midpoint = (max(x) - min(x))/2
	# 		}
	# 	}
	# 	if(is.null(index)) {
	# 		attr(dend_env$dend, "x") = pos
	# 	} else {
	# 		attr(dend_env$dend[[index]], "x") = pos
	# 		attr(dend_env$dend[[index]], "midpoint") = midpoint
	# 	}
	# }
	# update_midpoint()

	oe = try(dend_tmp <- as.dendrogram(as.hclust(dend)), silent = TRUE)

	if(!inherits(oe, "try-error")) {
		dend = edit_node(dend, function(d, ind) {
			if(length(ind) == 0) {
				attr(d, "midpoint") = attr(dend_tmp, "midpoint")
			} else {
				attr(d, "midpoint") = attr(dend_tmp[[ind]], "midpoint")
			}
			d
		})
	}

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
	dend
}

random_dend = function(n) {
	x = rnorm(n)
	dend = as.dendrogram(hclust(dist(1:n)))
	# set height to zero

	dendrapply(dend, function(x) {attr(x, "height") = 0; x})
}

zero_height_dend = function(n) {
	check_pkg("data.tree", bioc = FALSE)

	lt = data.tree::Node$new("foo")
	lt$node_height = 0
	for(i in 1:n) {
		lt$AddChildNode({
			node = data.tree::Node$new(paste0("foo", i))
			node$node_height = 0
			node
		})
	}
	dend = as.dendrogram(lt, heightAttribute = "node_height")
	
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
	as.numeric(names(which.min(tb)))
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
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
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
# get_classes(golub_cola_rh, depth = 2)
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
# data(golub_cola_rh)
# golub_cola_rh
setMethod(f = "show",
	signature = "HierarchicalPartition",
	definition = function(object) {

	qqcat("A 'HierarchicalPartition' object with '@{object@list[[1]]@top_value_method}:@{object@list[[1]]@partition_method}' method.\n")
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{ncol(object@.env$data)} columns.\n")
	qqcat("  Performed in total @{object@list[[1]]@n_partition*length(object@list)} partitions.\n")
	if(has_hierarchy(object)) {
		qqcat("  There are @{length(all_leaves(object))} groups.\n")
	}
	cat("\n")

	if(has_hierarchy(object)) {
		cat("Hierarchy of the partition:\n")
		hierarchy = object@hierarchy
		nodes = hierarchy[, 2]
		nc = nchar(nodes)
		names(nc) = nodes
		n = length(nc)

		parent = structure(hierarchy[, 1], names = hierarchy[, 2])
		all_leaves = all_leaves(object)
		n_columns = object@n_columns
		n_signatures = object@n_signatures

		lines = character(n)
		for(i in seq_len(n)) {
			
			lines[i] = paste0("  ", strrep("    ", nc[i] - 2), ifelse(grepl("0$", nodes[i]), "`-", "|-") ,"- ", nodes[i], qq(", @{n_columns[nodes[i]]} cols"), ifelse(is.na(n_signatures[nodes[i]]), "", qq(", @{n_signatures[nodes[i]]} signatures")))
			p = nodes[i]
			while(p != "0") {
				p = parent[p]
				if(!grepl("0$", p)) {
					substr(lines[i], (nc[p] - 2)*4+3, (nc[p] - 2)*4+3) = "|"
				}
			}
		}
		# substr(lines[1], 1, 1) = "+"
		lines = c(qq("  0, @{ncol(object@list[['0']])} cols"), lines)
		cat(lines, sep = "\n")
	} else {
		cat("No hierarchy found.\n")
	}
	
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
	definition = function(object, depth = max_depth(object),
	group_diff = cola_opt$group_diff,
	row_km = NULL,
	scale_rows = object[1]@scale_rows, 
	anno = get_anno(object), 
	anno_col = get_anno_col(object),
	show_column_names = FALSE, column_names_gp = gpar(fontsize = 8),
	verbose = TRUE, plot = TRUE, seed = 888,
	...) {

	if(depth <= 1) {
		stop_wrap("depth should be at least larger than 1.")
	}

	alf = all_leaves(object, depth = depth)
	ap = setdiff(all_nodes(object, depth = depth), all_leaves(object, depth = depth))

	sig_lt = list()
	for(p in ap) {
		best_k = suggest_best_k(object[[p]])
		if(verbose) qqcat("* get signatures at node @{p} with @{best_k} subgroups.\n")
		sig_tb = get_signatures(object[[p]], k = best_k, prefix = "  ", verbose = TRUE, plot = FALSE, simplify = TRUE, ...)
		
		sig_lt[[p]] = sig_tb
		# if(verbose) qqcat("  * find @{nrow(sig_tb)} signatures at node @{p}\n")
	}

	all_index = sort(unique(unlist(lapply(sig_lt, function(x) x[, 1]))))

	returned_df = data.frame(which_row = all_index)

	# filter by group_diff
	mat = object@.env$data[all_index, , drop = FALSE]
	class = get_classes(object, depth = depth)[, 1]

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

	dend = cluster_within_group(use_mat1, class)
	
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
# -method Method to visualize.
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
	definition = function(object, depth = max_depth(object), 
	method = c("euler", "upset"), verbose = interactive(), ...) {

	lt = object@list
	al = all_leaves(object)
	lt = lt[! names(lt) %in% al]
	lt = lt[nchar(names(lt)) <= depth]

	sig_list = lapply(lt, function(x) {
		tb = get_signatures(x, k = suggest_best_k(x), verbose = verbose, ..., plot = FALSE)
		if(is.null(tb)) {
			return(integer(0))
		} else {
			return(tb$which_row)
		}
	})

	l = sapply(sig_list, length) > 0
	if(any(l) && verbose) {
		qqcat("Following nodes have no signature found: \"@{paste(names(sig_list)[l], collapse=', ')}\"\n")
	}
	sig_list = sig_list[l]

	if(missing(method)) {
		if(length(sig_list) <= 3) {
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
		if(length(comb_size(m)) > 40) {
			m = m[order(comb_size(m), decreasing = TRUE)[1:40]]
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
	definition = function(object, depth = max_depth(object), 
	show_row_names = FALSE, row_names_gp = gpar(fontsize = 8),
	anno = get_anno(object[1]), anno_col = get_anno_col(object[1])) {

	cl = get_classes(object, depth = depth)[, 1]
	dend = calc_dend(object, depth = depth)

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
# Subset the HierarchicalPartition object
#
# == param
# -x A `HierarchicalPartition-class` object.
# -i Index. The value should be numeric or a node ID.
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
# data(golub_cola_rh)
# golub_cola_rh["01"]
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
# -x A `HierarchicalPartition-class` object
# -i Index. The value should be numeric or a node ID.
#
# == details
# On each node, there is a `ConsensusPartition-class` object.
#
# Note you cannot get a sub-hierarchy of the partition.
#
# == value
# A `ConsensusPartition-class` object.
#
"[[.HierarchicalPartition" = function(x, i) {
	x[i]
}


# == title
# Test correspondance between predicted classes and known factors
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
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
	colnames(m) = paste0(colnames(m), "(p)")
	df = cbind(n = nrow(class), m, n_class = apply(class, 2, function(x) length(unique(x))))
	return(df)
})

# == title
# Visualize columns after dimension reduction
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth of the hierarchy.
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
	definition = function(object,
	depth = max_depth(object), parent_node,
	top_n = NULL, method = c("PCA", "MDS", "t-SNE", "UMAP"),
	scale_rows = TRUE, verbose = TRUE, ...) {

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
# Max depth of the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
#
# == value
# A numeric value.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# max_depth(golub_cola_rh)
setMethod(f = "max_depth",
	signature = "HierarchicalPartition",
	definition = function(object) {

	if(has_hierarchy(object)) {
		max(nchar(object@hierarchy[, 2]))
	} else {
		1
	}
})

# == title
# All nodes in the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth in the hierarchy.
#
# == value
# A vector of node ID.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# all_nodes(golub_cola_rh)
setMethod(f = "all_nodes",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object)) {

	if(has_hierarchy(object)) {
		all_nodes = unique(as.vector(t(object@hierarchy)))
		if(!is.null(depth)) {
			all_nodes = all_nodes[nchar(all_nodes) <= depth]
		}
		return(all_nodes)
	} else {
		return(character(0))
	}
})

# == title
# All leaves in the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
# -depth Depth in the hierarchy.
#
# == value
# A vector of node ID.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# all_leaves(golub_cola_rh)
setMethod(f = "all_leaves",
	signature = "HierarchicalPartition",
	definition = function(object, depth = max_depth(object)) {

	if(has_hierarchy(object)) {
		hierarchy = unique(object@hierarchy)
		if(!is.null(depth)) {
			hierarchy = hierarchy[nchar(hierarchy[, 2]) <= depth, , drop = FALSE]
		}
		
		tb = table(hierarchy)
		tb = tb[tb <= 1]
		names(tb)	
	} else {
		"0"
	}
})

# == title
# Test whether a node is a leaf node
#
# == param
# -x A `HierarchicalPartition-class` object.
# -node A node ID.
#
is_leaf_node = function(x, node) {
	node %in% all_leaves(x)
}

get_children = function(object, node = "0") {
	hierarchy = unique(object@hierarchy)
	hierarchy[hierarchy[, 1] == node, 2]
}

# == title
# Get the original matrix
#
# == param
# -object A `HierarchicalPartition-class` object.
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
# -object A `HierarchicalPartition-class` object.
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
# -object A `HierarchicalPartition-class` object.
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
	for(i in seq_along(object@list)) {
		obj = object@list[[i]]
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
	}

	tb = data.frame(
		node = names(object@list),
		is_leaf = names(object@list) %in% all_leaves(object),
		best_k = best_k,
		"1-PAC" = stability,
		mean_silhouette = mean_silhouette,
		concordance = concordance,
		check.names = FALSE)
	tb$n_sample = sapply(tb$node, function(id) {
		ncol(object[[id]])
	})

	rntb = rownames(tb)
	l = tb$`1-PAC` >= 0.9 & !is.na(tb$best_k)

	tb = cbind(tb, ifelse(l, ifelse(tb$`1-PAC` <= 0.95, "*", "**"), ""), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""
	
	# l = tb[, 1] %in% all_leaves(object)
	# tb = cbind(tb, ifelse(l, "leaf", "node"), stringsAsFactors = FALSE)
	# colnames(tb)[ncol(tb)] = ""

	
	stop_reason = lapply(object@list, function(obj) {
		attr(obj, "stop_reason")
	})
	attr(tb, "stop_reason") = stop_reason

	class(tb) = c("hc_table_suggest_best_k", class(tb))
	tb
})

stop_reason_index = c(
	"Rand indices for all k were too high." = "z",
	"PAC score was too big." = "a",
	"Subgroup had too few columns." = "b",
	"There were too few signatures." = "c"
)

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
			return(stop_reason_index[x])
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
		cat("  ", a, ") ", names(which(stop_reason_index == a)), "\n", sep = "")
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
	title = qq("cola Report for Hierarchical Partitioning (@{object[1]@top_value_method}:@{object[1]@partition_method})"), 
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
# Number of rows in the matrix
#
# == param
# -x A `HierarchicalPartition-class` object.
#
setMethod(f = "nrow",
	signature = "HierarchicalPartition",
	definition = function(x) {
	nrow(x@.env$data)
})

# == title
# Number of columns in the matrix
#
# == param
# -x A `HierarchicalPartition-class` object.
#
setMethod(f = "ncol",
	signature = "HierarchicalPartition",
	definition = function(x) {
	ncol(x@.env$data)
})

# == title
# Row names of the matrix
#
# == param
# -x A `HierarchicalPartition-class` object.
#
setMethod(f = "rownames",
	signature = "HierarchicalPartition",
	definition = function(x) {
	rownames(x@.env$data)
})

# == title
# Column names of the matrix
#
# == param
# -x A `HierarchicalPartition-class` object.
#
setMethod(f = "colnames",
	signature = "HierarchicalPartition",
	definition = function(x) {
	colnames(x@.env$data)
})

# == title
# Dimension of the matrix
#
# == param
# -x A `HierarchicalPartition-class` object.
#
dim.HierarchicalPartition = function(x) {
	dim(x@.env$data)
}
