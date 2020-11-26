
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
        hierarchy = "matrix",
        subgroup = "character",
        subgroup_col = "character",
        subgroup_dend = "ANY",
        node_level = "list",
        call = "ANY",
        .env = "environment"
    )
)

# == title
# Hierarchical partition
#
# == param
# -data a numeric matrix where subgroups are found by columns.
# -top_n Number of rows with top values. Since the number of rows of sub-matrieces might be different,
#       here ``top_n`` can be set as a vector of values less than 1 which are treated as the fraction
#       of the rows from the input matrix.
# -top_value_method a single top-value method. Available methods are in `all_top_value_methods`.
# -partition_method a single partition method. Available methods are in `all_partition_methods`.
# -combination_method A list of combinations of top-value methods and partitioning methods. The value
#     can be a two-column data frame where the first column is the top-value methods and the second
#     column is the partitioning methods. Or it can be a vector of combination names in a form of
#     "top_value_method:partitioning_method".
# -anno A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation
#       to predicted subgroups will be tested.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
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
	top_n = seq(min(1000, round(nrow(data)*0.1)), 
		        min(3000, round(nrow(data)*0.3)), 
		        length.out = 3),
	top_value_method = c("SD", "ATC"), 
	partition_method = c("kmeans", "skmeans"),
	combination_method =  expand.grid(top_value_method, partition_method),
	anno = NULL, anno_col = NULL,
	PAC_cutoff = 0.2, min_samples = max(6, round(ncol(data)*0.01)), subset = Inf,
	min_n_signatures = round(nrow(data)*min_p_signatures), 
	min_p_signatures = 0.01,
	max_k = 4, verbose = TRUE, mc.cores = 1, help = TRUE, ...) {

	t1 = Sys.time()

	data = as.matrix(data)

	if(help) {
		if(identical(subset, Inf) && ncol(data) > 500) {
			qqcat_wrap("You have quite a lot of columns in the matrix. For reducing the runtime, you can set `subset` argument to a number less than the total number of columns or a subset of column indices. The classification of unselected columns are inferred from the classes of the selected columns. Set the argument 'help = FALSE' to turn off this message. Other tips: 1. set a single value for `top_value_method` and `partition_method`, 2. set a single value for `top_n`.\n")
			cat("\n")
		}
	}

	cl = match.call()

	if(min_samples < 3) {
		if(verbose) qqcat("! 'min_samples' was reset to 3.\n")
		min_samples = 3
	}

	check_pkg("Polychrome", bioc = FALSE)

	if(verbose) {
		qqcat("* hierarchical partition on a @{nrow(data)}x@{ncol(data)} matrix.\n")
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

	if(verbose) {
		if(length(combination_method) == 1) {
			qqcat("* running @{combination_method[[1]][1]}:@{combination_method[[1]][2]}.\n")
		} else {
			qqcat("* running @{length(combination_method)} combinations of top-value methods and partitioning methods.\n")
		}
	}

	if(!is.null(anno)) {
		if(is.atomic(anno)) {
			known_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = known_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = known_nm
			}
		}
		if(nrow(anno) != ncol(data)) {
			stop_wrap("nrow of `anno` should be the same as ncol of the matrix.")
		}
	}

	if(is.null(anno_col)) {
		anno_col = lapply(anno, ComplexHeatmap:::default_col)
	} else {
		if(ncol(anno) == 1 && is.atomic(anno_col)) {
			anno_col = list(anno_col)
			names(anno_col) = colnames(anno)
		} else if(is.null(names(anno_col))) {
			if(length(anno_col) == ncol(anno)) {
				names(anno_col) = colnames(anno)
			} else {
				anno_col = lapply(anno, ComplexHeatmap:::default_col)
			}
		}
		for(nm in names(anno)) {
			if(is.null(anno_col[[nm]])) {
				anno_col[[nm]] = ComplexHeatmap:::default_col(anno[[nm]])
			}
		}
	}
	if(is.null(anno)) {
		anno_col = NULL
	}

	# if(!multicore_supported()) {
	# 	if(mc.cores > 1) message("* mc.cores is reset to 1 because mclapply() is not supported on this OS.")
	# 	mc.cores = 1
	# }
	
	.hierarchical_partition = function(.env, column_index, node_id = '0', 
		min_samples = 6, max_k = 4, verbose = TRUE, mc.cores = 1, ...) {

		prefix = ""
		if(node_id != "0") {
			prefix = paste(rep("  ", nchar(node_id) - 1), collapse = "")
		}

		if(verbose) qqcat("@{prefix}================== node @{node_id} ============================\n")
		if(verbose) qqcat("@{prefix}* submatrix with @{length(column_index)} columns, node_id: @{node_id}.\n")

		if(length(column_index) < 2*min_samples) {
			if(verbose) qqcat("@{prefix}* number of samples is not enough to perform partitioning.\n")
			part = STOP_REASON["b"]
			# we need the following two values for other functions
			attr(part, "node_id") = node_id
			attr(part, "column_index") = column_index
			attr(part, "stop_reason") = STOP_REASON["b"]
			return(list(obj = part))
		}

		## all_top_value_list is only used in run_all_consensus_partition_methods(), we remove it here
	   	
	   	if(node_id == "0") {
	   		if(all(top_n <= 1)) {
				top_n = floor(top_n * nrow(data))
				top_n = top_n[top_n > 0]
			}
	   	}
		if(node_id != "0") {
			row_sd = rowSds(data[, column_index, drop = FALSE])
			qa = quantile(unique(row_sd[row_sd > 1e-10]), 0.95, na.rm = TRUE)*0.05
			l = row_sd > qa
			.env$row_index = which(l)
			if(verbose) qqcat("@{prefix}* @{sum(!l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")

			if(all(top_n <= 1)) {
				top_n = floor(top_n * length(.env$row_index))
				top_n = top_n[top_n > 0]
			}

			if(length(.env$row_index) < min(top_n)*1.2 || length(top_n) == 0) {
				if(verbose) qqcat("@{prefix}* number of rows is not enough to perform partitioning.\n")
				part = STOP_REASON["d"]
				# we need the following two values for other functions
				attr(part, "node_id") = node_id
				attr(part, "column_index") = column_index
				attr(part, "stop_reason") = STOP_REASON["d"]
				return(list(obj = part))
			}
		}
		if(mc.cores > 1 && verbose) {
			qqcat("@{prefix}* running consensus partitioning with @{mc.cores} cores.\n")
		}

		if(is.null(anno)) {
			anno2 = NULL
		} else {
			anno2 = anno[column_index, , drop = FALSE]
		}

		part_list = list(length(combination_method))
		if(length(column_index) <= subset) {
			if(node_id != "0") .env$all_top_value_list = NULL
		
			for(i in seq_along(combination_method)) {
				if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
				.env$column_index = column_index #note .env$column_index is only for passing to `consensus_partition()` function
				part_list[[i]] = consensus_partition(verbose = TRUE, .env = .env, max_k = max_k, prefix = prefix,
					top_n = top_n, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					mc.cores = mc.cores, anno = anno2, anno_col = anno_col, ...)
			}
		} else {
			.env$all_top_value_list = NULL

			# in consensus_partition_by_down_sampling(), .env$all_top_value_list is always reset to NULL
			# because the columns are randomly sampled and top_value changes. However, here we cache
			# the recent top_value for a top_value_method for downstream process
			all_top_value_list_ds = list()

			for(i in seq_along(combination_method)) {
				if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
				.env$column_index = column_index #note .env$column_index is only for passing to `consensus_partition()` function
				part_list[[i]] = consensus_partition_by_down_sampling(subset = subset, verbose = TRUE, .env = .env, max_k = max_k, prefix = prefix,
					top_n = top_n, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					mc.cores = mc.cores, .predict = FALSE, anno = anno2, anno_col = anno_col, ...)
				all_top_value_list_ds[[part_list[[i]]@top_value_method]] = .env$all_top_value_list[[part_list[[i]]@top_value_method]]
			}
		}
	   
		if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")

		stat_tb = sapply(part_list, function(part) {
			stat_df = get_stats(part, all_stats = TRUE)
			best_k = suggest_best_k(part)
			if(is.na(best_k)) {
				c(0, 0, 0, 0)
			} else {
				stat_df[stat_df[, "k"] == best_k, c("1-PAC", "mean_silhouette", "concordance", "area_increased")]
			}
		})
		stat_tb = as.data.frame(t(stat_tb))
		ind = do.call(order, -stat_tb)[1]
		part = part_list[[ind]]

		if(verbose) qqcat("@{prefix}* select @{part@top_value_method}:@{part@partition_method} with the most stability among methods.\n")

		dist_method = list(...)$dist_method
		if(is.null(dist_method)) dist_method = "euclidean"
		if(length(column_index) > subset) {

			part = convert_to_DownSamplingConsensusPartition(part, column_index, dist_method, verbose, prefix, mc.cores)
		}

		attr(part, "node_id") = node_id

		lt = list(obj = part)

	    best_k = suggest_best_k(part)
	    if(is.na(best_k)) {
	    	attr(lt$obj, "stop_reason") = STOP_REASON["z"]
	    	if(verbose) qqcat("@{prefix}* Rand index is too high, no meaningful subgroups, stop.\n")
	    	return(lt)
	    }
	    cl = get_classes(part, k = best_k)

	    # check PAC score
	    PAC_score = 1 - get_stats(part, k = best_k)[, "1-PAC"]
	    if(PAC_score > PAC_cutoff) {
	    	if(verbose) qqcat("@{prefix}* PAC score is too big (@{sprintf('%.2f', PAC_score)}), stop.\n")
	    	attr(lt$obj, "stop_reason") = STOP_REASON["a"]
	    	return(lt)
	    }

	    # check number of confident samples
	    # if(length(column_index) <= subset) {
	    # 	class_df = get_classes(part, k = best_k)
	    # 	p_confident_samples = sum(class_df[, "silhouette"] > 0.5)/nrow(class_df)
	    # } else {
	    # 	class_df = get_classes(part, k = best_k)
	    # 	p_confident_samples = sum(class_df[, "p"] < 0.05)/nrow(class_df)
	    # }
	    # if(p_confident_samples < 0.8) {
	    # 	if(verbose) qqcat("@{prefix}* Proportion of confident columns is too small (@{sprintf('%.2f', p_confident_samples)}), stop.\n")
	    # 	attr(lt$obj, "stop_reason") = STOP_REASON["e"]
	    # 	return(lt)
	    # }

	    dd = data[.env$row_index, column_index, drop = FALSE]
	    if(length(column_index) <= subset) {
	    	ri = order(.env$all_top_value_list[[part@top_value_method]], decreasing = TRUE)[1:max(get_param(part)[, "top_n"])]
	    } else {
	    	ri = order(all_top_value_list_ds[[part@top_value_method]], decreasing = TRUE)[1:max(get_param(part)[, "top_n"])]
		}
		if(length(ri) > 5000) ri = sample(ri, 5000)
		mean_dist = tapply(seq_len(ncol(dd)), cl$class, function(ind) {
			n = length(ind)
			if(n == 1) {
				return(Inf)
			}
			sum(dist(t(dd[ri, ind, drop = FALSE]))^2)/(n*(n-1)/2)
		})

	    kg1 = as.numeric(names(which.min(mean_dist)))

	    ## test the number of samples
	    sample_too_small = FALSE
	    if(best_k == 2) {

	    	kg2 = setdiff(1:2, kg1)

	    	set1 = which(cl$class %in% kg1)
	    	set2 = which(cl$class %in% kg2)

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
	    	for(ik in c(kg1, setdiff(sort(as.numeric(names(tb))), kg1))) {
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

	    .env$all_top_value_list = NULL

    	# if(sample_too_small) {
    	# 	if(verbose) qqcat("@{prefix}* some of the subgroups have too few columns (< @{min_samples}), won't split.\n")
    	# 	attr(lt$obj, "stop_reason") = STOP_REASON["b"]
    	# 	return(lt)
    	# }

    	# check the numbers of signatures
    	if(verbose) qqcat("@{prefix}* checking number of signatures in the best classification.\n")
    	if(length(column_index) <= subset) {
    		sig_df = get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE, simplify = TRUE, silhouette_cutoff = -Inf)
    	} else {
    		sig_df = get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE, simplify = TRUE, p_cutoff = Inf)
    	}
    	if(is.null(part@.env$signature_hash)) {
    		part@.env$signature_hash = list()
    	}
    	part@.env$signature_hash[[node_id]] = attr(sig_df, "hash")
    	
    	n_sig = nrow(sig_df)
    	p_sig = n_sig/nrow(part)
    	if(n_sig <= min_n_signatures && p_sig <= min_p_signatures) {
    		if(verbose) qqcat("@{prefix}* too few signatures (n = @{n_sig}, p = @{sprintf('%.3f', p_sig)}) under the best classification, stop.\n")
    		attr(lt$obj, "stop_reason") = STOP_REASON["c"]
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

	if(verbose) cat("* calculate top-values.\n")
	all_top_value_method = unique(sapply(combination_method, function(x) x[1]))
	all_top_value_list = lapply(all_top_value_method, function(tm) {
		if(verbose) qqcat("  - calculate @{tm} score for @{nrow(data)} rows.\n")
		all_top_value = get_top_value_method(tm)(data)
		all_top_value[is.na(all_top_value)] = -Inf
		return(all_top_value)
	})
	names(all_top_value_list) = all_top_value_method
	.env$all_top_value_list = all_top_value_list
	.env$node_0_top_value_list = .env$all_top_value_list
   	
	lt = .hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), min_samples = min_samples, 
		node_id = "0", max_k = min(max_k, ncol(data)-1), verbose = verbose, mc.cores = mc.cores, ...)

	qqcat("* formatting the results into a HierarchicalPartition object.\n")

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
	hp@node_level$best_k = sapply(.e$lt, function(x) {
		if(inherits(x, "ConsensusPartition")) {
			suggest_best_k(x)
		} else {
			NA
		}
	})
	if(nrow(hp@hierarchy) > 0) {
		tb = table(hp@hierarchy)
		leaves = names(tb[tb <= 1])
	} else {
		leaves = "0"
	}

	subgroup = rep("0", ncol(data))
	for(le in leaves) {
		if(inherits(.e$lt[[le]], "DownSamplingConsensusPartition")) {
			subgroup[ .e$lt[[le]]@full_column_index ] = le
		} else if(inherits(.e$lt[[le]], "ConsensusPartition")) {
			subgroup[ .e$lt[[le]]@column_index ] = le
		} else {
			subgroup[ attr(.e$lt[[le]], "column_index") ] = le
		}
	}
	hp@subgroup = subgroup
	names(hp@subgroup) = colnames(data)

	n = length(hp@list)
	n_columns = numeric(n); names(n_columns) = names(hp@list)
	n_signatures = rep(NA_real_, n); names(n_signatures) = names(hp@list)
	nodes = names(hp@list)
	for(i in seq_len(n)) {
		if(inherits(hp@list[[i]], "DownSamplingConsensusPartition")) {
			n_columns[i] = length(hp@list[[i]]@full_column_index)
		} else if(inherits(hp@list[[i]], "ConsensusPartition")) {
			n_columns[i] = length(hp@list[[i]]@column_index)
		} else {
			n_columns[i] = length(attr(hp@list[[i]], "column_index"))
		}
		if(nodes[i] %in% leaves) {
			n_signatures[i] = NA_real_
		} else {
			sig_tb = get_signatures(hp@list[[i]], k = suggest_best_k(hp@list[[i]]), verbose = FALSE, plot = FALSE, simplify = TRUE)
			n_signatures[i] = nrow(sig_tb)
		}
	}

	hp@node_level$n_columns = n_columns
	hp@node_level$n_signatures = n_signatures
	hp@node_level$p_signatures = n_signatures/nrow(data)

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

	if(has_hierarchy(hp)) {
		hp@subgroup_dend = subgroup_dend(hp)
		hp@node_level$node_height = structure(dend_node_apply(hp@subgroup_dend, function(d) attr(d, "height")),
			                            names = dend_node_apply(hp@subgroup_dend, function(d) attr(d, "node_id")))
	}

	hp@.env$combination_methods = combination_method

	t2 = Sys.time()
	if(verbose) cat("* totally used ", gsub("^ +", "", format(t2 - t1)), ".\n", sep = "")

	return(hp)
}


STOP_REASON = c(
	"z" = "Rand indices for all k were too high.",
	"a" = "PAC score was too big.",
	"b" = "Subgroup had too few columns.",
	"c" = "There were too few signatures.",
	"d" = "Matrix has too few rows (less than 'top_n').",
	"e" = "There were too few conficent columns to detect signatures."
)

STOP_REASON_INDEX = structure(names(STOP_REASON), names = unname(STOP_REASON))

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
		n_columns = object@node_level$n_columns
		n_signatures = object@node_level$n_signatures

		lines = character(n)

		si = NULL
		for(i in seq_len(n)) {
			
			lines[i] = paste0("  ", strrep("    ", nc[i] - 2), ifelse(grepl("0$", nodes[i]), "`-", "|-") ,"- ", nodes[i], qq(", @{n_columns[nodes[i]]} cols"))
			if(!is.na(n_signatures[nodes[i]])) {
				lines[i] = paste0(lines[i], qq(", @{n_signatures[nodes[i]]} signatures"))
			}
			stop_reason = attr(object@list[[ nodes[i] ]], "stop_reason")
			if(!is.null(stop_reason)) {
				lines[i] = paste0(lines[i], qq(" (@{STOP_REASON_INDEX[stop_reason]})"))
				si = c(si, STOP_REASON_INDEX[stop_reason])
			}

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

		if(length(si)) {
			si = sort(unique(si))
			cat("Stop reason:\n")
			for(s in si) {
				cat("  ", s, ") ", names(which(STOP_REASON_INDEX == s)), "\n", sep = "")
			}
		}
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

# == title
# Information on the nodes
#
# == param
# -object A `HierarchicalPartition-class` object.
#
# == details
# It returns the following node-level information:
#
# -id Node id.
# -n_columns Number of columns.
# -n_signatures Number of signatures.
# -p_signatures Percent of signatures.
# -node_height Height of the node in the dendrogram
# -is_leaf Whether the node is a leaf
#
setMethod(f = "node_info",
	signature = "HierarchicalPartition",
	definition = function(object) {

	df = do.call(cbind, object@node_level)
	best_method = sapply(object@list, function(x) {
		if(inherits(x, "ConsensusPartition")) {
			paste0(x@top_value_method, ":", x@partition_method)
		} else {
			"not applied"
		}
	})
	df = data.frame(id = rownames(df), 
		best_method = best_method,
		depth = nchar(rownames(df)), df)
	df$is_leaf = FALSE
	df$is_leaf[df$id %in% all_leaves(object)] = TRUE
	rownames(df) = NULL
	df
})


