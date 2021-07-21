
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
# -`functional_enrichment,HierarchicalPartition-method`: apply functional enrichment.
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
        param = "list",
        running_time = "ANY",
        call = "ANY",
        .env = "environment"
    )
)

# == title
# Hierarchical partition
#
# == param
# -data a numeric matrix where subgroups are found by columns.
# -top_n Number of rows with top values.
# -top_value_method a single or a vector of top-value methods. Available methods are in `all_top_value_methods`.
# -partition_method a single or a vector of partition methods. Available methods are in `all_partition_methods`.
# -combination_method A list of combinations of top-value methods and partitioning methods. The value
#     can be a two-column data frame where the first column is the top-value methods and the second
#     column is the partitioning methods. Or it can be a vector of combination names in a form of
#     "top_value_method:partitioning_method".
# -anno A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation
#       to predicted subgroups will be tested.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -mean_silhouette_cutoff The cutoff to test whether partition in current node is stable.
# -min_samples the cutoff of number of samples to determine whether to continue looking for subgroups.
# -group_diff Pass to `get_signatures,ConsensusPartition-method`.
# -fdr_cutoff Pass to `get_signatures,ConsensusPartition-method`.
# -subset Number of columns to randomly sample.
# -min_n_signatures Minimal number of signatures under the best classification.
# -filter_fun A self-defined function which filters the original matrix and returns a submatrix for partitioning.
# -max_k maximal number of partitions to try. The function will try ``2:max_k`` partitions. Note this is the number of
#        partitions that will be tried out on each node of the hierarchical partition. Since more subgroups will be found
#        in the whole partition hierarchy, on each node, ``max_k`` should not be set to a large value.
# -scale_rows Whether rows are scaled?
# -verbose whether print message.
# -mc.cores multiple cores to use. This argument will be removed in future versions.
# -cores Number of cores, or a ``cluster`` object returned by `parallel::makeCluster`.
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
	top_n = NULL,
	top_value_method = "ATC", 
	partition_method = "skmeans",
	combination_method =  expand.grid(top_value_method, partition_method),
	anno = NULL, anno_col = NULL,
	mean_silhouette_cutoff = 0.9, min_samples = max(6, round(ncol(data)*0.01)), subset = Inf,
	group_diff = ifelse(scale_rows, 0.5, 0), 
	fdr_cutoff = cola_opt$fdr_cutoff,
	min_n_signatures = NULL, 
	filter_fun = function(mat) {
		s = rowSds(mat)
		s > quantile(unique(s[s > 1e-10]), 0.05, na.rm = TRUE)
	},
	max_k = 4, scale_rows = TRUE, verbose = TRUE, mc.cores = 1, cores = mc.cores, help = TRUE, ...) {

	t1 = Sys.time()

	data = as.matrix(data)

	if(any(rowSds(data) == 0)) {
		stop_wrap("Deteched some rows have zero SD, please remove them.")
	}

	if(help) {
		if(identical(subset, Inf) && ncol(data) > 500) {
			qqcat_wrap("You have quite a lot of columns in the matrix. For reducing the runtime, you can set `subset` argument to a number less than the total number of columns or a subset of column indices. The classification of unselected columns are inferred from the classes of the selected columns. Set the argument 'help = FALSE' to turn off this message. Other tips: 1. set a single value for `top_value_method` and `partition_method`, 2. set a single value for `top_n`.\n")
			cat("\n")
		}
	}

	cl = match.call()
	rg = range(data)
	if(rg[1] >= 0 && rg[2] <= 1) {
		if( (!"top_value_method" %in% names(cl)) && (!"partition_method" %in% names(cl)) ) {
			top_value_method = "SD"
			partition_method = "kmeans"
			if(help) {
				qqcat_wrap("The range of the matrix is inside c(0, 1), thus, `top_value_method` is set to 'SD' and `partition_method` is set to 'kmeans'. You can stop this default behaviour by explictly setting values for these two arguments.")
				cat("\n")
			}
		}
	}

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

	combination_method = lapply(combination_method, unname)

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

	cores = get_nc(cores)

	.env = new.env(parent = emptyenv())
	.env$data = data
	total_n = ncol(data)

	if(verbose) cat("* calculate top-values.\n")
	all_top_value_method = unique(sapply(combination_method, function(x) x[1]))
	

	# all_top_value_list = lapply(all_top_value_method, function(tm) {
	# 	if(verbose) qqcat("  - calculate @{tm} score for @{nrow(data)} rows.\n")
	# 	if(tm == "ATC") {
	# 		if(identical(get_top_value_method("ATC"), ATC)) {
	# 			config_ATC(cores = cores)
	# 			if(verbose && cores > 1) qqcat("  - set @{cores} cores for ATC()\n")
	# 		}
	# 	}

	# 	all_top_value = get_top_value_method(tm)(data)
	# 	all_top_value[is.na(all_top_value)] = -Inf
	# 	return(all_top_value)
	# })
	# names(all_top_value_list) = all_top_value_method
	# .env$all_top_value_list = all_top_value_list
	.env$combination_method = combination_method

	.env$global_row_mean = rowMeans(data)
	.env$global_row_sd = rowSds(data)

	.env$group_diff = group_diff
	.env$fdr_cutoff = fdr_cutoff

	lt = .hierarchical_partition(.env = .env, column_index = seq_len(ncol(data)), subset = subset, anno = anno, anno_col = anno_col,
		min_samples = min_samples, node_id = "0", max_k = min(max_k, ncol(data)-1), verbose = verbose, cores = cores, mean_silhouette_cutoff = mean_silhouette_cutoff,
		top_n = top_n, min_n_signatures = min_n_signatures, group_diff = group_diff, fdr_cutoff = fdr_cutoff, scale_rows = scale_rows, 
		filter_fun = filter_fun, ...)

	if(verbose) qqcat("* formatting the results into a HierarchicalPartition object.\n")

	# reformat lt
	.e = new.env(parent = emptyenv())
	.e$hierarchy = matrix(nrow = 0, ncol = 2)
	.e$lt = list()
	reformat_lt = function(lt, .e) {
		nm = names(lt)
		parent_id = attr(lt$obj, "node_id")
		.e$lt[[parent_id]] = lt$obj

		for(nm in names(lt)) {
			if(grepl("^child\\d+$", nm)) {
				child_id = attr(lt[[nm]]$obj, "node_id")
				.e$hierarchy = rbind(.e$hierarchy, c(parent_id, child_id))
				reformat_lt(lt[[nm]], .e)
			}
		}
	}
	reformat_lt(lt, .e)

	hp = new("HierarchicalPartition")
	hp@hierarchy = .e$hierarchy
	hp@list = .e$lt
	hp@node_level$best_k = sapply(.e$lt, function(x) {
		if(inherits(x, "ConsensusPartition")) {
			attr(x, "best_k")
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
			if(attr(hp@list[[i]], "stop_reason") == STOP_REASON["c"]) {
				sig_tb = get_signatures(hp@list[[i]], k = attr(hp@list[[i]], "best_k"), verbose = FALSE, plot = FALSE, simplify = TRUE,
					group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = .env$global_row_mean, .scale_sd = .env$global_row_sd)
				n_signatures[i] = nrow(sig_tb)
			} else {
				n_signatures[i] = NA_real_
			}
		} else {
			sig_tb = get_signatures(hp@list[[i]], k = attr(hp@list[[i]], "best_k"), verbose = FALSE, plot = FALSE, simplify = TRUE,
				group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = .env$global_row_mean, .scale_sd = .env$global_row_sd)
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

	hp@.env$combination_methods = combination_method

	if(!is.null(.env$min_n_signatures)) min_n_signatures = .env$min_n_signatures

	hp@param = list(top_n = top_n, combination_method = combination_method, mean_silhouette_cutoff = mean_silhouette_cutoff, min_samples = min_samples,
		subset = subset, group_diff = group_diff, fdr_cutoff = fdr_cutoff, min_n_signatures = min_n_signatures, max_k = max_k, cores = cores)

	t2 = Sys.time()
	if(verbose) cat("* totally used ", gsub("^ +", "", format(t2 - t1)), ".\n", sep = "")

	hp@running_time = t2 - t1
	return(hp)
}

   
.hierarchical_partition = function(.env, column_index, node_id = '0', subset = Inf, anno = NULL, anno_col = anno_col,
	min_samples = 6, max_k = 4, verbose = TRUE, cores = 1, scale_rows = TRUE,
	filter_fun = function(mat) {
		s = rowSds(mat)
		qa = quantile(unique(s[s > 1e-10]), 0.05, na.rm = TRUE)
		s > qa
	},
	top_n = NULL, mean_silhouette_cutoff = 0.9, min_n_signatures = 100, group_diff = cola_opt$group_diff, fdr_cutoff = cola_opt$fdr_cutoff, ...) {

	prefix = ""
	if(node_id != "0") {
		prefix = paste(rep("  ", nchar(node_id) - 1), collapse = "")
	}

	if(verbose) qqcat(crayon::red("@{prefix}================== node @{node_id} ============================\n"))
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
   	data = .env$data
   	combination_method = .env$combination_method
   	global_row_mean = .env$global_row_mean
   	global_row_sd = .env$global_row_sd

   	top_n2 = top_n
	if(node_id != "0") {
		l = filter_fun(data[, column_index, drop = FALSE])
		if(is.logical(l)) {
			l[is.na(l)] = FALSE
			.env$row_index = which(l)
			if(verbose) qqcat("@{prefix}* @{sum(!l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")
		} else {
			.env$row_index = l
			if(verbose) qqcat("@{prefix}* @{nrow(data) - length(l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")
		}

		if(!is.null(top_n)) {
			if(length(top_n2) > 1) {
				min_n = min(top_n2)
				max_n = max(top_n2)
				top_n2 = (max_n - min_n) * length(column_index)/ncol(.env$data) + min_n
				top_n2 = unique(top_n2)
			}

			if(length(.env$row_index) < min(top_n2)*1.2 || length(top_n2) == 0) {
				if(verbose) qqcat("@{prefix}* number of rows is not enough to perform partitioning.\n")
				part = STOP_REASON["d"]
				# we need the following two values for other functions
				attr(part, "node_id") = node_id
				attr(part, "column_index") = column_index
				attr(part, "stop_reason") = STOP_REASON["d"]
				return(list(obj = part))
			}
		}
	} else {
		l = filter_fun(data)
		if(is.logical(l)) {
			l[is.na(l)] = FALSE
			.env$row_index = which(l)
			if(verbose) qqcat("@{prefix}* @{sum(!l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")
		} else {
			.env$row_index = l
			if(verbose) qqcat("@{prefix}* @{nrow(data) - length(l)}/@{nrow(data)} rows are removed for partitioning, due to very small variance.\n")
		}
	}
	if(cores > 1 && verbose) {
		qqcat("@{prefix}* running consensus partitioning with @{cores} cores.\n")
	}

	if(is.null(anno)) {
		anno2 = NULL
	} else {
		anno2 = anno[column_index, , drop = FALSE]
	}

	part_list = list()
	if(length(column_index) <= subset) {
		.env$all_top_value_list = NULL
	
		for(i in seq_along(combination_method)) {
			if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
			.env$column_index = column_index #note .env$column_index is only for passing to `consensus_partition()` function
			nm = qq("@{combination_method[[i]][1]}:@{combination_method[[i]][2]}-row")
			part_list[[nm]] = consensus_partition(verbose = verbose, .env = .env, max_k = max_k, prefix = prefix,
					top_n = top_n2, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					cores = cores, anno = anno2, anno_col = anno_col, scale_rows = scale_rows, ...)
			
			# if(verbose) qqcat("@{prefix}* .......................................................\n")
			# nm = qq("@{combination_method[[i]][1]}:@{combination_method[[i]][2]}-column")
			# part_list[[nm]] = consensus_partition(verbose = verbose, .env = .env, max_k = max_k, prefix = prefix,
			# 		top_n = top_n2, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
			# 		cores = cores, anno = anno2, anno_col = anno_col, sample_by = "column", ...)
		}
	} else {

		# in consensus_partition_by_down_sampling(), .env$all_top_value_list is always reset to NULL
		# because the columns are randomly sampled and top_value changes. However, here we cache
		# the recent top_value for a top_value_method for downstream process
		.env$all_top_value_list = NULL
		for(i in seq_along(combination_method)) {
			if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")
			.env$column_index = column_index #note .env$column_index is only for passing to `consensus_partition()` function
			nm = qq("@{combination_method[[i]][1]}:@{combination_method[[i]][2]}-row")

			part_list[[nm]] = consensus_partition_by_down_sampling(subset = subset, verbose = verbose, .env = .env, max_k = max_k, prefix = prefix,
					top_n = top_n2, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
					cores = cores, .predict = FALSE, anno = anno2, anno_col = anno_col, scale_rows = scale_rows, ...)

			# if(verbose) qqcat("@{prefix}* ........................................................\n")
			# nm = qq("@{combination_method[[i]][1]}:@{combination_method[[i]][2]}-column")
			# part_list[[nm]] = consensus_partition_by_down_sampling(subset = subset, verbose = verbose, .env = .env, max_k = max_k, prefix = prefix,
			# 		top_n = top_n2, top_value_method = combination_method[[i]][1], partition_method = combination_method[[i]][2], 
			# 		cores = cores, .predict = FALSE, anno = anno2, anno_col = anno_col, sample_by = "column", ...)
		}
	}

	if(verbose) qqcat("@{prefix}* -------------------------------------------------------\n")

	# find the best partitioning result
	stat_tb = lapply(names(part_list), function(nm) {
		part = part_list[[nm]]
		stat_df = get_stats(part, all_stats = TRUE)
		stat_df = as.data.frame(stat_df)
		stat_df$method = nm
		stat_df
	})
	stat_tb = do.call("rbind", stat_tb)

	stat_tb2 = stat_tb
	stat_tb = stat_tb[stat_tb$`mean_silhouette` >= 0.95, , drop = FALSE]
	if(nrow(stat_tb) == 1) {
		ind = do.call(order, -stat_tb[setdiff(colnames(stat_tb), "method")])[1]
		part = part_list[[ stat_tb[ind, "method"] ]]
		best_k = stat_tb[ind, "k"]

		if(verbose) qqcat("@{prefix}* select @{part@top_value_method}:@{part@partition_method} (@{best_k} groups) because this is the only stable partitioning result.\n")

	} else if(nrow(stat_tb) > 1) {
		stat_tb$n_signatures = -Inf
		for(i in seq_len(nrow(stat_tb))) {
			sig_tb = get_signatures(part_list[[stat_tb[i, "method"]]], k = stat_tb[i, "k"], plot = FALSE, verbose = FALSE, simplify = TRUE,
				group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = global_row_mean, .scale_sd = global_row_sd)
			stat_tb$n_signatures[i] = nrow(sig_tb)
		}

		stat_tb = do.call(rbind, tapply(1:nrow(stat_tb), stat_tb$method, function(ind) {
			if(length(ind) == 1) {
				return(stat_tb[ind, , drop = FALSE])
			} else {
				tb = stat_tb[ind, ]
				tb = tb[order(tb$k), ]
				j = 1
				for(i in seq_len(nrow(tb))[-1]) {
					if(tb$n_signatures[i]/tb$n_signatures[j] > 1.1^(i-j)) {
						j = i
					} else {
						tb$n_signatures[i] = -1
					}
				}
				tb = tb[tb$n_signatures > 0, , drop = FALSE]
			}
			tb
		}))

		ind = do.call(order, -stat_tb[, c("n_signatures", setdiff(colnames(stat_tb), c("n_signatures", "k", "method")))])[1]
		part = part_list[[ stat_tb[ind, "method"] ]]
		best_k = stat_tb[ind, "k"]

		if(verbose) qqcat("@{prefix}* select @{part@top_value_method}:@{part@partition_method} (@{best_k} groups) because it has the largest number of signatures among all @{nrow(stat_tb)} stable partitioning results.\n")
	} else {
		ind = do.call(order, -stat_tb2[c("mean_silhouette", setdiff(colnames(stat_tb2), "method"))])[1]
		part = part_list[[ stat_tb2[ind, "method"] ]]
		best_k = stat_tb2[ind, "k"]

		if(verbose) qqcat("@{prefix}* select @{part@top_value_method}:@{part@partition_method} (@{best_k} groups) as the best partitioning result.\n")
	}

	dist_method = list(...)$dist_method
	if(is.null(dist_method)) dist_method = "euclidean"
	if(length(column_index) > subset) {
		part = convert_to_DownSamplingConsensusPartition(part, column_index, dist_method, verbose, prefix, cores)
	}

	attr(part, "node_id") = node_id
	attr(part, "best_k") = best_k
	lt = list(obj = part)

	if(is.na(best_k)) {
    	attr(lt$obj, "stop_reason") = STOP_REASON["z"]
    	if(verbose) qqcat("@{prefix}* Jaccard index is too high, no meaningful subgroups, stop.\n")
    	return(lt)
    }

    if(best_k == 0) {
    	attr(lt$obj, "stop_reason") = STOP_REASON["z"]
    	if(verbose) qqcat("@{prefix}* Jaccard index is too high, no meaningful subgroups, stop.\n")
    	return(lt)
    }

    mean_silhouette = get_stats(part, k = best_k)[, "mean_silhouette"]
    if(mean_silhouette < mean_silhouette_cutoff) {
    	if(verbose) qqcat("@{prefix}* mean_silhouette score is too small (@{sprintf('%.2f', mean_silhouette)}), stop.\n")
    	attr(lt$obj, "stop_reason") = STOP_REASON["a"]
    	return(lt)
    }

    cl = get_classes(part, k = best_k)

    .env$all_top_value_list = NULL

	# check the numbers of signatures
	if(verbose) qqcat("@{prefix}* checking number of signatures in the best classification.\n")
	if(length(column_index) <= subset) {
		sig_df = get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE, simplify = TRUE,
			group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = global_row_mean, .scale_sd = global_row_sd)
	} else {
		sig_df = get_signatures(part, k = best_k, plot = FALSE, verbose = FALSE, simplify = TRUE,
			group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = global_row_mean, .scale_sd = global_row_sd)
	}
	if(is.null(part@.env$signature_hash)) {
		part@.env$signature_hash = list()
	}
	part@.env$signature_hash[[node_id]] = attr(sig_df, "hash")
	
	n_sig = nrow(sig_df)
	p_sig = n_sig/nrow(part)

	if(is.null(min_n_signatures)) {
		if(node_id == "0") {
			min_n_signatures = floor(n_sig*0.05)
			.env$min_n_signatures = min_n_signatures
		}
	}
	if(n_sig <= min_n_signatures && node_id != "0") {
		if(verbose) qqcat("@{prefix}* too few signatures (n = @{n_sig}, p = @{sprintf('%.3f', p_sig)}) under the best classification, stop.\n")
		attr(lt$obj, "stop_reason") = STOP_REASON["c"]
		return(lt)
	}

	if(verbose) qqcat("@{prefix}* best k = @{best_k}, partition into @{best_k} subgroups.\n")

    lt2 = lapply(sort(unique(cl$class)), function(i_class) {
    	set = cl$class == i_class
    	sub_node = paste0(node_id, i_class)
    	return(.hierarchical_partition(.env, column_index = column_index[set], node_id = sub_node, subset = subset, anno = anno, anno_col = anno_col,
    			min_samples = min_samples, max_k = min(max_k, length(set)-1), cores = cores, verbose = verbose, mean_silhouette_cutoff = mean_silhouette_cutoff,
    			top_n = top_n, min_n_signatures = min_n_signatures, group_diff = group_diff, fdr_cutoff = fdr_cutoff, scale_rows = scale_rows, 
    			filter_fun = filter_fun, ...))
    })

    for(i in seq_along(lt2)) {
    	if(!is.null(lt2[[i]])) {
    		lt[[qq("child@{i}")]] = lt2[[i]]
    	}
    }

    return(lt)
}

STOP_REASON = c(
	"z" = "Jaccard indices for all k were too high.",
	"a" = "Mean silhouette score was too small",
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

	if(length(object@param$combination_method) == 1) {
		qqcat("A 'HierarchicalPartition' object with '@{object@list[[1]]@top_value_method}:@{object@list[[1]]@partition_method}' method.\n")
	} else {
		qqcat("A 'HierarchicalPartition' object with @{length(object@param$combination_method)} combinations of top-value methods and partitioning methods.\n")
	}
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{ncol(object@.env$data)} columns.\n")
	qqcat("  Performed in total @{object@list[[1]]@n_partition*length(object@list)*length(object@param$combination_method)} partitions.\n")
	if(has_hierarchy(object)) {
		qqcat("  There are @{length(all_leaves(object))} groups under the following parameters:\n")
		qqcat("    - min_samples: @{object@param$min_samples}\n")
		qqcat("    - mean_silhouette_cutoff: @{object@param$mean_silhouette_cutoff}\n")
		qqcat("    - min_n_signatures: @{object@param$min_n_signatures} (signatures are selected based on:)\n")
		qqcat("      - fdr_cutoff: @{object@param$fdr_cutoff}\n")
		if(object@list[[1]]@scale_rows) {
			qqcat("      - group_diff (scaled values): @{object@param$group_diff}\n")
		} else {
			qqcat("      - group_diff: @{object@param$group_diff}\n")
		}
	}
	cat("\n")

	if(has_hierarchy(object)) {
		cat("Hierarchy of the partition:\n")
		hierarchy = object@hierarchy
		nodes = hierarchy[, 2]
		nc = nchar(nodes)
		names(nc) = nodes
		n = length(nc)
		each_node_max_k = tapply(names(nc), gsub("\\d$", "0", names(nc)), function(x) {
			max(as.numeric(substr(x, nchar(x), nchar(x))))
		})

		parent = structure(hierarchy[, 1], names = hierarchy[, 2])
		all_leaves = all_leaves(object)
		n_columns = object@node_level$n_columns
		n_signatures = object@node_level$n_signatures

		lines = character(n)

		si = NULL
		for(i in seq_len(n)) {
			k_end = each_node_max_k[as.character(gsub("\\d$", "0", names(nc)[i]))]
			lines[i] = paste0("  ", strrep("    ", nc[i] - 2), ifelse(grepl(qq("@{k_end}$"), nodes[i]), "`-", "|-") ,"- ", nodes[i], qq(", @{n_columns[nodes[i]]} cols"))
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
				k_end = each_node_max_k[as.character(gsub("\\d$", "0", p))]
				if((!grepl(qq("@{k_end}$"), p)) && (p != "0")) {
					substr(lines[i], (nc[p] - 2)*4+3, (nc[p] - 2)*4+3) = "|"
				}
			}
		}
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

# == title
# Information on the nodes
#
# == param
# -object A `HierarchicalPartition-class` object.
#
# == details
# It is the same as `node_info,HierarchicalPartition-method`.
setMethod(f = "node_level",
	signature = "HierarchicalPartition",
	definition = function(object) {
	node_info(object)
})

# == title
# Find the knee/elbow of a list of sorted points
#
# == param
# -x A numeric vector.
# -plot Whether to make the plot.
#
# == value
# A vector of two numeric values. One for the left knee and the second for the right knee.
#
# == example
# x = rnorm(1000)
# knee_finder2(x, plot = TRUE)
knee_finder2 = function(x, plot = FALSE) {
	
	y = sort(x)
	x = seq_along(y)

	n = length(x)

	a = (y[n] - y[1])/(x[n] - x[1])
	b = y[1] - a*x[1]
	d = a*x + b - y
	x1 = x[which.min(d)]
	x2 = x[which.max(d)]

	if(all(d >= 0)) x1 = NA
	if(all(d <= 0)) x2 = NA

	if(plot) {
		op = par(no.readonly = TRUE)
		par(mfrow = c(1, 2))
		plot(x, y, xlab = "index", ylab = "value")
		abline(a = b, b = a)
		abline(v = x1, col = "green")
		abline(v = x2, col = "red")

		plot(d, xlab = "inidex", ylab = "distance to diagonal line")
		abline(v = x1, col = "green")
		abline(v = x2, col = "red")
		par(op)
	}

	return(c(x1, x2))
}


find_top_n = function(x, plot = FALSE) {
	x = x[is.finite(x)]
	v = knee_finder2(x, plot = plot)
	length(x) - v[2] + 1
}

# == title 
# Merge node
#
# == param
# -object A `HierarchicalPartition-class` object.
# -node_id A vector of node IDs where each node is merged as a leaf node.
#
# == value
# A `HierarchicalPartition-class` object.
setMethod(f = "merge_node",
	signature = "HierarchicalPartition",
	definition = function(object, node_id) {

	if(length(node_id) > 1) {
		for(i in seq_along(node_id)) {
			object = merge_node(object, node_id[i])
		}
		return(object)
	}

	if(is_leaf_node(object, node_id)) {
		stop_wrap(qq("'@{node_id}' is a leaf node. Maybe you want to merge on node '@{gsub('.$', '', node_id)}'?"))
	}

	all_nodes = all_nodes(object)
	nodes_to_remove = all_nodes[ grepl(qq("^@{node_id}\\d+"), all_nodes) ]

	object@hierarchy = object@hierarchy[ !(object@hierarchy[, 1] %in% c(node_id, nodes_to_remove) | object@hierarchy[, 2] %in% nodes_to_remove), , drop = FALSE]

	object@node_level = lapply(object@node_level, function(x) {
		x[setdiff(names(x), nodes_to_remove)]
	})

	object@subgroup[ object@subgroup %in% nodes_to_remove ] = node_id
	object@subgroup_col = object@subgroup_col[setdiff(names(object@subgroup_col), nodes_to_remove)]

	if(has_hierarchy(object)) {
		object@subgroup_dend = subgroup_dend(object)
	}

	object@list = object@list[setdiff(names(object@list), nodes_to_remove)]

	object
})

# == title 
# Split node
#
# == param 
# -object A `HierarchicalPartition-class` object.
# -node_id A single ID of a node that is going to be split.
# -subset The same as in `hierarchical_partition`.
# -min_samples The same as in `hierarchical_partition`.
# -max_k max_k The same as in `hierarchical_partition`.
# -cores Number of cores.
# -verbose Whether to print messages.
# -top_n The same as in `hierarchical_partition`.
# -min_n_signatures The same as in `hierarchical_partition`.
# -group_diff The same as in `hierarchical_partition`.
# -fdr_cutoff The same as in `hierarchical_partition`.
#
# == details
# It applies hierarchical consensus partitioning on the specified node.
#
# == value
# A `HierarchicalPartition-class` object.
setMethod(f = "split_node",
	signature = "HierarchicalPartition",
	definition = function(object, node_id, 
	subset = object@param$subset,
	min_samples = object@param$min_samples, max_k = object@param$max_k, cores = object@param$cores, 
	verbose = TRUE, 
	top_n = object@param$top_n, min_n_signatures = object@param$min_n_signatures, 
	group_diff = object@param$group_diff, fdr_cutoff = object@param$fdr_cutoff) {

	if(length(node_id) > 1) {
		for(i in seq_along(node_id)) {
			object = split_node(object, node_id[i], subset = subset, min_samples = min_samples, max_k = max_k, cores = cores,
				verbose = verbose, top_n = top_n, min_n_signatures  = min_n_signatures, group_diff = group_diff, fdr_cutoff = fdr_cutoff)
		}
		return(object)
	}

	if(!is_leaf_node(object, node_id)) {
		if(verbose) qqcat("Remove node @{node_id} because it is not a leave node.\n")
		merge_node(object, node_id)
	}

	pnode = object@list[[node_id]]

	lt = .hierarchical_partition(object@.env, column_index = pnode@column_index, node_id = node_id, 
		subset = subset, anno = object@list[[1]]@anno, anno_col = object@list[[1]]@anno_col,
		min_samples = min_samples, max_k = max_k, cores = cores, verbose = verbose, 
		top_n = top_n, min_n_signatures = min_n_signatures, 
		group_diff = group_diff, fdr_cutoff = fdr_cutoff)

	# reformat lt
	.e = new.env(parent = emptyenv())
	.e$hierarchy = matrix(nrow = 0, ncol = 2)
	.e$lt = list()
	reformat_lt = function(lt, .e) {
		nm = names(lt)
		parent_id = attr(lt$obj, "node_id")
		.e$lt[[parent_id]] = lt$obj

		for(nm in names(lt)) {
			if(grepl("^child\\d+$", nm)) {
				child_id = attr(lt[[nm]]$obj, "node_id")
				.e$hierarchy = rbind(.e$hierarchy, c(parent_id, child_id))
				reformat_lt(lt[[nm]], .e)
			}
		}
	}
	reformat_lt(lt, .e)

	hierarchy = .e$hierarchy
	list = .e$lt

	if(nrow(hierarchy) == 0) {
		if(verbose) qqcat("- Node cannot be split because no hierarchy was generated.\n")
		return(object)
	}

	node_level = list()
	node_level$best_k = sapply(.e$lt, function(x) {
		if(inherits(x, "ConsensusPartition")) {
			attr(x, "best_k")
		} else {
			NA
		}
	})
	if(nrow(hierarchy) > 0) {
		tb = table(hierarchy)
		leaves = names(tb[tb <= 1])
	} else {
		leaves = "0"
	}

	subgroup = object@subgroup
	for(le in leaves) {
		if(inherits(.e$lt[[le]], "DownSamplingConsensusPartition")) {
			subgroup[ .e$lt[[le]]@full_column_index ] = le
		} else if(inherits(.e$lt[[le]], "ConsensusPartition")) {
			subgroup[ .e$lt[[le]]@column_index ] = le
		} else {
			subgroup[ attr(.e$lt[[le]], "column_index") ] = le
		}
	}
	object@subgroup = subgroup

	n = length(list)
	n_columns = numeric(n); names(n_columns) = names(list)
	n_signatures = rep(NA_real_, n); names(n_signatures) = names(list)
	nodes = names(list)
	for(i in seq_len(n)) {
		if(inherits(list[[i]], "DownSamplingConsensusPartition")) {
			n_columns[i] = length(list[[i]]@full_column_index)
		} else if(inherits(list[[i]], "ConsensusPartition")) {
			n_columns[i] = length(list[[i]]@column_index)
		} else {
			n_columns[i] = length(attr(list[[i]], "column_index"))
		}
		if(nodes[i] %in% leaves) {
			if(attr(list[[i]], "stop_reason") == STOP_REASON["c"]) {
				sig_tb = get_signatures(list[[i]], k = attr(list[[i]], "best_k"), verbose = FALSE, plot = FALSE, simplify = TRUE,
					group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = object@.env$global_row_mean, .scale_sd = object@.env$global_row_sd)
				n_signatures[i] = nrow(sig_tb)
			} else {
				n_signatures[i] = NA_real_
			}
		} else {
			sig_tb = get_signatures(list[[i]], k = attr(list[[i]], "best_k"), verbose = FALSE, plot = FALSE, simplify = TRUE,
				group_diff = group_diff, fdr_cutoff = fdr_cutoff, .scale_mean = object@.env$global_row_mean, .scale_sd = object@.env$global_row_sd)
			n_signatures[i] = nrow(sig_tb)
		}
	}

	node_level$n_columns = n_columns
	node_level$n_signatures = n_signatures
	node_level$p_signatures = n_signatures/nrow(object)

	subgroup_col = object@subgroup_col
	le = setdiff(as.vector(hierarchy), all_nodes(object))
	col_pal = Polychrome::kelly.colors(22)
	col_pal = col_pal[!(names(col_pal) %in% c("white", "black"))]
	col_pal = setdiff(col_pal, subgroup_col)
	if(length(le) <= length(col_pal)) {
		subgroup_col2 = structure(col_pal[seq_along(le)], names = le)
	} else {
		subgroup_col2 = structure(c(col_pal, rand_color(length(le) - length(col_pal))), names = le)
	}
	subgroup_col = c(subgroup_col, subgroup_col2)
	object@subgroup_col = subgroup_col

	for(nm in names(node_level)) {
		v = object@node_level[[nm]]
		cn = intersect(names(v), names(node_level[[nm]]))
		object@node_level[[nm]][cn] = node_level[[nm]][cn]
		cn2 = setdiff(names(node_level[[nm]]), names(v))
		object@node_level[[nm]] = c(object@node_level[[nm]], node_level[[nm]][cn2])
	}
	
	object@hierarchy = rbind(object@hierarchy, hierarchy)
	object@hierarchy = reorder_hierarchy(object@hierarchy)
	cn = intersect(names(object@list), names(list))
	object@list[cn] = list[cn]
	cn2 = setdiff(names(list), names(object@list))
	object@list = c(object@list, list[cn2])

	return(object)
})

reorder_hierarchy = function(hierarchy) {

	e = new.env()
	e$od = NULL
	look_child = function(hierarchy, node) {
		ind = which(hierarchy[, 1] == node)
		if(length(ind) == 0) {
			return(NULL)
		} else {
			for(i in ind) {
				e$od = c(e$od, i)
				look_child(hierarchy, hierarchy[i, 2])
			}
		}
	}
	look_child(hierarchy, "0")
	hierarchy[e$od, , drop = FALSE]

}
