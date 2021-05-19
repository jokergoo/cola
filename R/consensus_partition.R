
# == title
# Consensus partition
#
# == param
# -data A numeric matrix where subgroups are found by columns.
# -top_value_method A single top-value method. Available methods are in `all_top_value_methods`.
#                   Use `register_top_value_methods` to add a new top-value method.
# -top_n Number of rows with top values. The value can be a vector with length > 1. When n > 5000, 
#        the function only randomly sample 5000 rows from top n rows. If ``top_n`` is a vector, paritition
#        will be applied to every values in ``top_n`` and consensus partition is summarized from all partitions.
# -partition_method A single partitioning method. Available methods are in `all_partition_methods`.
#                   Use `register_partition_methods` to add a new partition method.
# -max_k Maximal number of subgroups to try. The function will try for ``2:max_k`` subgroups
# -k Alternatively, you can specify a vector k.
# -sample_by Should randomly sample the matrix by rows or by columns?
# -p_sampling Proportion of the submatrix which contains the top n rows to sample.
# -partition_repeat Number of repeats for the random sampling.
# -partition_param Parameters for the partition method which are passed to ``...`` in a registered partitioning method. See `register_partition_methods` for detail.
# -anno A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation
#       to predicted subgroups will be tested.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -scale_rows Whether to scale rows. If it is ``TRUE``, scaling method defined in `register_partition_methods` is used.
# -verbose Whether print messages.
# -mc.cores Multiple cores to use. This argument will be removed in future versions.
# -cores Number of cores, or a ``cluster`` object returned by `parallel::makeCluster`.
# -prefix Internally used.
# -.env An environment, internally used.
# -help Whether to print help messages.
#
# == details
# The function performs analysis in following steps:
#
# - calculate scores for rows by top-value method,
# - for each top_n value, take top n rows,
# - randomly sample ``p_sampling`` rows from the top_n-row matrix and perform partitioning for ``partition_repeats`` times,
# - collect partitions from all individual partitions and summarize a consensus partition.
#
# == return
# A `ConsensusPartition-class` object. Simply type object in the interactive R session
# to see which functions can be applied on it.
#
# == seealso
# `run_all_consensus_partition_methods` runs consensus partitioning with multiple top-value methods
# and multiple partitioning methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# m = cbind(rbind(matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
#           rbind(matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
#           rbind(matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
#                 matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20))
#          ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
# res = consensus_partition(m, partition_repeat = 10, top_n = c(10, 20, 50))
# res
consensus_partition = function(data,
	top_value_method = "ATC",
	top_n = NULL,
	partition_method = "skmeans",
	max_k = 6, 
	k = NULL,
	sample_by = "row",
	p_sampling = 0.8,
	partition_repeat = 50,
	partition_param = list(),
	anno = NULL,
	anno_col = NULL,
	scale_rows = NULL,
	verbose = TRUE,
	mc.cores = 1, cores = mc.cores,
	prefix = "",
	.env = NULL,
	help = cola_opt$help) {

	if(help && !missing(data)) {
		if(identical(subset, Inf) && ncol(data) > 500) {
			qqcat_wrap("You have quite a lot of columns in the matrix. For reducing the runtime, you can use the function `consensus_partition_by_down_sampling()` to apply to a subset of column. The classification of unselected columns are inferred from the classes of the selected columns. Set the argument 'help = FALSE' to turn off this message.\n")
		}
	}

	if(missing(data)) {
		data = .env$data
	} else {
		data = as.matrix(data)
	}

	if(max_k < 2) {
		stop_wrap("max_k should be no less than 2.")
	}

	if(max_k >= 10) {
		if(help) {
			qqcat_wrap("It is not recommended to set `max_k` larger than 10. Users are suggested to use `hierarchical_partition()` function to obtain more subgroups. Set the argument `help` to FALSE to turn off this message.\n")
		}
	}

	t = system.time(res <- .consensus_partition(
		data = data,
		top_value_method = top_value_method,
		top_n = top_n,
		partition_method = partition_method,
		k = if(is.null(k)) 2:max_k else k, 
		p_sampling = p_sampling,
		sample_by = sample_by,
		partition_repeat = partition_repeat,
		partition_param = partition_param,
		anno = anno,
		anno_col = anno_col,
		scale_rows = scale_rows,
		verbose = verbose,
		cores = cores,
		prefix = prefix,
		.env = .env))
	res@hash = digest(res)
	res@running_time = t[["elapsed"]]

	if(verbose) {
		tc = Sys.time()
		tf = format(tc + structure(t[["elapsed"]], units = "secs", class = "difftime") - tc)
		qqcat("@{prefix}* @{top_value_method}:@{partition_method} used @{tf}.\n")
	}

	# test if k = 2 is stable and correlates to column mean
	# if(help) {
	# 	if(is_stable_k(res, k = 2)) {
	# 		cl = get_classes(res, k = 2)[, "class"]
	# 		le = unique(cl)
	# 		cm_diff = rowMeans(get_matrix(res[, cl == le[1], drop = FALSE])) - colMeans(get_matrix(res[, cl == le[2]]))
	# 		tb = table(sign(cm_diff))
	# 		if(max(cm) - min(cm) > 1e-10) {
	# 			p <- t.test(cm[cl == le[1]], cm[cl == le[2]])$p.value
	# 			if(p < 0.01) {
	# 				message_wrap(qq("Classification with 2 groups is stable and it significantly (p-value  = @{p}) correlates to the column mean of the matrix. You might have system-level batch in the matrix. Consider to transform your matrix, such as by quantile normalization. Set `help = FALSE` to turn it off."))
	# 			}
	# 		}
	# 	}
	# }
	return(res)
}

.consensus_partition = function(data,
	top_value_method = "MAD",
	top_n = NULL,
	partition_method = "kmeans",
	k = 2:6, 
	p_sampling = 0.8,
	sample_by = c("row", "column"), 
	partition_repeat = 50,
	partition_param = list(),
	anno = NULL,
	anno_col = NULL,
	scale_rows = NULL,
	verbose = TRUE,
	cores = 1,
	prefix = "",
	.env = NULL) {

	# .env is used to store shared matrix because the matrix will be used multiple times in
	# run_all_consensus_partition_methods() and hierarchical_partition() and we don't want to
	# copy it multiple times
	# `column_index` is specifically for hierarchical_partition() where in each level only a subset
	# of columns in the original matrix is used
	# .env$column_index is only used to pass secret column_index parameter, 
	# the real column_index is stored as a slot in this object.
	if(is.null(.env)) {
		if(is.data.frame(data)) data = as.matrix(data)
		# if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		.env = new.env(parent = emptyenv())
		.env$data = data
		.env$column_index = seq_len(ncol(data))
	} else if(is.null(.env$data)) {
		if(is.data.frame(data)) data = as.matrix(data)
		# if(is.null(rownames(data))) rownames(data) = seq_len(nrow(data))

		.env$data = data
		.env$column_index = seq_len(ncol(data))
	} else if(is.null(.env$column_index)) {
		data = .env$data
		.env$column_index = seq_len(ncol(data))
	} else {
		data = .env$data
	}

	if(is.null(.env$row_index)) {
		.env$row_index = seq_len(nrow(data))
	}

	data = data[.env$row_index, .env$column_index, drop = FALSE]

	# process the annotation
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

		if(nrow(anno) != length(.env$column_index)) {
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

	if(verbose) qqcat("@{prefix}* run @{top_value_method}:@{partition_method} on a @{nrow(data)}x@{ncol(data)} matrix.\n")

	k = sort(k)
	l = k <= ncol(data)
	if(sum(l) != length(k)) {
		qqcat("@{prefix}* Following k (@{paste(k[!l], collapse=', ')}) are removed.\n")
	}
	k = k[l]
	if(length(k) == 0) {
		stop_wrap("There is no valid k.\n")
	}

	partition_fun = get_partition_method(partition_method, partition_param)

	get_top_value_fun = get_top_value_method(top_value_method)
	if(identical(get_top_value_method(top_value_method), ATC) && cores > 1) {
		config_ATC(cores = cores)
		if(verbose) qqcat("@{prefix}* set @{cores} cores for ATC()\n")
	}

	# also since one top value metric will be used for different partition methods,
	# we cache the top values for repetitive use
	# they are only used when column_index and row_index are not changed
	if(is.null(.env$all_top_value_list)) {
		if(verbose) qqcat("@{prefix}* calculating @{top_value_method} values.\n")
		all_top_value = get_top_value_fun(data)
		all_top_value[is.na(all_top_value)] = -Inf
		.env$all_top_value_list = list()
		.env$all_top_value_list[[top_value_method]] = all_top_value
	} else if(is.null(.env$all_top_value_list[[top_value_method]])) {
		if(verbose) qqcat("@{prefix}* calculating @{top_value_method} values.\n")
		all_top_value = get_top_value_fun(data)
		all_top_value[is.na(all_top_value)] = -Inf
		.env$all_top_value_list[[top_value_method]] = all_top_value
	} else {
		if(verbose) qqcat("@{prefix}* @{top_value_method} values have already been calculated. Get from cache.\n")
		all_top_value = .env$all_top_value_list[[top_value_method]]
	}

	if(is.null(top_n)) {
		top_n = find_top_n(all_top_value)
		top_n = min(top_n, round(length(all_top_value)*0.1))
	}
	top_n = round(top_n)
	l = top_n <= nrow(data)
	if(sum(l) != length(top_n)) {
		qqcat("@{prefix}* Following top_n (@{paste(top_n[!l], collapse = ', ')}) are removed.\n")
	}
	top_n = top_n[l]
	if(length(top_n) == 0) {
		stop_wrap("There is no valid top_n.\n")
	}

	if(top_n < 20) browser()

	if(is.null(scale_rows)) {
		scale_rows = TRUE
	}
	if(scale_rows) {
		scale_method = attr(partition_fun, "scale_method")
		if("z-score" %in% scale_method) {
			if(verbose) qqcat("@{prefix}* rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd\n")
			# data = t(scale(t(data)))
			row_mean = rowMeans(data)
			row_sd = rowSds(data)
			data = (data - row_mean)/row_sd
		} else if("min-max" %in% scale_method) {
			if(verbose) qqcat("@{prefix}* rows are scaled before sent to partition, method: 'min-max' (x - min)/(max - min)\n")
			row_min = rowMins(data)
			row_max = rowMaxs(data)
			row_range = row_max - row_min
			data = (data - row_min)/row_range
		} else {
			scale_rows = FALSE
		}
		attr(scale_rows, "scale_method") = scale_method
	}

	# in case NA is produced in scaling
	l = apply(data, 1, function(x) any(is.na(x)))
	if(any(l)) {
		if(verbose) qqcat("@{prefix}* remove @{sum(l)} rows with zero standard deviations.\n")
		data = data[!l, , drop = FALSE]
		all_top_value = all_top_value[!l]
		.env$row_index = .env$row_index[!l]
		l = top_n <= nrow(data)
		top_n = top_n[l]
		

		if(sum(l) != length(top_n)) {
			qqcat("@{prefix}* Following top_n (@{paste(top_n[!l], collapse = ', ')}) are removed.\n")
		}
		top_n = top_n[l]
		if(length(top_n) == 0) {
			stop_wrap("There is no valid top_n.\n")
		}
	}

	# now we do repetitive clustering
	param = data.frame(top_n = numeric(0), k = numeric(0), n_row = numeric(0))
	partition_list = list()
	n_cores = get_nc(cores)

	sample_by = match.arg(sample_by)[1]
	for(i in seq_len(length(top_n))) {
		ind = order(all_top_value, decreasing = TRUE)[ 1:top_n[i] ]

		if(length(ind) > 5000) {
			ind = sample(ind, 5000)
			if(verbose) qqcat("@{prefix}* get top @{top_n[i]} (randomly sampled 5000) rows by @{top_value_method} method\n")
		} else {
			if(verbose) qqcat("@{prefix}* get top @{top_n[i]} rows by @{top_value_method} method\n")
		}

		if(verbose && n_cores > 1) qqcat("@{prefix}  - @{partition_method} repeated for @{partition_repeat} times by @{sample_by}-sampling (p = @{p_sampling}) from top @{top_n[i]} rows (@{n_cores} cores).\n")

		registerDoParallel(cores)
		lt <- foreach(j = seq_len(partition_repeat)) %dopar% {

			param = data.frame(top_n = numeric(0), k = numeric(0), n_row = numeric(0))
			partition_list = list()
			
			if(sample_by == "row") {
				ind_sub = sample(ind, round(p_sampling*length(ind)))
				mat = data[ind_sub, , drop = FALSE]
				column_ind_sub = seq_len(ncol(data))
			} else if(sample_by == "column") {
				column_ind_sub = sample(ncol(data), round(p_sampling*ncol(data)))
				mat = data[ind, , drop = FALSE]
			}
			for(y in k) {
				if(interactive() && verbose && n_cores == 1) {
					msg = qq("@{prefix}  [k = @{y}] @{partition_method} repeated for @{j}@{ifelse(j %% 10 == 1, 'st', ifelse(j %% 10 == 2, 'nd', ifelse(j %% 10 == 3, 'rd', 'th')))} @{sample_by}-sampling (p = @{p_sampling}) from top @{top_n[i]} rows.")
					cat(strrep("\r", nchar(msg)))
					cat(msg)
				}
				partition_list = c(partition_list, list(list(partition_fun(mat, y, column_ind_sub))))
				param = rbind(param, data.frame(top_n = top_n[i], k = y, n_row = nrow(mat), n_col = ncol(mat), stringsAsFactors = FALSE))
			}

			return(list(param = param, partition_list = partition_list))
		}
		stopImplicitCluster()

		for(i in seq_along(lt)) {
			param = rbind(param, lt[[i]]$param)
			partition_list = c(partition_list, lt[[i]]$partition_list)
		}
		if(interactive() && verbose && n_cores == 1) cat("\n")
	}

	construct_consensus_object = function(param, partition_list, k, prefix = "  - ", verbose = TRUE) {

		partition_list = do.call("c", partition_list)

		partition_list = cl_ensemble(list = partition_list)
		if(verbose) qqcat("@{prefix}merge @{length(partition_list)} (@{partition_repeat}x@{length(top_n)}) partitions into a single ensemble object.\n")
		if(sample_by == "row") {
			partition_consensus = cl_consensus(partition_list)
		} else {
			partition_consensus = cl_consensus2(partition_list, k)
		}

		# note: number of class_ids may be less than k
		class_ids = as.vector(cl_class_ids(partition_consensus))
		# adjust the class labels according to the tightness of each subgroup
		if(verbose) qqcat("@{prefix}adjust the class labels according to the mean intra-group distance.\n")
		ri = order(all_top_value, decreasing = TRUE)[1:max(param[, "top_n"])]
		if(length(ri) > 5000) ri = sample(ri, 5000)
		mean_dist = tapply(seq_len(ncol(data)), class_ids, function(ind) {
			n = length(ind)
			if(n == 1) {
				return(Inf)
			}
			sum(dist(t(data[ri, ind, drop = FALSE]))^2)/(n*(n-1)/2)
		})
		if(length(mean_dist) < k) {
			mean_dist_foo = structure(rep(Inf, k - length(mean_dist)), names = setdiff(seq_len(k), class_ids))
			mean_dist = c(mean_dist, mean_dist_foo)
		}
		map = structure(names = names(mean_dist)[order(mean_dist)], names(mean_dist))
		unclassfied = setdiff(1:k, map)
		if(length(unclassfied)) {
			map = c(map, structure(unclassfied, names = unclassfied))
		}
		class_ids = as.numeric(map[as.character(class_ids)])

		if(verbose) qqcat("@{prefix}calculate global membership from all partitions.\n")
		membership_mat = cl_membership(partition_consensus)
		class(membership_mat) = "matrix"
		if(ncol(membership_mat) < k) {
			membership_mat = cbind(membership_mat, matrix(0, nrow = nrow(membership_mat), ncol = k - ncol(membership_mat)))
		}
		map2 = structure(names(map), names = map)
		membership_mat = membership_mat[, as.numeric(map2[as.character(1:k)])]

		colnames(membership_mat) = paste0("p", 1:ncol(membership_mat))
		rownames(membership_mat) = colnames(data)
		attr(membership_mat, "n_of_classes") = NULL
		attr(membership_mat, "is_cl_hard_partition") = NULL
		
		if(verbose) qqcat("@{prefix}adjust class labels for every single partition.\n")
		class_ids_by_top_n = tapply(seq_along(partition_list), param$top_n, function(ind) {
			if(sample_by == "row") {
				partition_consensus = cl_consensus(cl_ensemble(list = partition_list[ind]))
			} else {
				partition_consensus = cl_consensus2(cl_ensemble(list = partition_list[ind]), k)
			}
			ci = as.vector(cl_class_ids(partition_consensus))
			map = relabel_class(ci, class_ids, full_set = 1:k)
			as.numeric(map[as.character(ci)])
		})

		# adjust class labels in each membership matrix to fit to the consensus class labels
		membership_each = do.call("cbind", lapply(seq_along(partition_list), function(i) {
			x = partition_list[[i]]
			class_ids = class_ids_by_top_n[[as.character(param$top_n[i])]]
			class = as.vector(cl_class_ids(x))
			map = relabel_class(class, class_ids, full_set = 1:k)
			class = as.numeric(map[as.character(class)])
			as.integer(class)
		}))
		rownames(membership_each) = rownames(membership_mat)

		if(verbose) qqcat("@{prefix}calculate consensus matrix.\n")
		# consensus_mat = matrix(1, nrow = nrow(membership_mat), ncol = nrow(membership_mat))
		# for(i in seq_len(nrow(membership_each)-1)) {
		# 	for(j in (i+1):nrow(membership_each)) {
		# 		consensus_mat[i, j] = sum(membership_each[i, ] == membership_each[j, ] + 0)/ncol(membership_each)
		# 		consensus_mat[j, i] = consensus_mat[i, j]
		# 	}
	 # 	}
		membership_each2 = membership_each
		membership_each2[is.na(membership_each2)] = as.integer(0)
		consensus_mat = get_consensus_matrix(membership_each2)
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

		if(verbose) qqcat("@{prefix}calculate statistics for the consensus partition.\n")
		if(length(class_df$silhouette)*0.05 > 1) {
			l = class_df$silhouette >= quantile(class_df$silhouette, 0.05)
		} else {
			l = rep(TRUE, length(class_df$silhouette))
		}
		stat = list(
			ecdf = ecdf(consensus_mat[lower.tri(consensus_mat)]),
			"1-PAC" = 1 - PAC(consensus_mat[l, l, drop = FALSE]),
			mean_silhouette = mean(class_df$silhouette),
			concordance = concordance(membership_each, class_ids),
			cophcor = cophcor(consensus_mat),
			aPAC = aPAC(consensus_mat),
			FCC = FCC(consensus_mat[l, l, drop = FALSE])
		)

		# loc = cmdscale(dist(t(data)), k = 2)
		# stat$mean_group_dist_2PC = mean_group_dist(loc, class_ids)
		# loc = cmdscale(dist(t(data)), k = 3)
		# stat$mean_group_dist_3PC = mean_group_dist(loc, class_ids)
		
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
		
		top_n_level = unique(param[l, "top_n"])
		if(verbose) qqcat("@{prefix}* wrap results for k = @{y}\n")
		construct_consensus_object(param[l, ], partition_list[l], y, verbose = FALSE)
		
	})
	names(object_list) = as.character(k)

	rm(partition_list)
	gc(verbose = FALSE, reset = TRUE)

	## adjust class labels for different k should also be adjusted
	if(verbose) qqcat("@{prefix}* adjust class labels between different k.\n")
	reference_class = object_list[[1]]$class_df$class

	for(i in seq_along(k)[-1]) {
		class_df = object_list[[i]]$class_df
    	class = class_df[, "class"]

    	map = relabel_class(class, reference_class, full_set = 1:(k[i]))
    	l = which( (duplicated(map) | duplicated(map, fromLast = TRUE)) & map != names(map))
    	unmapped = setdiff(names(map), map)
    	if(any(l)) {
    		map[l] = unmapped
    	}
    	map2 = structure(names(map), names = map)
    	object_list[[i]]$class_df$class = as.numeric(map[as.character(class)])

    	# the class label for the global membership matrix needs to be adjusted
    	object_list[[i]]$membership = object_list[[i]]$membership[, as.numeric(map2[as.character(1:k[i])]) ]
		colnames(object_list[[i]]$membership) = paste0("p", 1:k[i])
			
		# the class label for the individual membership needs to be adjusted
		odim = dim(object_list[[i]]$membership_each)
		object_list[[i]]$membership_each = as.numeric(map[as.character(object_list[[i]]$membership_each)])
		dim(object_list[[i]]$membership_each) = odim

		reference_class = object_list[[i]]$class_df$class
	}

	# an additional metric for determine "best k"
	ak = sapply(object_list, function(obj) {
		f = obj$stat$ecdf
		x = seq(0, 1, length = 1000)
		n = length(x)
		sum((x[2:n] - x[1:(n-1)])*f(x[2:n]))
	})
	names(ak) = NULL
	delta_k = ak
	for(i in seq_along(k)[-1]) {
		delta_k[i] = (ak[i] - ak[i-1])/ak[i-1]
	}
	for(i in seq_along(object_list)) {
		object_list[[i]]$stat$area_increased = delta_k[i]
	}

	n_sample = ncol(data)
	# for(method in c("Rand", "cRand", "NMI", "KP", "FM", "Jaccard", "purity", "PS")) {
	for(method in c("Rand", "Jaccard")) {
		for(i in seq_along(k)) {
			if(i == 1) {
				cl1 = rep(1, n_sample)
			} else {
				cl1 = object_list[[i - 1]]$class_df$class
			}
			cl2 = object_list[[i]]$class_df$class
			object_list[[i]]$stat[[method]] = cl_agreement(as.cl_hard_partition(cl1), as.cl_hard_partition(cl2), method = method)[[1]]
		}
	}
	
	# process the annotations so it can be shared in `run_all_consensus_partition_methods()` and `hierarchical_partitions()`

	res = ConsensusPartition(object_list = object_list, k = k, n_partition = partition_repeat * length(top_n) * length(k),  
		partition_method = partition_method, top_value_method = top_value_method, top_n = top_n, top_value_list = all_top_value,
		anno = anno, anno_col = anno_col, scale_rows = scale_rows, sample_by = sample_by,
		column_index = .env$column_index, row_index = .env$row_index, .env = .env)

	.env$column_index = NULL
	return(res)
}



# == title
# Print the ConsensusPartition object
#
# == param
# -object A `ConsensusPartition-class` object.
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

	# fix older version where there was no sample_by slot
	error = try(object@sample_by, silent = TRUE)
	if(inherits(error, "try-error")) {
		object@sample_by = "row"
	}
	qqcat("A 'ConsensusPartition' object with k = @{paste(object@k, collapse = ', ')}.\n")
	qqcat("  On a matrix with @{length(object@row_index)} rows and @{length(object@column_index)} columns.\n")
	top_n_str = object@top_n
	qqcat("  Top rows (@{paste(top_n_str, collapse = ', ')}) are extracted by '@{object@top_value_method}' method.\n")
	qqcat("  Subgroups are detected by '@{object@partition_method}' method.\n")
	qqcat("  Performed in total @{object@n_partition} partitions by @{object@sample_by} resampling.\n")
	best_k = suggest_best_k(object, help = FALSE)
	if(is.na(best_k)) {
		qqcat("  There is no best k.\n")
	} else {
		qqcat("  Best k for subgroups seems to be @{best_k}.\n")
	}
	qqcat("\n")
	qqcat("Following methods can be applied to this 'ConsensusPartition' object:\n")
	txt = showMethods(classes = "ConsensusPartition", where = topenv(), printTo = FALSE)
	txt = grep("Function", txt, value = TRUE)
	fname = gsub("Function: (.*?) \\(package.*$", "\\1", txt)
	print(fname)
})

# == title
# Plot the empirical cumulative distribution (eCDF) curve of the consensus matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -... Other arguments.
#
# == details
# It plots eCDF curve for each k.
#
# This function is mainly used in `collect_plots` and `select_partition_number` functions.
#
# == value
# No value is returned.
#
# == seealso
# See `stats::ecdf` for a detailed explanation of the empirical cumulative distribution function.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# plot_ecdf(golub_cola["ATC", "skmeans"])
setMethod(f = "plot_ecdf",
	signature = "ConsensusPartition",
	definition = function(object, ...) {

	omar = par("mar")
	par(mar = c(4.1, 4.1, 1, 1))
	on.exit(par(mar = omar))
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "consensus value (x)", ylab = "P(X <= x)")
	for(i in seq_along(object@k)) {
		consensus_mat = get_consensus(object, k = object@k[i])
		f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		x = seq(0, 1, length = 100)
		y = f(x)
		x = c(0, x)
		y = c(0, y)
		lines(x, y, col = i)
	}
	legend("bottomright", pch = 15, legend = paste0("k = ", object@k), col = seq_along(object@k))
})

# == title
# Several plots for determining the optimized number of subgroups
#
# == param
# -object A `ConsensusPartition-class` object.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == details
# There are following plots made:
#
# - eCDF of the consensus matrix under each k, made by `plot_ecdf,ConsensusPartition-method`,
# - `PAC` score,
# - mean sihouette score,
# - the `concordance` for each partition to the consensus partition,
# - area increase of the area under the ECDF of consensus matrix with increasing k,
# - Rand index for current k compared to k - 1,
# - Jaccard coefficient for current k compared to k - 1,
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# select_partition_number(golub_cola["ATC", "skmeans"])
setMethod(f = "select_partition_number",
	signature = "ConsensusPartition",
	definition = function(object, mark_best = TRUE, all_stats = FALSE) {
	op = par(no.readonly = TRUE)

	m = get_stats(object, all_stats = all_stats)
	if(all_stats) STAT_USED = STAT_ALL
	m = m[, colnames(m) %in% c(STAT_USED, "area_increased", "Rand", "Jaccard"), drop = FALSE]
	nm = colnames(m)
	n = ncol(m)
	nr = floor(sqrt(n)+1)

	par(mfrow = c(nr, ceiling((n+1)/nr)), mar = c(4, 4, 1, 1))

	plot_ecdf(object, lwd = 1)

	best_k = suggest_best_k(object, help = FALSE)
	if(is.na(best_k)) best_k = -1

	for(i in seq_len(ncol(m))) {
		
		if(mark_best) {
			l = object@k == best_k
		} else {
			l = rep(FALSE, length(object@k))
		}
		plot(object@k, m[, i], type = "b", xlab = "k", ylab = nm[i])
		if(any(l)) {
			points(object@k, m[, i],
				pch = ifelse(l, 16, 1),
				col = ifelse(l, "red", "black"))
		}
	}

	par(xpd = NA, mar = c(4, 2, 1, 1))
	plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ann = FALSE)
# 	if(max(object@k) <= 6) {
# 		text(x = 0, y = 1,
# "Suggested rules:
# Jaccard index < 0.95;
# If 1-PAC >= 0.90,
#   take the maximum k;
# else take the k with higest votes of
#   1. max 1-PAC,
#   2. max mean silhouette,
#   3. max concordance.
# ", cex = 1.2, adj = c(0, 1))
# 	} else {
# 		text(x = 0, y = 1,
# "Suggested rules:
# Jaccard index < 0.95;
# take the k with higest votes of
#   1. max 1-PAC,
#   2. max mean silhouette,
#   3. max concordance.
# ", cex = 1.2, adj = c(0, 1))
# 	}

	if(mark_best) legend("topright", pch = 16, col = "Red", legend = "best k", cex = 1.5)

	par(op)
})


# == title
# Heatmap of the consensus matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -internal Used internally.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `consensus_partition` or `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_row_names Whether plot row names on the consensus heatmap (which are the column names in the original matrix)
# -row_names_gp Graphics parameters for row names.
# -simplify Internally used.
# -... other arguments
#
# == details
# For row i and column j in the consensus matrix, the value of corresponding x_ij
# is the probability of sample i and sample j being in a same group from all partitions.
#
# There are following heatmaps from left to right:
#
# - probability of the sample to stay in the corresponding group
# - silhouette scores which measure the distance of an item to the second closest subgroups.
# - predicted subgroups
# - consensus matrix.
# - more annotations if provided as ``anno``
#
# One thing that is very important to note is that since we already know the consensus subgroups from consensus
# partition, in the heatmap, only rows or columns within the group is clustered.
#
# == value
# No value is returned.
#
# == seealso
# `membership_heatmap,ConsensusPartition-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# consensus_heatmap(golub_cola["ATC", "skmeans"], k = 3)
setMethod(f = "consensus_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, internal = FALSE,
	anno = object@anno, anno_col = get_anno_col(object), 
	show_row_names = FALSE, row_names_gp = gpar(fontsize = 8),
	simplify = FALSE, ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

	if(inherits(object, "DownSamplingConsensusPartition")) {
		class_df = get_classes(object, k, reduce = TRUE)
	} else {
		class_df = get_classes(object, k)
	}
	class_ids = class_df$class

	consensus_mat = get_consensus(object, k)

	mat_col_od = column_order_by_group(factor(class_ids, levels = sort(unique(class_ids))), consensus_mat)

	membership_mat = get_membership(object, k)

	if(simplify) {
		ht_list = NULL
	} else {
		ht_list = Heatmap(membership_mat, name = "Prob", cluster_columns = FALSE, show_row_names = FALSE,
				width = unit(5, "mm")*k, col = colorRamp2(c(0, 1), c("white", "red")),
				show_column_names = !internal) + 
			Heatmap(class_df$silhouette, name = "Silhouette", width = unit(5, "mm"),
				show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "purple")),
				show_column_names = !internal)
	}
	ht_list = ht_list + Heatmap(class_ids, name = "Class", col = cola_opt$color_set_2,
		show_row_names = FALSE, width = unit(5, "mm"),
		show_column_names = !internal)
	
	ht_list = ht_list +	Heatmap(consensus_mat, name = "Consensus", show_row_names = FALSE, show_row_dend = FALSE,
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
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
			}
		}
		if(is.null(anno_col))
			ht_list = ht_list + rowAnnotation(df = anno, show_annotation_name = !internal)
		else {
			ht_list = ht_list + rowAnnotation(df = anno, col = anno_col, show_annotation_name = !internal)
		}
	}

	if(show_row_names && !is.null(rownames(consensus_mat))) {
		ht_list = ht_list + rowAnnotation(rn = anno_text(rownames(consensus_mat), gp = row_names_gp))
	}
	if(internal) {
		column_title = NULL
	} else {
		column_title = qq("consensus @{object@partition_method} with @{k} groups from @{object@n_partition/length(object@k)} partitions")
	}
	draw(ht_list, main_heatmap = "Consensus", column_title = column_title,
		show_heatmap_legend = !internal, show_annotation_legend = !internal)
})

# == title
# Heatmap of membership in each partition
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -internal Used internally.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `consensus_partition` or `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_column_names Whether show column names in the heatmap (which is the column name in the original matrix).
# -column_names_gp Graphics parameters for column names.
# -... Other arguments
#
# == details
# Each row in the heatmap is the membership in one single partition.
#
# Heatmap is split on rows by ``top_n``.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# membership_heatmap(golub_cola["ATC", "skmeans"], k = 3)
setMethod(f = "membership_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, internal = FALSE, 
	anno = object@anno, anno_col = get_anno_col(object),
	show_column_names = FALSE, column_names_gp = gpar(fontsize = 8), ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

	if(inherits(object, "DownSamplingConsensusPartition")) {
		class_df = get_classes(object, k, reduce = TRUE)
	} else {
		class_df = get_classes(object, k)
	}
	class_ids = class_df$class

	membership_mat = get_membership(object, k)
	col_fun = colorRamp2(c(0, 1), c("white", "red"))

	membership_each = get_membership(object, k, each = TRUE)
	membership_each = t(membership_each)
	mat_col_od = column_order_by_group(factor(class_ids, levels = sort(unique(class_ids))), membership_each)

	col = cola_opt$color_set_1[1:k]

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
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
			}
		}

		if(is.null(anno_col)) {
			bottom_anno = HeatmapAnnotation(df = anno,
				show_annotation_name = !internal, annotation_name_side = "right")
		} else {
			bottom_anno = HeatmapAnnotation(df = anno, col = anno_col,
				show_annotation_name = !internal, annotation_name_side = "right")
		}
	}

	param = get_param(object, k, unique = FALSE)

	top_n_level = unique(param$top_n)
	suppressWarnings(n_row_col <- structure(brewer.pal(length(top_n_level), "Accent"), names = top_n_level))
	
	ht = Heatmap(as.character(param$top_n), name = "top_n", col = n_row_col,
		width = unit(10, "mm"), show_row_names = FALSE, show_heatmap_legend = FALSE,
		show_column_names = FALSE) +
	Heatmap(membership_each, name = "Class", show_row_dend = FALSE, show_column_dend = FALSE, col = cola_opt$color_set_2,
		column_title = ifelse(internal, "", qq("membership heatmap, k = @{k}")), 
		column_order = mat_col_od, cluster_columns = FALSE,
		row_split = factor(param$top_n, levels = sort(unique(param$top_n))),
		cluster_row_slices = FALSE,
		top_annotation = HeatmapAnnotation(Prob = membership_mat,
			Class = class_ids, col = c(list(Class = cola_opt$color_set_2), Prob = col_fun),
			show_annotation_name = !internal),
		bottom_annotation = bottom_anno,
		show_column_names = show_column_names, column_names_gp = column_names_gp,
		row_title = NULL
	)
	if(internal) {
		row_title = NULL
	} else {
		row_title = qq("@{round(object@n_partition/length(object@k)/length(top_n_level))} x @{length(top_n_level)} random samplings")
	}
	draw(ht, main_heatmap = "Class", row_title = row_title,
		show_heatmap_legend = FALSE, show_annotation_legend = !internal)
	param2 = get_param(object, k)
	if(!internal) {
		for(i in seq_along(top_n_level)) {
			decorate_heatmap_body("top_n", slice = i, {
				grid.text(qq("top @{top_n_level[i]} rows"), rot = 90, gp = gpar(fontsize = 10))
			})
		}
	}
})
