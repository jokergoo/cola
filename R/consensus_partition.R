
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
# -partition_method A single partition method. Available methods are in `all_partition_methods`.
#                   Use `register_partition_methods` to add a new partition method.
# -max_k Maximal number of partitions to try. The function will try ``2:max_k`` partitions.
# -sample_by Should randomly sample the matrix by rows or by columns?
# -p_sampling Proportion of the submatrix which contains the top n rows to sample.
# -partition_repeat Number of repeats for the random sampling.
# -partition_param Parameters for the partition method which are passed to ``...`` in a registered partition method. See `register_partition_methods` for detail.
# -anno A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation
#       to predicted subgroups will be tested.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -scale_rows Whether to scale rows. If it is ``TRUE``, scaling method defined in `register_partition_methods` is used.
# -verbose Whether print messages.
# -mc.cores Multiple cores to use.
# -.env An environment, internally used.
#
# == details
# The function performs analysis in following steps:
#
# - calculate scores for rows by top-value method,
# - for each top_n value, take top n rows,
# - randomly sample ``p_sampling`` rows from the top_n-row matrix and perform partitioning for ``partition_repeats`` times,
# - collect partitions from all partitions and calculate consensus partitions.
#
# == return
# A `ConsensusPartition-class` object. Simply type object in the interactive R session
# to see which functions can be applied on it.
#
# == seealso
# `run_all_consensus_partition_methods` runs consensus partition with multiple top-value methods
# and multiple partition methods.
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
# cp = consensus_partition(m, partition_repeat = 10, top_n = c(10, 20, 50))
# cp
consensus_partition = function(data,
	top_value_method = "ATC",
	top_n = seq(min(1000, round(nrow(data)*0.1)), 
		        min(5000, round(nrow(data)*0.5)), 
		        length.out = 5),
	partition_method = "skmeans",
	max_k = 6, 
	sample_by = "row",
	p_sampling = 0.8,
	partition_repeat = 50,
	partition_param = list(),
	anno = NULL,
	anno_col = NULL,
	scale_rows = NULL,
	verbose = TRUE,
	mc.cores = 1,
	.env = NULL) {

	if(missing(data)) {
		data = .env$data
	}

	if(max_k < 2) {
		stop_wrap("max_k should be no less than 2.")
	}

	t = system.time(res <- .consensus_partition(
		data = data,
		top_value_method = top_value_method,
		top_n = top_n,
		partition_method = partition_method,
		k = 2:max_k, 
		p_sampling = p_sampling,
		sample_by = sample_by,
		partition_repeat = partition_repeat,
		partition_param = partition_param,
		anno = anno,
		anno_col = anno_col,
		scale_rows = scale_rows,
		verbose = verbose,
		mc.cores = mc.cores,
		.env = .env))
	res@hash = digest(res)
	res@running_time = t[["elapsed"]]

	if(verbose) {
		tc = Sys.time()
		tf = format(tc + structure(t[["elapsed"]], units = "secs", class = "difftime") - tc)
		qqcat("* @{top_value_method}:@{partition_method} used @{tf}.\n")
	}
	return(res)
}

.consensus_partition = function(data,
	top_value_method = "MAD",
	top_n = seq(min(1000, round(nrow(data)*0.1)), 
		        min(5000, round(nrow(data)*0.5)), 
		        length.out = 5),
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
	mc.cores = 1,
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

	data = data[, .env$column_index, drop = FALSE]

	if(verbose) qqcat("* on a @{nrow(data)}x@{ncol(data)} matrix.\n")

	k = sort(k)
	l = k <= ncol(data)
	if(sum(l) != length(k)) {
		qqcat("* Following k (@{paste(k[!l], collapse=', ')}) are removed.\n")
	}
	k = k[l]
	if(length(k) == 0) {
		stop_wrap("There is no valid k.\n")
	}
	
	top_n = round(top_n)
	l = top_n <= nrow(data)
	if(sum(l) != length(top_n)) {
		qqcat("* Following top_n (@{paste(top_n[!l], collapse = ', ')}) are removed.\n")
	}
	top_n = top_n[l]
	if(length(top_n) == 0) {
		stop_wrap("There is no valid top_n.\n")
	}

	partition_fun = get_partition_method(partition_method, partition_param)

	get_top_value_fun = get_top_value_method(top_value_method)

	# also since one top value metric will be used for different partition methods,
	# we cache the top values for repetitive use
	if(is.null(.env$all_top_value_list)) {
		if(verbose) qqcat("* calculating @{top_value_method} values.\n")
		all_top_value = get_top_value_fun(data)
		all_top_value[is.na(all_top_value)] = -Inf
		.env$all_top_value_list = list()
		.env$all_top_value_list[[top_value_method]] = all_top_value
	} else if(is.null(.env$all_top_value_list[[top_value_method]])) {
		if(verbose) qqcat("* calculating @{top_value_method} values.\n")
		all_top_value = get_top_value_fun(data)
		all_top_value[is.na(all_top_value)] = -Inf
		.env$all_top_value_list[[top_value_method]] = all_top_value
	} else {
		if(verbose) qqcat("* @{top_value_method} values have already been calculated. Get from cache.\n")
		all_top_value = .env$all_top_value_list[[top_value_method]]
	}

	if(is.null(scale_rows)) {
		scale_rows = TRUE
	}
	if(scale_rows) {
		scale_method = attr(partition_fun, "scale_method")
		if("z-score" %in% scale_method) {
			if(verbose) cat("* rows are scaled before sent to partition, method: 'z-score' (x - mean)/sd\n")
			data = t(scale(t(data)))
		} else if("min-max" %in% scale_method) {
			if(verbose) cat("* rows are scaled before sent to partition, method: 'min-max' (x - min)/(max - min)\n")
			row_min = rowMins(data)
			row_max = rowMaxs(data)
			row_range = row_max - row_min
			data = (data - row_min)/row_range
		} else {
			scale_rows = FALSE
		}
	}

	# in case NA is produced in scaling
	l = apply(data, 1, function(x) any(is.na(x)))
	if(any(l)) {
		if(verbose) qqcat("* remove @{sum(l)} rows with NA values.\n")
		data = data[!l, , drop = FALSE]
		all_top_value = all_top_value[!l]
		l = top_n <= nrow(data)
		top_n = top_n[l]

		if(sum(l) != length(top_n)) {
			qqcat("* Following top_n (@{paste(top_n[!l], collapse = ', ')}) are removed.\n")
		}
		top_n = top_n[l]
		if(length(top_n) == 0) {
			stop_wrap("There is no valid top_n.\n")
		}
	}

	# now we do repetitive clustering
	param = data.frame(top_n = numeric(0), k = numeric(0), n_row = numeric(0))
	partition_list = list()

	sample_by = match.arg(sample_by)[1]
	for(i in seq_len(length(top_n))) {
		ind = order(all_top_value, decreasing = TRUE)[ 1:top_n[i] ]

		if(length(ind) > 5000) {
			ind = sample(ind, 5000)
			if(verbose) qqcat("* get top @{top_n[i]} (randomly sampled 5000) rows by @{top_value_method} method\n")
		} else {
			if(verbose) qqcat("* get top @{top_n[i]} rows by @{top_value_method} method\n")
		}

		if(verbose && mc.cores > 1) qqcat("  - @{partition_method} repeated for @{partition_repeat} times by @{sample_by}-sampling (p = @{p_sampling}) from top @{top_n[i]} rows (@{mc.cores} cores).\n")

		lt = mclapply(seq_len(partition_repeat), function(j) {

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
				if(interactive() && verbose && mc.cores == 1) cat(strrep("\b", 100))
				if(interactive() && verbose && mc.cores == 1) qqcat("  [k = @{y}] @{partition_method} repeated for @{j}@{ifelse(j %% 10 == 1, 'st', ifelse(j %% 10 == 2, 'nd', ifelse(j %% 10 == 3, 'rd', 'th')))} @{sample_by}-sampling (p = @{p_sampling}) from top @{top_n[i]} rows.")
				partition_list = c(partition_list, list(list(partition_fun(mat, y, column_ind_sub))))
				param = rbind(param, data.frame(top_n = top_n[i], k = y, n_row = nrow(mat), n_col = ncol(mat), stringsAsFactors = FALSE))
			}

			return(list(param = param, partition_list = partition_list))
		}, mc.cores = mc.cores)
		for(i in seq_along(lt)) {
			param = rbind(param, lt[[i]]$param)
			partition_list = c(partition_list, lt[[i]]$partition_list)
		}
		if(interactive() && verbose && mc.cores == 1) cat("\n")
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
		if(verbose) qqcat("* wrap results for k = @{y}\n")
		construct_consensus_object(param[l, ], partition_list[l], y, verbose = FALSE)
		
	})
	names(object_list) = as.character(k)

	rm(partition_list)
	gc(verbose = FALSE, reset = TRUE)

	## adjust class labels for different k should also be adjusted
	if(verbose) qqcat("* adjust class labels between different k.\n")
	reference_class = object_list[[1]]$class_df$class
	for(i in seq_along(k)[-1]) {
		class_df = object_list[[i]]$class_df
    	class = class_df[, "class"]

    	map = relabel_class(class, reference_class, full_set = 1:(k[i]))
    	map2 = structure(names(map), names = map)
    	# unmapped = setdiff(as.character(1:k[i]), map)
    	# map = c(map, structure(unmapped, names = unmapped))
    	# map2 = c(map2, structure(unmapped, names = unmapped))
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
		anno = anno[.env$column_index, , drop = FALSE]
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

	res = ConsensusPartition(object_list = object_list, k = k, n_partition = partition_repeat * length(top_n) * length(k),  
		partition_method = partition_method, top_value_method = top_value_method, top_n = top_n,
		anno = anno, anno_col = anno_col, scale_rows = scale_rows, sample_by = sample_by,
		column_index = .env$column_index, .env = .env)

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
	qqcat("  On a matrix with @{nrow(object@.env$data)} rows and @{length(object@column_index)} columns.\n")
	top_n_str = object@top_n
	qqcat("  Top rows (@{paste(top_n_str, collapse = ', ')}) are extracted by '@{object@top_value_method}' method.\n")
	qqcat("  Subgroups are detected by '@{object@partition_method}' method.\n")
	qqcat("  Performed in total @{object@n_partition} partitions by @{object@sample_by} resampling.\n")
	best_k = suggest_best_k(object)
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
# Plot the empirical cumulative distribution curve (ECDF) of the consensus matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -... Other arguments.
#
# == details
# It plots ECDF curve for each k.
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
# data(cola_rl)
# plot_ecdf(cola_rl["sd", "hclust"])
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
# Several plots for determining the optimized number of partitions
#
# == param
# -object A `ConsensusPartition-class` object.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == details
# There are following plots made:
#
# - ECDF of the consensus matrix under each k, made by `plot_ecdf,ConsensusPartition-method`,
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
# data(cola_rl)
# select_partition_number(cola_rl["sd", "hclust"])
setMethod(f = "select_partition_number",
	signature = "ConsensusPartition",
	definition = function(object, all_stats = FALSE) {
	op = par(no.readonly = TRUE)

	m = get_stats(object, all_stats = all_stats)
	if(all_stats) STAT_USED = STAT_ALL
	m = m[, colnames(m) %in% c(STAT_USED, "area_increased", "Rand", "Jaccard"), drop = FALSE]
	nm = colnames(m)
	n = ncol(m)
	nr = floor(sqrt(n)+1)

	par(mfrow = c(nr, ceiling((n+1)/nr)), mar = c(4, 4, 1, 1))

	plot_ecdf(object, lwd = 1)

	best_k = suggest_best_k(object)
	if(is.na(best_k)) best_k = -1

	for(i in seq_len(ncol(m))) {
		l = object@k == best_k
		plot(object@k, m[, i], type = "b", xlab = "k", ylab = nm[i])
		if(any(l)) {
			points(object@k, m[, i],
				pch = ifelse(l, 16, 1),
				col = ifelse(l, "red", "black"))
		}
	}

	par(xpd = NA, mar = c(4, 2, 1, 1))
	plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ann = FALSE)
	text(x = 0, y = 1,
"Suggested rules:
Jaccard index < 0.95;
If 1-PAC >= 0.90,
  take the maximum k;
else take the k with higest votes of
  1. max 1-PAC,
  2. max mean_silhouette,
  3. max concordance.
", cex = 1.2, adj = c(0, 1))

	legend("topright", pch = 16, col = "Red", legend = "best k", cex = 1.5)

	par(op)
})


# == title
# Heatmap for the consensus matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions.
# -internal Used internally.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `consensus_partition` or `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_row_names Whether plot row names on the consensus heatmap (which are the column names in the original matrix)
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
# - predicted classes.
# - consensus matrix.
# - more annotations if provided as ``anno``
#
# One thing that is very important to note is that since we already know the consensus classes from consensus
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
# data(cola_rl)
# consensus_heatmap(cola_rl["sd", "hclust"], k = 3)
setMethod(f = "consensus_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, internal = FALSE,
	anno = get_anno(object), anno_col = get_anno_col(object), 
	show_row_names = FALSE, ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

	class_df = get_classes(object, k)
	class_ids = class_df$class

	consensus_mat = get_consensus(object, k)

	mat_col_od = column_order_by_group(factor(class_ids, levels = sort(unique(class_ids))), consensus_mat)

	membership_mat = get_membership(object, k)

	ht_list = Heatmap(membership_mat, name = "Prob", cluster_columns = FALSE, show_row_names = FALSE,
		width = unit(5, "mm")*k, col = colorRamp2(c(0, 1), c("white", "red")),
		show_column_names = !internal) + 
	Heatmap(class_df$silhouette, name = "Silhouette", width = unit(5, "mm"),
		show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "purple")),
		show_column_names = !internal) +
	Heatmap(class_ids, name = "Class", col = cola_opt$color_set_2,
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
		ht_list = ht_list + rowAnnotation(rn = anno_text(rownames(consensus_mat)))
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
# -k Number of partitions.
# -internal Used internally.
# -anno A data frame of annotations for the original matrix columns. 
#       By default it uses the annotations specified in `consensus_partition` or `run_all_consensus_partition_methods`.
# -anno_col A list of colors (color is defined as a named vector) for the annotations. If ``anno`` is a data frame,
#       ``anno_col`` should be a named list where names correspond to the column names in ``anno``.
# -show_column_names Whether show column names in the heatmap (which is the column name in the original matrix).
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
# data(cola_rl)
# membership_heatmap(cola_rl["sd", "hclust"], k = 3)
setMethod(f = "membership_heatmap",
	signature = "ConsensusPartition",
	definition = function(object, k, internal = FALSE, 
	anno = get_anno(object), anno_col = get_anno_col(object),
	show_column_names = FALSE, ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

	class_df = get_classes(object, k)
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
		show_column_names = show_column_names,
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

# == title
# Visualize column after dimension reduction
#
# == description
# Visualize samples (the matrix columns) after dimension reduction
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions.
# -top_n Top n rows to use. By default it uses all rows in the original matrix.
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -control A list of parameters for `Rtsne::Rtsne` or `umap::umap`.
# -internal Internally used.
# -nr If number of matrix rows is larger than this value, random ``nr`` rows are used.
# -silhouette_cutoff Cutoff of silhouette score. Data points with values less
#        than it will be mapped with cross symbols.
# -remove Whether to remove columns which have less silhouette scores than
#        the cutoff.
# -scale_rows Whether perform scaling on matrix rows.
# -verbose Whether print messages.
# -... Other arguments.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# dimension_reduction(cola_rl["sd", "kmeans"], k = 3)
setMethod(f = "dimension_reduction",
	signature = "ConsensusPartition",
	definition = function(object, k, top_n = NULL,
	method = c("PCA", "MDS", "t-SNE", "UMAP"), 
	control = list(),
	internal = FALSE, nr = 5000,
	silhouette_cutoff = 0.5, remove = FALSE,
	scale_rows = TRUE, verbose = TRUE, ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

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
	data = object@.env$data[, object@column_index, drop = FALSE]

	if(!is.null(top_n)) {
		top_n = min(c(top_n, nrow(data)))
		all_value = object@.env$all_top_value_list[[object@top_value_method]]
		ind = order(all_value)[1:top_n]
		if(length(ind) > 5000) ind = sample(ind, 5000)
		data = data[ind, , drop = FALSE]
	} else {
		top_n = nrow(data)
		if(nrow(data) > 5000) data = data[sample(1:nrow(data), 5000), , drop = FALSE]
	}

	class_df = get_classes(object, k)

	l = class_df$silhouette >= silhouette_cutoff
	
	op = par(c("mar", "xpd"))
	par(mar = c(4.1, 4.1, 4.1, 6), xpd = NA)
	on.exit(par(op))

	class_level = sort(as.character(unique(class_df$class)))
	n_class = length(class_level)
		
	if(remove) {
		dimension_reduction(data[, l], pch = 16, col = cola_opt$color_set_2[as.character(class_df$class[l])],
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (silhouette > @{silhouette_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
				pch = c(rep(16, n_class), 0),
				col = c(cola_opt$color_set_2[class_level], "white"), xjust = 0, yjust = 0.5,
				title = "Class", title.adj = 0.1, bty = "n",
				text.col = c(rep("black", n_class), "white"))
		}
	} else {
		dimension_reduction(data, pch = ifelse(l, 16, 4), col = cola_opt$color_set_2[as.character(class_df$class)],
			cex = 1, main = qq("@{method} on @{top_n} rows with highest @{object@top_value_method} scores@{ifelse(scale_rows, ', rows are scaled', '')}\n@{sum(l)}/@{length(l)} confident samples (silhouette > @{silhouette_cutoff})"),
			method = method, control = control, scale_rows = scale_rows, nr = nr, internal = internal, verbose = verbose, ...)
		if(!internal) {
			if(any(!l)) {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
						pch = c(rep(16, n_class), 4),
						col = c(cola_opt$color_set_2[class_level], "black"), xjust = 0, yjust = 0.5,
						title = "Class", title.adj = 0.1, bty = "n")
			} else {
				legend(x = par("usr")[2], y = mean(par("usr")[3:4]), legend = c(paste0("group", class_level), "ambiguous"), 
					pch = c(rep(16, n_class), 0),
					col = c(cola_opt$color_set_2[class_level], "white"), xjust = 0, yjust = 0.5,
					title = "Class", title.adj = 0.1, bty = "n",
					text.col = c(rep("black", n_class), "white"))
			}
		}
	}
})


# == title
# Visualize columns after dimension reduction
#
# == param
# -object A numeric matrix.
# -method Which method to reduce the dimension of the data. ``MDS`` uses `stats::cmdscale`,
#         ``PCA`` uses `stats::prcomp`. ``t-SNE`` uses `Rtsne::Rtsne`. ``UMAP`` uses
#         `umap::umap`.
# -pc Which two principle components to visualize
# -control A list of parameters for `Rtsne::Rtsne` or `umap::umap`.
# -pch Ahape of points.
# -col Color of points.
# -cex Aize of points.
# -main Title of the plot.
# -scale_rows Whether perform scaling on matrix rows.
# -nr If number of matrix rows is larger than this value, random ``nr`` rows are used.
# -internal Internally used.
# -verbose Whether print messages.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "dimension_reduction",
	signature = "matrix",
	definition = function(object, 
	pch = 16, col = "black", cex = 1, main = "",
	method = c("PCA", "MDS", "t-SNE", "UMAP"),
	pc = NULL, control = list(), 
	scale_rows = TRUE, nr = 5000,
	internal = FALSE, verbose = TRUE) {

	data = object
	if(nrow(data) > nr) {
		data = data[sample(1:nrow(data), nr), , drop = FALSE]
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

	method = match.arg(method)[1]

	if(is.null(pc)) {
		if(method == "PCA") {
			pc = 1:2
		} else if(method %in% c("UMAP", "t-SNE")) {
			pc = seq_len(min(c(10, ncol(data))))
		}
	} else if(length(pc) == 1) {
		pc = 1:pc
	}

	if(method %in% c("UMAP", "t-SNE")) {
		main = paste0(main, qq(", with @{length(pc)} PCs"))
	}

	if(scale_rows) data = t(scale(t(data)))
	l = apply(data, 1, function(x) any(is.na(x)))
	data = data[!l, ]

	if(internal) {
		omar = par("mar")
		par(mar = c(4.1, 4.1, 1, 1))
		on.exit(par(mar = omar))

		main = NULL
	}

	if(method == "MDS") {
		loc = cmdscale(dist(t(data)))
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "Coordinate 1", ylab = "Coordinate 2")
	} else if(method == "PCA") {
		fit = prcomp(t(data))
		sm = summary(fit)
		prop = sm$importance[2, pc]
		loc = fit$x[, pc]

		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = qq("PC@{pc[1]} (@{round(prop[1]*100)}%)"), ylab = qq("PC@{pc[2]} (@{round(prop[2]*100)}%)"))
	} else if(method == "t-SNE") {
		fit = prcomp(t(data))
		sm = summary(fit)
		loc = fit$x[, pc]

		param = list(X = loc)
		param = c(param, control)
		fit = do.call(Rtsne::Rtsne, param)
		loc = fit$Y
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "t-SNE 1", ylab = "t-SNE 2")
	} else if(method == "UMAP") {

		fit = prcomp(t(data))
		sm = summary(fit)
		loc = fit$x[, pc]

		param = list(d = loc)
		if(!"config" %in% names(control)) {
			# reset n_neighbors for small dataset
			control$config = umap::umap.defaults
			control$config$n_neighbors = min(umap::umap.defaults$n_neighbors, round(ncol(data)/2))
		} else {
			if(!"n_neighbors" %in% names(control$config)) {
				# reset n_neighbors for small dataset
				control$config$n_neighbors = min(umap::umap.defaults$n_neighbors, round(ncol(data)/2))
			}
		}
		param = c(param, control)
		fit = do.call(umap::umap, param)
		loc = fit$layout
		plot(loc, pch = pch, col = col, cex = cex, main = main, xlab = "UMAP 1", ylab = "UMAP 2")
	}
})

