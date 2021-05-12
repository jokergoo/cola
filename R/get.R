
# == title
# Get parameters
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -unique Whether to apply `base::unique` to rows of the returned data frame.
#
# == details
# It is mainly used internally.
#
# == value
# A data frame of parameters corresponding to the current k. In the data frame, each row corresponds
# to a partition run.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_param(obj)
# get_param(obj, k = 2)
# get_param(obj, unique = FALSE)
setMethod(f = "get_param",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k, unique = TRUE) {
	ind = which(object@k %in% k)
	if(unique) {
		df = unique(do.call("rbind", lapply(ind, function(i) object@object_list[[i]]$param)))
	} else {
		df = do.call("rbind", lapply(ind, function(i) object@object_list[[i]]$param))
	}
	rownames(df) = NULL
	return(df)
})

# == title
# Get consensus matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
#
# == details
# For row i and column j in the consensus matrix, the value of corresponding x_ij
# is the probability of sample i and sample j being in the same group from all partitions.
#
# == value
# A consensus matrix corresponding to the current k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_consensus(obj, k = 2)
setMethod(f = "get_consensus",
	signature = "ConsensusPartition",
	definition = function(object, k) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	i = which(object@k == k)
	object@object_list[[i]]$consensus
})

# == title
# Get membership matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -each Whether to return the percentage membership matrix which is summarized from all partitions
#       or the individual membership in every single partition run.
#
# == details
# If ``each == FALSE``, the value in the membership matrix is the probability
# to be in one subgroup, while if ``each == TRUE``, the membership matrix contains the 
# subgroup labels for every single partitions which are from randomly sampling from the original matrix.
#
# The percent membership matrix is calculated by `clue::cl_consensus`.
#
# == value
# - If ``each == FALSE``, it returns a membership matrix where rows correspond to the columns from the subgroups.
# - If ``each == TRUE``, it returns a membership matrix where rows correspond to the columns from the original matrix.
#
# == seealso
# `get_membership,ConsensusPartitionList-method` summarizes membership from partitions from all combinations
# of top-value methods and partitioning methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_membership(obj, k = 2)
# get_membership(obj, k = 2, each = TRUE)
setMethod(f = "get_membership",
	signature = "ConsensusPartition",
	definition = function(object, k, each = FALSE) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	i = which(object@k == k)
	if(each) {
		object@object_list[[i]]$membership_each
	} else {
		object@object_list[[i]]$membership
	}
})

# == title
# Get membership matrix
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of subgroups.
#
# == detail
# The membership matrix (the probability of each sample to be in one subgroup, if assuming columns represent samples) is inferred
# from the consensus partition of every combination of methods, weighted by the mean silhouette score of the partition
# for each method. So methods which give unstable partitions have lower weights 
# when summarizing membership matrix from all methods.
# 
# == value
# A membership matrix where rows correspond to the columns in the original matrix.
#
# == seealso
# `get_membership,ConsensusPartition-method` returns membership matrix for a single top-value method and partitioning method.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# get_membership(golub_cola, k = 2)
setMethod(f = "get_membership",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	lt = object@consensus_class[[as.character(k)]]
	lt$membership
})



# == title
# Get statistics
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups. The value can be a vector.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == details
# The statistics are:
#
# -1-PAC 1 - proportion of ambiguous clustering, calculated by `PAC`.
# -mean_silhouette The mean silhouette score. See https://en.wikipedia.org/wiki/Silhouette_(clustering) .
# -concordance The mean probability that each partition fits the consensus partition, calculated by `concordance`.
# -area_increased The increased area under eCDF (the empirical cumulative distribution function) curve to the previous k.
# -Rand This is the percent of pairs of samples that are both in a same cluster or both are not 
#       in a same cluster in the partition of ``k`` and ``k-1``. See https://en.wikipedia.org/wiki/Rand_index .
# -Jaccard The ratio of pairs of samples are both in a same cluster in the partition of ``k`` and ``k-1`` and the pairs
#          of samples are both in a same cluster in the partition ``k`` or ``k-1``.
#
# == value
# A matrix of partition statistics.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_stats(obj)
# get_stats(obj, k = 2)
setMethod(f = "get_stats",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k, all_stats = FALSE) {

	if(all_stats) STAT_USED = STAT_ALL

	m = matrix(nrow = length(object@k), ncol = length(STAT_USED) + 3)
	rownames(m) = object@k
	colnames(m) = c(STAT_USED, "area_increased", "Rand", "Jaccard")
	for(i in seq_along(object@k)) {
		m[i, ] = unlist(object@object_list[[i]]$stat[colnames(m)])
	}
	m = m[as.character(k), , drop = FALSE]
	# m = cbind(m, separation_rate = sapply(k, function(i) separation_rate(object, i)))
	m = cbind(k = k, m)
	return(m)
})


# == title
# Get statistics
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of subgroups. The value can only be a single value.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == value
# A matrix of partition statistics for a selected k. Rows in the 
# matrix correspond to combinations of top-value methods and partitioning methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# get_stats(golub_cola, k = 2)
setMethod(f = "get_stats",
	signature = "ConsensusPartitionList",
	definition = function(object, k, all_stats = FALSE) {

	if(all_stats) STAT_USED = STAT_ALL

	if(missing(k)) stop_wrap("k needs to be provided.")
	m = matrix(nrow = length(object@list), ncol = length(STAT_USED) + 3)
	rownames(m) = names(object@list)
	colnames(m) = c(STAT_USED, "area_increased", "Rand", "Jaccard")
	for(i in seq_along(object@list)) {
		ik = which(object@list[[i]]@k == k)
		m[i, ] = unlist(object@list[[i]]@object_list[[ik]]$stat[colnames(m)])
	}
	m = cbind(k = rep(k, nrow(m)), m)
	return(m)
})


# == title
# Get subgroup labels
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
#
# == return
# A data frame with subgroup labels and other columns which are entropy of the percent membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its group.
#
# If ``k`` is not specified, it returns a data frame with subgroup labels from all k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_classes(obj, k = 2)
# get_classes(obj)
setMethod(f = "get_classes",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k) {
	if(length(k) == 1) {
		i = which(object@k == k)
		object@object_list[[i]]$class_df
	} else {
		df = do.call("cbind", lapply(k, function(i) object@object_list[[as.character(i)]]$class_df$class))
		colnames(df) = paste0("k=", k)
		rownames(df) = rownames(object@object_list[[1]]$class_df)
		df
	}
})

# == title
# Get subgroup labels
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of subgroups.
# 
# == details 
# The subgroup labels are inferred by merging partitions from all methods
# by weighting the mean silhouette scores in each method. 
#
# == return
# A data frame with subgroup labels and other columns which are entropy of the percent membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its group.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# get_classes(golub_cola, k = 2)
setMethod(f = "get_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	lt = object@consensus_class[[as.character(k)]]
	df = lt$class_df
	rownames(df) = colnames(object)
	df
})

# == title
# Recalculate statistics in the ConsensusPartitionList object
#
# == param
# -rl A `ConsensusPartitionList-class` object.
#
# == details
# It updates the ``stat`` slot in the ConsensusPartitionList object, used internally.
#
recalc_stats = function(rl) {
	for(j in seq_along(rl@list)) {
		for(i in seq_along(rl@list[[j]]@k)) {
			rl@list[[j]]@object_list[[i]]$stat["PAC"] = NULL
			rl@list[[j]]@object_list[[i]]$stat["cophcor"] = NULL
			rl@list[[j]]@object_list[[i]]$stat["stability"] = NULL
			rl@list[[j]]@object_list[[i]]$stat["flatness"] = NULL
			rl@list[[j]]@object_list[[i]]$stat["FCC"] = NULL

			consensus_mat = rl@list[[j]]@object_list[[i]]$consensus

			if(length(rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"])*0.05 > 1) {
				l = rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"] >= quantile(rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"], 0.05)
			} else {
				l = rep(TRUE, length(rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"]))
			}

			rl@list[[j]]@object_list[[i]]$stat[["1-PAC"]] = 1 - PAC(consensus_mat[l, l, drop = FALSE])
			rl@list[[j]]@object_list[[i]]$stat[["cophcor"]] = cophcor(consensus_mat)
			rl@list[[j]]@object_list[[i]]$stat[["aPAC"]] = aPAC(consensus_mat)
			rl@list[[j]]@object_list[[i]]$stat[["FCC"]] = FCC(consensus_mat[l, l, drop = FALSE])
		}

	}
	return(rl)
}

# == title
# Suggest the best number of subgroups
#
# == param
# -object A `ConsensusPartition-class` object.
# -jaccard_index_cutoff The cutoff for Jaccard index for comparing to previous k.
# -help Whether to print help message.
#
# == details
# The best k is selected according to following rules:
#
# - All k with Jaccard index larger than 0.95 are removed because increasing k does not
#   provide enough extra information. If all k are removed, it is marked as no
#   subgroup is detected. 
# - For all k with 1-PAC score larger than 0.9, the
#   maximal k is taken as the best k, and other k are marked as optional k. 
# - If it does not fit the second rule. The k with the maximal vote of the highest
#   1-PAC score, highest mean silhouette, and highest concordance is taken as
#   the best k.
#
# Additionally, if 1-PAC for the best k is larger than 0.9 (10\% ambiguity for
# the partition), cola marks it as a stable partition. It should be noted that
# it is difficult to find the best k deterministically, we encourage users to
# compare results for all k and determine a proper one which best explain
# their studies.
#
# When k > 6, only the third rule is applied because 1-PAC does not work well for larger k.
#
# == see also
# The selection of the best k can be visualized by `select_partition_number`.
#
# == value
# The best k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# suggest_best_k(obj)
setMethod(f = "suggest_best_k",
	signature = "ConsensusPartition",
	definition = function(object, jaccard_index_cutoff = 0.95, stable_PAC = 0.1, help = TRUE) {

	# if(verbose) {
	# 	cat("This function only suggests the best k. It is recommended that users look ...")
	# }

	stat = get_stats(object)
	stat = stat[stat[, "Jaccard"] < jaccard_index_cutoff, , drop = FALSE]

	if(nrow(stat) == 0) {
		return(NA)
	}

	if(nrow(stat) == 1) {
		return(stat[, "k"])
	}

	if(help) {
		message_wrap("The best k suggested by this function might not reflect the real subgroups in the data (especially when you expect a large best k). It is recommended to directly look at the plots from select_partition_number() or other related plotting functions.")
	}

	if(min(stat[, "Jaccard"]) >= 0.75) {
		dec = c(which.max(stat[, "mean_silhouette"]))
		
		tb = table(dec)
		max_ind = order(tb, as.numeric(names(tb)), decreasing = TRUE)[1]
		x = rownames(stat)[as.numeric(names(tb[max_ind]))]
		
		return(as.numeric(x))
	}

	l = stat[, "1-PAC"] >= 1 - stable_PAC
	if(sum(l) == 1) {
		return(as.numeric(rownames(stat)[l]))
	} else if(sum(l) > 1) {
		best_k = as.numeric(rownames(stat)[l])
		best_k = structure(max(best_k), optional = setdiff(best_k, max(best_k)))
		return(best_k)
	}

	l = stat[-nrow(stat), "area_increased"] >= stat[-1, "area_increased"]
	if(sum(l) == 0) {
		return(NA)
	}
	stat = stat[l, , drop = FALSE]

	dec = c(which.max(stat[, "1-PAC"]),
		    which.max(stat[, "mean_silhouette"]),
		    which.max(stat[, "concordance"]))
	
	tb = table(dec)
	max_ind = order(tb, as.numeric(names(tb)), decreasing = TRUE)[1]
	x = rownames(stat)[as.numeric(names(tb[max_ind]))]
	
	as.numeric(x)
})


# == title
# Test whether the current k is the best/optional k
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -...  Pass to `suggest_best_k,ConsensusPartition-method`.
#
# == details
# Optional best k is also assigned as ``TRUE``.
#
# == value
# Logical scalar.
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# is_best_k(obj, k = 2)
# is_best_k(obj, k = 3)
setMethod(f = "is_best_k",
	signature = "ConsensusPartition",
	definition = function(object, k, ...) {

	best_k = suggest_best_k(object, ..., help = FALSE)
	if(is.na(best_k)) {
		return(FALSE)
	} else {
		if(k == best_k) {
			return(TRUE)
		} else {
			if(k %in% attr(best_k, "optional")) {
				return(TRUE)
			}
		}
	}
	return(FALSE)
})


# == title
# Test whether the current k corresponds to a stable partition
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of subgroups.
# -...  Pass to `suggest_best_k,ConsensusPartition-method`.
#
# == details
# if 1-PAC for the k is larger than 0.9 (10\% ambiguity for
# the partition), cola marks it as a stable partition.
#
# == value
# Logical scalar.
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# is_stable_k(obj, k = 2)
# is_stable_k(obj, k = 3)
setMethod(f = "is_stable_k",
	signature = "ConsensusPartition",
	definition = function(object, k, ...) {

	is_best_k(object, k, ...) & get_stats(object, k = k)[, "1-PAC"] >= 0.9
})

# == title
# Suggest the best number of subgroups
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -jaccard_index_cutoff The cutoff for Jaccard index for comparing to previous k.
#
# == details
# It basically gives the best k for each combination of top-value method and partitioning method by calling `suggest_best_k,ConsensusPartition-method`.
#
# 1-PAC score higher than 0.95 is treated as very stable partition (marked by ``**``) and higher than 0.9 is treated as stable partition (marked by ``*``).
#
# == value
# A data frame with the best k and other statistics for each combination of methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# suggest_best_k(golub_cola)
setMethod(f = "suggest_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object, jaccard_index_cutoff = 0.95) {

	best_k = NULL
	stability = NULL
	mean_silhouette = NULL
	concordance = NULL
	optional_k = NULL
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			k = suggest_best_k(obj, jaccard_index_cutoff, help = FALSE)
			best_k[nm] = k
			op = attr(k, "optional")
			if(is.null(op)) {
				optional_k[nm] = ""
			} else {
				optional_k[nm] = paste(op, collapse = ",")
			}
			if(is.na(best_k[nm])) {
				stability[nm] = NA
				mean_silhouette[nm] = NA
				concordance[nm] = NA
			} else {
				stat = get_stats(obj, k = best_k[nm])
				stability[nm] = stat[1, "1-PAC"]
				mean_silhouette[nm] = stat[1, "mean_silhouette"]
				concordance[nm] = stat[1, "concordance"]
			}
		}
	}

	tb = data.frame(best_k = best_k,
		"1-PAC" = stability,
		mean_silhouette = mean_silhouette,
		concordance = concordance,
		check.names = FALSE)

	rntb = rownames(tb)
	l = tb$`1-PAC` >= 0.9 & !is.na(tb$best_k)

	tb = cbind(tb, ifelse(l, ifelse(tb$`1-PAC` <= 0.95, "*", "**"), ""), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""
	if(any(optional_k != "")) {
		tb$optional_k = optional_k
	}
	return(tb[order(!is.na(tb[, "best_k"]) + 0, tb[, "1-PAC"], decreasing = TRUE), , drop = FALSE])
})

# == title
# Test whether the current k is the best/optional k
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of subgroups.
# -... Pass to `suggest_best_k,ConsensusPartitionList-method`.
#
# == details
# It tests on the partitions for every method.
#
# == value
# Logical vector.
#
# == example
# data(golub_cola)
# is_best_k(golub_cola, k = 3)
setMethod(f = "is_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object, k, ...) {

	x = sapply(object@list, function(x) is_best_k(x, k, ...))
	names(x) = names(object@list)
	return(x)
})

# == title
# Test whether the current k corresponds to a stable partition
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of subgroups.
# -... Pass to `suggest_best_k,ConsensusPartitionList-method`.
#
# == details
# It tests on the partitions for every method.
#
# == value
# Logical vector
#
# == example
# data(golub_cola)
# is_stable_k(golub_cola, k = 3)
setMethod(f = "is_stable_k",
	signature = "ConsensusPartitionList",
	definition = function(object, k, ...) {

	x = sapply(object@list, function(x) is_stable_k(x, k, ...))
	names(x) = names(object@list)
	return(x)
})

# == title
# Get the original matrix
#
# == param
# -object A `ConsensusPartitionList-class` object.
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# get_matrix(golub_cola)
setMethod(f = "get_matrix",
	signature = "ConsensusPartitionList",
	definition = function(object) {
	object@.env$data
})

# == title
# Get the original matrix
#
# == param
# -object A `ConsensusPartition-class` object.
# -full Whether to extract the complete original matrix.
# -include_all_rows Internally used.
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# obj = golub_cola["ATC", "skmeans"]
# get_matrix(obj)
setMethod(f = "get_matrix",
	signature = "ConsensusPartition",
	definition = function(object, full = FALSE, include_all_rows = FALSE) {
	if(!full) {
		if(include_all_rows) {
			object@.env$data[, object@column_index, drop = FALSE]
		} else {
			object@.env$data[object@row_index, object@column_index, drop = FALSE]
		}
	} else {
		object@.env$data
	}
})

# == title
# Get annotations
#
# == param
# -object A `ConsensusPartitionList-class` object.
#
# == value
# A data frame if ``anno`` was specified in `run_all_consensus_partition_methods`, or else ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno",
	signature = "ConsensusPartitionList",
	definition = function(object) {
	object@list[[1]]@anno
})

# == title
# Get annotations
#
# == param
# -object A `ConsensusPartition-class` object.
#
# == value
# A data frame if ``anno`` was specified in `run_all_consensus_partition_methods` or `consensus_partition`, or else ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno",
	signature = "ConsensusPartition",
	definition = function(object) {
	object@anno
})

# == title
# Get annotation colors
#
# == param
# -object A `ConsensusPartitionList-class` object.
#
# == value
# A list of color vectors or else ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno_col",
	signature = "ConsensusPartitionList",
	definition = function(object) {
	object@list[[1]]@anno_col
})

# == title
# Get annotation colors
#
# == param
# -object A `ConsensusPartition-class` object.
#
# == value
# A list of color vectors or else ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno_col",
	signature = "ConsensusPartition",
	definition = function(object) {
	object@anno_col
})
