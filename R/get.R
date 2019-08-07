
# == title
# Get parameters
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions.
# -unique Whether apply `base::unique` to rows of the returned data frame.
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
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
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
# -k Number of partitions.
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
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
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
# -k Number of partitions.
# -each Whether return the percentage membership matrix which is summarized from all partitions
#       or the individual membership in every random partition.
#
# == details
# If ``each == FALSE``, the value in the membership matrix is the probability
# to be in one class, while if ``each == TRUE``, the membership matrix contains the 
# class labels for every single partitions which are from randomly sampling subset
# of rows in the matrix.
#
# The percent membership matrix is calculated by `clue::cl_consensus`.
#
# == value
# If ``each == TRUE``, it returns a membership matrix where rows correspond to the columns in the original matrix.
#
# == seealso
# `get_membership,ConsensusPartitionList-method` summarizes membership from partitions from all combinations
# of top-value methods and partition methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
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
# -k Number of partitions.
#
# == detail
# The membership matrix (the probability of each sample to be in one group, if assuming columns represent samples) is inferred
# from the consensus partition of every combination of methods, weighted by the mean silhouette score of the partition
# for each method. So methods which give instable partitions have lower weights 
# when summarizing membership matrix from all methods.
# 
# == value
# A membership matrix where rows correspond to the columns in the original matrix.
#
# == seealso
# `get_membership,ConsensusPartition-method` returns membership matrix for a single top-value method and partition method.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# get_membership(cola_rl, k = 2)
setMethod(f = "get_membership",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	lt = object@consensus_class[[as.character(k)]]
	lt$membership
})



# == title
# Get statistics for the consensus partition
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions. The value can be a vector.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == details
# The statistics are:
#
# -PAC proportion of ambiguous clustering, calculated by `PAC`.
# -mean_silhouette the mean silhouette score. See https://en.wikipedia.org/wiki/Silhouette_(clustering) .
# -concordance the mean probability that each partition fits the consensus partition, calculated by `concordance`.
# -area_increased the increased area under ECDF (the empirical cumulative distribution function curve) to the previous k.
# -Rand the Rand index which is the percent of pairs of samples that are both in a same cluster or both are not 
#       in a same cluster in the partition of ``k`` and ``k-1``. See https://en.wikipedia.org/wiki/Rand_index .
# -Jaccard the ratio of pairs of samples are both in a same cluster in the partition of ``k`` and ``k-1`` and the pairs
#          of samples are both in a same cluster in the partition ``k`` or ``k-1``.
#
# == value
# A matrix of partition statistics.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
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
# Get statistics for consensus partitions from all methods
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of partitions. The value can only be a single value.
# -all_stats Whether to show all statistics that were calculated. Used internally.
#
# == value
# A matrix of partition statistics for a selected k. Rows in the 
# matrix correspond to combinations of top-value methods and partition methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# get_stats(cola_rl, k = 2)
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
# Get class IDs from the ConsensusPartition object
#
# == param
# -object A `ConsensusPartition-class` object.
# -k Number of partitions.
#
# == return
# A data frame with class IDs and other columns which are entropy of the percent membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its group.
#
# If ``k`` is not specified, it returns a data frame with class IDs from every k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
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
# Get class IDs from the ConsensusPartitionList object
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -k Number of partitions.
# 
# == details 
# The class IDs are inferred by merging partitions from all methods
# by weighting the mean silhouette scores in each method. 
#
# == return
# A data frame with class IDs and other columns which are entropy of the percent membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its group.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# get_classes(cola_rl, k = 2)
setMethod(f = "get_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	if(missing(k)) stop_wrap("k needs to be provided.")
	lt = object@consensus_class[[as.character(k)]]
	lt$class_df
})

# == title
# Recalculate statistics in the ConsensusPartitionList object
#
# == param
# -rl A `ConsensusPartitionList-class` object.
#
# == details
# It updates the statistics slot in the ConsensusPartitionList object, used internally.
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

			l = rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"] >= quantile(rl@list[[j]]@object_list[[i]]$class_df[, "silhouette"], 0.05)
			rl@list[[j]]@object_list[[i]]$stat[["1-PAC"]] = 1 - PAC(consensus_mat[l, l, drop = FALSE])
			rl@list[[j]]@object_list[[i]]$stat[["cophcor"]] = cophcor(consensus_mat)
			rl@list[[j]]@object_list[[i]]$stat[["aPAC"]] = aPAC(consensus_mat)
			rl@list[[j]]@object_list[[i]]$stat[["FCC"]] = FCC(consensus_mat[l, l, drop = FALSE])
		}

	}
	return(rl)
}

# == title
# Suggest the best number of partitions
#
# == param
# -object A `ConsensusPartition-class` object.
# -rand_index_cutoff The cutoff for Rand index compared to previous k.
#
# == details
# The best k is selected according to following rules:
#
# 1. k with rand index larger than ``rand_index_cutoff`` are removed. If all k are removed, the best k is defined as ``NA``.
# 2. If there are some k having ``1-PAC`` larger than 0.9, the largest k is selected as the best k.
# 3. If it does not fit rule 2, the k with highest vote of highest 1-PAC, mean_silhouette and concordance scores is
#    selected as the best k.
#
# `suggest_best_k` function only gives suggestion of selecting
# a reasonable best k. Users still need to look at the plots (e.g. by `select_partition_number` or `consensus_heatmap` functions), or even
# by checking whether the subgrouping gives a reasonable signatures by `get_signatures`, to pick a reasonable k that best explains their study. 
#
# The best k with 1-PAC larger than 0.9 is treated as a stable partition.
# 
# == value
# The best k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
# suggest_best_k(obj)
setMethod(f = "suggest_best_k",
	signature = "ConsensusPartition",
	definition = function(object, rand_index_cutoff = 0.95) {

	stat = get_stats(object)
	stat = stat[stat[, "Rand"] < rand_index_cutoff, , drop = FALSE]

	if(nrow(stat) == 0) {
		return(NA)
	}

	l = stat[, "1-PAC"] >= 0.9
	if(sum(l) == 1) {
		return(as.numeric(rownames(stat)[l]))
	} else if(sum(l) > 1) {
		best_k = as.numeric(rownames(stat)[l])
		best_k = structure(max(best_k), optional = setdiff(best_k, max(best_k)))
		return(best_k)
	}

	dec = c(which.max(stat[, "1-PAC"]),
		    which.max(stat[, "mean_silhouette"]),
		    which.max(stat[, "concordance"]))
	
	tb = table(dec)
	x = rownames(stat)[as.numeric(names(tb[which.max(tb)[1]]))]
	
	as.numeric(x)
})


# == title
# Suggest the best number of partitions
#
# == param
# -object A `ConsensusPartitionList-class` object.
# -rand_index_cutoff The cutoff for Rand index compared to previous k.
#
# == details
# It basically gives the best k for each combination of top-value method and partition method by calling `suggest_best_k,ConsensusPartition-method`.
#
# 1-PAC score higher than 0.95 is treated as very stable partition and higher than 0.9 is treated as stable partition.
#
# == value
# A data frame with the best k and other statistics for each combination of methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# suggest_best_k(cola_rl)
setMethod(f = "suggest_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object, rand_index_cutoff = 0.95) {

	best_k = NULL
	stability = NULL
	mean_silhouette = NULL
	concordance = NULL
	optional_k = NULL
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			k = suggest_best_k(obj, rand_index_cutoff)
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
# Get the original matrix
#
# == param
# -object A `ConsensusPartitionList-class` object
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# get_matrix(cola_rl)
setMethod(f = "get_matrix",
	signature = "ConsensusPartitionList",
	definition = function(object) {
	object@.env$data
})

# == title
# Get the original matrix
#
# == param
# -object A `ConsensusPartition-class` object
#
# == value
# A numeric matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
# get_matrix(obj)
setMethod(f = "get_matrix",
	signature = "ConsensusPartition",
	definition = function(object) {
	object@.env$data
})

# == title
# Get annotations
#
# == param
# -object A `ConsensusPartitionList-class` object
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
# -object A `ConsensusPartition-class` object
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
# -object A `ConsensusPartitionList-class` object
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
# -object A `ConsensusPartition-class` object
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
