
# == title
# Get parameters
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
# -unique whether apply `base::unique` to rows of the returned data frame.
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
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
#
# == details
# For row i and column j in the consensus matrix, the value of corresponding x_ij
# is the probability of sample i and sample j being in a same group from all partitions.
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
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
# -each whether return the percentage membership matrix which is summarized from all partitions
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
# -object a `ConsensusPartitionList-class` object.
# -k number of partitions.
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
# -object a `ConsensusPartition-class` object.
# -k number of partitions. The value can be a vector.
#
# == details
# The statistics are:
#
# -cophcor cophenetic correlation coefficient. It measures if hierarchical clustering is applied
#          on the consensus matrix, how good it correlates to the consensus matrix itself.
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
	definition = function(object, k = object@k) {

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
# -object a `ConsensusPartitionList-class` object.
# -k number of partitions. The value can only be a single value.
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
	definition = function(object, k) {
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
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
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
# -object a `ConsensusPartitionList-class` object.
# -k number of partitions.
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
			rl@list[[j]]@object_list[[i]]$stat[["1-PAC"]] = 1 - PAC_origin(consensus_mat[l, l, drop = FALSE])
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
# -object a `ConsensusPartition-class` object.
# -rand_index_cutoff the Rand index compared to previous k is larger than this value, it is filtered out.
#
# == details
# The best k is voted from 1) the k with the maximal cophcor value, 2) the k with the minimal PAC value,
# 3) the k with the maximal mean silhouette value and 4) the k with the maximal concordance value.
#
# There are scenarios that a better partition with k groups than k - 1 groups (e.g. for the sense of better sihouette score) 
# is only because of one tiny group of samples are separated and it is better to still put them back to the original group
# to improve the robustness of the subgrouping. For this, users can set the cutoff of Rand index by ``rand_index_cutoff`` to
# get rid of or reduce the effect of such cirsumstances.
#
# Honestly, it is hard or maybe impossible to say which k is the best one. `suggest_best_k` function only gives suggestion of selecting
# a reasonable k. Users still need to look at the plots (e.g. by `select_partition_number` or `consensus_heatmap` functions), or even
# by checking whether the subgrouping gives a reasonable signatures by `get_signatures`, to pick a reasonable k that best explains their study. 
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
	if(sum(l)) {
		return(max(as.numeric(rownames(stat)[l])))
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
# -object a `ConsensusPartitionList-class` object.
# -rand_index_cutoff the Rand index compared to previous k is larger than this, it is filtered out.
#
# == details
# It basically gives the best k for each combination of top-value method and partition method by calling `suggest_best_k,ConsensusPartition-method`.
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
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			best_k[nm] = suggest_best_k(obj, rand_index_cutoff)
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
	return(tb[order(!is.na(tb[, "best_k"]) + 0, tb[, "1-PAC"], decreasing = TRUE), , drop = FALSE])
})

# == title
# Get the original matrix
#
# == param
# -object a `ConsensusPartitionList-class` object
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
# -object a `ConsensusPartition-class` object
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
# -object a `ConsensusPartitionList-class` object
#
# == value
# A data frame if ``anno`` was specified in `run_all_consensus_partition_methods`, or ``NULL``.
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
# -object a `ConsensusPartition-class` object
#
# == value
# A data frame if ``anno`` was specified in `run_all_consensus_partition_methods` or `consensus_partition`, or ``NULL``.
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
# -object a `ConsensusPartitionList-class` object
#
# == value
# A list of color vectors or ``NULL``.
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
# -object a `ConsensusPartition-class` object
#
# == value
# A list of color vectors or ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_anno_col",
	signature = "ConsensusPartition",
	definition = function(object) {
	object@anno_col
})
