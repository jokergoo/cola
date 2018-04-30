
# == title
# Get parameters
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
# -unique whether apply `base::unique` to rows of the returned data frame.
#
# == value
# A data frame of parameters corresponding to the current k.
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
		unique(do.call("rbind", lapply(ind, function(i) object@object_list[[i]]$param)))
	} else {
		do.call("rbind", lapply(ind, function(i) object@object_list[[i]]$param))
	}
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
# is the probability of sample i and sample j being in a same subgroup from the repetitive 
# partitionings.
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
	i = which(object@k == k)
	object@object_list[[i]]$consensus
})

# == title
# Get membership matrix
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
# -each whether return the percentage membership matrix which is summarized from all repetitive partitionings
#       or the individual membership in every random partitioning.
#
# == details
# If ``each == TRUE``, the value in the membership matrix is the probability
# to be in one subgroup, while if ``each == FALSE``, the membership matrix contains the 
# class labels for the partitions in all repetitive partitionings with randomly sampling subset
# of rows in the matrix.
#
# The percent membership matrix is calculated by `clue::cl_consensus`.
#
# == value
# A membership matrix where rows correspond to the samples in the original matrix.
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
# The membership matrix (the probability of each sample to be in a subgroup) is inferred
# from the membership matrices of all combinations of methods, weighted by the mean silhouette score of the partitions
# for each method. So methods which give instable partitions have lower weights 
# when summarizing membership matrix from all methods.
# 
# == value
# A membership matrix where rows correspond to the samples in the original matrix.
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
	lt = object@consensus_class[[as.character(k)]]
	lt$membership
})



# == title
# Get statistics
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
# -mean_silhouette the mean silhouette score.
# -concordance the mean probability that each partition fits the consensus partition, calculated by `concordance`.
# -area_increased the increased area under ecdf to the previous k
# -Rand the Rand index which is the percent of pairs of samples that are both in a same cluster or both are not 
#       in a same cluster in the partition of ``k`` and ``k-1``.
# -Jaccard the ratio of pairs of samples are both in a same cluster in the partition of ``k`` and ``k-1`` and the pairs
#          of samples are both in a same cluster in the partition ``k`` or ``k-1``.
#
# == value
# A matrix of partition statistics for all k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
# get_stat(obj)
# get_stat(obj, k = 2)
setMethod(f = "get_stat",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k) {
	m = matrix(nrow = length(object@k), ncol = length(object@object_list[[1]]$stat) - 1)
	rownames(m) = object@k
	colnames(m) = setdiff(names(object@object_list[[1]]$stat), "ecdf")
	for(i in seq_along(object@k)) {
		m[i, ] = unlist(object@object_list[[i]]$stat[colnames(m)])
	}
	m = m[as.character(k), , drop = FALSE]
	# m = cbind(m, separation_rate = sapply(k, function(i) separation_rate(object, i)))
	return(m)
})


# == title
# Get statistics
#
# == param
# -object a `ConsensusPartitionList-class` object.
# -k number of partitions.
#
# == value
# A matrix of partition statistics for a selected k. Rows in the 
# matrix correspond to all combinations of top value methods and partition methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# get_stat(cola_rl, k = 2)
setMethod(f = "get_stat",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	m = matrix(nrow = length(object@list), ncol = length(object@list[[1]]@object_list[[1]]$stat) - 1)
	rownames(m) = names(object@list)
	colnames(m) = setdiff(names(object@list[[1]]@object_list[[1]]$stat), "ecdf")
	for(i in seq_along(object@list)) {
		ik = which(object@list[[i]]@k == k)
		m[i, ] = unlist(object@list[[i]]@object_list[[ik]]$stat[colnames(m)])
	}
	return(m)
})


# == title
# Get class from the ConsensusPartition object
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions.
#
# == return
# A data frame with class IDs and other columns which are entropy of the membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its group.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
# get_classes(obj, k = 2)
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
# Get class from the ConsensusPartitionList object
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
# A data frame with class IDs, membership, entropy and silhouette scores.
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
	lt = object@consensus_class[[as.character(k)]]
	lt$class_df
})

# == title
# Guess the best number of partitions
#
# == param
# -object a `ConsensusPartition-class` object.
# -rand_index_cutoff the Rand index compared to previous k is larger than this, it is filtered out.
#
# == details
# The best k is voted from 1) which k has the maximum cophcor value, 2) which k has the minimal PAC value,
# 3) which k has the maximum mean silhouette value and 4) which k has the maximum concordance value.
#
# == value
# The best k
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# obj = cola_rl["sd", "kmeans"]
# guess_best_k(obj)
setMethod(f = "guess_best_k",
	signature = "ConsensusPartition",
	definition = function(object, rand_index_cutoff = 0.9) {

	stat = get_stat(object)
	stat = stat[stat[, "Rand"] < rand_index_cutoff, , drop = FALSE]
	dec = c(which.max(stat[, "cophcor"]), 
		    which.min(stat[, "PAC"]),
		    which.max(stat[, "mean_silhouette"]),
		    which.max(stat[, "concordance"]))
	
	tb = table(dec)
	x = rownames(stat)[as.numeric(names(tb[which.max(tb)[1]]))]
	
	as.numeric(x)
})


# == title
# Get the best number of partitions
#
# == param
# -object a `ConsensusPartitionList-class` object.
# -rand_index_cutoff the Rand index compared to previous k is larger than this, it is filtered out.
#
# == details
# It basically gives best k for each combination of top value method and partition method.
#
# == value
# A data frame with best k for each combination of methods
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# guess_best_k(cola_rl)
setMethod(f = "guess_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object, rand_index_cutoff = 0.9) {

	best_k = NULL
	cophcor = NULL
	PAC = NULL
	mean_silhouette = NULL
	concordance = NULL
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			best_k[nm] = guess_best_k(obj, rand_index_cutoff)
			stat = get_stat(obj, k = best_k[nm])
			cophcor[nm] = stat[1, "cophcor"]
			PAC[nm] = stat[1, "PAC"]
			mean_silhouette[nm] = stat[1, "mean_silhouette"]
			concordance[nm] = stat[1, "concordance"]
		}
	}
	tb = data.frame(best_k = best_k,
		cophcor = cophcor,
		PAC = PAC,
		mean_silhouette = mean_silhouette,
		concordance = concordance)

	rntb = rownames(tb)
	l = tb$concordance > 0.9
	tb = cbind(tb, ifelse(l, ifelse(tb$concordance > 0.95, "**", "*"), ""), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""
	return(tb[order(tb[, "concordance"], decreasing = TRUE), , drop = FALSE])
})

# == title
# Get original matrix
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
# Get original matrix
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
