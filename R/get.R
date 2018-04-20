
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
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
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
# == value
# A consensus matrix corresponding to the current k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
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
# -each whether return the percentage membership matrix which is summarized from all random samplings
#       or the individual membership in every random sampling.
#
# == details
# If ``each == TRUE``, the value in the membership matrix is the probability
# to be in one subgroup, while if ``each == FALSE``, the membership matrix contains the 
# class labels for the partitions in all random samplings.
#
# == value
# A membership matrix where rows correspond to the samples in the original matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
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

setMethod(f = "get_membership",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	df = object@consensus_class[[as.character(k)]]
	df[, grep("^p\\d+$", colnames(df)), drop = FALSE]
})



# == title
# Get statistics
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions. The value can be a vector.
#
# == details
# The statistics are:
#
# -cophcor cophenetic correlation coefficient.
# -PAC proportion of ambiguous clustering, calculated by `PAC`.
# -mean_silhouette the mean silhouette score.
#
# == value
# A matrix of partition statistics for all k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
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
# -object a `ConsensusPartitionList-class` object
# -k number of partitions.
#
# == value
# A matrix of partition statistics for a selected k. Rows in the 
# matrix correspond to all combinations of top methods and partition methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# get_stat(rl, k = 2)
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
# -object a `ConsensusPartition-class` object
# -k number of partitions
#
# == return
# A data frame with class IDs and other columns which are entropy of the membership matrix
# and the silhouette scores which measure the stability of a sample to stay in its subgroup.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
# get_classes(obj, k = 2)
setMethod(f = "get_classes",
	signature = "ConsensusPartition",
	definition = function(object, k) {
	i = which(object@k == k)
	object@object_list[[i]]$class_df
})

# == title
# Get class from the ConsensusPartitionList object
#
# == param
# -object a `ConsensusPartitionList-class` object
# -k number of partitions
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
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# get_classes(rl, k = 2)
setMethod(f = "get_classes",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {
	df = object@consensus_class[[as.character(k)]]
	df[, !grepl("^p\\d+$", colnames(df)), drop = FALSE]
})

# == title
# Guess the best number of partitions
#
# == param
# -object a `ConsensusPartition-class` object
#
# == details
# It looks for the best k with highest cophenetic correlation coefficient
# or lowest PAC score or highest mean silhouette value.
#
# == value
# The best k
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
# guess_best_k(obj)
setMethod(f = "guess_best_k",
	signature = "ConsensusPartition",
	definition = function(object) {

	stat = get_stat(object)
	stat = stat[stat[, "separation_rate"] > 0.2, , drop = FALSE]
	dec = c(which.max(stat[, "cophcor"]), 
		    which.min(stat[, "iPAC"]),
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
# -object a `ConsensusPartitionList-class` object
#
# == value
# A data frame with best k for each combination of methods
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# guess_best_k(rl)
setMethod(f = "guess_best_k",
	signature = "ConsensusPartitionList",
	definition = function(object) {

	best_k = NULL
	cophcor = NULL
	iPAC = NULL
	mean_silhouette = NULL
	concordance = NULL
	for(tm in object@top_value_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			best_k[nm] = guess_best_k(obj)
			stat = get_stat(obj, k = best_k[nm])
			cophcor[nm] = stat[1, "cophcor"]
			iPAC[nm] = stat[1, "iPAC"]
			mean_silhouette[nm] = stat[1, "mean_silhouette"]
			concordance[nm] = stat[1, "concordance"]
		}
	}
	tb = data.frame(best_k = best_k,
		cophcor = cophcor,
		iPAC = iPAC,
		mean_silhouette = mean_silhouette,
		concordance = concordance)

	rntb = rownames(tb)
	l = tb$concordance > 0.9
	tb = cbind(tb, ifelse(l, ifelse(tb$concordance > 0.95, "**", "*"), ""), stringsAsFactors = FALSE)
	colnames(tb)[ncol(tb)] = ""
	return(tb)
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
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# get_matrix(rl)
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
# rl = readRDS(system.file("extdata/example.rds", package = "cola"))
# obj = rl["sd", "kmeans"]
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
