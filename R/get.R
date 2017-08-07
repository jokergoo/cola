
# == title
# Get parameters
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
# -unique whether apply `base::unique` to rows
#
# == value
# A data frame of parameters corresponding to the current k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "get_param",
	signature = "ConsensusPartition",
	definition = function(object, k, unique = TRUE) {
	i = which(object@k == k)
	if(unique) {
		unique(object@object_list[[i]]$param)
	} else {
		object@object_list[[i]]$param
	}
})

# == title
# Get consensus matrix
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
#
# == value
# A consensus matrix corresponding to the current k.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
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
# -object a `ConsensusPartition-class` object
# -k number of partitions
# -each whether return the percentage membership matrix or the membership
#       in every random sampling
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
setMethod(f = "get_stat",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k) {
	m = matrix(nrow = length(object@k), ncol = length(object@object_list[[1]]$stat) - 1)
	rownames(m) = object@k
	colnames(m) = setdiff(names(object@object_list[[1]]$stat), "ecdf")
	for(i in seq_along(object@k)) {
		m[i, ] = unlist(object@object_list[[i]]$stat[colnames(m)])
	}
	return(m[as.character(k), , drop = FALSE])
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
setMethod(f = "get_class",
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
setMethod(f = "get_class",
	signature = "ConsensusPartitionList",
	definition = function(object, k) {

	if(!is.null(object@consensus_class)) {
		return(object@consensus_class)
	}

	res = object
	partition_list = NULL
	mean_cophcor = NULL
	mean_silhouette = NULL
	reference_class = NULL
	for(tm in object@top_method) {
		for(pm in object@partition_method) {
			nm = paste0(tm, ":", pm)
			obj = object@list[[nm]]
			ik = which(obj@k == k)

			membership = get_membership(obj, k)
			if(is.null(reference_class)) {
	        	reference_class = get_class(obj, k)[, "class"]
	        } else {
	        	map = relabel_class(get_class(obj, k)[, "class"], reference_class)
	        	map2 = structure(names(map), names = map)
	        	membership = membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(membership) = paste0("p", 1:k)
			}

			partition_list = c(partition_list, list(as.cl_partition(membership)))
			mean_silhouette = c(mean_silhouette, mean(get_class(obj, k)[, "silhouette"]))
		}
	}
	mean_silhouette[mean_silhouette < 0] = 0
	consensus = cl_consensus(cl_ensemble(list = partition_list), weights = mean_silhouette)
	m = cl_membership(consensus)
	class(m) = "matrix"
	colnames(m) = paste0("p", 1:k)
	class = as.vector(cl_class_ids(consensus))
	df = cbind(as.data.frame(m), class = class)
	df$entropy = apply(m, 1, entropy)

	membership_each = do.call("cbind", lapply(partition_list, function(x) {
		as.vector(cl_class_ids(x))
	}))

	consensus_mat = matrix(1, nrow = nrow(m), ncol = nrow(m))
	for(i in 1:(nrow(membership_each)-1)) {
		for(j in (i+1):nrow(membership_each)) {
			consensus_mat[i, j] = sum(membership_each[i, ] == membership_each[j, ])/ncol(membership_each)
			consensus_mat[j, i] = consensus_mat[i, j]
		}
 	}
 
	df$silhouette = silhouette(class, dist(t(consensus_mat)))[, "sil_width"]

	return(df)
})
