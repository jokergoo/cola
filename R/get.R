
# == title
# Get parameters
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
# -unique whether apply `base::unique` to rows
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
# -each return the percentage membership matrix or the membership
#       in every random round
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
# Get parameters
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions. The value can be a vector.
#
setMethod(f = "get_stat",
	signature = "ConsensusPartition",
	definition = function(object, k = object@k) {
	m = matrix(nrow = length(k), ncol = length(object@object_list[[1]]$stat) - 1)
	rownames(m) = k
	colnames(m) = setdiff(names(object@object_list[[1]]$stat), "ecdf")
	for(i in seq_along(k)) {
		m[i, ] = unlist(object@object_list[[i]]$stat[colnames(m)])
	}
	return(m)
})


# == title
# Get parameters
#
# == param
# -object a `ConsensusPartitionList-class` object
# -k number of partitions.
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
# Get class from the consensus_partition object
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
#
# == return
# A data frame with class IDs and other columns.
#
setMethod(f = "get_class",
	signature = "ConsensusPartition",
	definition = function(object, k) {
	i = which(object@k == k)
	object@object_list[[i]]$class_df
})

# == title
# Get class from the consensus_partition_all_methods object
#
# == param
# -object a `ConsensusPartitionList-class` object
# -k number of partitions
# -... other arguments
# 
# == details 
# The class IDs is re-calculated by merging class IDs from all methods.
#
# == return
# A data frame with class IDs and other columns.
#
setMethod(f = "get_class",
	signature = "ConsensusPartitionList",
	definition = function(object, k, ...) {
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
	cbind(as.data.frame(m), class = class)
})
