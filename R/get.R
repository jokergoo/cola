
get_param = function(res, k) {
	i = which(res$k == k)
	res$object_list[[i]]$param
}

get_consensus = function(res, k) {
	i = which(res$k == k)
	res$object_list[[i]]$consensus
}

get_membership = function(res, k, each = FALSE) {
	i = which(res$k == k)
	if(each) {
		res$object_list[[i]]$membership_each
	} else {
		res$object_list[[i]]$membership
	}
}

get_stat = function(res, k = res$k) {
	m = matrix(nrow = length(k), ncol = length(res$object_list[[1]]$stat) - 1)
	rownames(m) = k
	colnames(m) = setdiff(names(res$object_list[[1]]$stat), "ecdf")
	for(i in seq_along(k)) {
		m[i, ] = unlist(res$object_list[[i]]$stat[colnames(m)])
	}
	m
}


#' General method for get_class
#'
#' @param x x
#' @param ... other arguments
#' 
#' @export
get_class = function(x, ...) {
	UseMethod("get_class", x)
}

#' Get class from the consensus_partition object
#'
#' @param x a `consensus_parittion` object
#' @param k number of partitions
#' @param ... other arguments.
#'
#' @return
#' A data frame with class IDs and other columns.
#' 
#' @export
get_class.consensus_partition = function(x, k, ...) {
	res = x
	i = which(res$k == k)
	
	res$object_list[[i]]$class_df
}

#' Get class from the consensus_partition_all_methods object
#'
#' @param x a `consensus_partition_all_methods` object
#' @param k number of partitions
#' @param ... other arguments
#' 
#' @details 
#' The class IDs is re-calculated by merging class IDs from all methods.
#'
#' @return
#' A data frame with class IDs and other columns.
#' 
#' @export
get_class.consensus_partition_all_methods = function(x, k, ...) {
	res = x
	partition_list = NULL
	mean_cophcor = NULL
	mean_silhouette = NULL
	reference_class = NULL
	for(tm in res$top_method) {
		for(pm in res$partition_method) {
			nm = paste0(tm, ":", pm)
			obj = res$list[[nm]]
			ik = which(obj$k == k)

			membership = obj$object_list[[ik]]$membership
			if(is.null(reference_class)) {
	        	reference_class = obj$object_list[[ik]]$classification$class
	        } else {
	        	map = relabel_class(obj$object_list[[ik]]$classification$class, reference_class, 1:k)
	        	map2 = structure(names(map), names = map)
	        	membership = membership[, as.numeric(map2[as.character(1:k)]) ]
				colnames(membership) = paste0("p", 1:k)
			}

			partition_list = c(partition_list, list(as.cl_partition(membership)))
			mean_cophcor = c(mean_cophcor, cophcor(obj$object_list[[ik]]$consensus))
			mean_silhouette = c(mean_silhouette, mean(obj$object_list[[ik]]$classification$silhouette))
		}
	}
	consensus = cl_consensus(cl_ensemble(list = partition_list), weights = 1)
	m = cl_membership(consensus)
	class(m) = "matrix"
	colnames(m) = paste0("p", 1:k)
	class = as.vector(cl_class_ids(consensus))
	cbind(as.data.frame(m), class = class)
}
