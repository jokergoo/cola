
.ENV = new.env()

get_top_value_fun = function(method) {
	.ENV$ALL_TOP_VALUE_FUN[[method]]
}

#' Register user-defined top value functions
#'
#' @param ... a named list of functions.
#' 
#' @details 
#' The user-defined function should only accept one argument which is the data
#' matrix and the scores are calculated by rows.
#' 
#' To remove a top method, use [remove_top_value_method()].
#'
#' @export
#' 
#' @examples 
#' ALL_TOP_VALUE_METHOD()
#' register_top_value_fun(mean = function(mat) rowMeans(mat),
#'                        median = function(mat) rowMedians(mat))
#' ALL_TOP_VALUE_METHOD()
register_top_value_fun = function(...) {
	lt = list(...)
	.ENV$ALL_TOP_VALUE_FUN = c(.ENV$ALL_TOP_VALUE_FUN, lt)
	.ENV$ALL_TOP_VALUE_METHOD = union(.ENV$ALL_TOP_VALUE_METHOD, names(lt))
}

#' All supported top methods
#'
#' @return a vector of supported top methods.
#' @export
#'
#' @examples
#' ALL_TOP_VALUE_METHOD()
ALL_TOP_VALUE_METHOD = function() {
	.ENV$ALL_TOP_VALUE_METHOD
}

#' @import matrixStats
register_top_value_fun(
	sd = function(mat) {
		rowSds(mat)
	},
	vc = function(mat) {
		s = rowMeans(mat)
		rowSds(mat)/(s + quantile(s, 0.1))
	},
	MAD = function(mat) {
		s = rowMeans(mat)
		rowMads(mat)/(rowMedians(mat) + quantile(s, 0.1))
	},
	AAC = function(mat) {
		AAC(t(mat))
	}
)

get_partition_fun = function(method, partition_param = list()) {
	fun = .ENV$ALL_PARTITION_FUN[[method]]
	fun2 = function(mat, k) {
		partition = do.call(fun, c(list(mat, k), partition_param))
		if(is.atomic(partition)) {
			as.cl_hard_partition(partition)
		} else {
			as.cl_hard_partition(cl_membership(partition))
		}
	}
	attr(fun2, "scale_method") = attr(fun, "scale_method")
	return(fun2)
}

#' Register user-defined partition functions
#'
#' @param ... a named list of functions
#' @param scale_method normally, data matrix are scaled by rows before sent to
#'        the partition function. The default scaling is apply by [base::scale()].
#'        However, some partition function may not accept negative values which may
#'        be produced by [base::scale()]. Here `scale_method` can be set to `rescale`
#'        which scale rows by `(x - min)/(max - min)`.
#' 
#' @details 
#' The user-defined function should only accept two arguments which are the data
#' matrix and the number of partitions. The function should return a vector of
#' partitions (or group classes).
#' 
#' The partition function is applied on rows.
#' 
#' To remove a partition method, use [remove_partition_method()].
#'
#' @export
#' 
#' @examples
#' ALL_PARTITION_METHOD()
#' register_top_value_fun(random = function(mat, k) sample(k, nrow(mat), replace = TRUE))
#' ALL_PARTITION_METHOD()
register_partition_fun = function(..., scale_method = c("standardization", "rescale", "none")) {
	
	scale_method = match.arg(scale_method)
	lt = list(...)
	if(length(scale_method) == 1) {
		scale_method = rep(scale_method, length(lt))
	}
	for(i in seq_along(lt)) {
		attr(lt[[i]], "scale_method") = scale_method[i]
	}
	.ENV$ALL_PARTITION_FUN = c(.ENV$ALL_PARTITION_FUN, lt)
	.ENV$ALL_PARTITION_METHOD = union(.ENV$ALL_PARTITION_METHOD, names(lt))
}

#' All supported partition methods
#'
#' @return a vector of supported partition methods
#' @export
#'
#' @examples
#' ALL_PARTITION_METHOD()
ALL_PARTITION_METHOD = function() {
	.ENV$ALL_PARTITION_METHOD
}

#' @import skmeans
#' @import mclust
#' @import cclust
#' @import cluster
#' @import kohonen
register_partition_fun(
	hclust = function(mat, k, ...) {
		hc = hclust(d = dist(t(mat)), ...)
		cutree(hc, k)
	},
	kmeans = function(mat, k, ...) {
		kmeans(t(mat), centers = k, ...)
	},
	# NMF = function(mat, k, ...) {
	# 	as.numeric(predict(nmf(x = mat, rank = k, ...)))
	# },
	skmeans = function(mat, k, ...) {
		skmeans(x = t(mat), k = k, ...)
	},
	Mclust = function(mat, k, ...) {
		Mclust(data = t(mat), G = k, ...)
	}, 
	# clara = function(mat, k, ...) {
	# 	clara(x = t(mat), k = k, ...)
	# },
	pam = function(mat, k, ...) {
		pam(t(mat), k = k, ...)
	},
	cclust = function(mat, k, ...) {
		cclust(x = t(mat), centers = k, ...)
	},
	som = function(mat, k, ...) {
		som(mat, grid = somgrid(1, k), ...)$unit.classif
	}
)

#' Remove top methods
#'
#' @param method name of the top methods to be removed.
#'
#' @export
remove_top_value_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_TOP_VALUE_METHOD, method)
	.ENV$ALL_TOP_VALUE_FUN = .ENV$ALL_TOP_VALUE_FUN[nm_keep]
	.ENV$ALL_TOP_VALUE_METHOD = nm_keep
}


#' Remove partition methods
#'
#' @param method name of the partition methods to be removed.
#'
#' @export
remove_partition_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_PARTITION_METHOD, method)
	.ENV$ALL_PARTITION_FUN = .ENV$ALL_PARTITION_FUN[nm_keep]
	.ENV$ALL_PARTITION_METHOD = nm_keep
}