
.ENV = new.env()

get_top_value_fun = function(method) {
	if(!method %in% .ENV$ALL_TOP_VALUE_METHODS) {
		stop(qq("top value method @{method} has not been defined yet."))
	}
	.ENV$ALL_TOP_VALUE_FUN[[method]]
}

# == title
# Register user-defined top value functions
#
# == param
# -... a named list of functions.
# 
# == details 
# The user-defined function should only accept one argument which is the data
# matrix and the scores are calculated by rows.
# 
# To remove a top method, use `remove_top_value_method`.
#
# == return
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == examples 
# all_top_value_methods()
# register_top_value_fun(AAC_spearman = function(mat) AAC(t(mat), cor_method = "spearman"),
#                        AAC_multicore = function(mat) AAC(t(mat), mc.cores = 2))
# all_top_value_methods()
# remove_top_value_method(c("AAC_spearman", "AAC_multicore"))
register_top_value_fun = function(...) {
	lt = list(...)
	.ENV$ALL_TOP_VALUE_FUN = c(.ENV$ALL_TOP_VALUE_FUN, lt)
	.ENV$ALL_TOP_VALUE_METHODS = union(.ENV$ALL_TOP_VALUE_METHODS, names(lt))
}

# == title
# All supported top methods
#
# == return 
# A vector of supported top methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# all_top_value_methods()
all_top_value_methods = function() {
	.ENV$ALL_TOP_VALUE_METHODS
}

register_top_value_fun(
	sd = function(mat) {
		rowSds(mat)
	},
	vc = function(mat) {
		s = rowMeans(mat)
		rowSds(mat)/(s + quantile(s, 0.1))
	},
	MAD = function(mat) {
		rowMads(mat)
	},
	AAC = function(mat) {
		AAC(t(mat))
	}
)

get_partition_fun = function(method, partition_param = list()) {
	if(!method %in% .ENV$ALL_PARTITION_METHODS) {
		stop(qq("partition method @{method} has not been defined yet."))
	}
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

# == title
# Register user-defined partition functions
#
# == param
# -... a named list of functions
# -scale_method normally, data matrix are scaled by rows before sent to
#        the partition function. The default scaling is apply by `base::scale`.
#        However, some partition function may not accept negative values which may
#        be produced by `base::scale`. Here ``scale_method`` can be set to ``rescale``
#        which scale rows by ``(x - min)/(max - min)``.
# 
# == details 
# The user-defined function should only accept three arguments. The first two arguments are the data
# matrix and the number of partitions. The third argument should always be `...` so that parameters
# for the partition function can be passed by ``partition_param`` from `consensus_partition` or `run_all_consensus_partition_methods`.
# The function should return a vector of partitions (or group classes).
# 
# The partition function is applied on rows.
# 
# To remove a partition method, use `remove_partition_method`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# all_partition_methods()
# register_partition_fun(random = function(mat, k) sample(k, nrow(mat), replace = TRUE))
# all_partition_methods()
# remove_partition_method("random")
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
	.ENV$ALL_PARTITION_METHODS = union(.ENV$ALL_PARTITION_METHODS, names(lt))
}

# == title
# All supported partition methods
#
# == return 
# A vector of supported partition methods
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# all_partition_methods()
all_partition_methods = function() {
	.ENV$ALL_PARTITION_METHODS
}

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
	hddc = function(mat, k, ...) {
		i = 0
		while(i < 50) {
			suppressWarnings(cl <- hddc(data = t(mat), K = k, show = FALSE, ...)$class)
			if(length(cl) != 0) {
				return(cl)
			}
			i = i + 1
		}
		return(rep(1, ncol(mat)))
	}, 
	pam = function(mat, k, ...) {
		pam(t(mat), k = k, ...)
	},
	cclust = function(mat, k, ...) {
		cclust(x = t(mat), centers = k, ...)
	},
	som = function(mat, k, ...) {
		i = 0
		while(i < 50) {
			cl = som(t(mat), grid = somgrid(1, k), ...)$unit.classif
			if(length(unique(cl)) > 1) {
				return(cl)
			}
			i = i + 1
		}
		return(rep(1, ncol(mat)))
	}
)

# == title
# Remove top methods
#
# == param
# -method name of the top methods to be removed.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
remove_top_value_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_TOP_VALUE_METHODS, method)
	.ENV$ALL_TOP_VALUE_FUN = .ENV$ALL_TOP_VALUE_FUN[nm_keep]
	.ENV$ALL_TOP_VALUE_METHODS = nm_keep
}


# == title
# Remove partition methods
#
# == param
# -method name of the partition methods to be removed.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
remove_partition_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_PARTITION_METHODS, method)
	.ENV$ALL_PARTITION_FUN = .ENV$ALL_PARTITION_FUN[nm_keep]
	.ENV$ALL_PARTITION_METHODS = nm_keep
}


brewer_pal_set1_col =  brewer.pal(9, "Set1")
brewer_pal_set2_col =  brewer.pal(8, "Set2")
