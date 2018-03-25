
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
# The registered top method will be used in `run_all_consensus_partition_methods`.
# 
# To remove a top method, use `remove_top_value_method`.
#
# == return
# No value is returned.
#
# == seealso
# `all_top_value_methods` lists all registered top methods.
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
	lt1 = lt[intersect(names(lt), .ENV$ALL_TOP_VALUE_METHODS)]
	lt2 = lt[setdiff(names(lt), .ENV$ALL_TOP_VALUE_METHODS)]
	if(length(lt1)) .ENV$ALL_TOP_VALUE_FUN[names(lt1)] = lt1
	if(length(lt2)) .ENV$ALL_TOP_VALUE_FUN = c(.ENV$ALL_TOP_VALUE_FUN, lt2)
	.ENV$ALL_TOP_VALUE_METHODS = names(.ENV$ALL_TOP_VALUE_FUN)
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
	cv = function(mat) {
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
		stop(qq("partition method '@{method}' has not been defined yet."))
	}
	fun = .ENV$ALL_PARTITION_FUN[[method]]
	fun2 = function(mat, k) {
		partition = do.call(fun, c(list(mat, k), partition_param))
		if(is.atomic(partition)) {
			x = as.cl_hard_partition(partition)
		} else {
			x = as.cl_hard_partition(cl_membership(partition))
		}

		nc = n_of_classes(x)
		# if(nc != k) {
		# 	cat(red(qq("!!! @{method}: number of classes (@{nc}) in the partition is not same as @{k} !!!\n")))
		# }

		return(x)
	}
	attr(fun2, "scale_method") = attr(fun, "scale_method")
	attr(fun2, "execution_scale") = attr(fun, "execution_scale")
	return(fun2)
}

# == title
# Register user-defined partition functions
#
# == param
# -... a named list of functions
# -scale_method normally, data matrix are scaled by rows before sent to
#        the partition function. The default scaling is apply by `base::scale`.
#        However, some partition function may not accept negative values which 
#        are produced by `base::scale`. Here ``scale_method`` can be set to ``rescale``
#        which scales rows by ``(x - min)/(max - min)``. Note here ``scale_method`` only means
#        the method to scale rows. When ``scale_rows`` is set to ``FALSE`` in `consensus_partition`
#        or `run_all_consensus_partition_methods`, there wil be no row scaling before doing partition.
#        The value for ``scale_method`` can be a vector if user specifies more than one partition function.
# 
# == details 
# The user-defined function should only accept three arguments. The first two arguments are the data
# matrix and the number of partitions. The third argument should always be `...` so that parameters
# for the partition function can be passed by ``partition_param`` from `consensus_partition` or `run_all_consensus_partition_methods`.
#
# The function should return a vector of partitions (or class labels).
# 
# The partition function should be applied on columns.
#
# The registered partition method will be used in `run_all_consensus_partition_methods`.
# 
# To remove a partition method, use `remove_partition_method`.
#
# == value
# No value is returned.
#
# == seealso
# `all_partition_methods` lists all registered partition methods.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# all_partition_methods()
# register_partition_fun(random = function(mat, k) sample(k, ncol(mat), replace = TRUE))
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

	m = matrix(rnorm(100), 10)
	for(i in seq_along(lt)) {
		t = microbenchmark(lt[[i]](m, 2), times = 10)
		attr(lt[[i]], "execution_scale") = mean(t$time)
	}

	lt1 = lt[intersect(names(lt), .ENV$ALL_PARTITION_METHODS)]
	lt2 = lt[setdiff(names(lt), .ENV$ALL_PARTITION_METHODS)]
	if(length(lt1)) .ENV$ALL_PARTITION_FUN[names(lt1)] = lt1
	if(length(lt2)) .ENV$ALL_PARTITION_FUN = c(.ENV$ALL_PARTITION_FUN, lt2)
	.ENV$ALL_PARTITION_METHODS = names(.ENV$ALL_PARTITION_FUN)
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
	skmeans = function(mat, k, ...) {
		skmeans(x = t(mat), k = k, ...)
	},
	# hddc = function(mat, k, ...) {
	# 	i = 0
	# 	while(i < 50) {
	# 		suppressWarnings(cl <- hddc(data = t(mat), K = k, show = FALSE, ...)$class)
	# 		if(length(cl) != 0) {
	# 			return(cl)
	# 		}
	# 		i = i + 1
	# 	}
	# 	return(rep(1, ncol(mat)))
	# }, 
	pam = function(mat, k, ...) {
		pam(t(mat), k = k, ...)
	},
	# cclust = function(mat, k, ...) {
	# 	cclust(x = t(mat), centers = k, ...)
	# },
	mclust = function(mat, k, ...) {
		pca = prcomp(t(mat))
		Mclust(pca$x[, 1:3], G = k, verbose = FALSE, ...)$classification
	},
	som = function(mat, k, ...) {
		kr = floor(sqrt(ncol(mat)))
		somfit = som(t(mat), grid = somgrid(kr, kr, "hexagonal"), ...)
		m = somfit$codes[[1]]
		m = m[seq_len(nrow(m)) %in% somfit$unit.classif, ]
		cl = cutree(hclust(dist(m)), k)
		group = numeric(ncol(mat))
		for(cl_unique in unique(cl)) {
			ind = as.numeric(gsub("V", "", names(cl)[which(cl == cl_unique)]))
			l = somfit$unit.classif %in% ind
			group[l] = cl_unique
		}
		group
	}
)

register_partition_fun(
	NMF = function(mat, k, ...) {
		fit = nnmf(A = mat, k = k, verbose = FALSE, ...)
		apply(fit$H, 2, which.max)
	}, scale_method = "rescale"
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


brewer_pal_set1_col = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
names(brewer_pal_set1_col) = 1:17
brewer_pal_set2_col = c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
names(brewer_pal_set2_col) = 1:16
