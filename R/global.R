
.ENV = new.env()

get_top_value_fun = function(method) {
	.ENV$ALL_TOP_VALUE_FUN[[method]]
}

register_top_value_fun = function(...) {
	lt = list(...)
	.ENV$ALL_TOP_VALUE_FUN = c(.ENV$ALL_TOP_VALUE_FUN, lt)
	.ENV$ALL_TOP_VALUE_METHOD = union(.ENV$ALL_TOP_VALUE_METHOD, names(lt))
}

ALL_TOP_VALUE_METHOD = function() {
	.ENV$ALL_TOP_VALUE_METHOD
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
		s = rowMeans(mat)
		rowMads(mat)/(rowMedians(mat) + quantile(s, 0.1))
	},
	AAC = function(mat) {
		AAC_cor(t(mat))
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

register_partition_fun = function(..., scale_method = c("normal", "positive", "none")) {
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

ALL_PARTITION_METHOD = function() {
	.ENV$ALL_PARTITION_METHOD
}

register_partition_fun(
	hclust = function(mat, k, ...) {
		hc = hclust(d = cor.dist(t(mat), abs = FALSE), ...)
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
	clara = function(mat, k, ...) {
		clara(x = t(mat), k = k, ...)
	},
	pam = function(mat, k, ...) {
		pam(t(mat), k = k, ...)
	},
	cclust = function(mat, k, ...) {
		cclust(x = t(mat), centers = k, ...)
	}
)

remove_top_value_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_TOP_VALUE_METHOD, method)
	.ENV$ALL_TOP_VALUE_FUN = .ENV$ALL_TOP_VALUE_FUN[nm_keep]
	.ENV$ALL_TOP_VALUE_METHOD = nm_keep
}

remove_partition_method = function(method) {
	nm_keep = setdiff(.ENV$ALL_PARTITION_METHOD, method)
	.ENV$ALL_PARTITION_FUN = .ENV$ALL_PARTITION_FUN[nm_keep]
	.ENV$ALL_PARTITION_METHOD = nm_keep
}