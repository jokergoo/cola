
strrep = function(x, times) {
	x = rep(x, times = times)
	return(paste0(x, collapse = ""))
}


# change the label of `class` to let `class` fits `ref` more, maximize sum(class ==ref)
# class = c(rep("a", 10), rep("b", 3))
# ref = c(rep("b", 4), rep("a", 9))
relabel_class = function(class, ref) {
	class = as.character(class)
	ref = as.character(ref)

	# if(length(intersect(class, ref)) == 0) {
	# 	stop("class and ref must be from same set of labels.")
	# }
	all_class = union(class, ref)
	n = length(all_class)

	m = matrix(0, nrow = n, ncol = n, dimnames = list(all_class, all_class))
	tb = table(class, ref)
	m[rownames(tb), colnames(tb)] = tb

	imap = solve_LSAP(m, maximum = TRUE)
	map = structure(rownames(m)[imap], names = rownames(m))

	df = data.frame(class = class, adjusted = map[class], ref = ref)
	attr(map, "df") = df
	return(map)
}

# columns only in one same level are clustered
column_order_by_group = function(factor, mat) {
	if(!is.factor(factor)) {
		factor = factor(factor, levels = unique(factor))
	}
	do.call("c", lapply(sort(levels(factor)), function(le) {
		m = mat[, factor == le, drop = FALSE]
		if(ncol(m) == 1) {
			which(factor == le)
		} else {
			hc1 = hclust(dist(t(m)))
			oe = try({ 
				hc1 = as.hclust(reorder(as.dendrogram(hc1), colSums(m)))
			}, silent = TRUE)
			col_od1 = hc1$order
			which(factor == le)[col_od1]
		}
	}))
}


# == title
# Adjust outliers
#
# == param
# -x a numeric vector
# -q percential to adjust
# 
# == details
# Vaules larger than percential ``1 - q`` are adjusted to ``1 - q`` and 
# values smaller than percential ``q`` are adjusted to ``q``.
#
# == value
# a numeric vector with same length as the original one.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# x = rnorm(10)
# x[1] = 100
# adjust_outlier(x)
adjust_outlier = function(x, q = 0.05) {
	qu = quantile(x, c(q, 1 - q), na.rm = TRUE)
	x[x < qu[1]] = qu[1]
	x[x > qu[2]] = qu[2]
	x
}

# == title
# Remove rows with low variance
#
# == param
# -m a numeric matrix
# -sd_quantile cutoff the quantile of standard variation
#
# == details
# The function first uses `adjust_outlier` to adjust outliers and 
# then removes rows with low standard variation.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
adjust_matrix = function(m, sd_quantile = 0.05) {
	if(is.data.frame(m)) {
		m = as.matrix(m)
	}
	m = t(apply(m, 1, adjust_outlier))
	row_sd = rowSds(m, na.rm = TRUE)
	l = abs(row_sd) < 1e-6
	m2 = m[!l, , drop = FALSE]
	row_sd = row_sd[!l]
	qqcat("@{sum(l)} rows have been removed for zero variance (sd < 1e-6).\n")

	qa = quantile(row_sd, sd_quantile, na.rm = TRUE)
	l = row_sd >= qa
	qqcat("@{sum(!l)} rows have been removed for too low variance (sd < @{sd_quantile} quantile)\n")
	m2[l, , drop = FALSE]
}
