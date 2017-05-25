
strrep = function(x, times) {
	x = rep(x, times = times)
	return(paste0(x, collapse = ""))
}


relabel_class = function(class, ref) {
	class = as.character(class)
	ref = as.character(ref)

	if(length(intersect(class, ref)) == 0) {
		stop("class and ref must be from same set of labels.")
	}
	all_class = union(class, ref)
	n = length(all_class)

	m = matrix(0, nrow = n, ncol = n, dimnames = list(all_class, all_class))
	tb = table(class, ref)
	m[rownames(tb), colnames(tb)] = tb

	imap = solve_LSAP(m, maximum = TRUE)
	map = structure(rownames(m)[imap], names = rownames(m))
	
	return(map)
}

column_order_by_group = function(factor, mat) {
	do.call("c", lapply(sort(unique(factor)), function(le) {
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
# -q quantile to adjust
# 
# == details
# Vaules larger than percential ``1 - q`` are adjusted to ``1 - q`` and 
# values smaller than percential ``q`` are adjusted to ``q``.
#
# == example
# x = rnorm(10)
# x[1] = 100
# adjust_outlier(x)
adjust_outlier = function(x, q = 0.05) {
	qu = quantile(x, c(q, 1 - q))
	x[x < qu[1]] = qu[1]
	x[x > qu[2]] = qu[2]
	x
}
