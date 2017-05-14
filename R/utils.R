
strrep = function(x, times) {
	x = rep(x, times = times)
	return(paste0(x, collapse = ""))
}

# == title
# Adjust class IDs
#
# == param
# -class a vector of class ID
# -reference a vector of reference class ID
# -full_set full set of class IDs
#
relabel_class = function(class, reference, full_set = unique(c(reference, class))) {
	class = as.character(class)
	reference = as.character(reference)
	
	tb = tapply(reference, class, table, simplify = FALSE)
	map = NULL
	while(length(tb)) {
		max_p = sapply(tb, function(x) max(x))
		which_max = sapply(tb, function(x) names(which.max(x)))
		max_i = which.max(max_p)
		map = c(map, structure(which_max[max_i], names = names(tb)[max_i]))
		tb = tb[-max_i]
		tb = lapply(tb, function(x) x[!(names(x) %in% map)])
		tb = tb[sapply(tb, length) > 0]
	}
	sdf = setdiff(full_set, map)
	if(length(sdf)) {
		unmapped = structure(sdf, names = setdiff(full_set, names(map)))
		map = c(map, unmapped)
	}
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
