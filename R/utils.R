
strrep = function(x, times) {
	x = rep(x, times = times)
	return(paste0(x, collapse = ""))
}


# assume class and ref have same set of labels
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

	imap = solve_LSAP(tb, maximum = TRUE)
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


#' Test whether known factors are correlated
#'
#' @param x a data frame or a vector which contains discrete or continuous variables.
#' @param y a data frame or a vector which contains discrete or continuous variables.
#' @param verbose whether print messages.
#' 
#' @details
#' Pairwise test is applied to every two columns in the data frame. Methods are:
#' 
#' - two numeric variables: correlation test by [stats::cor.test()] is applied;
#' - two character or factor variables: Chi-squared test by [stats::chisq.test()] is applied;
#' - one numeric variable and one character/factor variable: oneway ANOVA test by [stats::oneway.test()] is applied.
#' 
#' This function can be used to test the correlation between the predicted classes and other known factors.
#' 
#' @return a matrix of p-values.
#' @export
#' @importFrom stats cor.test oneway.test chisq.test
#' 
#' @examples
#' df = data.frame(v1 = rnorm(100), v2 = sample(letters[1:3], 100, replace = TRUE), 
#'     v3 = sample(LETTERS[5:6], 100, replace = TRUE))
#' test_between_known_factors(df)
#' x = runif(100)
#' test_between_known_factors(x, df)
test_between_known_factors = function(x, y = NULL, verbose = TRUE) {
	
	if(is.null(y)) {
		df = x
		nm = colnames(df)
		n = ncol(df)
		p.value = matrix(NA, nrow = n, ncol = n, dimnames = list(nm, nm))
		for(i in 1:(n-1)) {
			for(j in (i+1):n) {
				if(is.numeric(df[[i]]) && is.numeric(df[[j]])) {
					if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: correlation test\n")
					p.value[i, j] = cor.test(df[[i]], df[[j]])$p.value
				} else if(is.numeric(df[[i]]) && (is.character(df[[j]]) || is.factor(df[[j]]))) {
					if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: oneway ANOVA test\n")
					try({p.value[i, j] = oneway.test(df[[i]] ~ df[[j]])$p.value})
				} else if(is.numeric(df[[j]]) && (is.character(df[[i]]) || is.factor(df[[i]]))) {
					if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: oneway ANOVA test\n")
					try({p.value[i, j] = oneway.test(df[[j]] ~ df[[i]])$p.value})
				} else if ((is.character(df[[i]]) || is.factor(df[[i]])) && (is.character(df[[j]]) || is.factor(df[[j]]))) {
					if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: chiseq test\n")
					p.value[i, j] = chisq.test(df[[i]], df[[j]])$p.value
				}
			}
		}
	} else {
		if(is.atomic(x)) {
			df1 = data.frame(x)
			colnames(df1) = deparse(substitute(x))
		} else {
			df1 = x
		}
		if(is.atomic(y)) {
			df2 = data.frame(y)
			colnames(df2) = deparse(substitute(y))
		} else {
			df2 = y
		}
		nm1 = colnames(df1)
		nm2 = colnames(df2)
		n1 = ncol(df1)
		n2 = ncol(df2)
		p.value = matrix(NA, nrow = n1, ncol = n2, dimnames = list(nm1, nm2))
		for(i in 1:n1) {
			for(j in 1:n2) {
				if(is.numeric(df1[[i]]) && is.numeric(df2[[j]])) {
					if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: correlation test\n")
					p.value[i, j] = cor.test(df1[[i]], df2[[j]])$p.value
				} else if(is.numeric(df1[[i]]) && (is.character(df2[[j]]) || is.factor(df2[[j]]))) {
					if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: oneway ANOVA test\n")
					try({p.value[i, j] = oneway.test(df1[[i]] ~ df2[[j]])$p.value})
				} else if((is.character(df1[[i]]) || is.factor(df1[[i]])) && is.numeric(df2[[j]])) {
					if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: oneway ANOVA test\n")
					try({p.value[i, j] = oneway.test(df1[[i]] ~ df2[[j]])$p.value})
				} else if ((is.character(df1[[i]]) || is.factor(df1[[i]])) && (is.character(df2[[j]]) || is.factor(df2[[j]]))) {
					if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: chiseq test\n")
					p.value[i, j] = chisq.test(df1[[i]], df2[[j]])$p.value
				}
			}
		}
	}
	return(p.value)
}


adjust_outlier = function(x, q = 0.05) {
	qu = quantile(x, c(q, 1 - q))
	x[x < qu[1]] = qu[1]
	x[x > qu[2]] = qu[2]
	x
}
