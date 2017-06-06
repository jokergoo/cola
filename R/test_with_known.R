
# == title
# Test whether known factors are correlated
#
# == param
# -x a data frame or a vector which contains discrete or continuous variables.
# -y a data frame or a vector which contains discrete or continuous variables.
# -all_factors are all columns in ``x`` and ``y`` are enforced to be factors.
# -verbose whether print messages.
# 
# == details
# Pairwise test is applied to every two columns in the data frame. Methods are:
# 
# - two numeric variables: correlation test by `stats::cor.test` is applied;
# - two character or factor variables: Chi-squared test by `stats::fisher.test` is applied;
# - one numeric variable and one character/factor variable: oneway ANOVA test by `stats::fisher.test` is applied.
# 
# This function can be used to test the correlation between the predicted classes and other known factors.
# 
# == return 
# A matrix of p-values.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == examples
# df = data.frame(v1 = rnorm(100), v2 = sample(letters[1:3], 100, replace = TRUE), 
#     v3 = sample(LETTERS[5:6], 100, replace = TRUE))
# test_between_factors(df)
# x = runif(100)
# test_between_factors(x, df)
test_between_factors = function(x, y = NULL, all_factors = FALSE, verbose = TRUE) {
	
	if(is.null(y)) {
		df = x
		if(is.matrix(df)) {
			df = as.data.frame(df)
		}
		if(all_factors) {
			for(i in seq_len(ncol(df))) {
				df[[i]] = factor(df[[i]])
			}
		}
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
					if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: Fisher's exact test\n")
					try({p.value[i, j] = fisher.test(df[[i]], df[[j]], alternative = "greater")$p.value})
				}
			}
		}
	} else {
		if(is.matrix(x)) {
			x = as.data.frame(x)
		}
		if(is.matrix(y)) {
			y = as.data.frame(y)
		}
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
		if(all_factors) {
			for(i in seq_len(ncol(df1))) {
				df1[[i]] = factor(df1[[i]])
			}
			for(i in seq_len(ncol(df2))) {
				df2[[i]] = factor(df2[[i]])
			}
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
					try({p.value[i, j] = oneway.test(df2[[j]] ~ df1[[i]])$p.value})
				} else if ((is.character(df1[[i]]) || is.factor(df1[[i]])) && (is.character(df2[[j]]) || is.factor(df2[[j]]))) {
					if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: Fisher's exact test\n")
					try({p.value[i, j] = fisher.test(df1[[i]], df2[[j]], alternative = "greater")$p.value})
				}
			}
		}
	}
	return(p.value)
}

# == title
# Test correspondance between predicted and known classes
#
# == param
# -object a `ConsensusPartition-class` object
# -k number of partitions
# -known a vector or a data frame with known factors
#
# == value
# A matrix of p-values 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "test_to_known_factors",
	signature = "ConsensusPartition",
	definition = function(object, k, known = object@known_anno) {

	class = as.character(get_class(object, k)$class)

	m = test_between_factors(class, known, verbose = FALSE)
	rownames(m) = paste(object@top_method, object@partition_method, sep = ":")
	return(m)
})

# == title
# Test correspondance between predicted and known classes
#
# == param
# -object a `ConsensusPartitionList-class` object
# -k number of partitions
# -known a vector or a data frame with known factors
#
# == value
# A matrix of p-values 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "test_to_known_factors",
	signature = "ConsensusPartitionList",
	definition = function(object, k, known = object@list[[1]]@known_anno) {

	class_df = lapply(object@list, function(x) as.character(get_class(x, k)$class))
	class_df = as.data.frame(class_df)

	test_between_factors(class_df, known, verbose = FALSE)
})
