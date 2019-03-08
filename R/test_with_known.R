
# == title
# Test whether a list of factors are correlated
#
# == param
# -x a data frame or a vector which contains discrete or continuous variables.
#    if ``y`` is omit, pairwise testing for columns in ``x`` is performed.
# -y a data frame or a vector which contains discrete or continuous variables.
# -all_factors are all columns in ``x`` and ``y`` are enforced to be factors?
# -verbose whether to print messages.
# 
# == details
# Pairwise test is applied to every two columns in the data frames. Methods are:
# 
# - two numeric variables: correlation test by `stats::cor.test` is applied;
# - two character or factor variables: `stats::chisq.test` is applied;
# - one numeric variable and one character/factor variable: oneway ANOVA test by `stats::oneway.test` is applied.
# 
# This function can be used to test the correlation between the predicted classes and other known factors.
# 
# == return 
# A matrix of p-values. If there are NA values, basically it means there are no efficient data points to perform the test.
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
test_between_factors = function(x, y = NULL, all_factors = FALSE, verbose = FALSE) {
	
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
		for(i in seq_len(n-1)) {
			if(length(df[[i]]) < 2) {
				p.value[i, ] = NA
			} else {
				for(j in (i+1):n) {
					if(is.numeric(df[[i]]) && is.numeric(df[[j]])) {
						if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: correlation test\n")
						p.value[i, j] = cor.test(df[[i]], df[[j]])$p.value
					} else if(is.numeric(df[[i]]) && (is.character(df[[j]]) || is.factor(df[[j]]))) {
						if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: oneway ANOVA test\n")
						try({p.value[i, j] = oneway.test(df[[i]] ~ df[[j]])$p.value}, silent = TRUE)
					} else if(is.numeric(df[[j]]) && (is.character(df[[i]]) || is.factor(df[[i]]))) {
						if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: oneway ANOVA test\n")
						try({p.value[i, j] = oneway.test(df[[j]] ~ df[[i]])$p.value}, silent = TRUE)
					} else if ((is.character(df[[i]]) || is.factor(df[[i]])) && (is.character(df[[j]]) || is.factor(df[[j]]))) {
						if(verbose) qqcat("@{nm[i]} ~ @{nm[j]}: Chi-squared test\n")
						try({p.value[i, j] = chisq.test(df[[i]], df[[j]])$p.value}, silent = TRUE)
					}
					p.value[j, i] = p.value[i, j]
				}
			}
		}
		diag(p.value) = 0
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
		for(i in seq_len(n1)) {
			if(length(df1[[i]]) < 2) {
				p.value[i, ] = NA
			} else {
				for(j in seq_len(n2)) {
					if(is.numeric(df1[[i]]) && is.numeric(df2[[j]])) {
						if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: correlation test\n")
						p.value[i, j] = cor.test(df1[[i]], df2[[j]])$p.value
					} else if(is.numeric(df1[[i]]) && (is.character(df2[[j]]) || is.factor(df2[[j]]))) {
						if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: oneway ANOVA test\n")
						try({p.value[i, j] = oneway.test(df1[[i]] ~ df2[[j]])$p.value}, silent = TRUE)
					} else if((is.character(df1[[i]]) || is.factor(df1[[i]])) && is.numeric(df2[[j]])) {
						if(verbose) qqcat("@{nm2[j]} ~ @{nm1[i]}: oneway ANOVA test\n")
						try({p.value[i, j] = oneway.test(df2[[j]] ~ df1[[i]])$p.value}, silent = TRUE)
					} else if ((is.character(df1[[i]]) || is.factor(df1[[i]])) && (is.character(df2[[j]]) || is.factor(df2[[j]]))) {
						if(verbose) qqcat("@{nm1[i]} ~ @{nm2[j]}: Chi-squared test\n")
						try({p.value[i, j] = chisq.test(df1[[i]], df2[[j]])$p.value}, silent = TRUE)
					}
				}
			}
		}
	}
	return(p.value)
}

# == title
# Test correspondance between predicted classes and known factors
#
# == param
# -object a `ConsensusPartition-class` object.
# -k number of partitions. It uses all ``k`` if it is not set.
# -known a vector or a data frame with known factors. By default it is the annotation table set in `consensus_partition` or `run_all_consensus_partition_methods`.
# -silhouette_cutoff cutoff for sihouette scores. Samples with value less than this are omit.
# -verbose whether to print messages.
#
# == value
# A data frame with columns:
#
# - number of samples used to test after filtering by ``silhouette_cutoff``
# - p-values from the tests
# - number of partitions
#
# == seealso
# `test_between_factors`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# test_to_known_factors(cola_rl[1, 1], known = 1:40)
setMethod(f = "test_to_known_factors",
	signature = "ConsensusPartition",
	definition = function(object, k, known = get_anno(object), 
	silhouette_cutoff = 0.5, verbose = FALSE) {

	if(missing(k)) {
		all_k = object@k
		m = do.call("rbind", lapply(all_k, function(k) test_to_known_factors(object, k, known = known, 
			silhouette_cutoff = silhouette_cutoff, verbose = verbose)))
		return(m)
	}

	class_df = get_classes(object, k)
	l = class_df$silhouette >= silhouette_cutoff
	class = as.character(class_df$class)[l]

	if(is.null(known)) {
		stop_wrap("`known` should not be NULL.")
	}

	if(is.data.frame(known)) {
		known = known[l, , drop = FALSE]
	} else if(is.matrix(known)) {
		known = known[l, ,drop = FALSE]
	} else {
		known = known[l]
	}

	m = test_between_factors(class, known, verbose = verbose)
	rownames(m) = paste(object@top_value_method, object@partition_method, sep = ":")
	colnames(m) = paste0(colnames(m), "(p-value)")
	m = cbind(n = sum(l), m, k = k)
	return(m)
})

# == title
# Test correspondance between predicted classes and known factors
#
# == param
# -object a `ConsensusPartitionList-class` object.
# -k number of partitions. It uses all ``k`` if it is not set.
# -known a vector or a data frame with known factors. By default it is the annotation table set in `consensus_partition` or `run_all_consensus_partition_methods`.
# -silhouette_cutoff cutoff for sihouette scores. Samples with value less than this are omit.
# -verbose whether to print messages.
#
# == details
# The function basically sends each `ConsensusPartition-class` object to
# `test_to_known_factors,ConsensusPartition-method` and merges results afterwards.
#
# == value
# A data frame with columns:
#
# - number of samples used to test after filtering by ``silhouette_cutoff``
# - p-values from the tests
# - number of partitions
#
# == seealso
# `test_between_factors`, `test_to_known_factors,ConsensusPartition-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# test_to_known_factors(cola_rl, known = 1:40)
setMethod(f = "test_to_known_factors",
	signature = "ConsensusPartitionList",
	definition = function(object, k, known = get_anno(object), 
		silhouette_cutoff = 0.5, verbose = FALSE) {

	if(missing(k)) {
		m_list = lapply(object@list, function(x) {
			test_to_known_factors(x, known = known, silhouette_cutoff = silhouette_cutoff, verbose = verbose)
		})
	} else {
		m_list = lapply(object@list, function(x) {
			test_to_known_factors(x, k = k, known = known, silhouette_cutoff = silhouette_cutoff, verbose = verbose)
		})
	}
	do.call("rbind", m_list)
})
