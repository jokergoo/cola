

# == title
# AAC score
#
# == param
# -mat a numeric matrix. AAC score is calculated by rows.
# -cor_method pass to `stats::cor`.
# -min_cor minimal absolute correlation.
# -max_cor maximal absolute correlation.
# -mc.cores number of cores.
# -n_sampling when the number of columns are too big, to get the curmulative
#           distribution, actually we don't need to use all the rows, e.g. 1000
#           rows can already give a farely nice estimation for the distribution.
# -q_sd percential of the sd for the rows to ignore.
# 
# == details 
# For a given row in a matrix, the AAC score is the area above the curve of the curmulative density
# distribution of the absolute correlation to all other rows. Formally, if ``F_i(x)`` is the 
# CDF of the absolute correlation for row i, ``AAC_i = 1 - \\int_{min_cor}^{max_cor} F_i(x)``.
#
# == return 
# A vector of numeric values.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# set.seed(12345)
# nr1 = 100
# mat1 = matrix(rnorm(100*nr1), nrow = nr1)
# 
# nr2 = 10
# require(mvtnorm)
# sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
# mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))
# 
# nr3 = 50
# sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
# mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))
# 
# mat = rbind(mat1, mat2, mat3)
# AAC_score = AAC(mat)
# plot(AAC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
AAC = function(mat, cor_method = "pearson", min_cor = 0, max_cor = 1, 
	mc.cores = 1, n_sampling = 1000, q_sd = 0) {

	# internally we do it by columns to avoid too many t() callings
	mat = t(mat)

	col_sd = colSds(mat)
	l = col_sd >= quantile(col_sd, q_sd)
	v2 = numeric(length(col_sd))
	v2[!l] = -Inf

	mat = mat[, l, drop = FALSE]
	
	n = ncol(mat)
	if(mc.cores > 1) {
		le = cut(1:n, mc.cores)
		ind_list = split(1:n, le)
	} else {
		ind_list = list(1:n)
	}

	if(cor_method == "KNN_weighted") {
		mat = scale(mat)
	}

	v_list = mclapply(ind_list, function(ind) {
		v = numeric(length(ind))
		for(i in seq_along(ind)) {
			ind2 = seq_len(ncol(mat))[-ind[i]]
			if(length(ind2) > n_sampling) {
				ind2 = sample(ind2, n_sampling)
			}
			# if(cor_method == "KNN_weighted") {
			# 	cor_v = cor_KNN_weighted(mat, ind[i], ind2)
			# } else {
				suppressWarnings(cor_v <- abs(cor(mat[, ind[i], drop = FALSE], mat[, ind2, drop = FALSE], method = cor_method)))
			# }

			cor_v = cor_v[cor_v >= min_cor & cor_v <= max_cor]

			if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
				v[i] = 1
			} else {
				f = ecdf(cor_v)
				cor_v = seq(min_cor, max_cor, length = 100)
				n2 = length(cor_v)
				v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
			}
		}
		return(v)
	}, mc.cores = mc.cores)

	v = do.call("c", v_list)
	v = (max_cor - min_cor) - v
	names(v) = NULL

	v2[l] = v
	return(v2)
}

# cor_KNN_weighted = function(mat, ind1, ind2, k = 5) {

# 	nx = length(ind1)
# 	ny = length(ind2)
# 	cm = matrix(nrow = nx, ncol = ny)
# 	for(i in seq_along(ind1)) {
# 		for(j in seq_along(ind2)) {
# 			cm[i, j] = cor_KNN_weighted_with_two_vectors(mat[, c(ind1[i], ind2[j])], k = k)
# 		}
# 	}
	
# 	return(cm)
# }

# # m: 2 column matrix
# cor_KNN_weighted_with_two_vectors = function(m, k = 5) {
# 	dist_m = base::as.matrix(dist(m, diag = TRUE, upper = TRUE))
# 	wt = apply(dist_m, 1, function(x) {
# 		sum(.Internal(qsort(x, FALSE))[1:(k+1)])/k
# 	})
# 	wt = max(wt) - wt
# 	weightedCorr(m[, 1], m[, 2], weights = wt, method = "pearson")
# }


entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}

# == title
# PAC score
#
# == param
# -consensus_mat a consensus matrix.
# -x1 lower bound to define "ambiguous clustering". The value can be a vector.
# -x2 upper bound to define "ambihuous clustering". The value can be a vector.
# -trim percent of extreme values to trim if combinations of ``x1`` and ``x2`` are more than 10.
# 
# == details
# This a variant of the orignial PAC (proportion of ambiguous clustering) method.
# 
# For each ``x_1i`` in ``x1`` and ``x_2j`` in ``x2``, ``PAC_k = F(x_2j) - F(x_1i)``
# where ``F(x)`` is the ecdf of the consensus matrix (the lower triangle matrix without diagnals). 
# The final PAC is the mean of all ``PAC_k`` by removing top ``trim/2`` percent and bottom ``trim/2`` percent of all values.
# 
# == return 
# A single numeric score.
#
# == see also
# See https://www.nature.com/articles/srep06207 for explanation of PAC score.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == example
# data(cola_rl)
# PAC(get_consensus(cola_rl[1, 1], k = 2))
# PAC(get_consensus(cola_rl[1, 1], k = 3))
# PAC(get_consensus(cola_rl[1, 1], k = 4))
# PAC(get_consensus(cola_rl[1, 1], k = 5))
# PAC(get_consensus(cola_rl[1, 1], k = 6))
#
PAC = function(consensus_mat, x1 = seq(0.1, 0.3, by = 0.02),
	x2 = seq(0.7, 0.9, by = 0.02), trim = 0.2) {

	f = ecdf(consensus_mat[lower.tri(consensus_mat)])

	offset1 = x1
	offset2 = x2

	m_score = matrix(nrow = length(offset1), ncol = length(offset2))
	for(i in seq_along(offset1)) {
		for(j in seq_along(offset2)) {
			m_score[i, j] = f(offset2[j]) - f(offset1[i])
		}
	}
	if(length(m_score)*trim > 1) {
		mean(m_score, trim = trim)
	} else {
		mean(m_score)
	}
}


tot_withinss = function(class_id, mat) {
	s = tapply(seq_along(class_id), class_id, function(ind) {
		if(length(ind) == 1) {
			return(0)
		} else {
			m = mat[, ind, drop = FALSE]
			mean = rowMeans(m)
			m = cbind(mean, m)
			sum((dist(t(m))[1:length(ind)])^2)
		}
	})
	sum(s)
}


# pairwise_concordance_to_class = function(consensus_mat, class, reverse = FALSE) {
# 	class_level = unique(class)
# 	mat_logi = matrix(FALSE, nrow = nrow(consensus_mat), ncol = ncol(consensus_mat))
# 	for(cl in class_level) {
# 		l = class == cl
# 		mat_logi[l, l] = TRUE
# 	}
# 	diag(mat_logi) = FALSE
# 	if(reverse) {
# 		mat_logi2 = !mat_logi
# 		diag(mat_logi2) = FALSE
# 		mean(consensus_mat[mat_logi2])
# 	} else {
# 		mean(consensus_mat[mat_logi])
# 	}
# }

# concordance = function(consensus_mat, class) {
# 	pairwise_concordance_to_class(consensus_mat, class) - pairwise_concordance_to_class(consensus_mat, class, reverse = TRUE)
# }

# separation_rate = function(object, k = 3) {
	
# 	cl2 = get_classes(object, k)$class
# 	if(k == 2) {
# 		cl1 = rep(1, length(cl2))
# 	} else {
# 		cl1 = get_classes(object, k - 1)$class
# 	}

# 	m1 = outer(cl1, cl1, "==")
# 	m2 = outer(cl2, cl2, "==")

# 	m1[lower.tri(m1, diag = TRUE)] = FALSE
# 	m2[lower.tri(m2, diag = TRUE)] = FALSE

# 	sum(m1 & !m2)/sum(m1)
# }


cophcor = function(consensus_mat) {
	dis_consensus = as.dist(1 - consensus_mat)
	hc = hclust(dis_consensus, method = "average")
	cor(dis_consensus, cophenetic(hc))
}

# == title
# Concordance of partitions to the consensus partition
#
# == param
# -membership_each all repetetive parititions.
# -class consensus class ids.
#
# == details
# Class ids in ``membership_meach`` have already be adjusted to the consensus class ids
# to let ``sum(x1 == x_consensus)`` to get maximum.
#
# The concordance score is the mean probability of fiting the consensus class ids in all
# partitions.
#
# This function is used internally.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
concordance = function(membership_each, class) {
	mean(apply(membership_each, 2, function(x) sum(x == class, na.rm = TRUE)/length(x)))
}
