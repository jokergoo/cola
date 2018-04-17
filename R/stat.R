

# == title
# AAC score
#
# == param
# -mat a numeric matrix. AAC score is calculated by rows.
# -cor_method pass to `stats::cor`.
# -min_cor minimal absolute correlation.
# -max_cor maximal absolute correlation.
# -mc.cores number of cores.
# -n_sampling when the number of columns are too high, to get the curmulative
#           distribution, actually we don't need to use all the rows, e.g. 1000
#           rows can already give a farely nice estimation for the distribution.
# -q_sd percential of the sd to ignore
# 
# == details 
# AAC score for a given item is the area above the curve of the curmulative 
# distribution of the absolute correlation to all other items with ``x >= min_cor``.
#
# == return 
# A vector of AAC scores.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# set.seed(12345)
# require(matrixStats)
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
# -consensus_mat a consensus matrix
# 
# == details
# This a variant of the orignial PAC (proportion of ambiguous clustering) method.
#
# Assume x_1 in [0, 0.3] and x_2 in [0.7, 1], we calculate s = (\\int_{x_1}^{x_2} F(x) - F(x1)*(x_2 - x_1))/(x_2 - x_1) where F(x) is the CDF of the consensus matrix
# and the mean value is taken as the final value.
# 
# == return 
# A single PAC score.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
PAC = function(consensus_mat, original = FALSE) {
	f = ecdf(consensus_mat[lower.tri(consensus_mat)])

	offset1 = seq(0, 0.3, by = 0.05)
	offset2 = 1 - seq(0, 0.3, by = 0.05)

	m_score = matrix(nrow = length(offset1), ncol = length(offset2))
	for(i in seq_along(offset1)) {
		for(j in seq_along(offset2)) {
			if(original) {
				m_score[i, j] = f(offset2[j]) - f(offset1[i])
			} else {
				x = seq(offset1[i], offset2[j], length = 100)
				n2 = length(x)
				v = sum((x[2:n2] - x[1:(n2-1)])*f(x[-n2]))
				rg = offset2[j] - offset1[i]

				m_score[i, j] = abs(v - rg*f(offset1[i]))/rg
			}
		}
	}
	mean(m_score, trim = 0.2)
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


pairwise_concordance_to_class = function(consensus_mat, class, reverse = FALSE) {
	class_level = unique(class)
	mat_logi = matrix(FALSE, nrow = nrow(consensus_mat), ncol = ncol(consensus_mat))
	for(cl in class_level) {
		l = class == cl
		mat_logi[l, l] = TRUE
	}
	diag(mat_logi) = FALSE
	if(reverse) {
		mat_logi2 = !mat_logi
		diag(mat_logi2) = FALSE
		mean(consensus_mat[mat_logi2])
	} else {
		mean(consensus_mat[mat_logi])
	}
}

rand_index = function(object, k = 3) {
	if(k == 2) {
		return(NA)
	}

	membership_each1 = get_membership(object, k - 1, each = TRUE)
	membership_each2 = get_membership(object, k, each = TRUE)
	mat11 = matrix(1, nrow = nrow(membership_each1), ncol = nrow(membership_each1))
	mat01 = mat11
	mat10 = mat11
	mat00 = mat11
	for(i in 1:(nrow(membership_each1)-1)) {
		for(j in (i+1):nrow(membership_each1)) {
			mat11[i, j] = sum(membership_each1[i, ] == membership_each1[j, ] & membership_each2[i, ] == membership_each2[j, ])/ncol(membership_each1)
			mat11[j, i] = mat11[i, j]

			mat01[i, j] = sum(membership_each1[i, ] != membership_each1[j, ] & membership_each2[i, ] == membership_each2[j, ])/ncol(membership_each1)
			mat01[j, i] = mat01[i, j]

			mat10[i, j] = sum(membership_each1[i, ] == membership_each1[j, ] & membership_each2[i, ] != membership_each2[j, ])/ncol(membership_each1)
			mat10[j, i] = mat10[i, j]

			mat00[i, j] = sum(membership_each1[i, ] != membership_each1[j, ] & membership_each2[i, ] != membership_each2[j, ])/ncol(membership_each1)
			mat00[j, i] = mat00[i, j]
		}
 	}
 	a = mean(mat00[lower.tri(mat00)])
 	b = mean(mat01[lower.tri(mat00)])
 	c = mean(mat10[lower.tri(mat00)])
 	d = mean(mat11[lower.tri(mat00)])

 	a
}

separation_rate = function(object, k = 3) {
	
	cl2 = get_classes(object, k)$class
	if(k == 2) {
		cl1 = rep(1, length(cl2))
	} else {
		cl1 = get_classes(object, k - 1)$class
	}

	m1 = outer(cl1, cl1, "==")
	m2 = outer(cl2, cl2, "==")

	m1[lower.tri(m1, diag = TRUE)] = FALSE
	m2[lower.tri(m2, diag = TRUE)] = FALSE

	sum(m1 & !m2)/sum(m1)
}


cophcor = function(consensus_mat) {
	dis_consensus = as.dist(1 - consensus_mat)
	hc = hclust(dis_consensus, method = "average")
	cor(dis_consensus, cophenetic(hc))
}


