

# == title
# AAC score
#
# == param
# -mat a numeric matrix. AAC score is calculated by columns.
# -cor_method pass to `stats::cor`.
# -min_cor minimal absolute correlation.
# -max_cor maximal absolute correlation.
# -mc.cores number of cores.
# -n_sampling when the number of columns are too high, to get the curmulative
#           distribution, actually we don't need to use all the columns, e.g. 1000
#           columns can already give a farely nice estimation for the distribution.
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
# AAC_score = AAC(t(mat))
# plot(AAC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
AAC = function(mat, cor_method = "pearson", min_cor = 0, max_cor = 1, 
	mc.cores = 1, n_sampling = 1000, q_sd = 0) {

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
	v = (1 - min_cor) - v
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
PAC = function(consensus_mat) {
	f = ecdf(consensus_mat[lower.tri(consensus_mat)])

	offset1 = seq(0, 0.3, by = 0.05)
	offset2 = 1 - seq(0, 0.3, by = 0.05)

	m_score = matrix(nrow = length(offset1), ncol = length(offset2))
	for(i in seq_along(offset1)) {
		for(j in seq_along(offset2)) {
			x = seq(offset1[i], offset2[j], length = 100)
			n2 = length(x)
			v = sum((x[2:n2] - x[1:(n2-1)])*f(x[-n2]))
			rg = offset2[j] - offset1[i]

			m_score[i, j] = abs(v - rg*f(offset1[i]))/rg
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
