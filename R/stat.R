

# == title
# Ability to correlate other rows in the matrix (ATC score)
#
# == param
# -mat a numeric matrix. ATC score is calculated by rows.
# -cor_fun a function which calculates correlation.
# -min_cor minimal absolute correlation.
# -power power on the correlation values.
# -mc.cores number of cores.
# -n_sampling when there are too many rows in the matrix, to get the curmulative
#           distribution of how one row correlates other rows, actually we don't need to use all the rows in the matrix, e.g. 1000
#           rows can already give a farely nice estimation.
# -q_sd percentile of the standard deviation for the rows. Rows with values less than it are ignored.
# -... pass to ``cor_fun``, e.g. ``method = 'spearman'`` can be passed to ``cor_fun`` if the correlation function is `stats::cor`.
# 
# == details 
# For a given row in a matrix, the ATC score is the area above the curve of the curmulative density
# distribution of the absolute correlation to all other rows. Formally, if ``F_i(X)`` is the 
# cumulative distribution function of ``X`` where ``X`` is the absolute correlation for row i with power ``power`` (i.e. ``x = cor^power``),
# ``ATC_i = 1 - \\int_{min_cor}^1 F_i(X)``.
#
# By default the ATC scores are calculated by Pearson correlation, to use Spearman correlation, you can register
# the top-value method by:
#
#     register_top_value_methods(
#         "ATC_spearman" = function(m) ATC(m, method = "spearman")
#     )
#
# Similarly, to use a robust correlation method, e.g. `WGCNA::bicor` function, you can do like:
#
#     register_top_value_methods(
#         "ATC_bicor" = function(m) ATC(m, cor_fun = WGCNA::bicor)
#     )
#
# == return 
# A vector of numeric values with the same order as rows in the input matrix.
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
# ATC_score = ATC(mat)
# plot(ATC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
ATC = function(mat, cor_fun = stat::cor, min_cor = 0, power = 1,
	mc.cores = 1, n_sampling = 1000, q_sd = 0, ...) {

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

	v_list = mclapply(ind_list, function(ind) {
		v = numeric(length(ind))
		for(i in seq_along(ind)) {
			ind2 = seq_len(ncol(mat))[-ind[i]]
			if(length(ind2) > n_sampling) {
				ind2 = sample(ind2, n_sampling)
			}
			suppressWarnings(cor_v <- abs(cor(mat[, ind[i], drop = FALSE], mat[, ind2, drop = FALSE], ...)))
			cor_v = cor_v^power
			
			if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
				v[i] = 1
			} else {
				f = ecdf(cor_v)
				cor_v = seq(min_cor, 1, length = 1000)
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

entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}

# == title
# The proportion of ambiguous clustering (PAC score)
#
# == param
# -consensus_mat a consensus matrix.
# -x1 lower bound to define "ambiguous clustering". The value can be a vector.
# -x2 upper bound to define "ambihuous clustering". The value can be a vector.
# -trim percent of extreme values to trim if combinations of ``x1`` and ``x2`` are more than 10.
# 
# == details
# This a variant of the orignial PAC method.
# 
# For each ``x_1i`` in ``x1`` and ``x_2j`` in ``x2``, ``PAC_k = F(x_2j) - F(x_1i)``
# where ``F(x)`` is the cumulative distribution function of the consensus matrix (the lower triangle matrix without diagnals is only used). 
# The final PAC is the mean of all ``PAC_k`` by removing top ``trim/2`` percent and bottom ``trim/2`` percent of all values.
# 
# == return 
# A single numeric vaule.
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

	F = ecdf(consensus_mat[lower.tri(consensus_mat)])

	offset1 = x1
	offset2 = x2

	m_score = matrix(nrow = length(offset1), ncol = length(offset2))
	for(i in seq_along(offset1)) {
		for(j in seq_along(offset2)) {
			m_score[i, j] = F(offset2[j]) - F(offset1[i])
		}
	}
	if(length(m_score)*trim > 1) {
		mean(m_score, trim = trim)
	} else {
		mean(m_score)
	}
}

PAC2 = function(consensus_mat, x1 = seq(0.1, 0.3, by = 0.02),
	x2 = seq(0.7, 0.9, by = 0.02)) {

	F = ecdf(consensus_mat[lower.tri(consensus_mat)])

	n1 = length(x1)
	n2 = length(x2)

	if(n1 == 1) {
		a1 = 0
	} else {
		a1 = sapply(seq_len(n1 - 1), function(i) F(x1[i+1]) - F(x1[i]))
	}

	x2 = rev(x2)
	if(n2 == 1) {
		a2 = 0
	} else {
		a2 = sapply(seq_len(n2 - 1), function(i) F(x2[i]) - F(x2[i+1]))
	}

	a3 = F(x2[n2]) - F(x1[n1])

	sum(a1*(seq_len(n1 - 1)/n1)) + sum(a2*(seq_len(n2 - 1)/n2)) + a3
}


# == title
# The proportion of ambiguous clustering (PAC score)
#
# == param
# -consensus_mat a consensus matrix.
# -x1 lower bound to define "ambiguous clustering".
# -x2 upper bound to define "ambihuous clustering".
# 
# == details
# This is the original implementation of PAC method.
#
PAC_origin = function(consensus_mat, x1 = 0.1, x2 = 0.9) {
	F = ecdf(consensus_mat[lower.tri(consensus_mat)])
	F(x2) - F(x1)
}

stability = function(consensus_mat, x1 = 0.1, x2 = 0.9) {
	1 - PAC_origin(consensus_mat, x1, x2)
}

# == title
# Flatness of the CDF curve
#
# == param
# -consensus_mat a consensus matrix.
# -diff Difference of F(b) - F(a)
#
# == details
# For a in [0, 0.5] and for b in [0.5, 1], the flatness measures
# the flatness of the CDF curve of the consensus matrix, it is 
# calculated as the maximum width that fits F(b) - F(a) <= diff
#
# A flatness larger than 0.9 is treated as stable partitions.
#
FCC = function(consensus_mat, diff = 0.1) {
	F = ecdf(consensus_mat[lower.tri(consensus_mat)])
	a = seq(0, 0.5, length = 100)
	b = seq(0.5, 1, length = 100)

	grid = expand.grid(a, b)
	x = F(grid[, 2]) - F(grid[, 1])
	ind = which(x <= diff)
	max(grid[ind, 2] - grid[ind, 1])
}

# == title
# Adapted PAC scores
#
# == param
# -consensus_mat a consensus matrix.
#
# == details
# For the consensus values x, it is transformed to 1 - x if x < 0.5.
# After the transformation, for any pair of samples in the consensus matrix,
# If they are always in a same group or always in different groups, the 
# value x is both to 1. Thus, if the consensus matrix shows stable partitions,
# values x will be all close to 1. Reflecting to the CDF of x, the curve is 
# shifted to the right and the area under CDF curve should be very small.
#
# An aPAC value less than 0.05 is considered as stable partitions, which can
# be thought the proportion of abmiguous partitioning is less than 0.05.
#
aPAC = function(consensus_mat) {
	x = consensus_mat[lower.tri(consensus_mat)]
	x = ifelse(x < 0.5, 1 - x, x)
	x = (x - 0.5)*2
	F = ecdf(x)
	breaks = seq(0, 1, length = 1000)
	sum(F(breaks)*(breaks[2] - breaks[1]))
}

mean_group_dist = function(mat, factor) {
	le = unique(factor)
	nle = length(le)
	if(nle == 1) {
		return(0)
	}
	d1 = matrix(nrow = nle, ncol = nle)
	dm = as.matrix(dist(mat))

	for(i in 1:(nle-1)) {
		for(j in (i+1):nle) {
			l1 = factor == le[i]
			l2 = factor == le[j]
			d1[i, j] = mean(dm[l1, l2])
		}
	}
	mean(d1[!is.na(d1)])
}

cophcor = function(consensus_mat) {
	dis_consensus = as.dist(1 - consensus_mat)
	hc = hclust(dis_consensus, method = "average")
	cor(dis_consensus, cophenetic(hc))
}

# == title
# Concordance of partitions to the consensus partition
#
# == param
# -membership_each a matrix which contains partitions in every single runs.
# -class consensus class IDs.
#
# == details
# Class IDs in ``membership_each`` have already be adjusted to the consensus class IDs
# to let ``sum(x_single == x_consensus)`` reach maximum.
#
# The concordance score is the mean probability of fitting the consensus class IDs in all
# partitions.
#
# This function is used internally.
#
# == value
# A numeric value.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(cola_rl)
# membership_each = get_membership(cola_rl["sd", "kmeans"], each = TRUE, k = 3)
# consensus_classes = get_classes(cola_rl["sd", "kmeans"], k = 3)$class
# concordance(membership_each, consensus_classes)
concordance = function(membership_each, class) {
	mean(apply(membership_each, 2, function(x) sum(x == class, na.rm = TRUE)/length(x)))
}
