

# == title
# Ability to correlate to other rows - an approximate method
#
# == param
# -mat A numeric matrix. ATC score is calculated by rows.
# -cor_fun A function which calculates correlations.
# -min_cor Cutoff for the minimal absolute correlation.
# -power Power on the correlation values.
# -k_neighbours Nearest k neighbours. Note when this argument is set, there won't be subset sampling for calculating
#     correlations, whihc means, it will calculate correlation to all other rows.
# -mc.cores Number of cores. This argument will be removed in future versions.
# -cores Number of cores.
# -n_sampling When there are too many rows in the matrix, to get the curmulative
#           distribution of how one row correlates other rows, actually we don't need to use all the rows in the matrix, e.g. 1000
#           rows can already give a very nice estimation.
# -group A categorical variable. If it is specified, the correlation is only calculated for the rows in the same group as current row.
# -... Pass to ``cor_fun``.
#
# == details
# For a matrix with huge number of rows. It is not possible to calculate correlation to all other rows, thus the correlation is only
# calculated for a randomly sampled subset of othe rows.
#
ATC_approx = function(mat, cor_fun = stats::cor, min_cor = 0, power = 1, k_neighbours = NULL,
	mc.cores = 1, cores = mc.cores, n_sampling = c(1000, 500), 
	group = NULL, ...) {

	if(!is.null(group)) {
		if(length(group) != nrow(mat)) {
			stop_wrap("Length of `group` should be equal to nrow of the matrix.")
		}

		ind_list = split(1:nrow(mat), group)
		sl = lapply(ind_list, function(ind) {
			if(length(ind) == 1) {
				return(0)
			} else {
				ATC_approx(mat[ind, , drop = FALSE], cor_fun = cor_fun, min_cor = min_cor, power = power,
					cores = cores, n_sampling = n_sampling, group = NULL, ...)
			}
		})
		v = numeric(nrow(mat))
		for(i in seq_along(ind_list)) {
			v[ ind_list[[i]] ] = sl[[i]]
		}
		return(v)
	}

	if(length(n_sampling) == 1) n_sampling = c(n_sampling, Inf)

	# internally we do it by columns to avoid too many t() callings
	if(ncol(mat) > n_sampling[2]) {
		mat = mat[, sample(ncol(mat), n_sampling[2])]
	}
	mat = t(mat)

	n_cores = get_nc(cores)

	n = ncol(mat)
	if(n_cores > 1) {
		le = cut(1:n, n_cores)
		ind_list = split(1:n, le)
	} else {
		ind_list = list(1:n)
	}

	if(!is.null(k_neighbours)) k_neighbours = min(k_neighbours, ncol(mat)-1)

	registerDoParallel(cores)
	v_list <- foreach(ind = ind_list) %dopar% {
		v = numeric(length(ind))
		for(i in seq_along(ind)) {
			ind2 = seq_len(ncol(mat))[-ind[i]]
			if(length(ind2) > n_sampling[1] && is.null(k_neighbours)) {
				ind2 = sample(ind2, n_sampling[1])
			}
			suppressWarnings(cor_v <- abs(cor_fun(mat[, ind[i], drop = FALSE], mat[, ind2, drop = FALSE], ...)))
			if(!is.null(k_neighbours)) {
				# min_cor = cor_v[order(-cor_v)[k_neighbours]]
				cor_v = sort(cor_v, decreasing = TRUE)[seq_len(k_neighbours)]
				min_cor = 0
			}
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
	}
	stopImplicitCluster()

	v = do.call("c", v_list)
	v = (1 - min_cor) - v
	names(v) = NULL

	return(v)
}


# == title
# Ability to correlate to other rows
#
# == param
# -mat A numeric matrix. ATC score is calculated by rows.
# -cor_fun A function which calculates correlations.
# -min_cor Cutoff for the minimal absolute correlation.
# -power Power on the correlation values.
# -k_neighbours Nearest k neighbours.
# -mc.cores Number of cores. This argument will be removed in future versions.
# -cores Number of cores.
# -group A categorical variable. If it is specified, the correlation is only calculated for the rows in the same group as current row.
# -... Pass to ``cor_fun``.
# 
# == details 
# For a given row in a matrix, the ATC score is the area above the curve of the curmulative density
# distribution of the absolute correlation to all other rows. Formally, if ``F_i(X)`` is the 
# cumulative distribution function of ``X`` where ``X`` is the absolute correlation for row i with power ``power`` (i.e. ``x = cor^power``),
# ``ATC_i = 1 - \\int_{min_cor}^1 F_i(X)``.
#
# By default the ATC scores are calculated by Pearson correlation, to use Spearman correlation, you can register
# a new top-value method by:
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
# == seealso
# https://jokergoo.github.io/cola_supplementary/suppl_1_ATC/suppl_1_ATC.html
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
ATC = function(mat, cor_fun = stats::cor, min_cor = 0, power = 1, k_neighbours = -1, group = NULL, cores = 1, ...) {
	if(nrow(mat) > 30000) {
		ATC_approx(mat, cor_fun = cor_fun, min_cor = min_cor, power = power, k_neighbours = k_neighbours, cores = cores, ...)
	} else {

		if(!is.null(group)) {
			if(length(group) != nrow(mat)) {
				stop_wrap("Length of `group` should be equal to nrow of the matrix.")
			}

			ind_list = split(1:nrow(mat), group)
			sl = lapply(ind_list, function(ind) {
				if(length(ind) == 1) {
					return(0)
				} else {
					ATC(mat[ind, , drop = FALSE], cor_fun = cor_fun, min_cor = min_cor, power = power,
						cores = cores, group = NULL, ...)
				}
			})
			v = numeric(nrow(mat))
			for(i in seq_along(ind_list)) {
				v[ ind_list[[i]] ] = sl[[i]]
			}
			return(v)
		}

		mat = t(mat)
		n = ncol(mat)
		if(k_neighbours > 0) {
			k_neighbours = min(n-1, k_neighbours)
		}

		if(cores > 1) {
			le = cut(1:n, ceiling(sqrt(cores)))
			ind_list = split(1:n, le)

			ind_mat = expand.grid(seq_along(ind_list), seq_along(ind_list))

			corm = matrix(NA_real_, nrow = n, ncol = n)

			registerDoParallel(cores)
			cm_list <- foreach(i = 1:nrow(ind_mat)) %dopar% {
				block_i = ind_mat[i, 1]
				block_j = ind_mat[i, 2]

				ind1 = ind_list[[block_i]]
				ind2 = ind_list[[block_j]]

				if(block_i == block_j) {
					cm = abs(cor_fun(mat[, ind1, drop = FALSE], ...))
				} else {
					cm = abs(cor_fun(mat[, ind1, drop = FALSE], mat[, ind2, drop = FALSE], ...))
				}
				cm
			}
			stopImplicitCluster()

			for(i in 1:nrow(ind_mat)) {
				block_i = ind_mat[i, 1]
				block_j = ind_mat[i, 2]

				ind1 = ind_list[[block_i]]
				ind2 = ind_list[[block_j]]

				corm[ind1, ind2] = cm_list[[i]]
			}

			le = cut(1:n, cores)
			ind_list = split(1:n, le)
			registerDoParallel(cores)
			s_list <- foreach(ind = ind_list) %dopar% {
				rowATC(corm[ind, , drop = FALSE], min_cor = min_cor, power = power, k_neighbours = k_neighbours, self = ind)
			}
			stopImplicitCluster()

			unlist(s_list)

		} else {
			corm = abs(cor_fun(mat, ...))
			rowATC(corm, min_cor = min_cor, power = power, k_neighbours = k_neighbours, self = 1:n)
		}
	}
}

# only for testing purpose
ATC_single_test = function(x, power = 1, min_cor = 0) {
	x = x^power
	f = ecdf(x)
	cor_v = seq(min_cor, 1, length = 1000)
	n2 = length(cor_v)
	v = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
	1 - min_cor - v
}

entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}


PAC_old = function(consensus_mat, x1 = seq(0.1, 0.3, by = 0.02),
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
# -consensus_mat A consensus matrix.
# -x1 Lower bound to define "ambiguous clustering".
# -x2 Upper bound to define "ambihuous clustering".
# -class Subgroup labels. If it is provided, samples with silhouette score less than the 5^th percential are removed from PAC calculation.
# 
# == details 
# The PAC score is defined as F(x2) - F(x1) where F(x) is the CDF of the consensus matrix.
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
# data(golub_cola)
# PAC(get_consensus(golub_cola[1, 1], k = 2))
# PAC(get_consensus(golub_cola[1, 1], k = 3))
# PAC(get_consensus(golub_cola[1, 1], k = 4))
# PAC(get_consensus(golub_cola[1, 1], k = 5))
# PAC(get_consensus(golub_cola[1, 1], k = 6))
#
# # with specifying `class`
# PAC(get_consensus(golub_cola[1, 1], k = 2), 
#     class = get_classes(golub_cola[1, 1], k = 2)[, 1])
PAC = function(consensus_mat, x1 = 0.1, x2 = 0.9, class = NULL) {
	if(is.null(class)) {
		F = ecdf(consensus_mat[lower.tri(consensus_mat)])
		F(x2) - F(x1)
	} else {
		silhouette = cluster::silhouette(class, dist(t(consensus_mat)))[, "sil_width"]
		if(length(unique(class)) == 1) {
			PAC(consensus_mat, x1, x2)
		} else {
			l = silhouette >= quantile(silhouette, 0.05)
			PAC(consensus_mat[l, l, drop = FALSE], x1, x2)
		}
	}
}

stability = function(consensus_mat, x1 = 0.1, x2 = 0.9) {
	1 - PAC(consensus_mat, x1, x2)
}

# == title
# Flatness of the CDF curve
#
# == param
# -consensus_mat A consensus matrix.
# -diff Difference of F(b) - F(a).
#
# == details
# For a in [0, 0.5] and for b in [0.5, 1], the flatness measures
# the flatness of the CDF curve of the consensus matrix. It is 
# calculated as the maximum width that fits F(b) - F(a) <= diff
#
# == value
# A numeric value.
#
# == example
# data(golub_cola)
# FCC(get_consensus(golub_cola[1, 1], k = 2))
# FCC(get_consensus(golub_cola[1, 1], k = 3))
# FCC(get_consensus(golub_cola[1, 1], k = 4))
# FCC(get_consensus(golub_cola[1, 1], k = 5))
# FCC(get_consensus(golub_cola[1, 1], k = 6))
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
# -consensus_mat A consensus matrix.
#
# == details
# For the consensus values x, it is transformed to 1 - x if x < 0.5.
# After the transformation, for any pair of samples in the consensus matrix,
# If they are always in a same group or always in different groups, the 
# value x is both to 1. Thus, if the consensus matrix shows stable partitions,
# values x will be all close to 1. Reflected in the CDF of x, the curve is 
# shifted to the right and the area under CDF curve should be very small.
#
# An aPAC value less than 0.05 is considered as the stable partition, which can
# be thought the proportion of abmiguous partitioning is less than 0.05.
#
# == value
# A numeric value.
#
# == example
# data(golub_cola)
# aPAC(get_consensus(golub_cola[1, 1], k = 2))
# aPAC(get_consensus(golub_cola[1, 1], k = 3))
# aPAC(get_consensus(golub_cola[1, 1], k = 4))
# aPAC(get_consensus(golub_cola[1, 1], k = 5))
# aPAC(get_consensus(golub_cola[1, 1], k = 6))
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
# Concordance to the consensus partition
#
# == param
# -membership_each A matrix which contains partitions in every single runs where columns correspond to runs.
#          The object can be get from ``get_membership(..., each = TRUE)``.
# -class Consensus subgroup labels.
#
# == details
# Note subgroup labels in ``membership_each`` should already be adjusted to the consensus labels, i.e. by `relabel_class`.
#
# The concordance score is the mean proportion of samples having the same subgroup labels as the consensus labels
# among individual partition runs.
#
# == value
# A numeric value.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola)
# membership_each = get_membership(golub_cola["SD", "kmeans"], each = TRUE, k = 3)
# consensus_classes = get_classes(golub_cola["SD", "kmeans"], k = 3)$class
# concordance(membership_each, consensus_classes)
concordance = function(membership_each, class) {
	mean(apply(membership_each, 2, function(x) sum(x == class, na.rm = TRUE)/length(x[!is.na(x)])))
}
