

#' AAC score
#'
#' @param mat a numeric matrix. AAC score is calculated by columns.
#' @param cor_method pass to [stats::cor()].
#' @param min_cor minimal absolute correlation.
#' 
#' @details 
#' AAC score for a given item is the area above the curve of the curmulative 
#' distribution of the absolute correlation to all other items.
#'
#' @return A vector of AAC scores
#' @export
#' @import stats
#' @import parallel
#'
#' @examples
#' set.seed(12345)
#' require(matrixStats)
#' nr1 = 100
#' mat1 = matrix(rnorm(100*nr1), nrow = nr1)
#' 
#' nr2 = 10
#' require(mvtnorm)
#' sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
#' mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))
#' 
#' nr3 = 50
#' sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
#' mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))
#' 
#' mat = rbind(mat1, mat2, mat3)
#' AAC_score = AAC(t(mat))
#' plot(AAC_score, pch = 16, col = c(rep(1, nr1), rep(2, nr2), rep(3, nr3)))
AAC = function(mat, cor_method = "pearson", min_cor = 0.2, mc.cores = 1) {
	n = ncol(mat)
	
	if(mc.cores > 1) {
		le = cut(1:n, mc.cores)
		ind_list = split(n, le)
	} else {
		ind_list = list(1:n)
	}

	v_list = mclapply(ind_list, function(ind) {
		v = numeric(length(ind))
		for(i in ind) {
			suppressWarnings(cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE], method = cor_method)))
			if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
				v[i] = 1
			} else {
				f = ecdf(cor_v)
				cor_v = seq(min_cor, 1, length = 100)
				n2 = length(cor_v)
				v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
			}
		}
		return(v)
	}, mc.cores = mc.cores)

	v = do.call("c", v_list)
	v = (1 - min_cor) - v
	return(v)
}

entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}

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
	mean(m_score)
}

APN = function(class_id, membership_each) {
	ind_list = tapply(seq_along(class_id), class_id, function(x) x, simplify = FALSE)
	s = sum(sapply(ind_list, function(ind) {
		sum(apply(membership_each[ind, , drop = FALSE], 2, function(x) 1 - sum(x == class_id[ind])/length(ind)))
	}))
	s/ncol(membership_each)/length(ind_list)
}
