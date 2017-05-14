
# == title
# Calculate AAC value for each column in the matrix
#
# == param
# -mat a numeric matrix, AAC is calculated for each column
# -cor_method pass to `stats::cor`
# -min_cor minimal absolute correlation
#
AAC_cor = function(mat, cor_method = "pearson", min_cor = 0.2) {
	n = ncol(mat)
	v = numeric(n)
	for(i in 1:n) {
		suppressWarnings(cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE], method = cor_method)))
		if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
			v[i] = 1
		} else {
			f = ecdf(cor_v)
			cor_v = seq(min_cor, 1, length = 100)
			n2 = length(cor_v)
			v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
		}
		# if(interactive()) {
		# 	cat(strrep("\b", 100))
		# 	cat("AAC of the ecdf of correlations for ", i, "/", n, sep = "")
		# }
	}
	# if(interactive()) cat("\n")
	v = (1 - min_cor) - v
	return(v)
}

entropy = function(p) {
	n = length(p)
	p = p[p > 0]
	p2 = rep(1/n, n)
	-sum(p*log2(p))/abs(sum(p2*log2(p2)))
}

PAC = function(consensus_mat, x1 = 0.1, x2 = 0.9) {
	Fn = ecdf(consensus_mat[lower.tri(consensus_mat)])
	Fn(x2) - Fn(x1)
}


mean_cophcor = function(res_list, top_method = res_list$top_method, partition_method = res_list$partition_method) {
	
	df = data.frame(top_method = character(0), partition_method = character(0), k = numeric(0), mean_cophcor = numeric(0))
	for(i in seq_along(top_method)) {
	    for(j in seq_along(partition_method)) {    
	    	res = get_single_run(res_list, top_method = top_method[i], partition_method = partition_method[j])

	        for(ik in seq_along(res$k)) {
	        	k = res$k[ik]
	        	df = rbind(df, data.frame(top_method = top_method[i], partition_method = partition_method[j], k = k, 
	        		mean_cophcor = cophcor(res$object_list[[ik]]$consensus)))
	        }
	    }
	}

	gp = ggplot(df, aes(x = k, y = mean_cophcor)) + geom_line() + geom_point() + facet_grid(top_method ~ partition_method)
	print(gp)
}
