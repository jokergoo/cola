
ATC_definition = function() {
	cor_v = c(0.00296, 0.00323, 0.00493, 0.00511, 0.00513, 0.00653, 0.00838, 0.00859, 0.00923, 0.0116, 0.0141, 0.0142, 0.0144, 0.0161, 0.0186, 0.0195, 0.0205, 0.0222, 0.0231, 0.0237, 0.0242, 0.025, 0.0255, 0.0277, 0.028, 0.029, 0.029, 0.0314, 0.0322, 0.0355, 0.0357, 0.0363, 0.0379, 0.0392, 0.0411, 0.0423, 0.0426, 0.0441, 0.0449, 0.046, 0.0501, 0.0508, 0.0529, 0.0536, 0.054, 0.0541, 0.0549, 0.0567, 0.0576, 0.0581, 0.0605, 0.0608, 0.0639, 0.0671, 0.0682, 0.0698, 0.0717, 0.0737, 0.0747, 0.075, 0.0773, 0.0816, 0.0834, 0.0837, 0.0869, 0.0899, 0.0928, 0.0937, 0.0945, 0.0946, 0.0964, 0.0981, 0.099, 0.101, 0.102, 0.102, 0.103, 0.104, 0.106, 0.109, 0.111, 0.115, 0.115, 0.116, 0.116, 0.121, 0.124, 0.127, 0.127, 0.127, 0.133, 0.135, 0.135, 0.135, 0.138, 0.143, 0.144, 0.146, 0.156, 0.16, 0.164, 0.169, 0.173, 0.178, 0.179, 0.179, 0.194, 0.202, 0.211, 0.253, 0.396, 0.429, 0.447, 0.455, 0.457, 0.465, 0.466, 0.477, 0.487, 0.488, 0.498, 0.518, 0.519, 0.522, 0.524, 0.526, 0.528, 0.539, 0.54, 0.541, 0.544, 0.546, 0.548, 0.551, 0.552, 0.556, 0.558, 0.559, 0.565, 0.568, 0.568, 0.574, 0.574, 0.576, 0.576, 0.579, 0.58, 0.582, 0.585, 0.59, 0.591, 0.592, 0.594, 0.596, 0.597, 0.609, 0.614, 0.615, 0.628)

	cor_v = sort(cor_v)
	f = ecdf(cor_v)
	x = seq(0, 1, length = 100)
	y = f(x)

	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "P(X <= x)")
	x2 = seq(0, 1, length = 100)
	polygon(c(x2, rev(x2)), c(f(x2), rep(1, length(x2))), col = "#FF000040", border = NA)
	polygon(c(x2, rev(x2)), c(f(x2), rep(0, length(x2))), col = "#00FF0040", border = NA)
	lines(x, y)
	text(0.4, 0.9, "ATC")
}


ATC_simulation = function() {
	set.seed(12345)
	library(matrixStats)

	nr1 = 100
	mat1 = matrix(rnorm(100*nr1), nrow = nr1)

	nr2 = 10
	library(mvtnorm)
	sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
	mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))

	nr3 = 50
	sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
	mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))

	mat = t(rbind(mat1, mat2, mat3))

	# par(mfrow = c(2, 2))
	grid.newpage()
	par(mfrow = c(2, 2))
	plot.new()
	par(mar = c(3, 5, 1, 1), cex = 0.8)
	v = NULL
	min_cor = 0
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "P(X <= x)", xlab = "")
	for(i in 1:ncol(mat)) {
	    cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE]))
	    cor_v = sort(cor_v)
	    f = ecdf(cor_v)
	    # cor_v = c(0, cor_v, 1)
	    # n2 = length(cor_v)
	    # v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
	    cor_v = seq(min_cor, 1, length = 100)
	    n2 = length(cor_v)
	    v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
	    x = seq(0, 1, length = 100)
	    lines(x, f(x), col = ifelse(i <= nr1, "#00000080", ifelse(i <= nr1+nr2, "#FF000080", "#00FF0080")))
	}
	# abline(v = 0.2, lty = 2, col = "grey")
	legend("bottomright", lty = 1, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))
	v = (1 - min_cor) - v

	dx_list = NULL
	dy_list = NULL
	for(i in 1:ncol(mat)) {
	    cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE]))
	    d = density(cor_v)
	    dx_list = c(dx_list, list(d$x))
	    dy_list = c(dy_list, list(d$y))
	}
	max_y = max(unlist(dy_list))
	# plot(NULL, xlim = c(0, 1), ylim = c(0, max_y), ylab = "density", xlab = "")
	# for(i in seq_along(dx_list)) {
	# 	lines(dx_list[[i]], dy_list[[i]], col = ifelse(i <= nr1, "#00000080", ifelse(i <= nr1+nr2, "#FF000080", "#00FF0080")))
	# }
	# legend("topright", lty = 1, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))


	plot(v, pch = 16, xlab = "", ylab = "ATC", col = ifelse(seq_along(v) <= nr1, "#000000", ifelse(seq_along(v) <= nr1+nr2, "#FF0000", "#00FF00")))
	legend("topleft", pch = 16, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))


	plot(colSds(mat), pch = 16, xlab = "", ylab = "SD", col = ifelse(seq_along(v) <= nr1, "#000000", ifelse(seq_along(v) <= nr1+nr2, "#FF0000", "#00FF00")))

	pushViewport(viewport(x = 0, y = 1, width = 0.5, height = 0.5, just = c("left", "top")))
	library(ComplexHeatmap)
	split = factor(c(rep(1, 100), rep(2, 10), rep(3, 50)), levels = c(1, 2, 3))
	ht = Heatmap(t(mat), row_split = split, cluster_row_slices = FALSE, row_title = NULL,
		show_row_dend = FALSE, show_column_dend = FALSE, show_heatmap_legend = FALSE,
		left_annotation = rowAnnotation(split = split, col = list(split = c("1" = 1, "2" = 2, "3" = 3)),
			show_legend = FALSE))
	draw(ht, newpage = FALSE)
	popViewport()
}
