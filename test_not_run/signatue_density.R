
foo = function(n, k, r = seq(1, 2, by = 0.01)) {
	x = replicate(1000, sum(sample(c(TRUE, FALSE), n, replace = TRUE)))
	
	p = sapply(r, function(y) sum(x*y >= k)/length(x))
	fc = k/mean(x*y)
}

r = seq(1, 2, by = 0.01)
plot(r, foo(n, n/2 + 10, r), type = "n", ylim = c(0, 1))
n = c(10, 20, 50, 100, 500, 1000)
for(i in seq_along(n)) {
	lines(r, foo(n[i], n[i]/2 + 10, r), col = i)
}


foo = function(k) {
	cl = get_class(res, k = 3)$class
	gp = x$group

	gp = gp[gp == k]
	den = density(data[, 1])
	plot(den, type = "n", ylim = c(0, 1.5), main = qq("signature genes in subgroup @{k}"))
	pos = find_pos(data[names(gp), 1], den)
	points(pos, col = add_transparency(brewer_pal_set2_col[cl[i]], ifelse(as.numeric(gp) == cl[1], 0.5, 0.9)),
		pch = 16, cex = ifelse(as.numeric(gp) == cl[1], 1, 0.4))

	for(i in 2:ncol(data)) {
		den = density(data[, i])
		pos = find_pos(data[names(gp), i], den)
		points(pos, col = add_transparency(brewer_pal_set2_col[cl[i]], ifelse(as.numeric(gp) == cl[i], 0.5, 0.9)),
			pch = 16, cex = ifelse(as.numeric(gp) == cl[i], 1, 0.4))

	}
}


foo = function(k) {
	cl = get_class(res, k = 3)$class
	gp = x$group

	gp = gp[gp == k]
	plot(NULL, type = "n", xlim = c(0, 4), ylim = c(0, 1.5), main = qq("signature genes in subgroup @{k}"))
	lines(density(data[names(gp), 1]), col = brewer_pal_set2_col[cl[1]], lwd = ifelse(cl[1] == as.numeric(k), 2, 0.5))

	for(i in 2:ncol(data)) {
		lines(density(data[names(gp), i]), col = brewer_pal_set2_col[cl[i]], lwd = ifelse(cl[i] == as.numeric(k), 2, 0.5))

	}
}

par(mfrow = c(3, 1))
foo("1")
foo("2")
foo("3")


find_pos = function(x, den) {
	do.call("rbind", lapply(x, function(x2) {
		i = max(which(den$x <= x2))
		c(den$x[i], den$y[i] + runif(1, min = -0.05, max = 0.05))
	}))
}