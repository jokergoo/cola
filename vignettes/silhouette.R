library(circlize)
library(GetoptLong)

random_points = function(n, x = 0, y = 0, r = 1) {
	data.frame(x = runif(n, min = -0.5, max = 0.5)*r + x,
	           y = runif(n, min = -0.5, max = 0.5)*r + y)

}

mean_dist = function(x, y, df) {
	d = (nrow(df))
	for(i in seq_len(nrow(df))) {
		d[i] = sqrt((x - df[i, 1]) ^ 2 + (y - df[i, 2]) ^ 2)
	}
	mean(d)
}


make_silhouette_example = function() {
	set.seed(1234)
	df1 = random_points(5, x = -2, y = 0, r = 1)
	df2 = random_points(5, x = 2, y = 2)
	df3 = random_points(5, x = 1, y = -1)

	plot(rbind(df1, df2, df3), col = rep(2:4, each = 5), pch = 16, axes = FALSE, ann = FALSE, asp = 1)

	i = which.max(df1[, 1])
	segments(df1[i, 1], df1[i, 2], df1[-i, 1], df1[-i, 2], col = "grey")
	segments(df1[i, 1], df1[i, 2], df2[, 1], df2[, 2], col = "grey")
	segments(df1[i, 1], df1[i, 2], df3[, 1], df3[, 2], col = "grey")
	points(rbind(df1, df2, df3), col = rep(2:4, each = 5), pch = 16)
	draw.sector(center = c(-2, 0), rou1 = 1)
	draw.sector(center = c(2, 2), rou1 = 1)
	draw.sector(center = c(1, -1), rou1 = 1)
	text(df1[i, 1], df1[i, 2], 'x_i', adj = c(-0.1, -0.1))
	text(mean(df1[, 1]), mean(df1[, 2]), "group1")
	text(mean(df2[, 1]), mean(df2[, 2]), "group2")
	text(mean(df3[, 1]), mean(df3[, 2]), "group3")


	m1 = mean_dist(df1[i, 1], df1[i, 2], df1[-i, ])
	m2 = mean_dist(df1[i, 1], df1[i, 2], df2[, ])
	m3 = mean_dist(df1[i, 1], df1[i, 2], df3[, ])

	txt = qq("mean_1 = @{sprintf('%.2f', m1)}\nmean_2 = @{sprintf('%.2f', m2)}\nmean_3 = @{sprintf('%.2f', m3)}
	silhouette coefficient = 1 - mean_1/min(mean_2, mean_3) = @{sprintf('%.2f', 1 - m1/min(c(m2, m3)))}")
	par(xpd = NA)
	text(-3, -2.5, txt, adj = c(0, 0), cex = 0.7)
}

make_silhouette_example()
