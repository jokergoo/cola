
add_transparency = function (col, transparency = 0) {
    rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
}

# == title
# Correspond two rankings
#
# == param
# -x1 a vector of scores calculated by one metric.
# -x2 a vector of scores calculated by another metric.
# -name1 name of the first metric.
# -name2 name of the second metric.
# -col1 color of the lines for the first metric.
# -col2 color of the lines for the second metric.
# -top_n top n elements to visualize
# -transparency transparency of the connection lines.
# -pt_size size of the points, must be a `grid::unit` object
# -newpage whether to plot in a new graphic page.
# -ratio ratio of width of the left plot, connection lines and right plot.
# 
# == details
# In ``x1`` and ``x2``, the i^{th} element is the same object but with different 
# scores under different metrics.
# 
# ``x1`` and ``x2`` are sorted in the plot. Lines are connecting the same object
# for those with top values. Shadows on the two sides are the sorted values for ``x1``
# and ``x2``, and points corresponds to positions of the top objects but are randomly shifted
# on x-directions.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == examples
# require(matrixStats)
# mat = matrix(runif(1000), ncol = 10)
# x1 = rowSds(mat)
# x2 = rowMads(mat)
# correspond_between_two_rankings(x1, x2, name1 = "sd", name2 = "mad", top_n = 20)
correspond_between_two_rankings = function(x1, x2, name1 = "", name2 = "", 
	col1 = 1, col2 = 2, top_n = round(0.25*length(x1)), transparency = 0.9, 
	pt_size = unit(1, "mm"), newpage = FALSE, ratio = c(1, 1, 1)) {
	
	if(newpage) {
		grid.newpage()
	}

	if(length(x1) != length(x2)) {
		stop("Length of `x1` and `x2` should be the same.")
	}

	r1 = rank(x1, ties.method = "random")
	r2 = rank(x2, ties.method = "random")

	n = length(x1)
	text_height = grobHeight(textGrob("foo"))*2
	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3, widths = unit(ratio, "null")), 
		height = unit(1, "npc") - text_height - unit(1, "cm"), y = unit(1, "cm"), just = "bottom"))
	
	max_x1 = max(x1)
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
		xscale = c(0, max_x1), yscale = c(0, n + 1)))
	grid.segments(max_x1 - x1, r1, max_x1, r1, default.units = "native", gp = gpar(col = "#EFEFEF"))
	l = r2 >= n - top_n
	grid.points(unit(runif(sum(l)), "npc"), r1[l], default.units = "native", pch = 16, size = pt_size, gp = gpar(col = add_transparency(col2, 0.8)))
	grid.text(name1, x = 1, y = unit(n + 1, "native") + unit(1, "mm"), default.units = "npc", just = c("right", "bottom"))
	upViewport()

	max_x2 = max(x2)
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, 
		xscale = c(0, max_x2), yscale = c(0, n + 1)))
	grid.segments(0, r2, x2, r2, default.units = "native", gp = gpar(col = "#EFEFEF"))
	l = r1 >= n - top_n
	grid.points(unit(runif(sum(l)), "npc"), r2[l], default.units = "native", pch = 16, size = pt_size, gp = gpar(col = add_transparency(col1, 0.8)))
	grid.text(name2, x = 0, y = unit(n + 1, "native") + unit(1, "mm"), default.units = "native", just = c("left", "bottom"))
	upViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, xscale = c(0, 1), yscale = c(0, n + 1)))
	l = r1 >= n - top_n | r2 >= n - top_n
	# if(sum(!l)) grid.segments(0, r1[!l], 1, r2[!l], default.units = "native", gp = gpar(col = "#EEEEEE80"))
	if(sum(l)) grid.segments(0, r1[l], 1, r2[l], default.units = "native", gp = gpar(col = add_transparency("#000000", transparency)))
	grid.segments(c(0, 1), c(1, 1), c(0, 1), c(n - top_n, n - top_n), default.units = "native", gp = gpar(col = "#EEEEEE"))
	grid.segments(c(0, 1), c(n - top_n, n - top_n), c(0, 1), c(n, n), default.units = "native", gp = gpar(lwd = 4, col = c(col1, col2)))
	upViewport()

	upViewport()

	# add a venn diagram at the bottom
	n_intersect = length(intersect(order(x1, decreasing = TRUE)[1:top_n], order(x2, decreasing = TRUE)[1:top_n]))
	n_union = 2*top_n - n_intersect
	grid.roundrect(x = unit(0.25, "npc"), y = unit(0.4, "cm"), width = unit(0.5*top_n/n_union, "npc"), 
		height = unit(0.4, "cm"), gp = gpar(fill = add_transparency(col1, 0.5), col = NA), just = "left")
	grid.roundrect(x = unit(0.75, "npc"), y = unit(0.4, "cm"), width = unit(0.5*top_n/n_union, "npc"), 
		height = unit(0.4, "cm"), gp = gpar(fill = add_transparency(col2, 0.5), col = NA), just = "right")
	grid.text(qq("top @{top_n}"), x = unit(0.5, "npc"), y = unit(0.7, "cm"), just = "bottom", gp = gpar(fontsize = 8))

}

# == title
# Correspond between a list of rankings
#
# == param
# -lt a list of scores under different metrics.
# -top_n top n elements to visualize.
# -col a vector of colors for ``lt``.
# -... pass to `correspond_between_two_rankings`.
# 
# == details
# It makes plots for pairwise comparisons between different rankings.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
# == examples
# require(matrixStats)
# mat = matrix(runif(1000), ncol = 10)
# x1 = rowSds(mat)
# x2 = rowMads(mat)
# x3 = rowSds(mat)/rowMeans(mat)
# correspond_between_rankings(lt = list(sd = x1, mad = x2, vc = x3), 
#     top_n = 20)
correspond_between_rankings = function(lt, top_n = length(lt[[1]]), 
	col = brewer_pal_set1_col[1:length(lt)], ...) {

	nm = names(lt)
	n = length(lt)
	n_plots = n*(n-1)/2

	if(length(col) == 1) {
		col = rep(col, n)
	}

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = n_plots)))
	k = 0
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			k = k + 1
			pushViewport(viewport(layout.pos.row = 1, layout.pos.col = k))
			pushViewport(viewport(width = 0.9))
			correspond_between_two_rankings(lt[[i]], lt[[j]], nm[i], nm[j], col[i], col[j], top_n, ...)
			upViewport()
			upViewport()
		}
	}
	upViewport()
}
