
add_transparency = function (col, transparency = 0) {
    rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
}

correspond_between_two_rankings = function(x1, x2, name1 = "", name2 = "", col1 = 1, col2 = 2, top_n = length(x1)) {
	r1 = rank(x1, ties.method = "random")
	r2 = rank(x2, ties.method = "random")

	n = length(x1)
	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3), height = unit(1, "npc") - unit(2, "cm")))
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
		xscale = c(0, 1), yscale = c(0, n + 1)))
	l = r2 >= n - top_n
	grid.points(runif(sum(l)), r1[l], pch = 16, size = unit(1, "mm"), gp = gpar(col = add_transparency(col2, 0.8)))
	grid.text(name1, x = 1, y = unit(n + 1, "native") + unit(2, "mm"), default.units = "native", just = c("right", "bottom"))
	upViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, 
		xscale = c(0, 1), yscale = c(0, n + 1)))
	l = r1 >= n - top_n
	grid.points(runif(sum(l)), r2[l], pch = 16, size = unit(1, "mm"), gp = gpar(col = add_transparency(col1, 0.8)))
	grid.text(name2, x = 0, y = unit(n + 1, "native") + unit(2, "mm"), default.units = "native", just = c("left", "bottom"))
	upViewport()

	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, xscale = c(0, 1), yscale = c(0, n + 1)))
	l = r1 >= n - top_n | r2 >= n - top_n
	# if(sum(!l)) grid.segments(0, r1[!l], 1, r2[!l], default.units = "native", gp = gpar(col = "#EEEEEE80"))
	if(sum(l)) grid.segments(0, r1[l], 1, r2[l], default.units = "native", gp = gpar(col = "#00000004"))
	grid.segments(c(0, 1), c(1, 1), c(0, 1), c(n - top_n, n - top_n), default.units = "native", gp = gpar(col = "#EEEEEE"))
	grid.segments(c(0, 1), c(n - top_n, n - top_n), c(0, 1), c(n, n), default.units = "native", gp = gpar(lwd = 4, col = c(col1, col2)))
	upViewport()

	upViewport()
}

correspond_between_rankings = function(lt, top_n = length(lt[[1]]), 
	col = brewer.pal(length(lt), "Set1")) {

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
			correspond_between_two_rankings(lt[[i]], lt[[j]], nm[i], nm[j], col[i], col[j], top_n)
			upViewport()
			upViewport()
		}
	}
	upViewport()
}
