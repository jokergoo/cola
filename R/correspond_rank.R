
add_transparency = function (col, transparency = 0) {
    rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
}

#' Correspond two rankings
#'
#' @param x1 a vector of scores calculated by one metric.
#' @param x2 a vector of scores calculated by another metric.
#' @param name1 name of the first metric.
#' @param name2 name of the second metric.
#' @param col1 color of the lines for the first metric.
#' @param col2 color of the lines for the second metric.
#' @param top_n top n elements to visualize
#' 
#' @details
#' In `x1` and `x2`, the i^{th} element is the same object but with different 
#' scores under different metrics.
#'
#' @export
#' @import grid
#'
#' @examples
#' require(matrixStats)
#' mat = matrix(runif(1000), ncol = 10)
#' x1 = rowSds(mat)
#' x2 = rowMads(mat)
#' correspond_between_two_rankings(x1, x2, name1 = "sd", name2 = "mad", top_n = 20)
correspond_between_two_rankings = function(x1, x2, name1 = "", name2 = "", 
	col1 = 1, col2 = 2, top_n = length(x1)) {
	
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

#' Correspond between a list of rankings
#'
#' @param lt a list of scores under different metrics.
#' @param top_n top n elements to visualize.
#' @param col colors
#'
#' @export
#' @import grid
#' @import RColorBrewer
#'
#' @examples
#' require(matrixStats)
#' mat = matrix(runif(1000), ncol = 10)
#' x1 = rowSds(mat)
#' x2 = rowMads(mat)
#' x3 = rowSds(mat)/rowMeans(mat)
#' correspond_between_rankings(lt = list(sd = x1, mad = x2, vc = x3), 
#'     top_n = 20)
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
