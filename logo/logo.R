bg = c("#51B948", "#EE3F22", "#2B4099", "#FCAF24")

library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))

temp_file = tempfile()
rsvg::rsvg_svg("~/project/cola/logo/can_pink.svg", temp_file)
image = grImport2::readPicture(temp_file)
file.remove(temp_file)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.rect(gp = gpar(fill = bg[1], col = NA))
grImport2::grid.picture(image)
popViewport()

temp_file = tempfile()
rsvg::rsvg_svg("~/project/cola/logo/can_blue.svg", temp_file)
image = grImport2::readPicture(temp_file)
file.remove(temp_file)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.rect(gp = gpar(fill = bg[2], col = NA))
grImport2::grid.picture(image)
popViewport()

temp_file = tempfile()
rsvg::rsvg_svg("~/project/cola/logo/can_orange.svg", temp_file)
image = grImport2::readPicture(temp_file)
file.remove(temp_file)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.rect(gp = gpar(fill = bg[3], col = NA))
grImport2::grid.picture(image)
popViewport()

temp_file = tempfile()
rsvg::rsvg_svg("~/project/cola/logo/can_red.svg", temp_file)
image = grImport2::readPicture(temp_file)
file.remove(temp_file)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.rect(gp = gpar(fill = bg[4], col = NA))
grImport2::grid.picture(image)
popViewport()

popViewport()




# rsvg::rsvg_svg("~/project/cola/logo/can_pink.svg", temp_file)
# image = grImport2::readPicture(temp_file)
# file.remove(temp_file)

# n = length(image@content[[1]]@content)
# path = vector("list", n)
# for(i in seq_len(n)) {
# 	segments = image@content[[1]]@content[[i]]@d@segments

# 	path[[i]] = list(x = numeric(0), y = numeric(0))

# 	for(j in seq_along(segments)) {
# 		path[[i]]$x = c(path[[i]]$x, segments[[j]]@x)
# 		path[[i]]$y = c(path[[i]]$y, segments[[j]]@y)
# 	}
# }

# xscale = image@summary@xscale
# yscale = image@summary@yscale

# grid.newpage()
# pushViewport(viewport(xscale = xscale, yscale = yscale))
# for(i in seq_along(path)) {
# 	grid.polygon(x = path[[i]]$x, y = path[[i]]$y, default.units = "native")
# }
# popViewport()


