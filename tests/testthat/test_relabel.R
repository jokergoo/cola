

class = c(rep("a", 10), rep("b", 3))
ref = c(rep("b", 4), rep("a", 9))
map = relabel_class(class, ref)
attr(map, "df") = NULL

test_that("test relabel", {
	expect_equal(map, c("a" = "b", "b" = "a"))
})


adjusted = relabel_class(class, ref, return_map = FALSE)

test_that("test relabel", {
	expect_equal(adjusted, unname(map[class]))
})



class = c(rep(1, 10), rep(2, 3))
ref = c(rep(2, 4), rep(1, 9))
map = relabel_class(class, ref)
adjusted = relabel_class(class, ref, return_map = FALSE)

test_that("test relabel", {
	expect_equal(adjusted, unname(as.numeric(map[class])))
})
