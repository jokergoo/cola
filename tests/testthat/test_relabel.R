

class = c(rep("a", 10), rep("b", 3))
ref = c(rep("b", 4), rep("a", 9))
map = cola:::relabel_class(class, ref)
attr(map, "df") = NULL

test_that("test relabel", {
	expect_equal(map, c("a" = "b", "b" = "a"))
})
