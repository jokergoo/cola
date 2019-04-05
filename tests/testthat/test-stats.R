
m = matrix(1:18, nr = 9)
fa = c("a", "a", "a", "b", "b", "b", "c", "c", "c")
d1 = mean_group_dist(m, fa)

d2 = as.matrix(dist(m))
d2 = mean(c(mean(d2[1:3, 4:6]), mean(d2[1:3, 7:9]), mean(d2[4:6, 7:9])))

test_that("test mean_group_dist", {
	expect_equal(d1, d2)
})

fa = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
d1 = mean_group_dist(m, fa)
test_that("test mean_group_dist", {
	expect_equal(d1, 0)
})

fa = c("a", "a", "a", "b", "b", "b", "b", "b", "c")
d1 = mean_group_dist(m, fa)

d2 = as.matrix(dist(m))
d2 = mean(c(mean(d2[1:3, 4:8]), mean(d2[1:3, 9]), mean(d2[4:8, 9])))

test_that("test mean_group_dist", {
	expect_equal(d1, d2)
})
