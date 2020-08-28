
m1 = matrix(rnorm(10*2), nrow = 10)
m2 = matrix(rnorm(10*2), nrow = 10)

d1 = cola:::pdist(t(m1), t(m2))

d2 = as.matrix(dist(t(cbind(m1, m2))))[1:2, 3:4]
dimnames(d2) = NULL

test_that("test pdist", {
	expect_equal(all(abs(d1 - d2) < 1e-6), TRUE)
})

