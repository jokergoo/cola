
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


ATC_cgi_anno = function (mat, cgi_anno, min_cor = 0, mc.cores = 1) {
    mat = t(mat)
    n = ncol(mat)

    ind_list = split(seq_len(n), cgi_anno)

    v_list = mclapply(ind_list, function(ind) {
        v = numeric(length(ind))
        for (i in seq_along(ind)) {
            ind2 = ind[-i]
            if (length(ind2) > 1000) {
                ind2 = sample(ind2, 1000)
            }
            suppressWarnings(cor_v <- abs(cor(mat[, ind[i], drop = FALSE], mat[, ind2, drop = FALSE])))
            if(sum(is.na(cor_v))/length(cor_v) >= 0.75) {
                v[i] = 1
            } else {
                f = ecdf(cor_v)
                cor_v = seq(min_cor, 1, length = 1000)
                n2 = length(cor_v)
                v[i] = sum((cor_v[2:n2] - cor_v[1:(n2 - 1)]) * f(cor_v[-n2]))
            }
        }
        return(v)
    }, mc.cores = mc.cores)

    v = numeric(n)
    for(i in seq_along(v_list)) {
        v[ ind_list[[i]] ] = v_list[[i]]
    }
    v = 1 - min_cor - v
    names(v) = NULL
    return(v)
}

m = matrix(rnorm(100), 10)
group = c(rep("a", 5), rep("b", 5))
s1 = ATC(m, group = group, min_cor = 0.5)
s2 = ATC_cgi_anno(m, cgi_anno = group, min_cor = 0.5)

test_that("test ATC with groups", {
	expect_equal(s1, s2)
})
