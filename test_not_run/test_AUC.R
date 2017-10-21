
set.seed(12345)
library(matrixStats)

nr1 = 100
mat1 = matrix(rnorm(100*nr1), nrow = nr1)

nr2 = 10
library(mvtnorm)
sigma = matrix(0.8, nrow = nr2, ncol = nr2); diag(sigma) = 1
mat2 = t(rmvnorm(100, mean = rep(0, nr2), sigma = sigma))

nr3 = 50
library(mvtnorm)
sigma = matrix(0.5, nrow = nr3, ncol = nr3); diag(sigma) = 1
mat3 = t(rmvnorm(100, mean = rep(0, nr3), sigma = sigma))


mat = t(rbind(mat1, mat2, mat3))

par(mfrow = c(2, 2))
par(mar = c(3, 5, 1, 1), cex = 0.8)
v = NULL
min_cor = 0
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "P(X <= x)", xlab = "")
for(i in 1:ncol(mat)) {
    cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE]))
    cor_v = sort(cor_v)
    f = ecdf(cor_v)
    # cor_v = c(0, cor_v, 1)
    # n2 = length(cor_v)
    # v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
    cor_v = seq(min_cor, 1, length = 100)
    n2 = length(cor_v)
    v[i] = sum((cor_v[2:n2] - cor_v[1:(n2-1)])*f(cor_v[-n2]))
    x = seq(0, 1, length = 100)
    lines(x, f(x), col = ifelse(i <= nr1, "#00000080", ifelse(i <= nr1+nr2, "#FF000080", "#00FF0080")))
}
# abline(v = 0.2, lty = 2, col = "grey")
legend("bottomright", lty = 1, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))
v = (1 - min_cor) - v

dx_list = NULL
dy_list = NULL
for(i in 1:ncol(mat)) {
    cor_v <- abs(cor(mat[, i, drop = FALSE], mat[, -i, drop = FALSE]))
    d = density(cor_v)
    dx_list = c(dx_list, list(d$x))
    dy_list = c(dy_list, list(d$y))
}
max_y = max(unlist(dy_list))
plot(NULL, xlim = c(0, 1), ylim = c(0, max_y), ylab = "density", xlab = "")
for(i in seq_along(dx_list)) {
	lines(dx_list[[i]], dy_list[[i]], col = ifelse(i <= nr1, "#00000080", ifelse(i <= nr1+nr2, "#FF000080", "#00FF0080")))
}
legend("topright", lty = 1, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))


plot(v, pch = 16, xlab = "", ylab = "AAC (|cor| in c(0, 1))", col = ifelse(seq_along(v) <= nr1, "#000000", ifelse(seq_along(v) <= nr1+nr2, "#FF0000", "#00FF00")))
legend("topleft", pch = 16, col = 1:3, legend = c("nr = 100, cor = 0", "nr = 10, cor = 0.8", "nr = 50, cor = 0.5"))


plot(colSds(mat), pch = 16, xlab = "", ylab = "SD", col = ifelse(seq_along(v) <= nr1, "#000000", ifelse(seq_along(v) <= nr1+nr2, "#FF0000", "#00FF00")))



cor_v <- abs(cor(mat[, 160, drop = FALSE], mat[, -160, drop = FALSE]))
cor_v = sort(cor_v)
f = ecdf(cor_v)
x = seq(0, 1, length = 100)
y = f(x)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "P(X < x)")
x2 = seq(0, 1, length = 100)
polygon(c(x2, rev(x2)), c(f(x2), rep(1, length(x2))), col = "#FF000040", border = NA)
polygon(c(x2, rev(x2)), c(f(x2), rep(0, length(x2))), col = "#00FF0040", border = NA)
lines(x, y)
text(0.4, 0.9, "AAC")
text(0.6, 0.4, "AUC")

par(mfrow = c(1, 2))
res = get_single_run(res_list, "AAC", "kmeans")
m = get_consensus(res, k = 3)
f = ecdf(m)
plot(x, f(x), type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "P(X > x)")
for(i in 1:20) {
    x1 = runif(1, min = 0, max = 0.3)
    x2 = runif(1, min = 0.7, max = 1)
    x0 = seq(x1, x2, length = 100)
    polygon(c(x1, x2, rev(x0)), c(f(x1), f(x1), f(rev(x0))), col = "#FF000010", border = NA)
}
title("high PAC")

res = get_single_run(res_list, "AAC", "skmeans")
m = get_consensus(res, k = 3)
f = ecdf(m)
plot(x, f(x), type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "P(X > x)")
for(i in 1:20) {
    x1 = runif(1, min = 0, max = 0.3)
    x2 = runif(1, min = 0.7, max = 1)
    x0 = seq(x1, x2, length = 100)
    polygon(c(x1, x2, rev(x0)), c(f(x1), f(x1), f(rev(x0))), col = "#FF000010", border = NA)
}
title("low PAC")
