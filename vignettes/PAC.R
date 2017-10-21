
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))

x = seq(0, 1, length = 100)
m = mat_list[[1]]

f = ecdf(m)
plot(x, f(x), type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "P(X < x)")
for(i in 1:20) {
    x1 = runif(1, min = 0, max = 0.3)
    x2 = runif(1, min = 0.7, max = 1)
    x0 = seq(x1, x2, length = 100)
    polygon(c(x1, x2, rev(x0)), c(f(x1), f(x1), f(rev(x0))), col = "#FF000010", border = NA)
}

hc = hclust(dist(m))
od = hc$order
m = m[rev(od), od]
image(m, col = colorRampPalette(c("white", "blue"))(40), axes = FALSE, ann = FALSE)

m = mat_list[[2]]

f = ecdf(m)
plot(x, f(x), type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "P(X < x)")
for(i in 1:20) {
    x1 = runif(1, min = 0, max = 0.3)
    x2 = runif(1, min = 0.7, max = 1)
    x0 = seq(x1, x2, length = 100)
    polygon(c(x1, x2, rev(x0)), c(f(x1), f(x1), f(rev(x0))), col = "#FF000010", border = NA)
}


hc = hclust(dist(m))
od = hc$order
m = m[rev(od), od]
image(m, col = colorRampPalette(c("white", "blue"))(40), axes = FALSE, ann = FALSE)
