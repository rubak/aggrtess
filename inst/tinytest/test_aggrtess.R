# Make a subset of `frb` with `n_p` marked points
n_p <- 1000
X <- frb[floor(seq(1, npoints(frb), length.out = n_p))]
# Choose `n_c` centers for aggregation tess
n_c <- 10
cent <- X[floor(seq(1, npoints(X), length.out = n_c))]

a <- aggrtess(X, cent, target = c(50, 100))
expect_true(inherits(a, "aggrtess"))
expect_equal(a$tessel$n, n_c)

nam <- c("id", names(marks(X)), "min")
dev_all <- target_deviation(a, which = "all")
dev_bad <- target_deviation(a, which = "bad")
dev_worst <- target_deviation(a, which = "worst")
expect_true(inherits(dev_worst, "data.frame"))
expect_equal(names(dev_worst), nam)
expect_true(nrow(dev_worst) == 1)
expect_true(nrow(dev_all) == n_c)
expect_true(1 <= nrow(dev_bad) & nrow(dev_bad) <= n_c)
