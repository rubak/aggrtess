# Example with a single duplicated location
x <- data.frame(municipality = c("A", "B", "B"), x = c(1, 2, 2), y = c(1, 1, 1),
                households = c(1, 1, 2), individuals = c(4, 3, 8), stringsAsFactors = FALSE)

# Three ways to achieve same aggregated version:
y1 <- aggregate_duplicates(x, idcols = c("municipality", "x", "y"), sumcols = c("households", "individuals"))
y2 <- aggregate_duplicates(x, sumcols = c("households", "individuals"))
y3 <- aggregate_duplicates(x, idcols = c("municipality", "x", "y"))

# Expected result
y <- data.frame(municipality = c("A", "B"), x = c(1, 2), y = c(1, 1),
                households = c(1, 3), individuals = c(4, 11), stringsAsFactors = FALSE)
expect_identical(y, y1)
expect_identical(y, y2)
expect_identical(y, y3)

# Additional columns are dropped if not (possibly implicitly) referred to:
xx <- cbind(x, junk = 1:3)
yy <- aggregate_duplicates(xx, idcols = c("municipality", "x", "y"), sumcols = c("households", "individuals"))
expect_identical(y, yy)

# Check basic error handling:
expect_error(aggregate_duplicates(x)) # both NULL
expect_error(aggregate_duplicates(x, idcols = 1:3)) # Non-character
