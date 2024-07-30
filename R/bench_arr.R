microbenchmark::microbenchmark(
lapply(1:10, matrix),
lapply(1:10, as.matrix),
lapply(1:10, array)
)

diag(array(1:12, c(2,3)))
