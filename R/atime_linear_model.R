m <- atime.list$measurements

#https://stackoverflow.com/questions/43152845/expand-list-column-of-data-tables
m[, r:= as.character(.I)]
res <- m[,rbindlist(setNames(lapply(time, list), r), id="r")]

m <- res[m, on=.(r)]

m$time_obs <- m$V1
m$V1 <- NULL


summary(lm(time_obs~N + expr.name, data = m))
model.matrix(time_obs~N + expr.name, data = m)
m[,.(expr.name, .I)][expr.name == 'brock_speedup']
