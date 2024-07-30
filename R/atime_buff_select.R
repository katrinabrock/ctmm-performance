N_min <- 100
N_max <- 3527
N_length <- 3
times <- 1
seconds.limit <- 600

print(glue::glue('Max run time is {N_length*times*seconds.limit/3600} hours'))

tdir <- tempfile()
dir.create(tdir)
git2r::clone("https://github.com/katrinabrock/ctmm", tdir)
atime.list <- atime::atime_versions(
 pkg.path=tdir,
 seconds.limit = seconds.limit,
 times=times,
 N=round(exp(seq(log(N_min), log(N_max), length.out=N_length))),
 setup={
    data("buffalo")
    Cilla <- buffalo$Cilla[1:N, ]
 },
 expr=ctmm::ctmm.select(Cilla,ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour")),verbose=TRUE,cores=1),
 pre_speedup = "66d1f5180323ae5502417c908e4e686ba65d2334",
 post_speedup="37a1aa480ab326f524e4314731cbfb2107182180",
 memoise_no_speedup="cd390665b6076b113d27abad06f1f778f3bc66fe",
 memoise_with_speedup="a50c47178037483c76f521b60329395d1c538e00")
saveRDS(atime.list, '20240717_atime.rds')

atime.list <- readRDS('./data/20240717_atime.rds')
plot(atime.list)
