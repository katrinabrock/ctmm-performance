N_length <- 8
times <- 8
seconds.limit <- 150

print(glue::glue('Max run time is {N_length*times*seconds.limit/3600} hours'))

N_min <- 1000
N_max <- 3527
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
    m.iid <- ctmm::ctmm(sigma=23 %#% "km^2")
    m.ou <- ctmm::ctmm(sigma=23 %#% "km^2",tau=6 %#% "day")
    m.ouf <- ctmm::ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
 },
 expr={
    M.IID <- ctmm::ctmm.fit(Cilla,m.iid)
    M.OU <- ctmm::ctmm.fit(Cilla,m.ou)
    M.OUF <- ctmm::ctmm.fit(Cilla,m.ouf)
 },
 #pre_speedup = "66d1f5180323ae5502417c908e4e686ba65d2334",
 #chris_speedup="37a1aa480ab326f524e4314731cbfb2107182180",
 brock_speedup="334317e087cdc9b25d1e8486a877fd46e14451ac",
 vectorized="8e0d14e74549a5655a1823b186e9ab4899940de7"
)

saveRDS(atime.list, './data/20240729_atime.rds')

if(FALSE){
library(atime)
atime.list <- readRDS('./data/20240729_atime.rds')
plot(atime.list)
}
