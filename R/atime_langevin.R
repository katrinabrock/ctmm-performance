N_length <- 3
times <- 3
seconds.limit <- 30  

print(glue::glue('Max run time is {N_length*times*seconds.limit/60} minutes'))

N_min <- 100
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
    library(ctmm)
    data('buffalo')
    load(file.path('.', 'tests/testthat/testdata', 'test_Langevin.rda'))
    t <- (buffalo$Cilla$t/3600)[1:N]
    dt <- c(Inf, diff(t))
    dti <- sort(dt,index.return=TRUE)
    dtl <- unique(dti$x) # dt levels
    dti <- dti$ix # sort indices
 },
 expr={
    for(idx in 1:length(CTMMs)){
      ctmm:::Langevin(
          t = t,
          CTMM = c(
              CTMMs[[idx]],
              list(
                  dt = dt,
                  dti = dti,
                  dtl = dtl
              )
          ),
          DIM = 1
       )
    }
 },
 pre_speedup = "66d1f5180323ae5502417c908e4e686ba65d2334",
 chris_speedup="37a1aa480ab326f524e4314731cbfb2107182180",
 brock_speedup="caab83541300f475ae3b872862c31670b7e7a828",
 vectorized="58b7152c20367169d3d158a38e2e479292e7a787"
)
saveRDS(atime.list, './data/20240725_atime_tmp.rds')

library(atime)
atime.list <- readRDS('./data/20240725_atime_tmp.rds')
png()
plot(atime.list)
dev.off()
