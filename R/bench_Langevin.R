#devtools::load_all()
library(data.table)
source('../ctmm-performance/R/vectorizedLangevin.R')
source('../ctmm-performance/R/originalLangevin.R')

load_constructed_ctmms <- function(){
    inputs_dt <- fread('../ctmm-performance/data/langevin_in.csv')
    dt_in <- data.table::copy(inputs_dt)
    dt_in[,
       ctmm_id := .GRP,
       by = setdiff(names(dt_in), "dt")
    ]
    setorder(dt_in, col = 'ctmm_id')

    CTMMs <- list()
    for(focal_ctmm_id in 1:max(dt_in$ctmm_id)){
      dt_sub <- dt_in[ctmm_id == focal_ctmm_id]
      dt <- dt_sub$dt
      dti <- sort(dt,index.return=TRUE)
      dtl <- unique(dti$x) # dt levels
      dti <- dti$ix # sort indices
      K <- min(dt_sub[1, K], 1)
      CTMM <- list(
        dynamics = 'stationary',
        K = K,
        tau = dt_sub[1, c(tau_1, tau_2)][1:K],
        sigma = dt_sub[1, sigma],
        Omega2 = dt_sub[1, Omega2],
        f.nu = dt_sub[1, c(f, nu)],
        TfOmega2 = dt_sub[1, TT],
        dt = dt,
        dti = dti,
        dtl = dtl
      )
      # assumes dt is first
      CTMMs[[focal_ctmm_id]] <- CTMM
    }
    return(CTMMs)
}


load_saved_ctmms <- function(){
    data('buffalo')
    load(file.path('.', 'tests/testthat/testdata', 'test_Langevin.rda'))
    t <- (buffalo$Cilla$t/3600)
    dt <- c(Inf, diff(t))
    dti <- sort(dt,index.return=TRUE)
    dtl <- unique(dti$x) # dt levels
    dti <- dti$ix # sort indices
    return(lapply(CTMMs, function(x) c(x, list( dt = dt, dti = dti, dtl = dtl))))
}

runner <- function(fn, CTMMs){
    for(CTMM in CTMMs){
      fn(
          t = 1,
          CTMM = CTMM,
          DIM = 1
       )
    }
}

saved_ctmms <- load_saved_ctmms()
constructed_ctmms <- load_constructed_ctmms()

CTMMs <- constructed_ctmms
CTMMs <- saved_ctmms

TASK = 'profile'
if(TASK == 'profile'){
    Sys.setenv(RSTUDIO_PANDOC=pandoc::pandoc_locate())
    my_profile <- profvis::profvis({
        runner(o_Langevin, CTMMs)
        runner(v_Langevin, CTMMs)
    })
    htmlwidgets::saveWidget(my_profile, "../ctmm/profile_Langevin.html")
}

if(TASK == 'benchmark'){
    result <- microbenchmark::microbenchmark(
        runner(o_Langevin, CTMMs),
        runner(v_Langevin, CTMMs),
        check = 'identical'
    )
    print(result)
}