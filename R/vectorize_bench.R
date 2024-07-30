source('R/bench_lib.R')
source('R/vectorize_lib.R')

v_from_o <- function(args, dt_var) return(do.call(function(...) vect_from_orig(dt_var, ...), args))


vectorized_no_factory <- function(args, dt_var){
    n <- length(dt_var)
    K <- max(args[['K']],1)
    DIM <- args[['DIM']]

    LANGEVIN <- do.call(function(...) no_factory(dt_var, ...), args)
    preGreen <- LANGEVIN$Green
    preSigma <- LANGEVIN$Sigma

    Green <- array(0,c(n,K*DIM,K*DIM))
    Sigma <- array(0,c(n,K*DIM,K*DIM))
    for(i in 1:n){
        fixed <- fix_dim(preGreen[i,,], preSigma[i,,], args[['sigma']], DIM, K) 
        Green[i,,] <- fixed$Green
        Sigma[i,,] <- fixed$Sigma
    }
    return(list(Green = Green, Sigma = Sigma))
}

vectorized <- function(args, dt_var){
    n <- length(dt_var)
    K <- max(args[['K']],1)
    DIM <- args[['DIM']]
    custom_langevin_fn <- do.call(langevin_fn_factory, args)

    LANGEVIN <- custom_langevin_fn(dt_var, dt_var < Inf, args[['nu']]*dt_var>0.8813736)
    preGreen <- LANGEVIN$Green
    preSigma <- LANGEVIN$Sigma

    Green <- array(0,c(n,K*DIM,K*DIM))
    Sigma <- array(0,c(n,K*DIM,K*DIM))
    for(i in 1:n){
        fixed <- fix_dim(preGreen[i,,], preSigma[i,,], args[['sigma']], DIM, K) 
        Green[i,,] <- fixed$Green
        Sigma[i,,] <- fixed$Sigma
    }
    return(list(Green = Green, Sigma = Sigma))
}

orig_wrapped <- function(args, dt_var){
    n <- length(dt_var)
    K <- max(args[['K']],1)
    DIM <- args[['DIM']]
    Green <- array(0,c(n,K*DIM,K*DIM))
    Sigma <- array(0,c(n,K*DIM,K*DIM))
    for(i in 1:n){
        LANGEVIN <- do.call(function(...) orig(dt_var[i], ...), args)
        Green[i,,] <- LANGEVIN$Green
        Sigma[i,,] <- LANGEVIN$Sigma
    }
    return(list(Green = Green, Sigma = Sigma)) 
}

dti <- data.table::copy(inputs_dt)

dti[,
   ctmm_id := .GRP,
   by = setdiff(names(dti), "dt")
]

dti

setorder(dti, col = 'ctmm_id')
dt_vars <- list()
argss <- list()
for(focal_ctmm_id in 1:max(dti$ctmm_id)){
  dt_sub <- dti[ctmm_id == focal_ctmm_id]
  # assumes dt is first
  args <- as.list(dt_sub[1, ])
  args[['tau']] <- na.omit(c(args[['tau_1']], args[['tau_2']]))
  for(key in c('tau_1', 'tau_2', 'dt', 'ctmm_id')) args[[key]] <- NULL 
  argss[[focal_ctmm_id]] <- args
  dt_vars[[focal_ctmm_id]] <- dt_sub$dt
}

pm_Langevin <- function(fn){
    r <- list()
    for(focal_ctmm_id in 1:max(dti$ctmm_id)){
      r[[focal_ctmm_id]] <- fn(argss[[focal_ctmm_id]], dt_vars[[focal_ctmm_id]])
    }
    return(r)
}


TASK = 'benchmark'
if(TASK == 'profile'){
    Sys.setenv(RSTUDIO_PANDOC=pandoc::pandoc_locate())
    my_profile <- profvis::profvis({
        #lapply(1:1e3, function(x) pm_Langevin(orig_wrapped))
        #lapply(1:1e3, function(x) pm_Langevin(vectorized))
        lapply(1:1e4, function(x) pm_Langevin(vectorized_no_factory))
    })
    htmlwidgets::saveWidget(my_profile, "../ctmm/profile_vectorized.html")
}

if(TASK == 'benchmark'){
    result <- microbenchmark::microbenchmark(
        pm_Langevin(orig_wrapped),
        pm_Langevin(vectorized),
        pm_Langevin(vectorized_no_factory),
        pm_Langevin(v_from_o),
        check = 'identical'
    )
    print(result)
}

if(TASK == 'check'){
    options(error=recover)
    expected <- pm_Langevin(orig_wrapped)
    actual <- pm_Langevin(v_from_o)

    options(error=NULL)
    testthat::expect_identical(
        actual,
        expected
    )
}