# save(in_args, out, file=paste0('lang_test_data/', id, '.RData')
source('R/1.R')
devtools::load_all()

langevin_e2e_regression <- function(){
    dir_name <- './lang_test_data/'

    for(file_name in list.files(dir_name)){
        env <- new.env()
        load(file.path(dir_name, file_name), envir= env)
        if(FALSE){with(
            env$in_args$CTMM,{
            dt <- env$in_args$dt
            if(K==2 && dt < Inf) print(c(tau[1]>tau[2], f.nu[2]*dt>0.8813736))
            }
        )}
        #print(env$in_args$CTMM$K)
        in_args <- env$in_args
        id <- in_args[['id']]
        in_args[['id']] <- NULL
        #tryCatch(
            testthat::expect_identical(
            do.call(ctmm:::langevin, in_args),
           env$out
        )
       # , error = function(e){ cat(e$message); print('\n'); print(in_args$dt)})
    }
}
langevin_e2e_regression()




file_name <- list.files(dir_name)[2]
methods::getDataPart(env$in_args$CTMM$sigma)

#tau_tmp <- env$in_args$CTMM$tau
env$in_args$CTMM$tau[1] <- tau_tmp[2]
env$in_args$CTMM$tau[2] <- tau_tmp[1]
env$in_args[['id']] <- NULL
env$out <- do.call(ctmm:::langevin, env$in_args)

id <- paste0('dim2_', id)
#id <- 'tau1inf_zwkeyfgjvmst'
env$in_args$id <- id
save(in_args, out, file=paste0('lang_test_data/', id, '.RData'), envir = env)



cov_obj <- covr::function_coverage(ctmm:::langevin, langevin_e2e_regression())
cov_obj
#82.47%
#84.54%
#89.69%

covr::zero_coverage(cov_obj)
