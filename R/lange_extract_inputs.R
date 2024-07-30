library(ctmm)
dir_name <- './data/lang_test_data/'


my_dt <- data.table::rbindlist(lapply(list.files(dir_name), function(file_name) {
  env <- new.env()
  load(file.path(dir_name, file_name), envir= env)
  print(length(env$in_args$CTMM$tau))
  with(env$in_args,
      return(list(
        dt = dt,
        K = CTMM$K,
        tau_1 = CTMM$tau[1],
        tau_2 = CTMM$tau[2],
        sigma = methods::getDataPart(CTMM$sigma),
        Omega2 = CTMM$Omega2,
        f = CTMM$f.nu[1], # mean(f)
        nu = CTMM$f.nu[2], # nu || omega
        TT = CTMM$TfOmega2, # 2 f / Omega^2
        DIM = DIM
      ))
    )
}))

write.csv(my_dt, 'data/langevin_in.csv', row.names = FALSE)
