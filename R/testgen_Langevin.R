setwd('../ctmm')
source('../testgen/testgen.R')
#source('R/1.R')
devtools::load_all()

if(FALSE){
    for(idx in 1:length(args)){
        testthat::expect_equal(
            do.call(Langevin, args[[idx]]),
            expected[[idx]]
        )
    }
}


data <- generate_tests(
    ctmm:::Langevin,
    {
        data("buffalo")
        Cilla <- buffalo$Cilla
        if(FALSE){
        m.iid <- ctmm(sigma=23 %#% "km^2")
        M.IID <- ctmm.fit(Cilla,m.iid)
        }
        m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
        FITZ <- ctmm.select(Cilla,m.ouf,verbose=TRUE,cores=-1)
    }
)
stop()

#write_test_file(Langevin, data, file = '20240725_test_Langevin.R')

perms <- list()
name_collection <- unique(
    do.call(c,
      lapply(data$args, function(x) names(x$CTMM))
    )
)

unique(do.call(c, lapply(data$args, function(x)x$DIM)))
unique(do.call(c, lapply(data$args, function(x)x$CTMM$dynamics)))
unique(do.call(c, lapply(data$args, function(x)x$CTMM$timelink)))

length(unique(lapply(data$args, function(x)x$CTMM$dt)))

length(unique(lapply(data$args, function(x) x$t )))

unique(lapply(data$expected, function(x) dim(x$Green) ))

args <- data$args
expected<- data$expected

#t <- args[[1]]$t
#dt <- data$args[[1]]$CTMM$dt

#unique(dt == c(Inf, diff(t)))

# DIM 1 stationary
ctmm_attr_used <- c('dynamics','K', 'tau', 'sigma', 'Omega2', 'f.nu', 'TfOmega2')
for(i in 1:length(args)){
    args[[i]][['CTMM']] <- args[[i]][['CTMM']][ctmm_attr_used]
}

print(object.size(data$args), units='Kb')
print(object.size(args), units='Kb')
print(object.size(CTMMs), units='Kb')
args[[1]][['CTMM']][ctmm_attr_used[2:length(ctmm_attr_used)]]

Langevin_test_CTMMs <- lapply(args, function(x)x$CTMM)
usethis::use_data(Langevin_test_CTMMs, internal=TRUE, overwrite=TRUE)

write(put(CTMMs), 'tmp.R')

t <- buffalo$Cilla$t/3600

n <- length(t)
n_cases <- length(args)

Green <- array(dim = c(n_cases, n, 2, 2))
Sigma <- array(dim = c(n_cases, n, 2, 2))

for(i in 1:n_cases){
    k <- dim(expected[[i]]$Green)[2]
    #print(c(dim(Green[i,,1:k,1:k]),dim(expected[[i]]$Green)))
    Green[i,,1:k,1:k] <-  expected[[i]]$Green   
    Sigma[i,,1:k,1:k] <-  expected[[i]]$Sigma   
}

print(object.size(data$expected), units='Kb')
print(object.size(list(Green=Green,Sigma=Sigma)), units='Kb')

Langevin_expected<- data$expected
usethis::use_data(Langevin_expected, internal=TRUE, overwrite=TRUE)

dt <- c(Inf, diff(t))
dti <- sort(dt,index.return=TRUE)
dtl <- unique(dti$x) # dt levels
dti <- dti$ix # sort indices


for(idx in 1:length(args)){
    testthat::expect_identical(
        Langevin(
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
         ),
         expected = expected[[idx]]

          
    )
}
if(FALSE){
    #k <- max(CTMMs[[idx]]$K,1)
list(
            Green = array(Green[idx,,1:k,1:k], c(n,k,k)),
            Sigma = array(Sigma[idx,,1:k,1:k], c(n,k,k))
         )
}