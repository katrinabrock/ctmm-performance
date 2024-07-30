setwd('../ctmm')
source('../testgen/testgen.R')
#source('R/1.R')
devtools::load_all()
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
args <- data$args
expected<- data$expected
ctmm_attr_used <- c('dynamics','K', 'tau', 'sigma', 'Omega2', 'f.nu', 'TfOmega2')
for(i in 1:length(args)){
    args[[i]][['CTMM']] <- args[[i]][['CTMM']][ctmm_attr_used]
}
CTMMs <- lapply(args, function(x)x$CTMM)
save(CTMMs,expected,file='tests/testthat/testdata/test_Langevin.rda')
