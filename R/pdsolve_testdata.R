setwd('../ctmm')
source('../testgen/testgen.R')
#source('../ctmm/R/1.R')
devtools::load_all('../ctmm')
data <- generate_tests(
    ctmm:::pd.solve,
    {
        data("buffalo")
        Cilla <- buffalo$Cilla
        if(FALSE){
        m.iid <- ctmm(sigma=23 %#% "km^2")
        M.IID <- ctmm.fit(Cilla,m.iid)
        } else {
        m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
        FITZ <- ctmm.select(Cilla,m.ouf,verbose=TRUE,cores=-1)
        }
    }
)
args <- data$args
expected<- data$expected
save(args,expected,file='tests/testthat/testdata/test_pd_solve.rda')
