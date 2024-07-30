setwd('../ctmm')
source('../testgen/testgen.R')
#source('R/1.R')
devtools::load_all()
data <- generate_tests(
    ctmm:::Langevin,
    {
        data("wolf")
        Gamba <- wolf$Gamba
        PROTO <- ctmm(mean="periodic",period=c(24 %#% "hours",1 %#% "month"),circle=TRUE)
        SVF <- variogram(Gamba,res=3)
        GUESS <- ctmm.guess(Gamba,PROTO,variogram=SVF,interactive=FALSE)
        # CRAN policy limits to 2 processes (cores)
        FITS <- ctmm.select(Gamba,GUESS,verbose=TRUE,cores=-1)
    }
)
args <- data$args
expected<- data$expected
ctmm_attr_used <- c('dynamics','K', 'tau', 'sigma', 'Omega2', 'f.nu', 'TfOmega2')
for(i in 1:length(args)){
    args[[i]][['CTMM']] <- args[[i]][['CTMM']][ctmm_attr_used]
}
CTMMs <- lapply(args, function(x)x$CTMM)
save(CTMMs,expected,file='tests/testthat/testdata/test_Langevin_wolf.rda')
