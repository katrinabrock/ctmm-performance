setwd('../ctmm')
source('../testgen/testgen.R')
e <- new.env()
e$fn_data <- list(args = list(), expected = list())
devtools::load_all()

data("buffalo")
Cilla <- buffalo$Cilla
m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
FITZ <- ctmm.select(Cilla,m.ouf,verbose=TRUE,cores=-1)

expect_equal(do.call(Langevin, e$fn_data[['args']][[1]]),
e$fn_data[['expected']][[1]])

## stash observer code
source("R/1.R")
devtools::load_all()

my_ctmm <- e$fn_data[['args']][[1]][['CTMM']]

sapply(1:length(e$fn_data[['args']]), function(idx) e$fn_data[['args']][[idx]][['CTMM']][['dynamics']])


td <- readRDS('Ltestdata.rds')

sample_idx <- sample(1:length(td[['args']]), 20)
for(idx in sample_idx){
expect_equal(do.call(Langevin, td[['args']][[idx]]),
td[['expected']][[idx]])
}

unique(sapply(1:length(td[['args']]), function(idx) td[['args']][[idx]][['CTMM']][['dynamics']]))
