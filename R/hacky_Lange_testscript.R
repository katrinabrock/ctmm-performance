
td <- readRDS('Ltestdata.rds')

sample_idx <- 1:length(td[['args']]
#sample_idx <- sample(1:length(td[['args']]), 20) 

source("R/1.R")
devtools::load_all()

for(idx in sample_idx){
expect_equal(do.call(Langevin, td[['args']][[idx]]),
td[['expected']][[idx]])
}
