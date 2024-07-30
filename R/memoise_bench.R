### Setup ###
system('git checkout master')
devtools::load_all()
data("buffalo")

Cilla <- buffalo$Cilla

m.iid <- ctmm(sigma=23 %#% "km^2")
m.ou <- ctmm(sigma=23 %#% "km^2",tau=6 %#% "day")
m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))

TEST <- FALSE
for(i in 1:40){
    branch <- c('master', 'memoise')[(i %% 2) + 1]
    system(paste0('git checkout ', branch))
    source('R/1.R')
    devtools::load_all()
    if(TEST) ctmm.fit <- function(...){rnorm(10000)}
    Rprof(paste(branch, i, 'Rprof.out', sep='.'))
    M.IID <- ctmm.fit(Cilla,m.iid)
    M.OU <- ctmm.fit(Cilla,m.ou)
    M.OUF <- ctmm.fit(Cilla,m.ouf)
    Rprof(NULL)
    memoise::forget('ctmm:::langevin')
}

# Analyse 
library(data.table)

n <- 40
results <- data.table(
    idx = rep(0, n),
    branch = rep('', n),
    run_total_time = rep(0, n),
    pct_langevin = rep(0,n)
)
for(i in 1:n){
    branch <- c('master', 'memoise')[(i %% 2) + 1]
    filename <- paste(branch, i, 'Rprof.out', sep='.')
    summary <- summaryRprof(filename)
    results[i, "idx"] <- i
    results[i, "branch"] <- branch
    results[i, "run_total_time"] <- summary$sampling.time
    results[i, "pct_langevin"] <- summary$by.total['"langevin"', 'total.pct']
}

png('total_time.png')
boxplot(results$run_total_time ~ results$branch)
dev.off()

png('langevin_percent.png')
boxplot(results$pct_langevin ~ results$branch)
dev.off()
