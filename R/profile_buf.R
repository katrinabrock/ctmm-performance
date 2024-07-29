profile_filename <- "profile_postoptim.html"


library(profvis)
Sys.setenv(RSTUDIO_PANDOC=pandoc::pandoc_locate())
## Prep
#library(ctmm)
#source("R/1.R")
devtools::load_all()
data("buffalo")

Cilla <- buffalo$Cilla

m.iid <- ctmm(sigma=23 %#% "km^2")
m.ou <- ctmm(sigma=23 %#% "km^2",tau=6 %#% "day")
m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))

my_profile <- profvis({
    M.IID <- ctmm.fit(Cilla,m.iid)
    M.OU <- ctmm.fit(Cilla,m.ou)
    M.OUF <- ctmm.fit(Cilla,m.ouf)
})
htmlwidgets::saveWidget(my_profile, profile_filename)

Rprof(line.profiling = TRUE)
for(i in 1:10){ctmm.fit(Cilla,m.iid)}
Rprof(NULL)