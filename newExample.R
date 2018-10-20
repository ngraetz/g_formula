rm(list=ls())

library(tidyverse)
library(nnet)
library(parallel)
source("./gfit.R")

binAges <- function(x){
    cut(
        x, 
        breaks = c(15, 22, 29, Inf), 
        labels = paste0(c(15, 23, 30), "to", c(22, 29, Inf)))
}
# Clean the shit out of some of this data
# Since this is simulated data this should have been done pre script 
# but we are flexin on them with the tidyverse
DF <- read_csv("./ex.csv") %>%
    # remove the fat
    select(-l.partner.single, -l.partner.cohab, -l.partner.married) %>%
    select(-X1, -datapart, -fut, -region.A) %>%
    # turn the dummy variables into normal variables
    gather(partner, Count, `partner.single`:`partner.married`) %>%
    filter(Count >=1) %>%
    select(-Count) %>%
    select(-l.job.edu, -l.job.other, -l.job.full) %>%
    gather(job, Count, job.edu, job.full, job.other) %>%
    filter(Count >=1) %>%
    select(-Count) %>%
    # Clean up the new variables
    mutate(job=gsub("job\\.", "", job)) %>%
    mutate(partner=gsub("partner\\.", "", partner)) %>%
    mutate(agegroup=binAges(age)) %>%
    # do some more cleaning
    select(-(age1522:age30pl), -(age1519:age36pl)) %>%
    arrange(id, year) %>%
    group_by(id) %>%
    mutate(l.partner=lag(partner), l.job=lag(job)) %>%
    ungroup

# build up the model structure base
baseForm <- ~ parental.income + region.B +# time constants
    # time varying #
    (agegroup) * (l.birth + l.totalbirth + l.partner + l.job + l.totaledu)

# formula for models
formulas <- list(
    update(baseForm, birth ~ .),
    update(baseForm, job ~ .),
    update(baseForm, partner ~ .),
    censor ~ parental.income +
        region.B + 
        # time varying #
        (agegroup) * (birth + totalbirth + partner + job + totaledu)
)

families <- list(binomial, NULL, NULL, binomial) # families for models
functions <- list(glm, multinom, multinom, glm) # functions for results

# here are our lags that we want to keep updated
lags <- list(
    l.birth = function(DF) DF %>% mutate(l.birth=lag(birth)),
    l.totalbirth = function(DF) DF %>% mutate(l.totalbirth=lag(totalbirth)),
    l.job = function(DF) DF %>% mutate(l.job=lag(job)),
    l.totaledu = function(DF) DF %>% mutate(l.totaledu=lag(totaledu)),
    l.partner = function(DF) DF %>% mutate(l.partner=lag(partner))
)

# here are the deterministic and probabilistic rules
rules <- list(
    birth = function(DF, models) simPredict(DF, models, 1),
    job = function(DF, models) simPredict(DF, models, 2),
    partner = function(DF, models) simPredict(DF, models, 3),
    totaledu = function(DF, ...) DF$totaledu + (DF$job == "edu"),
    totalbirth = function(DF, ...) DF$totalbirth + DF$birth,
    censor = function(DF, models) simPredict(DF, models, 4)
)

boots <- 6 # number of bootstraps
replicationSize <- 6 # number of monte carlo replication

bootruns <- mclapply(1:boots, function(b){
    # sample individuals with replacement not rows
    sampleDF <- DF %>%
        select(id) %>%
        unique %>%
        sample_frac(replace=TRUE) %>%
        right_join(DF)
    
    # Fit the model
    gfitboot <- gfit.init(formulas, families, functions, data=sampleDF)
    
    # Run the "natural course" rules
    mcDF <- bind_rows(lapply(1:replicationSize, function(i) sampleDF))
    naturalDF <- progressSimulation(mcDF, lags, rules, gfitboot)
    list(natural=naturalDF, models=gfitboot)
}, mc.cores=6)
