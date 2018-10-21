rm(list=ls())

## Options needed for this run.
repo <- 'C:/Users/Nick/Documents/repos/g_formula/'
local_cores <- 4

## Set up packages and functions.
setwd(repo)
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
# If we want to do scenarios we would chnage these rules or the starting values.
rules <- list(
    agegroup = function(DF, ...) binAges(DF$age),
    birth = function(DF, models) simPredict(DF, models, 1),
    job = function(DF, models) simPredict(DF, models, 2),
    partner = function(DF, models) simPredict(DF, models, 3),
    totaledu = function(DF, ...) DF$totaledu + (DF$job == "edu"),
    totalbirth = function(DF, ...) DF$totalbirth + DF$birth,
    censor = function(DF, models) simPredict(DF, models, 4)
)

boots <- 15 # number of bootstraps 100 takes a while
replicationSize <- 6 # number of monte carlo replication
## I feel like this will be a time hog and we can parallelize it instead of rbinding... but maybe it doesn't slow it down that much.

set.seed(123)

if(!file.exists("./bootruns.Rds")){
  
    bootruns <- mclapply(1:boots, function(b) {
      
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
        list(natural=naturalDF)
        
    }, mc.cores=local_cores)
    
    saveRDS(bootruns, "./bootruns.Rds")
}

bootruns <- read_rds("./bootruns.Rds")

pSimsDF <- bind_rows(lapply(1:length(bootruns), function(i){
    bootruns[[i]]$natural %>%
        mutate(sim=i)})) %>%
    group_by(age, sim) %>%
    summarize(
        nbirths = sum(birth), 
        nmoms = n(),
        pedu = mean(job == "edu"),
        pother = mean(job == "other"),
        pfull = mean(job == "full"),
        psingle = mean(partner == "single"),
        pmarried = mean(partner == "married"),
        pcohab = mean(partner == "cohab")) %>%
    mutate(pbirths=nbirths/nmoms) %>%
    select(-nbirths, -nmoms)

actualDF <- DF %>%
    group_by(age) %>%
    summarize(
        nbirths = sum(birth), 
        nmoms = n(),
        pedu = mean(job == "edu"),
        pother = mean(job == "other"),
        pfull = mean(job == "full"),
        psingle = mean(partner == "single"),
        pmarried = mean(partner == "married"),
        pcohab = mean(partner == "cohab")) %>%
    mutate(pbirths=nbirths/nmoms) %>%
    select(-nbirths, -nmoms) %>%
    gather("measure", "mean", -age) %>%
    mutate(type="Observed")


# I dont feel great about these confidence intervals we should talk about this
# Also this hould be automated
pSimsDF %>%
    select(-sim) %>%
    summarise_all(
        list(
            mean = mean, 
            lwr = function(x) quantile(x, probs=.025),
            upr = function(x) quantile(x, probs=.975))) %>%
    ungroup %>%
    gather("Metric", "Value", -age) %>%
    mutate(measure=gsub("_[A-z ]*", "", Metric)) %>%
    mutate(statistic=gsub("[A-z ]*_", "", Metric)) %>%
    select(-Metric) %>%
    spread("statistic", "Value") %>%
    mutate(type="Estimate") %>%
    bind_rows(actualDF) %>%
    filter(measure != "pother") %>%
    ggplot(aes(x=age, y=mean, group=type, color=type, fill=type)) +
    geom_line(linetype=3) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.2) +
    theme_classic() +
    facet_wrap(~measure)

## ESTIMATING EFFECTS: calculate controlled direct effect (CDE) of education by editing rules list and rerunning.
## Key change: we need a new rule for when we want to update a variable by pulling from the natural course DF. So this requires
## having estimated the natural course already. Instead of simPredict(), we use simScenario() for everything except the outcome. We still 
## simulate probabilistically the outcome (and the censoring...?). 
##    simScenario() = draw stochastically from the natural course distribution (or intervention distribution) of the variable at time t.
##
## Notes: I'm confused now though... maybe we always need a reference scenario DF to compare to the natural course DF. Like Maarten intervenes on 
## partner in the input DF, simulates a new scenario DF, and then calculates effects by simulating and drawing from the intervention scenario
## for partner and drawing natural course for job. Soooo I guess if you want to know the effect of a variable, you need to frame it as the effect
## of some difference (like having college education instead of high school) so that you have a reference set. So I guess our two rules would be:
##    simPredict() = used to update all probabilistic things in the natural course (or any scenario where we edit the input DF).
##    simScenario() = used to update all probabilitic things in effect calculations by drawing stochastically from the target variable distribution 
##                    at time t for whatever DF you give it.
##
## So I thiiiiink the way we want to do this is always hand all rules the DF to update, the natural course DF, and the intervention DF. And then
## in the function the rule uses to actually update values, we manipulate which DF it should use based on what direct/indirect effect we are
## estimating.
rules <- list(
  agegroup = function(DF, ...) binAges(DF$age), # Still deterministic.
  birth = function(DF, natural_DF, intervention_DF, models, ...) simPredict(DF, models, 1), # Still probabilistic based on model.
  job = function(DF, natural_DF, intervention_DF, models) simScenario(natural_DF, models, 2), # Draw stochastically from natural course.
  partner = function(DF, natural_DF, intervention_DF, models) simScenario(intervention_DF, models, 3), # Draw stochastically from intervention.
  totaledu = function(DF, ...) DF$totaledu + (DF$job == "edu"), # Still deterministic.
  totalbirth = function(DF, ...) DF$totalbirth + DF$birth, # Still deterministic.
  censor = function(DF, natural_DF, intervention_DF, models, ...) simPredict(DF, models, 4) # Still probabilistic based on model (?).
)

