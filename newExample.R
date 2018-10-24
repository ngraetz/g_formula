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

# Here are the natural deterministic and probabilistic rules.
rules <- list(
    agegroup = function(DF, ...) binAges(DF$age),
    birth = function(DF, models) simPredict(DF, models, 1),
    job = function(DF, models) simPredict(DF, models, 2),
    partner = function(DF, models) simPredict(DF, models, 3),
    totaledu = function(DF, ...) DF$totaledu + (DF$job == "edu"),
    totalbirth = function(DF, ...) DF$totalbirth + DF$birth,
    censor = function(DF, models) simPredict(DF, models, 4)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
# Let's test replacing any partner='married' with partner='single'; that is, we will be calculating the effect of marriage in our natural course (Maarten's example).
intervention_rules <- list(
    partner = function(DF) DF %>% select(partner) %>% mutate(partner = replace(partner, partner == 'married', 'single'))
)

# Lastly, we need a set of rules for each specific effect to be simulated. Below are the rules for simulating the direct effect
# using the natural and intervention courses. These leverage simScenario(), which draws stochastically from the DF provided. 
# I need to add the rule set for the indirect effects. I think for whatever you want the effect for you draw from the intervention,
# and draw from the natural for everything else. The only thing that makes it a "direct" vs. "indirect" effect is whether or not it
# was the variable actually being intervened on in creating intervention_DF.
direct_effect_rules <- list(
  agegroup = function(DF, ...) binAges(DF$age), # Still deterministic.
  birth = function(DF, natural_DF, intervention_DF, models, ...) simPredict(DF, models, 1), # Still probabilistic based on model.
  job = function(DF, natural_DF, intervention_DF, models) simScenario(DF, natural_DF, models, 2), # Draw stochastically from natural course.
  partner = function(DF, natural_DF, intervention_DF, models) simScenario(DF, intervention_DF, models, 3), # Draw stochastically from intervention.
  totaledu = function(DF, ...) DF$totaledu + (DF$job == "edu"), # Still deterministic.
  totalbirth = function(DF, ...) DF$totalbirth + DF$birth, # Still deterministic.
  censor = function(DF, natural_DF, intervention_DF, models, ...) simPredict(DF, models, 4) # Still probabilistic based on model (?).
)

boots <- 15 # Number of bootstraps, 100 takes a while
replicationSize <- 6 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)

set.seed(123)

if(!file.exists("./bootruns.Rds")) {
  
    bootruns <- mclapply(1:boots, function(b) {
      
        # Sample individuals with replacement not rows
        sampleDF <- DF %>%
            select(id) %>%
            unique %>%
            sample_frac(replace=TRUE) %>%
            right_join(DF)
        
        # Fit the model
        gfitboot <- gfit.init(formulas, families, functions, data=sampleDF)
        
        # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
        mcDF <- bind_rows(lapply(1:replicationSize, function(i) sampleDF)) 
        
        # Run the "natural course" rules
        natural_DF <- progressSimulation(mcDF, lags, rules, gfitboot)
        
        # Run the "intervention course" rules
        intervention_DF <- progressSimulation(mcDF, lags, rules, gfitboot, intervention_rules)
        
        # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
        direct_effect_DF <- progressSimulation_dev(mcDF, lags, direct_effect_rules, gfitboot, natural_DF=natural_DF, intervention_DF=intervention_DF)
        
        # Simulate indirect effects
        # [to be added]
        
        # Return all courses simulated
        list(natural=natural_DF, intervention=intervention_DF, direct=direct_effect_DF)
        
    }, mc.cores=1)
    
    saveRDS(bootruns, "./bootruns.Rds")
    
}

bootruns <- read_rds("./bootruns.Rds")

natural_pSimsDF <- bind_rows(lapply(1:length(bootruns), function(i){
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
    select(-nbirths, -nmoms) %>%
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
      mutate(type="Natural")

intervention_pSimsDF <- bind_rows(lapply(1:length(bootruns), function(i){
  bootruns[[i]]$intervention %>%
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
  select(-nbirths, -nmoms) %>%
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
  mutate(type="Intervention")

direct_effect_pSimsDF <- bind_rows(lapply(1:length(bootruns), function(i){
  bootruns[[i]]$direct %>%
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
  select(-nbirths, -nmoms) %>%
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
  mutate(type="Direct effect")

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
natural_pSimsDF %>%
    bind_rows(intervention_pSimsDF) %>%
    bind_rows(direct_effect_pSimsDF) %>%
    bind_rows(actualDF) %>%
    filter(measure != "pother") %>%
    ggplot(aes(x=age, y=mean, group=type, color=type, fill=type)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
    theme_classic() +
    facet_wrap(~measure)


