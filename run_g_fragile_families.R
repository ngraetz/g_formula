rm(list=ls())

# real_data %>%
#   multiple_imputation() %>%
#   replicate_survey_weights() %>%
#   g_formula() %>%
#     fit_survey_glms %>%
#     pull_vcov %>%
#     bootstrap %>%
#   visual()
  
## Options needed for this run.
local_cores <- 1
use_image <- FALSE
image_path <- ''

## Set up packages and functions.
if(!use_image) {
library(session)
library(data.table)
library(tidyverse)
library(nnet)
library(parallel)
library(ggplot2)
library(haven)
library(naniar)
library(mice)
library(survey)
source("./gfit.R")

## Load dat clean imputed datur.
DF <- readRDS('./fragile_families_clean.RDS')

## Save image to use on secure remote server (not connected to internet)
save.image(paste0(repo, 'computing_image.RData'))
}
if(use_image) {
  load(paste0(repo, 'computing_image.RData'))
}

## COME BACK TO THIS: for now, override irregularly spaced age/year to be equal (1-5), and the put it back at the end. I don't know how to handle this.
DF[, age := wave]
DF[, year := wave]
DF <- DF %>% 
  mutate(age=wave) %>%
  mutate(year=wave)

## Indices: id, age, year [THESE ARE HARDCODED RIGHT NOW IN OUR GFIT]
## DV = c_ppvt
## Time-invariant IVs: f_age, f_immigrant, m_age, c_sex, c_lbw, m_living_parents_15, m_immigrant, m_race, m_cognitive
## Time-variant IVs: m_relation, m_kids, m_health, m_edu, m_poverty, m_f_in_jail, m_evicted, m_anxiety, m_depression, m_job_type
## Censoring variable: censor

## Lag it all.
tc_vars <- c('f_age', 'f_immigrant', 'm_age', 'c_sex', 'c_lbw', 'm_living_parents_15', 'm_immigrant', 'm_race', 'm_cognitive')
tv_vars <- c('c_ppvt','m_relation', 'm_kids', 'm_health', 'm_edu', 'm_poverty', 'm_f_in_jail', 'm_evicted', 'm_anxiety', 'm_depression', 'm_job_type')
DF <- as.data.table(DF)
DF[, (paste0('l.', tv_vars)) :=  shift(.SD), by='id', .SDcols=tv_vars]

# build up the model structure base
baseForm <- as.formula(paste0('~ ', paste(tc_vars, collapse = ' + '), ' + 
  (age) * (', paste(paste0('l.', tv_vars), collapse = ' + '), ')')) 

# formula for models
formulas <- list(
  update(baseForm, c_ppvt ~ .),
  update(baseForm, m_relation ~ .),
  update(baseForm, m_kids ~ .),
  update(baseForm, m_health ~ .),
  update(baseForm, m_edu ~ .),
  update(baseForm, m_poverty ~ .),
  update(baseForm, m_f_in_jail ~ .),
  update(baseForm, m_evicted ~ .),
  update(baseForm, m_anxiety ~ .),
  update(baseForm, m_depression ~ .),
  update(baseForm, m_job_type ~ .),
  as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
  (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
)
families <- list(gaussian, NULL, gaussian, NULL, NULL, NULL, binomial, binomial, binomial, binomial, NULL, binomial) # families for models
functions <- list(glm, multinom, glm, multinom, multinom, multinom, glm, glm, glm, glm, multinom, glm) # functions for results

# here are our lags that we want to keep updated
lags <- list(
  l.c_ppvt = function(DF) DF %>% mutate(l.c_ppvt=lag(c_ppvt)),
  l.m_relation = function(DF) DF %>% mutate(l.m_relation=lag(m_relation)),
  l.m_kids = function(DF) DF %>% mutate(l.m_kids=lag(m_kids)),
  l.m_health = function(DF) DF %>% mutate(l.m_health=lag(m_health)),
  l.m_edu = function(DF) DF %>% mutate(l.m_edu=lag(m_edu)),
  l.m_poverty = function(DF) DF %>% mutate(l.m_poverty=lag(m_poverty)),
  l.m_f_in_jail = function(DF) DF %>% mutate(l.m_f_in_jail=lag(m_f_in_jail)),
  l.m_evicted = function(DF) DF %>% mutate(l.m_evicted=lag(m_evicted)),
  l.m_anxiety = function(DF) DF %>% mutate(l.m_anxiety=lag(m_anxiety)),
  l.m_depression = function(DF) DF %>% mutate(l.m_depression=lag(m_depression)),
  l.m_job_type = function(DF) DF %>% mutate(l.m_job_type=lag(m_job_type))
)

# Here are the natural deterministic and probabilistic rules.
rules <- list(
  # Deterministic
    # None for now, but need to add "accumulating disadvantage" variables, like how many years lived in poor neighborhood, years in poverty, etc.
  # Probabilistic 
  c_ppvt = function(DF, models, ...) simPredict(DF, models, 1),
  m_relation = function(DF, models, ...) simPredict(DF, models, 2),
  m_kids = function(DF, models, ...) simPredict(DF, models, 3),
  m_health = function(DF, models, ...) simPredict(DF, models, 4),
  m_edu = function(DF, models, ...) simPredict(DF, models, 5),
  m_poverty = function(DF, models, ...) simPredict(DF, models, 6),
  m_f_in_jail = function(DF, models, ...) simPredict(DF, models, 7),
  m_evicted = function(DF, models, ...) simPredict(DF, models, 8),
  m_anxiety = function(DF, models, ...) simPredict(DF, models, 9),
  m_depression = function(DF, models, ...) simPredict(DF, models, 10),
  m_job_type = function(DF, models, ...) simPredict(DF, models, 11),
  censor = function(DF, models, ...) simPredict(DF, models, 12)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
# Let's test no HH poverty 
intervention_rules <- list(
  m_poverty = function(DF) DF %>% select(m_poverty) %>% mutate(m_poverty = replace(m_poverty, m_poverty %in% c("perc_0_49","perc_50_99"), "perc_100_199"))
)

# Lastly, we need a set of rules for each specific effect to be simulated. Below are the rules for simulating the direct effect
# using the natural and intervention courses. These leverage simScenario(), which draws stochastically from the DF provided. 
# I need to add the rule set for the indirect effects. I think for whatever you want the effect for you draw from the intervention,
# and draw from the natural for everything else. The only thing that makes it a "direct" vs. "indirect" effect is whether or not it
# was the variable actually being intervened on in creating intervention_DF. I think we could automate this to create
# a set of rules for each indirect effect implied by all the variables in the models, and automatically add a simulation step for
# each to the bootruns mclapply() below. I think even with a lot of variable this wouldn't add a ton of run time.
direct_effect_rules <- list(
  ## DV of interest
  c_ppvt = function(DF, models, ...) simPredict(DF, models, 1),
  ## Indirect pathways
  m_relation = function(DF, models, ...) simScenario(DF, natural_DF, models, 2),
  m_kids = function(DF, models, ...) simScenario(DF, natural_DF, models, 3),
  m_health = function(DF, models, ...) simScenario(DF, natural_DF, models, 4),
  m_edu = function(DF, models, ...) simScenario(DF, natural_DF, models, 5),
  ## Direct pathway
  m_poverty = function(DF, models, ...) simScenario(DF, intervention_DF, models, 6),
  ## Indirect pathways
  m_f_in_jail = function(DF, models, ...) simScenario(DF, natural_DF, models, 7),
  m_evicted = function(DF, models, ...) simScenario(DF, natural_DF, models, 8),
  m_anxiety = function(DF, models, ...) simScenario(DF, natural_DF, models, 9),
  m_depression = function(DF, models, ...) simScenario(DF, natural_DF, models, 10),
  m_job_type = function(DF, models, ...) simScenario(DF, natural_DF, models, 11),
  ## Censoring
  censor = function(DF, models, ...) simPredict(DF, models, 12)
)

boots <- 15 # Number of bootstraps, 100 takes a while
replicationSize <- 5 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
# maxit=1000 (might need to increase max iterations on multinomial models, I think most are just stopping at 100 because it is the default,
# not because they converge)

set.seed(80085)

#if(!file.exists("./ff_bootruns.Rds")) {
  
  bootruns <- mclapply(1:boots, function(b) {
    
    # Sample individuals with replacement not rows
    sampleDF <- DF %>%
      select(id) %>%
      unique %>%
      sample_frac(replace=TRUE) %>%
      left_join(DF)

    # Fit the model
    gfitboot <- gfit.init(formulas, families, functions, data=sampleDF)
    
    # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
    mcDF <- bind_rows(lapply(1:replicationSize, function(i) sampleDF)) 
    
    # Run the "natural course" rules
    natural_course_DF <- progressSimulation(mcDF, lags, rules, gfitboot)
    
    # Run the "intervention course" rules
    intervention_DF <- progressSimulation(mcDF, lags, rules, gfitboot, intervention_rules)
    
    # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
    direct_effect_DF <- progressSimulation(mcDF, lags, direct_effect_rules, gfitboot, natural_DF=natural_course_DF, intervention_DF=intervention_DF)
    
    # Simulate indirect effects [not tested]
    # indirect_job_effect_DF <- progressSimulation(mcDF, lags, indirect_job_effect_rules, gfitboot, natural_DF=natural_DF, intervention_DF=intervention_DF)
    
    # Return all courses simulated
    list(natural=natural_DF, intervention=intervention_DF, direct=direct_effect_DF)

  }, mc.cores=local_cores)
  
  saveRDS(bootruns, "./ff_bootruns_poverty_effects.Rds")
  
#}
  
bootruns <- read_rds("./bootruns.Rds")

natural_pSimsDF <- bind_rows(lapply(1:length(bootruns), function(i){
  bootruns[[i]]$natural %>%
    mutate(sim=i)})) %>%
  group_by(age, sim) %>%
  summarize(
    mppvt = mean(c_ppvt), 
    pcohab = mean(m_relation == "cohab"),
    pmarried = mean(m_relation == "married"),
    pnocontact = mean(m_relation == "no_contact"),
    pcontact = mean(m_relation == "contact"),
    mkids = mean(m_kids), 
    pexcellenthealth = mean(m_health == "Excellent"),
    ppoorhealth = mean(m_health == "Poor"),
    pcollege = mean(m_edu == "college"),
    psomecollege = mean(m_edu == "some_college"),
    phs = mean(m_edu == "hs"),
    plesshs = mean(m_edu == "less_hs"),
    phighincome = mean(m_poverty == "perc_300_plus"),
    ppov = mean(m_poverty == "perc_50_99"),
    pextpov = mean(m_poverty == "perc_0_49"),
    pfatherjail = mean(m_f_in_jail), 
    pevicted = mean(m_evicted), 
    pmotheranxiety = mean(m_anxiety), 
    pmotherdepression = mean(m_depression), 
    pfulltime = mean(m_job_type == "full_time"),
    punemployed = mean(m_job_type == "unemployed"),
    pparttime = mean(m_job_type == "part_time")
    ) %>%
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
    mppvt = mean(c_ppvt), 
    pcohab = mean(m_relation == "cohab"),
    pmarried = mean(m_relation == "married"),
    pnocontact = mean(m_relation == "no_contact"),
    pcontact = mean(m_relation == "contact"),
    mkids = mean(m_kids), 
    pexcellenthealth = mean(m_health == "Excellent"),
    ppoorhealth = mean(m_health == "Poor"),
    pcollege = mean(m_edu == "college"),
    psomecollege = mean(m_edu == "some_college"),
    phs = mean(m_edu == "hs"),
    plesshs = mean(m_edu == "less_hs"),
    phighincome = mean(m_poverty == "perc_300_plus"),
    ppov = mean(m_poverty == "perc_50_99"),
    pextpov = mean(m_poverty == "perc_0_49"),
    pfatherjail = mean(m_f_in_jail), 
    pevicted = mean(m_evicted), 
    pmotheranxiety = mean(m_anxiety), 
    pmotherdepression = mean(m_depression), 
    pfulltime = mean(m_job_type == "full_time"),
    punemployed = mean(m_job_type == "unemployed"),
    pparttime = mean(m_job_type == "part_time")
  ) %>%
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
    mppvt = mean(c_ppvt), 
    pcohab = mean(m_relation == "cohab"),
    pmarried = mean(m_relation == "married"),
    pnocontact = mean(m_relation == "no_contact"),
    pcontact = mean(m_relation == "contact"),
    mkids = mean(m_kids), 
    pexcellenthealth = mean(m_health == "Excellent"),
    ppoorhealth = mean(m_health == "Poor"),
    pcollege = mean(m_edu == "college"),
    psomecollege = mean(m_edu == "some_college"),
    phs = mean(m_edu == "hs"),
    plesshs = mean(m_edu == "less_hs"),
    phighincome = mean(m_poverty == "perc_300_plus"),
    ppov = mean(m_poverty == "perc_50_99"),
    pextpov = mean(m_poverty == "perc_0_49"),
    pfatherjail = mean(m_f_in_jail), 
    pevicted = mean(m_evicted), 
    pmotheranxiety = mean(m_anxiety), 
    pmotherdepression = mean(m_depression), 
    pfulltime = mean(m_job_type == "full_time"),
    punemployed = mean(m_job_type == "unemployed"),
    pparttime = mean(m_job_type == "part_time")
  ) %>%
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
    mppvt = mean(c_ppvt), 
    pcohab = mean(m_relation == "cohab"),
    pmarried = mean(m_relation == "married"),
    pnocontact = mean(m_relation == "no_contact"),
    pcontact = mean(m_relation == "contact"),
    mkids = mean(m_kids), 
    pexcellenthealth = mean(m_health == "Excellent"),
    ppoorhealth = mean(m_health == "Poor"),
    pcollege = mean(m_edu == "college"),
    psomecollege = mean(m_edu == "some_college"),
    phs = mean(m_edu == "hs"),
    plesshs = mean(m_edu == "less_hs"),
    phighincome = mean(m_poverty == "perc_300_plus"),
    ppov = mean(m_poverty == "perc_50_99"),
    pextpov = mean(m_poverty == "perc_0_49"),
    pfatherjail = mean(m_f_in_jail), 
    pevicted = mean(m_evicted), 
    pmotheranxiety = mean(m_anxiety), 
    pmotherdepression = mean(m_depression), 
    pfulltime = mean(m_job_type == "full_time"),
    punemployed = mean(m_job_type == "unemployed"),
    pparttime = mean(m_job_type == "part_time")
  ) %>%
  gather("measure", "mean", -age) %>%
  mutate(type="Observed")

# Put ages that aren't actually sequential back to their original observed values. 
redo_ages <- function(df) {
new_df <- as.data.table(df)
new_df[age==1, new_age := 0]
new_df[age==2, new_age := 1]
new_df[age==3, new_age := 3]
new_df[age==4, new_age := 5]
new_df[age==5, new_age := 9]
new_df[, age := new_age]
return(new_df)
}
natural_pSimsDF <- redo_ages(natural_pSimsDF)
intervention_pSimsDF <- redo_ages(intervention_pSimsDF)
direct_effect_pSimsDF <- redo_ages(direct_effect_pSimsDF)
actualDF <- redo_ages(actualDF)

## Plots
natural_pSimsDF %>%
  bind_rows(intervention_pSimsDF) %>%
  bind_rows(direct_effect_pSimsDF) %>%
  bind_rows(actualDF) %>%
  filter(!(measure %in% c('mkids',"mppvt"))) %>%
  ggplot(aes(x=age, y=mean, group=type, color=type, fill=type)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  theme_classic() +
  facet_wrap(~measure)

natural_pSimsDF %>%
  bind_rows(intervention_pSimsDF) %>%
  bind_rows(direct_effect_pSimsDF) %>%
  bind_rows(actualDF) %>%
  filter(measure %in% c("mppvt")) %>%
  ggplot(aes(x=age, y=mean, group=type, color=type, fill=type)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) +
  theme_classic() +
  facet_wrap(~measure)
