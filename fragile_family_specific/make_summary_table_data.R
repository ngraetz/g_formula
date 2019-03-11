## Make summary tables for paper.
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
repo <- 'C:/Users/ngraetz/Documents/repos/g_formula/'
local_cores <- 1
use_image <- FALSE
image_path <- ''

## Set up packages and functions.
if(!use_image) {
  setwd(repo)
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
  source("C:/Users/ngraetz/Documents/repos/fragile_familes/ff_functions.R")
  
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

## How many individuals at baseline have each number of total observations
DF %>%
  filter(censor!=1) %>%
  group_by(id) %>%
  count() %>%
  ungroup() %>%
  count(n)

## Indices: id, age, year [THESE ARE HARDCODED RIGHT NOW IN OUR GFIT]
## DV = c_ppvt
## Time-invariant IVs: f_age, f_immigrant, m_age, c_sex, c_lbw, m_living_parents_15, m_immigrant, m_race, m_cognitive
## Time-variant IVs: m_relation, m_kids, m_health, m_edu, m_poverty, m_f_in_jail, m_evicted, m_anxiety, m_depression, m_job_type
## Censoring variable: censor

## TABLE 1: wide by wave/race, long by variables
## All variables
tc_vars <- c('f_age', 'f_immigrant', 'm_age', 'c_sex', 'c_lbw', 'm_living_parents_15', 'm_immigrant', 'm_race', 'm_cognitive')
tv_vars <- c('c_ppvt','m_relation', 'm_kids', 'm_health', 'm_edu', 'm_poverty', 'm_f_in_jail', 'm_evicted', 'm_anxiety', 'm_depression', 'm_job_type')
DF <- as.data.table(DF)

## Add wave-specific censoring, neighborhood disadvantage, and then duration-weighted neighborhood advantage under last wave.
table_1 <- DF %>%
  filter(m_race %in% c('black','white')) %>%
  group_by(wave, m_race) %>%
  summarize(
    f_age = mean(f_age), 
    f_immigrant = mean(f_immigrant),
    m_age = mean(m_age),
    c_sex_male = mean(c_sex == 1),
    c_sex_female = mean(c_sex == 2),
    c_lbw = mean(c_lbw), 
    m_living_parents_15 = mean(m_living_parents_15),
    m_immigrant = mean(m_immigrant),
    m_cognitive = mean(m_cognitive),
    m_relation_cohab = mean(m_relation == "cohab"),
    m_relation_contact = mean(m_relation == "contact"),
    m_relation_married = mean(m_relation == "married"),
    m_relation_no_contact = mean(m_relation == "no_contact"),
    m_kids = mean(m_kids),
    m_health_excellent = mean(m_health == 'Excellent'), 
    m_health_very_good = mean(m_health == 'Very good'), 
    m_health_good = mean(m_health == 'Good'), 
    m_health_fair = mean(m_health == 'Fair'), 
    m_health_poor = mean(m_health == "Poor"),
    m_f_in_jail = mean(m_f_in_jail),
    m_evicted = mean(m_evicted),
    m_anxiety = mean(m_anxiety),
    m_depression = mean(m_depression),
    m_job_type_full_time = mean(m_job_type == 'full_time'),
    m_job_type_part_time = mean(m_job_type == 'part_time'),
    m_job_type_unemployed = mean(m_job_type == 'unemployed'),
    c_ppvt = mean(c_ppvt),
    censor = mean(censor)
  ) %>%
  ungroup() %>%
  gather(variable, value, -m_race, -wave) %>%
  mutate(value = replace(value, wave!=1 & variable %in% get_tc_vars(), '-')) %>%
  mutate(wave=as.character(wave)) %>% 
  unite(temp, m_race, wave) %>%
  spread(temp, value) %>%
  right_join(clean_cov_names()) %>%
  arrange(cov_sort) %>%
  select(-variable, -cov_sort) %>%
  select(name, everything())
  
saveRDS(table_1, 'C:/Users/ngraetz/Documents/repos/fragile_familes/output/table_1.RDS')
