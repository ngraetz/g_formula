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

## TABLE 1: wide by wave/race, long by variables
## All variables
tc_vars <- c('f_age', 'f_immigrant', 'm_age', 'c_sex', 'c_lbw', 'm_living_parents_15', 'm_immigrant', 'm_race', 'm_cognitive')
tv_vars <- c('c_ppvt','m_relation', 'm_kids', 'm_health', 'm_edu', 'm_poverty', 'm_f_in_jail', 'm_evicted', 'm_anxiety', 'm_depression', 'm_job_type')
DF <- as.data.table(DF)

## Add wave-specific censoring, neighborhood disadvantage, and then duration-weighted neighborhood advantage under last wave.
