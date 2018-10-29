repo <- 'C:/Users/Nick/Documents/repos/g_formula/'
local_cores <- 1

## Set up packages and functions.
setwd(repo)
library(tidyverse)
library(nnet)
library(parallel)

binAges <- function(x) {
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
saveRDS(DF, './data/maarten.RDS')
