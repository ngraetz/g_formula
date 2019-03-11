rm(list=ls())

## Options needed for this run.
repo <- 'E:/Data/Fragile_Families/11062018/'
local_cores <- 1
use_image <- TRUE
image_path <- 'E:/Data/Fragile_Families/11062018/'
image_path <- 'E:/Data/Fragile_Families/11142018/'
package_lib <- 'E:/Data/Fragile_Families/11072018/g_packages'
.libPaths(c(.libPaths(), package_lib))
library(session)

## Set up packages and functions.
if(!use_image) {
setwd(repo)
library(data.table)
library(tidyverse)
library(nnet)
library(parallel)
source("./gfit.R")

## Load dat clean imputed datur.
DF <- readRDS('./fragile_families_clean.RDS')

## Save image to use on secure remote server (not connected to internet)
save.image(paste0(repo, 'computing_image.RData'))
}

if(use_image) {
  setwd(repo)
  #load(paste0(image_path, 'computing_image.RData'))
  package_lib <- 'E:/Data/Fragile_Families/11072018/g_packages'
  .libPaths(c(.libPaths(), package_lib))
  library(session)
  restore.session(paste0(image_path, 'computing_image.RData'))
  source("./gfit.R")
  library(data.table, lib.loc = package_lib)
  library(tidyverse, lib.loc = package_lib)
  library(mice, lib.loc = package_lib)
  library(dplyr, lib.loc = package_lib)
  library(nnet, lib.loc = package_lib)
  library(parallel, lib.loc = package_lib)
}

## ATTACH CENSUS TRACT IDS AND CHARACTERISTICS.
## Make long by id, wave to merge to cleaned FF public-use.
geo <- fread("E:/Data/Fragile_Families/ffgeo.csv")
# poverty = pfbpl, welfare = ppuba, unemployment = puemp, female hhs = pfhhr, >HS = p25hs
covs <- c('pfbpl','ppuba','puemp','pfhhr','p25hs')
format_cov <- function(c, dt) {
  df <- dt[, c('idnum',grep(paste0('^tm.', c), names(dt), value=T)), with=F]
  df <- melt(df, id.vars='idnum', measure.vars = grep(paste0('^tm.', c), names(dt), value=T), value.name = c)
  for(w in 1:5) df[grep(paste0('tm',w), variable), wave := w]
  df[, (c) := as.numeric(get(c))]
  df[, variable := NULL]
  return(df)
}
all_covs <- lapply(covs, format_cov, dt=geo)
all_covs <- Reduce(merge, all_covs)
## Need to impute/interpolate all values.
sapply(all_covs, function(x) sum(is.na(x)))
# mi_list = mice(all_covs, method='pmm', m=5)
# imp_geo <- mice::complete(mi_list)
# imp_geo <- as.data.table(imp_geo)
# sapply(imp_geo, function(x) sum(is.na(x)))
## UNTIL I CAN LOAD MICE PACKAGE, drop missing neighborhood characteristics.
all_covs <- as.data.table(all_covs)
for(c in covs) all_covs <- all_covs[!is.na(get(c)),]
sapply(all_covs, function(x) sum(is.na(x)))

## Wodtke disadvantage index = PCA of poverty, unemployment, welfare, female-headed HHs,
## percent of residents age 25 and older without high school diploma, percent of residents
## age 25 and older in managerial or professional occupations.
## Divide into quintiles.
## Calculate duration-weighted exposure(?)
all_covs[, p25hs := 1 - p25hs]
pca <- prcomp(all_covs[, covs, with=F], center=T, scale.=T)
all_covs[, wodtke := pca$x[,1]]
all_covs[, id := idnum]

## Use first PC to create quintiles (used national distribution from Census/ACS?) For now just use raw first PC value.
## IV of interest is duration-weighted first PC (i.e. just average over all values experienced thus far in an individual's life).
## Have this be deterministic component in g-formula simulation. Do the same with other variables...?

## Also look at definition in Massey 2018 (standardize within years?)

## Merge to clean, imputed public-use dataset.
#DF <- DF %>% right_join(all_covs[, c('id','wave','wodtke')], by = c('wave','id'))
DF <- DF %>% left_join(all_covs[, c('id','wave','wodtke')], by = c('wave','id'))
DF <- DF %>% filter(!is.na(f_age)) ## For some reason, there are 47 obs in the neighborhood dataset that aren't in the imputed DF... check this later.

## COME BACK TO THIS: for now, override irregularly spaced age/year to be equal (1-5), and the put it back at the end. I don't know how to handle this.
DF <- DF %>% 
  mutate(age=wave) %>%
  mutate(year=wave) %>%
  mutate(c_sex_male=as.numeric(c_sex)) %>%
  mutate(c_sex_male=replace(c_sex_male, c_sex_male==2, 0))

DF %>% filter(m_race %in% c('white','black')) %>% group_by(m_race, id) %>% summarise(wodtke=mean(wodtke)) %>% ggplot(aes(x=wodtke, fill=m_race)) +
  geom_density(alpha=0.4) + labs(x='Wodtke index of neighborhood disadvantage (duration-weighted)',y='Proportion',fill='Mother race') + theme_minimal()

source('E:/Data/Fragile_Families/11052018/ff_functions.R')
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
    m_relation_cohab = mean(m_relation == 'cohab'),
    m_relation_contact = mean(m_relation == 'contact'),
    m_relation_married = mean(m_relation == 'married'),
    m_relation_no_contact = mean(m_relation == 'no_contact'),
    m_kids = mean(m_kids),
    m_health_excellent = mean(m_health == 'Excellent'),
    m_health_very_good = mean(m_health == 'Very good'),
    m_health_good = mean(m_health == 'Good'),
    m_health_fair = mean(m_health == 'Fair'),
    m_health_poor = mean(m_health == 'Poor'),
    m_f_in_jail = mean(m_f_in_jail),
    m_evicted = mean(m_evicted),
    m_anxiety = mean(m_anxiety),
    m_depression = mean(m_depression),
    m_job_type_full_time = mean(m_job_type == 'full_time'),
    m_job_type_part_time = mean(m_job_type == 'part_time'),
    m_job_type_unemployed = mean(m_job_type == 'unemployed'),
    m_wodtke = mean(wodtke, na.rm = T),
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
write.csv(table_1, 'E:/Data/Fragile_Families/table_1.csv', row.names = F)
DF <- DF %>% filter(!is.na(wodtke))

## Indices: id, age, year [THESE ARE HARDCODED RIGHT NOW IN OUR GFIT]
## DV = c_ppvt
## Time-invariant IVs: f_age, f_immigrant, m_age, c_sex, c_lbw, m_living_parents_15, m_immigrant, m_race, m_cognitive
## Time-variant IVs: m_relation, m_kids, m_health, m_edu, m_poverty, m_f_in_jail, m_evicted, m_anxiety, m_depression, m_job_type
## Censoring variable: censor

## Lag it all.
tc_vars <- c('f_age', 'f_immigrant', 'm_age', 'c_sex_male', 'c_lbw', 'm_living_parents_15', 'm_immigrant', 'm_cognitive', 'm_race')
tv_vars <- c('c_ppvt','m_relation', 'm_kids', 'm_health', 'm_edu', 'm_poverty', 'm_f_in_jail', 'm_evicted', 'm_anxiety', 'm_depression', 'm_job_type','wodtke')
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
  update(baseForm, wodtke ~ .),
  as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
  (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
)
families <- list(gaussian, NULL, gaussian, NULL, NULL, NULL, binomial, binomial, binomial, binomial, NULL, gaussian, binomial) # families for models
functions <- list(glm, multinom, glm, multinom, multinom, multinom, glm, glm, glm, glm, multinom, glm, glm) # functions for results

## VALIDATE MODELS
test_f <- as.formula(paste0('c_ppvt ~ ', paste(tc_vars, collapse = ' + '), ' + 
  (age) * (', paste(paste0('l.', tv_vars), collapse = ' + '), ')'))
m_num <- 1
m <- gfit.init(formulas[m_num], families[m_num], functions[m_num], data=DF)
summary(m[[1]])

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
  l.m_job_type = function(DF) DF %>% mutate(l.m_job_type=lag(m_job_type)),
  l.wodtke = function(DF) DF %>% mutate(l.wodtke=lag(wodtke))
)

# Here are the natural deterministic and probabilistic rules.
natural_rules <- list(
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
  wodtke = function(DF, models, ...) simPredict(DF, models, 12),
  censor = function(DF, models, ...) simPredict(DF, models, 13)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
# Let's test drawing neighborhood disadvantage for the black population from the natural course for the white population.
# What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
intervention_rules <- list(
  wodtke = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 12)
  #wodtke = function(DF, ...) DF %>% mutate(wodtke = replace(wodtke, m_race=='black', -1.5)) %>% select(wodtke) 
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
  m_relation = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 2),
  m_kids = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  m_health = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 4),
  m_edu = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 5),
  m_poverty = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 6),
  m_f_in_jail = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 7),
  m_evicted = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 8),
  m_anxiety = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 9),
  m_depression = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 10),
  m_job_type = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 11),
  ## Direct pathway
  #wodtke = function(DF, ...) DF %>% mutate(wodtke = replace(wodtke, m_race=='black', -1.5)) %>% select(wodtke), 
  wodtke = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 12),
  ## Censoring 
  censor = function(DF, models, ...) simPredict(DF, models, 13)
)

boots <- 15 # Number of bootstraps, 100 takes a while
replicationSize <- 5 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
# maxit=1000 (might need to increase max iterations on multinomial models, I think most are just stopping at 100 because it is the default,
# not because they converge)

set.seed(80085)
repo <- 'E:/Data/Fragile_Families/11062018/'
setwd(repo)


cor(DF %>% filter(m_race=='hispanic') %>% select(c_ppvt, wodtke), use='complete.obs')

#if(!file.exists("./ff_bootruns.Rds")) {
  
bootruns <- lapply(1:boots, function(b) {
    
    message(paste0('Bootstrap ', b))
    
    # Sample individuals with replacement not rows
    sampleDF_black <- DF %>%
      filter(m_race=='black') %>%
      select(id) %>%
      unique %>%
      sample_frac(replace=TRUE) %>%
      left_join(DF)
    
    # Fit the model
    gfitboot_black <- gfit.init(formulas, families, functions, data=sampleDF_black)
    
    # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
    mcDF_black <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_black)) 
    
    # Run the "natural course" rules
    natural_course_DF_black <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black)
  
    # Sample individuals with replacement not rows
    sampleDF_white <- DF %>%
      filter(m_race=='white') %>%
      select(id) %>%
      unique %>%
      sample_frac(replace=TRUE) %>%
      left_join(DF)

    # Fit the model
    gfitboot_white <- gfit.init(formulas, families, functions, data=sampleDF_white)

    # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
    mcDF_white <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_white))

    # Run the "natural course" rules
    natural_course_DF_white <- progressSimulation(mcDF_white, lags, natural_rules, gfitboot_white)
    
    # Run the "intervention course" rules for BLACK population: draw neighborhood disadvantage from WHITE natural course.
    # What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
    intervention_DF_black <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black, intervention_rules, natural_DF=natural_course_DF_white)
    # intervention_DF <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black, intervention_rules)
    
    # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
    direct_effect_DF <- progressSimulation(mcDF_black, lags, direct_effect_rules, gfitboot_black, natural_DF=natural_course_DF_black, intervention_DF=intervention_DF_black)

    # Simulate indirect effects [not tested]
    # indirect_job_effect_DF <- progressSimulation(mcDF, lags, indirect_job_effect_rules, gfitboot, natural_DF=natural_DF, intervention_DF=intervention_DF)
    
    # Return all courses simulated
    list(natural_black=natural_course_DF_black,
         natural_white=natural_course_DF_white,
         intervention=intervention_DF_black,
         direct=direct_effect_DF)

  })
  
saveRDS(bootruns, "./ff_bootruns_intervention.Rds")
  
#}

source("./gfit.R")
bootruns <- lapply(1:boots, function(b) {
  
  message(paste0('Bootstrap ', b))
  
  # Sample individuals with replacement not rows
  sampleDF_all <- DF %>%
    #filter(m_race=='black') %>%
    select(id) %>%
    unique %>%
    sample_frac(replace=TRUE) %>%
    left_join(DF)
  
  # Fit the model
  gfitboot_all <- gfit.init(formulas, families, functions, data=sampleDF_all)
  
  # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
  mcDF_all <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_all)) 
  
  # Run the "natural course" rules
  natural_course_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all)
  
  # Run the "intervention course" rules for BLACK population: draw neighborhood disadvantage from WHITE natural course.
  # What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
  # intervention_DF <- progressSimulation(mcDF_black, lags, rules, gfitboot_black, intervention_rules, natural_DF=natural_course_DF_white)
  # intervention_DF_black <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules)
  intervention_DF_black <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules, natural_DF=natural_course_DF_all %>% filter(m_race=='white'))
  
  # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
  direct_effect_DF_black <- progressSimulation(mcDF_all, lags, direct_effect_rules, gfitboot_all, natural_DF=natural_course_DF_all %>% filter(m_race=='black'), intervention_DF=intervention_DF_black)
  
  # Simulate indirect effects [not tested]
  # indirect_job_effect_DF <- progressSimulation(mcDF, lags, indirect_job_effect_rules, gfitboot, natural_DF=natural_DF, intervention_DF=intervention_DF)
  
  # Return all courses simulated
  list(natural_black=natural_course_DF_all %>% filter(m_race=='black'),
       natural_white=natural_course_DF_all %>% filter(m_race=='white'),
       intervention=intervention_DF_black %>% filter(m_race=='black'),
       direct=direct_effect_DF_black %>% filter(m_race=='black'))
  
})

saveRDS(bootruns, "./ff_bootruns_intervention.Rds")

#}
  
#bootruns <- read_rds("./bootruns.Rds")

## Aggregate
all_runs <- rbindlist(lapply(c('natural_black','natural_white','intervention','direct'), compile_sims, bs = bootruns))
all_runs <- redo_ages(all_runs)

actualDF <- DF %>%
  group_by(age, m_race) %>%
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
    pparttime = mean(m_job_type == "part_time"),
    mwodtke = mean(wodtke)
  ) %>%
  gather("measure", "mean", -m_race, -age)
actualDF <- redo_ages(actualDF)

## Plot 
cov_plot <- all_runs %>%
  #filter(!(measure %in% c('mkids',"mppvt","mwodtke"))) %>%
  mutate(type=replace(type, type=='natural_white', 'white')) %>%
  mutate(type=replace(type, type=='natural_black', 'black')) %>%
  filter((measure=='mppvt' & age>=3)|(measure!='mppvt')) %>%
  filter(measure!='censor') %>%
  mutate(measure=factor(measure, levels = c('mppvt','mwodtke', unique(actualDF$measure[!(actualDF$measure %in% c('mppvt','mwodtke'))]))))
cov_obs <- actualDF %>%
  #filter(!(measure %in% c('mkids',"mppvt","mwodtke")) & m_race %in% c('white','black')) %>%
  filter(m_race %in% c('white','black')) %>%
  mutate(type=m_race) %>%
  filter((measure=='mppvt' & age>=3)|(measure!='mppvt')) %>%
  filter(measure!='censor') %>%
  mutate(measure=factor(measure, levels = c('mppvt','mwodtke', unique(actualDF$measure[!(actualDF$measure %in% c('mppvt','mwodtke'))]))))
ggplot(data=cov_plot %>% filter(!(measure %in% c('mwodtke','mppvt'))), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  geom_point(data=cov_obs %>% filter(!(measure %in% c('mwodtke','mppvt'))), aes(fill=type), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated\ncourse'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')
ggplot(data=cov_plot %>% filter((measure %in% c('mwodtke','mppvt'))), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  geom_point(data=cov_obs %>% filter((measure %in% c('mwodtke','mppvt'))), aes(fill=type), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated\ncourse'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')



cov_plot <- all_runs %>%
  filter((measure %in% c("mppvt"))) %>%
  mutate(type=replace(type, type=='natural_white', 'white')) %>%
  mutate(type=replace(type, type=='natural_black', 'black'))
cov_obs <- actualDF %>%
  filter((measure %in% c("mppvt")) & m_race %in% c('white','black')) %>%
  mutate(type=m_race)
ggplot(data=cov_plot, aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=type), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure, scales='free_y')

cov_plot <- all_runs %>%
  filter(measure %in% c("mwodtke",'mppvt')) %>%
  mutate(type=replace(type, type=='natural_white', 'white')) %>%
  mutate(type=replace(type, type=='natural_black', 'black')) %>%
  filter((measure=='mppvt' & age>=3)|(measure=='mwodtke'))
cov_obs <- actualDF %>%
  filter((measure %in% c("mwodtke",'mppvt')) & m_race %in% c('white','black')) %>%
  mutate(type=m_race) %>%
  filter((measure=='mppvt' & age>=3)|(measure=='mwodtke'))
ggplot(data=cov_plot, aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=type), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure, scales = 'free_y')








natural_pSimsDF_black <- bind_rows(lapply(1:length(bootruns_black), function(i){
  bootruns_black[[i]]$natural %>%
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
    pparttime = mean(m_job_type == "part_time"),
    mwodtke = mean(wodtke)
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
  mutate(type="Natural") %>%
  mutate(m_race="black")

natural_pSimsDF_white <- bind_rows(lapply(1:length(bootruns_white), function(i){
  bootruns_white[[i]]$natural %>%
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
    pparttime = mean(m_job_type == "part_time"),
    mwodtke = mean(wodtke)
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
  mutate(type="Natural") %>%
  mutate(m_race="white")

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
  group_by(age, m_race) %>%
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
    pparttime = mean(m_job_type == "part_time"),
    mwodtke = mean(wodtke)
  ) %>%
  gather("measure", "mean", -m_race, -age) %>%
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
natural_pSimsDF_black <- redo_ages(natural_pSimsDF_black)
natural_pSimsDF_white <- redo_ages(natural_pSimsDF_white)
# intervention_pSimsDF <- redo_ages(intervention_pSimsDF)
# direct_effect_pSimsDF <- redo_ages(direct_effect_pSimsDF)
actualDF <- redo_ages(actualDF)

## Plots
cov_plot <- natural_pSimsDF_black %>%
  bind_rows(natural_pSimsDF_white) %>%
  # bind_rows(intervention_pSimsDF) %>%
  # bind_rows(direct_effect_pSimsDF) %>%
  # bind_rows(actualDF %>% filter(m_race=='black')) %>%
  filter(!(measure %in% c('mkids',"mppvt","mwodtke")))
cov_obs <- actualDF %>%
  filter(!(measure %in% c('mkids',"mppvt","mwodtke")) & m_race %in% c('white','black'))
ggplot(data=cov_plot, aes(x=age, y=mean)) +
  geom_line(aes(color=m_race), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=m_race), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=m_race), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure)

cov_plot <- natural_pSimsDF_black %>%
  bind_rows(natural_pSimsDF_white) %>%
  # bind_rows(intervention_pSimsDF) %>%
  # bind_rows(direct_effect_pSimsDF) %>%
  # bind_rows(actualDF %>% filter(m_race=='black')) %>%
  filter((measure %in% c("mwodtke")))
cov_obs <- actualDF %>%
  filter((measure %in% c("mwodtke")) & m_race %in% c('white','black'))
ggplot(data=cov_plot, aes(x=age, y=mean)) +
  geom_line(aes(color=m_race), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=m_race), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=m_race), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure)

cov_plot <- natural_pSimsDF_black %>%
  bind_rows(natural_pSimsDF_white) %>%
  # bind_rows(intervention_pSimsDF) %>%
  # bind_rows(direct_effect_pSimsDF) %>%
  # bind_rows(actualDF %>% filter(m_race=='black')) %>%
  filter((measure %in% c("mppvt")))
cov_obs <- actualDF %>%
  filter((measure %in% c("mppvt","mwodtke")) & m_race %in% c('white','black'))
ggplot(data=cov_plot, aes(x=age, y=mean)) +
  geom_line(aes(color=m_race), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=m_race), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=m_race), shape=21, color='black', size=5) + 
  theme_classic() +
  facet_wrap(~measure)

## GRID ARRANGE A 1*3 OF COVS (FACET), BIG IV, BIG DV
