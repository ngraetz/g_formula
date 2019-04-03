rm(list=ls())

## Options needed for this run.
repo <- 'E:/Data/Fragile_Families/'
local_cores <- 1
use_image <- TRUE
image_path <- 'E:/Data/Fragile_Families/11062018/'
image_path <- 'E:/Data/Fragile_Families/11142018/'
package_lib <- 'E:/Data/Fragile_Families/uploads/11072018/g_packages'
.libPaths(c(.libPaths(), package_lib))
library(session)
out_dir <- 'E:/Data/Fragile_Families/results'

## Set up packages and functions.
if(!use_image) {
  setwd(repo)
  library(data.table)
  library(tidyverse)
  library(nnet)
  library(parallel)
  source("./gfit.R")
  
  ## Load dat clean imputed datur.
  #old_DF <- readRDS('E:/Data/Fragile_Families/03122019/fragile_families_clean_weights.RDS')
  # ff.design <- svydesign(id=~id, weights=~mweight, data=old_DF)
  # test <- make_svy_ci('~c_ppvt')
  DF <- readRDS('E:/Data/Fragile_Families/uploads/03272019/noncensored_fragile_families_clean_weights.RDS')
  DF[wave==1, wave := 2]
  DF[wave==3, wave := 3]
  DF[wave==5, wave := 4]
  DF[wave==9, wave := 5]
  DF[wave==0, wave := 1]
  DF[wave==1, c_ppvt := NA]
  DF[wave==2, c_ppvt := NA]
  ff.design <- svydesign(id=~id, weights=~mweight, data=DF)
  test <- make_svy_ci('~c_ppvt')
  # DF <- DF[, c('id','wave','m_relation')]
  # setnames(DF, 'm_relation', 'new_m_relation')
  # DF <- merge(old_DF, DF, by=c('wave','id'))
  # DF[, m_relation := new_m_relation]
  
  ## Save image to use on secure remote server (not connected to internet)
  save.image(paste0(repo, 'computing_image.RData'))
}

if(use_image) {
  setwd(repo)
  #load(paste0(image_path, 'computing_image.RData'))
  package_lib <- 'E:/Data/Fragile_Families/uploads/11072018/g_packages'
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
  library(survey)
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
mi_list = mice(all_covs, method='pmm', m=5)
imp_geo <- mice::complete(mi_list)
imp_geo <- as.data.table(imp_geo)
sapply(imp_geo, function(x) sum(is.na(x)))
all_covs <- as.data.table(imp_geo)
## UNTIL I CAN LOAD MICE PACKAGE, drop missing neighborhood characteristics.
#all_covs <- as.data.table(all_covs)
#for(c in covs) all_covs <- all_covs[!is.na(get(c)),]
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

## COME BACK TO THIS: for now, override irregularly spaced age/year to be equal (1-5), and the put it back at the end. I don't know how to handle this
## for modelling lagged values... simulating forward 1 year at a time is easy enough, but we can't fit the models that way...
DF <- DF %>% 
  mutate(real_age=age) %>%
  mutate(age=wave) %>%
  mutate(year=wave) %>%
  mutate(c_sex_male=as.numeric(c_sex)) %>%
  mutate(c_sex_male=replace(c_sex_male, c_sex_male==2, 0))

source('E:/Data/Fragile_Families/11052018/ff_functions.R')
DF[, N := 1]
table_1_mean <- DF %>%
  #filter(m_race %in% c('black','white')) %>%
  group_by(wave) %>%
  summarize(
    c_sex_male = mean(c_sex == 1),
    c_lbw = mean(c_lbw),
    f_age = mean(f_age),
    f_immigrant = mean(f_immigrant),
    m_age = mean(m_age),
    m_living_parents_15 = mean(m_living_parents_15),
    m_immigrant = mean(m_immigrant),
    m_cognitive = mean(m_cognitive),
    m_racewhite = mean(m_race=='white'),
    m_raceblack = mean(m_race=='black'),
    m_racehispanic = mean(m_race=='hispanic'),
    m_raceother = mean(m_race=='other'),
    bin_unemployed = mean(m_job_type == 'unemployed'),
    bin_pov = mean(bin_pov),
    m_f_in_jail = mean(m_f_in_jail),
    m_depression = mean(m_depression),
    f_absent = mean(m_relation == 'no_contact'),
    m_health_poor = mean(m_health == 'Poor'),
    m_wodtke = mean(wodtke, na.rm = T),
    c_ppvt = mean(original_ppvt),
    N = sum(N)
  ) %>%
  ungroup() %>%
  mutate(percent = N / 4886) %>%
  gather(variable, value, -wave) %>%
  mutate(value = replace(value, wave!=1 & variable %in% get_tc_vars(), '-')) %>%
  mutate(wave=as.character(wave)) %>%
  unite(temp, wave) %>%
  spread(temp, value) %>%
  left_join(clean_cov_names()) %>%
  arrange(cov_sort) %>%
  select(-variable, -cov_sort) %>%
  select(name, everything())
# table_1_sd <- DF %>%
#   #filter(m_race %in% c('black','white')) %>%
#   group_by(wave) %>%
#   summarize(
#     f_age = sd(f_age),
#     m_age = sd(m_age),
#     m_cognitive = sd(m_cognitive),
#     m_kids = sd(m_kids),
#     m_wodtke = sd(wodtke, na.rm = T),
#     c_ppvt = sd(original_ppvt)
#   ) %>%
#   ungroup() %>%
#   gather(variable, value, -wave) %>%
#   mutate(value = replace(value, wave!=1 & variable %in% get_tc_vars(), '-')) %>%
#   mutate(wave=as.character(wave)) %>%
#   unite(temp, wave) %>%
#   spread(temp, value) %>%
#   right_join(clean_cov_names()) %>%
#   arrange(cov_sort) %>%
#   select(-variable, -cov_sort) %>%
#   select(name, everything())
# # setnames(table_1_mean, c('name','m1','m2','m3','m4','m5'))
# # setnames(table_1_sd, c('name','sd1','sd2','sd3','sd4','sd5'))
# # table1 <- merge(table_1_mean, table_1_sd, by='name')
# # calc_sd <- function(p) {
# #   sqrt(p * (1-p) / 4886)
# # }
# # table1 <- as.data.table(table1)
# # table1[, c('sd1','sd2','sd3','sd4','sd5','m1','m2','m3','m4','m5') := lapply(.SD,as.numeric), .SDcols = c('sd1','sd2','sd3','sd4','sd5','m1','m2','m3','m4','m5')]
# # for(i in c('1','2','3','4','5')) table1[is.na(get(paste0('sd',i))), (paste0('sd',i)) := calc_sd(get(paste0('m',i)))]
write.csv(table_1_mean, paste0(out_dir, '/table_1.csv'))
DF <- DF %>% filter(!is.na(wodtke))

## Make trend plot of crude cognition scores by neighborhood disadvantage.
ggplot() + 
  geom_point(data=DF,
             aes(x=c_ppvt,
                 y=d.wodtke))
test <- copy(as.data.table(DF))
test[d.wodtke<2, neighborhood_cat := 'Low Disadvantage']
test[d.wodtke>2, neighborhood_cat := 'High Disadvantage']
test <- test[, list(cognition=mean(f_absent)), by=c('age','neighborhood_cat')]
ggplot() +
  geom_line(data=test[age>=3,],
            aes(x=age,
                y=cognition,
                color=neighborhood_cat)) + 
  theme_minimal()
test <- make_svy_ci('~f_absent')

## Indices: id, age, year [THESE ARE HARDCODED RIGHT NOW IN OUR GFIT]
## DV = c_ppvt
## Time-invariant IVs: f_age, f_immigrant, m_age, c_sex, c_lbw, m_living_parents_15, m_immigrant, m_race, m_cognitive
## Time-variant IVs: m_relation, m_kids, m_health, m_edu, m_poverty, m_f_in_jail, m_evicted, m_anxiety, m_depression, m_job_type
## Censoring variable: censor

## Save imputed data
saveRDS(DF, './ff_clean_imputed.RDS')
DF <- readRDS('./ff_clean_imputed.RDS')

## For now, try a lot of simple binary indicators.
DF <- DF %>% mutate(bin_pov = as.numeric(m_poverty %in% c('perc_0_49','perc_50_99')))
DF <- DF %>% mutate(bin_unemployed = as.numeric(m_job_type %in% 'unemployed'))
DF <- DF %>% mutate(bin_married = as.numeric(m_relation %in% 'married'))
DF <- DF %>% mutate(bin_lesshs = as.numeric(m_edu %in% c('less_hs','hs')))
DF <- DF %>% mutate(f_absent = as.numeric(m_relation %in% 'no_contact'))
DF <- DF %>% mutate(m_poorhealth = as.numeric(m_health %in% c('Poor','Fair')))
DF <- as.data.table(DF)
ppvt_mean <- mean(DF[age>=3, c_ppvt], na.rm=T)
ppvt_sd <- sd(DF[age>=3, c_ppvt], na.rm=T)
DF[, original_ppvt := c_ppvt]
DF[, c_ppvt := (c_ppvt - ppvt_mean) / ppvt_sd]

## Lag it all.
tc_vars <- c('f_age', 'f_immigrant', 'm_age', 'm_living_parents_15', 'm_immigrant', 'm_cognitive', 'm_race','bin_lesshs')
tc_vars_child <- c('c_sex_male', 'c_lbw')
tv_vars <- c('c_ppvt','bin_unemployed','bin_pov','m_f_in_jail','m_depression','f_absent','m_poorhealth','wodtke')
DF <- as.data.table(DF)
DF[, (paste0('l.', tv_vars)) :=  shift(.SD), by='id', .SDcols=tv_vars]

## Create rolling averages within individuals for "accumulating disadvantage" variables,
## where we expect duration to be important (Wodtke calls these "duration-weighted" exposures).
## Obvious examples: exposure to poor neighborhood, exposure to extreme poverty, etc. 
## If you live in extreme poverty, its effect is not likely to disappear the second you leave poverty.
d_vars <- c('wodtke','bin_pov','m_f_in_jail')
DF[, (paste0('d.', d_vars)) :=  cumsum(.SD)/(1:.N), by='id', .SDcols=d_vars]
#DF <- DF %>% mutate(d.wodtke = ave(wodtke, id, FUN = function(x) cumsum(x) / seq_along(x)))

## PLOTS
DF[m_race=='white', m_racename := 'White']
DF[m_race=='black', m_racename := 'Black']
DF %>% filter(m_racename %in% c('White','Black')) %>% group_by(m_racename, id) %>% summarise(wodtke=mean(d.wodtke)) %>% 
  ggplot(aes(x=wodtke, fill=m_racename)) +
  geom_density(alpha=0.4) + 
  labs(x='Wodtke index of neighborhood disadvantage (duration-weighted)',y='Proportion',fill='Mother race') +
  geom_vline(xintercept = c(-2,2), size=3) + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) 

# build up the model structure base
# baseForm <- as.formula(paste0('~ ', paste(tc_vars, collapse = ' + '), ' + 
#                               (age) * (', paste(paste0('l.', tv_vars), collapse = ' + '), ')')) 
baseForm <- as.formula(paste0('~ ', paste(tc_vars, collapse = ' + '))) 
DF[, factor_age := as.character(age)]

# formula for models
formulas <- list(
  update(baseForm, c_ppvt ~ . + factor_age * (d.wodtke + d.bin_pov + l.bin_unemployed + d.m_f_in_jail + l.m_depression + l.f_absent) + c_sex_male + c_lbw),
  update(baseForm, f_absent ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_pov ~ . + factor_age * (l.wodtke + l.bin_pov  + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)), 
  update(baseForm, m_f_in_jail ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, m_depression ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_unemployed ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, wodtke ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent))
  # as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
  #                   (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
  )
formulas <- list(
  update(baseForm, c_ppvt ~ . + factor_age * d.wodtke),
  update(baseForm, f_absent ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_pov ~ . + factor_age * (l.wodtke + l.bin_pov  + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)), 
  update(baseForm, m_f_in_jail ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, m_depression ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_unemployed ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, wodtke ~ . + factor_age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent))
  # as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
  #                   (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
  )
formulas <- list(
  update(baseForm, c_ppvt ~ . + age * (l.c_ppvt + d.wodtke + d.bin_pov + l.bin_unemployed + d.m_f_in_jail + l.m_depression + l.f_absent) + c_sex_male + c_lbw),
  update(baseForm, f_absent ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  #update(baseForm, bin_lesshs ~ . + l.wodtke + l.bin_pov + l.bin_lesshs + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent),
  update(baseForm, bin_pov ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)), 
  update(baseForm, m_f_in_jail ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, m_depression ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_unemployed ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, wodtke ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  #update(baseForm, m_evicted ~ . + age * (l.wodtke + l.bin_pov + l.bin_lesshs + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent + l.m_evicted)),
  as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
                    (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
  )
formulas <- list(
  update(baseForm, c_ppvt ~ . + age * (d.wodtke + d.bin_pov + l.bin_unemployed + d.m_f_in_jail + l.m_depression + l.f_absent) + c_sex_male + c_lbw),
  update(baseForm, f_absent ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_pov ~ . + age * (l.wodtke + l.bin_pov  + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)), 
  update(baseForm, m_f_in_jail ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, m_depression ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, bin_unemployed ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  update(baseForm, wodtke ~ . + age * (l.wodtke + l.bin_pov + l.bin_unemployed + l.m_f_in_jail + l.m_depression + l.f_absent)),
  as.formula(paste0('censor ~ ', paste(tc_vars, collapse = ' + '), ' + 
                    (age) * (', paste(tv_vars, collapse = ' + '), ')')) 
  )
families <- list(gaussian, quasibinomial, quasibinomial, quasibinomial, quasibinomial, quasibinomial, gaussian) # families for models
functions <- list(svyglm, svyglm, svyglm, svyglm, svyglm, svyglm, svyglm) # functions for results
if(!(length(formulas) == length(functions) & length(formulas) == length(families))) message('FORMULAS, FAMILIES, FUNCTIONS MISALIGNED')

## RESCALE VARIABLES TO MEAN,SD of 0,1.
DF[, years_pov := d.bin_pov * real_age] ## How many years of life up until that point living in poverty? Versus l.bin_pov, which is was the individual just living in poverty last year?
#scale_vars <- c('c_ppvt','l.wodtke','d.wodtke','years_pov')
scale_vars <- 'c_ppvt_sd'
#DF <- DF %>% mutate(c_ppvt_sd = c_ppvt)
#DF <- as.data.table(DF %>% filter(age>=3) %>% mutate_at(funs(scale(.) %>% as.vector), .vars=scale_vars))
s_ppvt_mean <- svymean(DF[, c_ppvt], ff.design)
s_ppvt_sd <- sqrt(svyvar(DF[, c_ppvt], ff.design)[1])
ppvt_mean <- mean(DF[, c_ppvt])
ppvt_sd <- sd(DF[, c_ppvt])
DF[, s_normal_cog := (c_ppvt - s_ppvt_mean) / s_ppvt_sd]
DF[, c_normal_cog := (c_ppvt - ppvt_mean) / ppvt_sd]

## Test survey package
fs <- c_ppvt ~ f_age + f_immigrant + m_age + m_living_parents_15 + 
  m_immigrant + m_cognitive + m_race + age + d.wodtke + d.bin_pov + 
  l.bin_lesshs + l.bin_unemployed + d.m_f_in_jail + l.m_depression + 
  l.bin_married + c_sex_male + c_lbw + age:d.wodtke + age:d.bin_pov + 
  age:l.bin_lesshs + age:l.bin_unemployed + age:d.m_f_in_jail + 
  age:l.m_depression + age:l.bin_married

DF[, int_wodtke_bad := 2]
DF[, int_wodtke_good := -2]
d_vars <- c('int_wodtke_bad','int_wodtke_good')
DF[, (paste0('d.',d_vars)) :=  cumsum(.SD)/(1:.N), by='id', .SDcols=d_vars]

ff.design <- svydesign(id=~id, weights=~mweight, data=DF)
i <- 1
m1 <- svyglm(formulas[[i]],
             family=families[[i]], design=ff.design, data=DF)
summary(m1)
cor(DF[, c('m_depression','c_ppvt')])
sim <- rbinom(nrow(DF), 1, predict(m1, newdata=DF, type="response"))
## Validated models:
##  bin_pov (good)
##  bin_married (not really affected by anything time-varying...)
## Not validated models: bin_lesshs, m_f_in_jail, m_depression, bin_unemployed, wodtke

## VALIDATE MODELS
m_num <- 4
m <- gfit.init(formulas[m_num], families[m_num], functions[m_num], data=DF)
(m[[1]]$coefficients)

## Test survey-weighted random effects
library(lme4)
f <- formulas[m_num][[1]]
m1 <- glm(f, data=DF)
f <- update(f, c_ppvt ~ . + (1 | id))
m2 <- glmer(f, data=DF)
#m3 <- glmer(f, data=DF, weights = )
m1_coefs <- as.data.table(coef(m1))
m2_coefs <- as.data.table(coef(m2)$`id`)
m2_coefs <- m2_coefs[, lapply(.SD, mean, na.rm=T)]
m2_coefs <- melt(m2_coefs, value.name = 're')
m2_coefs[, no_re := m1_coefs$V1]
DF[, pred_re := predict(m2, DF)]
DF[, pred_nore := predict(m1, DF)]

# Here are our lags and duration-weighted variables that we want to keep updated throughout the simulation, after predicting new
lags <- list(
  
  ## Make age a factor (to get transition probabilities correct for each irregular age step)
  factor_age = function(DF) DF %>% mutate(factor_age=as.character(age)),
  
  ## Lags
  l.c_ppvt = function(DF) DF %>% mutate(l.c_ppvt=lag(c_ppvt)),
  l.bin_pov = function(DF) DF %>% mutate(l.bin_pov=lag(bin_pov)),
  l.m_f_in_jail = function(DF) DF %>% mutate(l.m_f_in_jail=lag(m_f_in_jail)),
  l.m_depression = function(DF) DF %>% mutate(l.m_depression=lag(m_depression)),
  l.bin_unemployed = function(DF) DF %>% mutate(l.bin_unemployed=lag(bin_unemployed)),
  l.f_absent = function(DF) DF %>% mutate(l.f_absent=lag(f_absent)),
  l.m_poorhealth = function(DF) DF %>% mutate(l.m_poorhealth=lag(m_poorhealth)),
  l.wodtke = function(DF) DF %>% mutate(l.wodtke=lag(wodtke)),
  
  ## Duration-weighted
  d.wodtke = function(DF) DF %>% mutate(d.wodtke = ave(wodtke, id, FUN = function(x) cumsum(x) / seq_along(x))),
  d.m_f_in_jail = function(DF) DF %>% mutate(d.m_f_in_jail = ave(m_f_in_jail, id, FUN = function(x) cumsum(x) / seq_along(x))),
  d.bin_pov = function(DF) DF %>% mutate(d.bin_pov = ave(bin_pov, id, FUN = function(x) cumsum(x) / seq_along(x)))
  #years_pov = function(DF) DF %>% mutate(years_pov = d.bin_pov * real_age) ## How many years of life up until that point living in poverty? Versus l.bin_pov, which is was the individual just living in poverty last year?
  
)

# Here are the natural deterministic and probabilistic rules.
natural_rules <- list(
  # Deterministic
  #wodtke = function(DF, ...) DF %>% mutate(wodtke = replace(wodtke, m_race=='black', -1.5)) %>% select(wodtke) 
  # None for now, but need to add "accumulating disadvantage" variables, like how many years lived in poor neighborhood, years in poverty, etc.
  # Probabilistic 
  c_ppvt = function(DF, models, ...) simPredict(DF, models, 1),
  f_absent = function(DF, models, ...) simPredict(DF, models, 2),
  #bin_lesshs = function(DF, models, ...) simPredict(DF, models, 3),
  bin_pov = function(DF, models, ...) simPredict(DF, models, 3),
  m_f_in_jail = function(DF, models, ...) simPredict(DF, models, 4),
  m_depression = function(DF, models, ...) simPredict(DF, models, 5),
  bin_unemployed = function(DF, models, ...) simPredict(DF, models, 6),
  wodtke = function(DF, models, ...) simPredict(DF, models, 7)
  #m_evicted = function(DF, models, ...) simPredict(DF, models, 9),
  #m_poorhealth = function(DF, models, ...) simPredict(DF, models, 10),
  #censor = function(DF, models, ...) simPredict(DF, models, 8)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
# Let's test drawing neighborhood disadvantage for the black population from the natural course for the white population.
# What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
intervention_rules_good <- list(
  #wodtke = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 8)
  #wodtke = function(DF, ...) DF %>% mutate(wodtke = replace(wodtke, m_race=='black', -2.5)) %>% select(wodtke) 
  #bin_married = function(DF, ...) DF %>% mutate(bin_married = replace(bin_married, bin_married==0, 1)) %>% select(bin_married)
  wodtke = function(DF, ...) DF %>% mutate(wodtke = -2) %>% select(wodtke)
)
intervention_rules_bad <- list(
  wodtke = function(DF, ...) DF %>% mutate(wodtke = 2) %>% select(wodtke)
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
  f_absent = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 2),
  #bin_lesshs = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  bin_pov = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  m_f_in_jail = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 4),
  m_depression = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 5),
  bin_unemployed = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 6),
  #m_evicted = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 9),
  #m_poorhealth = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 10),
  ## Direct pathway
  wodtke = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 7)
  ## Censoring 
  #censor = function(DF, models, ...) simPredict(DF, models, 8)
)

pov_indirect_effect_rules <- list(
  ## DV of interest
  c_ppvt = function(DF, models, ...) simPredict(DF, models, 1),
  ## Indirect pathways
  f_absent = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 2),
  #bin_lesshs = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  bin_pov = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 3),
  m_f_in_jail = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 4),
  m_depression = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 5),
  bin_unemployed = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 6),
  #m_evicted = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 9),
  #m_poorhealth = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 10),
  ## Direct pathway
  wodtke = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 7)
  ## Censoring 
  #censor = function(DF, models, ...) simPredict(DF, models, 8)
)

f_absent_indirect_effect_rules <- list(
  ## DV of interest
  c_ppvt = function(DF, models, ...) simPredict(DF, models, 1),
  ## Indirect pathways
  f_absent = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 2),
  #bin_lesshs = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  bin_pov = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 3),
  m_f_in_jail = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 4),
  m_depression = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 5),
  bin_unemployed = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 6),
  #m_evicted = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 9),
  #m_poorhealth = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 10),
  ## Direct pathway
  wodtke = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, natural_DF, 7)
  ## Censoring 
  #censor = function(DF, models, ...) simPredict(DF, models, 8)
)

boots <- 5 # Number of bootstraps, 100 takes a while
replicationSize <- 5 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
# maxit=1000 (might need to increase max iterations on multinomial models, I think most are just stopping at 100 because it is the default,
# not because they converge)

set.seed(80085)
repo <- 'E:/Data/Fragile_Families/11062018/'
setwd(repo)


cor(DF %>% filter(m_race=='white') %>% select(c_ppvt, wodtke), use='complete.obs')

#if(!file.exists("./ff_bootruns.Rds")) {

# bootruns <- lapply(1:boots, function(b) {
#   
#   message(paste0('Bootstrap ', b))
#   
#   # Sample individuals with replacement not rows
#   sampleDF_black <- DF %>%
#     filter(m_race=='black') %>%
#     select(id) %>%
#     unique %>%
#     sample_frac(replace=TRUE) %>%
#     left_join(DF)
#   
#   # Fit the model
#   gfitboot_black <- gfit.init(formulas, families, functions, data=sampleDF_black)
#   
#   # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
#   mcDF_black <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_black)) 
#   
#   # Run the "natural course" rules
#   natural_course_DF_black <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black)
#   
#   # Sample individuals with replacement not rows
#   sampleDF_white <- DF %>%
#     filter(m_race=='white') %>%
#     select(id) %>%
#     unique %>%
#     sample_frac(replace=TRUE) %>%
#     left_join(DF)
#   
#   # Fit the model
#   gfitboot_white <- gfit.init(formulas, families, functions, data=sampleDF_white)
#   
#   # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
#   mcDF_white <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_white))
#   
#   # Run the "natural course" rules
#   natural_course_DF_white <- progressSimulation(mcDF_white, lags, natural_rules, gfitboot_white)
#   
#   # Run the "intervention course" rules for BLACK population: draw neighborhood disadvantage from WHITE natural course.
#   # What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
#   intervention_DF_black <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black, intervention_rules, natural_DF=natural_course_DF_white)
#   # intervention_DF <- progressSimulation(mcDF_black, lags, natural_rules, gfitboot_black, intervention_rules)
#   
#   # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
#   direct_effect_DF <- progressSimulation(mcDF_black, lags, direct_effect_rules, gfitboot_black, natural_DF=natural_course_DF_black, intervention_DF=intervention_DF_black)
#   
#   # Simulate indirect effects [not tested]
#   # indirect_job_effect_DF <- progressSimulation(mcDF, lags, indirect_job_effect_rules, gfitboot, natural_DF=natural_DF, intervention_DF=intervention_DF)
#   
#   # Return all courses simulated
#   list(natural_black=natural_course_DF_black,
#        natural_white=natural_course_DF_white,
#        intervention=intervention_DF_black,
#        direct=direct_effect_DF)
#   
# })
# 
# saveRDS(bootruns, "./ff_bootruns_intervention.Rds")

#}

# source("./gfit.R")
# bootruns <- lapply(1:boots, function(b) {
#   
#   message(paste0('Bootstrap ', b))
#   
#   # Sample individuals with replacement not rows
#   sampleDF_all <- DF %>%
#     #filter(m_race=='black') %>%
#     select(id) %>%
#     unique %>%
#     sample_frac(replace=TRUE) %>%
#     left_join(DF)
#   
#   ## For testing a single bootstrap, just use the whole dataset so that I can see where the point estimate would be.
#   #sampleDF_all <- DF
#   
#   # Fit the model
#   gfitboot_all <- gfit.init(formulas, families, functions, data=sampleDF_all, survey_design=ff.design)
#   
#   # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
#   mcDF_all <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_all)) 
#   
#   # Run the "natural course" rules
#   natural_course_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all)
#   
#   # Run the "intervention course" rules for BLACK population: draw neighborhood disadvantage from WHITE natural course.
#   # What is the direct effect of neighborhood disparities on the total black-white disparity in the outcome?
#   # intervention_DF <- progressSimulation(mcDF_black, lags, rules, gfitboot_black, intervention_rules, natural_DF=natural_course_DF_white)
#   # intervention_DF_black <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules)
#   intervention_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules)
#   
#   # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
#   direct_effect_DF_all <- progressSimulation(mcDF_all, lags, direct_effect_rules, gfitboot_all, natural_DF=natural_course_DF_all, intervention_DF=intervention_DF_all)
#   
#   # Simulate indirect effects [not tested]
#   indirect_effect_DF_all <- progressSimulation(mcDF_all, lags, indirect_effect_rules, gfitboot_all, natural_DF=natural_course_DF_all, intervention_DF=intervention_DF_all)
#   
#   # Return all courses simulated
#   list(natural=natural_course_DF_all,
#        intervention=intervention_DF_all,
#        direct=direct_effect_DF_all,
#        indirect=indirect_effect_DF_all)
#   
# })
# 
# saveRDS(bootruns, "./ff_bootruns_intervention.Rds")

#}

source("E:/Data/Fragile_Families/gfit.R")
boots <- 10
replicationSize <- 1
ff.design <- svydesign(id=~id, weights=~mweight, data=DF)
bootruns <- lapply(1:boots, function(b) {
  
  message(paste0('Bootstrap ', b))
  
  # Sample individuals with replacement, not rows (this is how we propagate model error right now).
  sampleDF_all <- DF %>%
    select(id) %>%
    unique %>%
    sample_frac(replace=TRUE) %>%
    left_join(DF)
  
  ## For testing a single bootstrap, just use the whole dataset so that I can see where the point estimate would be.
  # sampleDF_all <- DF
  
  # Fit the model
  gfitboot_all <- gfit.init.survey(formulas, families, functions, data=sampleDF_all, survey_design=ff.design)

  # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
  mcDF_all <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_all)) 
  
  ## SURVEY CHANGE: Instead of simulating on a weighted dataset, we have to simulate on an unweighted nationally representative dataset. 
  ## The reason is the simScenario() gets all fucked up because weights in the observed data are individual-specific, but simScenario()
  ## draws from the provided courseDF randomly. We can either fix it this way so that random draws are fine, or don't draw at all and just
  ## use the exact values for each individual from the provided courseDF (I still can't think of a reason why this would be a bad thing to do).
  ## We can use sample_frac() with the survey weight to get an unweighted, nationally representative dataset of the same size as the survey dataset.
  # mcDF_all <- mcDF_all %>%
  #   filter(year==min(year)) %>%
  #   select(id,mweight) %>%
  #   unique %>%
  #   sample_frac(size=1, weight=mweight, replace=TRUE) %>%
  #   select(id) %>%
  #   left_join(mcDF_all)
  mcDF_all <- mcDF_all %>%
    filter(age==min(age)) %>%
    select(id,mweight) %>%
    unique %>%
    sample_frac(size=1, weight=mweight, replace=TRUE) %>%
    select(id) %>%
    left_join(mcDF_all)

  # Run the "natural course" rules
  message('NATURAL COURSE')
  natural_course_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all)

  # Intervention: good neighborhood (-2) vs. bad neighborhood (2)
  message('INTERVENTION COURSE 1')
  intervention_good <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules_good)
  message('INTERVENTION COURSE 2')
  intervention_bad <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules_bad)
  
  # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
  # CHANGE: do effects on difference between good and bad interventions, rather than between good intervention and natural course.
  message('DIRECT COURSE')
  direct_effect_DF_all <- progressSimulation(mcDF_all, lags, direct_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good)
  
  # Simulate indirect effects [not tested]
  message('INDIRECT COURSE 1')
  indirect_effect_pov <- progressSimulation(mcDF_all, lags, pov_indirect_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good)
  message('INDIRECT COURSE 2')
  indirect_effect_f_absent <- progressSimulation(mcDF_all, lags, f_absent_indirect_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good)
  
  # Return all courses simulated
  list(natural=natural_course_DF_all,
       intervention_good=intervention_good,
       intervention_bad=intervention_bad,
       direct=direct_effect_DF_all,
       indirect_effect_pov=indirect_effect_pov,
       indirect_effect_f_absent=indirect_effect_f_absent)
  
})

saveRDS(bootruns, "./ff_bootruns_intervention_survey_v3.Rds")

#bootruns <- read_rds("./bootruns.Rds")
#for(i in 1:boots) names(bootruns[[i]])[4] <- 'indirect'

## Make effect size table
all_effects <- lapply(c('natural','intervention_good','intervention_bad','direct','indirect_effect_pov','indirect_effect_f_absent'), compile_sims_table, bs=bootruns)
all_effects <- Reduce(merge, all_effects)
#for(c in c('mppvt','married','lesshs','pov','jail','depression','unemployed','censor','mwodtke')) {
for(c in c('mppvt')) {
  for(b in 1:boots) all_effects[, (paste0(c, '_total_', b)) := get(paste0(c, '_intervention_good_', b)) - get(paste0(c, '_intervention_bad_', b))]
  for(b in 1:boots) all_effects[, (paste0(c, '_direct_', b)) := get(paste0(c, '_direct_', b)) - get(paste0(c, '_intervention_bad_', b))]
  for(b in 1:boots) all_effects[, (paste0(c, '_pov_indirect_', b)) := get(paste0(c, '_indirect_effect_pov_', b)) - get(paste0(c, '_intervention_bad_', b))]
  for(b in 1:boots) all_effects[, (paste0(c, '_f_absent_indirect_', b)) := get(paste0(c, '_indirect_effect_f_absent_', b)) - get(paste0(c, '_intervention_bad_', b))]
  for(b in 1:boots) all_effects[, (paste0(c, '_p_direct_', b)) := get(paste0(c, '_direct_', b)) / get(paste0(c, '_total_', b)) * 100]
  for(b in 1:boots) all_effects[, (paste0(c, '_p_pov_indirect_', b)) := get(paste0(c, '_pov_indirect_', b)) / get(paste0(c, '_total_', b)) * 100]
  for(b in 1:boots) all_effects[, (paste0(c, '_p_f_absent_indirect_', b)) := get(paste0(c, '_f_absent_indirect_', b)) / get(paste0(c, '_total_', b)) * 100]
}
for(c in c('total','direct','pov_indirect','f_absent_indirect','p_direct','p_pov_indirect','p_f_absent_indirect')) {
  all_effects[, (paste0(c,'_lower')) := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep(paste0('^mppvt_', c), names(all_effects))]
  all_effects[, (paste0(c,'_mean')) := apply(.SD, 1, mean), .SDcols=grep(paste0('^mppvt_', c), names(all_effects))]
  all_effects[, (paste0(c,'_upper')) := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep(paste0('^mppvt_', c), names(all_effects))]
  all_effects[, (c) := paste0(round(get((paste0(c,'_mean'))),3), ' (', round(get((paste0(c,'_lower'))),3), ' - ', round(get((paste0(c,'_upper'))),3), ')')]
  
}
all_effects <- all_effects[, c('age','total','direct','pov_indirect','f_absent_indirect','p_direct','p_pov_indirect','p_f_absent_indirect')]
all_effects <- redo_ages(all_effects)
saveRDS(all_effects, paste0(out_dir, '/all_effects.RDS'))

## Make coefficient table with full DF
gfitboot_all <- gfit.init.survey(formulas, families, functions, data=DF, survey_design=ff.design)
saveRDS(gfitboot_all, paste0(out_dir, '/models.RDS'))

## Aggregate
all_runs <- rbindlist(lapply(c('natural','intervention_good','intervention_bad','direct','indirect_effect_pov','indirect_effect_f_absent'), compile_sims_simple, bs = bootruns))
all_runs <- redo_ages(all_runs)
# actualDF <- DF %>%
#   group_by(age) %>%
#   summarize(
#     mppvt = weighted.mean(c_ppvt, mweight), 
#     fabsent = weighted.mean(f_absent, mweight),
#     pov = weighted.mean(bin_pov, mweight),
#     jail = weighted.mean(m_f_in_jail, mweight),
#     depression = weighted.mean(m_depression, mweight),
#     unemployed = weighted.mean(bin_unemployed, mweight),
#     censor = weighted.mean(censor, mweight),
#     mwodtke = weighted.mean(wodtke, mweight),
#     mmaterial = weighted.mean(m_material, mweight),
#     mstress = weighted.mean(m_stress, mweight)
#   ) %>%
#   gather("measure", "mean", -age)
# actualDF <- redo_ages(actualDF)

make_svy_ci <- function(v) {
  cis <- svyby(as.formula(v), ~age, ff.design, svymean, ci=T)
  cis <- as.data.table(cis)
  setnames(cis, gsub('~','',v), 'mean')
  cis[, upr := mean + 1.96*se]
  cis[, lwr := mean - 1.96*se]
  cis[, measure := gsub('~','',v)]
  return(cis)
}
all_cis <- rbindlist(lapply(paste0('~', tv_vars), make_svy_ci))
all_cis <- redo_ages(all_cis)
all_cis[measure=='c_ppvt', measure := 'mppvt']
all_cis[measure=='wodtke', measure := 'mwodtke']
all_cis[measure=='bin_pov', measure := 'pov']
all_cis[measure=='f_absent', measure := 'fabsent']
all_cis[measure=='m_f_in_jail', measure := 'jail']
all_cis[measure=='bin_unemployed', measure := 'unemployed']
all_cis[measure=='m_depression', measure := 'depression']

## Plot 
good_col <- '#1f78b4'
bad_col <- '#e31a1c'
nat_col <- '#33a02c'
direct_col <- '#b2df8a'
cov_plot <- all_runs %>%
  filter(!(measure %in% c('censor','lesshs','years_pov','years'))) %>%
  mutate(measure=factor(measure, levels = c('mppvt','mwodtke','pov','fabsent','jail','unemployed','depression'))) %>%
  mutate(type=factor(type,levels=c('natural','intervention_good','intervention_bad','direct','indirect_effect_pov','indirect_effect_marriage'))) %>%
  as.data.table()
cov_plot[age<3 & measure=='mppvt', mean := NA]
cov_plot[age<3 & measure=='mppvt', lwr := NA]
cov_plot[age<3 & measure=='mppvt', upr := NA]
cov_plot[type=='intervention_good', type := 'Low Disadvantage']
cov_plot[type=='intervention_bad', type := 'High Disadvantage']
cov_plot[type=='direct', type := 'Direct effect']

cov_obs <- all_cis %>%
  mutate(type='observed') %>%
  filter(!(measure %in% c('censor','lesshs','years_pov','years'))) %>%
  mutate(measure=factor(measure, levels = c('mppvt','mwodtke','pov','fabsent','jail','unemployed','depression'))) %>%
  as.data.table()
cov_obs[age<3 & measure=='mppvt', mean := NA]
cov_obs[age<3 & measure=='mppvt', lwr := NA]
cov_obs[age<3 & measure=='mppvt', upr := NA]

## Make clean names/labels
clean_names <- data.table(measure=c('mppvt','mwodtke','pov','mmaterial','mstress','fabsent','jail','unemployed','depression'),
                          clean_measure=c('PPVT Score','Wodtke Index','Poverty','Material stress','Parenting stress','No contact father','Father incarcerated','Mother unemployed','Mother depressed'))
clean_names[, clean_measure := factor(clean_measure, levels=c('PPVT Score','Wodtke Index','Poverty','Material stress','Parenting stress','No contact father','Father incarcerated','Mother unemployed','Mother depressed'))]
cov_plot <- merge(cov_plot, clean_names, by='measure')
cov_obs <- merge(cov_obs, clean_names, by='measure')
## Figure 2: Natural course for all time-varying variables with data.
fig2_covs <- ggplot(data=cov_plot[type=='natural' & measure != 'dpov',], aes(x=age, y=mean)) +
  geom_line(color=nat_col, size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill=nat_col, alpha=.5) +
  geom_point(data=cov_obs[measure != 'dpov',], shape=21, fill=nat_col, color='black', size=5) + 
  geom_errorbar(data=cov_obs[measure != 'dpov',], aes(ymin=lwr,ymax=upr)) + 
  theme_bw() + 
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  facet_wrap(~clean_measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated course'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')
fig2_ppvt <- ggplot(data=cov_plot[type=='natural' & measure == 'mppvt' & age >= 3,], aes(x=age, y=mean)) +
  geom_line(color=nat_col, size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill=nat_col, alpha=.5) +
  geom_point(data=cov_obs[measure == 'mppvt' & age >= 3,], shape=21, fill=nat_col, color='black', size=5) + 
  geom_errorbar(data=cov_obs[measure == 'mppvt',], aes(ymin=lwr,ymax=upr)) + 
  theme_bw() +
  guides(fill=guide_legend(title='Simulated course'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')

## Figure 3: Intervention courses with direct course.
cols <- c('Low Disadvantage' = '#1f78b4', 'High Disadvantage' = '#e31a1c', 'Direct effect' = 'black', 'natural' = '#33a02c', 'indirect_effect_pov' = 'black', 'indirect_effect_marriage' = 'black')
cov_plot_courses <- copy(cov_plot)
# cov_plot_courses <- cov_plot_courses[type=='direct' & measure != 'mppvt', mean := NA]
# cov_plot_courses <- cov_plot_courses[type=='direct' & measure != 'mppvt', upr := NA]
fig3_covs <- ggplot(data=cov_plot_courses %>% filter(type %in% c('Low Disadvantage','High Disadvantage') & measure != 'dpov'), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols) + 
  #xlim(c(1,9)) + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  facet_wrap(~clean_measure, scales = 'free_y') +
  guides(fill=FALSE, color=guide_legend(title='Simulated course',override.aes = list(size=10))) + labs(x='Child age',y='Population average (proportion or mean)')
fig3_ppvt <- ggplot(data=cov_plot_courses %>% filter(type %in% c('Low Disadvantage','High Disadvantage') & measure == 'mppvt'), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols) + 
  xlim(c(3,9)) + 
  theme_bw() +
  guides(fill=FALSE, color=guide_legend(title='Simulated course',override.aes = list(size=10))) + labs(x='Child age',y='Population average (proportion or mean)')

fig4_covs <- ggplot(data=cov_plot_courses %>% filter(type %in% c('Low Disadvantage','High Disadvantage','Direct effect') & measure != 'dpov'), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols) + 
  #xlim(c(1,9)) + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  facet_wrap(~clean_measure, scales = 'free_y') +
  guides(fill=FALSE, color=guide_legend(title='Simulated course',override.aes = list(size=10))) + labs(x='Child age',y='Population average (proportion or mean)')
