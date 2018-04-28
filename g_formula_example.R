
## Code by Maarten J. Bijlsma, version April 20th 2018
## Max Planck Institute for Demographic Research
## Bijlsma@demogr.mpg.de

## An example application of the parametric g-formula
## using an artificially created fertility dataset

setwd('C:/Users/ngraetz/Documents/repos/g_formula')

# Load the data (see supplemental material)
ex.dat <- read.csv(file='ex.csv',head=T)

# First a description of the 'toy' example data:
names(ex.dat)
dim(ex.dat)
length(unique(ex.dat$id))

# the data is in the long format: we have multiple observations 
# per individual (id), and every observation is a row

# Variables are
# id: the identification number
# year: calendar year of observation
# age: years since birth of respondent
# age1522, age2329, age30pl are dummy variables for age category

# baseline variables (not time varying):
# region: regions's A, B, and C contain dummy indicators for region of origin
# parental.income: low, med (medium) and high contain dummy indicators for income level of parents

# time varying variables:
# birth: whether the respondent had a child in that year
# totalbirth: sum of births for that respondent up to that point in time
# job: full, edu, and other contain dummy variables for main activity of respondent
# totaledu: sum of years respondent has been in education during follow-up
# partner: single, cohab (cohabitation), and married contain dummy variables for partnership status
# l. indicates a 1-year lag of a time varying covariate
# cens is an indicator for censoring that occurred before 2008

# in our paper, the number of variables is larger
# but more variables are not needed to demonstrate the fundamentals
# of the application of the g-formula in R

# we want to model the relationships between these time-varying covariates
# in this example, all of them are either binary (birth) or
# multinomial (job and partner)
# I choose to model multinomial variables with sets of logistic regressions
# so I load predict functions that make multinomial predictions:

multinomial.predict.3var <- function(predict.x1,predict.x2,predict.x3,
                                     dataset) {
  
  # predict probabilities for each of the 3 outcome possibilities
  x1 <- predict(predict.x1,dataset,type='response')
  x2 <- predict(predict.x2,dataset,type='response')
  x3 <- predict(predict.x3,dataset,type='response')
  
  # the total of the three probabilities should be 1
  # but sometimes deviates a little from 1 due to
  # the numerical integration R uses, so we make sure it is 1:
  rescaler <- 1/(x1+x2+x3)
  
  # rescale constituent probabilities
  x1r <- x1*rescaler
  x2r <- x2*rescaler
  # x3r does not need to be determined
  
  # determine boundary values
  x1b <- x1r
  x2b <- x1r + x2r
  # x3b does not need to be determined
  
  # determine n
  n2 <- length(x1)
  
  # draw random variables
  decider <- runif(n2,0,1)
  
  # determine category
  x1c <- ifelse(decider <= x1b,1,0)
  x2c <- ifelse(decider > x1b & decider <= x2b,1,0)
  x3c <- ifelse(decider > x2b,1,0)
  
  return(cbind(x1c,x2c,x3c))
}

binomial.predict <- function(predict.x1,dataset) {
  
  x1 <- predict(predict.x1,dataset,type='response')
  
  x1c <- rbinom(length(x1),1,x1)
  
  return(x1c)
}

# Due to having repeated measurements per individual
# we want to bootstrap individuals and not individual observations
# so we also need a function that does that:

long.sample <- function(originaldata, originaldataid) {
  # select a bunch of IDs
  IDs <- unique(originaldataid)
  y <- sample(IDs,length(IDs),replace=T)
  z <- table(table(y))
  
  # from there, select a group once
  selectID <- sample(IDs,size=z[1],replace=F)
  newdata <- originaldata[which(originaldataid %in% selectID),]
  
  if(length(z) > 1) {
    
    for(i in 2:length(z)) {
      
      # select a new group of IDs that was not yet selected
      IDs2 <- setdiff(IDs,selectID)
      
      # from there, randomly select a group of people of the right size
      selectID2 <- sample(IDs2,size=z[i],replace=F)
      selectID <- c(selectID,selectID2) # so we don't re-select the newly
      # selected people either
      
      for(j in 1:i) {
        
        # copy the new dataset i number of times
        newdata <- rbind(newdata,originaldata[which(originaldataid %in% selectID2),])
        
      }
      
    }
    
    return(newdata)
  }
}

# Since we will also lag observations in the g-formula loop
# it is more efficient if we save the location of the time-varying columns
col.index.from <- NULL
col.index.from[1] <- grep(paste0('^','birth','$'), colnames(ex.dat))
col.index.from[2] <- grep(paste0('^','totalbirth','$'), colnames(ex.dat))
col.index.from[3] <- grep(paste0('^','job.full','$'), colnames(ex.dat))
col.index.from[4] <- grep(paste0('^','job.edu','$'), colnames(ex.dat))
col.index.from[5] <- grep(paste0('^','job.other','$'), colnames(ex.dat))
col.index.from[6] <- grep(paste0('^','totaledu','$'), colnames(ex.dat))
col.index.from[7] <- grep(paste0('^','partner.single','$'), colnames(ex.dat))
col.index.from[8] <- grep(paste0('^','partner.cohab','$'), colnames(ex.dat))
col.index.from[9] <- grep(paste0('^','partner.married','$'), colnames(ex.dat))

col.index.to <- NULL
col.index.to[1] <- grep(paste0('^','l.birth','$'), colnames(ex.dat))
col.index.to[2] <- grep(paste0('^','l.totalbirth','$'), colnames(ex.dat))
col.index.to[3] <- grep(paste0('^','l.job.full','$'), colnames(ex.dat))
col.index.to[4] <- grep(paste0('^','l.job.edu','$'), colnames(ex.dat))
col.index.to[5] <- grep(paste0('^','l.job.other','$'), colnames(ex.dat))
col.index.to[6] <- grep(paste0('^','l.totaledu','$'), colnames(ex.dat))
col.index.to[7] <- grep(paste0('^','l.partner.single','$'), colnames(ex.dat))
col.index.to[8] <- grep(paste0('^','l.partner.cohab','$'), colnames(ex.dat))
col.index.to[9] <- grep(paste0('^','l.partner.married','$'), colnames(ex.dat))

# the order of the above variables in the index.from and index.to
# must be the same!
# and of course, the length of both vectors must also be the same
length(col.index.from)
length(col.index.to)
intersect(col.index.from,col.index.to)

# we must specify the formulas for the regressions that we will use:
formula.outcome <- c("outcome ~ 
                     
                     # time constants #
                     
                     parental.income +
                     
                     region.B + # region.A is ref
                     
                     # time varying #
                     
                     (age1522 + age30pl) * # age2329 is ref
                     
                     (l.birth + 
                     l.totalbirth +
                     
                     l.partner.cohab + l.partner.married + # single is ref
                     
                     l.job.full + l.job.other + # edu is ref
                     l.totaledu) ")

# We use lagged versions of the time-varying covariates
# because then the causal (temporal) order is more certain
# but technically it is allowed to also use non-lagged terms
# Furthermore, we interact age with the time-varying covariates
# for example to allow the effect of the time-varying covariates
# to change as an individual ages

# assign formulas for each outcome:
formula.birth <- as.formula(sub('outcome','birth',formula.outcome))

formula.job.full <- as.formula(sub('outcome','job.full',formula.outcome))
formula.job.edu <- as.formula(sub('outcome','job.edu',formula.outcome))
formula.job.other <- as.formula(sub('outcome','job.other',formula.outcome))

formula.partner.married <- as.formula(sub('outcome','partner.married',formula.outcome))
formula.partner.cohab <- as.formula(sub('outcome','partner.cohab',formula.outcome))
formula.partner.single <- as.formula(sub('outcome','partner.single',formula.outcome))

# 1) We apply the same formula (that is, same list of predictors) to 
# every outcome. However, we could in theory custom-make a 
# model for each outcome
# 2) Note that we do not make models for the 'deterministically' determined
# outcomes, namely age (function of time)
# totaledu (function of time and edu, the latter of which is modelled)
# totalbirth (function of birth, which is modelled)
# since they are deterministic, they will also be taken care of
# deterministically in the g-formula, below
# 3) all of these models are for logistic or multinomial variables
# but if an outcome was count (e.g. Poisson distributed)
# or continuous (e.g. Gaussian distributed)
# the procedure would work largely the same: in the event
# of a Gaussian distribution, however, we would also want to save
# the standard deviation of the residuals of our linear regression models, 
# and include that in our predict function used in the g-formula below

# to model censoring, it makes sense to use non-lagged time-varying covariates
formula.outcome.c <- c("outcome ~ 
                       
                       # time constants #
                       
                       parental.income +
                       
                       region.B + # region.A is ref
                       
                       # time varying #
                       
                       (age1522 + age30pl) * # age2329 is ref
                       
                       (birth + 
                       totalbirth +
                       
                       partner.cohab + partner.married + # single is ref
                       
                       job.full + job.other + # edu is ref
                       totaledu)
                       
                       ")

formula.censor <- as.formula(sub('outcome','censor',formula.outcome.c))

# Finally, we want a vector, matrix or array in which we
# save the findings of our g-formula iterations:

# let's go with 250 iterations at first
bssize <- 3

# this will save info on birth and totalbirth
output.dat <- rep(NA,bssize*6*(38-15)) # row, col, 3rd dim
dim(output.dat) <- c(bssize,6,38-15)

# this will save info on parities
babytable <- babytable.2 <- babytable.3 <- rep(NA,bssize*15*(38-15))
dim(babytable) <- c(bssize,15,38-15)
dim(babytable.2) <- c(bssize,15,38-15)
dim(babytable.3) <- c(bssize,15,38-15)

# we could also make matrices like the above ones for
# our other time-varying covariates such as
# partnership, to see how this mediator changes when
# we intervene on education
# we omit this here

t1 <- Sys.time()

## let's start the actual G-formula: ##
for(bs in 1:bssize) {
  
  # sample individuals (with replacement) from ex.dat
  ex.sample <- long.sample(ex.dat,ex.dat$id)
  
  # (re)fit models to ex.sample
  # bsf = bootstrap fit
  # estimate relations
  fit.bsf.birth <- glm(formula.birth, family=binomial, data=ex.sample)
  
  fit.bsf.job.full <- glm(formula.job.full, family=binomial, data=ex.sample)
  fit.bsf.job.edu <- glm(formula.job.edu, family=binomial, data=ex.sample)
  fit.bsf.job.other <- glm(formula.job.other, family=binomial, data=ex.sample)
  
  fit.bsf.partner.married <- glm(formula.partner.married, family=binomial, data=ex.sample)
  fit.bsf.partner.cohab <- glm(formula.partner.cohab, family=binomial, data=ex.sample)
  fit.bsf.partner.single <- glm(formula.partner.single, family=binomial, data=ex.sample)
  
  fit.bsf.censor <- glm(formula.censor, family=binomial, data=ex.sample)
  # note: if the estimates from these models are desired
  # for example to construct a table similar to 'Table 2' in our paper
  # make a matrix similar to 'output.dat' where these estimates are stored
  # and then take the average over the iterations for each respective coefficient
  
  # take individuals at time 0 (in this dataset: year 1986)
  # and discard the other observations
  ex.sample <- ex.sample[ex.sample$year==1986,]
  # these observations will be taken as starting values
  
  # to reduce Monte Carlo error
  # in this example, we don't perform a Monte Carlo loop
  # but instead copy the dataset a number of times
  # (sampling variability, needed for correct standard errors
  # comes from the Bootstrap part of the procedure)
  ex.sample <- rbind(ex.sample,ex.sample,ex.sample,ex.sample,ex.sample,
                     ex.sample,ex.sample,ex.sample,ex.sample,ex.sample)
  
  ## Nick TO-DO: two options for assessing convergence on the direct/indirect effects. 1. Brute force, as currently set up: run with different numbers of copies above
  ## (5, 10, 20, 30) and show that estimates of effects converge at a certain number of copies. 2. Rewrite to actually do the MCMC loop, and keep track of when effects converge.
  ## I think we can also parallelize the shit out of this. 
  
  # since individuals were drawn with replacement
  # multiple 'individuals' in the data can have the same ID
  # therefore, give each individual a new ID, since
  # otherwise the sorting below will go wrong when
  # individuals are ordered by ID and time
  ex.sample$idnr <- 1:length(ex.sample$id)
  
  # make a copy of this, which can be used for the
  # intervention loop
  ex.sample.3 <- ex.sample.2 <- ex.sample
  
  ## NATURAL LOOP ##
  # start a loop that moves through the follow-up time units
  # this part of the g-formula tries to reproduce the empirical data
  # and is known as the 'natural course'
  # Nick TO-DO: right now just predicting t from t-1. We could make this an actual autoregressive model instead of simply lagged variables.
  # Really all of these models can change. The trick is just writing predict functions and wrapping everything up in the g-formula bootstrap/MC loop to run efficiently.
  # I think there's a ton of space to improve this.
  
  ## Nick TO-DO: remove a lot of this hard-coding around the dimensions of this specific dataset. Add as arguments to functions.
  for(t in 1987:2008) {
    
    message(paste0('Predicting natural course, ', t))
    # make a copy of the previous row and then update year and age
    # but don't make a copy if someone was censored in the previous year
    ex.sample.temp <- ex.sample[ex.sample$year==(t-1) & ex.sample$censor!=1,]
    ex.sample.temp$year <- ex.sample.temp$year+1
    ex.sample.temp$age <- ex.sample.temp$age+1
    
    # update the age categories
    ex.sample.temp$age1522 <- ifelse(ex.sample.temp$age <= 22,1,0)
    ex.sample.temp$age2329 <- ifelse(ex.sample.temp$age >= 23 & ex.sample.temp$age <= 29,1,0)
    ex.sample.temp$age30pl <- ifelse(ex.sample.temp$age >= 30,1,0)
    
    # lag values: since the new values for t have not yet been produced
    # I can actually take the column information from a row at t
    # and put it in the same row at t but then in a column meant for
    # lagged values
    # the first are the 'from' columns, and the second the 'to' columns
    # so: take values from the the 'from' column and put them in the 'to' column
    ex.sample.temp[,col.index.to] <- ex.sample.temp[,col.index.from]
    
    # join this newly created row with the old rows
    ex.sample <- rbind(ex.sample,ex.sample.temp)
    rm(ex.sample.temp)
    
    # order by ID variable
    ex.sample <- ex.sample[order(ex.sample$idnr,ex.sample$year),]
    
    ## predict the outcomes for year t, using the lagged observations
    # note: the order of predictions does not matter currently
    # because our predictons are based on lagged values
    # but if we allowed variables in year t to affect other variables
    # in year t, the order of the predictions must follow that causal ordering
    
    ## Nick TO-DO: I feel like we could just fit all the models within year t simultaneously to avoid having to care about order...
    ## But also possible there would just be identification problems.
    
    # predict birth
    ex.sample$birth[ex.sample$year==t] <- binomial.predict(fit.bsf.birth, ex.sample[ex.sample$year==t,])
    
    # update totalbirth
    # totalbirth is the totalbirth value that we had last year, plus
    # the birth that just happened (1) or not (0)
    ex.sample$totalbirth[ex.sample$year==t] <- ex.sample$l.totalbirth[ex.sample$year==t] +
      ex.sample$birth[ex.sample$year==t]
    
    # predict job (includes education)
    x <- multinomial.predict.3var(fit.bsf.job.full,fit.bsf.job.edu,fit.bsf.job.other,
                                  ex.sample[ex.sample$year==t,])
    ex.sample$job.full[ex.sample$year==t] <- x[,1]
    ex.sample$job.edu[ex.sample$year==t] <- x[,2]
    ex.sample$job.other[ex.sample$year==t] <- x[,3]
    rm(x)
    
    # update totaledu
    # if someone is in education this year, then we add one
    # to the total years spent in education
    ex.sample$totaledu[ex.sample$year==t] <- ex.sample$l.totaledu[ex.sample$year==t] +
      ex.sample$job.edu[ex.sample$year==t]
    
    # predict partnership
    x <- multinomial.predict.3var(fit.bsf.partner.married,fit.bsf.partner.cohab,fit.bsf.partner.single,
                                  ex.sample[ex.sample$year==t,])
    ex.sample$partner.married[ex.sample$year==t] <- x[,1]
    ex.sample$partner.cohab[ex.sample$year==t] <- x[,2]
    ex.sample$partner.single[ex.sample$year==t] <- x[,3]
    rm(x)
    
    # predict censoring
    # (has to be done last, since it depends on the values
    # of the (updated) covariates at time t, and not the lagged ones)
    ex.sample$censor[ex.sample$year==t] <- binomial.predict(fit.bsf.censor,ex.sample[ex.sample$year==t,])
    
  }
  
  # now let's create a counterfactual dataset for the TOTAL EFFECT
  # where we intervene on the data in some way
  # we will not allow marriage: we make these people single instead.
  
  for(t in 1987:2008) {
    
    message(paste0('Predicting counterfactual, ', t))
    # make a copy of the previous row and then update year and age
    # but don't make a copy if someone was censored in the previous year
    ex.sample.2.temp <- ex.sample.2[ex.sample.2$year==(t-1) & ex.sample.2$censor!=1,]
    ex.sample.2.temp$year <- ex.sample.2.temp$year+1
    ex.sample.2.temp$age <- ex.sample.2.temp$age+1
    
    # update the age categories
    ex.sample.2.temp$age1522 <- ifelse(ex.sample.2.temp$age <= 22,1,0)
    ex.sample.2.temp$age2329 <- ifelse(ex.sample.2.temp$age >= 23 & ex.sample.2.temp$age <= 29,1,0)
    ex.sample.2.temp$age30pl <- ifelse(ex.sample.2.temp$age >= 30,1,0)
    
    # lag values: since the new values for t have not yet been produced
    # I can actually take the column information from a row at t
    # and put it in the same row at t but then in a column meant for
    # lagged values
    # the first are the 'from' columns, and the second the 'to' columns
    # so: take values from the the 'from' column and put them in the 'to' column
    ex.sample.2.temp[,col.index.to] <- ex.sample.2.temp[,col.index.from]
    
    # join this newly created row with the old rows
    ex.sample.2 <- rbind(ex.sample.2,ex.sample.2.temp)
    rm(ex.sample.2.temp)
    
    # order by ID variable
    ex.sample.2 <- ex.sample.2[order(ex.sample.2$idnr,ex.sample.2$year),]
    
    ## predict the outcomes for year t, using the lagged observations
    # note: the order of predictions does not matter currently
    # because our predictons are based on lagged values
    # but if we allowed variables in year t to affect other variables
    # in year t, the order of the predictions must follow that causal ordering
    
    # predict birth
    ex.sample.2$birth[ex.sample.2$year==t] <- binomial.predict(fit.bsf.birth, ex.sample.2[ex.sample.2$year==t,])
    
    # update totalbirth
    # totalbirth is the totalbirth value that we had last year, plus
    # the birth that just happened (1) or not (0)
    ex.sample.2$totalbirth[ex.sample.2$year==t] <- ex.sample.2$l.totalbirth[ex.sample.2$year==t] +
      ex.sample.2$birth[ex.sample.2$year==t]
    
    # predict job (includes education)
    x <- multinomial.predict.3var(fit.bsf.job.full,fit.bsf.job.edu,fit.bsf.job.other,
                                  ex.sample.2[ex.sample.2$year==t,])
    ex.sample.2$job.full[ex.sample.2$year==t] <- x[,1]
    ex.sample.2$job.edu[ex.sample.2$year==t] <- x[,2]
    ex.sample.2$job.other[ex.sample.2$year==t] <- x[,3]
    rm(x)
    
    # update totaledu
    # if someone is in education this year, then we add one
    # to the total years spent in education
    ex.sample.2$totaledu[ex.sample.2$year==t] <- ex.sample.2$l.totaledu[ex.sample.2$year==t] +
      ex.sample.2$job.edu[ex.sample.2$year==t]
    
    # predict partnership
    x <- multinomial.predict.3var(fit.bsf.partner.married,fit.bsf.partner.cohab,fit.bsf.partner.single,
                                  ex.sample.2[ex.sample.2$year==t,])
    ex.sample.2$partner.married[ex.sample.2$year==t] <- x[,1]
    ex.sample.2$partner.cohab[ex.sample.2$year==t] <- x[,2]
    ex.sample.2$partner.single[ex.sample.2$year==t] <- x[,3]
    rm(x)
    
    # !! Intervention: don't allow marriage: make them single instead
    ex.sample.2$partner.single[ex.sample.2$partner.married==1 & ex.sample.2$year==t] <- 1
    ex.sample.2$partner.married[ex.sample.2$partner.married==1 & ex.sample.2$year==t] <- 0
    
    # predict censoring
    # (has to be done last, since it depends on the values
    # of the (updated) covariates at time t, and not the lagged ones)
    ex.sample.2$censor[ex.sample.2$year==t] <- binomial.predict(fit.bsf.censor, ex.sample.2[ex.sample.2$year==t,])
    
  }
  
  ## now let's create a counterfactual dataset for the NATURAL DIRECT EFFECT
  ## to achieve this, the mediator distribution must be taken from the intervention scenario
  # and the exposure from the natural course
  
  for(t in 1987:2008) {
    
    message(paste0('Predicting natural direct effect, ', t))
    # make a copy of the previous row and then update year and age
    # but don't make a copy if someone was censored in the previous year
    ex.sample.3.temp <- ex.sample.3[ex.sample.3$year==(t-1) & ex.sample.3$censor!=1,]
    ex.sample.3.temp$year <- ex.sample.3.temp$year+1
    ex.sample.3.temp$age <- ex.sample.3.temp$age+1
    
    # update the age categories
    ex.sample.3.temp$age1522 <- ifelse(ex.sample.3.temp$age <= 22,1,0)
    ex.sample.3.temp$age2329 <- ifelse(ex.sample.3.temp$age >= 23 & ex.sample.3.temp$age <= 29,1,0)
    ex.sample.3.temp$age30pl <- ifelse(ex.sample.3.temp$age >= 30,1,0)
    
    # lag values: since the new values for t have not yet been produced
    # I can actually take the column information from a row at t
    # and put it in the same row at t but then in a column meant for
    # lagged values
    # the first are the 'from' columns, and the second the 'to' columns
    # so: take values from the the 'from' column and put them in the 'to' column
    ex.sample.3.temp[,col.index.to] <- ex.sample.3.temp[,col.index.from]
    
    # join this newly created row with the old rows
    ex.sample.3 <- rbind(ex.sample.3,ex.sample.3.temp)
    rm(ex.sample.3.temp)
    
    # order by ID variable
    ex.sample.3 <- ex.sample.3[order(ex.sample.3$idnr,ex.sample.3$year),]
    
    ## predict the outcomes for year t, using the lagged observations
    # note: the order of predictions does not matter currently
    # because our predictons are based on lagged values
    # but if we allowed variables in year t to affect other variables
    # in year t, the order of the predictions must follow that causal ordering
    
    # predict birth
    ex.sample.3$birth[ex.sample.3$year==t] <- binomial.predict(fit.bsf.birth, ex.sample.3[ex.sample.3$year==t,])
    
    # update totalbirth
    # totalbirth is the totalbirth value that we had last year, plus
    # the birth that just happened (1) or not (0)
    ex.sample.3$totalbirth[ex.sample.3$year==t] <- ex.sample.3$l.totalbirth[ex.sample.3$year==t] +
      ex.sample.3$birth[ex.sample.3$year==t]
    
    # determine job
    # job is now no longer predicted using the covariates of the previous year
    # because that would allow our intervention on partnership to affect job
    # we want to eliminate the pathway via job to determine the direct effect of our intervention
    # there are various ways to 'control' job to get that direct effect (see Wang & Arah)
    # here I draw from the natural course distribution stochastically
    # specifically, I draw from the natural course distribution at the same time t

    # take probabilities from natural course
    pt1 <- mean(ex.sample$job.full[ex.sample$year==t],na.rm=T)
    pt2 <- mean(ex.sample$job.edu[ex.sample$year==t],na.rm=T)
    pt3 <- mean(ex.sample$job.other[ex.sample$year==t],na.rm=T)
    
    # take sample size from current dataset
    n <- length(ex.sample.3$job.full[ex.sample.3$year==t])
    
    # and then draw stochastically from multinomial distribution
    x <- rmultinom(n,1,prob=c(pt1,pt2,pt3))
    ex.sample.3$job.full[ex.sample.3$year==t] <- x[1,]
    ex.sample.3$job.edu[ex.sample.3$year==t] <- x[2,]
    ex.sample.3$job.other[ex.sample.3$year==t] <- x[3,]
    rm(x)
    # note, in the above I draw from the rows instead of the columns
    # due to how the rmultinom function works
    
    # update totaledu
    # if someone is in education this year, then we add one
    # to the total years spent in education
    ex.sample.3$totaledu[ex.sample.3$year==t] <- ex.sample.3$l.totaledu[ex.sample.3$year==t] +
      ex.sample.3$job.edu[ex.sample.3$year==t]
    
    ## determine partnership
    # take probabilities from the intervention scenario
    pt1 <- mean(ex.sample.2$partner.single[ex.sample.2$year==t])
    pt2 <- mean(ex.sample.2$partner.cohab[ex.sample.2$year==t])
    pt3 <- 0 # we know that marriage is not allowed in the intervention
    
    # take sample size from current dataset
    n <- length(ex.sample.3$partner.single[ex.sample.3$year==t])
    
    # and then draw stochastically from multinomial distribution
    x <- rmultinom(n,1,prob=c(pt1,pt2,pt3))
    ex.sample.3$partner.single[ex.sample.3$year==t] <- x[1,]
    ex.sample.3$partner.cohab[ex.sample.3$year==t] <- x[2,]
    ex.sample.3$partner.married[ex.sample.3$year==t] <- x[3,]
    rm(x)
    
    # predict censoring
    # (has to be done last, since it depends on the values
    # of the (updated) covariates at time t, and not the lagged ones)
    ex.sample.3$censor[ex.sample.3$year==t] <- binomial.predict(fit.bsf.censor, ex.sample.3[ex.sample.3$year==t,])
    
  }
  
  # after each iteration, we save the birth information from
  # both the natural course (ex.sample)
  # and from its counterfactual equivalent under the intervention (ex.sample.2)
  # total number of children in birth groups at various ages/years
  # and mean number of children in birth groups at various ages/years
  for(i in 1986:2008-1985) {
    output.dat[bs,1,i] <- sum(ex.sample$birth[ex.sample$year<=(i+1985)],na.rm=T)
    output.dat[bs,2,i] <- sum(ex.sample.2$birth[ex.sample.2$year<=(i+1985)],na.rm=T)
    output.dat[bs,3,i] <- sum(ex.sample.3$birth[ex.sample.3$year<=(i+1985)],na.rm=T)
    
    output.dat[bs,4,i] <- mean(ex.sample$totalbirth[ex.sample$year==(i+1985)],na.rm=T)
    output.dat[bs,5,i] <- mean(ex.sample.2$totalbirth[ex.sample.2$year==(i+1985)],na.rm=T)
    output.dat[bs,6,i] <- mean(ex.sample.3$totalbirth[ex.sample.3$year==(i+1985)],na.rm=T)
  }
  
  ## frequency table of children at various ages/years
  for(i in 1986:2008-1985) {
    # natural course
    temp.table <- table(ex.sample$totalbirth[ex.sample$year==(i+1985)])
    
    if(length(temp.table) <= 15) {
      
      extra <- 15-length(temp.table)
      
      temp.table <- c(temp.table,rep(0,extra))  
    }
    babytable[bs,,i] <- temp.table
    
    # intervention
    temp.table.2 <- table(ex.sample.2$totalbirth[ex.sample.2$year==(i+1985)])
    
    if(length(temp.table.2) <= 15) {
      
      extra <- 15-length(temp.table.2)
      
      temp.table.2 <- c(temp.table.2,rep(0,extra))  
    }
    
    if(length(temp.table.2) > 15) {
      
      temp.table.2 <- temp.table.2[1:15]  
    }
    babytable.2[bs,,i] <- temp.table.2
    
    # intervention direct effect
    temp.table.3 <- table(ex.sample.3$totalbirth[ex.sample.3$year==(i+1985)])
    
    if(length(temp.table.3) <= 15) {
      
      extra <- 15-length(temp.table.3)
      
      temp.table.3 <- c(temp.table.3,rep(0,extra))  
    }
    
    if(length(temp.table.3) > 15) {
      
      temp.table.3 <- temp.table.3[1:15]  
    }
    babytable.3[bs,,i] <- temp.table.3
    
    
  }
  
  # as indicated above, we could here also save information
  # on the mediators, especially if we are doing a mediation analysis
  
  print(bs)
}

# now we can compare
# the natural course with intervention
# and the empirical data with the natural course

babies.nc <- apply(output.dat[,1,],2,mean)
babies.int <- apply(output.dat[,2,],2,mean)
babies.int.de <- apply(output.dat[,3,],2,mean)

tfr.nc <- apply(output.dat[,4,],2,mean)
tfr.int <- apply(output.dat[,5,],2,mean)
tfr.int.de <- apply(output.dat[,6,],2,mean)

plot(16:38,babies.nc,type='l',lwd=2,xlab='Age',ylab='Number of babies born',
     main='Increase in Education Scenario')
lines(16:38,babies.int,lwd=2,lty=5)
lines(16:38,babies.int.de,lwd=2,lty=5,col='red')
legend(legend=c('Natural course', 'Scenario'),16,8000,lwd=2,
       lty=c(1,5))
# when all those in 'other' have a full-time job
# there is a clear difference between the number of babies born
# compared with the natural course

# measure the 'effect' on babies
TE <- apply(output.dat[,2,] - output.dat[,1,],2,mean); TE
apply(output.dat[,2,],2,mean) - apply(output.dat[,1,],2,mean)
# E[Yx* - Yx] = E[Yx*] - E[Yx]

# and the controlled direct effect
CDE <- apply(output.dat[,3,] - output.dat[,1,],2,mean); CDE

# so then the indirect effect is
TE - CDE

# %wise this is
(TE - CDE) / TE # about 6-7% (at the end) of the intervention on partnership goes via job

# we can also look at relative effects:
TE.rel <- apply(output.dat[,2,] / output.dat[,1,],2,mean); TE.rel
apply(output.dat[,2,],2,mean) / apply(output.dat[,1,],2,mean)
# E[Yx* / Yx] / E[Yx*] - E[Yx]

# but to determine significance of effects
# regardless of whether they are absolute or relative, we should take covariance into account
# let's go with a 5% significance level:
apply(output.dat[,2,] - output.dat[,1,],2,quantile,probs=c(0.025,0.975))
# note that this is NOT the same as determining a CI for the natural course
# and a CI for the intervention
# and then checking of those CI's overlap.
# comparing CI's for the scenarios is not an appropriate test since it does not take covariance into account
# and is therefore generally more conservative
# When we, instead, take a CI for the difference, this is a CI for the EFFECT, and therefore
# also functions as a test of the difference








## doing validity checks

## see also validity checks in Keil et al. and in my SSM paper





# in case there is a strong difference between
# the natural course and the empirical data
# we would have to make a more detailed model
# for example, we could add specific interaction terms between variables
# or perhaps we are missing important predictors of fertility in our model

# this is an artificially constucted dataset to show the steps of the
# g-formula in r. The relations in this data are not real
# Given that model calibration is straightforward
# we will not here demonstrate it further.
