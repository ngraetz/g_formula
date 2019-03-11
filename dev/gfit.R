# function for running multiple moel fits simulataneously
gfit.init <- function(formulas, families, functions, data=NULL, kwargs=NULL){
    if(class(formulas) != "list"){
        stop("formulas argument must be a list of formulas.")
    }
    if(class(families) != "list"){
        stop("families argument must be a list of family objects.")
    }
    if(class(functions) != "list"){
        stop("functions argument must be a list of functions.")
    }
    parLengths <- sapply(list(formulas, families, functions), length)
    if(length(unique(parLengths)) != 1){
        stop("families, formulas, and functions arguments must be same length.")
    }
    if(is.null(kwargs)){
        kwargs <- lapply(1:length(formulas), function(i) NULL)
    }
    if(length(kwargs) != length(formulas)){
        stop("kwargs must be a list of lists, equal in length to formulas.")
    }
    lapply(1:length(formulas), function(i){
        do.call(
            functions[[i]],
            c(
                list(formula=formulas[[i]], family=families[[i]], data=data),
                kwargs[[i]]
            )
        )
    })
}

# Helper function for multinomial and logistic
simPredict <- function(DF, model_list, model_index){
    model_ <- model_list[[model_index]]
    if(class(model_)[1] == "multinom"){
        probs <- predict(model_, type="probs", newdata=DF)
        opts <- colnames(probs)
        sim <- apply(probs, 1, function(p) sample(opts, 1, prob=p))
    }
    if(class(model_)[1] == "glm"){
        if(model_$family$family == 'binomial') sim <- rbinom(nrow(DF), 1, predict(model_, newdata=DF, type="response"))
        if(model_$family$family == 'gaussian') sim <- predict(model_, newdata=DF)
    }
    sim
}

# Helper function for drawing stochastically from natural course distribution (right now just multinomial and logistic).
# Nick note: I don't totally understand why we draw stochastically from the variable distribution instead of literally using the 
# values for each individual from the desired course... doesn't this kind of fuck up the covariance within individuals... 
# I guess that doesn't really matter though because it's not like we do anything with individual simulants after the fact.
simScenario <- function(DF, model_list, course_DF, model_index){
  model_ <- model_list[[model_index]]
  target_var_name <- all.vars(model_$call$formula)[1] ## Infer target variable from DV specified in model formula.
  if(class(model_)[1] == "multinom"){
    ## Get probability distribution from provided course simulation (natural course or intervention course).
    course_probs <- course_DF %>%              ## Use course DF to get probability distribution.
      filter(year == max(DF$year)) %>%         ## Sample over course distribution in the year currently being updated.
      select_(target = target_var_name) %>%
      count(target) %>%
      mutate(freq = n / sum(n)) %>% 
      select(target, freq) %>%
      spread(target, freq)
    probs <- do.call("rbind", replicate(nrow(DF), course_probs, simplify = FALSE))
    opts <- colnames(probs)
    sim <- apply(probs, 1, function(p) sample(opts, 1, prob=p))
  }
  if(class(model_)[1] == "glm"){
    if(model_$family$family == 'binomial') {
      course_prob <- course_DF %>% 
      filter(year == max(DF$year)) %>%
      select_(target = target_var_name) %>%
      summarise(mean(target, na.rm=T))
      probs <- rep(as.numeric(course_prob), nrow(DF))
      sim <- rbinom(nrow(DF), 1, probs)
    }
    if(model_$family$family == 'gaussian') {
      course_mean <- course_DF %>% 
        filter(year == max(DF$year)) %>%
        select_(target = target_var_name) %>%
        summarise(mean(target, na.rm=T))
      course_sd <- course_DF %>% 
        filter(year == max(DF$year)) %>%
        select_(target = target_var_name) %>%
        summarise(sd(target, na.rm=T))
      sim <- rnorm(nrow(DF), mean = as.numeric(course_mean), sd = as.numeric(course_sd))
    }
  }
  sim
}

# Helper function to enforce specified intervention to a dataset after it is updated in each step of progressSimulation().
assert_intervention_rules <- function(DF, rule) {
  DF <- DF %>%
    filter(year == max(DF$year)) %>%
    select_(target = rule[1]) %>%
    mutate(gsub(rule[2], rule[3], target))
}

# run the lag updates and the deterministic and probabilistic rules
progressSimulation <- function(data, lags, rules, models, intervention_rules=NULL, natural_DF=NULL, intervention_DF=NULL){
  ## Start at t=0.
  require(dplyr)
  times <- sort(unique(data$year))
  simDF <- filter(data, year == min(year)) %>%
    mutate(id=1:n())
  ## We have to assert the intervention rules at t=0 as well as in updating below.
  if(!is.null(intervention_rules)) {
    for(r in names(intervention_rules)){
      simDF[,r] <- intervention_rules[[r]](simDF, models, natural_DF)
    }
  }
  ## If we are running scenarios (we are passing an already-simulated natural course DF to draw from), we have to apply the scenario rules at t=0.
  ## NOTE: we can't predict at t=0, and I'm currently predicting the DV of interest and censoring in all scenarios. Does it make sense to draw both
  ## of these from the intervention course distribution at t=0? Right now the only difference between natural course and intervention course at t=0
  ## is in the variable we intervene on... 
  if(!is.null(intervention_DF)) {
    for(r in names(rules)){
      if(!grepl('simPredict', deparse(rules[[r]])[2])) {
        simDF[,r] <- rules[[r]](simDF, models, natural_DF, intervention_DF)
      }
      if(grepl('simPredict', deparse(rules[[r]])[2])) { # If it is something that needs to be predicted, instead draw from natural at t=0.
        new_rule <- deparse(rules[[r]])
        new_rule[2] <- gsub('models', 'models, natural_DF', new_rule[2])
        new_rule <- paste(new_rule, collapse = '')
        new_rule <- gsub('simPredict', 'simScenario', new_rule)
        new_rule <- gsub('function ', 'function', new_rule)
        new_rule <- eval(parse(text=new_rule))
        simDF[,r] <- new_rule(simDF, models, natural_DF, intervention_DF)
      }
    }
  }  
  ## Simulate forward year by year, updating all variables according to provided deterministic/probablistic rules.
  for(y in times[2:length(times)]){
    ## Progress forward one year and index lagged variables for this step.
    upDF <- simDF %>%
      filter(year == max(year)) %>%
      mutate(year=year+1) %>%
      mutate(age=age+1) %>%
      rbind(simDF) %>%
      arrange(id, year) %>%
      group_by(id)
    for(r in names(lags)){
      upDF <- lags[[r]](upDF)
    }
    upDF <- ungroup(upDF) %>%
      filter(year == max(year))
    ## Update this time step with rules for this simulation (probabilistic from models, draw from natural course, draw from intervention course).
    for(r in names(rules)){
      upDF[,r] <- rules[[r]](upDF, models, natural_DF, intervention_DF)
    }
    ## Assert any rules associated with this simulation (used to calculate direct/indirect effects).
    if(!is.null(intervention_rules)) {
      for(r in names(intervention_rules)){
        upDF[,r] <- intervention_rules[[r]](upDF, models, natural_DF)
      }
    }
    ## Add updated time step to simulated DF.
    simDF <- rbind(simDF, upDF)
  }
  simDF
}

compile_sims <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
    bs[[i]][[name]] %>%
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
      censor = mean(censor),
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
    mutate(type=name)
  return(df)
}

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


compile_sims_new <- function(name, bootruns) {
  df <- bind_rows(lapply(1:length(bootruns), function(i){
    bootruns[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    group_by(age) %>%
    select(-matches('^l.'), -id, -year, -wave) %>%
    mutate(
      pcohab = as.numeric(m_relation == "cohab"),
      pmarried = as.numeric(m_relation == "married"),
      pnocontact = as.numeric(m_relation == "no_contact"),
      pcontact = as.numeric(m_relation == "contact"),
      pexcellenthealth = as.numeric(m_health == "Excellent"),
      ppoorhealth = as.numeric(m_health == "Poor"),
      pcollege = as.numeric(m_edu == "college"),
      psomecollege = as.numeric(m_edu == "some_college"),
      phs = as.numeric(m_edu == "hs"),
      plesshs = as.numeric(m_edu == "less_hs"),
      phighincome = as.numeric(m_poverty == "perc_300_plus"),
      ppov = as.numeric(m_poverty == "perc_50_99"),
      pextpov = as.numeric(m_poverty == "perc_0_49"),
      pmotherdepression = as.numeric(m_depression), 
      pfulltime = as.numeric(m_job_type == "full_time"),
      punemployed = as.numeric(m_job_type == "unemployed"),
      pparttime = as.numeric(m_job_type == "part_time")
    ) %>%
    select(age, c_ppvt, m_kids, m_f_in_jail, m_evicted, m_anxiety, m_depression, wodtke, censor) %>%
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
    mutate(type=name)
  return(df)
}

compile_sims_simple <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
    bs[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    group_by(age, sim) %>%
    summarize(
      mppvt = mean(c_ppvt), 
      married = mean(bin_married),
      lesshs = mean(bin_lesshs),
      pov = mean(bin_pov),
      jail = mean(m_f_in_jail),
      depression = mean(m_depression),
      unemployed = mean(bin_unemployed),
      censor = mean(censor),
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
    mutate(type=name)
  return(df)
}


compile_sims_table <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
    bs[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    group_by(age, sim) %>%
    summarize(
      mppvt = mean(c_ppvt), 
      married = mean(bin_married),
      lesshs = mean(bin_lesshs),
      pov = mean(bin_pov),
      jail = mean(m_f_in_jail),
      depression = mean(m_depression),
      unemployed = mean(bin_unemployed),
      censor = mean(censor),
      mwodtke = mean(wodtke)
      ) %>%
    mutate(name=name)
  df <- dcast(as.data.table(df), age ~ name + sim, value.var=c('mppvt','married','lesshs','pov','jail','depression','unemployed','censor','mwodtke'))
  return(df)
}
