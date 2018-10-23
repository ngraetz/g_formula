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
    else{
        sim <- rbinom(nrow(DF), 1, predict(model_, newdata=DF, type="response"))
    }
    sim
}

# Helper function for drawing stochastically from natural course distribution (right now just multinomial and logistic).
# Nick note: I don't totally understand why we draw stochastically from the variable distribution instead of literally using the 
# values for each individual from the desired course... doesn't this kind of fuck up the covariance within individuals... 
simScenario <- function(DF, course_DF, model_list, model_index){
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
  else{
    course_prob <- course_DF %>% 
      filter(year == max(DF$year)) %>%
      select_(target = target_var_name) %>%
      summarise(mean(target, na.rm=T))
    probs <- rep(as.numeric(course_prob), nrow(DF))
    sim <- rbinom(nrow(DF), 1, probs)
  }
  sim
}

# run the lag updates and the deterministic and probabilistic rules
progressSimulation <- function(data, lags, rules, models){
    require(dplyr)
    times <- sort(unique(data$year))
    simDF <- filter(data, year == min(year)) %>%
        mutate(id=1:n())
    for(y in times[2:length(times)]){
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
        for(r in names(rules)){
            upDF[,r] <- rules[[r]](upDF, models)
        }
        simDF <- rbind(simDF, upDF)
    }
    simDF
}
