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
simNatural <- function(naturalDF, model_list, model_index){
  model_ <- model_list[[model_index]]
  if(class(model_)[1] == "multinom"){
    
  }
  else{

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
