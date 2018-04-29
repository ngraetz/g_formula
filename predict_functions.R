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
