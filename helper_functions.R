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
