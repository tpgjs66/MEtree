

predict.MEtree <- function(tree , MCMCglmm, newdata,formula, id=NULL, EstimateRandomEffects=TRUE,...){

  treePrediction <- predict.party(tree,newdata)

  # If we aren't estimating random effects, we
  # just use the tree for prediction.
  if(!EstimateRandomEffects){
    return(treePrediction)
  }
  # Get the group identifiers if necessary
  if(is.null(id)){
    id <- newdata[,as.character((MCMCglmm$Random$formula[[2]][[3]]))]
  }
  # Error-checking: the number of observations
  # in the dataset must match the sum of NumObs
  if(length(newdata[,id]) != dim(newdata)[1]){
    stop("number of observations in newdata does not match the length of the group identifiers")
  }
  ### Use the formula to get the target name
  TargetName <- formula[[2]]
  # Remove the name of the data frame if necessary
  if(length(TargetName)>1) TargetName <-TargetName[3]

  ActualTarget <- newdata[,toString(TargetName)]

  completePrediction <- treePrediction

  # Get the identities of the groups in the data
  # This will be slow - does LME have a faster way?
  uniqueID <- unique(id)

  # Get the random effects from the estimated MCMCglmm, in case there is overlap
  estRE <- ranef(object, use = ("mean"))

  for(i in 1:length(uniqueID)){
    # Identify the new group in the data
    nextID <- uniqueID[i]
    thisGroup <- id==nextID

    # If this group was in the original estimation, apply its random effect
    filter<-grepl(toString(uniqueID[i]),rownames(estRE))
    estEffect <- estRE[filter,]

    if(is.na(estEffect)){
      # Check for non-missing target
      nonMissing <- !is.na(ActualTarget[thisGroup])
      numAvailable <- sum(nonMissing)

      # If all the targets are missing, accept the
      # tree prediction; otherwise, estimate
      if(numAvailable>0) {
        R <- object$ErrorVariance * diag(numAvailable)
        D <- object$BetweenMatrix
        Z <- matrix(data=1,ncol=1, nrow=numAvailable)
        W <- solve(R + Z %*% D %*% t(Z))
        effect <- D %*% t(Z) %*% W %*%
          subset(ActualTarget[thisGroup] - treePrediction[thisGroup],subset=nonMissing)
        completePrediction[thisGroup] <- treePrediction[thisGroup]+effect
      }
    } else {
      completePrediction[thisGroup] <- treePrediction[thisGroup]+estEffect
    }

  }

  return(completePrediction)

}
