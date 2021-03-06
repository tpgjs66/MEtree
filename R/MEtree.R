

################################################################################
###                          Define MEtree function                          ###
################################################################################

MEtree5<-function(data,formula,random) {
  ErrorTolerance=10
  MaxIterations=200
  #parse formula
  Predictors<-paste(attr(terms(formula),"term.labels"),collapse="+")
  TargetName<-formula[[2]]
  Target<-data[,toString(TargetName)]
  #set up variables for loop
  ContinueCondition<-TRUE
  iterations<-0
  #set up the initial target
  OriginalTarget<-(Target)
  TargetLevels<-levels(OriginalTarget)
  oldDIC<- Inf

  # Make a new data frame to include all the new variables
  newdata <- data
  newdata[,TargetLevels] <- 0

  m.list<-list()
  tree.list<-list()
  while(ContinueCondition){

    # Count iterations
    iterations <- iterations+1
    print(paste("############### Main Iteration ",iterations,"###############"))

    # Target response will be updated from the previous result.
    if (iterations<2){
      newdata[,"OriginalTarget"] <- as.factor(OriginalTarget)
    }else {
      newdata[,"OriginalTarget"] <- as.factor(MCMCTarget)
    }

    # Build CHAID tree
    ctrl <- chaid_control(alpha2=0.05,alpha4=0.05,
                          minsplit = 2*floor(nrow(data)/100),
                          minbucket=floor(nrow(data)/100), minprob=1)
    tree <- chaid(formula(paste(c("OriginalTarget", Predictors),collapse = "~"))
                  ,data = newdata, control = ctrl)

    tree.list[[iterations]]<-tree

    # Get terminal node
    newdata[,"nodeInd"] <- 0
    newdata["nodeInd"] <-as.factor(predict(tree,newdata=newdata,type="node"))

    # Get variables (alternative-specific) that identify the node for
    # each observation
    predCHAID <- as.data.frame(predict(tree,newdata=newdata,type="prob"))

    for (i in TargetLevels) {
      newdata[,TargetLevels] <- predCHAID[,TargetLevels]
    }

    CHAIDTarget<-c()
    # Update adjusted target based on CHAID predicted probs.
    repeat{
      for(k in 1:length(OriginalTarget)){
        t<-TargetLevels
        # Draw a decision based on probs
        CHAIDTarget[k]<-sample(t,1,replace=FALSE,
                               prob=newdata[k,TargetLevels])
      }
      if ((length(table(OriginalTarget))==5)){break}
    }

    newdata[,"CHAIDTarget"] <- as.factor(CHAIDTarget)

    # Fit MCMCglmm
    k <- length(TargetLevels)
    I <- diag(k-1)
    J <- matrix(rep(1, (k-1)^2), c(k-1, k-1))

    prior <- list(
      G = list(G1 = list(V = diag(k-1), n = k-1)),
      R = list(fix=1,V= (1/k) * (I + J), n = k-1))

    m <- MCMCglmm(fixed = OriginalTarget ~ -1 + trait
                    +trait:(nodeInd+CHAIDTarget),

                  random = ~ idh(trait):HHID,# ~ idh(trait-1+nodeInd):ind.id ??
                  rcov = ~idh(trait):units,

                  prior = prior, # Add fix=1 if you want fix R-structure
                  burnin =1000,
                  nitt = 11000,
                  thin = 10,
                  # This option saves the posterior distribution of
                  # random effects in the Solution mcmc object:
                  pr = TRUE,
                  #pl = TRUE,
                  family = "categorical",
                  #saveX = TRUE,
                  #saveZ = TRUE,
                  #saveXL = TRUE,
                  data = newdata,
                  verbose = T,
                  #slice = T
                  #singular.ok = T
    )

    m.list[[iterations]]<-m

    #p <- predict(m,type="terms",interval="prediction")[,1]
    p <- (predict(m,type="terms",interval="none",posterior="all"))
    #p <- (predict(m,type="terms",interval="none",posterior="distribution"))
    #p <- (predict(m,type="terms",interval="none",posterior="mean"))
    #p <- (predict(m,type="terms",interval="none",posterior="mode"))

    # Predicted probability with marginalizing the random effect
    #p <- predict(m,type="terms", interval="none",posterior="mean",marginal=NULL)
    #p <- predict(m,type="terms", interval="none",
    #            posterior="mean",marginal=m$Random$formula)


    pred<-c()

    pred$b<-p[1:nrow(newdata)]
    pred$c<-p[(nrow(newdata)+1):(2*nrow(newdata))]
    pred$d<-p[(2*nrow(newdata)+1):(3*nrow(newdata))]
    pred$e<-p[(3*nrow(newdata)+1):(4*nrow(newdata))]

    pred<-as.data.frame(pred)
    pred<-cbind(a = 0, pred)

    colnames(pred)<-TargetLevels

    predMCMC<-c()
    predMCMC$a<-exp(pred[,1])/(1+exp(pred[,2])+exp(pred[,3])+exp(pred[,4])+exp(pred[,5]))
    predMCMC$b<-exp(pred[,2])/(1+exp(pred[,2])+exp(pred[,3])+exp(pred[,4])+exp(pred[,5]))
    predMCMC$c<-exp(pred[,3])/(1+exp(pred[,2])+exp(pred[,3])+exp(pred[,4])+exp(pred[,5]))
    predMCMC$d<-exp(pred[,4])/(1+exp(pred[,2])+exp(pred[,3])+exp(pred[,4])+exp(pred[,5]))
    predMCMC$e<-exp(pred[,5])/(1+exp(pred[,2])+exp(pred[,3])+exp(pred[,4])+exp(pred[,5]))

    predMCMC<-as.data.frame(predMCMC)
    colnames(predMCMC)<-TargetLevels

    # Get the DIC to check on convergence
    if(!(is.null(m))){
      newDIC <- m$DIC
      ContinueCondition <- (abs(oldDIC-newDIC)>ErrorTolerance &
                              iterations < MaxIterations)
      oldDIC <- newDIC
      print(paste("###### DIC : ", m$DIC, " ######"))

      # Update prob.
      for (i in TargetLevels){
        newdata[,TargetLevels] <- predMCMC[,TargetLevels]
      }

      # # Update adjusted target based on logit prob.
      # for(k in 1:length(AdjustedTarget)){
      #   AdjustedTarget[k]<-sum(cumsum(mlogitfit$probabilities[k,])<runif(1))+1
      #
      # }
      # newdata[,"AdjustedTarget"] <- AdjustedTarget


      # Update adjusted target based on MCMCglmm predicted probs.
      MCMCTarget<-c()

      repeat{
        for(k in 1:length(OriginalTarget)){
          t<-TargetLevels
          MCMCTarget[k]<-sample(t,1,replace=FALSE,
                                prob=newdata[k,TargetLevels])
        }
        if ((length(table(MCMCTarget))==5)){break}
      }
      newdata[,"MCMCTarget"] <- as.factor(MCMCTarget)

    }
    else{ ContinueCondition<-FALSE }
  }


  #return final model fits and convergence info.
  return(list(
    CHAID.tree=tree.list,
    MCMCglmm.fit=m.list,
    Conv.info=newDIC-oldDIC,
    n.iter=iterations
  ))
}

