#'@title  Monte Carlo simulation for one or two RR variables
#' 
#' @description Simulate and analyse bivariate data including either one or two RR variables. Useful for power analysis, parametric bootstraps or for testing the effects of noncompliance on the stability of estimates.
#' @param numRep number of replications
#' @param n sample size
#' @param pi true proportion of carriers of sensitive attribute (for 2 RR variables: \code{vector})
#' @param model either one or two RR model (as \code{vector}), see \code{\link{RRuni}}
#' @param p randomization probability (for 2 RR variables: a \code{list})
#' @param rho true correlation in population (also known as point-biserial/point-tetrachoric correlation, in case of one/two dichotomous RR variables)
#' @param complyRates vector with two values giving the proportions of participants who adhere to the instructions in the subset with or without the sensitive attribute, respectively (for 2 RR variables: a \code{list})
#' @param sysBias probability of responding 'yes' (coded as 1 in the RR variable) in case of non-compliance for carriers and non-carriers, respectively. See \code{\link{RRgen}}
#' @param method vector specifying which RR methods to be used in each replication. For a single RR variable, all three methods can be used. For 2 RR variables, only \code{\link{RRcor}} is available.
#' @param alpha significance threshold for testing the logistic regression parameter \code{beta}
#' @param groupRatio only for multiple group models: ratio of groups (for 2 RR variables: \code{vector})
#' @param MLest correct estimates of \code{RRuni} if pi is outside of [0,1]
#' @param nCPU integer: how many processors to use? (use 'max' for automatic detection on Windows)
#' @return an object \code{RRsimu} which contains the estimated parameters \code{parEsts} and a matrix \code{results} with mean parameters and standard errors across replications
#' @details In case of using only one RR variable, the second, directly measured variable is sampled from a normal distribution with shifted means, depending on the true state on the sensitive attribute (i.e., the true, underlying values on the RR variable). For dichotomous RR variables, this corresponds to the assumption of an ordinary t-test, where the dependent variable is normally distributed within groups with equal variance. The difference in means is chosen in a way, to obtain the point-biserial correlation defined by \code{rho}.
#' 
#' In case of two dichotomous RR variables, the true group membership of individuals is sampled from a 2x2 cross table. Within this table, probabilities are chosen in a way, to obtain the point-tetrachoric correlation defined by \code{rho}
#' 
#' Note, that for the FR model with multiple response categories (e.g., from 0 to 4), the specified \code{rho} is not the exact target of the sampling procedure. It assumes a normal distribution for each true state, with constant differences between the groups (i.e., it assumes an interval scaled variable).
#' @examples # Simulate data according to the Warner model
#' mcsim <-  RRsimu(numRep=5, n=200, pi=.3, model="Warner", p=2/12, rho=.6)
#' mcsim
#' plot(mcsim)
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
RRsimu <- function(numRep, n, pi, model, p, rho, complyRates=c(1,1), 
                   sysBias=c(0,0), method=c("RRuni","RRcor","RRlog"), 
                   alpha=0.05, groupRatio=0.5, MLest=TRUE, nCPU=1){
  try(stopCluster(cl.tmp), silent = T)
  
  modelNames <- c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","direct","CDM","CDMsym","SLD")
  model <- pmatch(model, modelNames, duplicates.ok=T )
  model <- modelNames[model]
  
  nn <- c("pi.true","pi.RRuni","pi2.true","pi2.RRuni","piSE.RRuni",
          "cor.true","cor.RRcor","beta.true","betaSE.true","beta.RRlog","betaSE.RRlog",
          "beta.deltaG2.RRlog","par2.true","par2.RRuni","par2SE.RRuni")
  #   allvars <- c(pi[1], pi[1], pi[2], pi[2] , NA, rho, rho, )
  nvar <- length(nn)
  
  # check: one or two RR variables
  model <- model[model!="direct"]
  if (length(model)<1 | length(model) > 2)
    stop("Definition of either one or two RR variables required (e.g., model='Warner').")
  twoRR <- length(model)==2
  if (!twoRR){
    # if (class(complyRates)!= "list") complyRates <- list(complyRates)
    # if only one RR variable included: transform correlation to mean difference
    diffMean <- sqrt(rho^2/(pi[1]*(1-pi[1])*(1-rho^2))) * sign(rho)
    if (model == "FR" & length(pi)>2){
      ppi <- pi
      if (min(pi)==0){
        ppi <- pi + .001
      }
      diffMean <-  sign(rho)*mean(sqrt(rho^2/(ppi*(1-ppi)*(1-rho^2)))) * sqrt(pi%*% pi)
    }
  }else{
    # two RR variables: think about data generation !
    if (missing(complyRates)) complyRates <- list(c(1,1),c(1,1))
    if (missing(groupRatio)) groupRatio <- c(.5,.5)
    if(length(pi)!=2 ||min(pi)<=0 ||max(pi)>=1) stop("Please provide a vector pi with two values in (0,1)")
  }
  
  ### CDM and CDMsym not available in RRcor
  if (any(model %in% c("CDM","CDMsym"))){
    method <- method[!(method == "RRcor")]
    warning("Models 'CDM' and 'CDMsym' are not available in RRcor because both define a third group of cheaters (see vignette('RRreg') ).")
  }
  
  ############################################### ONE RR   ###############################################
  getParEsts.oneRR <- function(replic){
    parEsts <- matrix(NA,nrow=replic,ncol=nvar)
    colnames(parEsts) <- nn
    parEsts <- as.data.frame(parEsts)
    for (i in 1:(replic)){
      # generate continuous and discrete variable (point-biserial correlation)
      dat <- RRgen(n=n,pi.true=pi,model=model,p=p,complyRates=complyRates,
                   sysBias=sysBias, groupRatio=groupRatio)
      n.true <- table(dat$true)
      #       nn <- names(n.true)
      #       if (length(nn) != length(p)){
      #         
      #       }
      dat$cov <- rep(0,n)
      
      # generate normal data dependent on true value on RR variable
      # k : loop through true states
      for (k in 0:(length(n.true)-1) ){
        #         for(l in 0:(length(pi)-1)){
        #           nsel <- names(n.true)[k+1]==paste(l)
        #           if(nsel){
        #             print(k*diffMean)
        dat[dat$true==k,]$cov <-  rnorm(n.true[k+1],k*diffMean,1)
        #           }
        #         }
      }
      if (is.null(dat$group)) 
        group <- rep(1,nrow(dat))
      else 
        group <- dat$group
      
      # analyse with RRcor and RRlog
      if ("RRcor" %in% method){
        cor1 <- RRcor(x=dat[,c("response","cov")],
                      models=c(model,"d"),p.list= list(p,1),
                      group=group )  
        parEsts[i,"cor.true"] <- cor(dat$true,dat$cov)
        parEsts[i,"cor.RRcor"] <- cor1$r["cov","response"]
      }
      
      uni1 <- RRuni(response=dat$response,model=model,p=p,group=group,MLest=MLest)  
      if ("RRuni" %in% method){
        parEsts[i,"pi.true"] <- table(dat$true)["1"]/n
        parEsts[i,"pi.RRuni"] <- uni1$pi[ifelse(model=="FR",2,1)]
        parEsts[i,"piSE.RRuni"] <- uni1$piSE[ifelse(model=="FR",2,1)]
        
        if (model %in% c("CDMsym","CDM")){
          parEsts[i,"pi.true"] <- table(dat[dat$comply==1,"true"])["1"]/n
          parEsts[i,"par2.true"] = 1- colMeans(dat)["comply"]
          parEsts[i,"par2.RRuni"] = uni1$gamma
          #         parEsts[i,"par2SE.RRuni"] = uni1$gammaSE
        }else if (model== "SLD"){
          parEsts[i,"par2.true"] = colMeans(dat[dat$true==1,])["comply"]
          parEsts[i,"par2.RRuni"] = uni1$t
          #         parEsts[i,"par2SE.RRuni"] = uni1$tSE
        } 
      }
      
      if ("RRlog" %in% method){
        parEsts[i,"beta.true"] <- NA
        # adjust dependent variable for CDM: having sensitive attribute and commiting in RR design!
        if (model %in% c("CDMsym","CDM")){
          dat$true <- dat$true & dat$comply
        }
        tryCatch({
          glm1 <- glm(true~cov,dat,family=binomial(link = "logit"))
          parEsts[i,"beta.true"] <- glm1$coefficients["cov"]
          parEsts[i,"betaSE.true"] <- summary(glm1)$coef[2,2]    
        },
        error=function(e) {})
        
        log1 <- RRlog(response~cov,data=dat,model=model,p=p,group=group, LR.test=T)
        parEsts[i,"beta.RRlog"] <- log1$coefficients["cov"]
        parEsts[i,"betaSE.RRlog"] <- sqrt(log1$vcov["cov","cov"])
        parEsts[i,"beta.deltaG2.RRlog"] <- -2*log1$deltaLogLik["cov"]
      } 
      
    }
    parEsts
  }
  
  ############################################### TWO RR   ###############################################
  getParEsts.twoRR <- function(replic){
    parEsts <- matrix(NA,nrow=replic,ncol=nvar)
    colnames(parEsts) <- nn
    parEsts <- as.data.frame(parEsts)
    for (i in 1:(replic)){
      # generate data for two discrete, dichotomous RR variables
      # for a fourfold table, with marginal proportions pi[1] and pi[2]
      # a b
      # c d
      cov <- rho * sqrt(pi[1]*(1-pi[1])*pi[2]*(1-pi[2]))
      a <- pi[1]*pi[2] + cov
      b <- pi[1]*(1-pi[2]) - cov
      c <- (1-pi[1])*pi[2] - cov
      d <- (1-pi[1])*(1-pi[2]) + cov
      if(any(c(a,b,c,d)<0))
        stop("The given correlation coefficient rho=",rho," (phi coefficient) is not observable with n=",n, " and pi=c(",pi[1],", ",pi[2],").")
      dat <- expand.grid(0:1,0:1)[sample(1:4, size=n, replace=TRUE, prob=c(a,b,c,d)),]
      colnames(dat) <- c("true1", "true2")
      RR1 <- RRgen(model =model[1], pi.true=.3, p=p[[1]], complyRates = complyRates[[1]],
                   sysBias = sysBias, groupRatio = groupRatio[1], trueState = dat$true1)
      colnames(RR1) <- paste0(colnames(RR1),1)
      RR2 <- RRgen(model =model[2], pi.true=.3, p=p[[2]], complyRates = complyRates[[2]],
                   sysBias = sysBias, groupRatio = groupRatio[2], trueState = dat$true2)
      colnames(RR2) <- paste0(colnames(RR2),2)
      dat <- cbind(RR1, RR2)
      group <- dat[,grep("^group",colnames(dat))]
      
      # analyse with RRcor and RRlog
      if ("RRcor" %in% method){
        cor1 <- RRcor(x=dat[,c("response1","response2")],
                      models=model,p.list= p,
                      group=group)
        parEsts[i,"cor.true"] <- cor(dat$true1,dat$true2)
        parEsts[i,"cor.RRcor"] <- cor1$r["response1","response2"]
      }
    }
    parEsts
  }
  
  ###############################################
  
  # multi core processing
  if (nCPU!=1){
#     require(doParallel, quietly=T)
    if (nCPU=="max"){
      try(nCPU <-  as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
      if (nCPU=="max") nCPU <- 2    
    }
    cl.tmp =  makeCluster(nCPU) 
    registerDoParallel(cl.tmp, cores=nCPU)
    #   try(rm(parEsts), silent=TRUE)
    if (!twoRR){
      parEsts <- foreach(k=1:nCPU , .combine= rbind,.packages='RRreg') %dopar% { 
        getParEsts.oneRR(ceiling(numRep/nCPU)) }
    }else{
      parEsts <- foreach(k=1:nCPU , .combine= rbind,.packages='RRreg') %dopar% { 
        getParEsts.twoRR(ceiling(numRep/nCPU)) }
    }
    stopCluster(cl.tmp)
  }else{
    if(!twoRR)
      parEsts <- getParEsts.oneRR(numRep)
    else
      parEsts <- getParEsts.twoRR(numRep)
  }
  
  # significance testing of beta coefficients
  #parEsts[,"beta.perctSignif.RRlog"] <- parEsts[,"beta.prob.RRlog"]  <alpha
  parEsts <- as.matrix(parEsts)
  
  # summarize simulation results
  NAs <-  colSums(is.na(parEsts))
  results <- data.frame(Mean=colMeans(parEsts,na.rm =T),
                        #                         SE=bs.se(parEsts,)
                        SE=apply(parEsts,2,sd,na.rm =T)
  )
  results <- cbind(results,NAs)
  results <- results[!is.na(results[,1]),]
  parEsts <- parEsts[,NAs!=numRep]
  
  
  sim <- list(parEsts = parEsts, results = results,
              model=model,p=p,n=n,numRep=numRep,NAs=NAs,
              complyRates=complyRates, sysBias=sysBias, groupRatio=groupRatio)
  class(sim) <- "RRsimu"
  sim
}


# Print results of RRsimu

#' @aliases RRsimu
#' @method print RRsimu
#' @export
print.RRsimu <- function(x,...){
  if (length(x$model) ==1)
    cat( paste0(x$model," ; n= ",x$n,"; randomization probability: ", gsub(", ",", ",toString(x$p)),"\n" ))
  if (length(x$model) ==2){
    cat( paste0("\n",x$model[1]," ; n= ",x$n,"; randomization probability: ",gsub(", ",", ",toString(x$p[[1]]))))
    cat( paste0("\n",x$model[2]," ; n= ",x$n,"; randomization probability: ",gsub(", ",", ",toString(x$p[[2]])),"\n" ))
  }  
  cat(paste0("\n(",x$numRep," replications)\n"))
  print( round(x$results,5))
}

# Plot simulation results
# Generic method to plot the results of the monte carlo simulation

#' @aliases RRsimu
#' @method plot RRsimu
#' @export
plot.RRsimu <- function(x,y,...){
  parEsts <- x$parEsts
  results <- x$results
  nvar <- ncol(parEsts)
  nn <- colnames(parEsts)
  
  # how many plots?
  root <- sqrt(nvar)
  mfrow <- c( floor(root),1+floor(nvar/floor(root))) 
  plot.new()
  par(mfrow=mfrow)
  
  for (i in 1:nvar){
    hist(x= parEsts[,i],freq=FALSE,nclass=30,
         main=colnames(parEsts)[i],col="grey",xlab=paste0("mean=",round(results[i,1],4),",sd=",round(results[i,2],4)))
    abline(v=results[i,1],col="red")
    if (nn[i] %in% c("pi.RRuni","pi.RRlog")){
      dd <- function(x) dnorm(x,mean=results[i,1],sd=results["piSE.RRuni",1])
      curve(dd,col="blue",n=500,add=T)
    }
    else if (nn[i] %in% c("par2.RRuni")){
      dd <- function(x) dnorm(x,mean=results[i,1],sd=results["par2SE.RRuni",1])
      curve(dd,col="blue",n=500,add=T)
    }
    else if(nn[i]=="beta.RRlog"){
      dd <- function(x) dnorm(x,mean=results[i,1],sd=results["betaSE.RRlog",1])
      curve(dd,col="blue",n=500,add=T)
    }
    # asymptotic distribution of beta in RRlog. should be chisq for diffMean=0
    else if (nn[i]=="beta.deltaG2.RRlog"){
      wr <- function(x) dchisq(x,1)
      curve(wr, col = "blue", add = TRUE,n=500)
    }
    #     else if (nn[i]=="beta.prob.RRlog"){
    #       curve(dunif, col = "blue", add = TRUE,n=500)
    #     }
  }
  par(mfrow=c(1,1))
}