#' Bivariate correlations including randomized response variables
#' 
#' \code{RRcor} calculates bivariate Pearson correlations of variables measured with or without RR.
#' @param x a numeric vector, matrix or data frame.
#' @param y \code{NULL} (default) or a vector, matrix or data frame with compatible dimensions to \code{x}. 
#' @param models a vector defining which RR design is used for each variable. Must be in the same order as variables appear in \code{x} and \code{y} (by columns). Available models: \code{Warner}, \code{Kuk}, \code{FR}, \code{Mangat}, \code{UQTknown}, \code{UQTunknown}, \code{Crosswise}, \code{SLD} and \code{direct} (i.e., no randomized response design)
#' @param p.list a \code{list} containing the randomization probabilities of the RR models defined in \code{models}. Either, all \code{direct}-variables (i.e., no randomized response) in \code{models} can be excluded in \code{p.list}; or, if specified, randomization probabilities \code{p} are ignored for \code{direct}-variables.
#' @param group a matrix defining the group membership of each participant (values 1 and 2) for all multiple group models(\code{SLD}, \code{UQTunknown}). If only one of these models is included in \code{models}, a vector can be used. For more than one model, each column should contain one grouping variable 
#' @param bs.n number of samples used to get bootstrapped standard errors
#' @param bs.type Choose between parametric (\code{p}) and/or nonparametric (\code{n}) bootstrap to obtain standard errors (note: \code{bs.n} has to be larger than 0). The parametric bootstrap is based on the assumption, that the continuous variable is normally distributed within each group defined by the true state of the RR variable. For polytomous forced response (FR) designs, the RR variable is assumed to have equally spaced distances between categories (i.e., that it is interval scaled)
#' @param nCPU number of CPUs used for the bootstrap
#' @return an object \code{RRcor}, can be analyzed by means of \code{\link{summary}}. Use \code{object$r} to get the correlation matrix only.
#' @details Correlations of RR variables are calculated by the method of Fox & Tracy (1984) by interpreting the variance induced by the RR procedure as uncorrelated measurement error. Since the error is independent, the correlation can be corrected to obtain an unbiased estimator.
#' @seealso \code{vignette('RRreg')} or \url{https://dl.dropboxusercontent.com/u/21456540/RRreg/index.html} for a detailed description of the RR models and the appropriate definition of \code{p} 
#' @references Fox, J. A., & Tracy, P. E. (1984). Measuring associations with randomized response. \emph{Social Science Research, 13}, 188-197.
#' @examples 
#' # generate first RR variable
#' n <-1000
#' p1 <- c(.4,.6)
#' gData <- RRgen(n,pi=.3,model="Kuk",p1)
#' n1 <- table(gData$true)[1]
#' gData[gData$true==0,4] <-  rnorm(n1)
#' gData[gData$true==1,4] <-  rnorm(n-n1,mean=0.7)
#' 
#' # generate second RR variable
#' p2 <- c(.7,.4)
#' temp <- RRgen(n1,pi=.2,model="UQTknown",p2)
#' gData$uqtResponse[gData$true==0] <- temp$response
#' gData$uqtTrue[gData$true==0] <- temp$true
#' temp <- RRgen(n-n1,pi=.7,model="Kuk",p2)
#' gData$uqtResponse[gData$true==1] <- temp$response
#' gData$uqtTrue[gData$true==1] <- temp$true
#' 
#' # show true correlation
#' cor(gData[,c("true","V4","uqtTrue")])
#' # estimated correlation
#' RRcor(x=gData[,c("response","V4","uqtTrue")],models=c("Kuk","d","UQTknown"),p.list= list(p1,p2) )
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
RRcor <- function(x, y=NULL, models, p.list, group=NULL, bs.n=0, bs.type=c("n","p"), nCPU=1){
  
  # input handling  (abkopiert von function 'cor' ; na.method funktioniert nicht
  #   na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
  #                              "everything", "na.or.complete"))   ## returns integer!
  #   if (is.na(na.method)) 
  #     stop("invalid 'use' argument")
  
  if (is.data.frame(y) || is.vector(y) || is.matrix(y) || is.numeric(y)){   
    y.name <- deparse(substitute(y))
    y <- as.matrix(y)
    my <- ncol(y)
    ynames <- colnames(y)
    if  ( is.null(ynames) & my==1){
      ynames <- y.name
    }else if( is.null(ynames)){
      ynames <- 1:my
      ynames <- paste( "y", ynames,sep="")
    }
  }
  if (is.data.frame(x)|| is.vector(x) || is.matrix(x) || is.numeric(y)) {
    x.name <- deparse(substitute(x))
    x <- as.matrix(x)
    mx <- ncol(x)
    xnames <- colnames(x)
    if  ( is.null(xnames) & mx==1){
      xnames <- x.name
    }else if( is.null(xnames)){
      xnames <- 1:mx
      xnames <- paste( "x", xnames,sep="")
    }
  }
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }
  
  # erzeuge gemeinsamen datensatz X
  if (!is.null(x) && !is.null(y)){
    if ( !is.null(mx) && !is.null(my) && mx+my <2){
      stop("Not enough variables specified")
    }
    X <- cbind(x,y)
    colnames(X) <- c(xnames,ynames)
    m <- mx+my
  }
  
  if (!is.null(x) && is.null(y)){
    X <- x
    colnames(X) <- xnames
    m <- mx
  }
  if (is.null(x) && !is.null(y)){
    X <- y
    colnames(X) <- ynames
    m <- my
  }
  modelNames <- c("Warner","UQTknown","UQTunknown","Mangat","Kuk","FR","Crosswise","direct","SLD")
  models <- pmatch(models, modelNames, duplicates.ok=T )
  models <- modelNames[models]
  n <- nrow(X)
  
  # adjust length and entries of p.list
  if (length(p.list)<length(models) && length(p.list) == sum(models!="direct") ){
    cnt <- 1
    p.list2 <- list()
    for (i in 1:length(models)){
      if(models[i]=="direct"){
        p.list2[[i]]  <- 1
      } else {
        p.list2[[i]] <-  p.list[[cnt]]
      }
      cnt <- ifelse(models[i]=="direct",cnt,cnt+1)
    }
    p.list <- p.list2
  }
  if (!missing(group) && !is.null(group)) group <- as.matrix(group)
  group.2g <- group
  RRcheck.cor(X,m,models,p.list,colnames(X),group)
  
  # generate matrix with groups and calculate additional parameters like gamma, t, piUQ
  group.mat <- matrix( 1, nrow=n,ncol=m)
  par2 <- rep(NA,m)   # additional parameters except pi
  cnt <- 1
  for (i in 1: m){
    if (models[i] %in% c("CDM","CDMsym","UQTunknown","SLD") ){
      group.mat[,i] <-  as.numeric(  group[,cnt])
      est <- RRuni(X[,i],model=models[i],p=p.list[[i]],group=group.mat[,i],MLest=F)
      cnt <- cnt+1
      switch (models[i],
              #               "CDM" = {par2[i] <- est$gamma},
              #               "CDMsym" = {par2[i] <- est$gamma},
              "UQTunknown" ={par2[i] <- est$piUQ},
              "SLD" = {par2[i] <- est$t})
    }        
  }
  #   if (ncol(group) != n) group <- group.mat
  group <- group.mat
  
  
  # x.est : estimated true values for respondents (dependent on group)
  x.est <- X
  # for SLD and UQTunknown: which group to select (always with larger p[i] !)
  group.sel <- rep(1,m)
  
  # used sample sizes for all pairwise comparisons
  n.mat <- diag(m) * n
#   g.list <- list()
#   quotient.list <- list()
#   cnt <- 0
  
  # quotient contains the scaling factors for the observed correlation
  # rows: different models; second column for 2-group models
  quotient <- rep(0,m)
  for (i in 1:m){
    # pairwise comparison to get adequate subgroups!
#     if (i != m){
#       for (j in (i+1):m){
#         # find all available subgroups
#         subgroups <- group[,c(i,j)][!duplicated(group[,c(i,j)]),,drop=F]
#         numsub <- nrow(subgroups)
#         # calculate information index per subgroup:
#         inf.idx <- rep(0, numsub)
#         for (ss in 1:numsub){
#           # get the subgroup sample size n[i,j]
#           subgroup.n <- sum(apply(group[,c(i,j)]  == subgroups[ss,], 1, any))
#           pi <- p.list[[i]][subgroups[ss,1]]
#           pj <- p.list[[j]][subgroups[ss,2]]
#           inf.idx[ss] <- subgroup.n * pi * pj
#         }
#         cnt <- cnt+1
#         # g.list contains the most informative subgroups (take care of the right order!)
#         g.list[[cnt]] <- subgroups[inf.idx == max(inf.idx),,drop=F]  
#       }
#     }
    
    p <- p.list[[i]]
    x.est[,i] <- get.x.est(X[,i],models[i],p,
                           par2[i],group[,i])
        
    # for each RR variable: get quotient for 2-group models: 2 quotients
#     quotient.list[[i]] <- 0
#     for (gg in 1:length(unique(group[,i]))){
#       sel <- group[,i] == gg
#       quotient.list[[i]][gg] <- get.quotient(X[sel,i],models[i],p, par2[i], group=gg)
#     }
    
    # for 2-group models: separately
    sel <- rep(T,n)
    if (models[i] %in% c("SLD","UQTunknown")){
      if (p[2] > p[1]){
        group.sel[i] <- 2
        sel <- group[,i] == 2
      }else{
        sel <- group[,i] == 1
      }
    }
    quotient[i] <- get.quotient(X[sel,i],models[i],p, par2[i], group=group.sel[i])
  }  
 
  r <-  diag(m) 
#   cnt <- 0
  # calculate correlation separately for each combination of variables
  for (i in 1: (m-1)){
    for (k in (i+1):m){
      # all possible group combinations
      #       gg <- unique(group[,c(i,k)])
      # for multiple group models
      #       for (g in 1:nrow(gg)){
      sel <- group[,i] == group.sel[i] & group[,k] == group.sel[k]
      
      # for multiple group models: go through all subgroups listed in g.list!
#       cnt <- cnt +1
#       for (gg in 1:nrow(g.list[[cnt]])){
#         # select subgroup
#         g.i <- g.list[[cnt]][gg,1]
#         g.j <- g.list[[cnt]][gg,2]
#         sel.i <- group[,i] == g.i
#         sel.j <- group[,j] == g.j
#         sel <- sel.i & sel.j
#         
#         r.obs <- cor(x.est[sel,i],x.est[sel,j])
# #         print( (1+ quotient.list[[i]][g.i])*(1+ quotient.list[[j]][g.j]) )
#         r[i,j] <- r[i,j] +  r.obs * sqrt( (1+ quotient.list[[i]][g.i])*(1+ quotient.list[[j]][g.j]) ) / nrow(g.list[[cnt]])
#         print(r.obs * sqrt( (1+ quotient.list[[i]][g.i])*(1+ quotient.list[[j]][g.j]) ))
#         r[i,j] <- r.sign( r[i,j], models[i], p.list[[i]], models[j], p.list[[j]]) 
#         n.mat[i,j] <- n.mat[i,j] + sum(sel)
#       }
      n.mat[i,k] <- sum(sel)

      #         sel <- group[,i] == gg[g,i] &  group[,k] == gg[g,k]
      r.obs <- cor(x.est[sel,i],x.est[sel,k])
      r[i,k] <- r[i,k] +  r.obs * sqrt( (1+ quotient[i])*(1+ quotient[k]) )
      r[i,k] <- r.sign( r[i,k],models[i],p.list[[i]] , models[k],p.list[[k]]) 
      #       }
    }
  }
  
  
  r[lower.tri(r)] <- t(r)[lower.tri(r)]
  n.mat[lower.tri(n.mat)] <- t(n.mat)[lower.tri(n.mat)]
  rownames(r) <- colnames(X)
  colnames(r) <- colnames(X)
  rownames(n.mat) <- colnames(X)
  colnames(n.mat) <- colnames(X)
  # character vector including randomization probabilities
  p.char <- character(m)
  for (i in 1:m){
    if (models[i]!="direct"){
      for (k in 1:length(p.list[[i]])){
        p.char[i] <- paste(p.char[i],(p.list[[i]][k]))
      }
    } 
  }
  
  rSE.p <- matrix(NA,nrow=m, ncol=m, dimnames=list(colnames(X), colnames(X))) 
  r.Prob <- rSE.p ; rSE.n <- rSE.p;
  # parametric bootstrap: SE
  if(bs.n>0 && "p" %in% bs.type){
    for (i in 1:(m-1)){
      for(j in (i+1):m){
#         cat("variables",1,"and",2,": ",paste( models[c(i,j)]))
        ci <- c(1,1); cj <- c(1,1);
        if (models[i] == "SLD") ci <- c(par2[i],1)
        if (models[j] == "SLD") cj <- c(par2[j],1)
        #         if (models[i] == "CDM") ci <- rep(1-par2[i], 2)
        #         if (models[j] == "CDM") cj <- rep(1-par2[j], 2)
        groupRatios <- 2 - colMeans(group.mat)
        
        # how many RR variables in pairwise bootstrap?
        if(sum(models[c(i,j)] == "direct") <2){
          if(models[i] != "direct")
            RRi <- RRuni(X[,i], model =models[i], p=p.list[[i]], group = group.mat[,i],MLest = T)
          if(models[j] != "direct")
            RRj <- RRuni(X[,j], model =models[j], p=p.list[[j]], group = group.mat[,j],MLest = T)
          # two RRs
          if(sum(models[c(i,j)] == "direct") ==0){
            mcsim <- RRsimu(numRep=bs.n, n=n, pi=c(RRi$pi, RRj$pi), model = models[c(i,j)],
                            p=p.list[c(i,j)], rho=0, complyRates=list(ci, cj), sysBias=c(0,0),
                            groupRatio=groupRatios[c(i,j)], method="RRcor", MLest=T, nCPU=nCPU)
          }else{
            # one RR, one nonRR variable
            sel <- grep("direct", models[c(i,j)])
            idx <- ifelse(sel==2,i,j)  # select RR variable
            mcsim <- RRsimu(numRep=bs.n, n=n, pi=ifelse(idx==i, RRi$pi, RRj$pi), model = models[idx],
                            p=unlist(p.list[idx]), rho=0, complyRates=ifelse(rep(idx,2)==i,ci, cj), sysBias=c(0,0),
                            groupRatio=groupRatios[idx], method="RRcor", MLest=T, nCPU=nCPU)
          }
          rSE.p[i,j] <- mcsim$results["cor.RRcor","SE"]
          #           r.Prob[i,j] <- sum(mcsim$parEsts[,"cor.RRcor"]>r[i,j])/bootstrap 
        }
        
      }
    }
    rSE.p[lower.tri(rSE.p)] <- t(rSE.p)[lower.tri(rSE.p)]
    #     r.Prob[lower.tri(r.Prob)] <- t(r.Prob)[lower.tri(r.Prob)]
  }
  
# nonparametric bootstrap
  if(bs.n>0 && "n" %in% bs.type){
    getBoot <- function(rep){
      bs.ests <- array(NA, dim=c(m,m,rep))
      for (i in 1: bs.n){
        sel <- sample(1:n, n,T)
        try(bs.ests[,,i] <- RRcor(x=X[sel,],models=models,p.list=p.list, group=group.2g[sel,] )$r)
      }
      bs.ests
    }
    if (nCPU == 1){
      bs.ests <- getBoot(rep=bs.n)   
    }else{
#       require(doParallel, quietly=T)
      if (nCPU=="max"){
        try(nCPU <-  as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')))
        if (nCPU=="max") nCPU <- 2
      }
      cl.tmp =  makeCluster(nCPU) 
      registerDoParallel(cl.tmp, cores=nCPU)
      rep <- ceiling(bs.n/nCPU)
      bs.ests <- foreach(k=1:nCPU , .packages='RRreg') %dopar% {getBoot(rep) }
      bs.ests <- array(unlist(bs.ests), dim = c(m,m, rep*nCPU))
      stopCluster(cl.tmp)
    }
#     print(apply(bs.ests,1:2,mean))
    bs.n.NA <- sum(is.na(bs.ests[2,1,]))
    rSE.n <- apply(bs.ests,1:2,sd, na.rm=T)
    diag(rSE.n) <- NA
    dimnames(rSE.n) <- dimnames(rSE.p)
  }
  
  res <- list(r = r, call = "Randomized response correlation",
              models = models, p.list = p.list , p.char = p.char, varNames = colnames(X), n.mat=n.mat,
              n=n, bs.n=bs.n, bs.type=bs.type, rSE.p=rSE.p, rSE.n=rSE.n)
  class(res) <- "RRcor"
  res
}


#' @aliases RRcor
#' @method print RRcor
#' @export
print.RRcor<-function(x,...){
  cat("Call: \n")
  write(x$call,"")
  cat("\nRandomized response variables:\n")
  TAB <- cbind(Variable = x$varNames,
               RRmodel = x$models,
               p = x$p.char)
  rownames(TAB)=1:length(x$models)
  print(TAB, quote = FALSE)
  if (min(x$n.mat) != x$n){
    cat(paste("\nSample size N:\n"))
    print( x$n.mat)
  }else{
    cat(paste("\nSample size N =",x$n,"\n"))
  }
  cat("\nEstimated correlation matrix:\n")
  print( round(x$r,6))
  
  if(x$bs.n >0 &&  "p" %in% x$bs.type){
    cat("\nStandard errors from parametric bootstrap:\n")
    print( round(x$rSE.p,6), quote = FALSE) 
  }
  if(x$bs.n >0 &&  "n" %in% x$bs.type){
    cat("\nStandard errors from nonparametric bootstrap:\n")
    print( round(x$rSE.n,6), quote = FALSE)
  }
  
}



# 
# 
# summary.RRcor<-function(object,...){
#   zval <- object$pi/object$rSE
#   TAB <- cbind(Estimate = object$r,
#                StdErr = object$rSE,
#                z.value=zval,
#                p.value=pnorm(-abs(zval)))
#   rownames(TAB) <- "r"
#   res <- list(call=object$call,n=object$n,
#               coefficients=TAB)
#   class(res) <- "summary.RRcor"
#   return(res)
# }
# 
# 
# #' @aliases RRcor
# #' @method print summary.RRcor
# #' @export
# print.summary.RRcor<-function(x,...){
#   cat("Call:\n")
#   write(x$call,"")
#   cat("Sample size: ")
#   write(x$n,"")
#   cat("\n")
#   printCoefmat(x$coefficients)
# }

