
# Main function to carry out boosting!
#---------------------------------------------------
# IN:
# dep: y (longitudinal outcome including death observation (last)) # RENAME to y!
# Xl: x for longitudinal predictor
# Xls: x for shared predictor
# last:  0 obs contributing to longitudinal model 1 obs for shared model
# delta: censoring indicator (person had event or not)
# id: 1xN vector of subjects
# time: all time points (included in Xls if betatimeind != 0)
# lambda: starting value baseline hazard
# alpha: starting value association parameter
# sigma2: starting value for sigma2^2  (in the longitudianal part) REMOVE!
# mstop_l: stopping iteration for long. predictor
# mstop_ls stopping iteration for shard predictor
# step.length: for the boosting update (default 0.1?) SET!
# betatimeind: indicating which coefficient is fixed time effect [number]
# -----------------------------------------------------
# Out:
#
#
# --------

JMboost <- function(y, Xl, Xls, last, delta, id, time, lambda = 1, alpha = 0.1,
                       mstop_l, mstop_ls = NULL, step.length = 0.1, betatimeind = 0){

  sigma2 <- var(y) # starting value for sigma2
  if(is.null(mstop_ls)){mstop_ls <- mstop_l} # use mstop_l if ls missing
  mstop <- max(mstop_l, mstop_ls) # general mstop, perhaps RENAME

  n <- length(y) # number of all obs
  N <- length(unique(id)) # number of subjects

  #---------------------------------------------
  # code to construct random effects
  # needs df2lambda from mboost
  Xran <- matrix(ncol=2*N, nrow = length(y), data=0) # two columns per subject (slope and int), one row for each obs
  unid <- order(unique(id))
  id <- rep(unid, as.vector(table(id)))
  for(i in 1:N){
    Xran[which(id==as.character(i)),i] <- 1
    Xran[which(id==as.character(i)),N+i] <- time[which(id==as.character(i))]
  }
  Xrana <- Xran[,1:N]
  #lambdaran <- mboost:::df2lambda(Xrana, 4, weights=1)[2]
  lambdaran <- mboost_intern(Xrana, 4, weights=1, fun = "df2lambda")[2]
  Xranat <- t(Xrana)
  Sa <- solve(Xranat%*%Xrana + lambdaran*diag(N))%*%Xranat
  Xranb <- Xran[,c((N+1):(2*N))]
  #lambdaran <- mboost:::df2lambda(Xranb, 4, weights=1)[2]
  lambdaran <- mboost_intern(Xranb, 4, weights = 1, fun = "df2lambda")[2]
  Xranbt <- t(Xranb)
  Sb <- solve(Xranbt%*%Xranb + lambdaran*diag(N))%*%Xranbt
  gamma0 <- rep(0,N)
  gamma1 <- rep(0,N)
  #---------------------------------------------------

  #
  # if no variables for longitudinal predictor, we use an intercept
  if(is.null(Xl)){
    X <- as.matrix(rep(1, length(y)))
  }else{X <- cbind(1, Xl)
  }
  # Initializing
  int <- 0 # starting value for intercept
  beta <- rep(0, ncol(X)) # initializing coefficients
  beta[1] <- int

  betals <- rep(0, ncol(Xls)) #initializing coefficients
  if(length(betals)==1){betals <- t(betals)} # R can not handle it otherwise

  # storing matrices and vectors
  gamma0_mat <- matrix(0,ncol=mstop,nrow=N)
  gamma1_mat <- matrix(0,ncol=mstop,nrow=N)
  beta_mat <- matrix(0, ncol=mstop, nrow=ncol(X))
  betals_mat <- matrix(0, ncol=mstop, nrow=ncol(Xls))
  alphavec <- rep(0, mstop)
  lambdavec <- rep(0, mstop)
  sigma2vec <- rep(0, mstop)


  m <- 0 # iteration counter

  #---------------------------------

  ### Outer loop
  for(m in 1:mstop){


    ###############################################################
    #### S1  ######################################################
    ###############################################################
    if(m <= mstop_l){
      if(ncol(X)>1){ # more than just intercept
        f <- X%*%beta  # this is the predictor we are actually boosting (long.)

        u <- (y - Xran%*%c(gamma0,gamma1) - Xls%*%betals - f)/sigma2 # neg
        RSS <- numeric(ncol(X) - 1) # all without intc.

        for(j in 2:ncol(X)){ # inner loop
          J <- (j-1)
          assign(paste0("fit",J),lm(u[last==0]~X[last==0,j]))
          RSS[J] <- mean(resid(get(paste0("fit",J)))^2) # mean RSS
        }
        best <- which.min(RSS) # select best base-learner
        beta[best+1] <- beta[best+1] + step.length*coef(get(paste0("fit",best)))[2]  # update
        int <- int + step.length*coef(get(paste0("fit",best)))[1] # update intc.
      }else{ # only intercept
        f <- X%*%beta
        u <- (y - Xran%*%c(gamma0,gamma1) - Xls%*%betals - f)/sigma2
        b0 <- mean(u)
        int <-  int + step.length*b0
      }
    }
    # storage
    beta[1] <-int
    beta_mat[,m] <- beta

    ###############################################################
    #### S2  ######################################################
    ###############################################################
    #
    #
    if(m <= mstop_ls){

      # predictor for last obs
      f_surv <- Xran[last==1,]%*%c(gamma0,gamma1) + Xls[last==1,]%*%betals

      # predictor for all obs
      f <-  Xran%*%c(gamma0,gamma1) + Xls%*%betals

      # here starts the magic gradient
      if(betatimeind >0){ # if there exists a fixed time effect

        if(sum(gamma1)!=0 | betals[betatimeind]!=0){ #  if the current time coefs are not 0
          # this is ngradient for the last obs
          u_surv <- delta*alpha - lambda*(exp(alpha*f_surv) - exp(alpha*(f_surv - (gamma1 + betals[betatimeind])*time[last==1])))/(gamma1 + betals[betatimeind])
        }else{#  if the current time coefs are  0
          u_surv <- delta*alpha - alpha*lambda*exp(alpha*f_surv)*time[last==1]
        }
      }else{ # if there is no fixed time effect.
        if(sum(gamma1)!=0){
        u_surv <- delta*alpha - lambda*(exp(alpha*f_surv) - exp(alpha(f_surv-gamma1*time[last==1])))/gamma1
        }else{
          u_surv <- delta*alpha - alpha*lambda*exp(alpha*f_surv)*time[last==1]
        }
      }

      # longitudinal contribution to ngradient
      uls <- (y - X%*%beta - f)*(1/sigma2)
      # contribution of last obs
      uls[which(last==1)] <-  u_surv

       # this is HEPP stuff
      RSS <- numeric(ncol(Xls) + 1)
      assign(paste0("fit",1),(rbind(Sa, Sb)%*%uls)) # this is the estimation of random effects
      RSS[1] <- mean((uls-(Xran%*%get(paste0("fit",1))))^2) # mean RSS of random effect base-learner

      for(j in 1:ncol(Xls)){ # inner loop
        J <- j+1
        assign(paste0("fit",J),lm(uls~Xls[,j])) # always only info at last obs
        RSS[J] <- mean(resid(get(paste0("fit",J)))^2)  # mean RSS
      }
      best <- which.min(RSS) # select best base-learner

      # update of random effects (if selected)
      if(best==1){
      gamma0 <- gamma0+step.length*get(paste0("fit",best))[c(1:length(gamma0))]
      gamma1 <- gamma1+step.length*get(paste0("fit",best))[c((length(gamma0)+1):(2*length(gamma0)))]
      }else{
        betals[best-1] <- betals[best-1] + step.length*coef(get(paste0("fit",best)))[2]
        int <- int + step.length*coef(get(paste0("fit",best)))[1]
        beta[1] <-int
      }


      # current version
      f_surv <- Xran[last==1,]%*%c(gamma0,gamma1) + Xls[last==1,]%*%betals

      # optimize lambda
      lambda <- optimize(riskJM, alpha = alpha, f_surv = f_surv, delta = delta, gamma1 = gamma1,
                         betals = betals, betatimeind = betatimeind, time = time[last==1],
                         interval=c(0,100))$minimum
      #
      if(m > 20){
        alpha <- optimize(riskJM, lambda = lambda, f_surv = f_surv,  delta = delta, gamma1 = gamma1,
                          betals = betals, betatimeind = betatimeind, time = time[last==1],
                          interval=c(-100,100))$minimum
      }else{
        alpha <- optimize(riskJM, lambda = lambda, f_surv = f_surv,  delta = delta, gamma1 = gamma1,
                          betals = betals, betatimeind = betatimeind, time = time[last==1],
                          interval=c(-1,1))$minimum}

    }

    # storage
    betals_mat[,m] <- betals
    gamma0_mat[,m] <- gamma0
    gamma1_mat[,m] <- gamma1



    ###############################################################
    #### S3  ######################################################
    ###############################################################
    # currently only sigma2

    f_ls <- Xls%*%betals + Xran%*%c(gamma0,gamma1)
    f_l <- X%*%beta

    risk2 <- function(y, f_l, f_ls, sigma2){
                     sum( - dnorm(y, mean = (f_l + f_ls), sd = sqrt(sigma2), log = TRUE))}
    sigma2 <- optimize(risk2, y=y[last==0], f_l=f_l[last==0], f_ls =f_ls[last==0], interval=c(0,100))$min

    # storage
    lambdavec[m] <- lambda
    alphavec[m] <- alpha
    sigma2vec[m] <- sigma2

  }


  structure(list(gamma1_mat=gamma1_mat,gamma0_mat=gamma0_mat, beta_mat= beta_mat, beta = beta,
                 betals_mat= betals_mat, betals = betals,
                 gamma0=gamma0, gamma1=gamma1, alpha= alpha, alphavec =alphavec,
                 lambda=lambda, lambdavec = lambdavec, sigma2vec=sigma2vec, sigma2=sigma2))
}


