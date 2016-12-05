## n: #{persons}
## n_i: #{observations person} (BEFORE censoring)
## beta: beta vector for the longitudinal predictor
## betals: beta vector for the shared predictor
## betatimeind: is there a time dependent effect in betals? Which one?
## alpha: association parameter
## lambda: constant baseline hazard

simJM <- function(n, n_i, alpha, beta, betals = 0, betatimeind = 0, lambda){

  ### generating id
  id <- rep(1:n, each=n_i)
  ### simulating time points
  ## mimicking yearly observations
  day_in_year <- sample(1:365, n*n_i, replace=TRUE)
  ## dann immer vom Jahresanfang aus aufsummieren?
  time <- rep(seq(0,(n_i-1)*365, 365), n)
  time <- time + day_in_year
  ## norming between 0 and 1
  time <- time/(n_i*365)
  ### generating eta and eta_ls
  ## eta_ls: random intercept + random Slope
  rint <- rnorm(n,0,.01)
  rslope <- rnorm(n,0,.01)
  R_mean <- cbind(rep(rint, each=n_i),rep(rslope, each=n_i))
  eta_ls_mean <- rowSums(cbind(1,time)*R_mean)

  if(all(betals==0)==FALSE){
    if(betatimeind!=0){
      betatime <- betals[betatimeind]
      betals <- betals[-betatimeind]
    }else{betatime <- 0}
    Xls <- matrix(nrow=n*n_i, ncol=length(betals))
    for(j in 1:ncol(Xls)){
      sa <- runif(n)
      Xls[,j] <- rep(sa, each=n_i)
      Xls[,j] <- (Xls[,j]-mean(Xls[,j]))/sd(Xls[,j])
    }
    ### is there a time dependent variable?
    eta_ls_mean = eta_ls_mean + Xls%*%betals + betatime*time
  }else{Xls=NULL
  betatime =0}
  ## eta_l: simple linear model
  X <- matrix(nrow=n*n_i, ncol=length(beta))
  X[,1] <- 1
  if(ncol(X)>1){
    for(j in 2:ncol(X)){
      sa <- sample(seq(0,3),1)
      X[,j] <- runif(n*n_i, 0,1)
      X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
    }}
  eta_l_mean <- X%*%beta
  ## longitudinal measurements
  y <- rnorm(n*n_i, eta_ls_mean + eta_l_mean, sqrt(.1))

  ### Event times
  ## Probabilities from event times
  pred_surv <- alpha*(eta_ls_mean)
  if(all(betals==0)){F_death <- 1-exp(-lambda*(exp(pred_surv)-exp(alpha*(R_mean[,1])))/(alpha*(R_mean[,2]+betatime)))}else{
    F_death <- 1-exp(-lambda*(exp(pred_surv)-exp(alpha*(R_mean[,1]+Xls%*%betals)))/(alpha*(R_mean[,2]+betatime)))}

  prob_death2 <-  matrix(nrow=n_i, ncol=n,
                         data=F_death)
  prob_death_ext <- rbind(0, prob_death2)
  prob_death_ext <- prob_death_ext[-nrow(prob_death_ext),]

  prob_death <- (prob_death2 - prob_death_ext)/(1-prob_death_ext)
  u <- matrix(nrow=nrow(prob_death), ncol=ncol(prob_death), runif(n*n_i))
  death <- matrix(nrow=n_i, ncol=n,
                  data=u<prob_death)
  time_mat <- matrix(nrow=n_i, data=time)
  helpmat <- (rbind(0,time_mat[-n_i,])-u*(rbind(0,time_mat[-n_i,])-time_mat)/prob_death)*death
  time_mat[helpmat!=0] <- helpmat[helpmat!=0]
  helpmat2 <- as.vector(helpmat)
  time[helpmat2!=0] <- helpmat2[helpmat2!=0]
  helpmat[helpmat==0] <- 1000
  minna <- function(x){min(x, na.rm=T)}
  death_time <- cbind(unique(id),apply(helpmat, 2, minna))
  c_i <- colSums(death)==0
  delta <- c_i==0
  time_mat_help <- time_mat*as.numeric(t(t(time_mat)<=(death_time[,2])))
  time_zero <- as.vector(time_mat_help)


  if(length(which(time_zero==0))>0){
    id <- id[-which(time_zero==0)]
    y <- y[-which(time_zero==0)]
    X <- X[-which(time_zero==0),]
    Xls <- Xls[-which(time_zero==0),]
    time <- time[-which(time_zero==0)]
    R_mean <- R_mean[-which(time_zero==0),]
    prob_death <- prob_death[-which(time_zero==0)]
    pred_surv <- pred_surv[-which(time_zero==0)]}
  c_i <- c_i
  last <- rep(FALSE, length(time))
  id_un <- unique(id)
  for(i in 1:length(id_un)){
    last[max(which(id==id_un[i]))] <- TRUE
  }

  if(betatimeind!=0){Xls <- cbind(Xls, time)}
  X <- X[,-1]


  return(list(
    #### longitudinal outcome
    "y" = y,
    #### longitudinal predictor fixed effect covariates
    "X" = X,
    #### shared predictor fixed effect covariates
    "Xls"=Xls,
    ### Values for the random effects
    "R_mean"= R_mean,
    #### id indicator
    "id"=id,
    #### Time vector
    "time"=time,
    #### Censoring indicator
    "delta" = delta,
    #### Event time (set to 1000, if the invidual is censored)
    "event_time"=death_time,
    ####Indicator per observation, if it is the last observation for the individual
    "last" = last
  )
  )
}
