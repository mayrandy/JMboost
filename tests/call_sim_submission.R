library(JMboost)
setup <- "small"
n <- 500
n2 <- 1000
n_i <- 5
alpha <- .5
beta <- c(2,1,-2)
betals <- c(1,-2,1)
lambda <- .1
like_cv <- function(delta, lambda, alpha, time, 
                    sigma, y, betatimeind, datals, betals, last, data, beta) {
  beta_t <- betals[betatimeind]
  etals_m <- (datals[,-betatimeind]%*%betals[-betatimeind])[last==1]
  etals <- (datals%*%betals)[last==0]
  etal <- (cbind(1,data)%*%beta)[last==0]
  if(beta_t!=0){int <- lambda*(1/(alpha*beta_t)*(exp(alpha*etals_m + alpha*beta_t*time) - exp(alpha*etals_m)))}else{
    int <- lambda*time*exp(alpha*etals_m)  
  }
  surv <- delta*log(lambda) + 
    delta*exp(alpha*etals_m + alpha*beta_t*time) - int
  
  long <- log(1/sqrt(2*pi*sigma)) - (y-etals - etal)^2/(2*sigma)
  like <- sum(surv) + sum(long)
  like
}

bom <- list()
likemat <- list()

for(k in 1:2){
set.seed(k)  
si <- simJM(n, n_i, alpha, beta, betals,  betatimeind = 3,lambda,4,4)
si2 <- simJM (n2, n_i, alpha, beta, betals,  betatimeind = 3,lambda, 4,4)
bo <-list()
mseq <- seq(30,300, length=10)
likemat[[k]] <- matrix(nrow=10, ncol=10)

for(m_akt_ls in 1:2){
  for(m_akt_l in 1:2){
    
    bo <- JMboost(y =si$y, X = si$X, Xls = si$Xls, last = si$last,
                  delta = si$delta, id = si$id, time = si$time, lambda = .1, alpha = 0.1,
                  mseq[m_akt_l],mstop_ls= mseq[m_akt_ls], step.length = .1, betatimeind = 3) 
    delta2 <- si2$delta[si2$id]
    like <- like_cv(delta2[si2$last==1], bo$lambda, bo$alpha, 
                    si2$time[si2$last==1], bo$sigma2, si2$y[si2$last==0], 
                    3, si2$Xls, bo$betals, si2$last, si2$X, bo$beta)      
    likemat[[k]][m_akt_l,m_akt_ls] <- like
   
  }
  
}

m <- which(likemat == max(likemat, na.rm=T), arr.ind = TRUE)[1,]
bom[[k]] <- JMboost(y =si$y, X = si$X, Xls = si$Xls, last = si$last,
               delta = si$delta, id = si$id, time = si$time, lambda = .1, alpha = 0.1,
               mseq[m[1]],mstop_ls= mseq[m[2]], step.length = .1, betatimeind = 3) 
  
}
bom[[k]]$beta
bom[[k]]$alpha
