beta <- matrix(nrow=100, ncol=length(bom[[1]]$beta))
betals <- matrix(nrow=100, ncol=length(bom[[1]]$betals))
alpha <- rep(NA, 100)
lambda <- rep(NA,100)
sigma <- rep(NA,100)
m <- matrix(nrow=100, ncol=2)
for(j in 1:2)
{
  beta[j,] <- bom[[j]]$beta
  betals[j,] <- bom[[j]]$betals
  lambda[j] <- bom[[j]]$lambda
  alpha[j] <- bom[[j]]$alpha
  sigma[j] <- bom[[j]]$sigma2
  m[j,] <-which(likemat[[j]] == max(likemat[[j]], na.rm=T), arr.ind = TRUE)[1,]
}

m_stop <- round(cbind(c(summary(mseq[m[,1]])[c(4)],
                        summary(mseq[m[,2]])[c(4)])))

paste(round(colMeans(beta[,1:3]),3),"(", round(apply(beta[,1:3], 2,sd),3),")", sep="")
paste(round(colMeans(betals[,1:3]),3),"(", round(apply(betals[,1:3], 2,sd),3),")", sep="")
beta_prop <- colSums(beta!=0)/100
paste(round(mean(beta_prop[-c(1:3)]),3), "(", round(sd(beta_prop[-c(1:3)]),3), ")", sep="")
beta_prop[1:3]
betals_prop <- colSums(betals!=0)/100
paste(round(mean(betals_prop[-c(1:3)]),3), "(", round(sd(betals_prop[-c(1:3)]),3), ")", sep="")
betals_prop[1:3]
round(mean(alpha), 3)
round(sd(alpha), 3)


if(ncol(beta_gr)<20){
beta_gr <- cbind(beta, betals)
boxplot(beta_gr, las=1, axes=FALSE, ylim=c(-2.5,2.5), col="lightgrey", border="grey",
        boxwex=c(rep(.85, 3), rep(.4,4), rep(.85,3), rep(.4,4)),
        at=c(1:3, seq(4,5.5,length=4),seq(6.5,8.5,1), seq(9.5,11,length=4)))
axis(2, las=1)
axis(4, las=1)
axis(1, labels = c("2","1", "-2"), at=c(1:3))
axis(1, labels = c("1","-2", "1"), at=seq(6.5,8.5,1))
mtext("0", side=1, line=1, cex=.8, at=c(4.75,10.25))
mtext(expression(paste(beta[l])), side=3, line=1, cex=1, at=3.5)
mtext(expression(paste(beta[ls])), side=3, line=1, cex=1, at=9)
segments(x0=0.5, x1=1.5, y0=2)
segments(x0=1.5, x1=2.5, y0= 1)
segments(x0=2.5, x1=3.5, y0=-2)
segments(x0=3.5, x1=5.975, y0=0)
segments(x0=6.025, x1=7, y0=1)
segments(x0=7, x1=8, y0= -2)
segments(x0=8, x1=9, y0=1)
segments(x0=9, x1=11.5, y0=0)
abline(v=c(seq(.5,3.5,1),5.98,6.02, seq(7,9,1),11.5), col="lightgrey", lty=2)
box()
}else{
  beta_gr <- cbind(beta[,1:3], rowMeans(beta[,-c(1:3)], na.rm=T), 
                        betals[,1:3], rowMeans(betals[,-c(1:3)], na.rm=T))
 boxplot(beta_gr, las=1, axes=FALSE, ylim=c(-2.5,2.5), col="lightgrey", border="grey",
          boxwex=c(rep(.8, 3), .8, rep(.8,3), .8),
          at=c(1:8))
  axis(2, las=1)
  axis(4, las=1)
  axis(1, labels = c("2","1", "-2", "0", "1","-2", "1", "0"), at=c(1:8))
  mtext(expression(paste(beta[l])), side=3, line=1, cex=1, at=2.5)
  mtext(expression(paste(beta[ls])), side=3, line=1, cex=1, at=6.5)
  segments(x0=0.5, x1=1.5, y0=2)
  segments(x0=1.5, x1=2.5, y0= 1)
  segments(x0=2.5, x1=3.5, y0=-2)
  segments(x0=3.5, x1=4.48, y0=0)
  segments(x0=4.52, x1=5.5, y0=1)
  segments(x0=5.5, x1=6.5, y0= -2)
  segments(x0=6.5, x1=7.5, y0=1)
  segments(x0=7.5, x1=8.5, y0=0)
  abline(v=c(seq(.5,3.5,1), 4.48,4.52,seq(5.5,8.5,1)), col="lightgrey", lty=2)
  box()
  
}

boxplot(alpha, las=1, ylim= c(0,1), axes=FALSE,
        main=expression(alpha), col="lightgrey", border="grey")
segments(x0=0.5, x1=1.5, y0=.5)
abline(v=seq(.5,3.5,1), 
       col="lightgrey", lty=2)
axis(2, las=1)
box()
