########################

# risk: Survival Likelihood for optimization of alpha and lambda
#---------------------------------------------------
# IN: 
# alpha: association parameter
# lamda: const. baseline hazard
# f_surv: shared predictor at last observation
# delta: censoring indicator 1xn
# gamma1: random slope
# betals: shared coefficients
# betatimeind: indicating which coefficient is fixed time effect [number]
# if betatimeind == 0 then we do not have any 
# time: last observation / time of death (1xn)
#---------------------------------------------------
# OUT:
# negative log likelihood (sum over obs)


riskJM <- function(alpha, lambda, f_surv, delta, gamma1, betals, betatimeind, time){
  if(betatimeind!=0){ # with fixed time effect ... 
    if(sum(gamma1)!=0 | betals[betatimeind]!=0) # do we actually have any time effect (random or fixed) 
      # at current iteration? YES
    { res <- -(sum(delta[delta==1]*log(lambda) + alpha*f_surv[delta==1]) -  
                 sum( lambda*(exp(alpha*f_surv) - exp(alpha*(f_surv - gamma1*time - betals[betatimeind]*time))) /
                        (alpha*(gamma1 + betals[betatimeind])) ))
    }else{  # at current iteration NO (random or fixed) time effect 
      res <- -(sum(delta[delta==1]*log(lambda) + alpha*f_surv[delta==1]) - sum(lambda*exp(alpha*f_surv)*time))
    }
  }else{ # needs to be fixed ! # 
    res <- -(sum(delta[delta==1]*log(lambda) +alpha*f_surv[delta==1]) -
               sum(lambda*(exp(alpha*f_surv)-exp(alpha*(f_surv-gamma1*time)))/(alpha*gamma1)))}
  return(res)
}


