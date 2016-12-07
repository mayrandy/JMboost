library(JMboost)
set.seed(123)
dat <- simJM(n = 400, n_i = 3, alpha = .5, beta = c(1,2,3), betals = c(2,3,1), betatimeind = 3,
            lambda = 0.6)

JMboost(y =dat$y, Xl = dat$Xl, Xls = dat$Xls, last = dat$last,
           delta = dat$delta, id = dat$id, time = dat$time, lambda = 1, alpha = 0.1,
           mstop_l = 10, mstop_ls = 10, step.length = .1, betatimeind = 3)
