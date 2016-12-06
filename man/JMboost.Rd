\name{JMboost}
\alias{JMboost}

\title{
Boosting joint models
}
\description{
Provides an implementation to boost joint models for longitudinal and time-to-event outcomes.
}
\usage{
JMboost(y, Xl, Xls, last, delta, id, time, lambda = 1, alpha = 0.1, m_stop_l, m_stop_ls = NULL, step.length, betatimeind = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{Xl}{
%%     ~~Describe \code{Xl} here~~
}
  \item{Xls}{
%%     ~~Describe \code{Xls} here~~
}
  \item{last}{
%%     ~~Describe \code{last} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{m_stop_l}{
%%     ~~Describe \code{m_stop_l} here~~
}
  \item{m_stop_ls}{
%%     ~~Describe \code{m_stop_ls} here~~
}
  \item{step.length}{
%%     ~~Describe \code{step.length} here~~
}
  \item{betatimeind}{
%%     ~~Describe \code{betatimeind} here~~
}
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Waldmann, E., Taylor-Robinson, D., Klein, N., Kneib, T., Pressler, T., Schmid, M., & Mayr, A. (2016). Boosting Joint Models for Longitudinal and Time-to-Event Data. arXiv preprint arXiv:1609.02686.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, Xl, Xls, last, delta, id, time, lambda = 1, alpha = 0.1,
    m_stop_l, m_stop_ls = NULL, step.length, betatimeind = 0)
{
    sigma2 <- var(y)
    if (is.null(m_stop_ls)) {
        m_stop_ls <- m_stop_l
    }
    m_stop <- max(m_stop_l, m_stop_ls)
    n <- length(y)
    N <- length(unique(id))
    Xran <- matrix(ncol = 2 * N, nrow = length(y), data = 0)
    unid <- order(unique(id))
    id <- rep(unid, as.vector(table(id)))
    for (i in 1:N) {
        Xran[which(id == as.character(i)), i] <- 1
        Xran[which(id == as.character(i)), N + i] <- time[which(id ==
            as.character(i))]
    }
    Xrana <- Xran[, 1:N]
    lambdaran <- mboost_intern(Xrana, 4, weights = 1, fun = "df2lambda")[2]
    Xranat <- t(Xrana)
    Sa <- solve(Xranat \%*\% Xrana + lambdaran * diag(N)) \%*\% Xranat
    Xranb <- Xran[, c((N + 1):(2 * N))]
    lambdaran <- mboost_intern(Xranb, 4, weights = 1, fun = "df2lambda")[2]
    Xranbt <- t(Xranb)
    Sb <- solve(Xranbt \%*\% Xranb + lambdaran * diag(N)) \%*\% Xranbt
    gamma0 <- rep(0, N)
    gamma1 <- rep(0, N)
    if (is.null(Xl)) {
        X <- as.matrix(rep(1, length(y)))
    }
    else {
        X <- cbind(1, Xl)
    }
    int <- 0
    beta <- rep(0, ncol(X))
    beta[1] <- int
    betals <- rep(0, ncol(Xls))
    if (length(betals) == 1) {
        betals <- t(betals)
    }
    gamma0_mat <- matrix(0, ncol = m_stop, nrow = N)
    gamma1_mat <- matrix(0, ncol = m_stop, nrow = N)
    beta_mat <- matrix(0, ncol = m_stop, nrow = ncol(X))
    betals_mat <- matrix(0, ncol = m_stop, nrow = ncol(Xls))
    alphavec <- rep(0, m_stop)
    lambdavec <- rep(0, m_stop)
    sigma2vec <- rep(0, m_stop)
    m <- 0
    for (m in 1:m_stop) {
        if (m <= m_stop_l) {
            if (ncol(X) > 1) {
                f <- X \%*\% beta
                u <- (y - Xran \%*\% c(gamma0, gamma1) - Xls \%*\%
                  betals - f)/sigma2
                RSS <- numeric(ncol(X) - 1)
                for (j in 2:ncol(X)) {
                  J <- (j - 1)
                  assign(paste0("fit", J), lm(u[last == 0] ~
                    X[last == 0, j]))
                  RSS[J] <- mean(resid(get(paste0("fit", J)))^2)
                }
                best <- which.min(RSS)
                beta[best + 1] <- beta[best + 1] + step.length *
                  coef(get(paste0("fit", best)))[2]
                int <- int + step.length * coef(get(paste0("fit",
                  best)))[1]
            }
            else {
                f <- X \%*\% beta
                u <- (y - Xran \%*\% c(gamma0, gamma1) - Xls \%*\%
                  betals - f)/sigma2
                b0 <- mean(u)
                int <- int + step.length * b0
            }
        }
        beta[1] <- int
        beta_mat[, m] <- beta
        if (m <= m_stop_ls) {
            f_surv <- Xran[last == 1, ] \%*\% c(gamma0, gamma1) +
                Xls[last == 1, ] \%*\% betals
            f <- Xran \%*\% c(gamma0, gamma1) + Xls \%*\% betals
            if (betatimeind > 0) {
                if (sum(gamma1) != 0 | betals[betatimeind] !=
                  0) {
                  u_surv <- delta * alpha - lambda * (exp(alpha *
                    f_surv) - exp(alpha * (f_surv - (gamma1 +
                    betals[betatimeind]) * time[last == 1])))/(gamma1 +
                    betals[betatimeind])
                }
                else {
                  u_surv <- delta * alpha - alpha * lambda *
                    exp(alpha * f_surv) * time[last == 1]
                }
            }
            else {
                if (sum(gamma1) != 0) {
                  u_surv <- delta * alpha - lambda * (exp(alpha *
                    f_surv) - exp(alpha(f_surv - gamma1 * time[last ==
                    1])))/gamma1
                }
                else {
                  u_surv <- delta * alpha - alpha * lambda *
                    exp(alpha * f_surv) * time[last == 1]
                }
            }
            uls <- (y - X \%*\% beta - f) * (1/sigma2)
            uls[which(last == 1)] <- u_surv
            RSS <- numeric(ncol(Xls) + 1)
            assign(paste0("fit", 1), (rbind(Sa, Sb) \%*\% uls))
            RSS[1] <- mean((uls - (Xran \%*\% get(paste0("fit",
                1))))^2)
            for (j in 1:ncol(Xls)) {
                J <- j + 1
                assign(paste0("fit", J), lm(uls ~ Xls[, j]))
                RSS[J] <- mean(resid(get(paste0("fit", J)))^2)
            }
            best <- which.min(RSS)
            if (best == 1) {
                gamma0 <- gamma0 + step.length * get(paste0("fit",
                  best))[c(1:length(gamma0))]
                gamma1 <- gamma1 + step.length * get(paste0("fit",
                  best))[c((length(gamma0) + 1):(2 * length(gamma0)))]
            }
            else {
                betals[best - 1] <- betals[best - 1] + step.length *
                  coef(get(paste0("fit", best)))[2]
                int <- int + step.length * coef(get(paste0("fit",
                  best)))[1]
                beta[1] <- int
            }
            f_surv <- Xran[last == 1, ] \%*\% c(gamma0, gamma1) +
                Xls[last == 1, ] \%*\% betals
            lambda <- optimize(riskJM, alpha = alpha, f_surv = f_surv,
                delta = delta, gamma1 = gamma1, betals = betals,
                betatimeind = betatimeind, time = time[last ==
                  1], interval = c(0, 100))$minimum
            if (m > 20) {
                alpha <- optimize(riskJM, lambda = lambda, f_surv = f_surv,
                  delta = delta, gamma1 = gamma1, betals = betals,
                  betatimeind = betatimeind, time = time[last ==
                    1], interval = c(-100, 100))$minimum
            }
            else {
                alpha <- optimize(riskJM, lambda = lambda, f_surv = f_surv,
                  delta = delta, gamma1 = gamma1, betals = betals,
                  betatimeind = betatimeind, time = time[last ==
                    1], interval = c(-1, 1))$minimum
            }
        }
        betals_mat[, m] <- betals
        gamma0_mat[, m] <- gamma0
        gamma1_mat[, m] <- gamma1
        f_ls <- Xls \%*\% betals + Xran \%*\% c(gamma0, gamma1)
        f_l <- X \%*\% beta
        risk2 <- function(y, f_l, f_ls, sigma2) {
            sum(-dnorm(y, mean = (f_l + f_ls), sd = sqrt(sigma2),
                log = TRUE))
        }
        sigma2 <- optimize(risk2, y = y[last == 0], f_l = f_l[last ==
            0], f_ls = f_ls[last == 0], interval = c(0, 100))$min
        lambdavec[m] <- lambda
        alphavec[m] <- alpha
        sigma2vec[m] <- sigma2
    }
    structure(list(gamma1_mat = gamma1_mat, gamma0_mat = gamma0_mat,
        beta_mat = beta_mat, beta = beta, betals_mat = betals_mat,
        betals = betals, gamma0 = gamma0, gamma1 = gamma1, alpha = alpha,
        alphavec = alphavec, lambda = lambda, lambdavec = lambdavec,
        sigma2vec = sigma2vec, sigma2 = sigma2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
