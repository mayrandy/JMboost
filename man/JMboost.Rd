\name{JMboost}
\alias{JMboost}

\title{
Boosting joint models
}
\description{
Provides an implementation to boost joint models for longitudinal and time-to-event outcomes.
}
\usage{
JMboost(y, Xl, Xls, last, delta, id, time,
        lambda = 1, alpha = 0.1, mstop_l,
        mstop_ls = NULL, step.length = 0.1,
        betatimeind = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Longitudinal outcome including observation at or before event (last).
}
  \item{Xl}{
covariate matrix for longitudinal predictor
}
  \item{Xls}{
covariate matrix for shared predictor
}
  \item{last}{
vector indicating if observation refers to longitduinal (FALSE) or shared predictor (TRUE); i.e. if it is the patient's last observation before event.
}
  \item{delta}{
censoring indicator (person had event or not)
}
  \item{id}{
1xN vector of subjects
}
  \item{time}{
all time points (included in Xls if betatimeind != 0)
}
  \item{lambda}{
starting value baseline hazard
}
  \item{alpha}{
starting value association parameter
}
  \item{mstop_l}{
Stopping iteration longitudinal predictor.
}
  \item{mstop_ls}{
Stopping iteration shared predictor.
}
  \item{step.length}{
Step length for boosting updates, typical value is 0.1 (default).
}
  \item{betatimeind}{
  indicating which coefficient in shared predictor is time variable (number),
}
}
\details{
Implements a gradient boosting algorithn to fit joint models considering a longitudinal predictor (w.r.t outcome y) and a shared predictor (affecting the longitudinal outcome and the time-to-event outcome). The relation between the two outcomes is described via the association parameter alpha. For more details, see Waldmann et al. (2016).

}
\references{
Waldmann, E., Taylor-Robinson, D., Klein, N., Kneib, T., Pressler, T., Schmid, M., & Mayr, A. (2016). Boosting Joint Models for Longitudinal and Time-to-Event Data. arXiv preprint arXiv:1609.02686.}
\examples{
set.seed(123)
dat <- simJM(n = 400, n_i = 3, alpha = .5,
            beta = c(1,2,3), betals = c(2,3,1),
            betatimeind = 3, lambda = 0.6)

j1 <- JMboost(y = dat$y, Xl = dat$X, Xls = dat$Xls,
              last = dat$last, delta = dat$delta,
              id = dat$id, time = dat$time, lambda = 1, alpha = 0.1,
              mstop_l = 100, mstop_ls = 100, step.length = 0.1,
              betatimeind = 3)
}
