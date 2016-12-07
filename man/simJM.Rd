\name{simJM}
\alias{simJM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data frames for joint modelling of longitudinal and time-to-event outcome.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simJM(n = 100, n_i = 5, alpha = 0.5,
      beta, betals = 0, betatimeind = 0, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of individuals.
}
  \item{n_i}{
number of observations per individual (before censoring).
}
  \item{alpha}{
association parameter.
}
  \item{beta}{
coefficient beta vector for longitudinal outcome.
}
  \item{betals}{
coefficient beta vector for shared predictor.
}
  \item{betatimeind}{
index of time-variable in shared predictor.
}
  \item{lambda}{
constant bazeline hazard.
}
}
\details{
Simulates a data frame for joint models, considering a longitudinal and a shared predictor. For more details, see the Appendix of Waldmann et al. (2016).

}
\value{
A data frame.
}\references{
Waldmann, E., Taylor-Robinson, D., Klein, N., Kneib, T., Pressler, T., Schmid, M., & Mayr, A. (2016). Boosting Joint Models for Longitudinal and Time-to-Event Data. arXiv preprint arXiv:1609.02686.}
\examples{
set.seed(123)
dat <- simJM(n = 400, n_i = 3, alpha = .5,
            beta = c(1,2,3), betals = c(2,3,1),
            betatimeind = 3, lambda = 0.6)

j1 <- JMboost(y = dat$y, Xl = dat$Xl, Xls = dat$Xls,
              last = dat$last, delta = dat$delta,
              id = dat$id, time = dat$time, lambda = 1, alpha = 0.1,
              mstop_l = 100, mstop_ls = 100, step.length = 0.1,
              betatimeind = 3)
}
