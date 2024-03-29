% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hnlmix.r
\name{hnlmix}
\alias{hnlmix}
\alias{deviance.hnlm}
\alias{residuals.hnlm}
\alias{fitted.hnlm}
\alias{print.hnlm}
\title{Generalized Nonlinear Regression using h-likelihood for a Random Parameter}
\usage{
hnlmix(
  y = NULL,
  distribution = "normal",
  mixture = "normal",
  random = NULL,
  nest = NULL,
  mu = NULL,
  shape = NULL,
  linear = NULL,
  pmu = NULL,
  pshape = NULL,
  pmix = NULL,
  prandom = NULL,
  delta = 1,
  common = FALSE,
  envir = parent.frame(),
  print.level = 0,
  typsize = abs(p),
  ndigit = 10,
  gradtol = 1e-05,
  stepmax = 10 * sqrt(p \%*\% p),
  steptol = 1e-05,
  iterlim = 100,
  fscale = 1,
  eps = 1e-04,
  points = 5
)
}
\arguments{
\item{y}{A response vector of uncensored data, a two column matrix for
binomial data or censored data, with the second column being the censoring
indicator (1: uncensored, 0: right censored, -1: left censored), or an
object of class, \code{response} (created by
\code{\link[rmutil]{restovec}}) or \code{repeated} (created by
\code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
\code{repeated} data object contains more than one response variable, give
that object in \code{envir} and give the name of the response variable to
be used here.}

\item{distribution}{The distribution for the response: binomial, beta
binomial, double binomial, mult(iplicative) binomial, Poisson, negative
binomial, double Poisson, mult(iplicative) Poisson, gamma count, Consul
generalized Poisson, logarithmic series, geometric, normal, inverse Gauss,
logistic, exponential, gamma, Weibull, extreme value, Cauchy, Pareto,
Laplace, Levy, beta, simplex, or two-sided power. (For definitions of
distributions, see the corresponding [dpqr]distribution help.)}

\item{mixture}{The mixing distribution for the random parameter (whose
initial values are supplied in \code{prandom}): normal, logistic, inverse
Gauss, gamma, inverse gamma, Weibull, or beta. The first two have zero
location parameter, the next three have unit location parameter, and the
last one has location parameter set to 0.5.}

\item{random}{The name of the random parameter in the \code{mu} formula.}

\item{nest}{The cluster variable classifying observations by the unit upon
which they were observed. Ignored if \code{y} or \code{envir} has class,
\code{response} or \code{repeated}.}

\item{mu}{A user-specified formula containing named unknown parameters,
giving the regression equation for the location parameter. This may contain
the keyword, \code{linear} referring to a linear part.}

\item{shape}{A user-specified formula containing named unknown parameters,
giving the regression equation for the shape parameter. This may contain
the keyword, \code{linear} referring to a linear part. If nothing is
supplied, this parameter is taken to be constant. This parameter is the
logarithm of the usual one.}

\item{linear}{A formula beginning with ~ in W&R notation, specifying the
linear part of the regression function for the location parameter or list
of two such expressions for the location and/or shape parameters.}

\item{pmu}{Vector of initial estimates for the location parameters. These
must be supplied either in their order of appearance in the formula or in a
named list.}

\item{pshape}{Vector of initial estimates for the shape parameters. These
must be supplied either in their order of appearance in the expression or
in a named list.}

\item{pmix}{If NULL, this parameter is estimated from the variances. If a
value is given, it is taken as fixed.}

\item{prandom}{Either one estimate of the random effects or one for each
cluster (see \code{nest}), in which case the last value is not used. If the
location parameter of the mixing distribution is zero, the last value is
recalculated so that their sum is zero; if it is unity, they must all be
positive and the last value is recalculated so that the sum of their
logarithms is zero; if it is 0.5, they must all lie in (0,1) and the last
value is recalculated so that the sum of their logits is zero.}

\item{delta}{Scalar or vector giving the unit of measurement (always one
for discrete data) for each response value, set to unity by default. For
example, if a response is measured to two decimals, \code{delta=0.01}. If
the response is transformed, this must be multiplied by the Jacobian. The
transformation cannot contain unknown parameters. For example, with a log
transformation, \code{delta=1/y}. (The delta values for the censored
response are ignored.)}

\item{common}{If TRUE, the formulae with unknowns for the location and
shape have names in common. All parameter estimates must be supplied in
\code{pmu}.}

\item{envir}{Environment in which model formulae are to be interpreted or a
data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
name of the response variable should be given in \code{y}. If \code{y} has
class \code{repeated}, it is used as the environment.}

\item{print.level}{Arguments for nlm.}

\item{typsize}{Arguments for nlm.}

\item{ndigit}{Arguments for nlm.}

\item{gradtol}{Arguments for nlm.}

\item{stepmax}{Arguments for nlm.}

\item{steptol}{Arguments for nlm.}

\item{iterlim}{Arguments for nlm.}

\item{fscale}{Arguments for nlm.}

\item{eps}{Arguments for nlm.}

\item{points}{Arguments for nlm.}
}
\value{
A list of class \code{hnlmix} is returned that contains all of the
relevant information calculated, including error codes.

The two variances and shrinkage estimates of the random effects are
provided.
}
\description{
\code{hnlmix} fits user-specified nonlinear regression equations to one or
both parameters of the common one and two parameter distributions. One
parameter of the location regression is random with some specified mixing
distribution.
}
\details{
It is recommended that initial estimates for \code{pmu} and \code{pshape}
be obtained from \code{gnlr}.

These nonlinear regression models must be supplied as formulae where
parameters are unknowns. (See \code{\link[rmutil]{finterp}}.)
}
\examples{

dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
#y <- rgamma(20,2+0.3*dose,scale=2)+rep(rnorm(4,0,4),rep(5,4))
y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
      11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
      20.972267, 17.178012)
resp <- restovec(matrix(y, nrow=4, byrow=TRUE), name="y")
reps <- rmna(resp, tvcov=tvctomat(matrix(dose, nrow=4, byrow=TRUE), name="dose"))

# same linear normal model with random normal intercept fitted four ways
# compare with growth::elliptic(reps, model=~dose, preg=c(0,0.6), pre=4)
glmm(y~dose, nest=individuals, data=reps)
gnlmm(reps, mu=~dose, pmu=c(8.7,0.25), psh=3.5, psd=3)
hnlmix(reps, mu=~a+b*dose+rand, random="rand", pmu=c(8.7,0.25),
	pshape=3.44, prandom=0)

# gamma model with log link and random normal intercept fitted three ways
glmm(y~dose, family=Gamma(link=log), nest=individuals, data=reps, points=8)
gnlmm(reps, distribution="gamma", mu=~exp(a+b*dose), pmu=c(2,0.03),
	psh=1, psd=0.3)
hnlmix(reps, distribution="gamma", mu=~exp(a+b*dose+rand), random="rand",
	pmu=c(2,0.04), pshape=1, prandom=0)

# gamma model with log link and random gamma mixtures
hnlmix(reps, distribution="gamma", mixture="gamma",
	mu=~exp(a*rand+b*dose), random="rand", pmu=c(2,0.04),
	pshape=1.24, prandom=1)
hnlmix(reps, distribution="gamma", mixture="gamma",
	mu=~exp(a+b*dose)*rand, random="rand", pmu=c(2,0.04),
	pshape=1.24, prandom=1)

}
\author{
J.K. Lindsey
}
\keyword{models}
