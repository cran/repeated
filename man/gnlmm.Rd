\name{gnlmm}
\title{Fit a Generalized Nonlinear Mixed Model}
\usage{
gnlmm(y, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
	nest=NULL, pmu=NULL, pshape=NULL, psd=NULL, exact=F, wt=1,
	delta=1, shfn=F, scale=NULL, points=10,
	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p\%*\%p),
	steptol=0.00001, iterlim=100, fscale=1)
}
\alias{gnlmm}
\arguments{
\item{y}{A response vector for uncensored data, a two column matrix
for binomial data or censored data, with the second column being the
censoring indicator (1: uncensored, 0: right censored, -1: left
censored), or an object of class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{distribution}{Either a character string containing the name of the
distribution or a function giving the -log likelihood and calling the
location and shape functions.}
\item{mu}{A user-specified function of \code{pmu}, and possibly
\code{linear}, giving the regression equation for the location. This
may contain a linear part as the second argument to the function. It
may also be a formula beginning with ~, specifying a either linear
regression function for the location parameter in the Wilkinson and
Rogers notation or a general function with named unknown
parameters. If none is supplied, the location is taken to be constant
unless the linear argument is given.}
\item{shape}{A user-specified function of \code{pshape}, and possibly
\code{linear} and/or \code{mu}, giving the regression equation for the
dispersion or shape parameter. This may contain a linear part as the
second argument to the function and the location function as last
argument (in which case \code{shfn} must be set to TRUE). It may also
be a formula beginning with ~, specifying either a linear regression
function for the shape parameter in the Wilkinson and Rogers notation
or a general function with named unknown parameters. If none is
supplied, this parameter is taken to be constant unless the linear
argument is given. This parameter is the logarithm of the usual one.}
\item{linear}{A formula beginning with ~, specifying the linear part
of the regression function for the location parameter or list of two
such expressions for the location and/or shape parameters.}
\item{nest}{The variable classifying observations by the unit upon
which they were observed. Ignored if y has class, response.}
\item{pmu}{Vector of initial estimates for the location parameters.
If \code{mu} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pshape}{Vector of initial estimates for the shape parameters.
If \code{shape} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{psd}{Initial estimate of the standard deviation of the normal
mixing distribution.}
\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}
\item{wt}{Weight vector.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default. Ignored if y has class, response. For example, if a response
is measured to two decimals, delta=0.01. If the response is
transformed, this must be multiplied by the Jacobian. The
transformation cannot contain unknown parameters. For example, with a
log transformation, \code{delta=1/y}. (The delta values for the
censored response are ignored.)}
\item{shfn}{If true, the supplied shape function depends on the
location (function). The name of this location function must be the
last argument of the shape function.}
\item{scale}{The scale on which the random effect is applied:
identity, log, logit, reciprocal, or exp.}
\item{points}{The number of points for Gauss-Hermite integration of
the random effect.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{gnlmm} fits user-specified nonlinear regression equations to one
or both parameters of the common one and two parameter distributions
(binomial, beta binomial, double binomial, mult(iplicative) binomial,
Poisson, negative binomial, double Poisson, mult(iplicative) Poisson,
gamma count, Consul generalized Poisson, logarithmic series,
geometric, normal, inverse Gauss, logistic, exponential, gamma,
Weibull, extreme value, Cauchy, Pareto, Laplace, and Levy; all but the
binomial-based distributions may be right and/or left censored).

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)

The intercept of the location regression has a normally-distributed
random effect. This normal mixing distribution is computed by
Gauss-Hermite integration.

The \code{scale} of the random effect is the link function to be
applied. For example, if it is \code{log}, the supplied mean function,
\code{mu}, is transformed as log(exp(mu)+sd), where sd is the random
effect parameter.

It is recommended that initial estimates for \code{pmu} and
\code{pshape} be obtained from \code{gnlr}.
}
\value{
A list of class gnlr is returned.
The printed output includes the -log likelihood (not the deviance),
the corresponding AIC, the maximum likelihood estimates, standard
errors, and correlations. A list is returned that contains all of the
relevant information calculated, including error codes.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{finterp}}, \code{\link{fmr}}, \code{\link{glm}},
\code{\link{glmm}}, \code{\link{gnlr}}, \code{\link{gnlr3}},
\code{\link{lm}}.
}
\examples{
# data objects
sex <- c(0,1)
sx <- tcctomat(sex)
dose <- matrix(rpois(20,10),nrow=2)
dd <- tvctomat(dose)
# vectors for functions
dose <- as.vector(t(dose))
sex <- c(rep(0,10),rep(1,10))
nest <- rbind(rep(1,10),rep(2,10))
y <- exp(0.2+0.3*dose+0.5*sex)+rgamma(20,2,5)+rep(rnorm(10),rep(2,10))
y <- restovec(matrix(y, nrow=2), nest=nest)
reps <- rmna(y, ccov=sx, tvcov=dd)
#
# log linear regression with gamma distribution
mu <- function(p) exp(p[1]+p[2]*sex+p[3]*dose)
print(z <- gnlr(y, dist="gamma", mu=mu, pmu=c(1,0,0), pshape=1))
gnlmm(y, dist="gamma", mu=mu, nest=nest, pmu=z$coef[1:3],
	pshape=z$coef[4], psd=0.1, points=3)
# or equivalently
gnlmm(y, dist="gamma", mu=~exp(b0+b1*sex+b2*dose), nest=nest,
	pmu=z$coef[1:3], pshape=z$coef[4], psd=0.1, points=3, envir=reps)
# or with identity link
print(z <- gnlr(y, dist="gamma", mu=~sex+dose, pmu=c(1,0,0), pshape=1))
gnlmm(y, dist="gamma", mu=~sex+dose, nest=nest, pmu=z$coef[1:3],
	pshape=z$coef[4], psd=0.1, points=3)
# or
gnlmm(y, dist="gamma", mu=~b0+b1*sex+b2*dose, nest=nest, pmu=z$coef[1:3],
	pshape=z$coef[4], psd=0.1, points=3, envir=reps)
#
# nonlinear regression with gamma distribution
mu <- function(p) p[4]+exp(p[1]+p[2]*sex+p[3]*dose)
print(z <- gnlr(y, dist="gamma", mu=mu, pmu=c(1,0,0,1), pshape=1))
gnlmm(y, dist="gamma", mu=mu, nest=nest, pmu=z$coef[1:4],
	pshape=z$coef[5], psd=0.1, points=3)
# or
gnlmm(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*dose), nest=nest,
	pmu=c(z$coef[4],z$coef[1:3]), pshape=z$coef[5], psd=0.1,
	points=3, envir=reps)
#
# include regression for the shape parameter with same mu function
shape <- function(p) p[1]+p[2]*sex
print(z <- gnlr(y, dist="gamma", mu=mu, shape=shape, pmu=c(1,0,0,1),
	pshape=rep(1,2)))
gnlmm(y, dist="gamma", mu=mu, shape=shape, nest=nest,
	pmu=z$coef[1:4], pshape=z$coef[5:6], psd=0.1, points=3)
# or
gnlmm(y, dist="gamma", mu=mu, shape=shape, nest=nest, pmu=z$coef[1:4],
	pshape=z$coef[5:6], psd=0.1, points=3, envir=reps)
# or
gnlmm(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*dose), shape=shape,
	nest=nest, pmu=c(z$coef[4],z$coef[1:3]),
	pshape=z$coef[5:6], psd=0.1, points=3, envir=reps)
}
