\name{gar}
\title{Generalized autoregression}
\usage{
gar(response, distribution="normal", times=NULL, totals=NULL, censor=NULL,
	delta=NULL, mu=NULL, shape=NULL, shfn=F, common=F, preg=NULL,
	pdepend=NULL, pshape=NULL, transform="identity", link="identity",
	autocorr="exponential", order=1, envir=sys.frame(sys.parent()),
	print.level=0, ndigit=10, gradtol=0.00001,steptol=0.00001,
	fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\alias{gar}
\alias{coefficients.gar}
\alias{deviance.gar}
\alias{fitted.gar}
\alias{residuals.gar}
\alias{print.gar}
\arguments{
\item{response}{A list of two or three column matrices with responses,
corresponding times, and possibly a censor indicator, for each
individual, one matrix or dataframe of responses, or an object of
class, response (created by \code{\link{restovec}}) or repeated
(created by \code{\link{rmna}}).}
\item{distribution}{The distribution to be fitted: Bernoulli, Poisson,
exponential, negative binomial, mult Poisson, double Poisson, beta
binomial, mult binomial, double binomial, normal, inverse Gauss,
logistic, gamma, Weibull, Cauchy, Laplace, Levy, Pareto, gen(eralized)
gamma, gen(eralized) logistic, Hjorth, Burr, gen(eralized) Weibull,
gen(eralized) extreme value, gen(eralized) inverse Gauss, or power
exponential.}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced. Ignored if response has class,
response or repeated.}
\item{totals}{An appropriate scalar, vector, or matrix of binomial
totals (only applicable for binomial, beta binomial, mult binomial,
double binomial). Ignored if response has class, response or repeated.}
\item{censor}{If response is a matrix, a matrix of the same size
containing the censor indicator: 1=uncensored, 0=right-censored,
-1=left-censored. Ignored if response has class, response or repeated.}
\item{delta}{Scalar or vector giving the unit of measurement for each
response value, set to unity by default. For example, if a response is
measured to two decimals, delta=0.01. If the response has been
pretransformed, this must be multiplied by the Jacobian. This
transformation cannot contain unknown parameters. For example, with a
log transformation, \code{delta=1/y}. (The delta values for the
censored response are ignored.) The jacobian is calculated
automatically for the transform option. Ignored if response has class,
response or repeated.}
\item{mu}{A user-specified function of \code{pmu} giving the
regression equation for the location. It may also be a formula
beginning with ~, specifying either a linear regression function for
the location parameter in the Wilkinson and Rogers notation or a
general function with named unknown parameters. It must yield a  value
for each observation on each individual.}
\item{shape}{An optional user-specified shape regression function;
this may depend on the location (function) through its second
argument, in which case, \code{shfn} must be set to TRUE. It may also
be a formula beginning with ~, specifying either a linear regression
function for the shape parameter in the Wilkinson and Rogers notation
or a general function with named unknown parameters.}
\item{shfn}{If TRUE, the supplied shape function depends on the
location function. The name of this location function must be the
last argument of the shape function.}
\item{common}{If TRUE, \code{mu} and \code{shape} must both be
functions with, as argument, a vector of parameters having some or all
elements in common between them so that indexing is in common
between them; all parameter estimates must be supplied in \code{preg}.
If FALSE, parameters are distinct between the two functions and
indexing starts at one in each function.}
\item{preg}{The initial parameter estimates for the location
regression function. If \code{mu} is a formula with unknown
parameters, their estimates must be supplied either in their order of
appearance in the expression or in a named list.}
\item{pdepend}{One or two estimates of the dependence parameters for
the Kalman update. With one, it is Markovian and, with two, it is
nonstationary. For the latter, the \code{order} must be one.}
\item{pshape}{Zero to two estimates for the shape parameters of the
distribution if \code{shape} is not a function; otherwise, estimates
for the parameters in this function, with one extra at the end for
three-parameter distributions. If \code{shape} is a formula with
unknown parameters, their estimates must be supplied either in their
order of appearance in the expression or in a named list.}
\item{transform}{Transformation of the response variable: \code{identity},
\code{exp}, \code{square}, \code{sqrt}, or \code{log}.}
\item{link}{Link function for the mean: \code{identity}, \code{exp},
\code{square}, \code{sqrt}, \code{log}, \code{logit}, or
\code{cloglog} (last two only for binary data).}
\item{autocorr}{The form of the (second if two) dependence function:
\code{exponential} is the usual rho^|t_i-t_j|; \code{gaussian} is
rho^((t_i-t_j)^2); \code{cauchy} is 1/(1+rho(t_i-t_j)^2);
\code{spherical} is ((|t_i-t_j|rho)^3-3|t_i-t_j|rho+2)/2 for
|t_i-t_j|<=1/rho and zero otherwise; \code{IOU} is the integrated
Ornstein-Uhlenbeck process, (2rho min(t_i,t_j)+exp(-rho t_i)
+exp(-rho t_j)-1 -exp(rho|ti-t_j|))/2rho^3.}
\item{order}{First- or second-order stationary autoregression.}
\item{envir}{Environment in which model formulae are to be
interpreted or a data object of class, repeated, tccov, or tvcov.
If \code{response} has class \code{repeated}, it is used as the
environment.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{gar} fits a first- or second-order generalized autoregression,
possibly with Kalman update over time (first-order only).

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)

Marginal and individual profiles can be plotted using
\code{\link{profile}} and \code{\link{iprofile}} and
residuals with \code{\link{plot.residuals}}.
}
\references{
Lindsey, J.K. (1997) Applying Generalized Linear Models. Springer,
pp.\ 93--101

Lambert, P. (1996) Statistics in Medicine 15, 1695-1708
}
\value{
A list of classes \code{gar} and \code{recursive} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{finterp}},
\code{\link{gnlmm}}, \code{\link{gnlr}}, \code{\link{iprofile}},
\code{\link{kalcount}}, \code{\link{kalseries}},
\code{\link{kalsurv}}, \code{\link{plot.residuals}},
\code{\link{profile}}, \code{\link{read.list}},
\code{\link{restovec}}, \code{\link{rmna}}, \code{\link{tcctomat}},
\code{\link{tvctomat}}.
}
\examples{
# first-order one-compartment model
# data objects for formulae
dose <- c(2,5)
dd <- tcctomat(dose)
times <- matrix(rep(1:20,2), nrow=2, byrow=T)
tt <- tvctomat(times)
# vector covariates for functions
dose <- c(rep(2,20),rep(5,20))
times <- rep(1:20,2)
# functions
mu <- function(p) exp(p[1]-p[3])*(dose/(exp(p[1])-exp(p[2]))*
	(exp(-exp(p[2])*times)-exp(-exp(p[1])*times)))
shape <- function(p) exp(p[1]-p[2])*times*dose*exp(-exp(p[1])*times)
# response
conc <- matrix(rgamma(40,shape(log(c(0.1,0.4))),mu(log(c(1,0.3,0.2)))),
	ncol=20,byrow=T)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=T)[,1:19])
conc <- restovec(ifelse(conc>0,conc,0.01))
reps <- rmna(conc, ccov=dd, tvcov=tt)
# constant shape parameter
gar(conc, dist="gamma", times=1:20, mu=mu,
	preg=log(c(1,0.4,0.1)), pdepend=0.5, pshape=1)
# or
gar(conc, dist="gamma", times=1:20, mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	preg=list(absorption=0,elimination=log(0.4),volume=log(0.1)),
	pdepend=0.5, pshape=1, envir=reps)
# (if the covariates contained NAs, reps would have to be used as
# response instead of conc)
#
# time dependent shape parameter
gar(conc, dist="gamma", times=1:20, mu=mu, shape=shape,
	preg=log(c(1,0.4,0.1)), pdepend=0.5, pshape=log(c(1,0.2)))
# or
gar(conc, dist="gamma", times=1:20, mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	shape=~exp(b1-b2)*times*dose*exp(-exp(b1)*times),
	preg=list(absorption=0,elimination=log(0.4),volume=log(0.1)),
	pdepend=0.5, pshape=list(b1=0,b2=log(0.2)), envir=reps)
}
