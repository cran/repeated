name{kalseries}
\title{Repeated Measurements Models for Continuous Variables in the
Presence of Time-varying Covariates and with Frailty or Serial
Dependence}
\usage{
kalseries(response, times=NULL, intensity="exponential",
	depend="independence", mu=NULL, shape=NULL, density=F, ccov=NULL,
	tvcov=NULL, torder=0, interaction=NULL, preg=NULL, ptvc=NULL,
	pintercept=NULL, pshape=1, pinitial=1, pdepend=NULL, delta=NULL,
	transform="identity", link="identity",
	envir=sys.frame(sys.parent()), print.level=0,
	ndigit=10, gradtol=0.00001, steptol=0.00001, iterlim=100,
	fscale=1, typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\alias{kalseries}
\alias{coefficients.kalseries}
\alias{deviance.kalseries}
\alias{fitted.kalseries}
\alias{residuals.kalseries}
\alias{print.kalseries}
\arguments{
\item{response}{A list of two column matrices with responses and
corresponding times for each individual, one matrix or dataframe of
response values, or an object of class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced. Ignored if response has class,
response or repeated.}
\item{intensity}{The form of function to be put in the Pareto distribution.
Choices are exponential, Weibull, gamma, normal, logistic, Cauchy,
Student, log normal, log logistic, log Cauchy, log Student, and
gen(eralized) logistic.}
\item{depend}{Type of dependence. Choices are independence, Markov,
serial, and frailty.}
\item{mu}{A regression function for the location parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. If there are only time-constant
covariates, give the initial estimates in preg; if any covariates are
time-varying, give all initial estimates in ptvc.}
\item{shape}{A regression function for the shape parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. It must yield one value per observation.}
\item{density}{If TRUE, the density of the function specified in
\code{intensity} is used instead of the intensity.}
\item{ccov}{A vector or matrix containing time-constant baseline
covariates with one row per individual, a model formula using
vectors of the same size, or an object of class, tccov (created by
\code{\link{tcctomat}}). Ignored if response has class, repeated.}
\item{tvcov}{A list of matrices with time-varying covariate values,
observed at the event times in \code{response}, for each individual
(one column per variable), one matrix or dataframe of such
covariate values, or an object of class, tvcov (created by
\code{\link{tvctomat}}). If a time-varying covariate is observed at
arbitrary time, \code{\link{gettvc}} can be used to find the most
recent values for each response and create a suitable list. Ignored if
response has class, repeated.}
\item{torder}{The order of the polynomial in time to be fitted.}
\item{interaction}{Vector of length equal to the number of
time-constant covariates, giving the levels of interactions between
them and the polynomial in time in the \code{linear model}.}
\item{preg}{Initial parameter estimates for the regression model:
intercept, one for each covariate in \code{ccov}, and \code{torder} plus
sum(\code{interaction}). If a location function (\code{mu}) is
supplied that contains time-varying covariates, all initial estimates
must be given in ptvc. If \code{mu} is a formula with unknown
parameters, their estimates must be supplied either in their order of
appearance in the expression or in a named list.}
\item{ptvc}{Initial parameter estimates for the coefficients of the
time-varying covariates, as many as in \code{tvcov}. If a location
function (mu) is supplied that contains time-varying covariates, all
initial estimates must be given here.}
\item{pintercept}{The initial estimate of the intercept for the
generalized logistic intensity.}
\item{pshape}{An initial estimate for the shape parameter of the
intensity function (except exponential intensity). If \code{shape} is
a function or formula, the corresponding initial estimates. If
\code{shape} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pinitial}{An initial estimate for the initial parameter. (With
frailty dependence, this is the frailty parameter.)}
\item{pdepend}{An initial estimate for the serial dependence parameter.}
\item{delta}{Scalar or vector giving the unit of measurement for each
response value, set to unity by default. For example, if a response is
measured to two decimals, delta=0.01. If the response has been
pretransformed, this must be multiplied by the Jacobian. This
transformation cannot contain unknown parameters. For example, with a
log transformation, \code{delta=1/y}. The jacobian is calculated
automatically for the transform option. Ignored if response has class,
response or repeated.}
\item{transform}{Transformation of the response variable: \code{identity},
\code{exp}, \code{square}, \code{sqrt}, or \code{log}.}
\item{link}{Link function for the mean: \code{identity}, \code{exp},
\code{square}, \code{sqrt}, or \code{log}.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{kalseries} is designed to handle repeated measurements models with
time-varying covariates. The distributions have two extra parameters
as compared to the functions specified by \code{intensity} and are
generally longer tailed than those distributions. Dependence among
observations on a unit can be through frailty (a type of random effect)
or serial dependence over time.

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)

Marginal and individual profiles can be plotted using
\code{\link{plot.profile}} and \code{\link{plot.iprofile}} and
residuals with \code{\link{plot.residuals}}.

If the responses on a unit are clustered, not longitudinal, use the
failty dependence with the default exponential intensity.
}
\value{
A list of classes \code{kalseries} and \code{recursive} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{finterp}},
\code{\link{gar}}, \code{\link{gettvc}}, \code{\link{gnlmm}},
\code{\link{gnlr}}, \code{\link{kalcount}}, \code{\link{kalsurv}},
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{rmna}},
\code{\link{tcctomat}}, \code{\link{tvctomat}}.
}
\examples{
treat <- c(0,0,1,1)
tr <- tcctomat(treat)
dose <- matrix(rpois(20,10), ncol=5)
dd <- tvctomat(dose)
y <- restovec(matrix(rnorm(20), ncol=5))
reps <- rmna(y, ccov=tr, tvcov=dd)
#
# normal intensity, independence model
kalseries(y, intensity="normal", dep="independence", preg=1, pshape=5)
# random effects
kalseries(y, intensity="normal", dep="frailty", preg=1, pdep=0.1, psh=5)
# serial dependence
kalseries(y, intensity="normal", dep="serial", preg=1, pinitial=1,
	pdep=0.1, psh=5)
#
# add time-constant variable
kalseries(y, intensity="normal", dep="serial", pinitial=1,
	pdep=0.1, psh=5, preg=c(1,0), ccov=treat)
# or equivalently
kalseries(y, intensity="normal", mu=~treat, dep="serial", pinitial=1,
	pdep=0.1, psh=5, preg=c(1,0))
# or
kalseries(y, intensity="normal", mu=~b0+b1*treat, dep="serial",
	pinitial=1, pdep=0.1, psh=5, preg=c(1,0), envir=reps)
#
# add time-varying variable
kalseries(y, intensity="normal", dep="serial", pinitial=1, pdep=0.1,
	psh=5, preg=c(1,0), ccov=treat, ptvc=0, tvc=dose)
# or equivalently, from the environment
kalseries(y, intensity="normal",
	mu=~b0+b1*rep(treat,rep(5,4))+b2*as.vector(t(dose)),
	dep="serial", pinitial=1, pdep=0.1, psh=5, ptvc=c(1,0,0))
# or from the reps data object
kalseries(y, intensity="normal", mu=~b0+b1*treat+b2*dose,
	dep="serial", pinitial=1, pdep=0.1, psh=5,
	ptvc=c(1,0,0), envir=reps)
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
conc <- matrix(rgamma(40,shape(log(c(0.01,1))),mu(log(c(1,0.3,0.2)))),
	ncol=20,byrow=T)
conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
	ncol=20,byrow=T)[,1:19])
conc <- restovec(ifelse(conc>0,conc,0.01))
reps <- rmna(conc, ccov=dd, tvcov=tt)
#
# constant shape parameter
kalseries(reps, intensity="gamma", dep="independence", mu=mu,
	ptvc=c(-1,-1.1,-1), pshape=1.5, envir=reps)
# or
kalseries(reps, intensity="gamma", dep="independence",
	mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
	pshape=1.2, envir=reps)
# add serial dependence
kalseries(reps, intensity="gamma", dep="serial", pdep=0.9,
	mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
	pshape=1.2, envir=reps)
# time dependent shape parameter
kalseries(reps, intensity="gamma", dep="independence", mu=mu,
	shape=shape, ptvc=c(-1,-1.1,-1), pshape=c(-3,0),
	envir=reps)
# or
kalseries(reps, intensity="gamma", dep="independence",
	mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
	shape=~exp(b1-b2)*times*dose*exp(-exp(b1)*times),
	pshape=list(b1=-3,b2=0), envir=reps)
# add serial dependence
kalseries(reps, intensity="gamma", dep="serial", pdep=0.5,
	mu=~exp(absorption-volume)*
	dose/(exp(absorption)-exp(elimination))*
	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
	shape=~exp(b1-b2)*times*dose*exp(-exp(b1)*times),
	pshape=list(b1=-3,b2=0), envir=reps)
}
