\name{kalcount}
\title{Repeated Measurements Models for Counts in the Presence of
Time-varying Covariates and with Frailty or Serial Dependence}
\usage{
kalcount(response, times=NULL, origin=0, intensity="exponential",
	depend="independence", update="Markov", mu=NULL, shape=NULL,
	density=F, ccov=NULL, tvcov=NULL, preg=NULL, ptvc=NULL,
	pbirth=NULL, pintercept=NULL, pshape=1, pinitial=1, pdepend=NULL,
	envir=sys.frame(sys.parent()), print.level=0, ndigit=10,
	gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
	typsiz=abs(p), stepmax=10*sqrt(p\%*\%p))
}
\arguments{
\item{response}{A list of two column matrices with counts and
corresponding times for each individual, one matrix or dataframe of
counts, or an object of class, \code{response} (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).
The time origin is taken to be zero and the given times to be the ends
of periods (since the previous time given) in which the counts occurred.}
\item{times}{When response is a matrix, a vector of possibly unequally
spaced times when they are the same for all individuals or a matrix of
times. Not necessary if equally spaced. Ignored if response has class,
\code{response} or \code{repeated}.}
\item{origin}{If the time origin is to be before the start of
observations, a positive constant to be added to all times.}
\item{intensity}{The form of function to be put in the Pareto distribution.
Choices are exponential, Weibull, gamma, log normal, log logistic,
log Cauchy, log Student, and gen(eralized) logistic.}
\item{depend}{Type of dependence. Choices are independence, frailty,
and serial.}
\item{update}{Type of for serial dependence. Choices are Markov, serial,
event, cumulated, count, and kalman. With frailty dependence,
weighting by length of observation time may be specified by setting
update to \code{time}.}
\item{mu}{A regression function for the location parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation (a log link is assumed)
or a general function with named unknown parameters. If there are only
time-constant covariates, give the initial estimates in preg; if any
covariates are time-varying, give all initial estimates in ptvc.}
\item{shape}{A regression function for the shape parameter or a
formula beginning with ~, specifying either a linear regression
function in the Wilkinson and Rogers notation or a general function
with named unknown parameters. It must yield one value per observation.}
\item{density}{If TRUE, the density of the function specified in
\code{intensity} is used instead of the intensity.}
\item{ccov}{A vector or matrix containing time-constant baseline
covariates with one row per individual, a model formula using
vectors of the same size, or an object of class, \code{tccov} (created
by \code{\link{tcctomat}}). If response has class, \code{repeated},
the covariates must be supplied as a Wilkinson and Rogers formula
unless none are to be used or \code{mu} is given.}
\item{tvcov}{A list of matrices with time-varying covariate values,
observed in the time periods in \code{response}, for each individual
(one column per variable), one matrix or dataframe of such covariate
values, or an object of class, \code{tvcov} (created by
\code{\link{tvctomat}}). If response has class, \code{repeated}, the
covariates must be supplied as a Wilkinson and Rogers formula unless
none are to be used or \code{mu} is given.}
\item{preg}{Initial parameter estimates for the regression model:
intercept plus one for each covariate in \code{ccov}. If a location
function (mu) is supplied that contains time-varying covariates, all
initial estimates must be given in ptvc. If \code{mu} is a formula
with unknown parameters, their estimates must be supplied either in
their order of appearance in the expression or in a named list.}
\item{ptvc}{Initial parameter estimates for the coefficients of the
time-varying covariates, as many as in \code{tvcov}. If a location
function (mu) is supplied that contains time-varying covariates, all
initial estimates must be given here.}
\item{pbirth}{If supplied, this is the initial estimate for the
coefficient of the birth model.}
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
\item{pdepend}{An initial estimate for the serial dependence
parameter.}
\item{envir}{Environment in which model formulae are to be
interpreted or a data object of class, \code{repeated}, \code{tccov},
or \code{tvcov}.
If \code{response} has class \code{repeated}, it is used as the
environment.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{kalcount} is designed to handle repeated measurements models
with time-varying covariates. The distributions have two extra
parameters as compared to the functions specified by \code{intensity}
and are generally longer tailed than those distributions. Dependence
among observations on a unit can be through frailty (a type of random
effect) or serial dependence over time.

Here, the variance, with exponential intensity, is a quadratic
function of the mean, whereas, for \code{\link{nbkal}}, it is
proportional to the mean function.

If the counts on a unit are clustered, not longitudinal, use the
failty dependence with the default exponential intensity.

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)

Marginal and individual profiles can be plotted using
\code{\link{profile}} and \code{\link{iprofile}} and
residuals with \code{\link{plot.residuals}}.
}
\value{
A list of classes \code{kalcount} and \code{recursive} is returned.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{carma}}, \code{\link{elliptic}}, \code{\link{finterp}},
\code{\link{gar}}, \code{\link{gettvc}}, \code{\link{gnlmm}},
\code{\link{gnlr}}, \code{\link{iprofile}}, \code{\link{kalseries}},
\code{\link{kalsurv}}, \code{\link{nbkal}}, \code{\link{profile}},
\code{\link{read.list}}, \code{\link{restovec}}, \code{\link{rmna}},
\code{\link{tcctomat}}, \code{\link{tvctomat}}.
}
\examples{
treat <- c(0,0,1,1)
tr <- tcctomat(treat)
dose <- # matrix(rpois(20,10),ncol=5)
	matrix(c(9,13,16,7,12,6,9,10,11,9,10,10,7,9,9,9,8,10,15,4),
		ncol=5,byrow=T)
dd <- tvctomat(dose)
y <- # matrix(rpois(20,1+3*rep(treat,5)),ncol=5)
	restovec(matrix(c(1,1,1,1,0,1,0,1,0,5,3,3,4,1,4,4,2,3,2,5),
		ncol=5,byrow=T))
reps <- rmna(y, ccov=tr, tvcov=dd)
#
# log normal intensity, independence model
kalcount(y, intensity="log normal", dep="independence", preg=1,
	pshape=0.1)
# random effects
kalcount(y, intensity="log normal", dep="frailty", pdep=0.1, preg=1,
	psh=0.1)
# serial dependence
kalcount(y, intensity="log normal", dep="serial", pinitial=0.1,
	preg=1, pdep=0.01, psh=0.1)
# add time-constant variable
kalcount(y, intensity="log normal", pinitial=0.1, psh=0.1,
	preg=c(1,0), ccov=treat)
# or equivalently
kalcount(y, intensity="log normal", mu=~treat, pinitial=0.1,
	psh=0.1, preg=c(1,0))
# or
kalcount(y, intensity="log normal", mu=~b0+b1*treat,
	pinitial=0.1, psh=0.1, preg=c(1,0), envir=reps)
# add time-varying variable
kalcount(y, intensity="log normal", pinitial=0.1, psh=0.1,
	preg=c(1,0), ccov=treat, ptvc=0, tvc=dose)
# or equivalently, from the environment
kalcount(y, intensity="log normal",
	mu=~b0+b1*rep(treat,rep(5,4))+b2*as.vector(t(dose)),
	pinitial=0.1, psh=0.1, ptvc=c(1,0,0))
# or from the reps data object
kalcount(y, intensity="log normal", mu=~b0+b1*treat+b2*dose,
	pinitial=0.1, psh=0.1, ptvc=c(1,0,0), envir=reps)
}
