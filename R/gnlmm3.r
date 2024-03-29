#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 1998, 1999, 2000, 2001, 2004 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     gnlmm3(y=NULL, distribution="normal", mu=NULL, shape=NULL,
#	nest=NULL, family=NULL, linear=NULL, pmu=NULL, pshape=NULL,
#	pfamily=NULL, psd=NULL, exact=FALSE, wt=1, scale=NULL, points=10,
#	common=FALSE, delta=1, envir=parent.frame(), print.level=0,
#	typsize=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit Generalized nonlinear mixed models with a variety of
# three parameter distributions and a normal random effect.




##' Generalized Nonlinear Mixed Models for Three-parameter Distributions
##' 
##' \code{gnlmm3} fits user-specified nonlinear regression equations to one or
##' more parameters of the common three parameter distributions. The intercept
##' of the location regression has a normally-distributed random effect. This
##' normal mixing distribution is computed by Gauss-Hermite integration.
##' 
##' The \code{scale} of the random effect is the link function to be applied.
##' For example, if it is \code{log}, the supplied mean function, \code{mu}, is
##' transformed as exp(log(mu)+sd), where sd is the random effect parameter.
##' 
##' It is recommended that initial estimates for \code{pmu}, \code{pshape}, and
##' \code{pfamily} be obtained from \code{gnlr3}.
##' 
##' Nonlinear regression models can be supplied as formulae where parameters
##' are unknowns in which case factor variables cannot be used and parameters
##' must be scalars. (See \code{\link[rmutil]{finterp}}.)
##' 
##' The printed output includes the -log likelihood (not the deviance), the
##' corresponding AIC, the maximum likelihood estimates, standard errors, and
##' correlations.
##' 
##' 
##' @param y A response vector for uncensored data, a two column matrix for
##' binomial data or censored data, with the second column being the censoring
##' indicator (1: uncensored, 0: right censored, -1: left censored), or an
##' object of class, \code{response} (created by
##' \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
##' \code{\link[rmutil]{rmna}}) or \code{\link[rmutil]{lvna}}). If the
##' \code{repeated} data object contains more than one response variable, give
##' that object in \code{envir} and give the name of the response variable to
##' be used here.
##' @param distribution Either a character string containing the name of the
##' distribution or a function giving the -log likelihood and calling the
##' location, shape, and family functions. Distributions are Box-Cox
##' transformed normal, generalized inverse Gauss, generalized logistic,
##' Hjorth, generalized gamma, Burr, generalized Weibull, power exponential,
##' Student t, generalized extreme value, power variance function Poisson, and
##' skew Laplace. (For definitions of distributions, see the corresponding
##' [dpqr]distribution help.)
##' @param mu A user-specified function of \code{pmu}, and possibly
##' \code{linear}, giving the regression equation for the location. This may
##' contain a linear part as the second argument to the function. It may also
##' be a formula beginning with ~, specifying a either linear regression
##' function for the location parameter in the Wilkinson and Rogers notation or
##' a general function with named unknown parameters. If it contains unknown
##' parameters, the keyword \code{linear} may be used to specify a linear part.
##' If nothing is supplied, the location is taken to be constant unless the
##' linear argument is given.
##' @param shape A user-specified function of \code{pshape}, and possibly
##' \code{linear} and/or \code{mu}, giving the regression equation for the
##' dispersion or shape parameter. This may contain a linear part as the second
##' argument to the function and the location function as last argument (in
##' which case \code{shfn} must be set to TRUE). It may also be a formula
##' beginning with ~, specifying either a linear regression function for the
##' shape parameter in the Wilkinson and Rogers notation or a general function
##' with named unknown parameters. If it contains unknown parameters, the
##' keyword \code{linear} may be used to specify a linear part and the keyword
##' \code{mu} to specify a function of the location parameter. If nothing is
##' supplied, this parameter is taken to be constant unless the linear argument
##' is given. This parameter is the logarithm of the usual one.
##' @param family A user-specified function of \code{pfamily}, and possibly
##' \code{linear}, for the regression equation of the third (family) parameter
##' of the distribution. This may contain a linear part that is the second
##' argument to the function. It may also be a formula beginning with ~,
##' specifying either a linear regression function for the family parameter in
##' the Wilkinson and Rogers notation or a general function with named unknown
##' parameters. If neither is supplied, this parameter is taken to be constant
##' unless the linear argument is given. In most cases, this parameter is the
##' logarithm of the usual one.
##' @param linear A formula beginning with ~ in W&R notation, specifying the
##' linear part of the regression function for the location parameter or list
##' of two such expressions for the location and/or shape parameters.
##' @param nest The variable classifying observations by the unit upon which
##' they were observed. Ignored if \code{y} or \code{envir} has class,
##' response.
##' @param pmu Vector of initial estimates for the location parameters. If
##' \code{mu} is a formula with unknown parameters, their estimates must be
##' supplied either in their order of appearance in the expression or in a
##' named list.
##' @param pshape Vector of initial estimates for the shape parameters. If
##' \code{shape} is a formula with unknown parameters, their estimates must be
##' supplied either in their order of appearance in the expression or in a
##' named list.
##' @param pfamily Vector of initial estimates for the family parameters. If
##' \code{family} is a formula with unknown parameters, their estimates must be
##' supplied either in their order of appearance in the expression or in a
##' named list.
##' @param psd Initial estimate of the standard deviation of the normal mixing
##' distribution.
##' @param exact If TRUE, fits the exact likelihood function for continuous
##' data by integration over intervals of observation, i.e. interval censoring.
##' @param wt Weight vector.
##' @param delta Scalar or vector giving the unit of measurement (always one
##' for discrete data) for each response value, set to unity by default.
##' Ignored if y has class, response. For example, if a response is measured to
##' two decimals, \code{delta=0.01}. If the response is transformed, this must
##' be multiplied by the Jacobian. The transformation cannot contain unknown
##' parameters. For example, with a log transformation, \code{delta=1/y}. (The
##' delta values for the censored response are ignored.)
##' @param scale The scale on which the random effect is applied:
##' \code{identity}, \code{log}, \code{logit}, \code{reciprocal}, or
##' \code{exp}.
##' @param points The number of points for Gauss-Hermite integration of the
##' random effect.
##' @param common If TRUE, at least two of \code{mu}, \code{shape}, and
##' \code{family} must both be either functions with, as argument, a vector of
##' parameters having some or all elements in common between them so that
##' indexing is in common between them or formulae with unknowns. All parameter
##' estimates must be supplied in \code{pmu}. If FALSE, parameters are distinct
##' between the two functions and indexing starts at one in each function.
##' @param envir Environment in which model formulae are to be interpreted or a
##' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
##' name of the response variable should be given in \code{y}. If \code{y} has
##' class \code{repeated}, it is used as the environment.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of class \code{gnlm} is returned that contains all of the
##' relevant information calculated, including error codes.
##' @author J.K. Lindsey
### @seealso \code{\link[rmutil]{finterp}}, \code{\link[gnlm]{fmr}},
### \code{\link{glm}}, \code{\link[repeated]{gnlmix}},
### \code{\link[repeated]{glmm}}, \code{\link[gnlm]{gnlr}},
### \code{\link[gnlm]{gnlr3}}, \code{\link[repeated]{hnlmix}},
### \code{\link{lm}}, \code{\link[gnlm]{nlr}}, \code{\link[stats]{nls}}.
##' @keywords models
##' @examples
##' 
##' # data objects
##' sex <- c(0,1,1)
##' sx <- tcctomat(sex)
##' #dose <- matrix(rpois(30,10),nrow=3)
##' dose <- matrix(c(8,9,11,9,11,11,7,8,7,12,8,8,9,10,15,10,9,9,20,14,4,7,
##' 	4,13,10,13,6,13,11,17),nrow=3)
##' dd <- tvctomat(dose)
##' # vectors for functions
##' dose <- as.vector(t(dose))
##' sex <- c(rep(0,10),rep(1,20))
##' nest <- rbind(rep(1,10),rep(2,10),rep(3,10))
##' #y <- (rt(30,5)+exp(0.2+0.3*dose+0.5*sex+rep(rnorm(3),rep(10,3))))*3
##' y <- c(62.39712552,196.94419614,2224.74940087,269.56691601,12.86079662,
##' 	14.96743546, 47.45765042,156.51381687,508.68804438,281.11065302,
##' 	92.32443655, 81.88000484, 40.26357733, 13.04433670, 15.58490237,
##' 	63.62154867, 23.69677549, 53.52885894, 88.02507682, 34.04302506,
##' 	44.28232323,116.80732423,106.72564484, 25.09749055, 12.61839145,
##' 	-0.04060996,153.32670123, 63.25866087, 17.79852591,930.52558064)
##' y <- restovec(matrix(y, nrow=3), nest=nest, name="y")
##' reps <- rmna(y, ccov=sx, tvcov=dd)
##' #
##' # log linear regression with Student t distribution
##' mu <- function(p) exp(p[1]+p[2]*sex+p[3]*dose)
##' ## print(z <- gnlm::gnlr3(y, dist="Student", mu=mu, pmu=c(0,0,0), pshape=1, pfamily=1))
##' ## starting values for pmu and pshape from z$coef[1:3] and z$coef[4] respectively
##' ## starting value for pfamily in z$coef[5]
##' gnlmm3(y, dist="Student", mu=mu, nest=nest, pmu=c(3.69,-1.19, 0.039),
##'        pshape=5, pfamily=0, psd=40, points=3)
##' # or equivalently
##' gnlmm3(y, dist="Student", mu=~exp(b0+b1*sex+b2*dose), nest=nest,
##'        pmu=c(3.69,-1.19, 0.039), pshape=5, pfamily=0, psd=40,
##'        points=3, envir=reps)
##' \dontrun{
##' # or with identity link
##' print(z <- gnlm::gnlr3(y, dist="Student", mu=~sex+dose, pmu=c(0.1,0,0), pshape=1,
##' 	pfamily=1))
##' gnlmm3(y, dist="Student", mu=~sex+dose, nest=nest, pmu=z$coef[1:3],
##' 	pshape=z$coef[4], pfamily=z$coef[5], psd=50, points=3)
##' # or
##' gnlmm3(y, dist="Student", mu=~b0+b1*sex+b2*dose, nest=nest, pmu=z$coef[1:3],
##' 	pshape=z$coef[4], pfamily=z$coef[5], psd=50, points=3, envir=reps)
##' #
##' # nonlinear regression with Student t distribution
##' mu <- function(p) p[1]+exp(p[2]+p[3]*sex+p[4]*dose)
##' print(z <- gnlm::gnlr3(y, dist="Student", mu=mu, pmu=c(1,1,0,0), pshape=1,
##' 	pfamily=1))
##' gnlmm3(y, dist="Student", mu=mu, nest=nest, pmu=z$coef[1:4],
##' 	pshape=z$coef[5], pfamily=z$coef[6], psd=50, points=3)
##' # or
##' mu2 <- function(p, linear) p[1]+exp(linear)
##' gnlmm3(y, dist="Student", mu=mu2, linear=~sex+dose, nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5], pfamily=z$coef[6], psd=50,
##' 	points=3)
##' # or
##' gnlmm3(y, dist="Student", mu=~a+exp(linear), linear=~sex+dose, nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5], pfamily=z$coef[6], psd=50,
##' 	points=3)
##' # or
##' gnlmm3(y, dist="Student", mu=~b4+exp(b0+b1*sex+b2*dose), nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5], pfamily=z$coef[6], psd=50,
##' 	points=3, envir=reps)
##' #
##' # include regression for the shape parameter with same mu function
##' shape <- function(p) p[1]+p[2]*sex
##' print(z <- gnlm::gnlr3(y, dist="Student", mu=mu, shape=shape, pmu=z$coef[1:4],
##' 	pshape=c(z$coef[5],0), pfamily=z$coef[6]))
##' gnlmm3(y, dist="Student", mu=mu, shape=shape, nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5:6], pfamily=z$coef[7],
##' 	psd=5, points=3)
##' # or
##' gnlmm3(y, dist="Student", mu=mu, shape=shape, nest=nest, pmu=z$coef[1:4],
##' 	pshape=z$coef[5:6], pfamily=z$coef[7], psd=5, points=3,
##' 	envir=reps)
##' # or
##' gnlmm3(y, dist="Student", mu=~b4+exp(b0+b1*sex+b2*dose), shape=~a1+a2*sex,
##' 	nest=nest, pmu=z$coef[1:4], pshape=z$coef[5:6],
##' 	pfamily=z$coef[7], psd=5, points=3, envir=reps)
##' }
##' @export gnlmm3
gnlmm3 <- function(y=NULL, distribution="normal", mu=NULL, shape=NULL,
	nest=NULL, family=NULL, linear=NULL, pmu=NULL, pshape=NULL,
	pfamily=NULL, psd=NULL, exact=FALSE, wt=1, scale=NULL, points=10,
	common=FALSE, delta=1, envir=parent.frame(), print.level=0,
	typsize=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
#
# Burr cdf
#
pburr <- function(q, m, s, f) 1-(1+(q/m)^s)^-f
#
# generalized logistic cdf
#
pglogis <- function(y, m, s, f) (1+exp(-sqrt(3)*(y-m)/(s*pi)))^-f
#
# generalized Weibull cdf
#
pgweibull <- function(y, s, m, f) (1-exp(-(y/m)^s))^f
#
# Hjorth cdf
#
phjorth <- function(y, m, s, f) 1-(1+s*y)^(-f/s)*exp(-(y/m)^2/2)
#
# generalized inverse Gauss cdf
#
pginvgauss <- function(y, m, s, f)
	.C("pginvgauss_c",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		#DUP=FALSE,
		PACKAGE="repeated")$res
#
# power exponential cdf
#
ppowexp <- function(y, m, s, f){
	z <- .C("ppowexp_c",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		#DUP=FALSE,
		PACKAGE="repeated")$res
	ifelse(y-m>0,0.5+z,0.5-z)}
#
# power variance function Poisson density
#
dpvfpois <- function(y, m, s, f)
	.C("dpvfp_c",
		as.integer(y),
		as.double(m),
		as.double(s),
		as.double(f),
		as.integer(length(y)),
		as.double(rep(1,length(y))),
		res=double(length(y)),
		#DUP=FALSE,
		PACKAGE="repeated")$res
#
# skew Laplace cdf
#
pskewlaplace <- function(q,m,s,f){
u <- (q-m)/s
ifelse(u>0,1-exp(-f*abs(u))/(1+f^2),f^2*exp(-abs(u)/f)/(1+f^2))}
call <- sys.call()
#
# check distribution
#
distribution <- match.arg(distribution,c("normal","inverse Gauss",
	"logistic","Hjorth","gamma","Burr","Weibull","extreme value",
	"Student t","power exponential","power variance function Poisson",
	"skew Laplace"))
#
# check for parameters common to location and shape functions
#
if(common){
	if(sum(is.function(mu)+is.function(shape)+is.function(family))<2&&sum(inherits(mu,"formula")+inherits(shape,"formula")+inherits(family,"formula"))<2)
		stop("with common parameters, at least two of mu, shape, and family must be functions or formulae")
	if((!is.function(mu)&&!inherits(mu,"formula")&&!is.null(mu))||(!is.function(shape)&&!inherits(shape,"formula")&&!is.null(shape))||(!is.function(family)&&!inherits(family,"formula")&&!is.null(family)))
		stop("with common parameters, mu, shape, and family must either be functions, formulae, or NULL")
	if(!is.null(linear))stop("linear cannot be used with common parameters")}
if(!is.null(scale))scale <- match.arg(scale,c("identity","log",
	"reciprocal","exp"))
#
# count number of parameters
#
npl <- length(pmu)
nps <- length(pshape)
npf <- length(pfamily)
if(is.null(psd))stop("An initial value of psd must be supplied")
np <- npl+nps+npf+1
#
# find number of observations now for creating null functions
#
n <- if(inherits(envir,"repeated")||inherits(envir,"response"))sum(nobs(envir))
	else if(inherits(envir,"data.frame"))dim(envir)[1]
	else if(is.vector(y,mode="numeric"))length(y)
	else if(is.matrix(y))dim(y)[1]
	else sum(nobs(y))
if(n==0)stop(paste(deparse(substitute(y)),"not found or of incorrect type"))
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("gnlr3 only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("gnlr3 does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(inherits(envir,"repeated")||inherits(envir,"response"))
		deparse(substitute(envir))
	else NULL
#
# find linear part of each regression and save model for printing
#
lin1 <- lin2 <- lin3 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]
	lin3 <- linear[[3]]}
else lin1 <- linear
if(inherits(lin1,"formula")&&is.null(mu)){
	mu <- lin1
	lin1 <- NULL}
if(inherits(lin2,"formula")&&is.null(shape)){
	shape <- lin2
	lin2 <- NULL}
if(inherits(lin3,"formula")&&is.null(family)){
	family <- lin3
	lin3 <- NULL}
if(inherits(lin1,"formula")){
	lin1model <- if(respenv){
		if(!is.null(attr(finterp(lin1,.envir=y,.name=envname),"parameters")))
			attr(finterp(lin1,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(lin1,.envir=envir,.name=envname),"parameters")))
			attr(finterp(lin1,.envir=envir,.name=envname),"model")}}
else lin1model <- NULL
if(inherits(lin2,"formula")){
	lin2model <- if(respenv){
		if(!is.null(attr(finterp(lin2,.envir=y,.name=envname),"parameters")))
			attr(finterp(lin2,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(lin2,.envir=envir,.name=envname),"parameters")))
			attr(finterp(lin2,.envir=envir,.name=envname),"model")}}
else lin2model <- NULL
if(inherits(lin3,"formula")){
	lin3model <- if(respenv){
		if(!is.null(attr(finterp(lin3,.envir=y,.name=envname),"parameters")))
			attr(finterp(lin3,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(lin3,.envir=envir,.name=envname),"parameters")))
			attr(finterp(lin3,.envir=envir,.name=envname),"model")}}
else lin3model <- NULL
#
# check if linear contains W&R formula
#
if(inherits(lin1,"formula")){
	tmp <- attributes(if(respenv)finterp(lin1,.envir=y,.name=envname)
		else finterp(lin1,.envir=envir,.name=envname))
	lf1 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	if(length(tmp$model)==1){
		if(is.null(mu))mu <- ~1
		else stop("linear must contain covariates")}
	rm(tmp)}
else lf1 <- 0
if(inherits(lin2,"formula")){
	tmp <- attributes(if(respenv)finterp(lin2,.envir=y,.name=envname)
		else finterp(lin2,.envir=envir,.name=envname))
	lf2 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	if(length(tmp$model)==1){
		if(is.null(shape))shape <- ~1
		else stop("linear must contain covariates")}
	rm(tmp)}
else lf2 <- 0
if(inherits(lin3,"formula")){
	tmp <- attributes(if(respenv)finterp(lin3,.envir=y,.name=envname)
		else finterp(lin3,.envir=envir,.name=envname))
	lf3 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	if(length(tmp$model)==1){
		if(is.null(family))family <- ~1
		else stop("linear must contain covariates")}
	rm(tmp)}
else lf3 <- 0
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu2 <- sh2 <- fa2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	# modify formulae
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	if(inherits(shape,"formula")){
		sh2 <- if(respenv)finterp(shape,.envir=y,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
	if(inherits(family,"formula")){
		fa2 <- if(respenv)finterp(family,.envir=y,.name=envname)
			else finterp(family,.envir=envir,.name=envname)}
	# modify functions
	if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
			tmp <- parse(text=deparse(mu)[-1])
			mu <- if(respenv)fnenvir(mu,.envir=y,.name=envname)
				else fnenvir(mu,.envir=envir,.name=envname)
			mu2 <- mu
			attr(mu2,"model") <- tmp}
		else mu2 <- mu}
	if(is.function(shape)){
		if(is.null(attr(shape,"model"))){
		        tmp <- parse(text=deparse(shape)[-1])
		        shape <- if(respenv)fnenvir(shape,.envir=y,.name=envname)
		        	else fnenvir(shape,.envir=envir,.name=envname)
		        sh2 <- shape
		        attr(sh2,"model") <- tmp}
		else sh2 <- shape}
	if(is.function(family)){
		if(is.null(attr(family,"model"))){
		        tmp <- parse(text=deparse(family)[-1])
		        family <- if(respenv)fnenvir(family,.envir=y,.name=envname)
		        	else fnenvir(family,.envir=envir,.name=envname)
		        fa2 <- family
		        attr(fa2,"model") <- tmp}
		else fa2 <- family}}
else {
     if(is.function(mu)&&is.null(attr(mu,"model")))mu <- fnenvir(mu)
     if(is.function(shape)&&is.null(attr(shape,"model")))
		shape <- fnenvir(shape)
     if(is.function(family)&&is.null(attr(family,"model")))
		family <- fnenvir(family)}
#
# transform location formula to function and check number of parameters
#
if(inherits(mu,"formula")){
	if(npl==0)stop("formula for mu cannot be used if no parameters are estimated")
	linarg <- if(lf1>0) "linear" else NULL
	mu3 <- if(respenv)finterp(mu,.envir=y,.name=envname,.args=linarg)
		else finterp(mu,.envir=envir,.name=envname,.args=linarg)
	npt1 <- length(attr(mu3,"parameters"))
	if(is.character(attr(mu3,"model"))){
	# W&R formula
		if(length(attr(mu3,"model"))==1){
		# intercept model
			tmp <- attributes(mu3)
			mu3 <- function(p) p[1]*rep(1,n)
			attributes(mu3) <- tmp}}
	else {
	# formula with unknowns
		if(npl!=npt1&&!common&&lf1==0){
			cat("\nParameters are ")
			cat(attr(mu3,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu3,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}}
else if(!is.function(mu)){
	mu3 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	mu3 <- mu
	npt1 <- length(attr(mu3,"parameters"))-(lf1>0)}
#
# if linear part, modify location function appropriately
#
if(lf1>0){
	if(is.character(attr(mu3,"model")))
		stop("mu cannot be a W&R formula if linear is supplied")
	dm1 <- if(respenv)wr(lin1,data=y)$design
		else wr(lin1,data=envir)$design
	if(is.null(mu2))mu2 <- mu3
	mu1 <- function(p)mu3(p,dm1%*%p[(npt1+1):(npt1+lf1)])}
else {
	if(lf1==0&&length(mu3(pmu))==1){
		mu1 <- function(p) mu3(p)*rep(1,n)
		attributes(mu1) <- attributes(mu3)}
	else {
		mu1 <- mu3
		rm(mu3)}}
#
# give appropriate attributes to mu1 for printing
#
if(is.null(attr(mu1,"parameters"))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=y))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
		else attributes(fnenvir(mu1))}
#
# check that correct number of estimates was supplied
#
nlp <- npt1+lf1
if(!common&&nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
npl <- if(common) 1 else npl+1
npl1 <- if(common&&!inherits(lin2,"formula")) 1 else nlp+2
np1 <- npl+nps
#
# transform shape formula to function and check number of parameters
#
if(inherits(shape,"formula")){
	if(nps==0&&!common)
		stop("formula for shape cannot be used if no parameters are estimated")
	old <- if(common)mu1 else NULL
	linarg <- if(lf2>0) "linear" else NULL
	sh3 <- if(respenv)finterp(shape,.envir=y,.start=npl1,.name=envname,.old=old,.args=linarg)
		else finterp(shape,.envir=envir,.start=npl1,.name=envname,.old=old,.args=linarg)
	npt2 <- length(attr(sh3,"parameters"))
	if(is.character(attr(sh3,"model"))){
	# W&R formula
		if(length(attr(sh3,"model"))==1){
		# intercept model
			tmp <- attributes(sh3)
			sh3 <- function(p) p[npl1]*rep(1,n)
			sh2 <- fnenvir(function(p) p[1]*rep(1,n))
			attributes(sh3) <- tmp}}
	else {
	# formula with unknowns
		if(nps!=npt2&&!common&&lf2==0){
			cat("\nParameters are ")
			cat(attr(sh3,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh3,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}}
else if(!is.function(shape)){
	sh3 <- function(p) p[npl1]*rep(1,n)
	sh2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
else {
	sh3 <- function(p) shape(p[npl1:np])
	attributes(sh3) <- attributes(shape)
	npt2 <- length(attr(sh3,"parameters"))-(lf2>0)}
#
# if linear part, modify shape function appropriately
#
if(lf2>0){
	if(is.character(attr(sh3,"model")))
		stop("shape cannot be a W&R formula if linear is supplied")
	dm2 <- if(respenv)wr(lin2,data=y)$design
		else wr(lin2,data=envir)$design
	if(is.null(sh2))sh2 <- sh3
	sh1 <- sh3(p,dm2%*%p[(npl1+lf2-1):np])}
else {
	sh1 <- sh3
	rm(sh3)}
#
# give appropriate attributes to sh1 for printing
#
if(is.null(attr(sh1,"parameters"))){
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn")){
			if(respenv)attributes(fnenvir(shape,.envir=y))
			else attributes(fnenvir(shape,.envir=envir))}
		else attributes(shape)}
		else attributes(fnenvir(sh1))}
#
# check that correct number of estimates was supplied
#
nlp <- npt2+lf2
if(!common&&nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))
nps1 <- if(common&&!inherits(family,"formula")) 1
	else if(common&&inherits(family,"formula"))
		length(attr(mu1,"parameters"))+nlp+1
	else np1+1
#
# transform family formula to function and check number of parameters
#
if(inherits(family,"formula")){
	if(npf==0&&!common)
		stop("formula for family cannot be used if no parameters are estimated")
	old <- if(common)c(mu1,sh1) else NULL
	linarg <- if(lf3>0) "linear" else NULL
	fa3 <- if(respenv)finterp(family,.envir=y,.start=nps1,.name=envname,.old=old,.args=linarg)
		else finterp(family,.envir=envir,.start=nps1,.name=envname,.old=old,.args=linarg)
	npt3 <- length(attr(fa3,"parameters"))
	if(is.character(attr(fa3,"model"))){
	# W&R formula
		if(length(attr(fa3,"model"))==1){
		# intercept model
			tmp <- attributes(fa3)
			fa3 <- function(p) p[nps1]*rep(1,n)
			fa2 <- fnenvir(function(p) p[1]*rep(1,n))
			attributes(fa3) <- tmp}}
	else {
	# formula with unknowns
		if(npf!=npt3&&!common&&lf3==0){
			cat("\nParameters are ")
			cat(attr(fa3,"parameters"),"\n")
			stop(paste("pfamily should have",npt3,"estimates"))}
		if(is.list(pfamily)){
			if(!is.null(names(pfamily))){
				o <- match(attr(fa3,"parameters"),names(pfamily))
				pfamily <- unlist(pfamily)[o]
				if(sum(!is.na(o))!=length(pfamily))stop("invalid estimates for family - probably wrong names")}
			else pfamily <- unlist(pfamily)}}}
else if(!is.function(family)){
	fa3 <- function(p) p[nps1]*rep(1,n)
	fa2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt3 <- 1}
else {
	fa3 <- function(p) family(p[nps1:np])
	attributes(fa3) <- attributes(family)
	npt3 <- length(attr(fa3,"parameters"))-(lf3>0)}
#
# if linear part, modify family function appropriately
#
if(lf3>0){
	if(is.character(attr(fa3,"model")))
		stop("family cannot be a W&R formula if linear is supplied")
	dm3 <- if(respenv)wr(lin3,data=y)$design
		else wr(lin3,data=envir)$design
	if(is.null(fa2))fa2 <- fa3
	fa1 <- fa3(p,dm3%*%p[(nps1+lf3-1):np])}
else {
	fa1 <- fa3
	rm(fa3)}
#
# give appropriate attributes to fa1 for printing
#
if(is.null(attr(fa1,"parameters"))){
	attributes(fa1) <- if(is.function(family)){
		if(!inherits(family,"formulafn")){
			if(respenv)attributes(fnenvir(family,.envir=y))
			else attributes(fnenvir(family,.envir=envir))}
		else attributes(family)}
		else attributes(fnenvir(fa1))}
#
# check that correct number of estimates was supplied
#
nlp <- npt3+lf3
if(!common&&nlp!=npf)stop(paste("pfamily should have",nlp,"initial estimates"))
#
# when there are parameters common to location, shape, and family functions,
# check that correct number of estimates was supplied
#
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"),attr(fa1,"parameters"))))
	if(nlp!=npl)stop(paste("with a common parameter model, pmu should contain",nlp,"estimates"))}
pmu <- c(pmu,psd)
p <- c(pmu,pshape,pfamily)
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&any(!is.na(y$response$wt)))
		wt <- as.vector(y$response$wt)
	if(!is.null(y$response$delta))
		delta <- as.vector(y$response$delta)
	type <- y$response$type
	respname <- colnames(y$response$y)
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("gnlr3 does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	respname <- colnames(envir$response$y)[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$n)&&!all(is.na(envir$response$n[,col])))
		y <- cbind(y,envir$response$n[,col]-y)
	else if(!is.null(envir$response$censor)&&!all(is.na(envir$response$censor[,col])))
		y <- cbind(y,envir$response$censor[,col])
	if(!is.null(envir$response$wt))wt <- as.vector(envir$response$wt)
	if(!is.null(envir$response$delta))
		delta <- as.vector(envir$response$delta[,col])}
else if(inherits(envir,"data.frame")){
	respname <- deparse(substitute(y))
	y <- envir[[deparse(substitute(y))]]}
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("gnlr3 only handles univariate responses")
	if(!is.null(y$wt)&&any(!is.na(y$wt)))wt <- as.vector(y$wt)
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	type <- y$type
	respname <- colnames(y$y)
	y <- response(y)}
else respname <- deparse(substitute(y))
if(any(is.na(y)))stop("NAs in y - use rmna")
#
# if there is censoring, set up censoring indicators for likelihood
#
censor <- length(dim(y))==2&&dim(y)[2]==2
if(censor&&all(y[,2]==1)){
	y <- y[,1]
	censor <- FALSE}
if(censor){
	y[,2] <- as.integer(y[,2])
	if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
		stop("Censor indicator must be -1s, 0s, and 1s")
	cc <- ifelse(y[,2]==1,1,0)
	rc <- ifelse(y[,2]==0,1,ifelse(y[,2]==-1,-1,0))
	lc <- ifelse(y[,2]==-1,0,1)
	if(any(delta<=0&y[,2]==1))
		stop("All deltas for uncensored data must be positive")
	else {
		delta <- ifelse(delta<=0,0.000001,delta)
		delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001,delta)}}
else {
	if(!is.vector(y,mode="numeric"))stop("y must be a vector")
	if(min(delta)<=0)stop("All deltas for must be positive")}
#
# check that data are appropriate for distribution
#
if(distribution=="power variance function Poisson"){
	if(type!="unknown"&&type!="discrete")
		stop("discrete data required")
	if(censor)stop("censoring not allowed for power variance function Poisson")
	if(any(y<0))stop("All response values must be >= 0")}
else if(distribution!="logistic"&&distribution!="Student t"&&
	distribution!="power exponential"&&distribution!="skew Laplace"){
	if(type!="unknown"&&type!="duration"&&type!="continuous")
		stop("duration data required")
	if((censor&&any(y[,1]<=0))||(!censor&&any(y<=0)))
		stop("All response values must be > 0")}
else if(type!="unknown"&&type!="continuous"&&type!="duration")
	stop("continuous data required")
#
# prepare weights and unit of measurement
#
if(min(wt)<0)stop("All weights must be non-negative")
if(length(wt)==1)wt <- rep(wt,n)
if(length(delta)==1)delta <- rep(delta,n)
#
# check nesting vector
#
if(is.null(nest))stop("A nest vector must be supplied")
else if(length(nest)!=n)stop("nest must be the same length as the other variables")
if(is.factor(nest))nest <- as.numeric(nest)
nind <- length(unique(nest))
#
# prepare Gauss-Hermite integration
#
od <- length(nest)==nind
i <- rep(1:n,points)
ii <- rep(1:nind,points)
k <- NULL
for(j in 1:points)k <- c(k,nest+(j-1)*max(nest))
k <- as.integer(k)
#
# calculate quadrature points
#
quad <- gauss.hermite(points)
sd <- quad[rep(1:points,rep(n,points)),1]
qw <- quad[rep(1:points,rep(nind,points)),2]
#
# set default scale for random effect, if not specified
#
#
if(is.null(scale)){
	if(distribution=="normal"||distribution=="logistic"||
		distribution=="Student t"||distribution=="power exponential"||
		distribution=="skew Laplace")scale <- "identity"
	else scale <- "log"}
#
# modify location function to proper scale for random effect
#
mu4 <- if(scale=="identity") function(p) mu1(p)[i]+p[npl]*sd
	else if(scale=="log") function(p) exp(log(mu1(p))[i]+p[npl]*sd)
	else if(scale=="reciprocal") function(p) 1/(1/mu1(p)[i]+p[npl]*sd)
	else if(scale=="exp") function(p) log(exp(mu1(p))[i]+p[npl]*sd)
#
# check that functions returns appropriate values
#
if(any(is.na(mu1(pmu))))stop("The location regression returns NAs: probably invalid initial values")
if(any(is.na(sh1(p))))stop("The shape regression returns NAs: probably invalid initial values")
if(any(is.na(fa1(p))))stop("The family regression returns NAs: probably invalid initial values")
#
# create the appropriate likelihood function of one observation
#
if (!censor){
	ret <- switch(distribution,
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				y <- y^f/f
				jy <- y^(2*f-1)*delta/(2*f)
				norm <- sign(f)*pnorm(0,m,s)
				-wt*(log((pnorm(y+jy,m,s)-pnorm(y-jy,m,s)))
					-log(1-(f<0)-norm))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				norm <- sign(f)*pnorm(0,m,s)
				-wt*((f-1)*log(y)+log(dnorm(y^f/f,m,s))
					-log(1-(f<0)-norm))}
			const <- -wt*log(delta)}},
	"power exponential"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*log(ppowexp(y+delta/2,m,s)
					-ppowexp(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- 0.5*sh1(p)
				f <- exp(fa1(p))
				b <- 1+1/(2*f)
				wt*(t+(abs(y-mu4(p))/exp(t))^(2*f)/2+
					lgamma(b)+b*log(2))}
			const <- -wt*log(delta)}},
	"inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*log(pginvgauss(y+delta/2,m,s,f)
					-pginvgauss(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(-f*log(m)+(f-1)*log(y)-
					log(2*besselK(1/(s*m),abs(f)))-
					(1/y+y/m^2)/(2*s))}
			const <- -wt*log(delta)}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pglogis(y+delta/2,m,s,f)
					-pglogis(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				m <- (y-mu4(p))/exp(t)*sqrt(3)/pi
				wt*(-fa1(p)+m+t+(exp(fa1(p))+1)*
					log(1+exp(-m)))}
			const <- -wt*(log(delta*sqrt(3)/pi))}},
	Hjorth={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*log(phjorth(y+delta/2,m,s,f)-
					phjorth(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(-f*log(1+s*y)/s-(y/m)^2/2+
					log(y/m^2+f/(1+s*y)))}
			const <- -wt*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				u <- (m/s)^f
				-wt*log(pgamma((y+delta/2)^f,s,scale=u)
					-pgamma((y-delta/2)^f,s,scale=u))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				v <- s*f
				-wt*(v*(t-log(m))-(s*y/m)^f+u+(v-1)*log(y)
					-lgamma(s))}
			const <- -wt*log(delta)}},
	Burr={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pburr(y+delta/2,m,s,f)-
					pburr(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				y1 <- y/m
				-wt*(log(f*s/m)+(s-1)*log(y1)
					-(f+1)*log(1+y1^s))}
			const <- -wt*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pgweibull(y+delta/2,s,m,f)
					-pgweibull(y-delta/2,s,m,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				y1 <- (y/m)^s
				-wt*(t+u+(s-1)*log(y)-s*log(m)+
					(f-1)*log(1-exp(-y1))-y1)}
			const <- -wt*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*log(pt((y+delta/2-m)/s,f)-
					pt((y-delta/2-m)/s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(0.5*sh1(p))
				-wt*log(dt((y-mu4(p))/s,exp(fa1(p)))/s)}
			const <- -wt*(log(delta))}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				y1 <- y^f/f
				ey <- exp(y1)
				jey <- y^(f-1)*ey*delta/2
				norm <- sign(f)*exp(-m^-s)
				-wt*(log((pweibull(ey+jey,s,m)
					-pweibull(ey-jey,s,m))/
					(1-(f>0)+norm)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				f <- fa1(p)
				y1 <- y^f/f
				norm <- sign(f)*exp(-m^-s)
				-wt*(t+s*(y1-log(m))-(exp(y1)/m)^s
					+(f-1)*log(y)-log(1-(f>0)+norm))}
			const <- -wt*log(delta)}},
	"power variance function Poisson"={
			fcn <- function(p) {
				m <- mu4(p)
				-wt*log(dpvfpois(y,m,exp(sh1(p))/m,
					fa1(p)))}
			const <- 0},
        "skew Laplace"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pskewlaplace(y+delta/2,m,s,f)
					-pskewlaplace(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				t <- sh1(p)
				u <- fa1(p)
				f <- exp(u)
				-wt*(u+ifelse(y>m,-f*(y-m),(y-m)/f)/
					s-log(1+f^2)-t)}
			const <- -wt*log(delta)}})}
else {
	# censoring
	ret <- switch(distribution,
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				yy <- y[,1]^f/f
				jy <- y[,1]^(2*f-1)*delta/(2*f)
				norm <- sign(f)*pnorm(0,m,s)
				-wt*(cc*log((pnorm(yy+jy,m,s)-
					pnorm(yy-jy,m,s)))+log(lc-rc*(pnorm(yy,
					m,s)-(f>0)*norm)))/(1-(f<0)-norm)}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(0.5*t)
				f <- fa1(p)
				norm <- sign(f)*pnorm(0,m,s)
				-wt*(cc*(-(t+((y[,1]^f/f-m)/s)^2)/2+(f-1)*
					log(y[,1]))+log(lc-rc
					*(pnorm(y[,1]^f/f,m,s)
					-(f>0)*norm)))/(1-(f<0)-norm)}
			const <- wt*cc*(log(2*pi)/2-log(delta))}},
	"power exponential"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				-wt*(cc*log(ppowexp(y[,1]+delta/2,m,s,f)-
					ppowexp(y[,1]-delta/2,m,s,f))
					+log(lc-rc*ppowexp(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- 0.5*sh1(p)
				s <- exp(t)
				f <- exp(fa1(p))
				b <- 1+1/(2*f)
				-wt*(cc*(-t-(abs(y[,1]-mu4(p))/s)^(2*f)/2-
					lgamma(b)-b*log(2))+log(lc-rc
					*ppowexp(y[,1],m,s,f)))}
			const <- -wt*cc*(log(delta))}},
	"inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				f <- fa1(p)
				-wt*(cc*log(pginvgauss(y[,1]+delta/2,m,s,f)
					-pginvgauss((y[,1]-delta/2)^f/f,m,s,f))+
					log(lc-rc*pginvgauss(y[,1]^f/f,m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*(-f*log(m)+(f-1)*log(y[,1])-
					log(2*besselK(1/(s*m),abs(f)))-
					(1/y[,1]+y[,1]/m^2)/(2*s))+log(lc-rc
					*pginvgauss(y[,1],m,s,f)))}
			const <- -wt*cc*(log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				f <- exp(fa1(p))
				-wt*(cc*log(pglogis(y[,1]+delta/2,m,s,f)-
					pglogis(y[,1]-delta/2,m,s,f))
					+log(lc-rc*pglogis(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				y1 <- (y[,1]-m)/s
				u <- fa1(p)
				f <- exp(u)
				-wt*(cc*(u-y1-log(s)-(f+1)*log(1+exp(-y1)))
					+log(lc-rc*pglogis(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
	Hjorth={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*log(phjorth(y[,1]+delta/2,m,s,f)-
					phjorth(y[,1]-delta/2,m,s,f))
					+log(lc-rc*phjorth(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*(-f*log(1+s*y[,1])/s-(y[,1]/m)^2/2+
					log(y[,1]/m^2+f/(1+s*y[,1])))+
					log(lc-rc*phjorth(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				u <- (m/s)^f
				-wt*(cc*log(pgamma((y[,1]+delta/2)^f,s,
					scale=u)-pgamma((y[,1]-delta/2)^f,s,
					scale=u))+log(lc-rc*pgamma(y[,1]^f,s,
					scale=u)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				v <- s*f
				-wt*(cc*(v*(t-log(m))-(s*y[,1]/m)^f+u+(v-1)
					*log(y[,1])-lgamma(s))+log(lc-rc
					*pgamma(y[,1]^f,s,scale=(m/s)^f)))}
			const <- -wt*cc*log(delta)}},
	Burr={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pburr(y[,1]+delta/2,m,s,f)-
					pburr(y[,1]-delta/2,m,s,f))
					+log(lc-rc*pburr(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				y1 <- y[,1]/m
				-wt*(cc*(log(f*s/m)+(s-1)*log(y1)
					-(f+1)*log(1+y1^s))+
					log(lc-rc*pburr(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pgweibull(y[,1]+delta/2,s,m,f)-
					pgweibull(y[,1]-delta/2,s,m,f))
					+log(lc-rc*pgweibull(y[,1],s,m,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				y1 <- (y[,1]/m)^s
				-wt*(cc*(t+u+(s-1)*log(y[,1])-s*log(m)+
					(f-1)*log(1-exp(-y1))-y1)+log(lc-rc*
					pgweibull(y[,1],s,m,f)))}
			const <- -wt*cc*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pt((y[,1]+delta/2-m)/s,f)-
					pt((y[,1]-delta/2-m)/s,f))
					+log(lc-rc*pt((y[,1]-m)/s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(dt((y[,1]-m)/s,f)/s)
					+log(lc-rc*pt((y[,1]-m)/s,f)))}
			const <- -wt*cc*(log(delta))}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				y1 <- y[,1]^f/f
				ey <- exp(y1)
				jey <- y[,1]^(f-1)*ey*delta/2
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(cc*log(pweibull(ey+jey,s,m)-
					pweibull(ey-jey,s,m))
					+log(lc-rc*(pweibull(ey,s,m)-ind+
					(f>0)*norm))-log(1-ind+norm))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				f <- fa1(p)
				y1 <- y[,1]^f/f
				ey <- exp(y1)
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(cc*(t+s*(y1-log(m))-(ey/m)^s
					+(f-1)*log(y[,1]))+log(lc-rc*
					(pweibull(ey,s,m)-ind+(f>0)*norm))-
					log(1-ind+norm))}
			const <- -wt*cc*log(delta)}},
        "skew Laplace"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pskewlaplace(y[,1]+delta/2,m,s,f)
					-pskewlaplace(y[,1]-delta/2,m,s,f))
					+log(lc-rc*pskewlaplace(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				-wt*(cc*(u+ifelse(y>m,-f*(y-m),(y-m)/f)/
					s-log(1+f^2)-t)+log(lc-rc
					*pskewlaplace(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}})}
#
# put the individual elements of likelihood together, checking for underflow
#
fn <- function(p) {
	under  <- 0
	if(od)pr <- -fcn(p)
	else {
		pr <- NULL
		for(i in split(fcn(p),k))pr <- c(pr,-sum(i))}
	if(any(is.na(pr)))stop("NAs - unable to calculate probabilities.\n Try other initial values.")
	if(max(pr)-min(pr)>1400){
		if(print.level==2)cat("Log probabilities:\n",pr,"\n\n")
		stop("Product of probabilities is too small to calculate.\n Try fewer points.")}
	if(any(pr > 700))under <- 700-max(pr)
	else if(any(pr < -700))under <- -700-min(pr)
	tmp <- NULL
	for(i in split(qw*exp(pr+under),ii))tmp <- c(tmp,sum(i))
	-sum(log(tmp)-under)}
#
# check that the likelihood returns an appropriate value and optimize
#
if(fscale==1)fscale <- fn(p)
if(is.na(fn(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fn,p=p,hessian=TRUE,print.level=print.level,typsize=typsize,
	ndigit=ndigit,gradtol=gradtol,stepmax=stepmax,steptol=steptol,
	iterlim=iterlim,fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
#
# calculate fitted values and raw residuals
#
fitted.values <- as.vector(mu4(z0$estimate))
residuals <- y-fitted.values
#
# calculate se's
#
if(np==0)cov <- NULL
else if(np==1){
	cov <- 1/z0$hessian
	se <- as.vector(sqrt(cov))}
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)
	se <- sqrt(diag(cov))}
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
if(!is.null(sh2))sh1 <- sh2
if(!is.null(fa2))fa1 <- fa2
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	respname=respname,
	mu=mu1,
	shape=sh1,
	family=fa1,
	linear=list(lin1,lin2,lin3),
	linmodel=list(lin1model,lin2model,lin3model),
	common=common,
	scale=scale,
	points=points,
	prior.weights=wt,
	censor=censor,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=z0$minimum+np,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	npl=npl,
	npm=0,
	nps=nps,
	npf=npf,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlm"
return(z1)}
