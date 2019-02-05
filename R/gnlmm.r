#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
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
#     gnlmm(y=NULL, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
#	nest=NULL, pmu=NULL, pshape=NULL, psd=NULL, exact=FALSE, wt=1,
#	delta=1, shfn=FALSE, scale=NULL, points=10, common=FALSE,
#	envir=parent.frame(), print.level=0, typsize=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=sqrt(p%*%p)/10,
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit generalized nonlinear mixed models with normal
#  random effect



##' Generalized Nonlinear Mixed Models
##' 
##' \code{gnlmm} fits user-specified nonlinear regression equations to one or
##' both parameters of the common one and two parameter distributions. The
##' intercept of the location regression has a normally-distributed random
##' effect. This normal mixing distribution is computed by Gauss-Hermite
##' integration.
##' 
##' The \code{scale} of the random effect is the link function to be applied.
##' For example, if it is \code{log}, the supplied mean function, \code{mu}, is
##' transformed as exp(log(mu)+sd), where sd is the random effect parameter.
##' 
##' It is recommended that initial estimates for \code{pmu} and \code{pshape}
##' be obtained from \code{gnlr}.
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
##' be used here. The beta, simplex, and two-sided power distributions for
##' proportions do not allow censoring.
##' @param distribution Either a character string containing the name of the
##' distribution or a function giving the -log likelihood and calling the
##' location and shape functions. Distributions are binomial, beta binomial,
##' double binomial, mult(iplicative) binomial, Poisson, negative binomial,
##' double Poisson, mult(iplicative) Poisson, gamma count, Consul generalized
##' Poisson, logarithmic series, geometric, normal, inverse Gauss, logistic,
##' exponential, gamma, Weibull, extreme value, Cauchy, Pareto, Laplace, and
##' Levy, beta, simplex, and two-sided power. All but the binomial-based
##' distributions and the beta, simplex, and two-sided power may be right
##' and/or left censored. (For definitions of distributions, see the
##' corresponding [dpqr]distribution help.)
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
##' @param shfn If true, the supplied shape function depends on the location
##' (function). The name of this location function must be the last argument of
##' the shape function.
##' @param scale The scale on which the random effect is applied:
##' \code{identity}, \code{log}, \code{logit}, \code{reciprocal}, or
##' \code{exp}.
##' @param points The number of points for Gauss-Hermite integration of the
##' random effect.
##' @param common If TRUE, \code{mu} and \code{shape} must both be either
##' functions with, as argument, a vector of parameters having some or all
##' elements in common between them so that indexing is in common between them
##' or formulae with unknowns. All parameter estimates must be supplied in
##' \code{pmu}. If FALSE, parameters are distinct between the two functions and
##' indexing starts at one in each function.
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
##' dose <- matrix(rpois(30,10),nrow=3)
##' dd <- tvctomat(dose)
##' # vectors for functions
##' dose <- as.vector(t(dose))
##' sex <- c(rep(0,10),rep(1,20))
##' nest <- rbind(rep(1,10),rep(2,10),rep(3,10))
##' #y <- rgamma(30,2,scale=exp(0.2+0.1*dose+0.1*sex+rep(rnorm(3),rep(10,3)))/2)
##' y <- c(0.6490851,0.9313931,0.4765569,0.4188045,2.8339637,2.8158090,
##' 	2.6059975,2.9958184,2.7351583,3.2884980,1.1180961,0.9443986,1.7915571,
##' 	9.0013379,2.3969570,3.4227356,0.5045518,0.7452521,1.8712467,3.6814198,
##' 	0.1489849,1.0327552,0.6102406,1.1536620,2.9145237,9.2847798,5.6454605,
##' 	1.9759672,1.5798008,5.1024496)
##' y <- restovec(matrix(y, nrow=3), nest=nest, name="y")
##' reps <- rmna(y, ccov=sx, tvcov=dd)
##' #
##' # log linear regression with gamma distribution
##' mu <- function(p) exp(p[1]+p[2]*sex+p[3]*dose)
##' ## print(z <- gnlm::gnlr(y, dist="gamma", mu=mu, pmu=c(1,0,0), pshape=1))
##' ## starting values for pmu and pshape from z$coef[1:3] and z$coef[4] respectively
##' gnlmm(y, dist="gamma", mu=mu, nest=nest, pmu=c(1.101,0.326,-0.045),
##'      pshape=0.391, psd=0.1, points=3)
##' # or equivalently
##' gnlmm(y, dist="gamma", mu=~exp(b0+b1*sex+b2*dose), nest=nest,
##'      pmu=c(1.101,0.326,-0.045), pshape=0.391, psd=0.1, points=3, envir=reps)
##' \dontrun{
##' # or with identity link
##' print(z <- gnlm::gnlr(y, dist="gamma", mu=~sex+dose, pmu=c(0.1,0,0), pshape=1))
##' gnlmm(y, dist="gamma", mu=~sex+dose, nest=nest, pmu=z$coef[1:3],
##' 	pshape=z$coef[4], psd=0.1, points=3)
##' # or
##' gnlmm(y, dist="gamma", mu=~b0+b1*sex+b2*dose, nest=nest, pmu=z$coef[1:3],
##' 	pshape=z$coef[4], psd=0.1, points=3, envir=reps)
##' #
##' # nonlinear regression with gamma distribution
##' mu <- function(p) p[1]+exp(p[2]+p[3]*sex+p[4]*dose)
##' print(z <- gnlm::gnlr(y, dist="gamma", mu=mu, pmu=c(1,1,0,0), pshape=1))
##' gnlmm(y, dist="gamma", mu=mu, nest=nest, pmu=z$coef[1:4],
##' 	pshape=z$coef[5], psd=0.1, points=3)
##' # or
##' mu2 <- function(p, linear) p[1]+exp(linear)
##' gnlmm(y, dist="gamma", mu=mu2, linear=~sex+dose, nest=nest,
##' 	pmu=z$coef[1:4], pshape=1, psd=0.1, points=3)
##' # or
##' gnlmm(y, dist="gamma", mu=~a+exp(linear), linear=~sex+dose, nest=nest,
##' 	pmu=z$coef[1:4], pshape=1, psd=0.1, points=3)
##' # or
##' gnlmm(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*dose), nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5], psd=0.1,
##' 	points=3, envir=reps)
##' #
##' # include regression for the shape parameter with same mu function
##' shape <- function(p) p[1]+p[2]*sex
##' print(z <- gnlm::gnlr(y, dist="gamma", mu=mu, shape=shape, pmu=z$coef[1:4],
##' 	pshape=rep(1,2)))
##' gnlmm(y, dist="gamma", mu=mu, shape=shape, nest=nest,
##' 	pmu=z$coef[1:4], pshape=z$coef[5:6], psd=0.1, points=3)
##' # or
##' gnlmm(y, dist="gamma", mu=mu, shape=shape, nest=nest, pmu=z$coef[1:4],
##' 	pshape=z$coef[5:6], psd=0.1, points=3, envir=reps)
##' # or
##' gnlmm(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*dose), shape=~a1+a2*sex,
##' 	nest=nest, pmu=z$coef[1:4], pshape=z$coef[5:6], psd=0.1,
##' 	points=3, envir=reps)
##' }
##' @export gnlmm
gnlmm <- function(y=NULL, distribution="normal", mu=NULL, shape=NULL,
	linear=NULL, nest=NULL, pmu=NULL, pshape=NULL, psd=NULL, exact=FALSE,
	wt=1, delta=1, shfn=FALSE, scale=NULL, points=10, common=FALSE,
	envir=parent.frame(), print.level=0, typsize=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=sqrt(p%*%p)/10, steptol=0.00001,
	iterlim=100, fscale=1){
#
# inverse Gaussian cdf
#
pinvgauss <- function(y,m,s){
	t <- y/m
	v <- sqrt(y*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
#
# Laplace cdf
#
plaplace <- function(y,m,s){
	u <- (y-m)/s
	t <- exp(-abs(u))/2
	ifelse(u<0,t,1-t)}
#
# Levy cdf (a pure R version is imported from rmutil)
#
# plevy <- function(y, m, s)
# 	.C("plevy_c",
# 		as.double(y),
# 		as.double(m),
# 		as.double(s),
# 		as.double(1),
# 		len=as.integer(n),
# 		eps=as.double(1.0e-6),
# 		pts=as.integer(5),
# 		max=as.integer(16),
# 		err=integer(1),
# 		res=double(n),
# 		#DUP=FALSE,
# 		PACKAGE="repeated")$res
#
# simplex cdf
#
psimplex <- function(y, m, s)
        z <- .C("psimplex_c",
        	as.double(y),
        	as.double(m),
        	as.double(s),
        	as.double(1),
        	len=as.integer(n),
        	eps=as.double(1.0e-6),
        	pts=as.integer(5),
        	max=as.integer(16),
        	err=integer(1),
        	res=double(n),
        	#DUP=FALSE,
		PACKAGE="repeated")$res

call <- sys.call()
#
# check distribution
#
distribution <- match.arg(distribution,c("binomial","beta binomial",
	"double binomial","mult binomial","Poisson","negative binomial",
	"double Poisson","mult Poisson","gamma count","Consul","logarithmic",
	"geometric","normal","inverse Gauss","logistic","exponential","gamma",
	"Weibull","extreme value","Pareto","Cauchy","Laplace","Levy","beta",
	"simplex","two-sided power"))
shp <- distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"
if(!shp){
	pshape <- NULL
	shfn <- common <- FALSE}
#
# check for parameters common to location and shape functions
#
if(common){
	if(!is.function(mu)&&!inherits(mu,"formula"))
		stop("with common parameters, mu must be a function or formula")
	if(!is.function(shape)&&!inherits(shape,"formula"))
		stop("with common parameters, shape must be a function or formula")
	if(!is.null(linear))stop("linear cannot be used with common parameters")}
if(!is.null(scale))scale <- match.arg(scale,c("identity","log","logit",
	"reciprocal","exp"))
#
# count number of parameters
#
npl <- length(pmu)
nps <- length(pshape)
if(is.null(psd))stop("An initial value of psd must be supplied")
np <- npl+nps+1
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
		stop("gnlmm only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("gnlmm does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(inherits(envir,"repeated")||inherits(envir,"response"))
		deparse(substitute(envir))
	else NULL
#
# find linear part of each regression and save model for printing
#
lin1 <- lin2 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]}
else lin1 <- linear
if(inherits(lin1,"formula")&&is.null(mu)){
	mu <- lin1
	lin1 <- NULL}
if(inherits(lin2,"formula")&&is.null(shape)){
	shape <- lin2
	lin2 <- NULL}
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
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu2 <- sh2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	# modify formulae
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	if(inherits(shape,"formula")){
		sh2 <- if(respenv)finterp(shape,.envir=y,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
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
		else sh2 <- shape}}
else {
     if(is.function(mu)&&is.null(attr(mu,"model")))mu <- fnenvir(mu)
     if(is.function(shape)&&is.null(attr(shape,"model")))
		shape <- fnenvir(shape)}
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
#
# transform shape formula to function and check number of parameters
#
if(inherits(shape,"formula")){
	if(nps==0&&!common)
		stop("formula for shape cannot be used if no parameters are estimated")
	old <- if(common)mu1 else NULL
	mufn <- if(shfn)"mu" else NULL
	mufn <- c(mufn,if(lf2>0) "linear" else NULL)
	sh4 <- if(respenv)finterp(shape,.envir=y,.start=npl1,.name=envname,.old=old,.args=mufn)
		else finterp(shape,.envir=envir,.start=npl1,.name=envname,.old=old,.args=mufn)
	tmp <- attributes(sh4)
	sh3 <- if(shfn)function(p) sh4(p,mu1(p)) else sh4
	attributes(sh3) <- tmp
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
		if(nps!=npt2&&!common){
			cat("\nParameters are ")
			cat(attr(sh3,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh3,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}}
else if(!is.function(shape)&&shp){
	sh3 <- function(p) p[npl1]*rep(1,n)
	sh2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
else if(shp){
	sh3 <- if(shfn)function(p) shape(p[npl1:np], mu1(p))
		else function(p) shape(p[npl1:np])
	attributes(sh3) <- attributes(shape)
	npt2 <- length(attr(sh3,"parameters"))-(lf2>0)-shfn}
else sh3 <- NULL
#
# if linear part, modify shape function appropriately
#
if(lf2>0){
	if(is.character(attr(sh3,"model")))
		stop("shape cannot be a W&R formula if linear is supplied")
	dm2 <- if(respenv)wr(lin2,data=y)$design
		else wr(lin2,data=envir)$design
	if(is.null(sh2))sh2 <- sh3
	sh1 <- if(shfn)function(p)sh3(p,dm2%*%p[(npl1+lf2-1):np],mu1(p))
		else sh3(p,dm2%*%p[(npl1+lf2-1):np])}
else {
	sh1 <- sh3
	rm(sh3)}
#
# if distribution has a shape parameter, give appropriate attributes to
# sh1 for printing and check that correct number of estimates was supplied
#
if(shp){
	if(is.null(attr(sh1,"parameters"))){
		attributes(sh1) <- if(is.function(shape)){
			if(!inherits(shape,"formulafn")){
				if(respenv)attributes(fnenvir(shape,.envir=y))
				else attributes(fnenvir(shape,.envir=envir))}
			else attributes(shape)}
			else attributes(fnenvir(sh1))}
	nlp <- npt2+lf2
	if(!common&&nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))}
#
# when there are parameters common to location and shape functions,
# check that correct number of estimates was supplied
#
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"))))-shfn
	if(nlp!=npl)stop(paste("with a common parameter model, pmu should contain",nlp,"estimates"))}
pmu <- c(pmu,psd)
p <- c(pmu,pshape)
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
	if(length(nobs(y))==1&&nobs(y)==1)
		nest <- 1:length(y$response$y)
	else nest <- covind(y)
	type <- y$response$type
	respname <- colnames(y$response$y)
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("gnlmm does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	respname <- colnames(envir$response$y)[col]
	y <- envir$response$y[,col]
	if(length(nobs(envir))==1&&nobs(envir)==1)nest <- 1:length(y)
	else nest <- covind(envir)
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
	if(dim(y$y)[2]>1)stop("gnlmm only handles univariate responses")
	if(!is.null(y$wt)&&any(!is.na(y$wt)))wt <- as.vector(y$wt)
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	nest <- covind(y)
	type <- y$type
	respname <- colnames(y$y)
	y <- response(y)}
else respname <- deparse(substitute(y))
if(any(is.na(y)))stop("NAs in y - use rmna")
#
# check that data are appropriate for distribution
#
if(distribution=="binomial"||distribution=="double binomial"||
	distribution=="beta binomial"||distribution=="mult binomial"){
	# binomial data
	if(type!="unknown"&&type!="nominal")stop("nominal data required")
	if(distribution=="binomial"&&(is.vector(y)||(length(dim(y))==2&&
		dim(y)[2]==1))&&all(y==0|y==1))y <- cbind(y,1-y)
	if(length(dim(y))!=2||dim(y)[2]!=2)
		stop(paste("Two column matrix required for response: successes and failures"))
	if(any(y<0))stop("All response values must be positive")
	nn <- y[,1]+y[,2]
	censor <- FALSE}
else {
	# censoring present?
	censor <- length(dim(y))==2&&dim(y)[2]==2
	if(censor&&all(y[,2]==1)){
		y <- y[,1]
		censor <- FALSE}
	if(!censor)if(!is.vector(y,mode="numeric"))stop("y must be a vector")
	if(censor&&(distribution=="beta"||distribution=="simplex"||
		distribution=="two-sided power"||distribution=="gamma count"||
		distribution=="gamma count"||distribution=="logarithmic"))
		stop("Censoring not allowed for this distribution")
	if(distribution=="double Poisson"||distribution=="mult Poisson")
		my <- if(censor)3*max(y[,1]) else 3*max(y)}
if(distribution=="inverse Gauss"||distribution=="exponential"||
	distribution=="gamma"||distribution=="Weibull"||
	distribution=="extreme value"){
	if(type!="unknown"&&type!="duration"&&type!="continuous")
		stop("duration data required")
	if((censor&&any(y[,1]<=0))||(!censor&&any(y<=0)))
		stop("All response values must be > 0")}
else if(distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="gamma count"||distribution=="double Poisson"||
	distribution=="mult Poisson"){
	if(type!="unknown"&&type!="discrete")stop("discrete data required")
	if(any(y<0))stop("All response values must be >= 0")}
else if(distribution=="logarithmic"){
	if(type!="unknown"&&type!="discrete")stop("discrete data required")
	if(any(y[wt>0]<1))stop("All response values must be integers > 0")}
else if(distribution=="beta"||distribution=="simplex"||
	distribution=="two-sided power"){
	if(type!="unknown"&&type!="continuous")stop("continuous data required")
	if(any(y<=0)||any(y>=1))
		stop("All response values must lie between 0 and 1")
	if(distribution=="two-sided power"){
		y1 <- y+delta/2
		y2 <- y-delta/2}}
else if(distribution!="binomial"&&distribution!="double binomial"&&
	distribution!="beta binomial"&&distribution!="mult binomial"&&
	type!="unknown"&&type!="continuous"&&type!="duration")
	stop("continuous data required")
#
# if there is censoring, set up censoring indicators for likelihood
#
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
	if(min(delta)<=0)stop("All deltas for must be positive")}
#
# prepare weights and unit of measurement
#
if(length(wt)==1)wt <- rep(wt,n)
else if(length(wt)!=n)stop("wt must be the same length as the other variables")
if(min(wt)<0)stop("All weights must be non-negative")
if(length(delta)==1)delta <- rep(delta,n)
else if(length(delta)!=n)stop("delta must be the same length as the other variables")
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
if(is.null(scale)){
	if(distribution=="binomial"||distribution=="beta binomial"||
		distribution=="double binomial"||distribution=="mult binomial")
		scale <- "logit"
	else if(distribution=="normal"||distribution=="logistic"||
		distribution=="Cauchy"||distribution=="Laplace")
		scale <- "identity"
	else scale <- "log"}
#
# modify location function to proper scale for random effect
#
mu4 <- if(scale=="logit") function(p){
		pp <- exp(log(mu1(p)/(1-mu1(p)))[i]+p[npl]*sd)
		pp/(1+pp)}
	else if(scale=="identity") function(p) mu1(p)[i]+p[npl]*sd
	else if(scale=="log") function(p) exp(log(mu1(p))[i]+p[npl]*sd)
	else if(scale=="reciprocal") function(p) 1/(1/mu1(p)[i]+p[npl]*sd)
	else if(scale=="exp") function(p) log(exp(mu1(p))[i]+p[npl]*sd)
#
# check that location function returns appropriate values
#
if(any(is.na(mu4(pmu))))stop("The location regression returns NAs")
if(distribution=="Levy"&&((!censor&&any(y<=mu1(pmu)))||(censor&&any(y[,1]<=mu1(pmu)))))
	stop("location parameter must be strictly less than corresponding observation")
#
# check that shape function returns appropriate values
#
if(distribution=="Pareto"&&exp(sh1(p))<=1)stop("shape parameters must be > 0")
if(distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"){
	if(any(is.na(sh1(p))))stop("The shape regression returns NAs")
	if(od)stop("Some individuals must have more than one observation")}
#
# create the appropriate likelihood function of one observation
#
if (!censor){
	ret <- switch(distribution,
	binomial={
		fcn <- function(p) {
			m <- mu4(p)
			-wt*(y[,1]*log(m)+y[,2]*log(1-m))}
		const <- -wt*lchoose(nn,y[,1])},
	"beta binomial"={
		fcn <- function(p) {
			m <- mu4(p)
			s <- exp(sh1(p))
			t <- s*m
			u <- s*(1-m)
			-wt*(lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u))}
		const <- -wt*lchoose(nn,y[,1])},
	"double binomial"={
		fcn <- function(p) {
			-.C("ddb_c",as.integer(y[,1]),as.integer(nn),
				as.double(mu4(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),
				#DUP=FALSE,
				PACKAGE="repeated")$res}
		const <- 0},
	"mult binomial"={
		fcn <- function(p) {
			-.C("dmb_c",as.integer(y[,1]),as.integer(nn),
				as.double(mu4(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),
				#DUP=FALSE,
				PACKAGE="repeated")$res}
		const <- 0},
	Poisson={
		fcn <- function(p) {
			m <- mu4(p)
			-wt*(-m+y*log(m))}
		const <- wt*lgamma(y+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu4(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(lgamma(y+s)-lgamma(s)+s*t+y*log(m)
				-(y+s)*log(s+m))}
		const <- wt*lgamma(y+1)},
	"double Poisson"={
		fcn <- function(p) {
			-.C("ddp_c",as.integer(y),as.integer(my),
				as.double(mu4(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),
				#DUP=FALSE,
				PACKAGE="repeated")$res}
		const <- 0},
	"mult Poisson"={
		fcn <- function(p) {
			-.C("dmp_c",as.integer(y),as.integer(my),
				as.double(mu4(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),
				#DUP=FALSE,
				PACKAGE="repeated")$res}
		const <- 0},
	"gamma count"={
		fcn <- function(p) {
			m <- mu4(p)
			s <- exp(sh1(p))
			u <- m*s
			-wt*log(ifelse(y==0,1-pgamma(u,(y+1)*s,1),
				pgamma(u,y*s+(y==0),1)-
				pgamma(u,(y+1)*s,1)))}
		const <- 0},
	Consul={
		fcn <- function(p) {
			m <- mu4(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-y*t)}
		const <- wt*lgamma(y+1)},
	logarithmic={
		fcn <- function(p) {
			m <- exp(mu4(p))
			m <- m/(1+m)
			-wt*(y*log(m)-log(y)-log(-log(1-m)))}
		const <- 0},
	geometric={
		fcn <- function(p) {
			m <- mu4(p)
			-wt*(y*log(m)-(y+1)*log(1+m))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				-wt*log(pnorm(y+delta/2,m,s)
					-pnorm(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(t+(y-mu4(p))^2/exp(t))/2}
			const <- wt*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*log(pinvgauss(y+delta/2,m,s)-
					pinvgauss(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				wt*(t+(y-m)^2/(y*exp(t)*m^2))/2}
			const <- wt*(log(2*pi*y^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*log(plogis(y+delta/2,m,s)
					-plogis(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)*sqrt(3)/pi
				wt*((y-m)/s+t+2*log(1+exp(-(y-m)/s)))}
			const <- -wt*(log(pi/sqrt(3))+log(delta))}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				-wt*log(pcauchy(y+delta/2,m,s)
					-pcauchy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				wt*log(s*(1+(y-m)^2/s^2))}
			const <- -wt*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*log(plaplace((y+delta/2-m)/s)
					-plaplace((y-delta/2-m)/s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(abs(y-mu4(p))/exp(t)+t)}
			const <- -wt*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*log(plevy(y+delta/2,m,s)
					-plevy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*(0.5*log(s/(2*pi))-1.5*log(y-m)-
					s/(2*(y-m)))}
			const <- -wt*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu4(p)*(s-1))
				-wt*log((1+(y-delta/2)*t)^-s
					-(1+(y+delta/2)*t)^-s)}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu4(p)*(s-1))
				-wt*(log(s*t)-(s+1)*log(1+y*t))}
			const <- -wt*log(delta)}},
        exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				-wt*log(-exp(-(y+delta/2)/m)
					+exp(-(y-delta/2)/m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				wt*(log(m)+y/m)}
			const <- -wt*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*log(pgamma(y+delta/2,s,scale=u)
					-pgamma(y-delta/2,s,scale=u))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(s*(t-log(m)-y/m)+(s-1)*log(y)-lgamma(s))}
			const <- -wt*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*log(pweibull(y+delta/2,s,m)
					-pweibull(y-delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+(s-1)*log(y)-s*log(m)-(y/m)^s)}
			const <- -wt*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				-wt*log(pweibull(ey+ey*delta/2,s,m)
					-pweibull(ey-ey*delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+s*y-s*log(m)-(exp(y)/m)^s)}
			const <- -wt*log(delta)}},
        beta={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				m <- mu4(p)*s
				s <- s-m
				-sum(wt*log(pbeta(y+delta/2,m,s)
					-pbeta(y-delta/2,m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				m <- mu4(p)*s
				s <- s-m
				-sum(wt*dbeta(y,m,s,0,TRUE))}
			const <- -wt*log(delta)}},
        simplex={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-sum(wt*log(psimplex(y+delta/2,m,s)
					-psimplex(y-delta/2,m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				sum(wt*(((y-m)/(m*(1-m)))^2/(y*(1-y)*s)+
					t+3*log(y*(1-y)))/2)}
			const <- wt*(log(2*pi)/2-log(delta))}},
        "two-sided power"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-sum(wt*log(ifelse(y1<m,m*(y1/m)^s,
					1-(1-m)*((1-y1)/(1-m))^s)
					-ifelse(y2<m,m*(y2/m)^s,
					1-(1-m)*((1-y2)/(1-m))^s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-sum(wt*(t+ifelse(y<m,(s-1)*log(y/m),
					(s-1)*log((1-y)/(1-m)))))}
			const <- -wt*log(delta)}})}
else {
	# censored data
	ret <- switch(distribution,
	Poisson={
		fcn <- function(p) {
			m <- mu4(p)
			-wt*(cc*(-m+y[,1]*log(m))+
				log(lc-rc*ppois(y[,1],m)))}
		const <- wt*cc*lgamma(y[,1]+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu4(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(cc*(lgamma(y[,1]+s)-lgamma(s)
				+s*t+y[,1]*log(m)-(y[,1]+s)*log(s+m))+
				log(lc-rc*pnbinom(y[,1],s,1/(1+m/s))))}
		const <- wt*cc*lgamma(y[,1]+1)},
	geometric={
		fcn <- function(p) {
			m <- mu4(p)
			-wt*(cc*(y[,1]*log(m)-(y[,1]+1)*log(1+m))+
				log(lc-rc*pgeom(y[,1],1/(1+m))))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pnorm(y[,1]+delta/2,m,s)-
					pnorm(y[,1]-delta/2,m,s))
					+log(lc-rc*pnorm(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/s)/2)+log(lc-rc
					*pnorm(y[,1],m,sqrt(s))))}
			const <- wt*cc*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pinvgauss(y[,1]+delta/2,m,s)-
					pinvgauss(y[,1]-delta/2,m,s))
					+log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/(y[,1]*s*m^2))/2)+
					log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- wt*cc*(log(2*pi*y[,1]^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*(cc*log(plogis(y[,1]+delta/2,m,s)-
					plogis(y[,1]-delta/2,m,s))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				y1 <- (y[,1]-m)/s
				-wt*(cc*(-y1-log(s)-2*log(1+exp(-y1)))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- -wt*cc*log(delta)}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pcauchy(y[,1]+delta/2,m,s)-
					pcauchy(y[,1]-delta/2,m,s))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p)/2)
				-wt*(-cc*log(s*(1+(y[,1]-m)^2/s^2))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plaplace((y[,1]+delta/2-m)/s)-
					plaplace((y[,1]-delta/2-m)/s))
					+log(lc-rc*plaplace((y[,1]-m)/s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-abs(y[,1]-m)/s-t)+log(lc-rc
					*plaplace((y[,1]-m)/s)))}
			const <- -wt*cc*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plevy(y[,1]+delta/2,m,s)-
					plevy(y[,1]-delta/2,m,s))
					+log(lc-rc*plevy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(0.5*log(s/(2*pi))-1.5*log(y[,1]-m)-
					s/(2*(y[,1]-m)))+log(lc-rc
					*plevy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu4(p)*(s-1))
				-wt*(cc*log((1+(y[,1]-delta/2)*t)^-s-
					(1+(y[,1]+delta/2)*t)^-s)
					+log(lc-rc*(-(1+(y[,1])*t)^-s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu4(p)*(s-1))
				-wt*(cc*(log(s*t)-(s+1)*log(1+y[,1]*t))
					+log(lc-rc*(1-(1+y[,1]*t)^-s)))}
			const <- -wt*cc*log(delta)}},
	exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				-wt*(cc*log(-exp(-(y[,1]+delta/2)/m)
					+exp(-(y[,1]-delta/2)/m))+
					log(lc-rc*(1-exp(-y[,1]/m))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				-wt*(cc*(-log(m)-y[,1]/m)+log(lc-rc*
					(1-exp(-y[,1]/m))))}
			const <- -wt*cc*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*(cc*log(pgamma(y[,1]+delta/2,s,scale=u)-
					pgamma(y[,1]-delta/2,s,scale=u))
					+log(lc-rc*pgamma(y[,1],s,scale=u)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(s*(t-log(m)-y[,1]/m)+(s-1)*log(y[,1])
					-lgamma(s))+log(lc-rc
					*pgamma(y[,1],s,scale=m/s)))}
			const <- -wt*cc*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pweibull(y[,1]+delta/2,s,m)-
					pweibull(y[,1]-delta/2,s,m))
					+log(lc-rc*pweibull(y[,1],s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(t+(s-1)*log(y[,1])-s*log(m)
					-(y[,1]/m)^s)+log(lc-rc*
					pweibull(y[,1],s,m)))}
			const <- -wt*cc*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu4(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				-wt*(cc*log(pweibull(ey+ey*delta/2,s,m)-
					pweibull(ey-ey*delta/2,s,m))
					+log(lc-rc*pweibull(ey,s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu4(p)
				t <- sh1(p)
				s <- exp(t)
				ey <- exp(y[,1])
				-wt*(cc*(t+s*y[,1]-s*log(m)-(ey/m)^s)+log(lc-
					rc*pweibull(ey,s,m)))}
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
z0 <- nlm(fn, p=p, hessian=TRUE, print.level=print.level, typsize=typsize,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
#
# calculate fitted values and raw residuals
#
fitted.values <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial")
		as.vector((y[,1]+y[,2])*mu4(z0$estimate))
	else as.vector(mu4(z0$estimate))
residuals <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial"||censor)
		y[,1]-fitted.values
	else y-fitted.values
#
# calculate se's
#
if(npl+nps==1){
	cov <- 1/z0$hessian
	se <- sqrt(cov)}
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==npl+nps)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=npl+nps,nrow=npl+nps)
	se <- sqrt(diag(cov))}
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
if(!is.null(sh2))sh1 <- sh2
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	respname=respname,
	mu=mu1,
	shape=sh1,
	linear=list(lin1,lin2),
	linmodel=list(lin1model,lin2model),
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
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlm"
return(z1)}
