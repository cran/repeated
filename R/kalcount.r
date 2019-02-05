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
#     kalcount(response=NULL, times=NULL, origin=0, intensity="exponential",
#	depend="independence", update="Markov", mu=NULL, shape=NULL,
#	density=FALSE, ccov=NULL, tvcov=NULL, preg=NULL, ptvc=NULL,
#	pbirth=NULL, pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
#	pfamily=NULL, envir=parent.frame(), print.level=0, ndigit=10,
#	gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
#	typsize=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit various distributions inserted into a Pareto
#  distribution with serial dependence or gamma frailties using
#  Kalman-type update for longitudinal count data.



##' Repeated Measurements Models for Counts with Frailty or Serial Dependence
##' 
##' \code{kalcount} is designed to handle repeated measurements models with
##' time-varying covariates. The distributions have two extra parameters as
##' compared to the functions specified by \code{intensity} and are generally
##' longer tailed than those distributions. Dependence among observations on a
##' unit can be through gamma or power variance family frailties (a type of
##' random effect), with or without autoregression, or serial dependence over
##' time.
##' 
##' By default, a gamma mixture of the distribution specified in
##' \code{intensity} is used, as the conditional distribution in the
##' \code{serial} dependence models, and as a symmetric multivariate (random
##' effect) model for \code{frailty} dependence.
##' 
##' Unless specified otherwise, the time origin is taken to be zero. The given
##' times are the ends of the periods in which the counts occurred.
##' 
##' Here, the variance, with exponential intensity, is a quadratic function of
##' the mean, whereas, for \code{\link[repeated]{nbkal}}, it is proportional to
##' the mean function.
##' 
##' If the counts on a unit are clustered, not longitudinal, use the failty
##' dependence with the default exponential intensity, yielding a multivariate
##' negative binomial distribution.
##' 
##' Nonlinear regression models can be supplied as formulae where parameters
##' are unknowns in which case factor variables cannot be used and parameters
##' must be scalars. (See \code{\link[rmutil]{finterp}}.)
##' 
##' Marginal and individual profiles can be plotted using
##' \code{\link[rmutil]{mprofile}} and \code{\link[rmutil]{iprofile}} and
##' residuals with \code{\link[rmutil]{plot.residuals}}.
##' 
##' 
##' @param response A list of two column matrices with counts and corresponding
##' times for each individual, one matrix or dataframe of counts, or an object
##' of class, \code{response} (created by \code{\link[rmutil]{restovec}}) or
##' \code{repeated} (created by \code{\link[rmutil]{rmna}} or
##' \code{\link[rmutil]{lvna}}). If the \code{repeated} data object contains
##' more than one response variable, give that object in \code{envir} and give
##' the name of the response variable to be used here.
##' @param times When response is a matrix, a vector of possibly unequally
##' spaced times when they are the same for all individuals or a matrix of
##' times. Not necessary if equally spaced. Ignored if response has class,
##' \code{response} or \code{repeated}.
##' @param origin If the time origin is to be before the start of observations,
##' a positive constant to be added to all times.
##' @param intensity The form of function to be put in the Pareto distribution.
##' Choices are exponential, Weibull, gamma, log normal, log logistic, log
##' Cauchy, log Student, and gen(eralized) logistic.
##' @param depend Type of dependence. Choices are \code{independence},
##' \code{frailty}, and \code{serial}.
##' @param update Type of for serial dependence. Choices are \code{Markov},
##' \code{serial}, \code{event}, \code{cumulated}, \code{count}, and
##' \code{kalman}. With \code{frailty} dependence, weighting by length of
##' observation time may be specified by setting update to \code{time}.
##' @param mu A regression function for the location parameter or a formula
##' beginning with ~, specifying either a linear regression function in the
##' Wilkinson and Rogers notation (a log link is assumed) or a general function
##' with named unknown parameters. Give the initial estimates in \code{preg} if
##' there are no time-varying covariates and in \code{ptvc} if there are.
##' @param shape A regression function for the shape parameter or a formula
##' beginning with ~, specifying either a linear regression function in the
##' Wilkinson and Rogers notation or a general function with named unknown
##' parameters. It must yield one value per observation.
##' @param density If TRUE, the density of the function specified in
##' \code{intensity} is used instead of the intensity.
##' @param ccov A vector or matrix containing time-constant baseline covariates
##' with one row per individual, a model formula using vectors of the same
##' size, or an object of class, \code{tccov} (created by
##' \code{\link[rmutil]{tcctomat}}). If response has class, \code{repeated},
##' the covariates must be supplied as a Wilkinson and Rogers formula unless
##' none are to be used or \code{mu} is given.
##' @param tvcov A list of matrices with time-varying covariate values,
##' observed in the time periods in \code{response}, for each individual (one
##' column per variable), one matrix or dataframe of such covariate values, or
##' an object of class, \code{tvcov} (created by
##' \code{\link[rmutil]{tvctomat}}). If response has class, \code{repeated},
##' the covariates must be supplied as a Wilkinson and Rogers formula unless
##' none are to be used or \code{mu} is given.
##' @param preg Initial parameter estimates for the regression model: intercept
##' plus one for each covariate in \code{ccov}. If \code{mu} is a formula or
##' function, the parameter estimates must be given here only if there are no
##' time-varying covariates. If \code{mu} is a formula with unknown parameters,
##' their estimates must be supplied either in their order of appearance in the
##' expression or in a named list.
##' @param ptvc Initial parameter estimates for the coefficients of the
##' time-varying covariates, as many as in \code{tvcov}. If \code{mu} is a
##' formula or function, the parameter estimates must be given here if there
##' are time-varying covariates present.
##' @param pbirth If supplied, this is the initial estimate for the coefficient
##' of the birth model.
##' @param pintercept The initial estimate of the intercept for the generalized
##' logistic intensity.
##' @param pshape An initial estimate for the shape parameter of the intensity
##' function (except exponential intensity). If \code{shape} is a function or
##' formula, the corresponding initial estimates. If \code{shape} is a formula
##' with unknown parameters, their estimates must be supplied either in their
##' order of appearance in the expression or in a named list.
##' @param pinitial An initial estimate for the initial parameter. With
##' \code{frailty} dependence, this is the frailty parameter.
##' @param pdepend An initial estimate for the serial dependence parameter. For
##' \code{frailty} dependence, if a value is given here, an autoregression is
##' fitted as well as the frailty.
##' @param pfamily An optional initial estimate for the second parameter of a
##' two-parameter power variance family mixture instead of the default gamma
##' mixture. This yields a gamma mixture as \code{family -> 0}, an inverse
##' Gauss mixture for \code{family = 0.5}, and a compound distribution of a
##' Poisson-distributed number of gamma distributions for \code{-1 < family <
##' 0}.
##' @param envir Environment in which model formulae are to be interpreted or a
##' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
##' name of the response variable should be given in \code{response}. If
##' \code{response} has class \code{repeated}, it is used as the environment.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of classes \code{kalcount} and \code{recursive} is returned.
##' @author J.K. Lindsey
### @seealso \code{\link[growth]{carma}}, \code{\link[growth]{elliptic}},
### \code{\link[rmutil]{finterp}}, \code{\link[repeated]{gar}},
### \code{\link[rmutil]{gettvc}}, \code{\link[repeated]{gnlmm}},
### \code{\link[gnlm]{gnlr}}, \code{\link[rmutil]{iprofile}},
### \code{\link[repeated]{kalseries}}, \code{\link[event]{kalsurv}},
### \code{\link[rmutil]{mprofile}}, \code{\link[repeated]{nbkal}},
### \code{\link[rmutil]{read.list}}, \code{\link[rmutil]{restovec}},
### \code{\link[rmutil]{rmna}}, \code{\link[rmutil]{tcctomat}},
### \code{\link[rmutil]{tvctomat}}.
##' @keywords models
##' @examples
##' 
##' treat <- c(0,0,1,1)
##' tr <- tcctomat(treat)
##' dose <- # matrix(rpois(20,10),ncol=5)
##' 	matrix(c(9,13,16,7,12,6,9,10,11,9,10,10,7,9,9,9,8,10,15,4),
##' 		ncol=5,byrow=TRUE)
##' dd <- tvctomat(dose)
##' y <- # matrix(rpois(20,1+3*rep(treat,5)),ncol=5)
##' 	restovec(matrix(c(1,1,1,1,0,1,0,1,0,5,3,3,4,1,4,4,2,3,2,5),
##' 		ncol=5,byrow=TRUE))
##' reps <- rmna(y, ccov=tr, tvcov=dd)
##' #
##' # log normal intensity, independence model
##' kalcount(y, intensity="log normal", dep="independence", preg=1,
##' 	pshape=0.1)
##' # random effect
##' kalcount(y, intensity="log normal", dep="frailty", pinitial=0.1,
##' 	preg=1, psh=0.1)
##' # serial dependence
##' kalcount(y, intensity="log normal", dep="serial", pinitial=0.1,
##' 	preg=1, pdep=0.01, psh=0.1)
##' # random effect and autoregression
##' kalcount(y, intensity="log normal", dep="frailty", pinitial=0.1,
##' 	pdep=0.1, preg=1, psh=0.1)
##' # add time-constant variable
##' kalcount(y, intensity="log normal", pinitial=0.1, psh=0.1,
##' 	preg=c(1,0), ccov=treat)
##' # or equivalently
##' kalcount(y, intensity="log normal", mu=~treat, pinitial=0.1,
##' 	psh=0.1, preg=c(1,0))
##' # or
##' kalcount(y, intensity="log normal", mu=~b0+b1*treat,
##' 	pinitial=0.1, psh=0.1, preg=c(1,0), envir=reps)
##' # add time-varying variable
##' kalcount(y, intensity="log normal", pinitial=0.1, psh=0.1,
##' 	preg=c(1,0), ccov=treat, ptvc=0, tvc=dose)
##' # or equivalently, from the environment
##' dosev <- as.vector(t(dose))
##' kalcount(y, intensity="log normal", mu=~b0+b1*rep(treat,rep(5,4))+b2*dosev,
##' 	pinitial=0.1, psh=0.1, ptvc=c(1,0,0))
##' # or from the reps data object
##' kalcount(y, intensity="log normal", mu=~b0+b1*treat+b2*dose,
##' 	pinitial=0.1, psh=0.1, ptvc=c(1,0,0), envir=reps)
##' # try power variance family
##' kalcount(y, intensity="log normal", mu=~b0+b1*treat+b2*dose,
##' 	pinitial=0.001, psh=14, ptvc=c(5,-1,0.1), envir=reps,
##' 	pfamily=0.8)
##' 
##' @aliases kalcount deviance.kalcount residuals.kalcount fitted.kalcount print.kalcount 
##' @export kalcount
kalcount <- function(response=NULL, times=NULL, origin=0,
	intensity="exponential", depend="independence", update="Markov",
	mu=NULL, shape=NULL, density=FALSE, ccov=NULL, tvcov=NULL, preg=NULL,
	ptvc=NULL, pbirth=NULL, pintercept=NULL, pshape=NULL, pinitial=1,
	pdepend=NULL, pfamily=NULL, envir=parent.frame(), print.level=0,
	ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
	typsize=abs(p), stepmax=10*sqrt(p%*%p)){
#
# Pareto (beta) likelihood function for serial dependence
#
kcountb <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("kcountb_c",
		p=as.double(p),
		y=as.double(resp$response$times),
		origin=as.double(origin),
		c=as.integer(resp$response$y),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		model=as.integer(mdl),
		density=as.integer(density),
		pfamily=as.integer(!is.null(pfamily)),
		dep=as.integer(dep),
		birth=as.integer(birth),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")
	z$like}
#
# Pareto (beta) likelihood function for frailty
#
countfb <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("countfb_c",
		p=as.double(p),
		y=as.double(resp$response$times),
		c=as.integer(resp$response$y),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		model=as.integer(mdl),
		density=as.integer(density),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		frser=as.integer(frser),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")
	z$like}
call <- sys.call()
#
# check model
#
tmp <- c("exponential", "Weibull","gamma","gen logistic","log normal",
	"log logistic","log Cauchy","log Laplace")
mdl <- match(intensity <- match.arg(intensity,tmp),tmp)
depend <- match.arg(depend,c("independence","serial","frailty"))
tmp <- c("Markov","serial","event","cumulated","count","kalman","time")
dep <- match(update <- match.arg(update,tmp),tmp)
rf <- !is.null(mu)
sf <- !is.null(shape)
#
# check if a data object is being supplied
#
type <- "unknown"
respenv <- exists(deparse(substitute(response)),envir=parent.frame())&&
	inherits(response,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(response$response$y)[2]>1)
		stop("kalcount only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("kalcount does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)deparse(substitute(response))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# if envir, remove extra (multivariate) responses
#
if(!respenv&&inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("kalcount does not handle data with NAs")
	cn <- deparse(substitute(response))
	if(length(grep("\"",cn))>0)cn <- response
	if(length(cn)>1)stop("only one response variable allowed")
	response <- envir
	col <- match(cn,colnames(response$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- response$response$type[col]
	if(dim(response$response$y)[2]>1){
		response$response$y <- response$response$y[,col,drop=FALSE]}}
#
# find the covariates
#
if(respenv||inherits(envir,"repeated")){
	if(rf)resp <- rmna(response$response)
	else {
		resp <- response
		# time-constant covariates
		if(is.null(ccov))resp$ccov <- NULL
		else if(inherits(ccov,"formula")){
			if(any(is.na(match(rownames(attr(terms(ccov,data=response),"factors")),colnames(response$ccov$ccov)))))
				stop("ccov covariate(s) not found")
			tmp <- wr(ccov,data=response,expand=FALSE)$design
			resp$ccov$ccov <- tmp[,-1,drop=FALSE]
			rm(tmp)}
		else stop("ccov must be a W&R formula")
		# time-varying covariates
		if(is.null(tvcov))resp$tvcov <- NULL
		else if(inherits(tvcov,"formula")){
			if(any(is.na(match(rownames(attr(terms(tvcov,data=response),"factors")),colnames(response$tvcov$tvcov)))))
				stop("tvcov covariate(s) not found")
			tmp <- wr(tvcov,data=response)$design
			resp$tvcov$tvcov <- tmp[,-1,drop=FALSE]
			rm(tmp)}
		else stop("tvcov must be a W&R formula")}
	nccov <- if(rf||is.null(resp$ccov$ccov)) 0
		 else  dim(resp$ccov$ccov)[2]
	ttvc <- if(rf||is.null(resp$tvcov$tvcov)) 0
		 else  dim(resp$tvcov$tvcov)[2]}
else {
	# make response object
	if(inherits(response,"response")){
		resp <- response
		type <- response$type}
	else resp <- restovec(response,times)
	# time-constant covariates
	if(is.null(ccov))nccov <- 0
	else {
		if(!inherits(ccov,"tccov")){
			ccname <- deparse(substitute(ccov))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- deparse(substitute(ccov))
				if(dim(ccov)[2]>1){
					tmp <- NULL
					for(i in 1:dim(ccov)[2])tmp <- c(tmp,paste(ccname,i,sep=""))
					ccname <- tmp}}
			ccov <- tcctomat(ccov,names=ccname)}
		nccov <- if(rf) 0 else dim(ccov$ccov)[2]}
	# time-varying covariates
	if(is.null(tvcov))ttvc <- 0
	else {
		if(!inherits(tvcov,"tvcov")){
			tvcname <- deparse(substitute(tvcov))
			if(is.list(tvcov)&&dim(tvcov[[1]])[2]>1){
				if(is.null(colnames(tvcov[[1]]))){
					tvcname <- deparse(substitute(tvcov))
					tmp <- NULL
					for(i in 1:dim(tvcov[[1]])[2])tmp <- c(tmp,paste(tvcname,i,sep=""))
					tvcname <- tmp}
				else tvcname <- colnames(tvcov[[1]])}
			tvcov <- tvctomat(tvcov, names=tvcname)}
		ttvc <- if(rf) 0 else dim(tvcov$tvcov)[2]}
	resp <- rmna(response=resp, tvcov=tvcov, ccov=ccov)
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
n <- dim(resp$response$y)[1]
nobs <- nobs(resp)
nind <- length(nobs)
if((inherits(envir,"repeated")&&(length(nobs)!=length(nobs(envir))||
	any(nobs!=nobs(envir))))||(inherits(envir,"tvcov")&&
	(length(nobs)!=length(envir$tvcov$nobs)||any(nobs!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
if(is.null(resp$response$times))
	stop("These models cannot be fitted without times.")
if(!is.null(resp$response$nest)&&depend=="serial")
	stop("kalcount does not handle two levels of nesting with serial dependence")
if(type!="unknown"&&type!="discrete")stop("discrete data required")
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu3 <- sh3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	if(inherits(mu,"formula")){
		mu3 <- if(respenv)finterp(mu,.envir=response,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	else if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
		        tmp <- parse(text=deparse(mu)[-1])
		        mu <- if(respenv)fnenvir(mu,.envir=response,.name=envname)
		        	else fnenvir(mu,.envir=envir,.name=envname)
		        mu3 <- mu
		        attr(mu3,"model") <- tmp}
		else mu3 <- mu}
	if(inherits(shape,"formula")){
		sh3 <- if(respenv)finterp(shape,.envir=response,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
	else if(is.function(shape)){
		if(is.null(attr(shape,"model"))){
		        tmp <- parse(text=deparse(shape)[-1])
		        shape <- if(respenv)fnenvir(shape,.envir=response,.name=envname)
		        	else fnenvir(shape,.envir=envir,.name=envname)
		        sh3 <- shape
		        attr(sh3,"model") <- tmp}
		else sh3 <- shape}}
#
# transform location formula to function and check number of parameters
#
npreg <- length(preg)
mu1 <- sh1 <- v <- b <- NULL
if(inherits(mu,"formula")){
	pr <- if(npreg>0)preg else ptvc
	npr <- length(pr)
	mu2 <- if(respenv)finterp(mu,.envir=response,.name=envname,.expand=is.null(preg))
		else finterp(mu,.envir=envir,.name=envname,.expand=is.null(preg))
	npt1 <- length(attr(mu2,"parameters"))
	if(is.character(attr(mu2,"model"))){
	# W&R formula
		if(length(attr(mu2,"model"))==1){
			mu1 <- function(p) exp(p[1]*rep(1,n))
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}
		else {
			mu1 <- function(p) exp(mu2(p))
			attributes(mu1) <- attributes(mu2)}}
	else {
	# formula with unknowns
		if(npr!=npt1&&length(ptvc)!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg or ptvc should have",npt1,"estimates"))}
		if(is.list(pr)){
			if(!is.null(names(pr))){
				o <- match(attr(mu2,"parameters"),names(pr))
				pr <- unlist(pr)[o]
				if(sum(!is.na(o))!=length(pr))
					stop("invalid estimates for mu - probably wrong names")}
			else pr <- unlist(pr)
			if(npreg>0)preg <- pr else ptvc <- pr}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else if(!is.character(attr(mu2,"model"))){
			mu1 <- mu2
			rm(mu2)}}}
else if(is.function(mu))mu1 <- mu
#
# give appropriate attributes to mu1 for printing
#
if(!is.null(mu1)&&is.null(attr(mu1,"parameters"))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=response))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,.envir=response))
			else attributes(fnenvir(mu1,.envir=envir))}}
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)){
	if(is.null(ptvc)&&nlp!=npreg)
		stop(paste("preg should have",nlp,"initial estimates"))
	else if(!is.null(ptvc)&&length(ptvc)!=nlp)
		stop(paste("ptvc should have",nlp,"initial estimates"))}
#
# transform shape formula to function and check number of parameters
#
nps <- length(pshape)
if(inherits(shape,"formula")){
	sh2 <- if(respenv)finterp(shape,.envir=response,.name=envname)
		else finterp(shape,.envir=envir,.name=envname)
	npt2 <- length(attr(sh2,"parameters"))
	# W&R formula
	common <- FALSE
	npl1 <- if(common&&!inherits(shape,"formula")) 1 else nlp+1
	if(is.character(attr(sh2,"model"))){
		if(length(attr(sh2,"model"))==1){
			sh1 <- function(p) p[npl1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
	else {
	# formula with unknowns
		if(nps!=npt2){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))
					stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(is.function(shape))sh1 <- shape
#
# give appropriate attributes to sh1 for printing
#
if(!is.null(sh1)&&is.null(attr(sh1,"parameters")))
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn")){
			if(respenv)attributes(fnenvir(shape,.envir=response))
			else attributes(fnenvir(shape,.envir=envir))}
		else attributes(shape)}
		else {
			if(respenv)attributes(fnenvir(sh1,.envir=response))
			else attributes(fnenvir(sh1,.envir=envir))}
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(shape))length(attr(sh1,"parameters"))
	else if(is.null(shape))NULL
	else npt2
if(!is.null(nlp)&&nlp!=nps)
	stop(paste("pshape should have",nlp,"initial estimates"))
#
# check that appropriate functions supplied
#
if(rf&&!is.function(mu1))stop("mu must be a formula or function")
if(sf&&!is.function(sh1))stop("shape must be a formula or function")
if(origin<0)stop("Origin must be positive")
birth <- !is.null(pbirth)
tvc <- length(ptvc)
if(!rf&&(ttvc>0&&tvc!=ttvc||ttvc==0&&tvc>0))
	stop(paste(ttvc,"initial estimates of coefficients for time-varying covariates must be supplied"))
if(rf&&birth)stop("Birth models cannot be fitted with a mean function")
if(intensity=="exponential"){
	sf <- FALSE
	pshape <- NULL}
else {
	if(is.null(pshape))
		stop("Initial value of the shape parameter must be supplied")
	if(!sf){
		if(pshape<=0)stop("shape must be positive")
		pshape <- log(pshape)}}
if(intensity=="gen logistic"){
	if(is.null(pintercept))stop("Initial value of the intercept parameter must be supplied")}
else pintercept <- NULL
if(pinitial<=0)stop("Estimate of initial parameter must greater than 0")
else pinitial <- log(pinitial)
#
# check what dependence is required
#
frser <- FALSE
if(depend=="independence"){
	pdepend <- NULL
	dep <- 0}
else if(depend=="serial"){
	if(update=="time")stop("time update can only be used with frailty")
	if(is.null(pdepend))
		stop("An estimate of the dependence parameter must be supplied")
	else if(pdepend<=0|pdepend>=1)
		stop("Dependence parameter must be between 0 and 1")
	else pdepend <- log(pdepend/(1-pdepend))}
else if(depend=="frailty"){
	frser <- !is.null(pdepend)
	if(frser){
		if(pdepend<=0)stop("Dependence parameter must be positive")
		pdepend <- log(pdepend)}
	if(update=="time")dep <- 1
	else {
		dep <- 0
		update <- "no"}}
#
# check number of parameter estimates
#
if(!is.null(pfamily)&&depend=="frailty")
		stop("pfamily cannot be used with frailty dependence")
if(rf&&npreg>0)nccov <- npreg-1
if(!rf&&1+nccov!=npreg)
	stop(paste(1+nccov,"regression estimates must be supplied"))
#
# set up for location function and make sure it gives correct number
# of values
#
if(rf){
	if(tvc>0&&nccov>0)
		stop("With a mean function, initial estimates must be supplied either in preg or in ptvc")
	if(tvc>0){
		if(length(mu1(ptvc))!=n)
			stop("The mu function or formula must provide an estimate for each observation")
		tvc <- tvc-1}
	else {
		lp <- length(mu1(preg))
		if(lp==1){
			if(nccov==0)mu1 <- function(p) rep(p[1],nind)
			else stop("Number of estimates does not correspond to mu function")}
		else if(lp!=nind)
		     stop("The mu function or formula must provide an estimate for each individual")}}
if(sf&&length(sh1(pshape))!=n)
	stop("The shape function must provide an estimate for each observation")
#
# check that the likelihood function gives appropriate values and call nlm
#
np <- 1+nccov+tvc+1+(depend=="serial")+birth+nps+(!is.null(pintercept))+frser+
	(!is.null(pfamily))
nps1 <- np-nps+1
p <- c(preg,pbirth,ptvc,pinitial,pdepend,pfamily,pintercept,pshape)
if(depend=="frailty")count <- countfb
else count <- kcountb
if(fscale==1)fscale <- count(p)
if(is.na(count(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(count, p=p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
#
# calculate se's
#
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
#
# calculate recursive fitted values
#
z <- if(depend=="frailty"){
	if(rf)b <- mu1(z0$estimate)
	if(sf)v <- sh1(z0$estimate[nps1:np])
	.C("countfb_c",
        	p=as.double(z0$estimate),
        	y=as.double(resp$response$times),
        	c=as.integer(resp$response$y),
        	x=as.double(resp$ccov$ccov),
        	nind=as.integer(nind),
        	nobs=as.integer(nobs),
        	nbs=as.integer(n),
        	nccov=as.integer(nccov),
        	model=as.integer(mdl),
        	density=as.integer(density),
        	tvc=as.integer(tvc),
        	tvcov=as.double(resp$tvcov$tvcov),
        	fit=as.integer(1),
        	pred=double(n),
		rpred=double(n),
        	rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		frser=as.integer(frser),
        	like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")}
else {
	if(rf)b <- mu1(z0$estimate)
	if(sf)v <- sh1(z0$estimate[nps1:np])
	.C("kcountb_c",
        	p=as.double(z0$estimate),
        	y=as.double(resp$response$times),
        	origin=as.double(origin),
        	c=as.integer(resp$response$y),
        	x=as.double(resp$ccov$ccov),
        	nind=as.integer(nind),
        	nobs=as.integer(nobs),
        	nbs=as.integer(n),
        	nccov=as.integer(nccov),
        	model=as.integer(mdl),
        	density=as.integer(density),
		pfamily=as.integer(!is.null(pfamily)),
        	dep=as.integer(dep),
        	birth=as.integer(birth),
        	tvc=as.integer(tvc),
        	tvcov=as.double(resp$tvcov$tvcov),
        	fit=as.integer(1),
        	pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")}
#
# return appropriate attributes on functions
#
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
z <- list(
	call=call,
	intensity=intensity,
	pfamily=!is.null(pfamily),
	mdl=mdl,
	frser=frser,
	mu=mu1,
	npr=1+nccov+tvc+birth,
	shape=sh1,
	nps=nps,
	density=density,
	depend=depend,
	update=update,
	birth=birth,
	response=resp$response,
	pred=z$pred,
	rpred=z$rpred,
	ccov=resp$ccov,
	tvcov=resp$tvcov,
	maxlike=z0$minimum,
	aic=z0$minimum+np,
	df=n-np,
	npt=np,
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=corr,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- c("kalcount","recursive")
return(z)}

### standard methods
###
##' @export
deviance.kalcount <- function(object, ...) 2*object$maxlike
##' @export
fitted.kalcount <- function(object, recursive=TRUE, ...)
if(recursive) object$rpred else object$pred

##' @export
residuals.kalcount <- function(object, type = "response", recursive=TRUE, ...){
if(type=="response") object$response$y-object$rpred
else (object$response$y-object$rpred)/sqrt(object$rpred)}

### print method
##' @export
print.kalcount <- function(x,digits=max(3,.Options$digits-3),correlation=TRUE,...){
  z<-x
if(!is.null(z$ccov))nccov <- dim(z$ccov$ccov)[2]
else nccov <- 0
expm <- z$intensity!="exponential"&&!is.function(z$shape)
glm <- z$intensity=="gen logistic"
npt <- if(is.function(z$shape)) z$npt-z$nps else z$npt
deppar <- (z$depend=="serial"||z$depend=="Markov")
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
if(z$density)cat(z$intensity," density",sep="")
else cat(z$intensity," intensity",sep="")
if(z$depend=="independence")cat(" with independence\n")
else if(z$depend=="frailty")
	cat(" with",z$depend,"dependence",if(!z$frser)"and",z$update,
		"weight",if(z$frser)"and AR","\n")
else cat(" with ",z$update," update\n",sep="")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Location parameters\n")
if(!is.null(attr(z$mu,"formula")))
	cat(deparse(attr(z$mu,"formula")),sep="\n")
else if(!is.null(attr(z$mu,"model"))){
	t <- deparse(attr(z$mu,"model"))
	t[1] <- sub("expression\\(","",t[1])
	t[length(t)] <- sub("\\)$","",t[length(t)])
	cat(t,sep="\n")}
coef.table <- cbind(z$coef[1:z$npr],z$se[1:z$npr])
if(inherits(z$mu,"formulafn"))
	cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else attr(z$mu,"parameters")
else {
	cname <- "(Intercept)"
	if(nccov)cname <- c(cname,colnames(z$ccov$ccov))
	if(z$birth)cname <- c(cname,"birth")
	if(!is.null(z$tvcov))cname <- c(cname,colnames(z$tvcov$tvcov))}
dimnames(coef.table) <- list(cname, c("estimate","se"))
print.default(coef.table, digits=digits, print.gap=2)
if(is.function(z$shape))cat("\nDependence parameters\n")
else cat("\nNonlinear parameters\n")
coef <- exp(z$coef[(npt-deppar-z$frser-z$pfamily-glm-expm):npt])
cname <- "initial"
if(deppar){
	coef[2] <- coef[2]/(1+coef[2])
	cname <- c(cname,"depend")}
if(z$frser)cname <- c(cname,"AR")
if(z$pfamily){
	cname <- c(cname,"family")
	coef[length(coef)-glm-expm] <- z$coef[npt-glm-expm]}
if(glm){
	cname <- c(cname,"intercept")
	coef[length(coef)-expm] <- z$coef[npt-expm]
	if(expm){
		cname <- c(cname,"asymptote")
		coef[length(coef)] <- 1/coef[length(coef)]}}
else if(expm)cname <- c(cname,"shape")
coef.table <- cbind(z$coef[(npt-deppar-z$frser-z$pfamily-glm-expm):npt],
	z$se[(npt-deppar-z$frser-z$pfamily-glm-expm):npt],coef)
dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
print.default(coef.table, digits=digits, print.gap=2)
if(inherits(z$shape,"formulafn")){
	cat("\nShape parameters\n")
	if(!is.null(attr(z$shape,"formula")))
		cat(deparse(attr(z$shape,"formula")),sep="\n")
	else if(!is.null(attr(z$shape,"model"))){
		t <- deparse(attr(z$shape,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$shape,"model")))
		attr(z$shape,"model")
		else attr(z$shape,"parameters")
	coef.table <- cbind(z$coef[(z$npt-z$nps+1):z$npt],
		z$se[(z$npt-z$nps+1):z$npt])
	dimnames(coef.table) <- list(cname, c("estimate","se"))
	print.default(coef.table, digits=digits, print.gap=2)}
if(correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}}
