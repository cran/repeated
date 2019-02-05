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
#     kalseries(response=NULL, times=NULL, intensity="exponential",
#	depend="independence", mu=NULL, shape=NULL, density=FALSE, ccov=NULL,
#	tvcov=NULL, torder=0, interaction=NULL, preg=NULL, ptvc=NULL,
#	pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
#	pfamily=NULL, delta=NULL, transform="identity", link="identity",
#	envir=parent.frame(), print.level=0, ndigit=10, gradtol=0.00001,
#	steptol=0.00001, fscale=1, iterlim=100, typsize=abs(p),
#	stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit various distributions inserted into a Pareto
#  distribution with serial dependence or gamma frailties using
#  Kalman-type update for continuous longitudinal data.

# .First.lib <- function(lib, pkg)
# 	library.dynam("repeated", pkg, lib)
# require(rmutil)



##' Repeated Measurements Models for Continuous Variables with Frailty or
##' Serial Dependence
##' 
##' \code{kalseries} is designed to handle repeated measurements models with
##' time-varying covariates. The distributions have two extra parameters as
##' compared to the functions specified by \code{intensity} and are generally
##' longer tailed than those distributions. Dependence among observations on a
##' unit can be through gamma or power variance family frailties (a type of
##' random effect), with or without autoregression, or one of two types of
##' serial dependence over time.
##' 
##' By default, a gamma mixture of the distribution specified in
##' \code{intensity} is used, as the conditional distribution in the
##' \code{Markov} and \code{serial} dependence models, and as a symmetric
##' multivariate (random effect) model for \code{frailty} dependence. For
##' example, with a Weibull \code{intensity} and \code{frailty} dependence,
##' this yields a multivariate Burr distribution and with \code{Markov} or
##' \code{serial} dependence, univariate Burr conditional distributions.
##' 
##' If a value for \code{pfamily} is used, the gamma mixture is replaced by a
##' power variance family mixture.
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
##' @param response A list of two column matrices with responses and
##' corresponding times for each individual, one matrix or dataframe of
##' response values, or an object of class, \code{response} (created by
##' \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
##' \code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
##' \code{repeated} data object contains more than one response variable, give
##' that object in \code{envir} and give the name of the response variable to
##' be used here.
##' @param times When response is a matrix, a vector of possibly unequally
##' spaced times when they are the same for all individuals or a matrix of
##' times. Not necessary if equally spaced. Ignored if response has class,
##' \code{response} or \code{repeated}.
##' @param intensity The form of function to be put in the Pareto distribution.
##' Choices are exponential, Weibull, gamma, normal, logistic, Cauchy, log
##' normal, log logistic, log Cauchy, log Student, inverse Gauss, and
##' gen(eralized) logistic. (For definitions of distributions, see the
##' corresponding [dpqr]distribution help.)
##' @param depend Type of dependence. Choices are \code{independence},
##' \code{Markov}, \code{serial}, and \code{frailty}.
##' @param mu A regression function for the location parameter or a formula
##' beginning with ~, specifying either a linear regression function in the
##' Wilkinson and Rogers notation or a general function with named unknown
##' parameters. Give the initial estimates in \code{preg} if there are no
##' time-varying covariates and in \code{ptvc} if there are.
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
##' observed at the event times in \code{response}, for each individual (one
##' column per variable), one matrix or dataframe of such covariate values, or
##' an object of class, \code{tvcov} (created by
##' \code{\link[rmutil]{tvctomat}}). If a time-varying covariate is observed at
##' arbitrary time, \code{\link[rmutil]{gettvc}} can be used to find the most
##' recent values for each response and create a suitable list. If response has
##' class, \code{repeated}, the covariates must be supplied as a Wilkinson and
##' Rogers formula unless none are to be used or \code{mu} is given.
##' @param torder The order of the polynomial in time to be fitted.
##' @param interaction Vector of length equal to the number of time-constant
##' covariates, giving the levels of interactions between them and the
##' polynomial in time in the \code{linear model}.
##' @param preg Initial parameter estimates for the regression model:
##' intercept, one for each covariate in \code{ccov}, and \code{torder} plus
##' sum(\code{interaction}). If \code{mu} is a formula or function, the
##' parameter estimates must be given here only if there are no time-varying
##' covariates. If \code{mu} is a formula with unknown parameters, their
##' estimates must be supplied either in their order of appearance in the
##' expression or in a named list.
##' @param ptvc Initial parameter estimates for the coefficients of the
##' time-varying covariates, as many as in \code{tvcov}. If \code{mu} is a
##' formula or function, the parameter estimates must be given here if there
##' are time-varying covariates present.
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
##' @param delta Scalar or vector giving the unit of measurement for each
##' response value, set to unity by default. For example, if a response is
##' measured to two decimals, delta=0.01. If the response has been
##' pretransformed, this must be multiplied by the Jacobian. This
##' transformation cannot contain unknown parameters. For example, with a log
##' transformation, \code{delta=1/y}. The jacobian is calculated automatically
##' for the transform option. Ignored if response has class, \code{response} or
##' \code{repeated}.
##' @param transform Transformation of the response variable: \code{identity},
##' \code{exp}, \code{square}, \code{sqrt}, or \code{log}.
##' @param link Link function for the mean: \code{identity}, \code{exp},
##' \code{square}, \code{sqrt}, or \code{log}.
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
##' @return A list of classes \code{kalseries} and \code{recursive} is
##' returned.
##' @author J.K. Lindsey
### @seealso \code{\link[growth]{carma}}, \code{\link[growth]{elliptic}},
### \code{\link[rmutil]{finterp}}, \code{\link[repeated]{gar}},
### \code{\link[rmutil]{gettvc}}, \code{\link[repeated]{gnlmm}},
### \code{\link[gnlm]{gnlr}}, \code{\link[rmutil]{iprofile}},
### \code{\link[repeated]{kalcount}}, \code{\link[event]{kalsurv}},
### \code{\link[rmutil]{mprofile}}, \code{\link[rmutil]{read.list}},
### \code{\link[rmutil]{restovec}}, \code{\link[rmutil]{rmna}},
### \code{\link[rmutil]{tcctomat}}, \code{\link[rmutil]{tvctomat}}.
##' @keywords models
##' @examples
##' 
##' treat <- c(0,0,1,1)
##' tr <- tcctomat(treat)
##' dose <- matrix(rpois(20,10), ncol=5)
##' dd <- tvctomat(dose)
##' y <- restovec(matrix(rnorm(20), ncol=5), name="y")
##' reps <- rmna(y, ccov=tr, tvcov=dd)
##' #
##' # normal intensity, independence model
##' kalseries(y, intensity="normal", dep="independence", preg=1, pshape=5)
##' \dontrun{
##' # random effect
##' kalseries(y, intensity="normal", dep="frailty", preg=1, pinitial=1, psh=5)
##' # serial dependence
##' kalseries(y, intensity="normal", dep="serial", preg=1, pinitial=1,
##' 	pdep=0.1, psh=5)
##' # random effect and autoregression
##' kalseries(y, intensity="normal", dep="frailty", preg=1, pinitial=1,
##' 	pdep=0.1, psh=5)
##' #
##' # add time-constant variable
##' kalseries(y, intensity="normal", dep="serial", pinitial=1,
##' 	pdep=0.1, psh=5, preg=c(1,0), ccov=treat)
##' # or equivalently
##' kalseries(y, intensity="normal", mu=~treat, dep="serial", pinitial=1,
##' 	pdep=0.1, psh=5, preg=c(1,0))
##' # or
##' kalseries(y, intensity="normal", mu=~b0+b1*treat, dep="serial",
##' 	pinitial=1, pdep=0.1, psh=5, preg=c(1,0), envir=reps)
##' #
##' # add time-varying variable
##' kalseries(y, intensity="normal", dep="serial", pinitial=1, pdep=0.1,
##' 	psh=5, preg=c(1,0), ccov=treat, ptvc=0, tvc=dose)
##' # or equivalently, from the environment
##' dosev <- as.vector(t(dose))
##' kalseries(y, intensity="normal",
##' 	mu=~b0+b1*rep(treat,rep(5,4))+b2*dosev,
##' 	dep="serial", pinitial=1, pdep=0.1, psh=5, ptvc=c(1,0,0))
##' # or from the reps data object
##' kalseries(y, intensity="normal", mu=~b0+b1*treat+b2*dose,
##' 	dep="serial", pinitial=1, pdep=0.1, psh=5, ptvc=c(1,0,0),
##' 	envir=reps)
##' # try power variance family instead of gamma distribution for mixture
##' kalseries(y, intensity="normal", mu=~b0+b1*treat+b2*dose,
##' 	dep="serial", pinitial=1, pdep=0.1, psh=5, ptvc=c(1,0,0),
##' 	pfamily=0.1, envir=reps)
##' # first-order one-compartment model
##' # data objects for formulae
##' dose <- c(2,5)
##' dd <- tcctomat(dose)
##' times <- matrix(rep(1:20,2), nrow=2, byrow=TRUE)
##' tt <- tvctomat(times)
##' # vector covariates for functions
##' dose <- c(rep(2,20),rep(5,20))
##' times <- rep(1:20,2)
##' # functions
##' mu <- function(p) exp(p[1]-p[3])*(dose/(exp(p[1])-exp(p[2]))*
##' 	(exp(-exp(p[2])*times)-exp(-exp(p[1])*times)))
##' shape <- function(p) exp(p[1]-p[2])*times*dose*exp(-exp(p[1])*times)
##' # response
##' conc <- matrix(rgamma(40,shape(log(c(0.01,1))),
##' 	scale=mu(log(c(1,0.3,0.2))))/shape(log(c(0.1,0.4))),ncol=20,byrow=TRUE)
##' conc[,2:20] <- conc[,2:20]+0.5*(conc[,1:19]-matrix(mu(log(c(1,0.3,0.2))),
##' 	ncol=20,byrow=TRUE)[,1:19])
##' conc <- restovec(ifelse(conc>0,conc,0.01))
##' reps <- rmna(conc, ccov=dd, tvcov=tt)
##' #
##' # constant shape parameter
##' kalseries(reps, intensity="gamma", dep="independence", mu=mu,
##' 	ptvc=c(-1,-1.1,-1), pshape=1.5)
##' # or
##' kalseries(reps, intensity="gamma", dep="independence",
##' 	mu=~exp(absorption-volume)*
##' 	dose/(exp(absorption)-exp(elimination))*
##' 	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
##' 	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
##' 	pshape=1.2)
##' # add serial dependence
##' kalseries(reps, intensity="gamma", dep="serial", pdep=0.9,
##' 	mu=~exp(absorption-volume)*
##' 	dose/(exp(absorption)-exp(elimination))*
##' 	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
##' 	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
##' 	pshape=0.2)
##' # time dependent shape parameter
##' kalseries(reps, intensity="gamma", dep="independence", mu=mu,
##' 	shape=shape, ptvc=c(-1,-1.1,-1), pshape=c(-3,0))
##' # or
##' kalseries(reps, intensity="gamma", dep="independence",
##' 	mu=~exp(absorption-volume)*
##' 	dose/(exp(absorption)-exp(elimination))*
##' 	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
##' 	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
##' 	shape=~exp(b1-b2)*times*dose*exp(-exp(b1)*times),
##' 	pshape=list(b1=-3,b2=0))
##' # add serial dependence
##' kalseries(reps, intensity="gamma", dep="serial", pdep=0.5,
##' 	mu=~exp(absorption-volume)*
##' 	dose/(exp(absorption)-exp(elimination))*
##' 	(exp(-exp(elimination)*times)-exp(-exp(absorption)*times)),
##' 	ptvc=list(absorption=-1,elimination=-1.1,volume=-1),
##' 	shape=~exp(b1-b2)*times*dose*exp(-exp(b1)*times),
##' 	pshape=list(b1=-3,b2=0))
##' }
##' @aliases kalseries deviance.kalseries residuals.kalseries fitted.kalseries print.kalseries 
##' @export kalseries
kalseries <- function(response=NULL, times=NULL, intensity="exponential",
	depend="independence", mu=NULL, shape=NULL, density=FALSE, ccov=NULL,
	tvcov=NULL, torder=0, interaction=NULL, preg=NULL, ptvc=NULL,
	pintercept=NULL, pshape=NULL, pinitial=1, pdepend=NULL,
	pfamily=NULL, delta=NULL, transform="identity", link="identity",
	envir=parent.frame(), print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, fscale=1, iterlim=100, typsize=abs(p),
	stepmax=10*sqrt(p%*%p)){
#
# Pareto (beta) likelihood function for serial dependence
#
series <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("kserie_c",
		p=as.double(p),
		y=as.double(y),
		t=as.double(times),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		npv=as.integer(npv),
		model=as.integer(mdl),
		link=as.integer(lnk),
		density=as.integer(density),
		pfamily=as.integer(!is.null(pfamily)),
		dep=as.integer(dep),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bbb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")
	z$like}
#
# Pareto (beta) likelihood function for frailty
#
serief <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("krand_c",
		p=as.double(p),
		y=as.double(y),
		t=as.double(times),
		x=as.double(resp$ccov$ccov),
		nind=as.integer(nind),
		nobs=as.integer(nobs),
		nbs=as.integer(n),
		nccov=as.integer(nccov),
		npv=as.integer(npv),
		model=as.integer(mdl),
		link=as.integer(lnk),
		density=as.integer(density),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(n),
		rpred=double(n),
		rf=as.integer(rf),
		bbb=as.double(b),
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
tmp <- c("exponential","Weibull","gamma","gen logistic","normal",
	"logistic","Cauchy","Laplace","log normal","log logistic",
	"log Cauchy","log Laplace","inverse Gauss")
mdl <- match(intensity <- match.arg(intensity,tmp),tmp)
tmp <- c("independence","Markov","serial","frailty")
dep <- match(depend <- match.arg(depend,tmp),tmp)-1
transform <- match.arg(transform,c("identity","exp","square","sqrt","log"))
tmp <- c("identity","exp","square","sqrt","log")
lnk <- match(link <- match.arg(link,tmp),tmp)
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
		stop("kalseries only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("kalseries does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)deparse(substitute(response))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# if envir, remove extra (multivariate) responses
#
if(!respenv&&inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("kalseries does not handle data with NAs")
	cn <- deparse(substitute(response))
	if(length(grep("\"",cn))>0)cn <- response
	if(length(cn)>1)stop("only one response variable allowed")
	response <- envir
	col <- match(cn,colnames(response$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- response$response$type[col]
	if(dim(response$response$y)[2]>1){
		response$response$y <- response$response$y[,col,drop=FALSE]
		if(!is.null(response$response$delta)){
			response$response$delta <- response$response$delta[,col,drop=FALSE]
			if(all(response$response$delta==1)||all(is.na(response$response$delta)))response$response$delta <- NULL}}}
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
	else resp <- restovec(response,times,delta=delta)
	# time-constant covariates
	if(is.null(ccov))nccov <- 0
	else {
		if(!inherits(ccov,"tccov")){
			ccname <- deparse(substitute(ccov))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- deparse(substitute(ccov))
				if(dim(ccov)[2]>1){
					tmp <- NULL
					for(i in 1:dim(ccov)[2])
						tmp <- c(tmp,paste(ccname,i,sep=""))
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
					for(i in 1:dim(tvcov[[1]])[2])
						tmp <- c(tmp,paste(tvcname,i,sep=""))
					tvcname <- tmp}
				else tvcname <- colnames(tvcov[[1]])}
			tvcov <- tvctomat(tvcov,names=tvcname)}
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
if((intensity=="exponential"||intensity=="Weibull"||intensity=="gamma"
	||intensity=="log normal"||intensity=="log logistic"
	||intensity=="log Cauchy"||intensity=="log Laplace"
	||intensity=="inverse Gauss")&&type!="unknown"&&type!="duration"
	&&type!="continuous")stop("duration data required")
if((intensity=="gen logistic"||intensity=="normal"||intensity=="logistic"
	||intensity=="Cauchy"||intensity=="Laplace")&&type!="unknown"
	&&type!="continuous"&&type!="duration")stop("continuous data required")
if(!is.null(resp$response$censor))stop("kalseries does not handle censoring")
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
	mu2 <- if(respenv)
		finterp(mu,.envir=response,.name=envname,.expand=is.null(preg))
		else finterp(mu,.envir=envir,.name=envname,.expand=is.null(preg))
	npt1 <- length(attr(mu2,"parameters"))
	if(is.character(attr(mu2,"model"))){
	# W&R formula
		if(length(attr(mu2,"model"))==1){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
	# formula with unknowns
		if(npr!=npt1){
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
		else {
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
			attributes(fnenvir(sh1,.envir=envir))}
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
tvc <- length(ptvc)
if(!rf&&(ttvc>0&&tvc!=ttvc||ttvc==0&&tvc>0))
	stop(paste(ttvc,"initial estimates of coefficients for time-varying covariates must be supplied"))
if(intensity=="exponential"){
	sf <- FALSE
	pshape <- NULL}
else {
	if(is.null(pshape))
		stop("Initial value of the shape parameter must be supplied")
	if(!sf){
		if(pshape<=0)stop("Shape must be positive")
		else pshape <- log(pshape)}}
if(intensity=="gen logistic"){
	if(is.null(pintercept))
		stop("Initial value of the intercept parameter must be supplied")}
else pintercept <- NULL
if(pinitial<=0)stop("Estimate of initial parameter must greater than 0")
else pinitial <- log(pinitial)
#
# check what dependence is required
#
frser <- FALSE
if(depend=="independence")pdepend <- NULL
else if(depend=="frailty"){
	frser <- !is.null(pdepend)
	if(frser){
		if(pdepend<=0)stop("Dependence parameter must be positive")
		pdepend <- log(pdepend)}}
else {
	if(is.null(pdepend))
		stop("An estimate of the dependence parameter must be supplied")
	else if(pdepend<=0||pdepend>=1)
		stop("Dependence parameter must be between 0 and 1")
	else pdepend <- log(pdepend/(1-pdepend))}
#
# check times
#
if(is.null(resp$response$times)){
	if(depend=="serial"||depend=="Markov")
		stop("No times. Serial and Markov dependence cannot be fitted.")
	ave <- times <- 0}
else {
	if(torder>length(unique(resp$response$times)))
		stop("torder is too large for the number of distinct times")
	ave <- mean(resp$response$times)
	times <- resp$response$times-ave}
#
# check interactions with time
#
if(!is.null(interaction)){
	if(length(interaction)!=nccov)
		stop(paste(nccov,"interactions with time must be specified"))
	else if(any(interaction>torder))
		stop(paste("Interactions can be at most of order ",torder))}
else interaction <- rep(0,nccov)
#
# check number of parameter estimates
#
if(!is.null(pfamily)&&depend=="frailty")
		stop("pfamily cannot be used with frailty dependence")
npv <- torder+sum(interaction)
if(rf&&npreg>0)nccov <- npreg-1
if(!rf&&1+nccov+npv!=npreg)
	stop(paste(1+nccov+npv,"regression estimates must be supplied"))
y <- resp$response$y
#
# check Jacobian and transformations
if(!is.null(resp$response$delta))jacob <- if(length(resp$response$delta)==1)
		-n*log(resp$response$delta)
 	else -sum(log(resp$response$delta))
else jacob <- 0
if((mdl<=3||mdl>=8)&&any(y<=0))stop("All responses must be positive")
if(transform=="exp"){
	jacob <- jacob-sum(y)
	y <- exp(y)}
else if(transform!="identity"){
	if(any(y==0))stop("Zero response values: invalid transformation")
	if(transform=="square"){
		jacob <- jacob-sum(log(abs(y)))-n*log(2)
		y  <- y^2}
	if(any(y<0))stop("Nonpositive response values: invalid transformation")
	if(transform=="sqrt"){
		jacob <- jacob+sum(log(y))/2+n*log(2)
		y <- sqrt(y)}
	else if(transform=="log"){
		jacob <- jacob+sum(log(y))
		y <- log(y)}}
#
# set up for location function and make sure it gives correct number
# of values
#
if(rf){
	if(tvc>0&&nccov>0)
		stop("With a mean function or formula, initial estimates must be supplied either in preg or in ptvc")
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
np <- 1+nccov+npv+tvc+1+(depend=="serial"||depend=="Markov")+
	nps+(!is.null(pintercept))+frser+(!is.null(pfamily))
nps1 <- np-nps+1
p <- c(preg,ptvc,pinitial,pdepend,pfamily,pintercept,pshape)
if(dep==3)serie <- serief
else serie <- series
if(fscale==1)fscale <- serie(p)
if(is.na(serie(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(serie, p=p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
p <- z0$estimate
like <- z0$minimum+jacob
#
# calculate se's
#
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
#
# calculate recursive fitted values
#
if(mdl==4)z <- list()
else {
	z <- if(depend=="frailty"){
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		.C("krand_c",
			p=as.double(p),
			y=as.double(y),
			t=as.double(times),
			x=as.double(resp$ccov$ccov),
			nind=as.integer(nind),
			nobs=as.integer(nobs),
			nbs=as.integer(n),
			nccov=as.integer(nccov),
			npv=as.integer(npv),
			model=as.integer(mdl),
			link=as.integer(lnk),
			density=as.integer(density),
			torder=as.integer(torder),
			inter=as.integer(interaction),
			tvc=as.integer(tvc),
			tvcov=as.double(resp$tvcov$tvcov),
			fit=as.integer(1),
			pred=double(n),
			rpred=double(n),
			rf=as.integer(rf),
			bbb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			frser=as.integer(frser),
			like=double(1),
			#DUP=FALSE,
			PACKAGE="repeated")}
	else {
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		z <- .C("kserie_c",
			p=as.double(p),
			y=as.double(y),
			t=as.double(times),
			x=as.double(resp$ccov$ccov),
			nind=as.integer(nind),
			nobs=as.integer(nobs),
			nbs=as.integer(n),
			nccov=as.integer(nccov),
			npv=as.integer(npv),
			model=as.integer(mdl),
			link=as.integer(lnk),
			density=as.integer(density),
			pfamily=as.integer(!is.null(pfamily)),
			dep=as.integer(dep),
			torder=as.integer(torder),
			inter=as.integer(interaction),
			tvc=as.integer(tvc),
			tvcov=as.double(resp$tvcov$tvcov),
			fit=as.integer(1),
			pred=double(n),
			rpred=double(n),
			rf=as.integer(rf),
			bbb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			like=double(1),
			#DUP=FALSE,
			PACKAGE="repeated")}}
if(rf&&tvc>0){
	nccov <- tvc
	tvc <- 0}
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
	npr=1+nccov+tvc+torder+sum(interaction),
	shape=sh1,
	nps=nps,
	density=density,
	depend=depend,
	torder=torder,
	interaction=interaction,
	response=resp$response,
	pred=z$pred,
	rpred=z$rpred,
	transform=transform,
	ccov=resp$ccov,
	tvcov=resp$tvcov,
	link=link,
	maxlike=like,
	aic=like+np,
	df=n-np,
	npt=np,
	npv=npv,
	coefficients=p,
	se=se,
	cov=cov,
	corr=corr,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- if(mdl==4)"kalseries" else c("kalseries","recursive")
return(z)}

### standard methods
###
##' @export
deviance.kalseries <- function(object, ...) 2*object$maxlike
##' @export
fitted.kalseries <- function(object, recursive=TRUE, ...)
if(recursive) object$rpred else object$pred
##' @export
residuals.kalseries <- function(object, recursive=TRUE, ...){
if(object$transform=="exp")object$response$y <- exp(object$response$y)
else if(object$transform=="square")object$response$y  <- object$response$y^2
else if(object$transform=="sqrt")object$response$y <- sqrt(object$response$y)
else if(object$transform=="log")object$response$y <- log(object$response$y)
res <- if(recursive)object$response$y-object$rpred else object$response$y-object$pred
class(res) <- "residuals"
res}

### print method
###
##' @export
print.kalseries <- function(x,digits=max(3,.Options$digits-3),correlation=TRUE, ...){
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
if(!is.null(z$response$time))
	cat("Times centred at      ",mean(z$response$time),"\n")
cat("Transformation        ",z$trans,"\n")
cat("Link function         ",z$link,"\n\n")
if(z$density)cat(z$intensity," density",sep="")
else cat(z$intensity," intensity",sep="")
if(z$depend=="independence")cat(" with independence\n")
else cat(" with",z$depend,"dependence",if(z$frser)"and AR","\n")
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
	if(z$torder){
		cname <- c(cname,paste("t^",1:z$torder,sep=""))
		if(length(z$interaction)>0){
			for(i in 1:nccov)if(z$interaction[i]>0){
				cname <- c(cname,paste(colnames(z$ccov$ccov)[i],".t^",1:z$interaction[i],sep=""))}}}
	if(!is.null(z$tvcov))cname <- c(cname,colnames(z$tvcov$tvcov))}
dimnames(coef.table) <- list(cname, c("estimate","se"))
print.default(coef.table, digits=digits, print.gap=2)
if(is.function(z$shape))cat("\nDependence parameters\n")
else cat("\nNonlinear parameters\n")
coef <- exp(z$coef[(npt-deppar-z$frser-z$pfamily-glm-expm):npt])
cname <- "initial"
if(deppar){
	cname <- c(cname,"depend")
	coef[2] <- coef[2]/(1+coef[2])}
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
if(z$depend=="frailty"){
	tmp <- trigamma(exp(-z$coef[npt-deppar-expm-z$frser]))
	cat("Correlation =",tmp/(tmp+trigamma(1)),"\n")}
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
