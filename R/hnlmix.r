#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 2002 J.K. Lindsey
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
#     hnlmix(y=NULL, distribution="normal", mixture="normal",
#	random=NULL, nest=NULL, mu=NULL, shape=NULL, linear=NULL,
#	pmu=NULL, pshape=NULL, pmix=NULL, prandom=NULL, delta=1, common=FALSE,
#	envir=parent.frame(), print.level=0, typsize=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
#	iterlim=100, fscale=1, eps=1.0e-4)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with one arbitrary
#  random parameter using h-likelihood



##' Generalized Nonlinear Regression using h-likelihood for a Random Parameter
##' 
##' \code{hnlmix} fits user-specified nonlinear regression equations to one or
##' both parameters of the common one and two parameter distributions. One
##' parameter of the location regression is random with some specified mixing
##' distribution.
##' 
##' It is recommended that initial estimates for \code{pmu} and \code{pshape}
##' be obtained from \code{gnlr}.
##' 
##' These nonlinear regression models must be supplied as formulae where
##' parameters are unknowns. (See \code{\link[rmutil]{finterp}}.)
##' 
##' 
##' @param y A response vector of uncensored data, a two column matrix for
##' binomial data or censored data, with the second column being the censoring
##' indicator (1: uncensored, 0: right censored, -1: left censored), or an
##' object of class, \code{response} (created by
##' \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
##' \code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
##' \code{repeated} data object contains more than one response variable, give
##' that object in \code{envir} and give the name of the response variable to
##' be used here.
##' @param distribution The distribution for the response: binomial, beta
##' binomial, double binomial, mult(iplicative) binomial, Poisson, negative
##' binomial, double Poisson, mult(iplicative) Poisson, gamma count, Consul
##' generalized Poisson, logarithmic series, geometric, normal, inverse Gauss,
##' logistic, exponential, gamma, Weibull, extreme value, Cauchy, Pareto,
##' Laplace, Levy, beta, simplex, or two-sided power. (For definitions of
##' distributions, see the corresponding [dpqr]distribution help.)
##' @param mixture The mixing distribution for the random parameter (whose
##' initial values are supplied in \code{prandom}): normal, logistic, inverse
##' Gauss, gamma, inverse gamma, Weibull, or beta. The first two have zero
##' location parameter, the next three have unit location parameter, and the
##' last one has location parameter set to 0.5.
##' @param random The name of the random parameter in the \code{mu} formula.
##' @param nest The cluster variable classifying observations by the unit upon
##' which they were observed. Ignored if \code{y} or \code{envir} has class,
##' \code{response} or \code{repeated}.
##' @param mu A user-specified formula containing named unknown parameters,
##' giving the regression equation for the location parameter. This may contain
##' the keyword, \code{linear} referring to a linear part.
##' @param shape A user-specified formula containing named unknown parameters,
##' giving the regression equation for the shape parameter. This may contain
##' the keyword, \code{linear} referring to a linear part. If nothing is
##' supplied, this parameter is taken to be constant. This parameter is the
##' logarithm of the usual one.
##' @param linear A formula beginning with ~ in W&R notation, specifying the
##' linear part of the regression function for the location parameter or list
##' of two such expressions for the location and/or shape parameters.
##' @param pmu Vector of initial estimates for the location parameters. These
##' must be supplied either in their order of appearance in the formula or in a
##' named list.
##' @param pshape Vector of initial estimates for the shape parameters. These
##' must be supplied either in their order of appearance in the expression or
##' in a named list.
##' @param pmix If NULL, this parameter is estimated from the variances. If a
##' value is given, it is taken as fixed.
##' @param prandom Either one estimate of the random effects or one for each
##' cluster (see \code{nest}), in which case the last value is not used. If the
##' location parameter of the mixing distribution is zero, the last value is
##' recalculated so that their sum is zero; if it is unity, they must all be
##' positive and the last value is recalculated so that the sum of their
##' logarithms is zero; if it is 0.5, they must all lie in (0,1) and the last
##' value is recalculated so that the sum of their logits is zero.
##' @param delta Scalar or vector giving the unit of measurement (always one
##' for discrete data) for each response value, set to unity by default. For
##' example, if a response is measured to two decimals, \code{delta=0.01}. If
##' the response is transformed, this must be multiplied by the Jacobian. The
##' transformation cannot contain unknown parameters. For example, with a log
##' transformation, \code{delta=1/y}. (The delta values for the censored
##' response are ignored.)
##' @param common If TRUE, the formulae with unknowns for the location and
##' shape have names in common. All parameter estimates must be supplied in
##' \code{pmu}.
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
##' @param eps Arguments for nlm.
##' @param points Arguments for nlm.
##' @return A list of class \code{hnlmix} is returned that contains all of the
##' relevant information calculated, including error codes.
##' 
##' The two variances and shrinkage estimates of the random effects are
##' provided.
##' @author J.K. Lindsey
### @seealso \code{\link[growth]{carma}}, \code{\link[rmutil]{finterp}},
### \code{\link[growth]{elliptic}}, \code{\link[repeated]{glmm}},
### \code{\link[repeated]{gnlmix}}, \code{\link[repeated]{gnlmm}},
### \code{\link[gnlm]{gnlr}}, \code{\link[repeated]{kalseries}},
### \code{\link[gnlm]{nlr}}, \code{\link[stats]{nls}}.
##' @keywords models
##' @aliases hnlmix deviance.hnlm residuals.hnlm fitted.hnlm print.hnlm 
##' @examples
##' 
##' dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
##' #y <- rgamma(20,2+0.3*dose,scale=2)+rep(rnorm(4,0,4),rep(5,4))
##' y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
##'        1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
##'       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
##'       20.972267, 17.178012)
##' resp <- restovec(matrix(y, nrow=4, byrow=TRUE), name="y")
##' reps <- rmna(resp, tvcov=tvctomat(matrix(dose, nrow=4, byrow=TRUE), name="dose"))
##' 
##' # same linear normal model with random normal intercept fitted four ways
##' # compare with growth::elliptic(reps, model=~dose, preg=c(0,0.6), pre=4)
##' glmm(y~dose, nest=individuals, data=reps)
##' gnlmm(reps, mu=~dose, pmu=c(8.7,0.25), psh=3.5, psd=3)
##' hnlmix(reps, mu=~a+b*dose+rand, random="rand", pmu=c(8.7,0.25),
##' 	pshape=3.44, prandom=0)
##' 
##' # gamma model with log link and random normal intercept fitted three ways
##' glmm(y~dose, family=Gamma(link=log), nest=individuals, data=reps, points=8)
##' gnlmm(reps, distribution="gamma", mu=~exp(a+b*dose), pmu=c(2,0.03),
##' 	psh=1, psd=0.3)
##' hnlmix(reps, distribution="gamma", mu=~exp(a+b*dose+rand), random="rand",
##' 	pmu=c(2,0.04), pshape=1, prandom=0)
##' 
##' # gamma model with log link and random gamma mixtures
##' hnlmix(reps, distribution="gamma", mixture="gamma",
##' 	mu=~exp(a*rand+b*dose), random="rand", pmu=c(2,0.04),
##' 	pshape=1.24, prandom=1)
##' hnlmix(reps, distribution="gamma", mixture="gamma",
##' 	mu=~exp(a+b*dose)*rand, random="rand", pmu=c(2,0.04),
##' 	pshape=1.24, prandom=1)
##' 
##' @export hnlmix
hnlmix <- function(y=NULL, distribution="normal", mixture="normal",
	random=NULL, nest=NULL, mu=NULL, shape=NULL, linear=NULL,
	pmu=NULL, pshape=NULL, pmix=NULL, prandom=NULL, delta=1, common=FALSE,
	envir=parent.frame(), print.level=0, typsize=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
	iterlim=100, fscale=1, eps=1.0e-4, points=5){
#
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
bindata <- distribution=="binomial"||distribution=="double binomial"||
	distribution=="beta binomial"||distribution=="mult binomial"
#
# check mixture
#
#mixture <- match.arg(mixture,c("normal","logistic","Cauchy","Laplace",
#	"gamma","inverse Gauss","Weibull","beta"))
mixture <- match.arg(mixture,c("normal","logistic","gamma","inverse gamma",
	"inverse Gauss","Weibull","beta"))
mean0 <- if(mixture=="normal"||mixture=="logistic"||mixture=="Cauchy"||
	 mixture=="Laplace")0 else if(mixture=="beta")2 else 1
fixed <- !is.null(pmix)
if(fixed&&pmix<0)stop("pmix must be positive")
#
# check random parameters
#
if(is.null(random))stop("name of random parameter must be supplied")
if(!is.character(random))stop("random must be the name of a parameter")
if(length(random)>1)stop("only one random parameter allowed")
#
# check for parameters common to location and shape functions
#
if(common&&!is.null(linear))
	stop("linear cannot be used with common parameters")
#
# count number of parameters
#
npl <- length(pmu)
nps <- length(pshape)
np1 <- npl+nps
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("hnlmix only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("hnlmix does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# check that formulae with unknowns are supplied
#
if(!inherits(mu,"formula"))stop("mu must be a formula")
if(shp&&!is.null(shape)&&!inherits(shape,"formula"))
	stop("shape must be a formula")
#
# find linear part of each regression
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
# if a data object was supplied, modify formulae or functions to read from it
#
lin1a <- lin2a <- mu2 <- sh2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	# modify formulae
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname,.args=random)
			else finterp(mu,.envir=envir,.name=envname,.args=random)}
	if(inherits(shape,"formula")){
		sh2 <- if(respenv)finterp(shape,.envir=y,.name=envname)
			else finterp(shape,.envir=envir,.name=envname)}
	if(inherits(lin1,"formula")){
		lin1a <- if(respenv)finterp(lin1,.envir=y,.name=envname)
			else finterp(lin1,.envir=envir,.name=envname)}
	if(inherits(lin2,"formula")){
		lin2a <- if(respenv)finterp(lin2,.envir=y,.name=envname)
			else finterp(lin2,.envir=envir,.name=envname)}
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
# check if linear contains W&R formula
#
if(inherits(lin1,"formula")){
	tmp <- attributes(if(respenv)finterp(lin1,.envir=y,.name=envname)
		else finterp(lin1,.envir=envir,.name=envname))
	lf1 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	else if(length(tmp$model)==1)stop("linear must contain covariates")
	rm(tmp)}
else lf1 <- 0
if(inherits(lin2,"formula")){
	tmp <- attributes(if(respenv)finterp(lin2,.envir=y,.name=envname)
		else finterp(lin2,.envir=envir,.name=envname))
	lf2 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	else if(length(tmp$model)==1)stop("linear must contain covariates")
	rm(tmp)}
else lf2 <- 0
#
# transform location formula to function and check number of parameters
#
if(lf1>0)random <- c(random,"linear")
mu3 <- if(respenv)finterp(mu,.envir=y,.name=envname,.args=random)
	else finterp(mu,.envir=envir,.name=envname,.args=random)
npt1 <- length(attr(mu3,"parameters"))
if(is.character(attr(mu3,"model")))stop("mu cannot be a W&R formula")
if(npl!=npt1&&!common&&lf1==0){
	cat("\nParameters are ")
	cat(attr(mu3,"parameters"),"\n")
	stop(paste("pmu should have",npt1,"estimates"))}
if(is.list(pmu)){
	if(!is.null(names(pmu))){
		o <- match(attr(mu3,"parameters"),names(pmu))
		pmu <- unlist(pmu)[o]
		if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
	else pmu <- unlist(pmu)}
#
# if linear part, modify location function appropriately
#
if(lf1>0){
	dm1 <- if(respenv)wr(lin1,data=y)$design
		else wr(lin1,data=envir)$design
	if(is.null(mu2))mu2 <- mu3
	mu1 <- function(p,random)mu3(p,random,dm1%*%p[(npt1+1):(npt1+lf1)])}
else {
	mu1 <- mu3
	rm(mu3)}
#
# check that correct number of estimates was supplied
#
nlp <- npt1+lf1
if(!common&&nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
npl1 <- if(common&&!inherits(shape,"formula")) 1 else nlp+1
#
# transform shape formula to function and check number of parameters
#
sh3 <- NULL
if(inherits(shape,"formula")){
        old <- if(common)mu1 else NULL
        mufn <- if(lf2>0) "linear" else NULL
        sh3 <- if(respenv)finterp(shape,.envir=y,.start=npl1,.name=envname,.old=old,.args=mufn)
        	else finterp(shape,.envir=envir,.start=npl1,.name=envname,.old=old,.args=mufn)
        npt2 <- length(attr(sh3,"parameters"))
        if(is.character(attr(sh3,"model")))stop("shape cannot be a W&R formula")
        if(nps!=npt2&&!common&&lf2==0){
        	cat("\nParameters are ")
        	cat(attr(sh3,"parameters"),"\n")
        	stop(paste("pshape should have",npt2,"estimates"))}
        if(is.list(pshape)){
        	if(!is.null(names(pshape))){
        		o <- match(attr(sh3,"parameters"),names(pshape))
        		pshape <- unlist(pshape)[o]
        		if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
        	else pshape <- unlist(pshape)}}
else if(is.null(shape)&&shp){
	sh3 <- function(p) p[npl1]*rep(1,n)
	sh2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
#
# if linear part, modify shape function appropriately
#
if(lf2>0){
	dm2 <- if(respenv)wr(lin2,data=y)$design
		else wr(lin2,data=envir)$design
	if(is.null(sh2))sh2 <- sh3
	sh1 <- sh3(p,dm2%*%p[(npl1+lf2-1):np])}
else {
	sh1 <- sh3
	rm(sh3)}
#
# if distribution has a shape parameter, check that correct number of
# estimates was supplied
#
if(shp){
	nlp <- npt2+lf2
	if(!common&&nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))}
#
# when there are parameters common to location and shape functions,
# check that correct number of estimates was supplied
#
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"))))
	if(nlp!=npl)stop(paste("with a common parameter model, pmu should contain",nlp,"estimates"))}
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(!is.null(y$response$delta))
		delta <- as.vector(y$response$delta)
	nest <- covind(y)
	type <- y$response$type
	envir <- y$response
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("hnlmix does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	nest <- covind(envir)
	type <- envir$response$type[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$n)&&!all(is.na(envir$response$n[,col])))
		y <- cbind(y,envir$response$n[,col]-y)
	else if(!is.null(envir$response$censor)&&!all(is.na(envir$response$censor[,col])))
		y <- cbind(y,envir$response$censor[,col])
	if(!is.null(envir$response$delta))
		delta <- as.vector(envir$response$delta[,col])
	envir <- envir$response
	envir$y <- envir$y[,col,drop=FALSE]
	if(!is.null(envir$n))envir$n <- envir$n[,col,drop=FALSE]
	else if(!is.null(envir$censor))envir$censor <- envir$censor[,col,drop=FALSE]
	if(!is.null(envir$delta))envir$delta <- envir$delta[,col,drop=FALSE]}
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("hnlmix only handles univariate responses")
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	nest <- covind(y)
	type <- y$type
	envir <- y
	y <- response(y)}
unest <- unique(nest)
nnest <- length(unest)
unest <- unest[-nnest]
if(is.null(nest)||nnest==1)stop("appropriate nest indicator required")
if(is.null(prandom))
	stop(paste("one or",nnest,"values must be supplied for prandom"))
if(length(prandom)==1)prandom <- rep(prandom,nnest)
else if(length(prandom)!=nnest)
	stop(paste(nnest,"values must be supplied for prandom"))
if(mean0==1){
	if(any(prandom<=0))stop("prandom must all be positive")
	prandom <- log(prandom)}
else if(mean0==2){
	if(any(prandom<=0|prandom>=1))stop("prandom must all be in (0,1)")
	prandom <- log(prandom/(1-prandom))}
p <- c(pmu,pshape,prandom[-nnest])
np <- length(p)
if(any(is.na(y)))stop("NAs in y - use rmna")
#
# check that data are appropriate for distribution
#
if(bindata){
	# binomial data
	if(type!="unknown"&&type!="nominal")stop("nominal data required")
	if((is.vector(y)||(length(dim(y))==2&&dim(y)[2]==1))&&all(y==0|y==1))
		y <- cbind(y,1-y)
	if(length(dim(y))!=2||dim(y)[2]!=2)
		stop(paste("Two column matrix required for response: successes and failures"))
	if(any(y<0))stop("All response values must be positive")
	n <- dim(y)[1]
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
		my <- if(censor)3*max(y[,1]) else 3*max(y)
	n <- if(length(dim(y))==2)dim(y)[1] else length(y)}
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
	if(any(y<1))stop("All response values must be integers > 0")}
else if(distribution=="beta"||distribution=="simplex"||
	distribution=="two-sided power"){
	if(type!="unknown"&&type!="continuous")stop("continuous data required")
	if(any(y<=0)||any(y>=1))
		stop("All response values must lie between 0 and 1")}
else if(!bindata&&type!="unknown"&&type!="continuous"&&type!="duration")
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
	lc <- ifelse(y[,2]==-1,0,1)}
else cc <- 1
#
# prepare weights and unit of measurement
#
wt <- rep(1,n)
if(length(delta)==1)delta <- rep(delta,n)
else if(length(delta)!=n)
	stop("delta must be the same length as the other variables")
delta2 <- mean(delta)
#
# check that location function returns appropriate values
#
if(any(is.na(mu1(pmu,0))))stop("The location model returns NAs: probably invalid initial values")
if(distribution=="Levy"&&any(y<=mu1(p)))
	stop("location parameter must be strictly less than corresponding observation")
#
# check that shape function returns appropriate values
#
if(shp&&any(is.na(sh1(p))))
	stop("The shape model returns NAs: probably invalid initial values")
if(distribution=="Pareto"&&exp(sh1(p))<=1)stop("shape parameters must be > 0")
#
# set up distribution
#
if(!censor)fcn <- switch(distribution,
	binomial=function(m,p) dbinom(y[,1],nn,m,TRUE),
	"beta binomial"=function(m,p){
		s <- exp(sh1(p))
		t <- s*m
		u <- s*(1-m)
		lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u)+lchoose(nn,y[,1])},
	"double binomial"=function(m,p)
		.C("ddb",as.integer(y[,1]),as.integer(nn),
			as.double(m),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res,
	"mult binomial"=function(m,p)
		.C("dmb",as.integer(y[,1]),as.integer(nn),
			as.double(m),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res,
	Poisson=function(m,p)dpois(y,m,TRUE),
	"negative binomial"=function(m,p)dnbinom(y,exp(sh1(p)),mu=m,log=TRUE),
	"double Poisson"=function(m,p)
		.C("ddp",as.integer(y),as.integer(my),as.double(m),
			as.double(exp(sh1(p))),as.integer(n),as.double(wt),
			res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res,
	"mult Poisson"=function(m,p)
		.C("dmp",as.integer(y),as.integer(my),
			as.double(m),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res,
	"gamma count"=function(m,p){
		s <- exp(sh1(p))
		u <- m*s
		ifelse(y==0,pgamma(u,(y+1)*s,1,lower.tail=FALSE,log.p=TRUE),
			log(pgamma(u,y*s+(y==0),1)-pgamma(u,(y+1)*s,1)))},
	Consul=function(m,p){
		s <- exp(sh1(p))
		log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-y*log(s)-
			lgamma(y+1)},
	logarithmic=function(m,p){
		m <- 1/(1+exp(-m))
		y*log(m)-log(y)-log(-log(1-m))},
	geometric=function(m,p)y*log(m)-(y+1)*log(1+m),
	normal=function(m,p)dnorm(y,m,exp(sh1(p)/2),TRUE),
	"inverse Gauss"=function(m,p){
		t <- sh1(p)
		s <- exp(t)
		-(t+(y-m)^2/(y*s*m^2)+log(2*pi*y^3))/2},
	logistic=function(m,p)dlogis(y,m,exp(sh1(p))*sqrt(3)/pi,TRUE),
	Cauchy=function(m,p)dcauchy(y,m,exp(sh1(p)/2),TRUE),
	Laplace=function(m,p){
		t <- sh1(p)
		s <- exp(t)
		-abs(y-m)/s-t-log(2)},
	Pareto=function(m,p){
		s <- exp(sh1(p))
		t <- 1/(m*(s-1))
		log(s*t)-(s+1)*log(1+y*t)},
	exponential=function(m,p)dexp(y,1/m,TRUE),
	gamma=function(m,p){
		s <- exp(sh1(p))
		dgamma(y,s,scale=m/s,log=TRUE)},
	Weibull=function(m,p)dweibull(y,exp(sh1(p)),m,TRUE),
	"extreme value"=function(m,p)y+dweibull(exp(y),exp(sh1(p)),m,TRUE),
	beta=function(m,p){
		s <- exp(sh1(p))
		m <- m*s
		s <- s-m
		dbeta(y,m,s,log=TRUE)},
	simplex=function(m,p){
		t <- sh1(p)
		s <- exp(t)
		-(((y-m)/(m*(1-m)))^2/(y*(1-y)*s)+t+3*log(y*(1-y))+
			log(2*pi))/2},
	"two-sided power"=function(m,p){
		t <- sh1(p)
		s <- exp(t)
		t+(s-1)*ifelse(y<m,log(y/m),log((1-y)/(1-m)))})
else fcn <- switch(distribution,
	Poisson=function(m,p)cc*dpois(y[,1],m,TRUE)+log(lc-rc*ppois(y[,1],m)),
	"negative binomial"=function(m,p){
		s <- exp(sh1(p))
		cc*dnbinom(y[,1],s,mu=m,log=TRUE)+log(lc-rc*pnbinom(y[,1],s,m))},
	geometric=function(m,p)
		cc*(y[,1]*log(m)-(y[,1]+1)*log(1+m))+
			log(lc-rc*pgeom(y[,1],1/(1+m))),
	normal=function(m,p){
		s <- exp(sh1(p)/2)
		cc*dnorm(y[,1],m,s,TRUE)+log(lc-rc*pnorm(y[,1],m,s))},
	"inverse Gauss"=function(m,p){
		s <- exp(sh1(p))
		v <- sqrt(s*y[,1]/2)
		-cc*(log(s)+(y[,1]-m)^2/(y[,1]*s*m^2)+log(2*pi*y[,1]^3))/2+
			log(lc-rc*(pnorm((y[,1]/m-1)/v)+
			exp(2/(m*s))*pnorm(-(y[,1]/m+1)/v)))},
	logistic=function(m,p){
		s <- exp(sh1(p))
		cc*dlogis(y[,1],m,s*sqrt(3)/pi,TRUE)+
			log(lc-rc*plogis(y[,1],m,s*sqrt(3)/pi))},
	Cauchy=function(m,p){
		s <- exp(sh1(p)/2)
		cc*dcauchy(y[,1],m,s,TRUE)+log(lc-rc*pcauchy(y[,1],m,s))},
	Laplace=function(m,p){
		v <- sh1(p)
		s <- exp(v)
		u <- abs(y[,1]-m)/s
		t <- exp(-u)/2
		-cc*(u+v+log(2))+log(lc-rc*(ifelse(u<0,t,1-t)))},
	Pareto=function(m,p){
		s <- exp(sh1(p))
		t <- 1/(m*(s-1))
		cc*(log(s*t)-(s+1)*log(1+y[,1]*t))+
			log(lc-rc*((1+y[,1]/(m*(s-1)))^(-s)))},
	exponential=function(m,p)
		cc*dexp(y[,1],1/m,TRUE)+log(lc-rc*pexp(y[,1],1/m)),
	gamma=function(m,p){
		s <- exp(sh1(p))
		t <- m/s
		cc*dgamma(y[,1],s,scale=t,log=TRUE)+
			log(lc-rc*pgamma(y[,1],s,scale=t))},
	Weibull=function(m,p){
		s <- exp(sh1(p))
		cc*dweibull(y[,1],s,m,TRUE)+log(lc-rc*pweibull(y[,1],s,m))},
	"extreme value"=function(m,p){
		s <- exp(sh1(p))
		yy <- exp(y[,1])
		cc*(y[,1]+dweibull(yy,s,m,TRUE))+log(lc-rc*pweibull(yy,s,m))})
#
# set up mixing distribution
#
mix <- switch(mixture,
        normal=function(r,ss)dnorm(r,0,sqrt(ss),log=TRUE),
	logistic=function(r,ss)dlogis(r,0,sqrt(3*ss)/pi,log=TRUE),
	Cauchy=function(r,ss){
#		s <- ?? # bruce edit s to ss in dcauchy
		dcauchy(r,0,ss,log=TRUE)},
        Laplace=function(r,ss){
#		s <- ?? # bruce edit s to ss in Laplace
		-abs(r)/ss-log(ss)-log(2)},
	gamma=function(r,ss)dgamma(r,1/ss,scale=ss,log=TRUE),
        "inverse gamma"=function(r,ss)dgamma(1/r,1/ss,scale=ss,log=TRUE)/r^2,
        "inverse Gauss"=function(r,ss)
		-(log(ss)+(r-1)^2/(r*ss)+log(2*pi*r^3))/2,
	Weibull=function(r,ss){
		fn <- function(z)gamma(1+2/z)-gamma(1+1/z)^2-ss
		s <- uniroot(fn,c(0.02,20))$root
		dweibull(r,s,1,log=TRUE)},
	beta=function(r,ss){
#		ss <- 0.25*ss/(1+ss) # variance <= 0.25
#		s  <- 0.125/ss-0.5
		s <- 0.5/ss
		dbeta(r,s,s,log=TRUE)})
#
# combine to create the appropriate likelihood function
#
like <- if(mean0==0)function(p){
		r <- c(p[np1+unest],-sum(p[np1+unest]))
		ss <- if(!fixed)sum(r^2)/nnest else pmix
		m <- mu1(p,r[nest])
		var0 <- if(censor)sum((y[,1]-m)^2)/n
			else if(bindata)sum((y[,1]-m*nn)^2)/n
			else sum((y-m)^2)/n
		-sum(fcn(m,p))-sum(mix(r,var0+ss))}
	else if(mean0==1) function(p){
		r <- exp(p[np1+unest])
		r <- c(r,nnest-sum(r))
		ss <- if(!fixed)sum((r-1)^2)/nnest else pmix
		m <- mu1(p,r[nest])
		var0 <- if(censor)sum((y[,1]-m)^2)/n
			else if(bindata)sum((y[,1]-m*nn)^2)/n
			else sum((y-m)^2)/n
		-sum(fcn(m,p))-sum(mix(r,var0+ss))}
	else function(p){
		r <- 1/(1+exp(-p[np1+unest]))
		r <- c(r,0.5*nnest-sum(r))
		ss <- if(!fixed)sum((r-0.5)^2)/nnest else pmix
		m <- mu1(p,r[nest])
		var0 <- if(censor)sum((y[,1]-m)^2)/n
			else if(bindata)sum((y[,1]-m*nn)^2)/n
			else sum((y-m)^2)/n
		-sum(fcn(m,p))-sum(mix(r,var0+ss))}
tlike <- function(p){
	if(mean0==0)r <- c(p[np1+unest],-sum(p[np1+unest]))
	else if(mean0==1){
		r <- exp(p[np1+unest])
		r <- c(r,nnest-sum(r))}
	else if(mean0==2){
		r <- 1/(1+exp(-p[np1+unest]))
		r <- c(r,0.5*nnest-sum(r))}
	m <- mu1(p,r[nest])
	-sum(fcn(m,p)+cc*log(delta))}
#
# check that the likelihood returns an appropriate value and optimize
#
tmp <- like(p)
if(is.na(tmp)||abs(tmp)==Inf)
	stop("Likelihood returns Inf or NA: invalid initial values, wrong model, or probabilities too small to calculate")
if(fscale==1)fscale <- tmp
z0 <- nlm(like,p=p,hessian=TRUE,print.level=print.level,typsize=typsize,
	ndigit=ndigit,gradtol=gradtol,stepmax=stepmax,steptol=steptol,
	iterlim=iterlim,fscale=fscale)
z0$minimum <- z0$minimum-sum(log(delta))
#
# calculate se's
#
if(np==0)cov <- NULL
else if(np==1)cov <- 1/z0$hessian
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)}
se <- sqrt(diag(cov))
maxlike <- sum(tlike(z0$estimate))
#
# calculate random effects, variances, fitted values, and
# mixing shape parameter
#
if(mean0==0){
	r <- c(z0$estimate[np1+unest],-sum(z0$estimate[np1+unest]))
	var <- sum(r^2)/nnest}
else if(mean0==1){
	r <- exp(z0$estimate[np1+unest])
	r <- c(r,nnest-sum(r))
	if(mean0==1)var <- sum((r-1)^2)/nnest}
else if(mean0==2){
	r <- 1/(1+exp(-z0$estimate[np1+unest]))
	r <- c(r,0.5*nnest-sum(r))
	var <- sum((r-0.5)^2)/nnest}
pred <- mu1(z0$estimate,if(mean0==0)0 else if(mean0==1)1 else 0.5)
if(bindata)pred <- nn*pred
rpred <- mu1(z0$estimate,r[nest])
if(bindata)rpred <- nn*rpred
var0 <- if(censor||bindata)sum((y[,1]-rpred)^2)/n else sum((y-rpred)^2)/n
if(!fixed)pmix <- switch(mixture,
	normal=var,
	logistic=sqrt(3*var)/pi,
	gamma=1/var,
	"inverse gamma"=1/var,
	"inverse Gauss"=var,
	Weibull={
		fn <- function(z)gamma(1+2/z)-gamma(1+1/z)^2-var
		uniroot(fn,c(0.02,20))$root},
	beta=0.125/var-0.5)
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
if(!is.null(sh2))sh1 <- sh2
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
z1 <- list(
	call=call,
	response=envir,
	delta=delta,
	distribution=distribution,
	mixture=mixture,
	mixvar=c(var0,var),
	mu=mu1,
	shape=sh1,
	linear=list(lin1,lin2),
	linmodel=list(lin1model,lin2model),
	common=common,
	maxlike=maxlike,
	penalty=z0$minimum-maxlike,
	pred=pred,
	rpred=rpred,
	aic=maxlike+np-1,
	df=n-np+1,
	coefficients=z0$estimate,
	npl=npl,
	npr=nnest-1,
	nps=nps,
	pmix=pmix,
	fixed=fixed,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- c("hnlm","recursive")
return(z1)}

### standard methods
###
##' @export
deviance.hnlm <- function(object, ...) 2*object$maxlike
##' @export
fitted.hnlm <- function(object, recursive=TRUE, ...) if(recursive) object$rpred else object$pred
##' @export
residuals.hnlm <- function(object, recursive=TRUE, ...)
	if(recursive) object$response$y-object$rpred else object$response$y-object$pred

### print method
### 
##' @export
print.hnlm <- function(x, correlation=TRUE, ...) {
  z<-x
sht <- z$nps>0||!is.null(z$shape)
npl <- z$npl
np1 <- z$npl+1
np1a <- z$npl+1
np2 <- z$npl+z$nps
np3 <- np2
np4 <- np3+1
np <- z$npl+z$nps+z$npr
mean0 <- if(z$mixture=="normal"||z$mixture=="logistic"||z$mixture=="Cauchy"||
	 z$mixture=="Laplace")0 else if(z$mixture=="beta")2 else 1
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
if(!is.null(z$dist))cat(z$dist,"distribution\n\n")
cat(z$mixture,"mixing distribution\n\n")
if(z$npl>0||!is.null(z$mu)){
	cat("Location function:\n")
	if(!is.null(attr(z$mu,"formula")))
		cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(!is.null(z$linear[[1]])){
		cat("Linear part:\n")
		print(z$linear[[1]])}}
if(sht){
	cat("\nLog shape function:\n")
	if(!is.null(attr(z$shape,"formula")))
		cat(deparse(attr(z$shape,"formula")),sep="\n")
	else if(!is.null(attr(z$shape,"model"))){
		t <- deparse(attr(z$shape,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(!is.null(z$linear[[2]])){
		cat("Linear part:\n")
		print(z$linear[[2]])}
	if(!is.null(z$family)){
		cat("\n(Log) family function:\n")
		if(!is.null(attr(z$family,"formula")))
			cat(deparse(attr(z$family,"formula")),sep="\n")
		else if(!is.null(attr(z$family,"model"))){
			t <- deparse(attr(z$family,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(!is.null(z$linear[[3]])){
			cat("Linear part:\n")
			print(z$linear[[3]])}}}
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Penalty           ",z$penalty,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
if(npl>0){
	if(z$common)cat("Common parameters:\n")
	else cat("Location parameters:\n")
	cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else if(length(grep("linear",attr(z$mu,"parameters")))>0)
		attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
		else attr(z$mu,"parameters")
	if(!is.null(z$linmodel[[1]]))cname <- c(cname,z$linmodel[[1]])
	coef.table <- cbind(z$coefficients[1:npl],z$se[1:npl])
	if(!z$common){
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table,digits=4,print.gap=2)
		cname <- coef.table <- NULL}}
if(z$common||z$nps>0){
	if(!is.null(z$shape))cname <- c(cname,
			if(is.character(attr(z$shape,"model")))
			attr(z$shape,"model")
		else if(length(grep("linear",attr(z$shape,"parameters")))>0||
			length(grep("mu",attr(z$shape,"parameters")))>0)
			attr(z$shape,"parameters")[grep("\\[",attr(z$shape,"parameters"))]
		else attr(z$shape,"parameters"))
	if(!is.null(z$linmodel[[2]]))cname <- c(cname,z$linmodel[[2]])
	if(!z$common)coef.table <- cbind(z$coefficients[np1a:np2],
		z$se[np1a:np2])
	if(z$common&&is.null(z$family)){
		dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
		print.default(coef.table,digits=4,print.gap=2)}
	if(is.null(z$shape)&&z$nps==1){
		coef.table <- cbind(z$coefficients[np2],z$se[np2])
		cname <- " "}}
if(z$nps>0&&!z$common){
	cat("\nShape parameters:\n")
	dimnames(coef.table) <- list(cname,c("estimate","se"))
	print.default(coef.table,digits=4,print.gap=2)
	cname <- coef.table <- NULL}
if(!is.null(z$family)){
	if(!z$common)cat("\nFamily parameters:\n")
	cname <- c(cname,if(is.character(attr(z$family,"model")))
		attr(z$family,"model")
		else if(length(grep("linear",attr(z$family,"parameters")))>0)
		attr(z$family,"parameters")[grep("\\[",attr(z$family,"parameters"))]
		else attr(z$family,"parameters"))
	if(!is.null(z$linmodel[[3]]))cname <- c(cname,z$linmodel[[3]])
	if(z$common){
		dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
		print.default(coef.table,digits=4,print.gap=2)}
	else {
		coef.table <- cbind(z$coefficients[np3:np],z$se[np3:np])
		dimnames(coef.table) <- list(cname,c("estimate","se"))
		print.default(coef.table,digits=4,print.gap=2)}}
if(z$fixed)
	cat("\nFixed mixing shape parameter:",z$pmix,"\n")
else
	cat("\nMixing shape parameter:",z$pmix,"\n")
cat("\nVariances: conditional = ",z$mixvar[1],", mixing = ",z$mixvar[2],"\n",
	sep="")
cat("\nRandom effect parameters:\n")
r <- c(z$coefficients[np4:np],-sum(z$coefficients[np4:np]))
if(mean0==1){
	rr <- exp(z$coefficients[np4:np])
	rr <- c(rr,z$npr+1-sum(rr))
	r[z$npr+1] <- log(rr[z$npr+1])}
else if(mean0==2){
	rr <- 1/(1+exp(-z$coefficients[np4:np]))
	rr <- c(rr,0.5*z$npr+0.5-sum(rr))
	r[z$npr+1] <- log(rr[z$npr+1]/(1-rr[z$npr+1]))}
coef.table <- cbind(r,c(z$se[np4:np],NA))
if(mean0)coef.table <- cbind(coef.table,rr)
dimnames(coef.table) <- list(1:(z$npr+1),
	c(if(mean0)"estimate"else"effect","se",if(mean0)"effect"))
print.default(coef.table,digits=4,print.gap=2)
if(correlation){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr,digits=4)}
invisible(z)}
