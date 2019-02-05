#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 2000, 2001, 2002 J.K. Lindsey
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
#     gnlmix(y=NULL, distribution="normal", mixture="normal",
#	random=NULL, nest=NULL, mu=NULL, shape=NULL, linear=NULL,
#	pmu=NULL, pshape=NULL, pmix=NULL, delta=1, common=FALSE,
#	envir=parent.frame(), print.level=0, typsize=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
#	iterlim=100, fscale=1, eps=1.0e-4, points=5, steps=10)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with one arbitrary
#  random parameter



##' Generalized Nonlinear Regression with a Random Parameter
##' 
##' \code{gnlmix} fits user-specified nonlinear regression equations to one or
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
##' binomial data, or an object of class, \code{response} (created by
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
##' @param mixture The mixing distribution for the random parameter: normal,
##' Cauchy, logistic, Laplace, inverse Gauss, gamma, inverse gamma, Weibull,
##' beta, simplex, or two-sided power. The first four have zero location
##' parameter, the next three have unit location parameter, and the last two
##' have location parameter set to 0.5.
##' @param random The name of the random parameter in the \code{mu} formula.
##' @param nest The variable classifying observations by the unit upon which
##' they were observed. Ignored if \code{y} or \code{envir} has class,
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
##' @param pmix Initial estimate for the logarithm of the dispersion parameter
##' of the mixing distribution.
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
##' @param eps Precision of the Romberg integration.
##' @param steps For the Romberg integration, the maximum number of steps, by
##' default set to 10.
##' @param points For the Romberg integration, the number of extrapolation
##' points so that 2*points is the order of integration, by default set to 5;
##' points=2 is Simpson's rule.
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
### @seealso \code{\link[growth]{carma}}, \code{\link[rmutil]{finterp}},
### \code{\link[growth]{elliptic}}, \code{\link[repeated]{glmm}},
### \code{\link[repeated]{gnlmm}}, \code{\link[gnlm]{gnlr}},
### \code{\link[repeated]{hnlmix}}, \code{\link[repeated]{kalseries}},
### \code{\link[gnlm]{nlr}}, \code{\link[stats]{nls}}.
##' @keywords models
##' @examples
##' 
##' dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
##' #y <- rgamma(20,shape=2+0.3*dose,scale=2)+rep(rnorm(4,0,4),rep(5,4))
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
##' gnlmix(reps, mu=~a+b*dose+rand, random="rand", pmu=c(8.7,0.25),
##' 	pshape=3.44, pmix=2.3)
##' \dontrun{
##' # gamma model with log link and random normal intercept fitted three ways
##' glmm(y~dose, family=Gamma(link=log), nest=individuals, data=reps, points=8)
##' gnlmm(reps, distribution="gamma", mu=~exp(a+b*dose), pmu=c(2,0.03),
##' 	psh=1, psd=0.3)
##' gnlmix(reps, distribution="gamma", mu=~exp(a+b*dose+rand), random="rand",
##' 	pmu=c(2,0.04), pshape=1, pmix=-2)
##' 
##' # gamma model with log link and random gamma mixtures
##' gnlmix(reps, distribution="gamma", mixture="gamma",
##' 	mu=~exp(a*rand+b*dose), random="rand", pmu=c(2,0.04),
##' 	pshape=1.24, pmix=3.5)
##' gnlmix(reps, distribution="gamma", mixture="gamma",
##' 	mu=~exp(a+b*dose)*rand, random="rand", pmu=c(2,0.04),
##' 	pshape=1.24, pmix=2.5)
##' }
##' @export gnlmix
##' @importFrom graphics par plot 
##' @importFrom stats as.formula dbeta dbinom dcauchy dexp dgamma dlogis dnbinom dnorm dpois dt dweibull gaussian glm glm.control model.frame model.matrix na.fail nlm pbeta pcauchy pexp pgamma pgeom plogis pnbinom pnorm ppois pt pweibull qnorm summary.glm terms uniroot update.formula
##' @import rmutil
##' @useDynLib repeated, .registration=TRUE
gnlmix <- function(y=NULL, distribution="normal", mixture="normal",
	random=NULL, nest=NULL, mu=NULL, shape=NULL, linear=NULL,
	pmu=NULL, pshape=NULL, pmix=NULL, delta=1, common=FALSE,
	envir=parent.frame(), print.level=0, typsize=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
	iterlim=100, fscale=1, eps=1.0e-4, points=5, steps=10){

int1 <- function(ff, aa, bb)
	.C("romberg_c",
		ff,
		as.double(aa),
		as.double(bb),
		len=as.integer(nnest),
		eps=as.double(eps),
		pts=as.integer(points),
		max=as.integer(steps),
		err=integer(1),
		res=double(nnest),
		PACKAGE="repeated")$res
inta <- function(f){
	ff <- function(x) f(1/x)/(x*x)
	int1(ff,neg1,zero)+int1(f,neg1,pos1)+int1(ff,zero,pos1)}
intb <- function(f){
	ff <- function(x) f(1/x)/(x*x)
	int1(f,zero,pos1)+int1(ff,zero,pos1)}
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
#
# check mixture
#
mixture <- match.arg(mixture,c("normal","logistic","Cauchy","Laplace",
	"gamma","inverse gamma","inverse Gauss","Weibull","Levy","beta",
	"simplex","two-sided power"))
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
np <- npl+nps+1
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("gnlmix only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("gnlmix does not handle data with NAs")}
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
if(is.null(pmix))stop("a value must be supplied for pmix")
p <- c(pmu,pshape,pmix)
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
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("gnlmix does not handle data with NAs")
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
		delta <- as.vector(envir$response$delta[,col])}
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("gnlmix only handles univariate responses")
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	nest <- covind(y)
	type <- y$type
	y <- response(y)}
if(is.null(nest)||length(unique(nest))==1)
	stop("appropriate nest indicator required")
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
	lc <- ifelse(y[,2]==-1,0,1)}
else cc <- 1
#
# fixed values for integration routines
#
nnest <- length(unique(nest))
neg1 <- rep(-1,nnest)
pos1 <- rep(1,nnest)
zero <- rep(0,nnest)
#
# prepare weights and unit of measurement
#
wt <- rep(1,n)
if(length(delta)==1)delta <- rep(delta,n)
else if(length(delta)!=n)
	stop("delta must be the same length as the other variables")
#
# check that location function returns appropriate values
#
if(any(is.na(mu1(pmu,0))))stop("The location model returns NAs: probably invalid initial values")
if(distribution=="Levy"&&any(y<=mu1(p)))
	stop("location parameter must be strictly less than corresponding observation")
#
# check that shape function returns appropriate values
#
if(distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"&&any(is.na(sh1(p))))
	stop("The shape model returns NAs: probably invalid initial values")
if(distribution=="Pareto"&&exp(sh1(p))<=1)stop("shape parameters must be > 0")
#
# set up distribution
#
if(!censor)fcn <- switch(distribution,
	binomial=function(p,r) dbinom(y[,1],y[,1]+y[,2],mu1(p,r)),
	"beta binomial"=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))
		t <- s*m
		u <- s*(1-m)
		exp(lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u)+lchoose(nn,y[,1]))},
	"double binomial"=function(p,r)
		exp(.C("ddb_c",as.integer(y[,1]),as.integer(nn),
			as.double(mu1(p,r)),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res),
	"mult binomial"=function(p,r)
		exp(.C("dmb_c",as.integer(y[,1]),as.integer(nn),
			as.double(mu1(p,r)),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res),
	Poisson=function(p,r)dpois(y,mu1(p,r)),
	"negative binomial"=function(p,r)dnbinom(y,exp(sh1(p)),mu1(p,r)),
	"double Poisson"=function(p,r)
		exp(.C("ddp_c",as.integer(y),as.integer(my),
			as.double(mu1(p,r)),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res),
	"mult Poisson"=function(p,r)
		exp(.C("dmp_c",as.integer(y),as.integer(my),
			as.double(mu1(p,r)),as.double(exp(sh1(p))),
			as.integer(n),as.double(wt),res=double(n),
			#DUP=FALSE,
			PACKAGE="repeated")$res),
	"gamma count"=function(p,r) {
		m <- mu1(p,r)
		s <- exp(sh1(p))
		u <- m*s
		ifelse(y==0,pgamma(u,(y+1)*s,1,lower.tail=FALSE),
			pgamma(u,y*s+(y==0),1)-pgamma(u,(y+1)*s,1))},
	Consul=function(p,r) {
		m <- mu1(p,r)
		t <- sh1(p)
		s <- exp(t)
		exp(log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-
			y*t-lgamma(y+1))},
	logarithmic=function(p,r) {
		m <- 1/(1+exp(-mu1(p,r)))
		exp(y*log(m)-log(y)-log(-log(1-m)))},
	geometric=function(p,r) {
		m <- mu1(p,r)
		exp(y*log(m)-(y+1)*log(1+m))},
	normal=function(p,r) dnorm(y,mu1(p,r),exp(sh1(p)/2)),
	"inverse Gauss"=function(p,r) {
		m <- mu1(p,r)
		t <- sh1(p)
		exp(-(t+(y-m)^2/(y*exp(t)*m^2)+log(2*pi*y^3))/2)},
	logistic=function(p,r)dlogis(y,mu1(p,r),exp(sh1(p))*sqrt(3)/pi),
	Cauchy=function(p,r) dcauchy(y,mu1(p,r),exp(sh1(p)/2)),
	Laplace=function(p,r) {
		t <- sh1(p)
		exp(-(abs(y-mu1(p,r))/exp(t)+t))/2},
	Levy=function(p,r) {
		m <- mu1(p,r)
		s <- exp(sh1(p))
		exp(0.5*log(s/(2*pi))-1.5*log(y-m)-s/(2*(y-m)))/2},
	Pareto=function(p,r) {
		s <- exp(sh1(p))
		t <- 1/(mu1(p,r)*(s-1))
		exp(log(s*t)-(s+1)*log(1+y*t))},
	exponential=function(p,r) dexp(y,1/mu1(p,r)),
	gamma=function(p,r){
		s <- exp(sh1(p))
		dgamma(y,s,scale=mu1(p,r)/s)},
	Weibull=function(p,r) dweibull(y,exp(sh1(p)),mu1(p,r)),
	"extreme value"=function(p,r)
		exp(y)*dweibull(exp(y),exp(sh1(p)),mu1(p,r)),
	beta=function(p,r) {
		s <- exp(sh1(p))
		m <- mu1(p,r)*s
		s <- s-m
		dbeta(y,m,s)},
	simplex=function(p,r) {
		m <- mu1(p,r)
		t <- sh1(p)
		s <- exp(t)
		exp(-(((y-m)/(m*(1-m)))^2/(y*(1-y)*s)+
			t+3*log(y*(1-y))/2+log(2*pi)/2))},
	"two-sided power"=function(p,r){
		m <- mu1(p,r)
		t <- sh1(p)
		s <- exp(t)
		exp(t+(s-1)*ifelse(y<m,log(y/m),log((1-y)/(1-m))))})
else fcn <- switch(distribution,
	Poisson=function(p,r){
		m <- mu1(p,r)
		dpois(y[,1],m)^cc*(lc-rc*ppois(y[,1],m))},
	"negative binomial"=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))
		dnbinom(y[,1],s,m)^cc*(lc-rc*pnbinom(y[,1],s,m))},
	geometric=function(p,r){
		m <- mu1(p,r)
		(y[,1]*log(m)-(y[,1]+1)*log(1+m))^cc*
			(lc-rc*pgeom(y[,1],1/(1+m)))},
	normal=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p)/2)
		dnorm(y[,1],m,s)^cc*(lc-rc*pnorm(y[,1],m,s))},
	"inverse Gauss"=function(p,r){
		m <- mu1(p,r)
		t <- sh1(p)
		s <- exp(t)
		v <- sqrt(s*y[,1]/2)
		exp(-cc*(t+(y[,1]-m)^2/(y[,1]*s*m^2)+log(2*pi*y[,1]^3))/2)*
			(lc-rc*(pnorm((y[,1]/m-1)/v)+
			exp(2/(m*s))*pnorm(-(y[,1]/m+1)/v)))},
	logistic=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))*sqrt(3)/pi
		dlogis(y[,1],m,s)^cc*(lc-rc*plogis(y[,1],m,s))},
	Cauchy=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p)/2)
		dcauchy(y[,1],m,s)^cc*(lc-rc*pcauchy(y[,1],m,s))},
	Laplace=function(p,r){
		u <- abs(y[,1]-mu1(p,r))/s
		t <- exp(-u)/2
		v <- sh1(p)
		s <- exp(v)
		-exp(cc*(u+v+log(2)))*(lc-rc*(ifelse(u<0,t,1-t)))},
	Pareto=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))
		t <- 1/(m*(s-1))
		exp(cc*(log(s*t)-(s+1)*log(1+y[,1]*t)))*
			(lc-rc*((1+y[,1]/(m*(s-1)))^(-s)))},
	exponential=function(p,r){
		m <- mu1(p,r)
		dexp(y[,1],1/m)^cc*(lc-rc*pexp(y[,1],1/m))},
	gamma=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))
		dgamma(y[,1],s,scale=m/s)^cc*
			(lc-rc*pgamma(y[,1],s,scale=m/s))},
	Weibull=function(p,r){
		m <- mu1(p,r)
		s <- exp(sh1(p))
		dweibull(y[,1],s,m)^cc*(lc-rc*pweibull(y[,1],s,m))},
	"extreme value"=function(p,r){
		yy <- exp(y[,1])
		m <- mu1(p,r)
		s <- exp(sh1(p))
		(yy*dweibull(yy,s,m))^cc*(lc-rc*pweibull(yy,s,m))})
#
# set up mixing distribution
#
mix <- switch(mixture,
        normal=function(p,r) dnorm(r,0,exp(p/2)),
	logistic=function(p,r) dlogis(r,0,exp(p)*sqrt(3)/pi),
	Cauchy=function(p,r) dcauchy(r,0,exp(p)),
        Laplace=function(p,r) {
		tmp <- exp(p)
		exp(-abs(r)/tmp)/(2*tmp)},
	gamma=function(p,r) {
		tmp <- exp(p)
		dgamma(r,tmp,scale=1/tmp)},
	"inverse gamma"=function(p,r) {
		tmp <- exp(p)
		dgamma(1/r,tmp,scale=1/tmp)/r^2},
        "inverse Gauss"=function(p,r)
        	exp(-((r-1)^2/(r*exp(p))+p+log(2*pi*r^3))/2),
	Weibull=function(p,r) dweibull(r,exp(p),1),
	beta=function(p,r) {
		tmp <- 0.5*exp(p)
		dbeta(r,tmp,tmp)},
        simplex=function(p,r) exp(-(((r-0.5)/0.25)^2/(r*(1-r)*exp(p))+
			p+3*log(r*(1-r))/2+log(2*pi)/2)),
	"two-sided power"=function(p,r)
		p*exp((p-1)*ifelse(r<0.5,log(r/0.5),log((1-r)/0.5))))
#
# combine to create the appropriate likelihood function
#
if(mixture=="normal"||mixture=="logistic"||mixture=="Cauchy"||mixture=="Laplace")
	like <- function(p){
		fn <- function(r)
			mix(p[np],r)*capply(fcn(p,r[nest])*delta^cc,nest,prod)
		-sum(log(inta(fn)))}
else if(mixture=="gamma"||mixture=="inverse gamma"||mixture=="inverse Gauss"||
	mixture=="Weibull")
	like <- function(p){
		fn <- function(r)
			mix(p[np],r)*capply(fcn(p,r[nest])*delta^cc,nest,prod)
		-sum(log(intb(fn)))}
else like <- function(p){
		fn <- function(r)
			mix(p[np],r)*capply(fcn(p,r[nest])*delta^cc,nest,prod)
		-sum(log(int1(fn,0,1)))}
#
# check that the likelihood returns an appropriate value and optimize
#
tmp <- like(p)
if(is.na(tmp)||abs(tmp)==Inf)
	stop("Likelihood returns Inf or NAs: invalid initial values, wrong model, or probabilities too small to calculate")
if(fscale==1)fscale <- tmp
z0 <- nlm(like,p=p,hessian=TRUE,print.level=print.level,typsize=typsize,
	ndigit=ndigit,gradtol=gradtol,stepmax=stepmax,steptol=steptol,
	iterlim=iterlim,fscale=fscale)
#
# calculate fitted values and raw residuals
#
fitted.values <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial")
		as.vector((y[,1]+y[,2])*mu1(z0$estimate,0))
	else as.vector(mu1(z0$estimate,0))
residuals <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial")
		y[,1]-fitted.values
	else y-fitted.values
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
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
if(!is.null(sh2))sh1 <- sh2
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	mixture=mixture,
	likefn=like,
	mu=mu1,
	shape=sh1,
	mix=NULL,
	censor=censor,
	linear=list(lin1,lin2),
	linmodel=list(lin1model,lin2model),
	common=common,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=z0$minimum+np,
	df=n-np,
	coefficients=z0$estimate,
	npl=npl,
	npm=1,
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
