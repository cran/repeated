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
#     glmm(formula, family=gaussian, data=list(), weights=NULL,
#	offset=NULL, nest, delta=1, maxiter=20, points=10, print.level=0,
#	control=glm.control(epsilon=0.0001,maxit=10,trace=FALSE))
#
#  DESCRIPTION
#
#    A function to fit generalized linear mixed models with normal
#  random effect

# does not work: na.omit


##' Generalized Linear Mixed Models
##' 
##' \code{glmm} fits a generalized linear mixed model with a random intercept
##' using a normal mixing distribution computed by Gauss-Hermite integration.
##' For the normal, gamma, and inverse Gaussian distributions, the deviances
##' supplied are -2 log likelihood, not the usual \code{\link{glm}} deviance;
##' the degrees of freedom take into account estimation of the dispersion
##' parameter.
##' 
##' If weights and/or offset are to be used or the formula transforms some
##' variables, all of the data must be supplied in a dataframe. Because the
##' \code{\link{glm}} function is such a hack, if this is not done, weird error
##' messages will result.
##' 
##' na.omit is not allowed.
##' 
##' 
##' @param formula A symbolic description of the model to be fitted. If it
##' contains transformations of the data, including cbind for binomial data, a
##' dataframe must be supplied.
##' @param family A description of the error distribution and link function to
##' be used in the model; see \code{\link{family}} for details.
##' @param data A dataframe containing the variables in the model, that is
##' optional in simple cases, but required in certain situations as specified
##' elsewhere in this help page.
##' @param weights An optional weight vector. If this is used, data must be
##' supplied in a data.frame.
##' @param offset The known component in the linear predictor. If this is used,
##' data must be supplied in a data.frame. An offset cannot be specified in the
##' model formula.
##' @param nest The variable classifying observations by the unit (cluster)
##' upon which they were observed.
##' @param delta If the response variable has been transformed, this is the
##' Jacobian of that transformation, so that AICs are comparable.
##' @param maxiter The maximum number of iterations of the outer loop for
##' numerical integration.
##' @param points The number of points for Gauss-Hermite integration of the
##' random effect.
##' @param print.level If set equal to 2, the log probabilities are printed out
##' when the underflow error is given.
##' @param control A list of parameters for controlling the fitting process.
##' @return \code{glmm} returns a list of class \code{glmm}
##' @author J.K. Lindsey
### @seealso \code{\link{family}}, \code{\link[gnlm]{fmr}}, \code{\link{glm}},
### \code{\link{glm.control}}, \code{\link[repeated]{gnlmix}},
### \code{\link[repeated]{gnlmm}}, \code{\link[gnlm]{gnlr}},
### \code{\link[gnlm]{gnlr3}}, \code{\link[repeated]{hnlmix}},
### \code{\link[stats]{nls}}.
##' @keywords models
##' @examples
##' 
##' # Poisson counts
##' nest <- gl(5,4)
##' y <- rpois(20,5+2*as.integer(nest))
##' # overdispersion model
##' glmm(y~1, family=poisson, nest=gl(20,1), points=3)
##' # clustered model
##' glmm(y~1, family=poisson, nest=nest, points=3)
##' #
##' # binomial data with model for overdispersion
##' df <- data.frame(r=rbinom(10,10,0.5), n=rep(10,10), x=c(rep(0,5),
##' 	rep(1,5)), nest=1:10)
##' glmm(cbind(r,n-r)~x, family=binomial, nest=nest, data=df)
##' 
##' @export glmm
glmm <- function(formula, family=gaussian, data=list(), weights=NULL,
	offset=NULL, nest, delta=1, maxiter=20, points=10, print.level=0,
	control=glm.control(epsilon=0.0001,maxit=10,trace=FALSE)){
call <- sys.call()
#
# find the design matrix
#
md <- missing(data)
if(missing(data))data <- parent.frame()
else if(inherits(data,"repeated"))data <- as.data.frame(data)
mf <- model.frame(terms(formula,data=data),data,na.action=na.fail)
#
# find the distribution
#
if(is.character(family))family <- get(family)
if(is.function(family))family <- family()
if(is.vector(mf[,1])||(is.matrix(mf[,1])&&dim(mf[,1])[2]==1)){
	y <- as.vector(mf[,1])
	if(family$family=="binomial")y2 <- rep(1,length(y))-y}
else {
	y <- as.vector(mf[,1][,1])
	y2 <- as.vector(mf[,1][,2])}
#
# check the nesting variable
#
slen <- length(y)
if(!md)nest <- model.frame(terms(as.formula(paste("~",deparse(substitute(nest)),sep="")),data=data),data,na.action=na.fail)[[1]]
if(length(nest)!=slen)
	stop("the nesting variable is not the same length as the other variables")
if(is.factor(nest))nest <- as.numeric(nest)
else nest <- as.numeric(factor(nest))
if(length(unique(nest))!=max(nest))
	stop("(codes of) the nesting variable must be consecutively numbered")
nind <- length(unique(nest))
if(length(nest)==nind&&family$family!="binomial"&&family$family!="poisson")
	stop("Some individuals must have more than one observation")
#
# set up Gauss-Hermite quadrature points
#
i <- rep(1:slen,points)
ii <- rep(1:nind,points)
k <- NULL
for(j in 1:points)k <- c(k,nest+(j-1)*max(nest))
k <- as.integer(k)
quad <- gauss.hermite(points)
sd <- quad[rep(1:points,rep(slen,points)),1]
qw <- quad[rep(1:points,rep(nind,points)),2]
#
# modify formula to contain standard deviation of mixture
#
nmodel <- update.formula(formula,.~.+sd)
nnmf <- if(md)mf[i,,drop=FALSE] else data[i,,drop=FALSE]
#
# set up adjusted weights and offset
#
if(missing(weights))lwt <- rep(1,slen)
else {
	if(!md){
		col <- match(deparse(substitute(weights)),names(data))
		if(is.na(col))stop("weights not found")
		lwt <- data[[col]]}
	else lwt <- as.vector(weights)}
nobs <- sum(lwt)
if(missing(offset))off <- rep(0,slen)
else if(md)
	stop("offset can only be used when a data.frame is specified")
else {
	col <- match(deparse(substitute(offset)),names(data))
	if(is.na(col))stop("offset not found")
	off <- data[[col]]}
#
# put everything in dataframe and run glm to find initial estimates
#
if(!md){
	data <- data.frame(data,lwt,off)
	zz <- glm(formula,family=family,data=data,control=control,
		weights=lwt,offset=off)}
else {
	if(!is.null(weights))
		stop("weights can only be used when a data.frame is specified")
	zz <- glm(formula,family=family,data=data,control=control)}
#
# expand offset and find deviance and d.f.
#
offset <- off[i]
ndev <- zz$deviance
rdf <- zz$df.res-1
ndf <- zz$df.null
fv <- family$linkinv(zz$linear[i]+sd)
sc <- ndev/nobs
#
# define log likelihood function
#
fpw <- switch(family$family,
	binomial= function()
		y*log(fv+0.0001)+(y2)*log(1-fv+0.0001),
	poisson= function() -fv+y*log(fv),
	Gamma= function() (log(y/fv)-y/fv)/sc,
	gaussian= function() -(y-fv)^2/sc/2,
	inverse.gaussian= function() -(y-fv)^2/(y*fv^2)/(2*sc))
#
# iterate, calling glm and checking for underflow
#
for(j in 1:maxiter){
	under  <- 0
	odev <- ndev
	ppr <- NULL
	for(ij in split(lwt*fpw(),k))ppr <- c(ppr,sum(ij))
	if(any(is.na(ppr)))stop("NAs - try another link function")
	if(max(ppr)-min(ppr)>1410){
		if(print.level==2)cat("Log probabilities:\n",ppr,"\n\n")
		stop("Product of probabilities is too small to calculate.\n Try fewer points.")}
	if(any(ppr > 705))under <- 705-max(ppr)
	else if(any(ppr < -705))under <- -705-min(ppr)
	pw <- qw*exp(ppr+under)
	pr <- NULL
	for(ij in split(pw,ii))pr <- c(pr,sum(ij))
	pw <- lwt*(pw/pr[ii])[k]
	nmf <- data.frame(nnmf,sd,pw,offset)
	z <- glm(nmodel,family=family,data=nmf,weights=pw,
		offset=offset,control=control)
	fv <- family$linkinv(z$linear)
	ndev <- -sum(log(pr)-under)
	sc <- z$dev/nobs
	if((odev-ndev)^2<0.00001)break}
#
# prepare object to return
#
z$deviance <- -2*sum(log(pr)-under)
formula <- update.formula(formula,.~1)
class(formula) <- "formula"
if(!md)
	zz <- glm(formula,family=family,data=data,control=control,
		weights=lwt,offset=off)
else zz <- glm(formula,family=family,data=data,control=control)
#
# calculate corrected AIC and deviance
#
switch(family$family,
	binomial={
		sc <- NULL
		z$aic <- z$deviance-2*sum(lwt*lchoose(y+y2,y))
		z$deviance <- z$deviance+2*sum(lwt*(y*
			log(ifelse(y,y,1))+y2*log(ifelse(y2,y2,1))
			-(y+y2)*log(ifelse(y+y2,y+y2,1))))},
	poisson={
		sc <- NULL
		z$aic <- z$deviance+2*sum(lwt*lgamma(y+1))
		z$deviance <- z$deviance+2*sum(lwt*(-y+
			y*log(ifelse(y,y,1))))},
	Gamma={
		sc1 <- zz$dev/sum(lwt)
		z$null.deviance <- 2*sum(lwt*(log(sc1)/sc1-(log(y/zz$fit)-y/zz$fit)/sc1+log(y)+lgamma(1/sc1)-log(delta)))
		z$deviance <- z$deviance+2*sum(lwt*(log(sc)/sc+log(y)+lgamma(1/sc)-log(delta)))
		z$aic <- z$deviance+2},
	gaussian={
		z$null.deviance <- sum(lwt)*(log(2*pi*zz$dev/sum(lwt))+1)-2*sum(lwt*log(delta))
		z$deviance <- z$deviance+sum(lwt)*log(2*pi*sc)-2*sum(lwt*log(delta))
		z$aic <- z$deviance+2},
	inverse.gaussian={
		z$null.deviance <- sum(lwt*(log(2*pi*zz$dev/sum(lwt)*y^3))+1-2*log(delta))
		z$deviance <- z$deviance+sum(lwt*(log(2*pi*sc*y^3)-2*log(delta)))
		z$aic <- z$deviance+2})
z$call <- call
z$aic <- z$aic+2*z$qr$rank
z$df.null <- ndf-!is.null(sc)
z$df.residual <- rdf-!is.null(sc)
z$iter <- j
z$scale <- sc
class(z) <- c("glmm",class(z))
z}

### standard print methods
###
print.glmm <- function(z,...){
#print.glm(z,...)
  class(z) <- "glm"
  print(z)
if(!is.null(z$scale)){
	cat("\nModel deviances are -2 log likelihood\n")
	cat("Model dispersion:      ",z$scale,"\n")}
cat("Normal mixing variance:",z$coef[names(z$coef)=="sd"]^2,"\n")}

print.summary.glmm <- function(z,...){
#print.summary.glm(z,...)
class(z) <- "summary.glm"
print(z)
if(!is.null(z$scale)){
	cat("Model deviances are -2 log likelihood\n")
	cat("Model dispersion:      ",z$scale,"\n")}
cat("Normal mixing variance:",z$coef[rownames(z$coef)=="sd",1]^2,"\n")}

summary.glmm <- function(z,...){
zz <- summary.glm(z,...)
class(zz) <- c("summary.glmm",class(zz))
if(!is.null(z$scale))zz$scale <- z$scale
zz}
