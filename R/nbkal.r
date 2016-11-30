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
#     nbkal(response, times, mu, preg, pdepend, kalman=TRUE,
#	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
#	fscale=1, iterlim=100, typsize=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit a correlated negative binomial model with Kalman update
#  Adapted from Gauss code written by P. Lambert



##' Negative Binomial Models with Kalman Update
##' 
##' \code{nbkal} fits a negative binomial regression with Kalman update over
##' time. The variance is proportional to the mean function, whereas, for
##' \code{\link[repeated]{kalcount}} with exponential intensity, it is a
##' quadratic function of the mean.
##' 
##' Marginal and individual profiles can be plotted using
##' \code{\link[rmutil]{mprofile}} and \code{\link[rmutil]{iprofile}} and
##' residuals with \code{\link[rmutil]{plot.residuals}}.
##' 
##' 
##' @param response A list of two column matrices with counts and corresponding
##' times for each individual, one matrix or dataframe of counts, or an object
##' of class, response (created by \code{\link[rmutil]{restovec}}) or repeated
##' (created by \code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}).
##' @param times When response is a matrix, a vector of possibly unequally
##' spaced times when they are the same for all individuals or a matrix of
##' times. Not necessary if equally spaced. Ignored if response has class,
##' response or repeated.
##' @param mu The mean function.
##' @param preg The initial parameter estimates for the mean function.
##' @param pdepend The estimates for the dependence parameters, either one or
##' three.
##' @param kalman If TRUE, fits the kalman update model, otherwise, a standard
##' negative binomial distribution.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of classes \code{nbkal} and \code{recursive} is returned.
##' @author P. Lambert and J.K. Lindsey
### @seealso \code{\link[repeated]{gar}}, \code{\link[repeated]{gnlmm}},
### \code{\link[gnlm]{gnlr}}, \code{\link[rmutil]{iprofile}}
### \code{\link[repeated]{kalcount}}, \code{\link[rmutil]{mprofile}},
### \code{\link[rmutil]{read.list}}, \code{\link[rmutil]{rmna}},
### \code{\link[rmutil]{restovec}}, \code{\link[rmutil]{tcctomat}},
### \code{\link[rmutil]{tvctomat}}.
##' @references Lambert, P. (1996) Applied Statistics 45, 31-38.
##' 
##' Lambert, P. (1996) Biometrics 52, 50-55.
##' @keywords models
##' @examples
##' 
##' y <- matrix(rnbinom(20,5,0.5), ncol=5)
##' times <- matrix(rep(seq(10,50,by=10),4), ncol=5, byrow=TRUE)
##' y0 <- matrix(rep(rnbinom(5,5,0.5),4), ncol=5, byrow=TRUE)
##' mu <- function(p) p[1]*log(y0)+(times<30)*p[2]*
##' 	(times-30)+(times>30)*p[3]*(times-30)
##' nbkal(y, preg=c(1.3,0.008,-0.05), times=times, pdep=1.2, mu=mu)
##' 
##' @aliases nbkal print.nbkal
##' @export nbkal
nbkal <- function(response, times, mu, preg, pdepend, kalman=TRUE,
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	fscale=1, iterlim=100, typsize=abs(p), stepmax=10*sqrt(p%*%p)){
#
# negative binomial likelihood function
#
likenb <- function(p){
	nu <- abs(p[length(p)])
	eta <- nu*exp(mu(p))
	sum(lgamma(response$response$y+1)+lgamma(eta)
		-lgamma(response$response$y+eta)+
		-eta*log(nu)+(response$response$y+eta)*log(1+nu))}
#
# dynamic negative binomial likelihood function
#
likekal <- function(p){
	eta <- exp(mu(p))
	phi <- exp(-abs(p[length(p)-np1+1]))
	if(np1>1){
		delta <- abs(p[length(p)-np1+2])
		alpha <- exp(p[length(p)-np1+3])
		alpha <- alpha/(1+alpha)}
	else delta <- alpha <- 1
	nm <- like <- 0
	for(ii in 1:nind){
		kapred <- 1
		upspred <- 1/0.000001
		for(jj in 1:nobs(response)[ii]){
			nm <- nm+1
			if(jj>1){
				kapred <- kap
				upspred <- exp(-phi*(response$response$times[nm]-
					response$response$times[nm-1]))*ups}
			ups <- upspred+alpha*eta[nm]/delta
			kap <- kapred+alpha*(response$response$y[nm]-kapred*
				eta[nm])/(upspred*delta+eta[nm])
			tp3 <- upspred*delta
			tp4 <- kapred*tp3
			like <- like+(tp4-delta+1)*log(tp3)-
				(tp4-delta+1+response$response$y[nm])*
				log(tp3+eta[nm])+
				response$response$y[nm]*log(eta[nm])-
				log(response$response$y[nm]+tp4-delta+1)-
				lbeta(tp4-delta+1,response$response$y[nm]+1)}}
	-like}
#
# negative binomial likelihood function for recursive prediction
#
likepred <- function(p){
	eta <- exp(mu(p))
	phi <- exp(-abs(p[length(p)-np1+1]))
	if(np1>1){
		delta <- abs(p[length(p)-np1+2])
		alpha <- exp(p[length(p)-np1+3])
		alpha <- alpha/(1+alpha)}
	else delta <- alpha <- 1
	nm <- 0
	rpred <- NULL
	for(ii in 1:nind){
		kapred <- 1
		upspred <- 1/0.000001
		for(jj in 1:nobs(response)[ii]){
			nm <- nm+1
			if(jj>1){
				kapred <- kap
				upspred <- exp(-phi*(response$response$times[nm]-
					response$response$times[nm-1]))*ups}
			ups <- upspred+alpha*eta[nm]/delta
			kap <- kapred+alpha*(response$response$y[nm]-kapred*
				eta[nm])/(upspred*delta+eta[nm])
			rpred <- c(rpred,eta[nm]*(kap*ups*delta-delta+1)/
				(ups*delta))}}
	rpred}
call <- sys.call()
#
# set up response if not a data object
#
if(!inherits(response,"repeated")){
	if(!inherits(response,"response"))response <- restovec(response,times)
	response <- rmna(response=response)}
yy <- ifelse(response$response$y==0,1,response$response$y)
#
# optimize using nlm
#
nind <- length(nobs(response))
p <- c(preg,pdepend)
np <- length(p)
np1 <- length(pdepend)
if(!kalman)z0 <- nlm(likenb, p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
else z0 <- nlm(likekal, p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
rpred <- if(kalman)likepred(z0$estimate) else NULL
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
z <- list(
	call=call,
	mu=mu,
	response=response$response,
	maxlike=z0$minimum,
	aic=z0$minimum+np,
	df=dim(response$response$y)[1]-np,
	np=np,
	kalman=kalman,
	full=length(pdepend)==3,
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=corr,
	pred=exp(mu(z0$estimate)),
	rpred=rpred,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
if(kalman)class(z) <- c("nbkal","recursive")
else class(z) <- "nbkal"
return(z)}

### standard methods
###
##' @export
deviance.nbkal <- function(object, ...) 2*object$maxlike
##' @export
fitted.nbkal <- function(object, recursive=TRUE, ...) if(recursive) object$rpred else object$pred
##' @export
residuals.nbkal <- function(object, recursive=TRUE, ...)
	if(recursive) object$response$y-object$rpred else object$response$y-object$pred

### print method
###
##' @export
print.nbkal <- function(x, digits=max(3,.Options$digits-3), ...){
  z<-x
np1 <- ifelse(z$full&z$kalman,3,1)
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
t <- deparse(z$mu)
cat("Location function:",t[2:length(t)],sep="\n")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Location parameters\n")
coef.table <- cbind(z$coef[1:(z$np-np1)],z$se[1:(z$np-np1)])
cname <- NULL
for(i in 1:dim(coef.table)[1])cname <- c(cname,paste("p",i,sep=""))
dimnames(coef.table) <- list(cname, c("estimate","se"))
print.default(coef.table, digits=digits, print.gap=2)
cat("\nNonlinear parameters\n")
cname <- ifelse(z$kalman,"phi","nu")
tmp <- exp(-exp(-z$coef[(z$np-np1+1)]))
if(np1>1){
	cname <- c(cname,"delta","alpha")
	tmp <- c(tmp,z$coef[(z$np-np1+2)],exp(z$coef[(z$np-np1+3)])
		/(1+exp(z$coef[(z$np-np1+2)])))}
coef.table <- cbind(z$coef[(z$np-np1+1):z$np],
	z$se[(z$np-np1+1):z$np],tmp)
dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
print.default(coef.table, digits=digits, print.gap=2)
cat("\nCorrelation matrix\n")
print.default(z$corr, digits=digits)}
