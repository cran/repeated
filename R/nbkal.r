#
#  repeated : A Library of Repeated Measurements Models
#  Copyright (C) 1998 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     nbkal(response, times, mu, preg, pdepend, kalman=TRUE,
#	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
#	fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit a correlated negative binomial model with Kalman update
#  Adapted from Gauss code written by P. Lambert

nbkal <- function(response, times, mu, preg, pdepend, kalman=TRUE,
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p%*%p)){
likenb <- function(p){
	nu <- abs(p[length(p)])
	eta <- nu*exp(mu(p))
	sum(lgamma(response$response$y+1)+lgamma(eta)
		-lgamma(response$response$y+eta)+
		-eta*log(nu)+(response$response$y+eta)*log(1+nu))}
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
		for(jj in 1:response$response$nobs[ii]){
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
		for(jj in 1:response$response$nobs[ii]){
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
if(!inherits(response,"repeated")){
	if(!inherits(response,"response"))response <- restovec(response,times)
	response <- rmna(response=response)}
yy <- ifelse(response$response$y==0,1,response$response$y)
nind <- length(response$response$nobs)
p <- c(preg,pdepend)
np <- length(p)
np1 <- length(pdepend)
if(!kalman)z0 <- nlm(likenb, p, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
else z0 <- nlm(likekal, p, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
rpred <- if(kalman)likepred(z0$estimate) else NULL
if(any(is.na(z0$hessian)))a <- 0
else a <- qr(z0$hessian)$rank
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
	df=length(response$response$y)-np,
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

coefficients.nbkal <- function(z) z$coefficients
deviance.nbkal <- function(z) 2*z$maxlike
fitted.nbkal <- function(z, recursive=TRUE)
	if(recursive) z$rpred else z$pred
residuals.nbkal <- function(z, recursive=TRUE)
	if(recursive) z$response$y-z$rpred else z$response$y-z$pred

print.nbkal <- function(z, digits = max(3, .Options$digits - 3)) {
	np1 <- ifelse(z$full&z$kalman,3,1)
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$response$nobs),"\n")
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
	for(i in 1:nrow(coef.table))cname <- c(cname,paste("p",i,sep=""))
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
