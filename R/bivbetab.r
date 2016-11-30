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
#     biv.betab(freq, x=NULL, p, depend=TRUE, print.level=0,
#	typsize=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit a bivariate beta-binomial regression



##' Bivariate Beta-binomial Regression Models
##' 
##' \code{biv.betab} fits dependent (logit) linear regression models to a
##' bivariate beta-binomial distribution.
##' 
##' 
##' @param freq A matrix containing four columns corresponding to 00, 01, 10,
##' and 11 responses.
##' @param x A matrix of explanatory variables, containing pairs of columns,
##' one for each response, and the same number of rows as freq.
##' @param p Initial parameter estimates: intercept, dependence (if depend is
##' TRUE, and one for each pair of columns of x.
##' @param depend If FALSE, the independence (logistic) model is fitted.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of class \code{bivbetab} is returned.
##' @author J.K. Lindsey
##' @keywords models
##' @examples
##' 
##' y <- matrix(  c( 2, 1, 1,13,
##' 		 4, 1, 3, 5,
##' 		 3, 3, 1, 4,
##' 		15, 8, 1, 6),ncol=4,byrow=TRUE)
##' first <- c(0,0,1,1)
##' second <- c(0,1,0,1)
##' self <- cbind(first,second)
##' other <- cbind(second,first)
##' biv.betab(y,cbind(self,other),p=c(-1,2,1,1))
##' # independence
##' biv.betab(y,cbind(self,other),p=c(-1,1,1),dep=FALSE)
##' 
##' @export biv.betab
biv.betab <- function(freq, x=NULL, p, depend=TRUE, print.level=0,
	typsize=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
#
# set up likelihood function
#
like <- function(p){
	dll <- 0
	for(kk in 1:n){
		tt1 <- p[1]*(freq[kk,2]+freq[kk,3]+2*freq[kk,4])
		if(depend)tt1 <- tt1+p[2]*freq[kk,4]
		if(full)for(i in 1:np)
			tt1 <- tt1+p[i+pn]*(x[kk,i]*(freq[kk,2]+freq[kk,4])+x[kk,i+np]*(freq[kk,3]+freq[kk,4]))
		dll <- dll+tt1
		tt2 <- 1
		t1 <- t2 <- p[1]
		if(full)for(i in 1:np){
			t1 <- t1+p[i+pn]*x[kk,i]
			t2 <- t2+p[i+pn]*x[kk,i+np]}
		tt2 <- tt2+exp(t1)+exp(t2)
		t1 <- 2*p[1]
		if(depend)t1 <- t1+p[2]
		if(full)for(i in 1:np)t1 <- t1+p[i+pn]*(x[kk,i]+x[kk,i+np])
		tt2 <- tt2+exp(t1)
		dll <- dll-(freq[kk,1]+freq[kk,2]+freq[kk,3]+freq[kk,4])*log(tt2)}
	-dll}
call <- sys.call()
#
# check that correct data were supplied
#
if(!is.matrix(freq))stop("freq must be a matrix")
else n <- dim(freq)[1]
if(missing(x))np <- 0
else if(!is.matrix(x))stop("x must be a matrix")
else  {
	np <- dim(x)[2]/2
	if(trunc(dim(x)[2]/2)!=np)
		stop("x must contain an even number of columns")}
if(!dim(freq)[2]==4)stop("freq must have four columns")
if(!is.null(x)&&!dim(freq)[1]==dim(x)[1])
	stop("x and freq must have the same number of rows")
full <- np>0
if(length(p)!=np+1+depend)
	stop(paste(np+1+depend,"parameter estimates must be supplied"))
#
# optimize with nlm
#
pn <- depend+1
z0 <- nlm(like, p=p, hessian=TRUE, print.level=print.level, typsize=typsize,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
np <- length(p)
#
# calculate se's
#
a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
	else qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
z1 <- list(
	call=call,
	maxlike=z0$minimum,
	aic=z0$minimum+length(p),
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "bevbetab"
return(z1)}

### print method
###
print.bevbetab <- function(z) {
np <- length(z$coef)
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("-Log likelihood   ",z$maxlike,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Parameter estimates:\n")
coef.table <- cbind(z$coef,z$se)
dimnames(coef.table) <- list(seq(1,np), c("estimate", "se"))
print.default(coef.table, digits=4, print.gap=2)
if(np>1){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr, digits=4)}}
