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
#     biv.betab(freq, x=NULL, p, depend=T, print.level=0,
#	typsiz=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit a bivariate beta-binomial regression

biv.betab <- function(freq, x=NULL, p, depend=T, print.level=0,
	typsiz=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
like <- function(p){
	dll <- 0
	for(kk in 1:n){
		tt1 <- p[1]*(y[kk,2]+y[kk,3]+2*y[kk,4])
		if(depend)tt1 <- tt1+p[2]*y[kk,4]
		if(full)for(i in 1:np)
			tt1 <- tt1+p[i+pn]*(x[kk,i]*(y[kk,2]+y[kk,4])+x[kk,i+np]*(y[kk,3]+y[kk,4]))
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
		dll <- dll-(y[kk,1]+y[kk,2]+y[kk,3]+y[kk,4])*log(tt2)}
	-dll}
call <- sys.call()
if(!is.matrix(y))stop("y must be a matrix")
else n <- dim(y)[1]
if(missing(x))np <- 0
else if(!is.matrix(x))stop("x must be a matrix")
else  {
	np <- dim(x)[2]/2
	if(trunc(dim(x)[2]/2)!=np)stop("x must contain an even number of columns")}
if(!dim(y)[2]==4)stop("y must have four columns")
if(!missing(x)&&!dim(y)[1]==dim(x)[1])stop("x and y must have the same number of rows")
full <- np>0
if(length(p)!=np+1+depend)stop(paste(np+1+depend,"parameter estimates must be supplied"))
pn <- depend+1
z0 <- nlm(like, p=p, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
np <- length(p)
a <- qr(z0$hessian)
if(a$rank==np)cov <- solve(z0$hessian)
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
		print.default(z$corr, digits=4)}
	invisible(z)}
