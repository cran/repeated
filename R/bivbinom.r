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
#     biv.binom(freq, marg1=~1, marg2=~1, interaction=~1, pmarg1=1,
#	pmarg2=1,pinteraction=1, print.level=0, typsize=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit a marginal bivariate binomial regression



##' Marginal Bivariate Binomial Regression Models
##' 
##' \code{biv.binom} fits (logit) linear regression models to a marginal
##' bivariate binomial distribution. The covariates must be of length K, that
##' is the number of 2x2 tables.
##' 
##' 
##' @param freq A four-column matrix containing K 2x2 frequency tables.
##' @param marg1 The model formula for the first margin.
##' @param marg2 The model formula for the second margin.
##' @param interaction The model formula for the interaction.
##' @param pmarg1 Initial parameter estimates for the first margin regression.
##' @param pmarg2 Initial parameter estimates for the second margin regression.
##' @param pinteraction Initial parameter estimates for the interaction
##' regression.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of class \code{bivbinom} is returned.
##' @author J.K. Lindsey
##' @keywords models
##' @examples
##' 
##' # 5 2x2 tables
##' Freq <- matrix(rpois(20,10),ncol=4)
##' x <- c(6,8,10,12,14)
##' print(z <- biv.binom(Freq,marg1=~x,marg2=~x,inter=~x,pmarg1=c(-2,0.08),
##' 	pmarg2=c(-2,0.1),pinter=c(3,0)))
##' 
##' @export biv.binom
biv.binom <- function(freq, marg1=~1, marg2=~1, interaction=~1, pmarg1=1,
	pmarg2=1,pinteraction=1, print.level=0, typsize=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
#
# set up likelihood function
#
like1 <- function(p){
	pia <- 1/(1+exp(dm1 %*% p[1:npt1]))
	pib <- 1/(1+exp(dm2 %*% p[(npt1+1):(npt1+npt2)]))
	psi <- exp(dm3 %*% p[(npt1+npt2+1):(npt1+npt2+npt3)])
	k1 <- 0.5/(psi-1)
	k2 <- 1+(pia+pib)*(psi-1)
	k3 <- 4*psi*(1-psi)*pia*pib
	k4 <- k1*(k2-sqrt(k2*k2+k3))
	s11 <- (psi!=1)*k4+(psi==1)*pia*pib
	s12 <- pia-s11
	s21 <- pib-s11
	s22 <- 1-s11-s12-s21
	ss[seq(1,n-3,by=4)] <- s11
	ss[seq(2,n-2,by=4)] <- s12
	ss[seq(3,n-1,by=4)] <- s21
	ss[seq(4,n,by=4)] <- s22
	list(like=-sum(freq*log(ss*(ss>0)+0.0001)),fitted=tot*ss,ss=ss)}
call <- sys.call()
#
# check input and prepare constants
#
if(!is.matrix(freq))stop("freq must be a matrix")
if(!dim(freq)[2]==4)stop("freq must have 4 columns")
freq <- as.vector(t(freq))
n <- length(freq)
nn <- n/4
tot <- rep(capply(freq,as.integer(gl(nn,4,n))),rep(4,nn))
ss <- rep(0,n)
#
# change formula for margins and interaction into matrices
#
p <- c(pmarg1,pmarg2,pinteraction)
if(inherits(marg1,"formula")){
	mt <- terms(marg1)
	if(is.numeric(mt[[2]])){
		dm1 <- matrix(rep(1,nn),ncol=1)
		colnames(dm1) <- "(Intercept)"
		npt1 <- 1}
	else {
		mf <- model.frame(mt,parent.frame(),na.action=na.fail)
		dm1 <- model.matrix(mt, mf)
		npt1 <- dim(dm1)[2]}
	np <- npt1
	nam1 <- colnames(dm1)}
else stop("marg1 must be a model formula")
if(npt1!=length(pmarg1))
	stop(paste(npt1,"parameter estimates must be supplied for marg1"))
if(inherits(marg2,"formula")){
	mt <- terms(marg2)
	if(is.numeric(mt[[2]])){
		dm2 <- matrix(rep(1,nn),ncol=1)
		colnames(dm2) <- "(Intercept)"
		npt2 <- 1}
	else {
		mf <- model.frame(mt,parent.frame(),na.action=na.fail)
		dm2 <- model.matrix(mt, mf)
		npt2 <- dim(dm2)[2]}
	np <- np+npt2
	nam2 <- colnames(dm2)}
else stop("marg2 must be a model formula")
if(npt2!=length(pmarg2))
	stop(paste(npt2,"parameter estimates must be supplied for marg2"))
if(inherits(interaction,"formula")){
	mt <- terms(interaction)
	if(is.numeric(mt[[2]])){
		dm3 <- matrix(rep(1,nn),ncol=1)
		colnames(dm3) <- "(Intercept)"
		npt3 <- 1}
	else {
		mf <- model.frame(mt,parent.frame(),na.action=na.fail)
		dm3 <- model.matrix(mt, mf)
		npt3 <- dim(dm3)[2]}
	np <- np+npt3
	nam3 <- colnames(dm3)}
else stop("interaction must be a model formula")
if(npt3!=length(pinteraction))
	stop(paste(npt3,"parameter estimates must be supplied for interaction"))
#
# call nlm to optimize
#
like1a <- function(p) like1(p)$like
z1 <- nlm(like1a, p=p, hessian=TRUE, print.level=print.level, typsize=typsize,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
fit <- like1(z1$estimate)$fitted
maxlike <- sum(fit-freq*log(fit)+lgamma(freq+1))
aic <- maxlike+length(p)+nn
#
# calculate se's
#
a <- if(any(is.na(z1$hessian))||any(abs(z1$hessian)==Inf))0
	else qr(z1$hessian)$rank
if(a==np)cov <- solve(z1$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
z <- list(call=call,
	maxlike=maxlike,
	aic=aic,
	df=n-length(p)-nn,
	fitted.values=fit,
	coefficients=z1$estimate,
	npt=list(npt1,npt2,npt3),
	vnames=list(nam1,nam2,nam3),
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z1$gradient,
	iterations=z1$iterations,
	code=z1$code)
class(z) <- "bivbinom"
return(z)}

### print method
#' @export
print.bivbinom <- function(x, ...){
  z <- x
	npt1 <- z$npt[[1]]
	npt2 <- z$npt[[2]]
	npt3 <- z$npt[[3]]
	np <- npt1+npt2+npt3
	cat("Bivariate binomial marginal regression\n")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	cat("Margin one parameters:\n")
	coef.table <- cbind(z$coefficients[1:npt1], z$se[1:npt1])
	dimnames(coef.table) <- list(z$vnames[[1]], c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	cat("\nMargin two parameters:\n")
	coef.table <- cbind(z$coefficients[(npt1+1):(npt1+npt2)], z$se[(npt1+1):(npt1+npt2)])
	dimnames(coef.table) <- list(z$vnames[[2]], c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	cat("\nInteraction parameters:\n")
	coef.table <- cbind(z$coefficients[(npt1+npt2+1):(npt1+npt2+npt3)], z$se[(npt1+npt2+1):(npt1+npt2+npt3)])
	dimnames(coef.table) <- list(z$vnames[[3]], c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr, digits=4)}
