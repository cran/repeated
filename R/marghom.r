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
#     marg.hom(freq,marg1,marg2)
#
#  DESCRIPTION
#
#    A function to fit the marginal homogeneity model to a square
#  contingency table

marg.hom <- function(freq,marg1,marg2){
	call <- sys.call()
	n <- length(unique(marg1))
	if(length(unique(marg2))!=n)stop("a square contingency table must be supplied")
	fit <- res <- freq
	dv <- 0
	test <- T
	a <- matrix(0,ncol=n-1,nrow=length(freq))
	while(test){
		pw <- 1/fit
		for(i in 1:(n-1))a[,i] <- ((marg1==i)-(marg2==i))*fit
		z0 <- glm(freq~a-1,weight=pw)
		fit <- freq-z0$fitted
		test <- (dv-z0$dev)^2>0.0001
		dv <- z0$dev}
	z <- list(call=call,
		model=z0,
		deviance=2*sum(freq*log(freq/fit)),
		df=n-1,
		aic=sum(fit-freq*log(fit)+lgamma(freq+1))+length(freq)-n+1,
		fitted=fit,
		residuals=(freq-fit)/sqrt(fit))
	class(z) <- "marginal"
	return(z)}

print.marginal <- function(z){
	cat("Marginal homogeneity model\n")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	cat("Deviance          ",z$deviance,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n\n")
	cat("Parameter values\n")
	cat(z$model$coef,"\n")}
