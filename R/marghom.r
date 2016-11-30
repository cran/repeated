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
#     marg.hom(freq,marg1,marg2)
#
#  DESCRIPTION
#
#    A function to fit the marginal homogeneity model to a square
#  contingency table



##' Marginal Homogeneity Models
##' 
##' \code{marg.hom} fits a marginal homogeneity model to a contingency table
##' that has two margins of equal size.
##' 
##' 
##' @param freq Vector of frequencies
##' @param marg1 Factor variable for the first margin
##' @param marg2 Factor variable for the second margin
##' @return A list containing the call, the model, the deviance, the degrees of
##' freedom, the aic, the fitted values, and the residuals is returned.
##' @author J.K. Lindsey
##' @keywords models
##' @examples
##' 
##' # 4x4x2 table in Freq, with margins indexed by Left and Right
##' Freq <- rpois(32,10)
##' Left <- gl(4,1,32)
##' Right <- gl(4,4,32)
##' marg.hom(Freq, Left, Right)
##' 
##' @aliases marg.hom print.marginal
##' @export marg.hom
marg.hom <- function(freq,marg1,marg2){
call <- sys.call()
#
# check input and initialize for iteration
#
n <- length(unique(marg1))
if(length(unique(marg2))!=n)stop("a square contingency table must be supplied")
fit <- res <- freq
dv <- 0
test <- TRUE
lev <- levels(marg1)
a <- matrix(0,ncol=n-1,nrow=length(freq))
#
# iterate, calling glm with least squares
#
while(test){
	pw <- ifelse(fit>0,1/fit,0)
	for(i in 1:(n-1))a[,i] <- ((marg1==lev[i])-(marg2==lev[i]))*fit
	z0 <- glm(freq~a-1,weights=pw)
	fit <- freq-z0$fitted
	test <- (dv-z0$dev)^2>0.0001
	dv <- z0$dev}
z <- list(call=call,
	model=z0,
	deviance=2*sum(freq*log(ifelse(freq>0&fit>0,freq/fit,1))),
	df=n-1,
	aic=-sum(dpois(freq,fit,TRUE))+length(freq)-n+1,
	fitted=fit,
	residuals=(freq-fit)/sqrt(fit))
class(z) <- "marginal"
return(z)}

### print method
###
##' @export
print.marginal <- function(x, ...){
  z <- x
cat("Marginal homogeneity model\n")
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
cat("Deviance          ",z$deviance,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n\n")
cat("Parameter values\n")
cat(z$model$coef,"\n")}
