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
#     catmiss(response, frequency, ccov=NULL)
#
#  DESCRIPTION
#
#    A function to calculate marginal probabilities from the fitted
#  values of a log linear model

# Woolson and Clarke (1984) calculation of marginal probabilities



##' Marginal Probabilities for Categorical Repeated Measurements with Missing
##' Data
##' 
##' \code{catmiss} calculates the marginal probabilities of repeated responses.
##' If there are missing values, it gives both the complete data estimates and
##' the estimates using all data. It is useful, for example, when a log linear
##' model is fitted; the resulting fitted values can be supplied to
##' \code{catmiss} to obtain the estimates of the marginal probabilities for
##' the model. (Note however that the standard errors do not take into account
##' the fitting of the model.)
##' 
##' 
##' @param response A matrix with one column for each of the repeated measures
##' and one row for each possible combination of responses, including the
##' missing values, indicated by NAs.
##' @param frequency A vector containing the frequencies. Its length must be a
##' multiple of the number of rows of \code{response}. Responses are arranged
##' in blocks corresponding to the various possible combinations of values of
##' the explanatory variables.
##' @param ccov An optional matrix containing the explanatory variables
##' (time-constant covariates) as columns, with one line per block of responses
##' in \code{frequency}. Thus, the number of rows of response times the number
##' of rows of \code{ccov} equals the length of \code{frequency}.
##' @return A matrix with the probabilities and their standard errors is
##' returned.
##' @author J.K. Lindsey
### @seealso \code{\link{glm}}, \code{\link[gnlm]{nordr}}
##' @keywords models
##' @examples
##' 
##' y <- rpois(27,15)
##' r1 <- gl(3,1,27)
##' r2 <- gl(3,3,27)
##' r3 <- gl(3,9)
##' # r1, r2, and r3 are factor variables with 3 indicating missing
##' # independence model with three binary repeated measures
##' # with missing values
##' print(z <- glm(y~r1+r2+r3, family=poisson))
##' # obtain marginal estimates (no observations with 3 missing values)
##' resp <- cbind(as.integer(r1), as.integer(r2), as.integer(r3))[1:26,]
##' resp <- ifelse(resp==3, NA, resp)
##' catmiss(resp, y[1:26])
##' 
##' @export catmiss
catmiss <- function(response, frequency, ccov=NULL){
#
# check that input is valid
#
if(!is.matrix(response))stop("response must be a matrix")
if(!is.vector(frequency,mode="numeric"))
	stop("frequency must be a vector")
if(!is.null(ccov)){
	if(is.vector(ccov,mode="numeric")||is.vector(ccov,mode="character"))ccov <- matrix(ccov,ncol=1)
	else if(!is.matrix(ccov))stop("ccov must be a vector or matrix")}
#
# calculate constants
#
nr <- dim(response)[1]
lf <- length(frequency)
ncr <- dim(response)[2]
response <- matrix(as.numeric(response),ncol=ncr)
res <- sort(unique(as.vector(response)),na.last=TRUE)
nc <- length(res)-any(is.na(response))
kron <- lf/nr
if(trunc(kron)!=kron)
	stop("length of frequency must be a multiple of number of rows of response")
if(!is.null(ccov)&&(dim(ccov)[1]*nr)!=lf)
	stop("ccov has incorrect number of rows")
#
# calculations for complete observations
#
pc <- resp <- NULL
for(i in 1:dim(response)[1])
	resp <- rbind(resp,if(any(is.na(response[i,])))rep(NA,ncr)
		else response[i,])
for(i in 1:ncr){
      tt <- as.matrix(tapply(frequency,list(rep(resp[,i],kron),gl(kron,nr,lf)),sum))
      if(i==1)tot <- matrix(rep(rep(1,nc)%*%tt,nc),ncol=kron,byrow=TRUE)
      pc <- c(pc,as.vector(t(tt/tot)))}
pc <- as.vector(matrix(pc,ncol=kron,byrow=TRUE))
coef.table <- cbind(pc,sqrt(pc*(1-pc)/rep(tot,rep(ncr,kron*nc))))
cnames <- c("complete","se")
#
# calculations for all observations
#
if(any(is.na(response))){
	total <- rep(capply(frequency,as.integer(gl(kron,nr,lf))),rep(nr,kron))
	prmat <- NULL
	b <- matrix(0,nrow=ncr*nc,ncol=ncr*(nc+1))
	for(i in 1:ncr){
	      jj <- 0
	      for(j in res[1:nc]){
		    jj <- jj+1
		    prmat <- rbind(prmat,as.numeric(!is.na(response[,i])&response[,i]==j))
		    b[jj+(i-1)*nc,jj+(i-1)*(nc+1)] <- 1
		    b[jj+(i-1)*nc,i*(nc+1)] <- -1}
	      prmat <- rbind(prmat,as.numeric(!is.na(response[,i])))}
	p1 <- diag(kron)%x%b
	p2 <- diag(kron)%x%prmat
	p3 <- frequency/total
	p4 <- as.vector(p2%*%p3)
	p <- as.vector(exp(p1%*%log(p4)))
	se <- diag(p)%*%p1%*%diag(1/p4)%*%p2
	coef.table <- cbind(coef.table,p,sqrt(diag(se%*%diag(p3*(1-p3)/total)%*%t(se))))
	cnames <- c(cnames,"all","se")}
#
# set up table to return
#
if(is.null(ccov))
	rnames <- paste(rep(1:nc,kron*ncr),rep(rep(1:ncr,rep(nc,ncr)),kron))
else {
	rnames <- NULL
	for(i in 1:dim(ccov)[1]){
		tmp <- NULL   
		for(j in 1:dim(ccov)[2])tmp <- paste(tmp,ccov[i,j])
		rnames <- c(rnames,tmp)}
	rnames <- paste(rep(1:nc,kron*ncr),rep(1:ncr,rep(nc,ncr)),rep(rnames,rep(ncr*nc,length(rnames))))}
dimnames(coef.table) <- list(rnames,cnames)
coef.table}
