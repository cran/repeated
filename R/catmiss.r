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
#     catmiss(response, frequency, ccov=NULL)
#
#  DESCRIPTION
#
#    A function to calculate marginal probabilities from the fitted
#  values of a log linear model

# Woolson and Clarke (1984) calculation of marginal probabilities

catmiss <- function(response, frequency, ccov=NULL){
	if(!is.matrix(response))stop("response must be a matrix")
	if(!is.vector(frequency,mode="numeric"))stop("frequency must be a vector")
	if(!is.null(ccov)){
		if(is.vector(ccov,mode="numeric")||is.vector(ccov,mode="character"))ccov <- matrix(ccov,ncol=1)
		else if(!is.matrix(ccov))stop("ccov must be a vector or matrix")}
	nr <- nrow(response)
	lf <- length(frequency)
	ncr <- ncol(response)
	response <- matrix(as.numeric(response),ncol=ncr)
	res <- sort(unique(response),na.last=T)
	nc <- length(res)-any(is.na(response))
	kron <- lf/nr
	if(trunc(kron)!=kron)stop("length of frequency must be a multiple of number of rows of response")
	if(!is.null(ccov)&&(nrow(ccov)*nr)!=lf)stop("ccov has incorrect number of rows")
# complete
	pc <- resp <- NULL
	for(i in 1:nrow(response))resp <-
		rbind(resp,if(any(is.na(response[i,])))rep(NA,ncr)
		else response[i,])
	for(i in 1:ncr){
	      tt <- as.matrix(tapply(frequency,list(rep(resp[,i],kron),gl(kron,nr,lf)),sum))
	      if(i==1)tot <- matrix(rep(rep(1,nc)%*%tt,nc),ncol=kron,byrow=T)
	      pc <- c(pc,as.vector(t(tt/tot)))}
	pc <- as.vector(matrix(pc,ncol=kron,byrow=T))
	coef.table <- cbind(pc,sqrt(pc*(1-pc)/rep(tot,rep(ncr,kron*nc))))
	cnames <- c("complete","se")
# all
if(any(is.na(response))){
	total <- rep(collapse(frequency,as.integer(gl(kron,nr,lf))),rep(nr,kron))
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
	if(is.null(ccov))rnames <- paste(rep(1:nc,kron*ncr),rep(rep(1:ncr,rep(nc,ncr)),kron))
	else {
	     rnames <- NULL
	     for(i in 1:nrow(ccov)){
		tmp <- NULL   
		for(j in 1:ncol(ccov))tmp <- paste(tmp,ccov[i,j])
		rnames <- c(rnames,tmp)}
	     rnames <- paste(rep(1:nc,kron*ncr),rep(1:ncr,rep(nc,ncr)),rep(rnames,rep(ncr*nc,length(rnames))))}
	dimnames(coef.table) <- list(rnames,cnames)
	coef.table
}
