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
#     logitord(y, id, out.ccov=NULL, drop.ccov=NULL, tvcov=NULL,
#	out.tvcov=!is.null(tvcov), drop.tvcov=!is.null(tvcov),
#	pout, pdrop, prand.out, prand.drop,
#	random.out.int=T, random.out.slope=!is.null(tvcov),
#	random.drop.int=T, random.drop.slope=!is.null(tvcov),
#	binom.mix=5, fcalls=900, eps=0.0001, print.level=0)
#
#  DESCRIPTION
#
#    A function to fit binary or ordinal random effects models with dropouts.

logitord <- function(y, id, out.ccov=NULL, drop.ccov=NULL, tvcov=NULL,
	out.tvcov=!is.null(tvcov), drop.tvcov=!is.null(tvcov),
	pout, pdrop, prand.out, prand.drop,
	random.out.int=T, random.out.slope=!is.null(tvcov),
	random.drop.int=T, random.drop.slope=!is.null(tvcov),
	binom.mix=5, fcalls=900, eps=0.0001, print.level=0){
# Fortran constants
	maxsub <- 5200
	maxcas <- 10
	maxbet <- 25
	maxsig <- 10

	call <- sys.call()
	cg <- length(unique(y[y>0]))
	n1 <- length(y)
	n2 <- sum(y>0)
	if(!is.null(out.ccov)&&!is.matrix(out.ccov)){
		if(inherits(out.ccov,"formula")){
			mt <- terms(out.ccov)
			out.ccov <- model.matrix(mt,model.frame(mt,sys.frame(sys.parent()),na.action=na.fail))
			if(colnames(out.ccov)[1]=="(Intercept)")out.ccov <- out.ccov[,-1,drop=F]}
		else if(is.vector(out.ccov,mode="numeric")){
			frname <- paste(deparse(substitute(out.ccov)))
			out.ccov <- matrix(out.ccov,ncol=1)
			colnames(out.ccov) <- frname}
		else stop("out.ccov must be a vector, matrix, or model formula")}
	if(!is.null(out.ccov)&&nrow(out.ccov)!=n1)stop("y and out.ccov must have the same number of observations")
	if(!is.null(drop.ccov)&&!is.matrix(drop.ccov)){
		if(inherits(drop.ccov,"formula")){
			mt <- terms(drop.ccov)
			drop.ccov <- as.matrix(model.matrix(mt,model.frame(mt,sys.frame(sys.parent()),na.action=na.fail)))
			if(colnames(drop.ccov)[1]=="(Intercept)")drop.ccov <- drop.ccov[,-1,drop=F]}
		else if(is.vector(drop.ccov,mode="numeric")){
			frname <- paste(deparse(substitute(drop.ccov)))
			drop.ccov <- matrix(drop.ccov,ncol=1)
			colnames(drop.ccov) <- frname}
		else stop("drop.ccov must be a vector, matrix, or model formula")}
	if(!is.null(drop.ccov)&&nrow(drop.ccov)!=n1)stop("y and drop.ccov must have the same number of observations")
	if(!is.null(tvcov)){
		if(!is.vector(tvcov,mode="numeric"))stop("tvcov must be a vector")
		if(length(tvcov)!=n1)stop("y and tvcov must have the same number of observations")}
	else out.tvcov <- drop.tvcov <- F
	if(length(id)!=n1)stop("Every response value must have an id")
	if(length(unique(id))>maxsub)stop("Too many subjects")
	if(max(table(id))>maxcas)stop("Too many measurements for some subjects")
	if(!random.out.int||!out.tvcov)random.out.slope <- F
	if(!random.drop.int||!out.tvcov)random.drop.slope <- F
	id <- matrix(id,ncol=1)
	rownames(id) <- 1:nrow(id)
	y1 <- cbind(id,ifelse(y==0,2,1),rep(1,n1),drop.ccov)
	if(drop.tvcov)y1 <- cbind(y1,tvcov)
	y2 <- cbind(cbind(id,y)[y>0,],rep(0,n2))
	n3d <- 0
	if(!is.null(drop.ccov))n3d <- n3d+ncol(drop.ccov)
	if(drop.tvcov)n3d <- n3d+1
	if(n3d>0)for(i in 1:n3d)y2 <- cbind(y2,rep(0,n2))
	n3o <- 0
	if(!is.null(out.ccov))n3o <- n3o+ncol(out.ccov)
	if(out.tvcov)n3o <- n3o+1
	if(n3o>0)for(i in 1:n3o)y1 <- cbind(y1,rep(0,n1))
	ccov1 <- if(!is.null(out.ccov)) out.ccov[y>0] else NULL
	tvcov1 <- if(out.tvcov) tvcov[y>0] else NULL
	y2 <- cbind(y2,ccov1,tvcov1)
	if(random.out.int){
		y1 <- cbind(y1,rep(1,n1))
		y2 <- cbind(y2,rep(1,n2))}
	if(random.drop.int){
		y1 <- cbind(y1,rep(1,n1))
		y2 <- cbind(y2,rep(0,n2))}
	if(random.out.slope)y1 <- cbind(y1,tvcov,tvcov)
	if(random.drop.slope)y1 <- cbind(y1,tvcov,tvcov)
	for(i in 1:2){
		if(random.out.slope)y2 <- cbind(y2,tvcov[y>0])
		if(random.drop.slope)y2 <- cbind(y2,rep(0,n2))}
	outname <- if(cg==2) "(Intercept)"
		else c(paste("(Intercept)",1:(cg-1),sep=""))
	if(!is.null(out.ccov)){
		if(!is.null(colnames(out.ccov)))outname <- c(outname,colnames(out.ccov))
		else outname <- c(outname,paste("p",1:ncol(out.ccov),sep=""))}
	if(out.tvcov)outname <- c(outname,paste(deparse(substitute(tvcov))))
	dropname <- "(Intercept)"
	if(!is.null(drop.ccov)){
		if(!is.null(colnames(drop.ccov)))dropname <- c(dropname,colnames(drop.ccov))
		else dropname <- c(dropname,paste("p",1:ncol(drop.ccov),sep=""))}
	if(drop.tvcov)dropname <- c(dropname,paste(deparse(substitute(tvcov))))
	data <- rbind(y1,y2)
	data <- data[order(data[,1],rownames(data)),]
	total1 <- cg+n3o+n3d
	if(total1>maxbet)stop("too many regression parameters")
	total2a <- random.out.int+random.drop.int+random.out.slope+
		random.drop.slope
	total2b <- random.out.slope+random.drop.slope
	if(total2a>maxsig||total2b>maxsig)stop("too many random parameters")
	nobs <- nrow(data)
	if(missing(pout)||length(pout)!=cg+n3o-1)stop(paste(cg+n3o-1,"pout estimates must be supplied"))
	if(missing(pdrop)||length(pdrop)!=n3d+1)stop(paste(n3d+1,"pdrop estimates must be supplied"))
	p <- c(pout,pdrop)
	if(total2a+total2b>0){
		p <- c(p,rep(0.5,total2a+total2b))
		if(!missing(prand.out)){
			if(length(prand.out)!=random.out.int+2*random.out.slope)stop(paste(random.out.int+2*random.out.slope,"prand.out estimates must be supplied"))
			if(random.out.int)p[seq(total1+1,total1+total2a+total2b,by=2)] <- prand.out}
		if(!missing(prand.drop)){
			if(length(prand.drop)!=random.drop.int+2*random.drop.slope)stop(paste(random.drop.int+2*random.drop.slope,"prand.drop estimates must be supplied"))
			if(random.drop.int)p[seq(total1+random.out.int+1,total1+total2a+total2b,by=2)] <- prand.drop}}
	total <- total1+total2a+total2b
	z <- .Fortran("logitord",
		y=as.double(data),
		upk=as.integer(binom.mix),
		eps=as.double(eps),
		fcalls=as.integer(fcalls),
		iout=as.integer(print.level),
		cg=as.integer(cg),
		total1=as.integer(total1),
		total2a=as.integer(total2a),
		total2b=as.integer(total2b),
		nobs=as.integer(nobs),
		p=as.double(p),
		x=double(total),
		ster=double(total),
		hess=double(total*total),
		hessinv=double(total*total),
		nflag=integer(1),
		iter=integer(1),
		ifun=integer(1),
		like=double(1),
		DUP=F)
	if(z$nflag>0)switch(as.character(z$nflag),
			"1"=warning("Maximum number of function evaluations has been used"),
			"2"=stop("Linear search failed to improve the function value. Either the function or the gradient is incorrectly coded"),
			"3"=stop("Search vector was not a descent direction. The convergence criterion may be too strict"))
	z$df <- nobs-length(p)
	z$nobs <- nobs
	o <- 1:total1
	if(total2a+total2b>0&&total2a+total2b<=4){
		if(total2b==0||random.drop.slope||!random.drop.int)
			o <- c(o,(total1+1):(total1+total2a+total2b))
		else o <- c(o,total1+1,total1+3,total1+4,total1+2)}
	else if(total2a+total2b==6)o <- c(o,seq(total1+1,total1+total2a+total2b,by=2),seq(total1+2,total1+total2a+total2b,by=2))
	z$x <- z$x[o]
	z$corr <- (-z$hessinv/(z$ster%o%z$ster))[o,o]
	z$ster <- z$ster[o]
	z$outname <- outname
	z$dropname <- dropname
	z$cg <- cg
	z$n3o <- n3o
	z$n3d <- n3d
	z$random.out.int <- random.out.int
	z$random.drop.int <- random.drop.int
	z$random.out.slope <- random.out.slope
	z$random.drop.slope <- random.drop.slope
	z$call <- call
	class(z) <- "logitord"
	z}

print.logitord <- function(z){
	np <- z$total1+z$total2a+z$total2b
	cat("\nLogit ordinal dropout model\n")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	cat("\n-Log likelihood     ",z$like,"\n")
	cat("Degrees of freedom  ",z$df,"\n")
	cat("AIC                 ",z$like+length(z$x),"\n")
	cat("Iterations          ",z$iter,"\n")
	cat("Function evaluations",z$ifun,"\n")
	cat("\nOutcome model\n")
	num <- 1:(z$cg-1+z$n3o)
	cat("Fixed effect parameters\n")
	coef.table <- cbind(z$x[num],z$ster[num])
	dimnames(coef.table) <- list(z$outname, c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	if(z$random.out.int){
		cat("\nRandom effect parameters\n")
		num <- (z$total1+1):(z$total1+z$random.out.int+2*z$random.out.slope)
		coef.table <- cbind(z$x[num],z$ster[num])
		cname <- "s11"
		if(z$random.out.slope)cname <- c(cname,"s12","s22")
		dimnames(coef.table) <- list(cname,c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	cat("\nDropout model\n")
	cat("Fixed effect parameters\n")
	num <- (z$cg+z$n3o):z$total1
	coef.table <- cbind(z$x[num],z$ster[num])
	dimnames(coef.table) <- list(z$dropname, c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	if(z$random.drop.int){
		cat("\nRandom effect parameters\n")
		num <- (z$total1+z$random.out.int+2*z$random.out.slope+1):(z$total1+z$random.out.int+2*z$random.out.slope+z$random.drop.int+2*z$random.drop.slope)
		coef.table <- cbind(z$x[num],z$ster[num])
		cname <- "d11"
		if(z$random.drop.slope)cname <- c(cname,"d12","d22")
		dimnames(coef.table) <- list(cname,c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	cat("\nCorrelations\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr, digits=4)
}
