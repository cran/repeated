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
#     glmm(formula, family=gaussian, data=list(), weights=NULL,
#	offset=NULL, nest, delta=1, maxiter=20, points=10, print.level=0,
#	control=glm.control(epsilon=0.0001,maxit=10,trace=FALSE))
#
#  DESCRIPTION
#
#    A function to fit generalized linear mixed models with normal
#  random effect

# does not work: na.omit
glmm <- function(formula, family=gaussian, data=list(), weights=NULL,
	offset=NULL, nest, delta=1, maxiter=20, points=10, print.level=0,
	control=glm.control(epsilon=0.0001,maxit=10,trace=FALSE)){
	call <- sys.call()
	md <- missing(data)
	if(missing(data))data <- sys.frame(sys.parent())
	mf <- model.frame(terms(formula,data=data),data,na.action=na.fail)
	slen <- nrow(mf)
	if(is.character(family))family <- get(family)
	if(is.function(family))family <- family()
	if(is.vector(mf[,1])){
		y <- as.vector(mf[,1])
		if(family$family=="binomial")y2 <- rep(1,length(y))-y}
	else {
		y <- as.vector(mf[,1][,1])
		y2 <- as.vector(mf[,1][,2])}
	slen <- length(y)
	if(length(nest)!=slen)stop("the nesting variable is not the same length as the other variables")
	if(is.factor(nest))nest <- codes(nest)
	nind <- length(unique(nest))
	if(length(nest)==nind&&family$family!="binomial"&&family$family!="poisson")stop("Some individuals must have more than one observation")
	i <- rep(1:slen,points)
	ii <- rep(1:nind,points)
	k <- NULL
	for(j in 1:points)k <- c(k,nest+(j-1)*max(nest))
	k <- as.integer(k)
	quad <- gauss.hermite(points)
	sd <- quad[rep(1:points,rep(slen,points)),1]
	qw <- quad[rep(1:points,rep(nind,points)),2]
	nmodel <- update.formula(formula,.~.+sd)
	nnmf <- if(md)mf[i,,drop=F] else data[i,,drop=F]
	if(is.null(weights))lwt <- rep(1,slen)
	else lwt <- as.vector(weights)
	nobs <- sum(lwt)
	if(is.null(offset))offset <- rep(0,slen)
	else if(!is.data.frame(data))stop("offset can only be used when a data.frame is specified")
	if(is.data.frame(data)){
		if(is.null(weights))weights <- rep(1,slen)
		data <- cbind(data,weights,offset)
		zz <- glm(formula,family=family,data=data,control=control,weights=weights)}
	else {
		if(!is.null(weights))stop("weights can only be used when a data.frame is specified")
		zz <- glm(formula,family=family,data=data,control=control)}
	offset <- offset[i]
	ndev <- zz$deviance
	rdf <- zz$df.res-1
	ndf <- zz$df.null
	fv <- family$linkinv(zz$linear[i]+sd)
	sc <- ndev/nobs
	fpw <- switch(family$family,
		binomial= function()
			y*log(fv+0.0001)+(y2)*log(1-fv+0.0001),
		poisson= function() -fv+y*log(fv),
		Gamma= function() (log(y/fv)-y/fv)/sc,
		gaussian= function() -(y-fv)^2/sc/2,
		inverse.gaussian= function() -(y-fv)^2/(y*fv^2)/(2*sc))
	for(j in 1:maxiter){
		under  <- 0
		odev <- ndev
		ppr <- NULL
		for(ij in split(lwt*fpw(),k))ppr <- c(ppr,sum(ij))
		if(any(is.na(ppr)))stop("NAs - try another link function")
		if(max(ppr)-min(ppr)>1410){
			if(print.level==2)cat("Log probabilities:\n",ppr,"\n\n")
			stop("Product of probabilities is too small to calculate.\n Try fewer points.")}
		if(any(ppr > 705))under <- 705-max(ppr)
		else if(any(ppr < -705))under <- -705-min(ppr)
		pw <- qw*exp(ppr+under)
		pr <- NULL
		for(ij in split(pw,ii))pr <- c(pr,sum(ij))
		pw <- lwt*(pw/pr[ii])[k]
		nmf <- data.frame(nnmf,sd,pw,offset)
		z <- glm(nmodel,family=family,data=nmf,weights=pw,
			offset=offset,control=control)
		fv <- family$linkinv(z$linear)
		ndev <- -sum(log(pr)-under)
		sc <- z$dev/nobs
		if((odev-ndev)^2<0.00001)break}
	z$deviance <- -2*sum(log(pr)-under)
	formula <- update.formula(formula,.~1)
	class(formula) <- "formula"
	if(is.data.frame(data))
		zz <- glm(formula,family=family,data=data,control=control,weights=weights)
	else zz <- glm(formula,family=family,data=data,control=control)
	switch(family$family,
		binomial={
			sc <- NULL
			z$aic <- z$deviance-2*sum(lwt*lchoose(y+y2,y))
			z$deviance <- z$deviance+2*sum(lwt*(y*
				log(ifelse(y,y,1))+y2*log(ifelse(y2,y2,1))
				-(y+y2)*log(ifelse(y+y2,y+y2,1))))},
		poisson={
			sc <- NULL
			z$aic <- z$deviance+2*sum(lwt*lgamma(y+1))
			z$deviance <- z$deviance+2*sum(lwt*(-y+
				y*log(ifelse(y,y,1))))},
		Gamma={
			sc1 <- zz$dev/sum(lwt)
			z$null.deviance <- 2*sum(lwt*(log(sc1)/sc1-(log(y/zz$fit)-y/zz$fit)/sc1+log(y)+lgamma(1/sc1)-log(delta)))
			z$deviance <- z$deviance+2*sum(lwt*(log(sc)/sc+log(y)+lgamma(1/sc)-log(delta)))
			z$aic <- z$deviance+2},
		gaussian={
			z$null.deviance <- sum(lwt)*(log(2*pi*zz$dev/sum(lwt))+1)-2*sum(lwt*log(delta))
			z$deviance <- z$deviance+sum(lwt)*log(2*pi*sc)-2*sum(lwt*log(delta))
			z$aic <- z$deviance+2},
		inverse.gaussian={
			z$null.deviance <- sum(lwt*(log(2*pi*zz$dev/sum(lwt)*y^3))+1-2*log(delta))
			z$deviance <- z$deviance+sum(lwt*(log(2*pi*sc*y^3)-2*log(delta)))
			z$aic <- z$deviance+2})
	z$call <- call
	z$aic <- z$aic+2*z$qr$rank
	z$df.null <- ndf-!is.null(sc)
	z$df.residual <- rdf-!is.null(sc)
	z$iter <- j
	z$scale <- sc
	class(z) <- c("glmm",class(z))
	z}

print.glmm <- function(z,...){
	print.glm(z,...)
	if(!is.null(z$scale)){
		cat("\nModel deviances are -2 log likelihood\n")
		cat("Model dispersion:      ",z$scale,"\n")}
	cat("Normal mixing variance:",z$coef[names(z$coef)=="sd"]^2,"\n")}

print.summary.glmm <- function(z,...){
	print.summary.glm(z,...)
	if(!is.null(z$scale)){
		cat("Model deviances are -2 log likelihood\n")
		cat("Model dispersion:      ",z$scale,"\n")}
	cat("Normal mixing variance:",z$coef[rownames(z$coef)=="sd",1]^2,"\n")}

summary.glmm <- function(z,...){
	zz <- summary.glm(z,...)
	class(zz) <- c("summary.glmm",class(zz))
	if(!is.null(z$scale))zz$scale <- z$scale
	zz}
