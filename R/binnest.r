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
#     binnnest(response, totals=NULL, nest=NULL, ccov=NULL, tvcov=NULL,
#	mu=~1, re1=~1, re2=~1, preg=NULL, pre1=NULL, pre2=NULL,
#	binom.mix=c(10,10), binom.prob=c(0.5,0.5), fcalls=900,
#	eps=0.01, print.level=0)
#
#  DESCRIPTION
#
#    A function to fit binary random effects models with two levels
#  of nesting.
#
binnest <- function(response, totals=NULL, nest=NULL, ccov=NULL, tvcov=NULL,
	mu=~1, re1=~1, re2=~1, preg=NULL, pre1=NULL, pre2=NULL,
	binom.mix=c(10,10), binom.prob=c(0.5,0.5), fcalls=900,
	eps=0.01, print.level=0){
# Fortran constants
maxt1 <- maxt2 <- maxt3 <- 10

call <- sys.call()
if(length(binom.prob)==1)binom.prob <- c(binom.prob,binom.prob)
if(any(binom.prob<=0)||any(binom.prob>=1))stop("binom.prob parameters must be between zero and one")
if(length(binom.mix)==1)binom.mix <- c(binom.mix,binom.mix)
total1 <- length(preg)
total2 <- length(pre1)
total3 <- length(pre2)
total <- total1+total2+total3
dimw <- total*(total+7)/2
if(!inherits(response,"repeated")){
	if(!inherits(response,"response")){
		if(is.matrix(response)||is.data.frame(response))response <- restovec(response,totals=totals,nest=nest,times=F)
		else if(is.list(response))response <- restovec(response,times=F)
		else stop("response must be a matrix, data.frame, list, or object of type repeated or response")}
	resp <- response$y
	if(is.null(response$n)){
		if(any(response$y!=0&response$y!=1))stop("if binomial totals are not supplied, all responses must be 0 or 1")
		else resp <- cbind(resp,rep(1,length(response$y)))}
	else resp <- cbind(resp,response$n)
	if(!is.null(ccov)){
		if(!inherits(ccov,"tccov"))ccov <- tcctomat(ccov)
		resp <- cbind(resp,ccov$ccov[covind(response),,drop=F])}
	if(!is.null(tvcov)){
		if(!inherits(tvcov,"tvcov"))tvcov <- tvctomat(tvcov)
		resp <- cbind(resp,tvcov$tvcov)}
	nind <- length(response$nobs)
	numsubj <- length(response$y)
	nest <- response$nest
	regname <- if(ncol(resp)>2)c("(Intercept)",colnames(resp[,3:ncol(resp)]))
	else "(Intercept)"
	if(total2==1)re1name <- "(Intercept)"
	else if(!is.null(total2))stop("pre1 should supply 1 initial estimate")
	if(total3==1)re2name <- "(Intercept)"
	else if(!is.null(total3))stop("pre2 should supply 1 initial estimate")
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
else {
	envname <- paste(deparse(substitute(response)))
	resp <- response$response$y
	if(is.null(response$response$n)){
		if(any(response$response$y!=0&response$response$y!=1))stop("if binomial totals are not supplied, all responses must be 0 or 1")
		else resp <- cbind(resp,rep(1,length(response$response$y)))}
	else resp <- cbind(resp,response$response$n)
	if(inherits(mu,"formula")){
		if(as.character(mu)[2]!="1"){
			tmp <- attr(finterp(mu,envir=response,name=envname),"model")
			if(!is.matrix(tmp))stop("mu must be a W&R formula")
			if(ncol(tmp)!=total1)stop(paste("preg should contain",ncol(tmp),"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=F])
			regname <- gsub("\\[.i]","",colnames(tmp))}
		else {
			regname <- "(Intercept)"
			if(total1!=1)stop("preg should contain 1 initial estimate")}}
	else stop("mu must be a W&R formula")
	if(inherits(re1,"formula")&&!is.null(pre1)){
		if(as.character(re1)[2]!="1"){
			tmp <- attr(finterp(re1,envir=response,name=envname),"model")
			if(!is.matrix(tmp))stop("re1 must be a W&R formula")
			if(ncol(tmp)!=total2)stop(paste("pre1 should contain",ncol(tmp),"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=F])
			re1name <- gsub("\\[.i]","",colnames(tmp))}
		else {
			re1name <- "(Intercept)"
			if(total2!=1)stop("pre1 should contain 1 initial estimate")}}
	else pre1 <- re1name <- NULL
	if(inherits(re2,"formula")&&!is.null(pre2)){
		if(as.character(re2)[2]!="1"){
			tmp <- attr(finterp(re2,envir=response,name=envname),"model")
			if(!is.matrix(tmp))stop("re2 must be a W&R formula")
			if(ncol(tmp)!=total3)stop(paste("pre2 should contain",ncol(tmp),"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=F])
			re2name <- gsub("\\[.i]","",colnames(tmp))}
		else {
			re2name <- "(Intercept)"
			if(total3!=1)stop("pre2 should contain 1 initial estimate")}}
	else pre2 <- re2name <- NULL
	nind <- length(response$response$nobs)
	numsubj <- length(response$response$y)
	nest <- response$response$nest}
nobs1 <- nobs2 <- NULL
if(!is.null(nest))for(i in 1:nind){
	nobs1 <- c(nobs1,length(unique(nest[covind(response)==i])))
	nobs2 <- c(nobs2,as.vector(table(nest[covind(response)==i])))}
else nobs1 <- nobs2 <- rep(1,nind)
rm(response)
maxmother <- max(nobs1)
maxkid <- max(nobs2)
p <- c(preg,pre1,pre2)
z0 <- .Fortran("binnest",
        Fvalue=double(1),
        res=integer(3),  #  Iter_N,Fun_N,flag
        x=double(total),
        g=double(total),
        hess=double(total*total),
        p=as.double(p),
        numcase1=as.integer(nind),
	numcase2=as.integer(length(nobs2)),
	numsubj=as.integer(numsubj),
	maxkid=as.integer(maxkid),
	maxmother=as.integer(maxmother),
        totcol=as.integer(ncol(resp)),
        total1=as.integer(total1),
	total2=as.integer(total2),
	total3=as.integer(total3),
	uph1in=as.integer(binom.mix[1]),
	uph2in=as.integer(binom.mix[2]),
	fcalls=as.integer(fcalls),
	dimw=as.integer(dimw),
        par=as.double(c(eps,binom.prob[1:2])),
	case1=as.integer(nobs1),
        case2=as.integer(nobs2),
        subject=as.double(resp),  # numsubj*(total+2)
	iout=as.integer(print.level),
        hab=double(maxmother*total1),
        hac=double(maxmother*total2),
        had=double(maxmother*total3),
        ha=double(maxmother),
        v1=double(binom.mix[1]),
        v2=double(binom.mix[2]),
        h1choo=double(binom.mix[1]),
        h2choo=double(binom.mix[2]),
        hn=double(binom.mix[2]),
        h1=double(binom.mix[1]),
        h2=double(binom.mix[2]),
        betakk=double(maxmother*maxkid),
        sig1kk=double(maxmother*maxkid),
        sig2kk=double(maxmother*maxkid),
        mother=integer(maxmother),
        rr=double(nind*maxmother*maxkid),
        r=double(nind*maxmother*maxkid),
        sn=double(nind*maxmother*maxkid),
        z=double(nind*maxmother*maxkid*total1),
        uu1=double(nind*maxmother*maxkid*total2),
        uu2=double(nind*maxmother*maxkid*total3),
        w=double(dimw),
        habb=double(maxmother*total1*total1),
        habs1=double(maxmother*total1*total2),
        habs2=double(maxmother*total1*total3),
        has1s1=double(maxmother*total2*total2),
        has1s2=double(maxmother*total2*total3),
        has2s2=double(maxmother*total3*total3),
        ebb=double(maxmother*total1*total1),
        ebs1=double(maxmother*total1*total2),
        ebs2=double(maxmother*total1*total3),
        fs1s1=double(maxmother*total2*total2),
        fs1s2=double(maxmother*total2*total3),
        gs2s2=double(maxmother*total3*total3),
        e2bb=double(maxmother*maxmother*total1*total1),
        e2bs1=double(maxmother*maxmother*total1*total2),
        e2bs2=double(maxmother*maxmother*total1*total3),
        f2s1s1=double(maxmother*maxmother*total2*total2),
        f2s1s2=double(maxmother*maxmother*total2*total3),
        g2s2s2=double(maxmother*maxmother*total3*total3),
	dup=F)
if(z0$res[3]>0)switch(as.character(z0$res[3]),
		"1"=warning("Maximum number of function evaluations has been used"),
		"2"=stop("Linear search failed to improve the function value. Either the function or the gradient is incorrectly coded"),
		"3"=stop("Search vector was not a descent direction. The convergence criterion may be too strict"))
z0$hess <- matrix(-z0$hess,ncol=total)
if(any(is.na(z0$hess)))a <- 0
else a <- qr(z0$hess)$rank
if(a==total)cov <- solve(z0$hess)
else cov <- matrix(NA,ncol=total,nrow=total)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
z <- list(
        call=call,
	mu=mu,
	re1=re1,
	re2=re2,
        maxlike=z0$Fvalue,
        aic=z0$Fvalue+total,
        df=nrow(resp)-total,
	total1=total1,
	total2=total2,
	total3=total3,
        coefficients=z0$x,
	regname=regname,
	re1name=re1name,
	re2name=re2name,
        se=se,
        cov=cov,
        corr=corr,
        grad=z0$g,
        iterations=z0$res[1],
	ifun=z0$res[2],
        code=z0$res[3])
class(z) <- "binnest"
return(z)}

print.binnest <- function(z){
	np <- z$total1+z$total2+z$total3
	cat("\nNested binomial model\n")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n-Log likelihood     ",z$maxlike,"\n")
	cat("Degrees of freedom  ",z$df,"\n")
	cat("AIC                 ",z$aic,"\n")
	cat("Iterations          ",z$iter,"\n")
	cat("Function evaluations",z$ifun,"\n")
	cat("\nFixed effect parameters\n")
	coef.table <- cbind(z$coef[1:z$total1],z$se[1:z$total1])
	dimnames(coef.table) <- list(z$regname, c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	if(z$total2>0){
		cat("\nFirst level random effects parameters\n")
		num <- (z$total1+1):(z$total1+z$total2)
		coef.table <- cbind(z$coef[num],z$se[num])
		dimnames(coef.table) <- list(z$re1name,c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	if(z$total3>0){
		cat("\nSecond level random effects parameters\n")
		num <- (z$total1+z$total2+1):np
		coef.table <- cbind(z$coef[num],z$se[num])
		dimnames(coef.table) <- list(z$re2name,c("estimate", "se"))
		print.default(coef.table, digits=4, print.gap=2)}
	cat("\nCorrelations\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr, digits=4)
}
