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
#     hidden(response, totals=NULL, distribution="Bernoulli", pgamma, cmu=NULL,
#       tvmu=NULL, pcmu=NULL, ptvmu=NULL, pshape=NULL, pfamily=NULL,
#       delta=1, print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
#       fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit hidden Markov chain models in discrete time

hidden <- function(response, totals=NULL, distribution="Bernoulli", pgamma,
	cmu=NULL, tvmu=NULL, pcmu=NULL, ptvmu=NULL, pshape=NULL, pfamily=NULL,
	delta=1, print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	fscale=1, iterlim=100, typsiz=abs(p), stepmax=10*sqrt(p%*%p)){
likel <- function(p){
	if(is.function(cmu))pmu1 <- cmu(p[nm1:ncmu])
	if(is.function(tvmu))pmu2 <- tvmu(p[ncmu1:np1])
	if(np>np1){
		if(np3-np1==1)pshape <- rep(p[np3],states)
		else pshape <- p[np2:np3]}
	z <- .Fortran("hidden",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		nobs=as.integer(response$response$nobs),
		mobs=as.integer(mobs),
		s=as.double(response$response$y),
		n=as.integer(response$response$n),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		pos=as.integer(pos),
		gamma=double(states*states),
		model=as.integer(mdl),
		pcmu=as.double(pmu1),
		ptvmu=as.double(pmu2),
		pshape=as.double(pshape),
		pfam=as.double(p[np]),
		delta=double(states),
		nn=as.integer(nob),
		filter=double(states*nob),
		cf=as.logical(0),
		a=double(states),
		b=double(states*states),
		c=double(states),
		gmod=double(states*states),
		rhs=double(states),
		pivot=integer(states),
		qraux=double(states),
		work=double(2*states),
		like=double(1),
		DUP=F)
z$like}
like <- function(p){
	if(is.function(cmu))pmu1 <- cmu(p[nm1:ncmu])
	if(is.function(tvmu))pmu2 <- tvmu(p[ncmu1:np1])
	if(np>np1){
		if(np3-np1==1)pshape <- rep(p[np3],states)
		else pshape <- p[np2:np3]}
	z <- .Fortran("hidden",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		t=as.integer(response$response$nobs),
		mobs=as.integer(mobs),
		s=as.double(response$response$y),
		n=as.integer(response$response$n),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		pos=as.integer(pos),
		gamma=double(states*states),
		model=as.integer(mdl),
		pcmu=as.double(pmu1),
		ptvmu=as.double(pmu2),
		pshape=as.double(pshape),
		pfam=as.double(p[np]),
		delta=double(states),
		nn=as.integer(nob),
		filter=double(states*nob),
		cf=as.logical(1),
		a=double(states),
		b=double(states*states),
		c=double(states),
		gmod=double(states*states),
		rhs=double(states),
		pivot=integer(states),
		qraux=double(states),
		work=double(2*states),
		like=double(1),
		DUP=F)
z}
call <- sys.call()
tmp <- c("Bernoulli","Poisson","multinomial","binomial","exponential",
	"beta binomial","negative binomial","normal","inverse Gauss",
	"logistic","Cauchy","Laplace","Levy","Pareto","gamma","Weibull",
	"gen gamma","gen logistic","Hjorth","Burr","gen Weibull",
	"gen extreme value","gen inverse Gauss","power exponential")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
if(!is.matrix(pgamma)||ncol(pgamma)!=nrow(pgamma)){
	if(length(pgamma)==1&pgamma==1)pgamma <- as.matrix(pgamma)
	else stop("pgamma must be a square transition matrix")}
if(any(pgamma<0|pgamma>1))
     stop("All pgamma estimates must be between zero and one")
states <- nrow(pgamma)
if(mdl>5){
	if(is.null(pshape))stop("pshape estimate must be supplied")
	else if(length(pshape)!=1&&length(pshape)!=states)
	     stop(paste("pshape must have 1 or",states,"estimates"))}
if(mdl>16){
	if(is.null(pfamily))stop("pfamily estimate must be supplied")
	else if(length(pshape)!=1)
	     stop(paste("pshape must have one estimate"))}
if(states>1&any(diag(pgamma)==1))stop("Diagonal of pgamma cannot be 1")
pgamma <- ifelse(pgamma==0,1e-31,pgamma)
pg <- NULL
pos <- rep(1,states)
if(states>1)for(i in 1:states){
	if(sum(pgamma[i,])!=1)stop(paste("row",i,"of pgamma does not sum to 1"))
	for(j in 1:states){
		if(pgamma[i,j]>1e-30){
			pos[i] <- j
			pdg <- pgamma[i,j]
			break}}
	if(any(pgamma[i,]>=1e-30&pgamma[i,]!=1))
		pg <- c(pg,log(pgamma[i,(1:states)!=pos[i]&pgamma[i,]>=1e-30&pgamma[i,]!=1]/pdg))}
if(!inherits(response,"repeated")){
	if(!inherits(response,"response")){
		if(is.vector(response,mode="numeric"))response <- matrix(response,nrow=1)
		response <- restovec(response,totals=totals,delta=delta)}
	if(!is.null(response$delta))delta <- response$delta
	response <- rmna(response=response)}
if((distribution=="binomial"||distribution=="beta binomial")&&is.null(response$response$n))stop("totals must be supplied")
if(mdl<6){
	if(any(response$response$y<0))stop("all responses must be non-negative")}
else if((mdl!=8)&&(mdl!=10)&&(mdl!=11)&&(mdl!=12)&&(mdl!=13)&&(mdl!=18)&&(mdl!=24)&&(any(response$response$y<=0)))
	stop("all responses must be positive")
nosubj <- length(response$response$nobs)
mobs <- max(response$response$nobs)
nob <- length(response$response$y)
if(distribution=="multinomial"){
	if(min(response$response$y)<1)stop("multinomial categories must be numbered from 1")
	l <- max(response$response$y)-1
	if(l<1)stop("multinomial response must have at least 2 categories")}
else l <- 1
if(is.function(cmu)){
	if(is.null(pcmu))stop("Initial values of pcmu must be supplied")
	d <- dim(cmu(pcmu))
	if(length(d)<2)stop("cmu must return an array")
	if(d[2]!=states)stop(paste("cmu must return a",states,"column array"))
	if(d[1]!=nosubj)
		stop("cmu must return an array with one row per individual")
	if(distribution=="multinomial"){
		if(length(d)!=3)stop("cmu must return a 3 dimensional array")
		if(d[3]!=l)stop(paste("cmu must return an array with",l,"layers"))}}
else pmu1 <- array(0,c(states,nosubj,l))
if(is.function(tvmu)){
	if(is.null(ptvmu))stop("Initial values of ptvmu must be supplied")
	d <- dim(tvmu(ptvmu))
	if(length(d)<2)stop("tvmu must return an array")
	if(d[2]!=states)stop("tvmu must return a",states,"column array")
	else if(d[1]!=mobs)
	     stop("tvmu must return an array with one row per time point")
	if(distribution=="multinomial"){
		if(length(d)!=3)stop("tvmu must return a 3 dimensional array")
		if(d[3]!=l)stop(paste("tvmu must return an array with",l,"layers"))}}
else pmu2 <- array(0,c(states,mobs,l))
if(!is.function(cmu)&&!is.function(tvmu))
	stop("Either cmu or tvmu must be a function")
nm1 <- length(pg)+1
ncmu <- length(pg)+length(pcmu)
ncmu1 <- ncmu+1
ntvmu <- length(ptvmu)
p <- c(pg,pcmu,ptvmu,pshape,pfamily)
np1 <- length(p)-length(pshape)-length(pfamily)
np2 <- np1+1
np3 <- length(p)-length(pfamily)
np <- length(p)
pshape <- rep(0,states)
z <- nlm(likel,p, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
coef <- z$estimate[nm1:np1]
pshape <- if(np3>np1)z$estimate[np2:np3] else NULL
pfamily <- if(np>np3)z$estimate[np] else NULL
if(length(delta)==1)delta <- rep(delta,nob)
maxlike <- z$minimum - if(mdl>7) sum(log(delta)) else 0
z0 <- like(z$estimate)
if(any(is.na(z$hessian)))a <- 0
else a <- qr(z$hessian)$rank
if(a==np)cov <- solve(z$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
z1 <- list(
   call=call,
   distribution=distribution,
   response=response$response,
   maxlike=maxlike,
   aic=maxlike+np,
   states=states,
   cmu=cmu,
   tvmu=tvmu,
   coef=coef,
   pshape=pshape,
   pfamily=pfamily,
   ncmu=length(pcmu),
   ntvmu=length(ptvmu),
   gamma=matrix(z0$gamma,ncol=states),
   delta=z0$delta,
   cov=cov,
   corr=corr,
   se=se,
   filter=matrix(z0$filter,nrow=states),
   iterations=z$iter,
   code=z$code)
class(z1) <- "hidden"
z1}

print.hidden <- function(z, digits = max(3, .Options$digits - 3)){
	m <- z$states
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	cat(z$distribution)
	if(z$states>1)cat(" hidden Markov chain with",m,"states\n\n")
	else cat(" independence model\n\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$response$nobs),"\n")
	cat("Number of observations",length(z$response$y),"\n")
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n")
	if(z$ncmu>0){
		cat("\nTime-constant coefficients\n")
		t <- deparse(z$cmu)
		cat(t[2:length(t)],sep="\n")
		coef <- matrix(z$coef[1:z$ncmu],nrow=1)
		dimnames(coef) <- list(paste(1),paste(1:z$ncmu))
		print.default(coef, digits=digits, print.gap=2)}
	if(z$ntvmu>0){
		cat("\nTime-varying coefficients\n")
		t <- deparse(z$tvmu)
		cat(t[2:length(t)],sep="\n")
		coef <- matrix(z$coef[(z$ncmu+1):(z$ncmu+z$ntvmu)],nrow=1)
		dimnames(coef) <- list(paste(1),paste(1:z$ntvmu))
		print.default(coef, digits=digits, print.gap=2)}
	if(!is.null(z$pshape)){
		cat("\nDispersion parameters\n")
		coef <- matrix(z$pshape,nrow=1)
		dimnames(coef) <- list(paste(1),paste(1:length(z$pshape)))
		print.default(coef, digits=digits, print.gap=2)}
	if(!is.null(z$pfamily)){
		cat("\nFamily parameter\n")
		coef <- matrix(z$pfamily,nrow=1)
		dimnames(coef) <- list(paste(1),paste(1))
		print.default(coef, digits=digits, print.gap=2)}
	if(z$states>1){
		cat("\nTransition matrix\n")
		dimnames(z$gamma) <- list(paste(1:m),paste(1:m))
		print.default(z$gamma, digits=digits, print.gap=2)
		cat("\nStationary distribution\n")
		z$delta <- matrix(z$delta,nrow=1)
		dimnames(z$delta) <- list(paste(1),paste(1:m))
		print.default(z$delta, digits=digits, print.gap=2)}
}

plot.hidden <- function(z, nind=1, smooth=F, main=NULL, ylab=NULL,
	xlab="Time", xlim=NULL, ...){
	if(max(nind)>length(z$response$nobs))stop("no such individual")
	ns <- length(nind)
	ii <- covind(z$response)
	if(is.null(main)){
		main <- NULL
		for(i in nind)main <- c(main,paste("Individual ",i))}
	else if(length(main)!=ns){
		if(length(main==1))main <- rep(main,ns)
		else stop("main must have a name for each individual")}
	if(is.null(xlim))xlim <- c(0,max(z$resp$times))
	oldpar <- par(mfcol=c(z$states,ns))
	k <- 0
	for(i in nind){
		k <- k+1
		lenf <- sum(ii==i)
		for(j in 1:z$states){
			xl <- if(j==z$states)xlab else ""
			mn <- if(j==1)main[k] else ""
			yl <- if(k==1){
				if(is.null(ylab))paste("State ",j) else ylab}
			else ""
			if(smooth){
				y <- (z$filter[j,ii==i][1:(lenf-2)]+
					z$filter[j,ii==i][3:lenf])*0.25+
					z$filter[j,ii==i][2:(lenf-1)]*0.5
				y <- c(z$filter[j,ii==i][1],y,
					z$filter[j,ii==i][lenf])}
			else y <- z$filter[j,ii==i]
			plot(z$resp$times[ii==i],y,type="l", ylab=yl, xlab=xl,
				main=mn, ylim=c(0,1), xlim=xlim, ...)}}
	par(oldpar)}
