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
#     chidden(response, totals=NULL, times=NULL, distribution="Bernoulli",
#       pgamma, cmu=NULL, tvmu=NULL, pcmu=NULL, ptvmu=NULL, pshape=NULL,
#       pfamily=NULL, delta=1, print.level=0, ndigit=10, gradtol=0.00001,
#       steptol=0.00001, fscale=1, iterlim=100, typsiz=abs(p),
#	stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit hidden Markov chain models in continuous time

chidden <- function(response, totals=NULL, times=NULL,
	distribution="Bernoulli", pgamma, cmu=NULL,
	tvmu=NULL, pcmu=NULL, ptvmu=NULL, pshape=NULL,
	pfamily=NULL, delta=1, print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, fscale=1, iterlim=100, typsiz=abs(p),
	stepmax=10*sqrt(p%*%p)){
likel <- function(p){
	if(is.function(cmu))pmu1 <- cmu(p[nm1:ncmu])
	if(is.function(tvmu))pmu2 <- tvmu(p[ncmu1:np1])
	if(np>np1){
		if(np3-np1==1)pshape <- rep(p[np3],states)
		else pshape <- p[np2:np3]}
	z <- .Fortran("chidden",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		nobs=as.integer(response$response$nobs),
		mobs=as.integer(mobs),
		s=as.double(response$response$y),
		n=as.integer(response$response$n),
		times=as.double(response$response$times),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		gamma=double(states*states),
		val=double(states),
		vec=double(states*states),
		invec=double(states*states),
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
	z <- .Fortran("chidden",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		t=as.integer(response$response$nobs),
		mobs=as.integer(mobs),
		s=as.double(response$response$y),
		n=as.integer(response$response$n),
		times=as.double(response$response$times),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		gamma=double(states*states),
		val=double(states),
		vec=double(states*states),
		invec=double(states*states),
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
	z$gamma <- matrix(0,ncol=states,nrow=states)
	z$gamma[ipos] <- p[1:nm]
	diag(z$gamma) <- -z$gamma%*%rep(1,states)
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
states <- nrow(pgamma)
if(mdl>5){
	if(is.null(pshape))stop("pshape estimate must be supplied")
	else if(length(pshape)!=1&&length(pshape)!=states)
	     stop(paste("pshape must have 1 or",states,"estimates"))}
if(mdl>16){
	if(is.null(pfamily))stop("pfamily estimate must be supplied")
	else if(length(pshape)!=1)
	     stop(paste("pshape must have one estimate"))}
ipos <- pg <- NULL
if(states>1)for(i in 1:states){
	if(pgamma[i,i]!=0&sum(pgamma[i,])!=0)stop(paste("row",i,"of pgamma does not sum to 0"))
	for(j in 1:states)if(i!=j&&pgamma[i,j]!=0)
		ipos <- c(ipos,states*(j-1)+i)
	pg <- c(pg,pgamma[i,(1:states)!=i&pgamma[i,]!=0])}
if(!inherits(response,"repeated")){
	if(!inherits(response,"response")){
		if(is.vector(response,mode="numeric"))response <- matrix(response,nrow=1)
		response <- restovec(response,totals=totals,times=times,delta=delta)}
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
nm <- length(pg)
nm1 <-  nm+1
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
coef <- z$estimate[nm1:np]
pshape <- if(np3>np1)z$estimate[np2:np3] else NULL
pfamily <- if(np>np3)z$estimate[np] else NULL
maxlike <- z$minimum
delta <- if(!is.null(response$response$delta)) sum(log(response$response$delta)) else 0
maxlike <- z$minimum - if(mdl>7) delta else 0
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
