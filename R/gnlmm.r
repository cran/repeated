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
#     gnlmm(y, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
#	nest=NULL, pmu=NULL, pshape=NULL, psd=NULL, exact=F, wt=1,
#	delta=1, shfn=F, scale=NULL, points=10, envir=sys.frame(sys.parent()),
#	print.level=0, typsiz=abs(p), ndigit=10, gradtol=0.00001,
#	stepmax=sqrt(p%*%p)/10, steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit generalized nonlinear mixed models with normal
#  random effect

gnlmm <- function(y, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
	nest=NULL, pmu=NULL, pshape=NULL, psd=NULL, exact=F, wt=1, delta=1,
	shfn=F, scale=NULL, points=10, envir=sys.frame(sys.parent()),
	print.level=0, typsiz=abs(p), ndigit=10, gradtol=0.00001,
	stepmax=sqrt(p%*%p)/10, steptol=0.00001, iterlim=100, fscale=1){

pinvgauss <- function(y,m,s){
	t <- y/m
	v <- sqrt(y*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
plaplace <- function(y,m,s){
	u <- (y-m)/s
	t <- exp(-abs(u))/2
	ifelse(u<0,t,1-t)}
plevy <- function(y, m, s)
	.C("plevy",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(1),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		DUP=F)$res

call <- sys.call()
if(!missing(distribution)&&!is.function(distribution)){
	distribution <- match.arg(distribution,c("binomial","beta binomial",
	"double binomial","mult binomial","Poisson","negative binomial",
	"double Poisson","mult Poisson","gamma count","Consul","logarithmic",
	"geometric","normal","inverse Gauss","logistic","exponential","gamma",
	"Weibull","extreme value","Pareto","Cauchy","Laplace","Levy"))}
shp <- distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"
if(!missing(scale))scale <- match.arg(scale,c("identity","log","logit",
	"reciprocal","exp"))
if(!missing(pmu))npl <- length(pmu)
else npl <- 0
if(!missing(pshape))nps <- length(pshape)
else nps <- 0
if(missing(psd))stop("An initial value of psd must be supplied")
np <- npl+nps+1
if(np<1)stop("At least one parameter must be estimated")
if(is.function(distribution)){
	fcn <- distribution
	distribution <- "own"}
if(inherits(y,"repeated")){
	if(inherits(envir,"repeated")&&(length(y$response$nobs)!=length(envir$response$nobs)||any(y$response$nobs!=envir$response$nobs)))stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&!is.na(y$response$wt))wt <- y$response$wt
	if(!is.null(y$response$delta))delta <- y$response$delta
	if(length(y$response$nobs)==1&&y$response$nobs==1)
		nest <- 1:length(y$response$y)
	else nest <- covind(y)
	if(!is.null(y$response$censor)&&any(y$response$censor!=1))y <- cbind(y$response$y,y$response$censor)
	else if(!is.null(y$response$n))y <- cbind(y$response$y,y$response$n-y$response$y)
	else y <- y$response$y}
else if(inherits(y,"response")){
	if(inherits(envir,"repeated")&&(length(y$nobs)!=length(envir$response$nobs)||any(y$nobs!=envir$response$nobs)))stop("y and envir objects are incompatible")
	if(!is.null(y$wt)&&!is.na(y$wt))wt <- y$wt
	if(!is.null(y$delta))delta <- y$delta
	if(length(y$nobs)==1&&y$nobs==1)
		nest <- 1:length(y$y)
	else nest <- covind(y)
	if(!is.null(y$censor)&&any(y$censor!=1))y <- cbind(y$y,y$censor)
	else if(!is.null(y$n))y <- cbind(y$y,y$n-y$y)
	else y <- y$y}
if(distribution=="binomial"||distribution=="double binomial"||
	distribution=="beta binomial"||distribution=="mult binomial"){
	if(length(dim(y))!=2||ncol(y)!=2)
		stop(paste("Two column matrix required for response: successes and failures"))
	if(any(y<0))stop("All response values must be positive")
	n <- nrow(y)
	nn <- y[,1]+y[,2]
	censor <- F}
else {
	censor <- length(dim(y))==2&&ncol(y)==2
	if(censor&&all(y[,2]==1)){
		y <- y[,1]
		censor <- F}
	if(!censor){
		if(!is.vector(y,mode="numeric"))stop("y must be a vector")
		n <- length(y)
		if(distribution=="double Poisson"||distribution=="mult Poisson")
			my <- 3*max(y)}}
if((distribution=="inverse Gauss"||distribution=="exponential"||
	distribution=="gamma"||distribution=="Weibull"||
	distribution=="extreme value")&&((censor&&any(y[,1]<=0))||
	(!censor&&any(y<=0))))stop("All response values must be > 0")
if((distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="gamma count"||distribution=="double Poisson"||
	distribution=="mult Poisson")&&(any(y<0)))
	stop("All response values must be >= 0")
if(distribution=="logarithmic"&&any(y[wt>0]<1))
	stop("All response values must be integers > 0")
if(censor){
	n <- nrow(y)
	y[,2] <- as.integer(y[,2])
	if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
		stop("Censor indicator must be -1s, 0s, and 1s")
	cc <- ifelse(y[,2]==1,1,0)
	rc <- ifelse(y[,2]==0,1,ifelse(y[,2]==-1,-1,0))
	lc <- ifelse(y[,2]==-1,0,1)
	if(any(delta<=0&y[,2]==1))
		stop("All deltas for uncensored data must be positive")
	else {
		delta <- ifelse(delta<=0,0.000001,delta)
		delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001,delta)}}
else {
	if(min(delta)<=0)stop("All deltas for must be positive")}
if(length(wt)==1)wt <- rep(wt,n)
else if(length(wt)!=n)stop("wt must be the same length as the other variables")
if(min(wt)<0)stop("All weights must be non-negative")
if(length(delta)==1)delta <- rep(delta,n)
else if(length(delta)!=n)stop("delta must be the same length as the other variables")
if(is.null(nest))stop("A nest vector must be supplied")
else if(length(nest)!=n)stop("nest must be the same length as the other variables")
if(is.factor(nest))nest <- codes(nest)
nind <- length(unique(nest))
od <- length(nest)==nind
i <- rep(1:n,points)
ii <- rep(1:nind,points)
k <- NULL
for(j in 1:points)k <- c(k,nest+(j-1)*max(nest))
k <- as.integer(k)
quad <- gauss.hermite(points)
sd <- quad[rep(1:points,rep(n,points)),1]
qw <- quad[rep(1:points,rep(nind,points)),2]
lin1 <- lin2 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]}
else lin1 <- linear
if(inherits(mu,"formula"))lin1 <- mu
if(inherits(shape,"formula"))lin2 <- shape
lin1a <- lin2a <- mu3 <- sh3 <- name <- NULL
if(inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	type <- if(inherits(envir,"repeated"))"repeated"
		else if(inherits(envir,"tccov"))"tccov"
		else "tvcov"
	if(inherits(lin1,"formula")){
		lin1a <- finterp(lin1)
		class(lin1) <- c(class(lin1),type)}
	if(inherits(lin2,"formula")){
		lin2a <- finterp(lin2)
		class(lin2) <- c(class(lin2),type)}
	name <- paste(deparse(substitute(envir)))
	if(is.function(mu)){
		mu3 <- mu
		attributes(mu3) <- attributes(fnenvir(mu))
		class(mu) <- type
		mu <- fnenvir(mu,envir=envir,name=name)}
	if(is.function(shape)){
		sh3 <- shape
		attributes(sh3) <- attributes(fnenvir(shape))
		class(shape) <- type
		shape <- fnenvir(shape,envir=envir,name=name)}}
if(inherits(lin1,"formula")){
	mu2 <- finterp(lin1,envir=envir,name=name)
	npt1 <- length(attr(mu2,"parameters"))
	if(is.matrix(attr(mu2,"model"))){
		if(all(dim(attr(mu2,"model"))==1)){
			if(is.function(mu)){
				lin1 <- mu2
				mu1 <- function(p) mu(p,p[1]*rep(1,n))}
			else {
				mu1 <- function(p) p[1]*rep(1,n)
				attributes(mu1) <- attributes(mu2)}
			mu2 <- NULL}
		else {
			if(nrow(attr(mu2,"model"))!=n)stop("mu model matrix does not match number of response observations")
			if(is.function(mu)){
				dm1 <- attr(mu2,"model")
				lin1 <- mu2
				mu2 <- function(p) mu(p,dm1%*%p[1:npt1])}}}
	else {
		if(npl!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu2,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(!is.function(mu)){
	mu1 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	mu1 <- mu
	if(length(mu1(pmu))==1)mu1 <- function(p) mu(p)*rep(1,n)}
if(is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn"))attributes(fnenvir(mu))
		else attributes(mu)}
		else attributes(fnenvir(mu1))}
nlp <- if(is.function(mu)){
		if(is.null(lin1))length(attr(mu1,"parameters"))
		else length(attr(mu1,"parameters"))-1+npt1}
       else npt1
if(nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
npl <- npl+1
npl1 <- npl+1
if(inherits(lin2,"formula")){
	sh2 <- finterp(lin2,envir=envir,start=npl1,name=name)
	npt2 <- length(attr(sh2,"parameters"))
	if(is.matrix(attr(sh2,"model"))){
		if(all(dim(attr(sh2,"model"))==1)){
			if(is.function(shape)){
				lin2 <- sh2
				sh1 <- function(p) shape(p[npl1:np],p[npl1]*rep(1,n))}
			else {
				sh1 <- function(p) p[npl1]*rep(1,n)
				sh3 <- fnenvir(function(p) p[1]*rep(1,n))
				attributes(sh1) <- attributes(sh2)}
			sh2 <- NULL}
		else {
			if(nrow(attr(sh2,"model"))!=n)stop("shape model matrix does not match number of response observations")
			if(is.function(shape)){
				dm2 <- attr(sh2,"model")
				sh1 <- if(shfn)function(p) shape(p[npl1:np], dm2%*%p[npl1:(npl1+npt2-1)], mu1(p))
				else function(p) shape(p[npl1:np],dm2%*%p[npl1:(npl1+npt2-1)])}}}
	else {
		if(nps!=npt2){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(!is.function(shape)&&shp){
	sh1 <- function(p) p[npl1]*rep(1,n)
	sh3 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
else {
	sh1 <- if(shfn)function(p) shape(p[npl1:np], mu1(p))
		else function(p) shape(p[npl1:np])}
if(shp){
	if(is.null(attributes(sh1))){
		attributes(sh1) <- if(is.function(shape)){
			if(!inherits(shape,"formulafn"))attributes(fnenvir(shape))
			else attributes(shape)}
			else attributes(fnenvir(sh1))}
	nlp <- if(is.function(shape)){
			if(is.null(lin2))length(attr(sh1,"parameters"))-shfn
			else length(attr(sh1,"parameters"))-1+npt2-shfn}
		else npt2
	if(nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))}
pmu <- c(pmu,psd)
p <- c(pmu,pshape)
if(missing(scale)){
	if(distribution=="binomial"||distribution=="beta binomial"||
		distribution=="double binomial"||distribution=="mult binomial")
		scale <- "logit"
	else if(distribution=="normal"||distribution=="logistic"||
		distribution=="Cauchy"||distribution=="Laplace")
		scale <- "identity"
	else scale <- "log"}
mu2 <- if(scale=="logit") function(p){
		pp <- exp(log(mu1(p)/(1-mu1(p)))[i]+p[npl]*sd)
		pp/(1+pp)}
	else if(scale=="identity") function(p) mu1(p)[i]+p[npl]*sd
	else if(scale=="log") function(p) exp(log(mu1(p))[i]+p[npl]*sd)
	else if(scale=="reciprocal") function(p) 1/(1/mu1(p)[i]+p[npl]*sd)
	else if(scale=="exp") function(p) log(exp(mu1(p))[i]+p[npl]*sd)
if(any(is.na(mu2(pmu))))stop("The location model returns NAs")
if(distribution=="Levy"&&((!censor&&any(y<=mu1(pmu)))||(censor&&any(y[,1]<=mu1(pmu)))))
	stop("location parameter must be strictly less than corresponding observation")
if(distribution=="Pareto"&&exp(sh1(p))<=1)stop("shape parameters must be > 0")
if(distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"){
	if(any(is.na(sh1(p))))stop("The shape model returns NAs")
	if(od)stop("Some individuals must have more than one observation")}
if (!censor){
	ret <- switch(distribution,
	binomial={
		fcn <- function(p) {
			m <- mu2(p)
			-wt*(y[,1]*log(m)+y[,2]*log(1-m))}
		const <- -wt*lchoose(nn,y[,1])},
	"beta binomial"={
		fcn <- function(p) {
			m <- mu2(p)
			s <- exp(sh1(p))
			t <- s*m
			u <- s*(1-m)
			-wt*(lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u))}
		const <- -wt*lchoose(nn,y[,1])},
	"double binomial"={
		fcn <- function(p) {
			-.C("ddb",as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),
				DUP=F)$res}
		const <- 0},
	"mult binomial"={
		fcn <- function(p) {
			-.C("dmb",as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),
				DUP=F)$res}
		const <- 0},
	Poisson={
		fcn <- function(p) {
			m <- mu2(p)
			-wt*(-m+y*log(m))}
		const <- wt*lgamma(y+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu2(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(lgamma(y+s)-lgamma(s)+s*t+y*log(m)
				-(y+s)*log(s+m))}
		const <- wt*lgamma(y+1)},
	"double Poisson"={
		fcn <- function(p) {
			-.C("ddp",as.integer(y),as.integer(my),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),DUP=F)$res}
		const <- 0},
	"mult Poisson"={
		fcn <- function(p) {
			-.C("dmp",as.integer(y),as.integer(my),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),DUP=F)$res}
		const <- 0},
	"gamma count"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- exp(sh1(p))
			u <- m*s
			-wt*log(ifelse(y==0,1-pgamma(u,(y+1)*s,1),
				pgamma(u,y*s+(y==0),1)-
				pgamma(u,(y+1)*s,1)))}
		const <- 0},
	Consul={
		fcn <- function(p) {
			m <- mu2(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-y*t)}
		const <- wt*lgamma(y+1)},
	logarithmic={
		fcn <- function(p) {
			m <- exp(mu2(p))
			m <- m/(1+m)
			-wt*(y*log(m)-log(y)-log(-log(1-m)))}
		const <- 0},
	geometric={
		fcn <- function(p) {
			m <- mu2(p)
			-wt*(y*log(m)-(y+1)*log(1+m))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				-wt*log(pnorm(y+delta/2,m,s)
					-pnorm(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(t+(y-mu2(p))^2/exp(t))/2}
			const <- wt*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*log(pinvgauss(y+delta/2,m,s)-
					pinvgauss(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				wt*(t+(y-m)^2/(y*exp(t)*m^2))/2}
			const <- wt*(log(2*pi*y^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*log(plogis(y+delta/2,m,s)
					-plogis(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)*sqrt(3)/pi
				wt*((y-m)/s+t+2*log(1+exp(-(y-m)/s)))}
			const <- -wt*(log(pi/sqrt(3))+log(delta))}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				-wt*log(pcauchy(y+delta/2,m,s)
					-pcauchy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				wt*log(s*(1+(y-m)^2/s^2))}
			const <- -wt*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*log(plaplace((y+delta/2-m)/s)
					-plaplace((y-delta/2-m)/s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(abs(y-mu2(p))/exp(t)+t)}
			const <- -wt*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*log(plevy(y+delta/2,m,s)
					-plevy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(0.5*log(s/(2*pi))-1.5*log(y-m)-
					s/(2*(y-m)))}
			const <- -wt*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*log((1+(y-delta/2)*t)^-s
					-(1+(y+delta/2)*t)^-s)}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(log(s*t)-(s+1)*log(1+y*t))}
			const <- -wt*log(delta)}},
        exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				-wt*log(-exp(-(y+delta/2)/m)
					+exp(-(y-delta/2)/m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				wt*(log(m)+y/m)}
			const <- -wt*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*log(pgamma(y+delta/2,s,u)
					-pgamma(y-delta/2,s,u))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(s*(t-log(m)-y/m)+(s-1)*log(y)-lgamma(s))}
			const <- -wt*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*log(pweibull(y+delta/2,s,m)
					-pweibull(y-delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+(s-1)*log(y)-s*log(m)-(y/m)^s)}
			const <- -wt*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				-wt*log(pweibull(ey+ey*delta/2,s,m)
					-pweibull(ey-ey*delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+s*y-s*log(m)-(exp(y)/m)^s)}
			const <- -wt*log(delta)}},
	own={ const <- 0})}
else {
	ret <- switch(distribution,
	Poisson={
		fcn <- function(p) {
			m <- mu2(p)
			-wt*(cc*(-m+y[,1]*log(m))+
				log(lc-rc*ppois(y[,1],m)))}
		const <- wt*cc*lgamma(y[,1]+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu2(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(cc*(lgamma(y[,1]+s)-lgamma(s)
				+s*t+y[,1]*log(m)-(y[,1]+s)*log(s+m))+
				log(lc-rc*pnbinom(y[,1],s,1/(1+m/s))))}
		const <- wt*cc*lgamma(y[,1]+1)},
	geometric={
		fcn <- function(p) {
			m <- mu2(p)
			-wt*(cc*(y[,1]*log(m)-(y[,1]+1)*log(1+m))+
				log(lc-rc*pgeom(y[,1],1/(1+m))))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pnorm(y[,1]+delta/2,m,s)-
					pnorm(y[,1]-delta/2,m,s))
					+log(lc-rc*pnorm(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/s)/2)+log(lc-rc
					*pnorm(y[,1],m,sqrt(s))))}
			const <- wt*cc*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pinvgauss(y[,1]+delta/2,m,s)-
					pinvgauss(y[,1]-delta/2,m,s))
					+log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/(y[,1]*s*m^2))/2)+
					log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- wt*cc*(log(2*pi*y[,1]^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*(cc*log(plogis(y[,1]+delta/2,m,s)-
					plogis(y[,1]-delta/2,m,s))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				y1 <- (y[,1]-m)/s
				-wt*(cc*(-y1-log(s)-2*log(1+exp(-y1)))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- -wt*cc*log(delta)}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pcauchy(y[,1]+delta/2,m,s)-
					pcauchy(y[,1]-delta/2,m,s))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p)/2)
				-wt*(-cc*log(s*(1+(y[,1]-m)^2/s^2))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plaplace((y[,1]+delta/2-m)/s)-
					plaplace((y[,1]-delta/2-m)/s))
					+log(lc-rc*plaplace((y[,1]-m)/s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-abs(y[,1]-m)/s-t)+log(lc-rc
					*plaplace((y[,1]-m)/s)))}
			const <- -wt*cc*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plevy(y[,1]+delta/2,m,s)-
					plevy(y[,1]-delta/2,m,s))
					+log(lc-rc*plevy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(0.5*log(s/(2*pi))-1.5*log(y[,1]-m)-
					s/(2*(y[,1]-m)))+log(lc-rc
					*plevy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(cc*log((1+(y[,1]-delta/2)*t)^-s-
					(1+(y[,1]+delta/2)*t)^-s)
					+log(lc-rc*(-(1+(y[,1])*t)^-s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(cc*(log(s*t)-(s+1)*log(1+y[,1]*t))
					+log(lc-rc*(1-(1+y[,1]*t)^-s)))}
			const <- -wt*cc*log(delta)}},
	exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				-wt*(cc*log(-exp(-(y[,1]+delta/2)/m)
					+exp(-(y[,1]-delta/2)/m))+
					log(lc-rc*(1-exp(-y[,1]/m))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				-wt*(cc*(-log(m)-y[,1]/m)+log(lc-rc*
					(1-exp(-y[,1]/m))))}
			const <- -wt*cc*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*(cc*log(pgamma(y[,1]+delta/2,s,u)-
					pgamma(y[,1]-delta/2,s,u))
					+log(lc-rc*pgamma(y[,1],s,u)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(s*(t-log(m)-y[,1]/m)+(s-1)*log(y[,1])
					-lgamma(s))+log(lc-rc
					*pgamma(y[,1],s,m/s)))}
			const <- -wt*cc*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pweibull(y[,1]+delta/2,s,m)-
					pweibull(y[,1]-delta/2,s,m))
					+log(lc-rc*pweibull(y[,1],s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(t+(s-1)*log(y[,1])-s*log(m)
					-(y[,1]/m)^s)+log(lc-rc*
					pweibull(y[,1],s,m)))}
			const <- -wt*cc*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu2(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				-wt*(cc*log(pweibull(ey+ey*delta/2,s,m)-
					pweibull(ey-ey*delta/2,s,m))
					+log(lc-rc*pweibull(ey,s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu2(p)
				t <- sh1(p)
				s <- exp(t)
				ey <- exp(y[,1])
				-wt*(cc*(t+s*y[,1]-s*log(m)-(ey/m)^s)+log(lc-
					rc*pweibull(ey,s,m)))}
			const <- -wt*cc*log(delta)}},
	own={const <- 0})}
fn <- function(p) {
	under  <- 0
	if(od)pr <- -fcn(p)
	else {
		pr <- NULL
		for(i in split(fcn(p),k))pr <- c(pr,-sum(i))}
	if(any(is.na(pr)))stop("NAs - unable to calculate probabilities.\n Try other initial values.")
	if(max(pr)-min(pr)>1400){
		if(print.level==2)cat("Log probabilities:\n",pr,"\n\n")
		stop("Product of probabilities is too small to calculate.\n Try fewer points.")}
	if(any(pr > 700))under <- 700-max(pr)
	else if(any(pr < -700))under <- -700-min(pr)
	tmp <- NULL
	for(i in split(qw*exp(pr+under),ii))tmp <- c(tmp,sum(i))
	-sum(log(tmp)-under)}
if(fscale==1)fscale <- fn(p)
if(is.na(fn(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fn, p=p, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
if(inherits(mu1,"formulafn")&&!is.matrix(attr(mu1,"model")))
	cname <- attr(mu1,"parameters")
else if(inherits(lin1,"formula")){
	cname <- if(inherits(mu1,"formulafn"))colnames(attr(mu1,"model"))
		else colnames(dm1)
	if(is.null(cname)&&npl==2)cname <- "(Intercept)"
	if(is.function(mu)&&length(cname)<npl-1)
		cname <- c(cname,paste("p",(length(cname)+1):npl,sep=""))}
else cname <- paste("p",1:(npl-1),sep="")
if(inherits(sh1,"formulafn")&&!is.matrix(attr(sh1,"model")))
	sname <- attr(sh1,"parameters")
else if(inherits(lin2,"formula")){
	sname <- if(inherits(sh1,"formulafn"))colnames(attr(sh1,"model"))
		else colnames(dm2)
	if(is.null(sname)&&nps==1)sname <- "(Intercept)"
	if(is.function(shape)&&length(sname)<nps)
		sname <- c(sname,paste("p",(length(sname)+1):nps,sep=""))}
else sname <- paste("p",1:nps,sep="")
fitted.values <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial")
		as.vector((y[,1]+y[,2])*mu2(z0$estimate))
	else as.vector(mu2(z0$estimate))
residuals <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial"||censor)
		y[,1]-fitted.values
	else y-fitted.values
if(npl+nps==1){
	cov <- 1/z0$hessian
	se <- sqrt(cov)}
else {
	if(any(is.na(z0$hessian)))a <- 0
	else a <- qr(z0$hessian)$rank
	if(a==npl+nps)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=npl+nps,nrow=npl+nps)
	se <- sqrt(diag(cov))}
like.comp <- as.vector(fcn(z0$estimate)+const)
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	mu=mu1,
	shape=sh1,
	linear=list(lin1,lin2),
	scale=scale,
	points=points,
	prior.weights=wt,
	censor=censor,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	like.comp=like.comp,
	aic=z0$minimum+np,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	cname=cname,
	sname=sname,
	npl=npl,
	npm=0,
	nps=nps,
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlr"
return(z1)}
