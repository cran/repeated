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
#     gar(response, distribution="normal", times=NULL, totals=NULL,
#	censor=NULL, delta=NULL, mu=NULL, shape=NULL, shfn=F,
#	common=F, preg=NULL, pshape=NULL, pdepend=NULL,
#	transform="identity", link="identity", autocorr="exponential",
#	order=1, envir=sys.frame(sys.parent()), print.level=0, ndigit=10,
#	gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
#	typsiz=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit generalized nonlinear autoregression models with
#  various distributions

gar <- function(response, distribution="normal", times=NULL, totals=NULL,
	censor=NULL, delta=NULL, mu=NULL, shape=NULL, shfn=F, common=F,
	preg=NULL, pshape=NULL, pdepend=NULL, transform="identity",
	link="identity", autocorr="exponential", order=1,
	envir=sys.frame(sys.parent()),print.level=0, ndigit=10,
	gradtol=0.00001, steptol=0.00001,fscale=1, iterlim=100,
	typsiz=abs(p), stepmax=10*sqrt(p%*%p)){
likekal <- function(p){
	eta <- mu1(p)
	if(sh){
		shr <- sh1(p)
		theta <- c(p[npr1:nprd],exp(p[np]))}
	else theta <- p[npr1:np]
	z <- .C("gar",
		y=y,
		total=response$response$n,
		my=as.integer(3*max(y)),
		nobs=as.integer(response$response$nobs),
		nind=as.integer(nind),
		times=as.double(response$response$times),
		censor=as.integer(censor),
		cens=as.integer(!is.null(response$response$censor)),
		eta=eta,
		theta=theta,
		model=as.integer(mdl),
		thp=as.integer(thp),
		shape=shr,
		sh=as.integer(sh),
		link=as.integer(lnk),
		ar=as.integer(ar),
		order=as.integer(order),
		pred=double(n),
		rpred=double(n),
		like=double(1),
		DUP=F)
	z$like+jacob}
call <- sys.call()
tmp <- c("binomial","Poisson","exponential","negative binomial",
	"mult Poisson","double Poisson","beta binomial","mult binomial",
	"double binomial","normal","logistic","Cauchy", "Weibull","gamma",
	"Laplace","inverse Gauss","Pareto","Levy","gen gamma",
	"gen logistic","Hjorth","Burr","gen Weibull","gen extreme value",
	"gen inverse Gauss","power exponential")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
tmp <- c("exponential","gaussian","cauchy","spherical","IOU")
ar <- match(autocorr <- match.arg(autocorr,tmp),tmp)
transform <- match.arg(transform,c("identity","exp","square","sqrt","log"))
tmp <- c("identity","exp","square","sqrt","log","logit","cloglog")
lnk <- match(link <- match.arg(link,tmp),tmp)
if((link=="logit"||link=="cloglog")&&(mdl!=1&&mdl!=7&&mdl!=8&&mdl!=9))stop("logit and cloglog links can only be used with binary data")
if(is.null(times)&&is.matrix(response))times <- matrix(1,nrow=nrow(response),ncol=ncol(response))
respenv <- inherits(response,"repeated")
envname <- if(respenv)paste(deparse(substitute(response)))
	else NULL
if(!respenv){
	if(!inherits(response,"response")){
		response <- if(mdl==1||mdl==7||mdl==8||mdl==9)restovec(response, times=times, censor=censor, delta=delta, totals=totals)
		else restovec(response, times=times, censor=censor, delta=delta)}
	response <- rmna(response=response)}
if((inherits(envir,"repeated")&&
	(length(response$response$nobs)!=length(envir$response$nobs)||
	any(response$response$nobs!=envir$response$nobs)))||
	(inherits(envir,"tvcov")&&
	(length(response$response$nobs)!=length(envir$tvcov$nobs)||
	any(response$response$nobs!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
if(mdl==1||mdl==7||mdl==8||mdl==9){
	if(is.null(response$response$n)){
		if(any(response$response$y!=0&&response$response$y!=1))stop("responses must be binary if totals are not supplied")
		else response$response$n <- rep(1,length(response$response$y))}}
y <- response$response$y
n <- length(y)
mu3 <- sh3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(respenv||inherits(envir,"repeated"))"repeated"
		else if(inherits(envir,"tccov"))"tccov"
		else "tvcov"
	if(is.null(envname))envname <- paste(deparse(substitute(envir)))
	if(inherits(mu,"formula")){
		mu3 <- if(respenv)finterp(mu,envir=response,name=envname)
			else finterp(mu,envir=envir,name=envname)
		class(mu) <- c(class(mu),type)}
	else if(is.function(mu)){
		tmp <- parse(text=paste(deparse(mu))[-1])
		class(mu) <- type
		mu <- if(respenv)fnenvir(mu,envir=response,name=envname)
			else fnenvir(mu,envir=envir,name=envname)
		mu3 <- mu
		if(respenv)attr(mu3,"model") <- tmp}
	if(inherits(shape,"formula")){
		sh3 <- if(respenv)finterp(shape,envir=response,name=envname)
			else finterp(shape,envir=envir,name=envname)
		class(shape) <- c(class(shape),type)}
	else if(is.function(shape)){
		tmp <- parse(text=paste(deparse(shape))[-1])
		class(shape) <- type
		shape <- if(respenv)fnenvir(shape,envir=response,name=envname)
			else fnenvir(shape,envir=envir,name=envname)
		sh3 <- shape
		if(respenv)attr(sh3,"model") <- tmp}}
npr <- length(preg)
npr1 <- npr+1
nprd <- npr+length(pdepend)
nprd1 <- if(common) 1 else nprd+1
if(inherits(mu,"formula")){
	mu2 <- if(respenv)finterp(mu,envir=response,name=envname)
		else finterp(mu,envir=envir,name=envname)
	npt1 <- length(attr(mu2,"parameters"))
	if(is.matrix(attr(mu2,"model"))){
		if(all(dim(attr(mu2,"model"))==1)){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
		if(npr!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg should have",npt1,"estimates"))}
		if(is.list(preg)){
			if(!is.null(names(preg))){
				o <- match(attr(mu2,"parameters"),names(preg))
				preg <- unlist(preg)[o]
				if(sum(!is.na(o))!=length(preg))stop("invalid estimates for mu - probably wrong names")}
			else preg <- unlist(preg)}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.null(mu)){
	mu1 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else mu1 <- mu
if(!is.null(mu1)&&is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,envir=response))
			else attributes(fnenvir(mu,envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,envir=response))
			else attributes(fnenvir(mu1,envir=envir))}}
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)&&!common&&nlp!=npr)
	stop(paste("preg should have",nlp,"initial estimates"))
nps <- length(pshape)
if(inherits(shape,"formula")){
	sh2 <- if(respenv)finterp(shape,envir=response,start=nprd1,name=envname)
		else finterp(shape,envir=envir,start=nprd1,name=envname)
	npt2 <- length(attr(sh2,"parameters"))
	if(is.matrix(attr(sh2,"model"))){
		if(all(dim(attr(sh2,"model"))==1)){
			sh1 <- function(p) p[nprd1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
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
else if(!is.function(shape)&&distribution!="binomial"&&
	distribution!="Poisson"&&distribution!="exponential"&&
	distribution!="geometric"&&distribution!="logarithmic"){
	sh1 <- function(p) p[nprd1]*rep(1,n)
	if(length(pshape)!=1)stop("pshape must provide an estimate")}
else sh1 <- if(shfn)function(p) shape(p[nprd1:np], mu1(p))
	else function(p) shape(p[nprd1:np])
if(!is.null(sh1)&&is.null(attributes(sh1)))
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn")){
			if(respenv)attributes(fnenvir(shape,envir=response))
			else attributes(fnenvir(shape,envir=envir))}
		else attributes(shape)}
		else {
			if(respenv)attributes(fnenvir(sh1,envir=response))
			else attributes(fnenvir(sh1,envir=envir))}
nlp <- if(is.function(shape))length(attr(sh1,"parameters"))-shfn
	else if(is.null(shape))NULL
	else npt2
if(!is.null(nlp)&&!common&&nlp!=nps)
	stop(paste("pshape should have",nlp,"initial estimates"))
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"))))
	if(nlp!=npr)stop(paste("with a common parameter model, preg should contain",nlp,"estimates"))}
if(mdl<10){
	if(any(y<0))stop("all responses must be non-negative")}
else if((mdl!=10)&&(mdl!=11)&&(mdl!=12)&&(mdl!=15)&&(mdl!=18)&&(mdl!=20)&&(mdl!=26)&&(any(y<=0)))
	stop("all responses must be positive")
else if(distribution=="Levy"&&any(response$response$y<=mu1(preg)))
	stop("location function must be strictly less than corresponding observation")
if(distribution=="Pareto"&&pshape<=1)stop("shape parameter must be > 1")
censor <- response$response$censor
if(is.null(censor))censor <- rep(1,n)
else if(mdl==1||mdl==2||mdl==4||mdl==5||mdl==6||mdl==7||mdl==8||mdl==9)stop(paste("Censored data not allowed for",distribution,"distribution"))
nind <- length(response$response$nobs)
if(transform=="identity")jacob <- 0
else if(transform=="exp"){
	jacob <- -sum(y[censor==1])
	y <- exp(y)}
else {
	if(any(y<0))stop("Nonpositive response values: invalid transformation")
	else if(transform=="square"){
		jacob <- -sum(log(y[y>0&censor==1]))
		y  <- y^2}
	else if(transform=="sqrt"){
		jacob <- sum(log(y[y>0&censor==1]))/2
		y <- sqrt(y)}
	else if(any(y==0))stop("Zero response values: invalid transformation")
	else if(transform=="log"){
		jacob <- sum(log(y[censor==1]))
		y <- log(y)}}
if(!is.null(response$response$delta)){
	if(length(response$response$delta)==1)
		jacob <- jacob-length(y[censor==1])*log(response$response$delta)
	else jacob <- jacob-sum(log(response$response$delta[censor==1]))}
if(order!=1&&order!=2)stop("Autoregression must have order 1 or 2")
if(is.null(pdepend))
	stop("Initial estimates of the dependence parameters must be supplied")
if(order==2&&length(pdepend)!=2)stop("2 estimates of dependence parameters must be given")
else if(length(pdepend)!=1&&length(pdepend)!=2)
     stop("One or two estimates of dependence parameters must be given")
thp <- length(pdepend)==2&&order==1
sh <- is.function(shape)||inherits(shape,"formula")
if(length(mu1(preg))!=length(response$response$y))
	stop("The mu function must provide an estimate for each observation")
else if(any(is.na(mu1(preg))))
	stop("Non-numerical mu: probably invalid initial values")
if(any(pdepend<=0))stop("All dependence parameters must be positive")
p <- c(preg,-log(pdepend))
if(mdl>3){
	if(!sh)p <- c(p,log(pshape))
	else {
	     if(mdl>=19)p <- if(common)c(p,log(pshape[nps]))
			else c(p,pshape[1:(nps-1)],log(pshape[nps]))
	     else p <- c(p,pshape)}}
np <- length(p)
if(!sh){
	if((mdl<=3&&nps!=0)||(mdl>3&&mdl<19&&nps!=1)||(mdl>=19&&nps!=2))
		stop("Incorrect number of shape parameter estimates")
	else if(nps>0&&any(pshape<=0))
		stop("All shape parameters must be positive")
	shr <- rep(0,length(response$response$y))}
else {
	if(any(is.na(sh1(p))))stop("The shape model returns NAs: probably invalid initial values")
	if(length(sh1(p))!=length(response$response$y))stop("The shape function must provide an estimate for each observation")}
if(fscale==1)fscale <- likekal(p)
if(is.na(likekal(p)))stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(likekal, p, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
like <- z0$minimum
if(any(is.na(z0$hessian)))a <- 0
else a <- qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
eta <- mu1(z0$estimate)
if(sh){
	shr <- sh1(z0$estimate)
	theta <- c(z0$estimate[npr1:nprd],exp(z0$estimate[np]))}
else theta <- z0$estimate[npr1:np]
z <- .C("gar",
	y=y,
	total=response$response$n,
	my=as.integer(3*max(y)),
	nobs=as.integer(response$response$nobs),
	nind=as.integer(nind),
	times=as.double(response$response$times),
	censor=as.integer(censor),
	cens=as.integer(!is.null(response$response$censor)),
	eta=eta,
	theta=theta,
	model=as.integer(mdl),
	thp=as.integer(thp),
	shape=shr,
	sh=as.integer(sh),
	link=as.integer(lnk),
	ar=as.integer(ar),
	order=as.integer(order),
	pred=double(n),
	rpred=double(n),
	like=double(1),
	DUP=F)
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
z <- list(
	call=call,
	distribution=distribution,
	mu=mu1,
	formula=mu,
	shape=shape,
	sh1=sh1,
	shfn=shfn,
	common=common,
	response=response$response,
	link=link,
	order=order,
	autocorr=autocorr,
	transform=transform,
	maxlike=like,
	aic=like+np,
	df=length(response$response$y)-np,
	np=np,
	npr=npr,
	nps=nps,
	thp=thp,
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=corr,
	pred=z$pred,
	rpred=z$rpred,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- c("gar","recursive")
return(z)}

coefficients.gar <- function(z) z$coefficients
deviance.gar <- function(z) 2*z$maxlike
fitted.gar <- function(z, recursive=TRUE)
	if(recursive) z$rpred else z$pred
residuals.gar <- function(z, recursive=TRUE){
	if(z$transform=="exp")z$response$y <- exp(z$response$y)
	else if(z$transform=="square")z$response$y  <- z$response$y^2
	else if(z$transform=="sqrt")z$response$y <- sqrt(z$response$y)
	else if(z$transform=="log")z$response$y <- log(z$response$y)
	if(recursive) z$response$y-z$rpred else z$response$y-z$pred}

print.gar <- function(z, digits = max(3, .Options$digits - 3)) {
	np1 <- if(z$distribution=="binomial"||z$distribution=="exponential"
			||z$distribution=="Poisson") 0
		else if(z$distribution=="gen gamma"
			||z$distribution=="gen logistic"
			||z$distribution=="Hjorth"||z$distribution=="Burr"
			||z$distribution=="gen Weibull"
			||z$distribution=="gen extreme value"
			||z$distribution=="gen inverse Gauss"
			||z$distribution=="power exponential") 2
		else 1
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$response$nobs),"\n")
	cat("Number of observations",length(z$response$y),"\n")
	cat("Transformation        ",z$trans,"\n")
	cat("Link function         ",z$link,"\n\n")
	cat(z$distribution,"distribution\n")
	if(z$order==1)cat("First order ")
	else cat("Second order ")
	cat(z$autocorr,"dependence\n")
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	if(z$common)cat("Location model\n")
	else cat("Location parameters\n")
	if(!is.null(attr(z$mu,"formula")))
		cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.matrix(attr(z$mu,"model")))
		colnames(attr(z$mu,"model"))
		else attr(z$mu,"parameters")
	coef.table <- cbind(z$coef[1:z$npr],z$se[1:z$npr])
	colname <- c("estimate","se")
	if(!z$common){
		dimnames(coef.table) <- list(cname,colname)
		print.default(coef.table, digits=digits, print.gap=2)}
	else {
		cat("\nShape model\n")
	        if(!is.null(attr(z$sh1,"formula")))
	        	cat(deparse(attr(z$sh1,"formula")),sep="\n")
	        else if(!is.null(attr(z$sh1,"model"))){
	        	t <- deparse(attr(z$sh1,"model"))
	        	t[1] <- sub("expression\\(","",t[1])
	        	t[length(t)] <- sub("\\)$","",t[length(t)])
	        	cat(t,sep="\n")}
	        cname <- c(cname,if(is.matrix(attr(z$sh1,"model")))
				colnames(attr(z$sh1,"model"))
			else attr(z$sh1,"parameters")[1:(length(attr(z$sh1,"parameters"))-z$shfn)])
		cname <- unique(cname)
		cat("\nCommon parameters\n")
		dimnames(coef.table) <- list(cname,colname)
		print.default(coef.table, digits=digits, print.gap=2)}
	if(z$thp||z$order==2){
		cat("\nDependence parameters\n")
		if(z$thp)cname <- c("phi","rho")
		else cname <- c("rho1","rho2")
		coef.table <- cbind(z$coef[(z$npr+1):(z$npr+2)],
			z$se[(z$npr+1):(z$npr+2)],
			exp(-z$coef[(z$npr+1):(z$npr+2)]))}
	else {
		cat("\nDependence parameter\n")
		cname <- "rho"
		coef.table <- cbind(z$coef[z$npr+1],
			z$se[z$npr+1],exp(-z$coef[z$npr+1]))}
	dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
	print.default(coef.table, digits=digits, print.gap=2)
	if(np1>0&&!z$common){
		cat("\nShape parameters\n")
		if(is.null(z$shape)){
			cname <- "shape"
			if(np1==2)cname <- c(cname,"psi")
			coef.table <- cbind(z$coef[(z$np-np1+1):z$np],
				z$se[(z$np-np1+1):z$np],
				exp(z$coef[(z$np-np1+1):z$np]))
				dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))}
		else {
	                if(!is.null(attr(z$sh1,"formula")))
	                	cat(deparse(attr(z$sh1,"formula")),sep="\n")
	                else if(!is.null(attr(z$sh1,"model"))){
	                	t <- deparse(attr(z$sh1,"model"))
	                	t[1] <- sub("expression\\(","",t[1])
	                	t[length(t)] <- sub("\\)$","",t[length(t)])
	                	cat(t,sep="\n")}
	                cname <- if(is.matrix(attr(z$sh1,"model")))
					colnames(attr(z$sh1,"model"))
				else attr(z$sh1,"parameters")[1:(length(attr(z$sh1,"parameters"))-z$shfn)]
			np2 <- length(cname)
			coef.table <- cbind(z$coef[(z$np-np2+1):z$np],
				z$se[(z$np-np2+1):z$np])
			dimnames(coef.table) <- list(cname, c("estimate","se"))
	                if(np1==2){
	                	cname[length(cname)] <- "psi"
	                	coef.table <- cbind(coef.table,c(rep(NA,nrow(coef.table)-1),exp(z$coef[z$np])))
	                	colname <- c(colname,"parameter")}}
		print.default(coef.table,digits=digits,print.gap=2)}
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}
