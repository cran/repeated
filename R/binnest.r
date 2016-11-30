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


##' Binary Random Effects Models with Two Levels of Nesting
##' 
##' \code{binnest} is designed to handle binary and binomial data with two
##' levels of nesting. The first level is the individual and the second will
##' consist of clusters within individuals.
##' 
##' The variance components at the two levels can only depend on the covariates
##' if \code{response} has class, \code{repeated}.
##' 
##' 
##' @param response A list of three column matrices with counts, corresponding
##' totals (not necessary if the response is binary), and (second-level)
##' nesting indicator for each individual, one matrix or dataframe of such
##' counts, or an object of class, response (created by
##' \code{\link[rmutil]{restovec}}) or repeated (created by
##' \code{\link[rmutil]{rmna}}).
##' @param totals If \code{response} is a matrix or dataframe, a corresponding
##' matrix or dataframe of totals (not necessary if the response is binary).
##' Ignored otherwise.
##' @param nest If \code{response} is a matrix or dataframe, a corresponding
##' matrix or dataframe of nesting indicators. Ignored otherwise.
##' @param ccov If \code{response} is a matrix, dataframe, list, or object of
##' class, \code{response}, a matrix of time-constant covariates or an object
##' of class, \code{tccov} (created by \code{\link[rmutil]{tcctomat}}). All of
##' these covariates are used in the fixed effects part of the model. Ignored
##' if response has class, \code{repeated}.
##' @param tvcov If \code{response} is a matrix, dataframe, list, or object of
##' class, \code{response}, an object of class, \code{tvcov} (created by
##' \code{\link[rmutil]{tvctomat}}). All of these covariates are used in the
##' fixed effects part of the model. Ignored if response has class,
##' \code{repeated}.
##' @param mu If \code{response} has class, \code{repeated}, a formula
##' beginning with ~, specifying a linear regression function for the fixed
##' effects, in the Wilkinson and Rogers notation, containing selected
##' covariates in the response object. (A logit link is assumed.)
##' @param re1 If \code{response} has class, \code{repeated}, a formula
##' beginning with ~, specifying a linear regression function for the variance
##' of the first level of nesting, in the Wilkinson and Rogers notation,
##' containing selected covariates in the response object. If NULL, a random
##' effect is not fitted at this level. (A log link is assumed.)
##' @param re2 If \code{response} has class, \code{repeated}, a formula
##' beginning with ~, specifying a linear regression function for the variance
##' of the second level of nesting, in the Wilkinson and Rogers notation,
##' containing selected covariates in the response object. If NULL, a random
##' effect is not fitted at this level. (A log link is assumed.)
##' @param preg Initial parameter estimates for the fixed effect regression
##' model: either the model specified by \code{mu} or else the intercept plus
##' one for each covariate in \code{ccov} and \code{tvcov}.
##' @param pre1 Initial parameter estimates for the first level of nesting
##' variance model: either the model specified by \code{re1} or just the
##' intercept. If NULL, a random effect is not fitted at this level.
##' @param pre2 Initial parameter estimates for the second level of nesting
##' variance model: either the model specified by \code{re1} or just the
##' intercept. If NULL, a random effect is not fitted at this level.
##' @param binom.mix A vector of two values giving the totals for the binomial
##' distributions used as the mixing distributions at the two levels of
##' nesting.
##' @param binom.prob A vector of two values giving the probabilities in the
##' binomial distributions used as the mixing distributions at the two levels
##' of nesting. If they are 0.5, the mixing distributions approximate normal
##' mixing distributions; otherwise, they are skewed.
##' @param fcalls Number of function calls allowed.
##' @param eps Convergence criterion.
##' @param print.level If 1, the iterations are printed out.
##' @return A list of classes \code{binnest} is returned.
##' @author T.R. Ten Have and J.K. Lindsey
##' @seealso \code{\link[repeated]{gar}}, \code{\link[rmutil]{read.list}},
##' \code{\link[rmutil]{restovec}}, \code{\link[rmutil]{rmna}},
##' \code{\link[rmutil]{tcctomat}}, \code{\link[rmutil]{tvctomat}}.
##' @references Ten Have, T.R., Kunselman, A.R., and Tran, L. (1999) Statistics
##' in Medicine 18, 947-960.
##' @keywords models
##' @examples
##' 
##' #y <- rbind(matrix(rbinom(20,1,0.6), ncol=4),
##' #	matrix(rbinom(20,1,0.4), ncol=4))
##' y <- matrix(c(1,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,0,1,0,
##' 	1,1,1,1,1,1,1,1,0,1,1,0),nrow=10,ncol=4,byrow=TRUE)
##' resp <- restovec(y, nest=1:4, times=FALSE)
##' ccov <- tcctomat(c(rep(0,5),rep(1,5)), name="treatment")
##' reps <- rmna(resp, ccov=ccov)
##' 
##' # two random effects
##' binnest(reps, mu=~treatment, preg=c(1,0), pre1=1, pre2=1)
##' 
##' # first level random effect only
##' binnest(reps, mu=~treatment, preg=c(1,-1), pre1=1)
##' 
##' # second level random effect only
##' binnest(reps, mu=~treatment, preg=c(1,-1), pre2=1)
##' 
##' @export binnest
binnest <- function(response, totals=NULL, nest=NULL, ccov=NULL, tvcov=NULL,
	mu=~1, re1=~1, re2=~1, preg=NULL, pre1=NULL, pre2=NULL,
	binom.mix=c(10,10), binom.prob=c(0.5,0.5), fcalls=900,
	eps=0.01, print.level=0){
#
# Fortran constants
#
maxt1 <- maxt2 <- maxt3 <- 10
call <- sys.call()
#
# check initial estimates
if(length(binom.prob)==1)binom.prob <- c(binom.prob,binom.prob)
if(any(binom.prob<=0)||any(binom.prob>=1))
	stop("binom.prob parameters must be between zero and one")
if(length(binom.mix)==1)binom.mix <- c(binom.mix,binom.mix)
#
# calculate constants
#
total1 <- length(preg)
total2 <- length(pre1)
t2 <- if(total2>0)total2 else 1
total3 <- length(pre2)
t3 <- if(total3>0)total3 else 1
total <- total1+total2+total3
dimw <- total*(total+7)/2
#
# set up response and regression functions
#
type <- "unknown"
if(inherits(response,"repeated")){
	if(dim(response$response$y)[2]>1)
		stop("binnest only handles univariate responses")
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("binnest does not handle data with NAs")
	type <- response$response$type
	envname <- paste(deparse(substitute(response)))
	# prepare response
	resp <- response$response$y
	if(is.null(response$response$n)){
		if(any(response$response$y!=0&response$response$y!=1))
			stop("if binomial totals are not supplied, all responses must be 0 or 1")
		else resp <- cbind(resp,rep(1,dim(response$response$y)[1]))}
	else resp <- cbind(resp,response$response$n)
	# location function or formula
	if(inherits(mu,"formula")){
		if(as.character(mu)[2]!="1"&&is.character(attr(finterp(mu,.envir=response,.name=envname),"model"))){
			tmp <- wr(mu,data=response)$design
			if(dim(tmp)[2]!=total1)
				stop(paste("preg should contain",dim(tmp)[2],"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=FALSE])
			regname <- gsub("\\[.i]","",colnames(tmp))}
		else {
			regname <- "(Intercept)"
			if(total1!=1)
				stop("preg should contain 1 initial estimate")}}
	else stop("mu must be a W&R formula")
	# function or formula for first level of random effect
	if(inherits(re1,"formula")&&!is.null(pre1)){
		if(as.character(re1)[2]!="1"&&is.character(attr(finterp(re1,.envir=response,.name=envname),"model"))){
			tmp <- wr(re1,data=response)$design
			if(dim(tmp)[2]!=total2)
				stop(paste("pre1 should contain",dim(tmp)[2],"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=FALSE])
			re1name <- gsub("\\[.i]","",colnames(tmp))}
		else {
			re1name <- "(Intercept)"
			if(total2!=1)
				stop("pre1 should contain 1 initial estimate")}}
	else pre1 <- re1name <- NULL
	# function or formula for second level of random effect
	if(inherits(re2,"formula")&&!is.null(pre2)){
		if(as.character(re2)[2]!="1"&&is.character(attr(finterp(re2,.envir=response,.name=envname),"model"))){
			tmp <- wr(re2,data=response)$design
			if(dim(tmp)[2]!=total3)
				stop(paste("pre2 should contain",dim(tmp)[2],"initial estimates"))
			resp <- cbind(resp,tmp[,-1,drop=FALSE])
			re2name <- gsub("\\[.i]","",colnames(tmp))}
		else {
			re2name <- "(Intercept)"
			if(total3!=1)
				stop("pre2 should contain 1 initial estimate")}}
	else pre2 <- re2name <- NULL
	nind <- length(nobs(response))
	numsubj <- dim(response$response$y)[1]
	nest <- response$response$nest}
else {
	# not a data object, so make one for the response
	if(inherits(response,"response")){
		if(dim(response$y)[2]>1)
			stop("binnest only handles univariate responses")}
	else {
		if(is.matrix(response)||is.data.frame(response))
			response <- restovec(response,totals=totals,nest=nest,times=FALSE)
		else if(is.list(response))
			response <- restovec(response,times=FALSE)
		else stop("response must be a matrix, data.frame, list, or object of type repeated or response")}
	resp <- response$y
	type <- response$type
	if(is.null(response$n)){
		if(any(response$y!=0&response$y!=1))
			stop("if binomial totals are not supplied, all responses must be 0 or 1")
		else resp <- cbind(resp,rep(1,dim(response$y)[1]))}
	else resp <- cbind(resp,response$n)
	# prepare covariates by making data objects
	if(!is.null(ccov)){
		if(!inherits(ccov,"tccov"))ccov <- tcctomat(ccov)
		resp <- cbind(resp,ccov$ccov[covind(response),,drop=FALSE])}
	if(!is.null(tvcov)){
		if(!inherits(tvcov,"tvcov"))tvcov <- tvctomat(tvcov)
		resp <- cbind(resp,tvcov$tvcov)}
	nind <- length(nobs(response))
	numsubj <- dim(response$y)[1]
	nest <- response$nest
	regname <- if(dim(resp)[2]>2)
		c("(Intercept)",colnames(resp[,3:dim(resp)])[2])
		else "(Intercept)"
	if(total2==1)re1name <- "(Intercept)"
	else if(!is.null(total2))stop("pre1 should supply 1 initial estimate")
	if(total3==1)re2name <- "(Intercept)"
	else if(!is.null(total3))stop("pre2 should supply 1 initial estimate")
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
if(type!="unknown"&&type!="nominal")stop("nominal data required")
#
# calculate numbers of observations at two levels
#
nobs1 <- nobs2 <- NULL
if(!is.null(nest))for(i in 1:nind){
	nobs1 <- c(nobs1,length(unique(nest[covind(response)==i])))
	nobs2 <- c(nobs2,as.vector(table(nest[covind(response)==i])))}
else nobs1 <- nobs2 <- rep(1,nind)
rm(response)
maxmother <- max(nobs1)
maxkid <- max(nobs2)
p <- c(preg,pre1,pre2)
#
# call Fortran to optimize likelihood
#
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
	totcol=as.integer(dim(resp)[2]),
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
	hac=double(maxmother*t2),
	had=double(maxmother*t3),
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
	uu1=double(nind*maxmother*maxkid*t2),
	uu2=double(nind*maxmother*maxkid*t3),
	w=double(dimw),
	habb=double(maxmother*total1*total1),
	habs1=double(maxmother*total1*t2),
	habs2=double(maxmother*total1*t3),
	has1s1=double(maxmother*t2*t2),
	has1s2=double(maxmother*t2*t3),
	has2s2=double(maxmother*t3*t3),
	ebb=double(maxmother*total1*total1),
	ebs1=double(maxmother*total1*t2),
	ebs2=double(maxmother*total1*t3),
	fs1s1=double(maxmother*t2*t2),
	fs1s2=double(maxmother*t2*t3),
	gs2s2=double(maxmother*t3*t3),
	e2bb=double(maxmother*maxmother*total1*total1),
	e2bs1=double(maxmother*maxmother*total1*t2),
	e2bs2=double(maxmother*maxmother*total1*t3),
	f2s1s1=double(maxmother*maxmother*t2*t2),
	f2s1s2=double(maxmother*maxmother*t2*t3),
	g2s2s2=double(maxmother*maxmother*t3*t3),
	DUP=TRUE,
	PACKAGE="repeated")
#
# warnings if convergence problems
#
if(z0$res[3]>0)switch(as.character(z0$res[3]),
		"1"=warning("Maximum number of function evaluations has been used"),
		"2"=stop("Linear search failed to improve the function value. Either the function or the gradient is incorrectly coded"),
		"3"=stop("Search vector was not a descent direction. The convergence criterion may be too strict"))
#
# calculate se's
#
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
        df=dim(resp)[1]-total,
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

### print method
###
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
	print.default(z$corr, digits=4)}
