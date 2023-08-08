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
#     hidden(response=NULL, totals=NULL, distribution="Bernoulli",
#	mu=NULL, cmu=NULL, tvmu=NULL, pgamma, pmu=NULL, pcmu=NULL,
#       ptvmu=NULL, pshape=NULL, pfamily=NULL, par=NULL, delta=NULL,
#       envir=parent.frame(), print.level=0,
#	ndigit=10, gradtol=0.00001, steptol=0.00001,
#       fscale=1, iterlim=100, typsize=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit hidden Markov chain models in discrete time
# Authors: J.K. Lindsey and P.J. Lindsey for the ordinal models



##' Discrete-time Hidden Markov Chain Models
##' 
##' \code{hidden} fits a two or more state hidden Markov chain model with a
##' variety of distributions. All series on different individuals are assumed
##' to start at the same time point. Time points are equal, discrete steps.
##' 
##' To fit an `observed' Markov chain, as well, with Bernoulli or multinomial
##' responses, use the lagged response as a time-varying covariate. For
##' quantitative responses, specifying \code{par} allows an `observed"
##' autoregression to be fitted as well as the hidden Markov chain.
##' 
##' All functions and formulae for the location parameter are on the
##' (generalized) logit scale for the Bernoulli, binomial, and multinomial
##' distributions.
##' 
##' If \code{cmu} and \code{tvmu} are used, these two mean functions are
##' additive so that interactions between time-constant and time-varying
##' variables are not possible.
##' 
##' The object returned can be plotted to give the probabilities of being in
##' each hidden state at each time point. For distributions other than the
##' multinomial, proportional odds, and continuation ratio, the (recursive)
##' predicted values can be plotted using \code{\link[rmutil]{mprofile}} and
##' \code{\link[rmutil]{iprofile}}.
##' 
##' See MacDonald, I.L. and Zucchini, W. (1997) Hidden Markov and Other Models
##' for Discrete-valued Time Series. Chapman and Hall.
##' 
##' 
##' @aliases hidden plot.hidden
##' @param response A list of two or three column matrices with counts or
##' category indicators, times, and possibly totals (if the distribution is
##' binomial), for each individual, one matrix or dataframe of counts, or an
##' object of class, \code{response} (created by
##' \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
##' \code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
##' \code{repeated} data object contains more than one response variable, give
##' that object in \code{envir} and give the name of the response variable to
##' be used here. If there is only one series, a vector of responses may be
##' supplied instead.
##' 
##' Multinomial and ordinal categories must be integers numbered from 0.
##' @param totals If response is a matrix, a corresponding matrix of totals if
##' the distribution is binomial. Ignored if response has class,
##' \code{response} or \code{repeated}.
##' @param distribution Bernoulli, Poisson, multinomial, proportional odds,
##' continuation ratio, binomial, exponential, beta binomial, negative
##' binomial, normal, inverse Gauss, logistic, gamma, Weibull, Cauchy, Laplace,
##' Levy, Pareto, gen(eralized) gamma, gen(eralized) logistic, Hjorth, Burr,
##' gen(eralized) Weibull, gen(eralized) extreme value, gen(eralized) inverse
##' Gauss, power exponential, skew Laplace, or Student t. (For definitions of
##' distributions, see the corresponding [dpqr]distribution help.)
##' @param mu A general location function with two possibilities: (1) a list of
##' formulae (with parameters having different names) or functions (with one
##' parameter vector numbering for all of them) each returning one value per
##' observation; or (2) a single formula or function which will be used for all
##' states (and all categories if multinomial) but with different parameter
##' values in each so that pmu must be a vector of length the number of
##' unknowns in the function or formula times the number of states (times the
##' number of categories minus one if multinomial).
##' @param cmu A time-constant location function with three possibilities: (1)
##' a list of formulae (with parameters having different names) or functions
##' (with one parameter vector numbering for all of them) each returning one
##' value per individual; (2) a single formula or function which will be used
##' for all states (and all categories if multinomial) but with different
##' parameter values in each so that pcmu must be a vector of length the number
##' of unknowns in the function or formula times the number of states (times
##' the number of categories minus one if multinomial); or (3) a function
##' returning an array with one row for each individual, one column for each
##' state of the hidden Markov chain, and, if multinomial, one layer for each
##' category but the last. If used, this function or formula should contain the
##' intercept. Ignored if \code{mu} is supplied.
##' @param tvmu A time-varying location function with three possibilities: (1)
##' a list of formulae (with parameters having different names) or functions
##' (with one parameter vector numbering for all of them) each returning one
##' value per time point; (2) a single formula or function which will be used
##' for all states (and all categories if multinomial) but with different
##' parameter values in each so that ptvmu must be a vector of length the
##' number of unknowns in the function or formula times the number of states
##' (times the number of categories minus one if multinomial); or (3) a
##' function returning an array with one row for each time point, one column
##' for each state of the hidden Markov chain, and, if multinomial, one layer
##' for each category but the last. This function or formula is usually a
##' function of time; it is the same for all individuals. It only contains the
##' intercept if \code{cmu} does not. Ignored if \code{mu} is supplied.
##' @param pgamma A square \code{mxm} matrix of initial estimates of the hidden
##' Markov transition matrix, where \code{m} is the number of hidden states.
##' Rows must sum to one. If the matrix contains zeroes or ones, these are
##' fixed and not estimated. (Ones cannot appear on the diagonal.) If a
##' \code{1x1} matrix or a scalar value of 1 is given, the independence model
##' is fitted.
##' @param pmu Initial estimates of the unknown parameters in \code{mu}.
##' @param pcmu Initial estimates of the unknown parameters in \code{cmu}.
##' @param ptvmu Initial estimates of the unknown parameters in \code{tvmu}.
##' @param pshape Initial estimate(s) of the dispersion parameter, for those
##' distributions having one. This can be one value or a vector with a
##' different value for each state.
##' @param pfamily Initial estimate of the family parameter, for those
##' distributions having one.
##' @param par Initial estimate of the autoregression parameter.
##' @param pintercept For multinomial, proportional odds, and continuation
##' ratio models, \code{p-2} initial estimates for intercept contrasts from the
##' first intercept, where \code{p} is the number of categories.
##' @param delta Scalar or vector giving the unit of measurement (always one
##' for discrete data) for each response value, set to unity by default. For
##' example, if a response is measured to two decimals, delta=0.01. If the
##' response is transformed, this must be multiplied by the Jacobian. For
##' example, with a log transformation, \code{delta=1/response}. Ignored if
##' response has class, \code{response} or \code{repeated}.
##' @param envir Environment in which model formulae are to be interpreted or a
##' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
##' name of the response variable should be given in \code{response}. If
##' \code{response} has class \code{repeated}, it is used as the environment.
##' @param print.level Arguments for nlm.
##' @param typsize Arguments for nlm.
##' @param ndigit Arguments for nlm.
##' @param gradtol Arguments for nlm.
##' @param stepmax Arguments for nlm.
##' @param steptol Arguments for nlm.
##' @param iterlim Arguments for nlm.
##' @param fscale Arguments for nlm.
##' @return A list of classes \code{hidden} and \code{recursive} (unless
##' multinomial, proportional odds, or continuation ratio) is returned that
##' contains all of the relevant information calculated, including error codes.
##' @author J.K. Lindsey and P.J. Lindsey
### @seealso \code{\link[repeated]{chidden}}, \code{\link[repeated]{gar}},
### \code{\link[repeated]{gnlmm}}, \code{\link[rmutil]{iprofile}},
### \code{\link[repeated]{kalcount}}, \code{\link[rmutil]{mprofile}},
### \code{\link[repeated]{nbkal}}, \code{\link[rmutil]{read.list}},
### \code{\link[rmutil]{restovec}}, \code{\link[rmutil]{rmna}}.
##' @references MacDonald, I.L. and Zucchini, W. (1997) Hidden Markov and other
##' Models for Discrete-valued Time Series. Chapman & Hall.
##' @keywords models
##' @examples
##' 
##' # generate two random Poisson sequences with change-points
##' set.seed(8)
##' y <- rbind(c(rpois(5,1), rpois(15,5)), c(rpois(15,1), rpois(5,5)))
##' print(z <- hidden(y,dist="Poisson", cmu=~1, pcmu=c(1,5),
##'                   pgamma=matrix(c(0.95,0.05,0.05,0.95),ncol=2)))
##' # or equivalently
##' mu <- function(p) array(rep(p[1:2],rep(2,2)), c(2,2))
##' print(z <- hidden(y,dist="Poisson", cmu=mu, pcmu=c(1,5),
##'                   pgamma=matrix(c(0.95,0.05,0.05,0.95),ncol=2)))
##' # or
##' # param nind For plotting: numbers of individuals to plot.
##' # param state For plotting: states to plot.
##' print(z <- hidden(y,dist="Poisson", mu=~rep(a,40), pmu=c(1,5),
##'                   pgamma=matrix(c(0.95,0.05,0.05,0.95),ncol=2)))
##' par(mfrow=c(3,2))
##' plot(z, nind=1:2)
##' plot(z, nind=1:2, smooth=TRUE)
##' plot(iprofile(z), lty=2)
##' plot(mprofile(z), add=TRUE)
##' plot(iprofile(z), nind=2, lty=2)
##' plot(mprofile(z), nind=2, add=TRUE)
##' 
##' 
##' @export
hidden <- function(response=NULL, totals=NULL, distribution="Bernoulli",
	mu=NULL, cmu=NULL, tvmu=NULL, pgamma, pmu=NULL, pcmu=NULL, ptvmu=NULL,
	pshape=NULL, pfamily=NULL, par=NULL, pintercept=NULL, delta=NULL,
	envir=parent.frame(), print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, fscale=1, iterlim=100, typsize=abs(p),
	stepmax=10*sqrt(p%*%p)){
#
# likelihood function for nlm
#
likel <- function(p){
	if(ord)for(i in 1:states){
		if(!all(diff(p[np4+l*(i-1)+(1:l)])>0))
			p[np4+l*(i-1)+(1:l)] <- p[np4+l*(i-1)+1]*seq(1,1.1,length=l)}
	if(is.function(muf))pmu0 <- muf(p[nm1:nmu])
	if(is.function(cmuf))pmu1 <- cmuf(p[nm1:ncmu])
	if(is.function(tvmuf))pmu2 <- tvmuf(p[ncmu1:np])
	if(np>np1){
		if(np3-np1==1)pshape <- rep(p[np3],states)
		else pshape <- p[np2:np3]}
	z <- .Fortran("hidden_f",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		nobs=as.integer(nobs(response)),
		mobs=as.integer(mobs),
		s=as.double(if(mdl==4)y else response$response$y),
		n=as.integer(response$response$n),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		pos=as.integer(pos),
		gamma=double(states*states),
		model=as.integer(mdl),
		lgam=as.double(lgam),
		ismu=as.logical(ismu),
		pmu=as.double(pmu0),
		pcmu=as.double(pmu1),
		ptvmu=as.double(pmu2),
		pshape=as.double(pshape),
		pfam=as.double(p[np5]),
		ppar=as.logical(!is.null(par)),
		par=as.double(p[np]),
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
		work2=double(states),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")
z$like}
#
# likelihood function returning parameter values
#
like <- function(p){
	if(is.function(muf))pmu0 <- muf(p[nm1:nmu])
	if(is.function(cmuf))pmu1 <- cmuf(p[nm1:ncmu])
	if(is.function(tvmuf))pmu2 <- tvmuf(p[ncmu1:np])
	if(np>np1){
		if(np3-np1==1)pshape <- rep(p[np3],states)
		else pshape <- p[np2:np3]}
	z <- .Fortran("hidden_f",
		x=as.double(p),
		as.integer(states),
		iq=as.integer(nosubj),
		nobs=as.integer(nobs(response)),
		mobs=as.integer(mobs),
		s=as.double(if(mdl==4)y else response$response$y),
		n=as.integer(response$response$n),
		l=as.integer(l),
		pgamma=as.double(pgamma),
		pos=as.integer(pos),
		gamma=double(states*states),
		model=as.integer(mdl),
		lgam=as.double(lgam),
		ismu=as.logical(ismu),
		pmu=as.double(pmu0),
		pcmu=as.double(pmu1),
		ptvmu=as.double(pmu2),
		pshape=as.double(pshape),
		pfam=as.double(p[np5]),
		ppar=as.logical(!is.null(par)),
		par=as.double(p[np]),
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
		work2=double(states),
		like=double(1),
		#DUP=FALSE,
		PACKAGE="repeated")
z}
call <- sys.call()
#
# check distribution
#
tmp <- c("Bernoulli","Poisson","multinomial","continuation ratio",
	 "proportional odds","binomial","exponential","beta binomial",
	 "negative binomial","normal","inverse Gauss","logistic",
	 "Cauchy","Laplace","Levy","Pareto","gamma","Weibull",
	 "gen gamma","gen logistic","Hjorth","Burr","gen Weibull",
	 "gen extreme value","gen inverse Gauss","power exponential",
	 "skew Laplace","Student t")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
if(mdl>3)mdl <- mdl+1
ord <- distribution=="continuation ratio"||distribution=="proportional odds"
ord2 <- !(!is.null(cmu)&&!is.null(tvmu))
if(!is.null(par)&&(ord||distribution=="Bernoulli"||distribution=="multinomial"))
	stop("par not allowed with multinomial and ordinal responses")
#
# check estimates of transition matrix
#
if(!is.matrix(pgamma)||dim(pgamma)[2]!=dim(pgamma)[1]){
	if(length(pgamma)==1&&pgamma==1)pgamma <- as.matrix(pgamma)
	else stop("pgamma must be a square transition matrix")}
if(any(pgamma<0)||any(pgamma>1))
	stop("All pgamma estimates must be between zero and one")
states <- dim(pgamma)[1]
if(states>1&&any(diag(pgamma)==1))stop("Diagonal of pgamma cannot be 1")
pgamma <- ifelse(pgamma==0,1e-31,pgamma)
pg <- NULL
pos <- rep(1,states)
if(states>1)for(i in 1:states){
	if(round(sum(pgamma[i,]),2)!=1)
		stop(paste("row",i,"of pgamma does not sum to 1"))
	for(j in 1:states){
		if(pgamma[i,j]>1e-30){
			pos[i] <- j
			pdg <- pgamma[i,j]
			break}}
	if(any(pgamma[i,]>=1e-30&pgamma[i,]!=1))
		pg <- c(pg,log(pgamma[i,(1:states)!=pos[i]&pgamma[i,]>=1e-30&pgamma[i,]!=1]/pdg))}
#
# check extra parameters in distribution
#
if(mdl>8){
	if(is.null(pshape))stop("pshape estimate must be supplied")
	else if(length(pshape)!=1&&length(pshape)!=states)
		stop(paste("pshape must have 1 or",states,"estimates"))}
else pshape <- NULL
if(mdl>19){
	if(is.null(pfamily))stop("pfamily estimate must be supplied")
	else if(length(pshape)!=1)
		stop(paste("pshape must have one estimate"))}
else pfamily <- NULL
#
# check if a data object is being supplied
#
type <- "unknown"
respenv <- exists(deparse(substitute(response)),envir=parent.frame())&&
	inherits(response,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(response$response$y)[2]>1){
		if(mdl==3)mdl <- 4
		else stop("hidden only handles univariate responses")}
	if(!is.null(response$NAs)&&any(response$NAs))
		stop("hidden does not handle data with NAs")
	type <- response$response$type}
envname <- if(respenv)deparse(substitute(response))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
#if response is not a data object, make one
#
if(!respenv){
	if(inherits(envir,"repeated")){
		if(!is.null(envir$NAs)&&any(envir$NAs))
			stop("hidden does not handle data with NAs")
		cn <- deparse(substitute(response))
		if(length(grep("\"",cn))>0)cn <- response
		if(length(cn)>1){
			if(mdl==3)mdl <- 4
			else stop("only one response variable allowed")}
		response <- envir
		col <- match(cn,colnames(response$response$y))
		if(is.na(col))stop(paste("response variable",cn,"not found"))
		type <- response$response$type[col]
		if(dim(response$response$y)[2]>1){
			response$response$y <- response$response$y[,col,drop=FALSE]
			if(!is.null(response$response$delta)){
				response$response$delta <- response$response$delta[,col,drop=FALSE]
				if(all(response$response$delta==1)||all(is.na(response$response$delta)))response$response$delta <- NULL}}}
	else {
        	if(!inherits(response,"response")){
        		if(is.vector(response,mode="numeric")){
        			times <- 1:length(response)
        			response <- matrix(response,nrow=1)}
        		else if(is.matrix(response))times <- 1:dim(response)[2]
        		else if(is.array(response)&&length(dim(response))==3&&mdl==3){
				times <- 1:dim(response)[2]
				mdl <- 4}
			else stop("response must be a vector, matrix, or array")
			response <- restovec(response,times=times,totals=totals,delta=delta)
			rm(times)}
		type <- response$type
		response <- rmna(response=response)}}
#
# check if response is appropriate for distribution
#
if(is.null(response$response$times))stop("no times available")
if(distribution=="binomial"||distribution=="beta binomial"){
	if(type!="unknown"&&type!="nominal")stop("nominal data required")
	if(is.null(response$response$n))stop("totals must be supplied")}
if(mdl<11){
	if((distribution=="Poisson"||distribution=="negative binomial")
		&&type!="unknown"&&type!="discrete")
		stop("discrete data required")
	if(distribution=="Bernoulli"){
		if(type!="unknown"&&type!="nominal")
			stop("nominal data required")
		if(any(response$response$y!=0&response$response$y!=1))
			stop("responses must be binary for Bernoulli")}
	if(distribution=="multinomial"
		&&type!="unknown"&&type!="nominal"&&type!="ordinal")
		stop("nominal or ordinal data required")
	if((distribution=="continuation ratio"||
		distribution=="proportional odds")
		&&type!="unknown"&&type!="ordinal")
		stop("ordinal data required")
	if(any(response$response$y<0))
		stop("all responses must be non-negative")}
else if(mdl!=11&&mdl!=13&&mdl!=14&&mdl!=15&&mdl!=16&&mdl!=21&&mdl!=27
	&&mdl!=28&&mdl!=29){
	if(type!="unknown"&&type!="duration"&&type!="continuous")
		stop("duration data required")
	if(any(response$response$y<=0))stop("all responses must be positive")}
else if(type!="unknown"&&type!="continuous"&&type!="duration")
	stop("continuous data required")
lgam <- NULL
if(mdl==3||ord){
	l <- length(unique(response$response$y))-1
	if(min(response$response$y)!=0||max(response$response$y)!=l)
		stop(paste("multinomial and ordinal categories must be numbered from 0 to",l))
	if(l<1)stop("multinomial and ordinal responses must have at least 2 categories")}
else if(mdl==4){
	l <- dim(response$response$y)[2]-1
	lgam <- lgamma(apply(response$response$y,1,sum)+1)-
		apply(lgamma(response$response$y+1),1,sum)
	y <- as.vector(t(response$response$y))}
else l <- 1
if(ord){
	if(length(pintercept)!=l*states)
		stop(paste("pintercept should contain",l*states,"estimates"))
	for(i in 1:states){
		if(!all(diff(pintercept[(i-1)*l+(1:l)])>0))
			stop("intercepts must be monotone increasing for each state")}}
else pintercept <- NULL
if(!is.null(response$response$nest)&&!is.null(par))
	stop("hidden does not handle two levels of nesting with an AR")
if(!is.null(response$response$censor))stop("hidden does not handle censoring")
nosubj <- length(nobs(response))
mobs <- max(nobs(response))
nob <- dim(response$response$y)[1]
ismu <- !is.null(mu)
if(ismu)cmu <- tvmu <- pcmu <- ptvmu <- NULL
else pmu <- NULL
#
# if a data object was supplied, modify formulae or functions to read from it
#
muf <- mu1 <- mu3 <- cmuf <- cmu1 <- cmu3 <- tvmuf <- tvmu1 <- tvmu3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	if(!is.null(mu)){
		cl <- is.list(mu)
		if(!cl)mu <- list(mu)
		tmp2 <- mu3 <- list()
		for(i in 1:length(mu)){
			if(inherits(mu[[i]],"formula")){
				mu3 <- c(mu3,list(if(respenv)finterp(mu[[i]],.envir=response,.name=envname,.intercept=!ord)
					else finterp(mu[[i]],.envir=envir,.name=envname,.intercept=!ord)))
				tmp2 <- c(tmp2,list(mu[[i]]))}
			else if(is.function(mu[[i]])){
				tmp <- parse(text=deparse(mu[[i]])[-1])
				tmp1 <- if(respenv)fnenvir(mu[[i]],.envir=response,.name=envname)
					else fnenvir(mu[[i]],.envir=envir,.name=envname)
				tmp2 <- c(tmp2,list(tmp1))
				attr(tmp1,"model") <- tmp
				mu3 <- c(mu3,list(tmp1))}}
		mu <- tmp2
		rm(tmp2)
		if(!cl){
			mu <- mu[[1]]
			mu3 <- mu3[[1]]}}
	if(!is.null(cmu)){
		cl <- is.list(cmu)
		if(!cl)cmu <- list(cmu)
		tmp2 <- cmu3 <- list()
		for(i in 1:length(cmu)){
			if(inherits(cmu[[i]],"formula")){
				cmu3 <- c(cmu3,list(if(respenv)finterp(cmu[[i]],.envir=response,.name=envname,.intercept=!ord)
					else finterp(cmu[[i]],.envir=envir,.name=envname,.intercept=!ord)))
			tmp2 <- c(tmp2,list(cmu[[i]]))}
			else if(is.function(cmu[[i]])){
				tmp <- parse(text=deparse(cmu[[i]])[-1])
				tmp1 <- if(respenv)fnenvir(cmu[[i]],.envir=response,.name=envname,.expand=FALSE)
					else fnenvir(cmu[[i]],.envir=envir,.name=envname,.expand=FALSE)
				tmp2 <- c(tmp2,list(tmp1))
				attr(tmp1,"model") <- tmp
				cmu3 <- c(cmu3,list(tmp1))}}
		cmu <- tmp2
		rm(tmp2)
		if(!cl){
			cmu <- cmu[[1]]
			cmu3 <- cmu3[[1]]}}
	if(!is.null(tvmu)){
		cl <- is.list(tvmu)
		if(!cl)tvmu <- list(tvmu)
		tmp2 <- tvmu3 <- list()
		for(i in 1:length(tvmu)){
			if(inherits(tvmu[[i]],"formula")){
				tvmu3 <- c(tvmu3,list(if(respenv)finterp(tvmu[[i]],.envir=response,.name=envname,.intercept=!ord)
					else finterp(tvmu[[i]],.envir=envir,.name=envname,.intercept=!ord)))
			tmp2 <- c(tmp2,list(tvmu[[i]]))}
			else if(is.function(tvmu[[i]])){
				tmp <- parse(text=deparse(tvmu[[i]])[-1])
				tmp1 <- if(respenv)fnenvir(tvmu[[i]],.envir=response,.name=envname)
					else fnenvir(tvmu[[i]],.envir=envir,.name=envname)
				tmp2 <- c(tmp2,list(tmp1))
				attr(tmp1,"model") <- tmp
				tvmu3 <- c(tvmu3,list(tmp1))}}
		tvmu <- tmp2
		rm(tmp2)
		if(!cl){
			tvmu <- tvmu[[1]]
			tvmu3 <- tvmu3[[1]]}}}
#
# transform general formulae to functions and check number of parameters
#
npr0 <- if(ord)length(pmu)/states else length(pmu)/states/l
start <- NULL
if(inherits(mu,"formula")){
	mu1 <- if(respenv)finterp(mu,.envir=response,.name=envname,.intercept=!ord)
		else finterp(mu,.envir=envir,.name=envname,.intercept=!ord)
	npt0 <- length(attr(mu1,"parameters"))
	if(is.character(attr(mu1,"model"))){
	# W&R formula
		if(ord&&length(attr(mu1,"model"))==0)mu1 <- NULL}
	else {
	# formula with unknowns
		if(npr0!=npt0){
			cat("\nParameters are ")
			cat(attr(mu1,"parameters"),"\n")
			stop(paste("pmu should have",if(ord) npt0*states else npt0*states*l,"estimates"))}
		if(is.list(pmu))pmu <- unlist(pmu)}
	if(is.null(mu1))pmu <- NULL}
else if(is.list(mu)){
	if(!ord&&length(mu)!=states*l)
		stop(paste("list of mu formulae or functions must have length",states*l))
	if(ord&&length(mu)!=states)
		stop(paste("list of mu formulae or functions must have length",states))
	# only count parameters if list of formulae
	if(inherits(mu[[1]],"formula"))start <- 1
	tmp <- list()
	for(i in mu)if(inherits(i,"formula")){
		mu1 <- if(respenv)finterp(i,.envir=response,.name=envname,.start=start,.intercept=!ord)
			else finterp(i,.envir=envir,.name=envname,.start=start,.intercept=!ord)
		npt0 <- length(attr(mu1,"parameters"))
		if(ord&&is.null(attr(mu1,"model")))
			stop("all formulae must contain covariates")
		if(is.list(pmu))pmu <- unlist(pmu)
		# make new list and count parameters
		tmp <- c(tmp,list(mu1))
		start <- start+npt0}
	else if(is.function(i)){
		if(inherits(i,"formulafn"))start <- c(start,attr(i,"parameters"))
		tmp <- c(tmp,list(i))}
	if(!is.null(start)){
		if(is.numeric(start))start <- start-1
		else start <- length(unique(start))}
	mu1 <- tmp
	rm(tmp)}
else if(!is.null(mu))mu1 <- mu
#
# give appropriate attributes to mu1 for printing
#
if(!is.null(mu1)&&is.null(attr(mu1,"parameters"))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=response))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,.envir=response))
			else attributes(fnenvir(mu1,.envir=envir))}}
#
# if possible, check that correct number of estimates was supplied
#
if(!is.list(mu)){
	nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
		else if(is.null(mu1))NULL
		else npt0
	if(!is.null(nlp)&&nlp!=npr0)
		stop(paste("pmu should have",if(ord)nlp*states else nlp*states*l,"initial estimates"))}
if(!is.null(start)&&start!=length(pmu))
	stop(paste("pmu should have",start,"initial estimates"))
#
# transform time-constant formulae to functions and check number of parameters
#
npr1 <- if(ord)length(pcmu)/states else length(pcmu)/states/l
start <- NULL
if(inherits(cmu,"formula")){
	cmu1 <- if(respenv)finterp(cmu,.envir=response,.expand=FALSE,.name=envname,.intercept=!ord)
		else finterp(cmu,.envir=envir,.expand=FALSE,.name=envname,.intercept=!ord)
	npt1 <- length(attr(cmu1,"parameters"))
	if(is.character(attr(cmu1,"model"))){
	# W&R formula
		if(length(attr(cmu1,"model"))==1&&length(cmu1(pcmu))!=nosubj){
			tmp <- function(p) rep(p[1],nosubj)
			attributes(tmp) <- attributes(cmu1)
			cmu1 <- tmp
			rm(tmp)}
		if(ord&&length(attr(cmu1,"model"))==0)cmu1 <- NULL}
	else {
	# formula with unknowns
		if(npr1!=npt1){
			cat("\nParameters are ")
			cat(attr(cmu1,"parameters"),"\n")
			stop(paste("pcmu should have",if(ord)npt1*states
				else npt1*states*l,"estimates"))}
		if(is.list(pcmu))pcmu <- unlist(pcmu)}
	if(is.null(cmu1))pcmu <- NULL}
else if(is.list(cmu)){
	if(!ord&&length(cmu)!=states*l)
		stop(paste("list of cmu formulae or functions must have length",states*l))
	if(ord&&length(cmu)!=states)
		stop(paste("list of cmu formulae or functions must have length",states))
	# only count parameters if list of formulae
	if(inherits(cmu[[1]],"formula"))start <- 1
	tmp <- list()
	for(i in cmu)if(inherits(i,"formula")){
		cmu1 <- if(respenv)finterp(i,.envir=response,.expand=FALSE,.name=envname,.start=start,.intercept=!ord)
			else finterp(i,.envir=envir,.expand=FALSE,.name=envname,.start=start,.intercept=!ord)
		npt1 <- length(attr(cmu1,"parameters"))
		if(is.character(attr(cmu1,"model"))){
		# W&R formula
			if(ord&&length(attr(cmu1,"model"))==1)cm1 <- NULL}
		else if(is.list(pcmu))pcmu <- unlist(pcmu)
		# make new list and count parameters
		tmp <- c(tmp,list(cmu1))
		start <- start+length(attr(cmu1,"parameters"))}
	else if(is.function(i)){
		if(inherits(i,"formulafn"))start <- c(start,attr(i,"parameters"))
		tmp <- c(tmp,list(i))}
	if(!is.null(start)){
		if(is.numeric(start))start <- start-1
		else start <- length(unique(start))}
	cmu1 <- tmp
	rm(tmp)}
else if(!is.null(cmu))cmu1 <- cmu
#
# give appropriate attributes to cmu1 for printing
#
if(!is.null(cmu1)&&is.null(attr(cmu1,"parameters"))){
	attributes(cmu1) <- if(is.function(cmu)){
		if(!inherits(cmu,"formulafn")){
			if(respenv)attributes(fnenvir(cmu,.envir=response))
			else attributes(fnenvir(cmu,.envir=envir))}
		else attributes(cmu)}
		else {
			if(respenv)attributes(fnenvir(cmu1,.envir=response))
			else attributes(fnenvir(cmu1,.envir=envir))}}
#
# if possible, check that correct number of estimates was supplied
#
if(!is.function(cmu)&&!is.list(cmu)){
	nlp <- if(is.function(cmu1))length(attr(cmu1,"parameters"))
		else if(is.null(cmu1))NULL
		else npt1
	if(!is.null(nlp)&&nlp!=npr1)
		stop(paste("pcmu should have",if(ord)nlp*states else nlp*states*l,"initial estimates"))}
if(!is.null(start)&&start!=length(pcmu))
	stop(paste("pcmu should have",start,"initial estimates"))
#
# transform time-varying formulae to functions and check number of parameters
#
npr2 <- if(ord)length(ptvmu)/states else length(ptvmu)/states/l
start <- NULL
if(inherits(tvmu,"formula")){
	tvmu1 <- if(respenv)finterp(tvmu,.envir=response,.name=envname,.intercept=!ord)
		else finterp(tvmu,.envir=envir,.name=envname,.intercept=!ord)
	npt2 <- length(attr(tvmu1,"parameters"))
	if(is.character(attr(tvmu1,"model"))){
	# W&R formula
		if(length(attr(tvmu1,"model"))==1&&length(tvmu1(ptvmu))!=mobs){
			tmp <- function(p) rep(p[1],mobs)
			attributes(tmp) <- attributes(tvmu1)
			tvmu1 <- tmp
			rm(tmp)}
		else if(ord&&length(attr(tvmu1,"model"))==0)tvmu1 <- NULL}
	else {
	# formula with unknowns
		if(npr2!=npt2){
			cat("\nParameters are ")
			cat(attr(tvmu1,"parameters"),"\n")
			stop(paste("ptvmu should have",if(ord) npt2*states else npt2*states*l,"estimates"))}
		if(is.list(ptvmu))ptvmu <- unlist(ptvmu)}
	if(is.null(tvmu1))ptvmu <- NULL}
else if(is.list(tvmu)){
	if(!ord&&length(tvmu)!=states*l)
		stop(paste("list of tvmu formulae or functions must have length",states*l))
	if(ord&&length(tvmu)!=states)
		stop(paste("list of tvmu formulae or functions must have length",states))
	# only count parameters if list of formulae
	if(inherits(tvmu[[1]],"formula"))start <- 1
	tmp <- list()
	for(i in tvmu)if(inherits(i,"formula")){
		tvmu1 <- if(respenv)finterp(i,.envir=response,.expand=FALSE,.name=envname,.start=start,.intercept=!ord)
			else finterp(i,.envir=envir,.expand=FALSE,.name=envname,.start=start,.intercept=!ord)
		npt2 <- length(attr(tvmu1,"parameters"))
		if(is.character(attr(tvmu1,"model"))){
		# W&R formula (don't know how to change index of p[i])
			if(ord&&length(attr(tvmu1,"model"))==1)tvm1 <- NULL}
		else if(is.list(ptvmu))ptvmu <- unlist(ptvmu)
		# make new list and count parameters
		tmp <- c(tmp,list(tvmu1))
		start <- start+length(attr(tvmu1,"parameters"))}
	else if(is.function(i)){
		if(inherits(i,"formulafn"))start <- c(start,attr(i,"parameters"))
		tmp <- c(tmp,list(i))}
	if(!is.null(start)){
		if(is.numeric(start))start <- start-1
		else start <- length(unique(start))}
	tvmu1 <- tmp
	rm(tmp)}
else if(!is.null(tvmu))tvmu1 <- tvmu
#
# give appropriate attributes to tvmu1 for printing
#
if(!is.null(tvmu1)&&is.null(attr(tvmu1,"parameters"))){
	attributes(tvmu1) <- if(is.function(tvmu)){
		if(!inherits(tvmu,"formulafn")){
			if(respenv)attributes(fnenvir(tvmu,.envir=response))
			else attributes(fnenvir(tvmu,.envir=envir))}
		else attributes(tvmu)}
		else {
			if(respenv)attributes(fnenvir(tvmu1,.envir=response))
			else attributes(fnenvir(tvmu1,.envir=envir))}}
#
# if possible, check that correct number of estimates was supplied
#
if(!is.function(tvmu)&&!is.list(tvmu)){
	nlp <- if(is.function(tvmu1))length(attr(tvmu1,"parameters"))
		else if(is.null(tvmu1))NULL
		else npt2
	if(!is.null(nlp)&&nlp!=npr2)
		stop(paste("ptvmu should have",if(ord)nlp*states else nlp*states*l,"initial estimates"))}
if(!is.null(start)&&start!=length(ptvmu))
	stop(paste("ptvmu should have",start,"initial estimates"))
#
# set up general regression function for likelihood
# if arg==FALSE, can print out general parameter estimates as matrices
#
arg <- TRUE
if(ord&&is.null(mu1))arg <- FALSE
if(is.function(mu1)){
	if(!ord&&is.null(pmu))stop("Initial values of pmu must be supplied")
	tmp <- mu1(pmu)
	if(length(tmp)!=nob)
		stop("mu function or formula must return a vector with one element per observation")
	if(any(is.na(tmp))||any(abs(tmp)==Inf))stop("mu returns NAs or Infs")
	rm(tmp)
	# will be able to print coefficients as matrices
	arg <- FALSE
	# create regression function by repeatedly calling
	# function or formula supplied
	if(distribution=="multinomial"){
		zm <- array(0,c(nob,states,l))
		muf <- function(p){
			k <- 0
			for(j in 1:l)for(i in 1:states){
				k <- k+1
				zm[,i,j] <- mu1(p[(1:npr0)+(k-1)*npr0])}
			zm}}
	else if(ord){
		zm <- array(0,c(nob,states,l))
		muf <- function(p){
			k <- 0
			for(i in 1:states){
				tmp <- mu1(p[(1:npr0)+(i-1)*npr0])
				for(j in 1:l){
					k <- k+1
					zm[,i,j] <- p[np0+k]+tmp}}
#				zm[,i,] <- rep(p[np0+(k+1):(k+l)],rep(nob,l))+
#					rep(mu1(p[(1:npr0)+(i-1)*npr0]),l)
#				k <- k+l}
			zm}}
	else {
		muf <- function(p){
			z <- NULL
			for(i in 1:states)
				z <- cbind(z,mu1(p[(1:npr0)+(i-1)*npr0]))
			z}}
	attributes(muf) <- attributes(mu)}
else if(is.list(mu1)){
	if(is.null(pmu))stop("Initial values of pmu must be supplied")
	for(i in mu1)if(length(i(pmu))!=nob)stop("each mu function or formula must return a vector with one element per observation")
	# create one regression function from list
	if(distribution=="multinomial"){
		zm <- array(0,c(nob,states,l))
		muf <- function(p){
			k <- 0
			for(j in 1:l)for(i in 1:states){
				k <- k+1
				zm[,i,j] <- mu1[[k]](p)}
			zm}}
	else if(ord){
		zm <- array(0,c(nob,states,l))
		muf <- function(p){
			k <- 0
			for(i in 1:states){
				tmp <- mu1[[i]](p)
				for(j in 1:l){
					k <- k+1
					zm[,i,j] <- p[np0+k]+tmp}}
			zm}}
	else {
		muf <- function(p){
			z <- NULL
			for(i in mu1)z <- cbind(z,i(p))
			z}}}
else pmu0 <- NULL
#
# set up time-constant regression function for likelihood
# if arc==FALSE, can print out time-constant parameter estimates as matrices
#
arc <- TRUE
if(is.function(cmu1)){
	if(!ord&&is.null(pcmu))stop("Initial values of pcmu must be supplied")
	# check if array returned or one function/formula for all states
	tmp <- cmu1(pcmu)
	d <- dim(tmp)
	if(length(d)==0){
	# one formula or function
		if(length(tmp)!=nosubj)
			stop("cmu function or formula must return a vector with one element per individual")
		if(any(is.na(tmp))||any(abs(tmp)==Inf))
			stop("cmu returns NAs or Infs")
		rm(tmp)
		# will be able to print coefficients as matrices
		arc <- FALSE
		# create regression function by repeatedly calling
		# function or formula supplied
		if(distribution=="multinomial"){
			zmc <- array(0,c(nosubj,states,l))
			cmuf <- function(p){
				k <- 0
				for(j in 1:l)for(i in 1:states){
					k <- k+1
					zmc[,i,j] <- cmu1(p[(1:npr1)+(k-1)*npr1])}
				zmc}}
		else if(ord){
			zmc <- array(0,c(nosubj,states,l))
			cmuf <- function(p){
				k <- 0
				for(i in 1:states){
					tmp <- cmu1(p[(1:npr1)+(i-1)*npr1])
					for(j in 1:l){
						k <- k+1
						zmc[,i,j] <- p[np0+k]+tmp}}
				zmc}}
		else {
			cmuf <- function(p){
				z <- NULL
				for(i in 1:states)
					z <- cbind(z,cmu1(p[(1:npr1)+(i-1)*npr1]))
				z}}
		attributes(cmuf) <- attributes(cmu)}
	else {
	# returns an array: check dimensions
		if(d[2]!=states)stop(paste("cmu must return a",states,"column array"))
		if(d[1]!=nosubj)
			stop("cmu must return an array with one row per individual")
		if(distribution=="multinomial"){
			if(length(d)!=3)
				stop("cmu must return a 3 dimensional array")
			if(d[3]!=l)
				stop(paste("cmu must return an array with",l,"layers"))}
		else if(length(d)!=2)
				stop("cmu must return a 2 dimensional array")
		arc <- FALSE
		cmuf <- cmu1}}
else if(is.list(cmu1)){
	if(is.null(pcmu))stop("Initial values of pcmu must be supplied")
	for(i in cmu1)if(length(i(pcmu))!=nosubj)stop("each cmu function or formula must return a vector with one element per individual")
	# create one regression function from list
	if(distribution=="multinomial"){
		zmc <- array(0,c(nosubj,states,l))
		cmuf <- function(p){
			k <- 0
			for(j in 1:l)for(i in 1:states){
				k <- k+1
				zmc[,i,j] <- cmu1[[k]](p)}
			zmc}}
	else if(ord){
		zmc <- array(0,c(nosubj,states,l))
		cmuf <- function(p){
			k <- 0
			for(i in 1:states){
				tmp <- cmu1[[i]](p)
				for(j in 1:l){
					k <- k+1
					zmc[,i,j] <- p[np0+k]+tmp}}
			zmc}}
	else {
		cmuf <- function(p){
			z <- NULL
			for(i in cmu1)z <- cbind(z,i(p))
			z}}}
else pmu1 <- array(0,c(states,nosubj,l))
#
# set up time-varying regression function for likelihood
# if arv==FALSE, can print out time-varying parameter estimates as matrices
#
arv <- TRUE
if(is.function(tvmu1)){
	if(!ord&&is.null(ptvmu))stop("Initial values of ptvmu must be supplied")
	# check if array returned or one function/formula for all states
	tmp <- tvmu1(ptvmu)
	d <- dim(tmp)
	if(length(d)==0){
	# one formula or function
		if(length(tmp)!=mobs)
			stop("tvmu function or formula must return a vector with one element per time point")
		if(any(is.na(tmp))||any(abs(tmp)==Inf))
			stop("tvmu returns NAs or Infs")
		rm(tmp)
		# will be able to print coefficients as matrices
		arv <- FALSE
		# create regression function by repeatedly calling
		# function or formula supplied
		if(distribution=="multinomial"){
			zmv <- array(0,c(mobs,states,l))
			tvmuf <- function(p){
				k <- 0
				for(j in 1:l)for(i in 1:states){
					k <- k+1
					zmv[,i,j] <- tvmu1(p[(1:npr2)+(k-1)*npr2])}
				zmv}}
		else if(ord){
			zmv <- array(0,c(mobs,states,l))
			tvmuf <- function(p){
				k <- 0
				for(i in 1:states){
					tmp <- tvmu1(p[(1:npr2)+(i-1)*npr2])
					for(j in 1:l){
						k <- k+1
						zmv[,i,j] <- ord2*p[np0+k]+tmp}}
				zmv}}
		else {
			tvmuf <- function(p){
				z <- NULL
				for(i in 1:states)
					z <- cbind(z,tvmu1(p[(1:npr2)+(i-1)*npr2]))
				z}}
		attributes(tvmuf) <- attributes(tvmu)}
	else {
	# returns an array: check dimensions
		if(d[2]!=states)stop(paste("tvmu must return a",states,"column array"))
		else if(d[1]!=mobs)
			stop("tvmu must return an array with one row per time point")
		if(distribution=="multinomial"){
			if(length(d)!=3)
				stop("tvmu must return a 3 dimensional array")
			if(d[3]!=l)
				stop(paste("tvmu must return an array with",l,"layers"))}
		else if(length(d)!=2)
				stop("tvmu must return a 2 dimensional array")
		arv <- FALSE
		tvmuf <- tvmu1}}
else if(is.list(tvmu1)){
	if(is.null(ptvmu))stop("Initial values of ptvmu must be supplied")
	for(i in tvmu1)if(length(i(ptvmu))!=mobs)stop("each tvmu function or formula must return a vector with one element per time point")
	# create one regression function from list
	if(distribution=="multinomial"){
		zmv <- array(0,c(mobs,states,l))
		tvmuf <- function(p){
			k <- 0
			for(j in 1:l)for(i in 1:states){
				k <- k+1
				zmv[,i,j] <- tvmu1[[k]](p)}
			zmv}}
	else if(ord){
		zmv <- array(0,c(mobs,states,l))
		tvmuf <- function(p){
			k <- 0
			for(i in 1:states){
				tmp <- tvmu1[[i]](p)
				for(j in 1:l){
					k <- k+1
					zmv[,i,j] <- ord2*p[np0+k]+tmp}}
			zmv}}
	else {
		tvmuf <- function(p){
			z <- NULL
			for(i in tvmu1)z <- cbind(z,i(p))
			z}}}
else pmu2 <- array(0,c(states,mobs,l))
if(!is.function(muf)&&!is.function(cmuf)&&!is.function(tvmuf)){
	if(ord){
		if(is.null(pintercept))
			stop(paste("pintercept must contain ",l," initial estimates",sep=""))
		else {
			ismu <- TRUE
			zm <- array(0,c(nob,states,l))
			muf <- function(p){
				k <- 0
				for(i in 1:states){
					tmp <- rep(1,nob)
					for(j in 1:l){
						k <- k+1
						zm[,i,j] <- p[np0+k]*tmp}}
				zm}}}
	else stop("Either mu, cmu or tvmu must be a list, function, or formula")}
#
# prepare to call nlm
#
nm1 <- length(pg)+1
nmu <- length(pg)+length(pmu)+length(pintercept)
ncmu <- length(pg)+length(pcmu)+length(ptvmu)+length(pintercept)
ncmu1 <- length(pg)+length(pcmu)+1
p <- c(pg,pmu,pcmu,ptvmu,pshape,pfamily,pintercept,par)
np <- length(p)
np0 <- np-length(pg)-length(pintercept)-length(par)
np1 <- np-length(pshape)-length(pfamily)-length(pintercept)-length(par)
np2 <- np1+1
np3 <- np-length(pfamily)-length(pintercept)-length(par)
np4 <- np-length(pintercept)-length(par)
np5 <- np-length(par)
pshape <- rep(0,states)
z <- nlm(likel,p, hessian=TRUE, print.level=print.level,
	typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
#
# prepare mle's
#
delta <- if(!is.null(response$response$delta))
	sum(log(response$response$delta)) else 0
maxlike <- z$minimum - delta
z0 <- like(z$estimate)
#
# calculate se's
#
a <- if(any(is.na(z$hessian))||any(abs(z$hessian)==Inf))0
	else qr(z$hessian)$rank
if(a==np)cov <- solve(z$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
#
# calculate predicted values
#
if(mdl<3||mdl>6){
	rpred <- 0
	if(is.function(muf)){
		rpred <- rpred+t(muf(z$estimate[nm1:nmu]))
		if(dim(rpred)[2]!=nob)rpred <- rpred[,rep(1:nosubj,nobs(response))]}
	if(is.function(cmuf))rpred <- rpred+t(cmuf(z$estimate[nm1:ncmu]))[,rep(1:nosubj,nobs(response))]
	if(is.function(tvmuf))rpred <- rpred+t(tvmuf(z$estimate[ncmu1:np]))[,sequence(nobs(response))]
	if(mdl==1||mdl==7)rpred <- 1/(1+exp(-rpred))
	pred <- rep(1,states)%*%(rpred*z0$delta)
	rpred <- rep(1,states)%*%(rpred*matrix(z0$filter,nrow=states))}
else pred <- rpred <- NULL
#
# return appropriate attributes on functions
#
if(!is.null(mu3))mu1 <- mu3
if(!is.null(cmu3))cmu1 <- cmu3
if(!is.null(tvmu3))tvmu1 <- tvmu3
z1 <- list(
	call=call,
	type="discrete",
	distribution=distribution,
	response=response$response,
	maxlike=maxlike,
	aic=maxlike+np,
	df=length(response$response$y)-np,
	states=states,
	levels=l,
	mu=mu1,
	arg=arg,
	cmu=cmu1,
	arc=arc,
	arv=arv,
	tvmu=tvmu1,
	coef=z$estimate[nm1:np1],
	pshape=if(np3>np1&&!ord)z$estimate[np2:np3] else NULL,
	pfamily=if(np>np3+!is.null(par)&&!ord)z$estimate[np] else NULL,
	pintercept=if(is.null(pintercept)) NULL else z$estimate[np2:np],
	nmu=length(pmu),
	ncmu=length(pcmu),
	ntvmu=length(ptvmu),
	gamma=matrix(z0$gamma,ncol=states),
	delta=z0$delta,
	par=if(!is.null(par))z0$par else NULL,
	np=np,
	cov=cov,
	corr=corr,
	se=se,
	filter=matrix(z0$filter,nrow=states),
	pred=pred,
	rpred=rpred,
	iterations=z$iter,
	code=z$code)
class(z1) <- if(is.null(rpred))"hidden" else c("hidden","recursive")
z1}

##' @export 
print.hidden <- function(x,digits=max(3,.Options$digits-3),correlation=FALSE, ...){
  z <- x
ord <- z$distribution=="continuation ratio"||
	z$distribution=="proportional odds"
m <- z$states
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(m>1)cat(z$type,"")
cat(z$distribution)
if(m>1)cat(" hidden Markov chain with",m,"states\n\n")
else cat(" independence model\n\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("Number of subjects    ",length(nobs(z)),"\n")
cat("Number of observations",length(z$response$y),"\n")
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n")
if(ord){
# intercept coefficients for ordinal models
	cat("\nIntercepts\n")
	# first check that intercepts are valid
	for(i in 1:m){
		tmp <- diff(z$pintercept[(i-1)*z$levels+(1:z$levels)])
		if(!all(tmp>0))cat(paste("intercepts not monotone increasing for state",i,"\n"))}
	coef <- matrix(z$pintercept,byrow=TRUE,nrow=m)
	dimnames(coef) <- list(paste("state ",1:m,sep=""),
		(1:z$levels)-1)
	print.default(coef,digits=digits,print.gap=2)}
if(z$nmu>0){
# general coefficients
	cat("\nRegression coefficients\n")
	if(!is.list(z$mu))z$mu <- list(z$mu)
	ll <- j <- l <- 0
	cname <- NULL
	for(i in z$mu){
	# print out functions/formulae used, by state if possible
		cname <- c(cname,if(is.character(attr(i,"model")))
			attr(i,"model") else attr(i,"parameters"))
		ll <- c(ll,length(cname))
		l <- l+1
		if(l==1&&length(z$mu)>1&&((!is.null(attr(i,"formula"))&&!is.null(attr(i,"model")))||is.function(i))){
			j <- j+1
			tmp <- paste("state",j)
			cat(tmp,"\n")}
		if(ord||l==z$levels)l <- 0
		if(!is.null(attr(i,"formula")))
			cat(deparse(attr(i,"formula")),sep="\n")
		else if(!is.null(attr(i,"model"))){
			t <- deparse(attr(i,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		else if(is.function(i))cat(deparse(i)[-1],sep="\n")}
	coef <- matrix(z$coef[1:z$nmu],nrow=1)
	if(z$arg){
		if(is.null(cname)||length(cname)!=z$nmu){
		# print out coefficients as one vector
			dimnames(coef) <- list(1,1:z$nmu)
			print.default(coef,digits=digits,print.gap=2)}
		else {
		# print out coefficients by state
			l2 <- if(ord)1 else z$levels
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:l2){
					coef1 <- coef[,(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1],drop=FALSE]
					dimnames(coef1) <- list(1,cname[(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1]])
					print.default(coef1,digits=digits,print.gap=2)}}}}
	else {
	# print out coefficients as matrices
		if(z$distribution=="multinomial"){
			npr0 <- z$nmu/m/z$levels
			if(is.null(cname)||length(cname)!=npr0)
				cname <- 1:npr0
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:z$levels)
					coef1 <- rbind(coef1,coef[(1:npr0)+(i-1)*npr0+(j-1)*npr0*m])
				dimnames(coef1) <- list(1:z$levels,cname)
				print.default(coef1,digits=digits,print.gap=2)}}
		else {
			npr0 <- z$nmu/m
			if(is.null(cname)||length(cname)!=npr0)
				cname <- 1:npr0
			coef1 <- matrix(coef,nrow=m,byrow=TRUE)
			dimnames(coef1) <- list(paste("state",1:m),cname)
			print.default(coef1,digits=digits,print.gap=2)}}}
if(z$ncmu>0){
# time-constant coefficients
	cat("\nTime-constant coefficients\n")
	if(!is.list(z$cmu))z$cmu <- list(z$cmu)
	ll <- j <- l <- 0
	cname <- NULL
	for(i in z$cmu){
	# print out functions/formulae used, by state if possible
		cname <- c(cname,if(is.character(attr(i,"model")))
			attr(i,"model") else attr(i,"parameters"))
		ll <- c(ll,length(cname))
		l <- l+1
		if(l==1&&length(z$cmu)>1&&((!is.null(attr(i,"formula"))&&!is.null(attr(i,"model")))||is.function(i))){
			j <- j+1
			tmp <- paste("state",j)
			cat(tmp,"\n")}
		if(ord||l==z$levels)l <- 0
		if(!is.null(attr(i,"formula")))
			cat(deparse(attr(i,"formula")),sep="\n")
		else if(!is.null(attr(i,"model"))){
			t <- deparse(attr(i,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		else if(is.function(i))cat(deparse(i)[-1],sep="\n")}
	coef <- matrix(z$coef[1:z$ncmu],nrow=1)
	if(z$arc){
		if(is.null(cname)||length(cname)!=z$ncmu){
		# print out coefficients as one vector
			dimnames(coef) <- list(1,1:z$ncmu)
			print.default(coef,digits=digits,print.gap=2)}
		else {
		# print out coefficients by state
			l2 <- if(ord)1 else z$levels
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:l2){
					coef1 <- coef[,(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1],drop=FALSE]
					dimnames(coef1) <- list(1,cname[(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1]])
					print.default(coef1,digits=digits,print.gap=2)}}}}
	else {
	# print out coefficients as matrices
		if(z$distribution=="multinomial"){
			npr1 <- z$ncmu/m/z$levels
			if(is.null(cname)||length(cname)!=npr1)
				cname <- 1:npr1
			k <- 0
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:z$levels){
					k <- k+1
					coef1 <- rbind(coef1,coef[(1:npr1)+(i-1)*npr1+(j-1)*npr1*m])}
				dimnames(coef1) <- list(1:z$levels,cname)
				print.default(coef1,digits=digits,print.gap=2)}}
		else {
			npr1 <- z$ncmu/m
			if(is.null(cname)||length(cname)!=npr1)
				cname <- 1:npr1
			coef1 <- matrix(coef,nrow=m,byrow=TRUE)
			dimnames(coef1) <- list(paste("state",1:m),cname)
			print.default(coef1,digits=digits,print.gap=2)}}}
if(z$ntvmu>0){
# time-varying coefficients
	cat("\nTime-varying coefficients\n")
	if(!is.list(z$tvmu))z$tvmu <- list(z$tvmu)
	ll <- j <- l <- 0
	cname <- NULL
	for(i in z$tvmu){
	# print out functions/formulae used, by state if possible
		cname <- c(cname,if(is.character(attr(i,"model")))
			attr(i,"model") else attr(i,"parameters"))
		ll <- c(ll,length(cname))
		l <- l+1
		if(l==1&&length(z$tvmu)>1&&((!is.null(attr(i,"formula"))&&!is.null(attr(i,"model")))||is.function(i))){
			j <- j+1
			tmp <- paste("state",j)
			cat(tmp,"\n")}
		if(ord||l==z$levels)l <- 0
		if(!is.null(attr(i,"formula")))
			cat(deparse(attr(i,"formula")),sep="\n")
		else if(!is.null(attr(i,"model"))){
			t <- deparse(attr(i,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		else if(is.function(i))cat(deparse(i)[-1],sep="\n")}
	coef <- matrix(z$coef[(z$ncmu+1):(z$ncmu+z$ntvmu)],nrow=1)
	if(z$arv){
		if(is.null(cname)||length(cname)!=z$ntvmu){
		# print out coefficients as one vector
			dimnames(coef) <- list(1,1:z$ntvmu)
			print.default(coef,digits=digits,print.gap=2)}
		else {
		# print out coefficients by state
			l2 <- if(ord)1 else z$levels
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:l2){
					coef1 <- coef[,(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1],drop=FALSE]
					dimnames(coef1) <- list(1,cname[(ll[i+(j-1)*m]+1):ll[i+(j-1)*m+1]])
					print.default(coef1,digits=digits,print.gap=2)}}}}
	else {
	# print out coefficients as matrices
		if(z$distribution=="multinomial"){
			npr2 <- z$ntvmu/m/z$levels
			if(is.null(cname)||length(cname)!=npr2)
				cname <- 1:npr2
			k <- 0
			for(i in 1:m){
				if(m>1)cat("state",i,"\n")
				coef1 <- NULL
				for(j in 1:z$levels){
					k <- k+1
					coef1 <- rbind(coef1,coef[(1:npr2)+(i-1)*npr2+(j-1)*npr2*m])}
				dimnames(coef1) <- list(1:z$levels,cname)
				print.default(coef1,digits=digits,print.gap=2)}}
		else {
			npr2 <- z$ntvmu/m
			if(is.null(cname)||length(cname)!=npr2)
				cname <- 1:npr2
			coef1 <- matrix(coef,nrow=m,byrow=TRUE)
			dimnames(coef1) <- list(paste("state",1:m),cname)
			print.default(coef1,digits=digits,print.gap=2)}}}
if(!is.null(z$pshape)){
	cat("\nDispersion parameters\n")
	coef <- matrix(z$pshape,nrow=1)
	dimnames(coef) <- list(1,1:length(z$pshape))
	print.default(coef, digits=digits, print.gap=2)}
if(!is.null(z$pfamily)){
	cat("\nFamily parameter\n")
	coef <- matrix(z$pfamily,nrow=1)
	dimnames(coef) <- list(1,1)
	print.default(coef, digits=digits, print.gap=2)}
if(z$states>1){
	cat("\nTransition matrix\n")
	dimnames(z$gamma) <- list(1:m,1:m)
	print.default(z$gamma, digits=digits, print.gap=2)
	if(is.null(z$changepoint)){
		cat("\nStationary distribution\n")
		z$delta <- matrix(z$delta,nrow=1)
		dimnames(z$delta) <- list(1,1:m)
		print.default(z$delta, digits=digits, print.gap=2)}}
if(!is.null(z$par)){
	cat("\nAutoregression parameter\n")
	if(z$type=="discrete"){
		coef <- matrix(z$par,nrow=1)
		dimnames(coef) <- list(1,1)}
	else {
		coef <- matrix(c(log(z$par),z$par),nrow=1)
		dimnames(coef) <- list(1,c("estimate","parameter"))}
	print.default(coef, digits=digits, print.gap=2)}
if(correlation){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,z$np),seq(1,z$np))
	print.default(z$corr,digits=digits)}}

# plot probability of being in each state at each time point
##' @export
plot.hidden <- function(x, nind=1, state=NULL, smooth=FALSE, main=NULL,
ylab=NULL, xlab="Time", xlim=NULL, ...){
  z <- x
if(max(nind)>length(nobs(z)))stop("no such individual")
ns <- length(nind)
ii <- covind(z$response)
if(is.null(main)){
	main <- NULL
	for(i in nind)main <- c(main,paste("Individual ",i))}
else if(length(main)!=ns){
	if(length(main==1))main <- rep(main,ns)
	else stop("main must have a name for each individual")}
if(is.null(xlim))xlim <- c(0,max(z$resp$times))
oldpar <- if(z$states>2&&is.null(state))par(mfcol=c(z$states,ns))
	else NULL
if(is.null(state)&&z$state==2)state <- 1
k <- 0
for(i in nind){
	k <- k+1
	lenf <- sum(ii==i)
	for(j in 1:z$states)if(is.null(state)||any(j==state)){
		xl <- if(j==z$states||z$states==2||!is.null(state))xlab else ""
		mn <- if(j==1)main[k] else ""
		yl <- if(k==1||z$states==2||!is.null(state)){
			if(is.null(ylab))paste("State ",j) else ylab}
		else ""
		if(smooth){
			y <- (z$filter[j,ii==i][1:(lenf-2)]+
				z$filter[j,ii==i][3:lenf])*0.25+
				z$filter[j,ii==i][2:(lenf-1)]*0.5
			y <- c(z$filter[j,ii==i][1],y,
				z$filter[j,ii==i][lenf])}
		else y <- z$filter[j,ii==i]
		plot(z$resp$times[ii==i],y,type="l",ylab=yl,
			xlab=xl,main=mn,xlim=xlim,ylim=c(0,1),...)}}
if(!is.null(oldpar))par(oldpar)}
