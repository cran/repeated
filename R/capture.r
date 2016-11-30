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
#     capture(z,n)
#
#  DESCRIPTION
#
#    Functions to fit capture-recapture models.

### function to calculate capture-recapture parameters from a Poisson glm
###


##' Capture-recapture Models
##' 
##' \code{capture} fits the Cormack capture-recapture model to \code{n} sample
##' periods. Set \code{n} to the appropriate value and type \code{eval(setup)}.
##' 
##' \code{n <- periods} # number of periods
##' 
##' \code{eval(setup)}
##' 
##' This produces the following variables -
##' 
##' \code{p[i]}: logit capture probabilities,
##' 
##' \code{pbd}: constant capture probability,
##' 
##' \code{d[i]}: death parameters,
##' 
##' \code{b[i]}: birth parameters,
##' 
##' \code{pw}: prior weights.
##' 
##' Then set up a Poisson model for log linear models:
##' 
##' \code{z <- glm(y~model, family=poisson, weights=pw)}
##' 
##' and call the function, \code{capture}.
##' 
##' If there is constant effort, then all estimates are correct. Otherwise,
##' \code{n[1]}, \code{p[1]}, \code{b[1]}, are correct only if there is no
##' birth in period 1.  \code{n[s]}, \code{p[s]}, are correct only if there is
##' no death in the last period.  \code{phi[s-1]} is correct only if effort is
##' constant in \code{(s-1, s)}.  \code{b[s-1]} is correct only if \code{n[s]}
##' and \code{phi[s-1]} both are.
##' 
##' 
##' @param z A Poisson generalized linear model object.
##' @param n The number of repeated observations.
##' @aliases setup
##' @return \code{capture} returns a matrix containing the estimates.
##' @author J.K. Lindsey
##' @keywords models
##' @examples
##' 
##' y <- c(0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,14,1,1,0,2,1,2,1,16,0,2,0,11,
##' 	2,13,10,0)
##' n <- 5
##' eval(setup)
##' # closed population
##' print(z0 <- glm(y~p1+p2+p3+p4+p5, family=poisson, weights=pw))
##' # deaths and emigration only
##' print(z1 <- update(z0, .~.+d1+d2+d3))
##' # immigration only
##' print(z2 <- update(z1, .~.-d1-d2-d3+b2+b3+b4))
##' # deaths, emigration, and immigration
##' print(z3 <- update(z2, .~.+d1+d2+d3))
##' # add trap dependence
##' print(z4 <- update(z3, .~.+i2+i3))
##' # constant capture probability over the three middle periods
##' print(z5 <- glm(y~p1+pbd+p5+d1+d2+d3+b2+b3+b4, family=poisson, weights=pw))
##' # print out estimates
##' capture(z5, n)
##' 
##' @export capture
capture <- function(z,n){
	aft <- bef <- rep(1,2^n-1)
	aa <- bb <- m <- rep(1,n)
	for(i in (1-n):-1){
		m[n+i] <- z$fit[1]+z$fit[2^(n+i-1)+1]
		aft <- aft*(2-rep(rep(1:2,rep(2^(n+i-1),2)),2^(-i))[1:(2^n-1)])
		aa[n+i] <- sum(aft*z$fit[1:(2^n-1)])
		bef <- bef*(2-rep(rep(1:2,rep(2^(-i),2)),2^(n+i-1))[1:(2^n-1)])
		bb[1-i] <- sum(bef*z$fit[1:(2^n-1)])}
	aa[n] <- bb[1] <- z$fit[1]
	cc <- aa[c(n,1:(n-1))]
	m[n] <- cc[n]
	ph <- z$fit[1]/m
	nest <- aa*m*bb/z$fit[1]^2
	phih <- aa*m/(cc*z$fit[1])
	if(max(ph[2:(n-1)])-min(ph[2:(n-1)])<=10^(-4)){
		nest <- nest*ph/ph[2]
		ph[1] <- ph[2]
		ph[n] <- ph[2]
		phih[n] <- z$fit[1]/(m[n]*ph[n])}
	phih[1] <- 1
	phih[n] <- phih[n]+(m[n-1]/m[n]-1)*10^(-4)*(phih[n-1]<=1-10^(-4))
	bh <- nest-phih*nest[c(n,1:(n-1))]
	bh[1] <- nest[1]
	zz <- cbind(1:n,nest,phih,ph,bh)
	colnames(zz) <- c("i","N(i)","Phi(i-1)","P(i)","B(i-1)")
	zz}

### produce the required variables for the number of periods, n, specified
### 
##' @export setup
setup <- expression({
	p1 <- as.numeric(gl(2,1,2^n))-1
	pbd <- p2 <- as.numeric(gl(2,2,2^n))-1
	p3 <- as.numeric(gl(2,4,2^n))-1
	b2 <- i1 <- p1*p2
	d1 <- i2 <- p2*p3
	pb <- pbd+p3
	if(n>=4){
		p4 <- as.numeric(gl(2,8,2^n))-1
		d2 <- i3 <- p3*p4
		b3 <- b2*p3
		d1 <- d2*p2
		pbd <- pbd+p3
		pb <- pbd+p4}
	if(n>=5){
		p5 <- as.numeric(gl(2,16,2^n))-1
		d3 <- i4 <- p4*p5
		b4 <- b3*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p4
		pb <- pbd+p5}
	if(n>=6){
		p6 <- as.numeric(gl(2,32,2^n))-1
		d4 <- i5 <- p5*p6
		b5 <- b4*p5
		d3 <- d4*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p5
		pb <- pbd+p6}
	if(n>=7){
		p7 <- as.numeric(gl(2,64,2^n))-1
		d5 <- i6 <- p6*p7
		b6 <- b5*p6
		d4 <- d5*p5
		d3 <- d4*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p6
		pb <- pbd+p7}
	pw <- rep(1,2^n)
	pw[2^n] <- 0
	pd <- pbd+p1
	pc <- pb+p1})

