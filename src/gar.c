/*
 *  repeated : A Library of Repeated Measurements Models
 *  Copyright (C) 1998 J.K. Lindsey
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 *  void gar(double y[], double total[], int *my, int nobs[], int *nind,
 *	 double times[], int censor[], int *cens, double eta[],
 *	 double theta[], int *model, int *thp, double shape[], int *sh,
 *	 int *link, int* ar, int* order, double pred[], double rpred[],
 *	 double *like)
 *
 *  DESCRIPTION
 *
 *    Function to compute the likelihood function for generalized nonlinear
 *  autoregression models with various distributions.
 *
 */

#include <math.h>
#include <stddef.h>
#include "Mathlib.h"

extern void dmp(int y[], int *my, double m[], double s[], int *nn,
       double wt[], double res[]);
extern void ddp(int y[], int *my, double m[], double s[], int *nn,
       double wt[], double res[]);
extern void dmb(int y[], int n[], double m[], double s[], int *nn,
       double wt[], double res[]);
extern void ddb(int y[], int n[], double m[], double s[], int *nn,
       double wt[], double res[]);
extern void plevy(double y[], double m[], double s[], double f[], int *len,
       double *eps, int *pts, int *max, int *err, double res[]);
extern void pginvgauss(double y[], double m[], double s[], double f[],
       int *len, double *eps, int *pts, int *max, int *err, double res[]);
extern void ppowexp(double y[], double m[], double s[], double f[], int *len,
	   double *eps, int *pts, int *max, int *err, double res[]);


extern double lgamma(double x);
extern double lbeta(double a, double b);
extern double lchoose(double n, double k);
extern double dbinom(double x, double n, double p);
extern double dpois(double x, double lambda);
extern double dexp(double x, double scale);
extern double pexp(double x, double scale);
extern double dweibull(double x, double shape, double scale);
extern double pweibull(double x, double shape, double scale);
extern double dgamma(double x, double shape, double scale);
extern double pgamma(double x, double shape, double scale);
extern double dnorm(double x, double mean, double sd);
extern double pnorm(double x, double mean, double sd);
extern double dlogis(double x, double location, double scale);
extern double plogis(double x, double location, double scale);
extern double dcauchy(double x, double location, double scale);
extern double pcauchy(double x, double location, double scale);
extern double bessel_k(double x, double alpha, double expo);

void gar(double y[], double total[], int *my, int nobs[], int *nind,
	 double times[], int censor[], int *cens, double eta[],
	 double theta[], int *model, int *thp, double shape[], int *sh,
	 int *link, int* ar, int* order, double pred[], double rpred[],
	 double *like){
  int nm,ii,jj,pos,bin,nn,iy,iy2,pts,max,err;
  double cres,cresp,cw,cwp,lmhu,lmhu2,diff,diff2,tmp,tmp2,lik,lambda,lambda2,
    th1,th2,th3,cres2,cresp2,wt,eps;

  nn=1;
  wt=1.0;
  pts=5;
  max=16;
  eps=1.0e-6;
  pos=(*model!=10)&&(*model!=11)&&(*model!=12)&&(*model!=15)&&(*model!=18)&&(*model!=20)&&(*model!=26);
  bin=(*model==1)||(*model==7)||(*model==8)||(*model==9);
  *like=0.;
  if(*thp){
    th1=theta[0];
    th2=theta[1];}
  else th2=theta[0];
  if(*order==2)th3=theta[1];
  if(*ar>2){
    th2=exp(-th2);
    if(*order==2)th3=exp(-th3);}
  if(*model>3&&!*sh)lambda=*model!=26?exp(theta[*thp+*order]):exp(0.5*theta[*thp+*order]);
  if(*model>16)lambda2=*sh?theta[*thp+*order]:exp(theta[*thp+*order+1]);
  if(*model==26)tmp=1.0+1.0/(2.0*lambda2);
  nm=0;
  for(ii=0;ii<*nind;ii++){
    cres=cres2=cw=lmhu=lmhu2=0.;
    for(jj=0;jj<nobs[ii];jj++){
      if(*sh&&*model>3)lambda=*model!=26?shape[nm]:sqrt(shape[nm]);
      if(jj>0){
	diff=times[nm]-times[nm-1];
	cresp=cres;
	cwp=cw;
	if(*thp){
	  cw=exp(-th1*diff);
	  cres=cw*cresp;
	  cw=cw*cwp;}
	else cres=cw=0.;
	switch(*ar){
	case 1: lmhu=exp(-th2*diff)*cresp/cwp; break;
	case 2: lmhu=exp(-th2*diff*diff)*cresp/cwp; break;
	case 3: lmhu=cresp/cwp/(1+th2*diff*diff); break;
	case 4: lmhu=diff<=1/th2?(pow(diff*th2,3)-3*th2*diff+2)*cresp/cwp/2:0;
	  break;
	case 5: lmhu=(2*th2*times[nm-1]+exp(-th2*times[nm])
		  +exp(-th2*times[nm-1])-1-exp(-th2*diff))/(2*pow(th2,3))
		  *cresp/cwp; break;}
	if(*order==2&&jj>1){
	  diff2=times[nm-1]-times[nm-2];
	  cresp2=cres2;
	  switch(*ar){
	  case 1: lmhu2=exp(-th3*diff2)*cresp2; break;
	  case 2: lmhu2=exp(-th3*diff2*diff2)*cresp2; break;
	  case 3: lmhu2=cresp2/(1+th3*diff2*diff2); break;
	  case 4: lmhu2=diff2<=1/th3?(pow(diff2*th3,3)-3*th3*diff2+2)*cresp2/2:0;
	    break;
	  case 5: lmhu2=(2*th3*times[nm-1]+exp(-th3*times[nm])
		    +exp(-th3*times[nm-1])-1-exp(-th3*diff2))/
		    (2*pow(th3,3))*cresp2; break;}}
	cres2=cresp;}
      cw++;
      switch(*link){
      case 1: pred[nm]=eta[nm]; break;
      case 2: pred[nm]=log(eta[nm]); break;
      case 3: pred[nm]=sqrt(eta[nm]); break;
      case 4: pred[nm]=eta[nm]*eta[nm]; break;
      case 5: pred[nm]=exp(eta[nm]); break;
      case 6:
	pred[nm]=exp(eta[nm]);
	pred[nm]=pred[nm]/(1+pred[nm]); break;
      case 7: pred[nm]=exp(-exp(eta[nm])); break;}
      lmhu+=pred[nm]+lmhu2;
      if(pos&&lmhu<=0.0)lmhu=0.01;
      if(*model==18&&lmhu>=y[nm])lmhu=y[nm]-0.01;
      if(bin){
	if(lmhu>=1.)lmhu=0.99;
	rpred[nm]=lmhu*total[nm];
	cres+=y[nm]-pred[nm]*total[nm];}
      else {
	rpred[nm]=lmhu;
	cres+=y[nm]-pred[nm];}
      if(!*cens||censor[nm]==1){
	switch(*model){
	case 1: /* binary distribution */
	  if(total[nm]==1)*like-=(int)y[nm]?log(lmhu):log(1-lmhu);
	  else *like-=log(dbinom(y[nm],total[nm],lmhu));
	  break;
	case 2: /* Poisson distribution */
	  *like-=log(dpois(y[nm],lmhu));
	  break;
	case 3: /* exponential distribution */
	  *like-=log(dexp(y[nm],lmhu));
	  break;
	case 4: /* negative binomial distribution */
	  lmhu*=lambda;
	  *like+=log(y[nm]+lmhu)+lbeta(y[nm]+1,lmhu)-lmhu*log(lambda)
	    +(y[nm]+lmhu)*log(1+lambda);
	  break;
	case 5: /* multiplicative Poisson distribution */
	  iy=y[nm];
	  dmp(&iy,my,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 6: /* double Poisson distribution */
	  iy=y[nm];
	  ddp(&iy,my,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 7: /* beta binomial distribution */
	  *like-=lbeta(y[nm]+lambda*lmhu,total[nm]-y[nm]+lambda*(1-lmhu))-lbeta(lambda*lmhu,lambda*(1-lmhu))+lchoose(total[nm],y[nm]);
	  break;
	case 8: /* multiplicative binomial distribution */
	  iy=y[nm]; iy2=total[nm];
	  dmb(&iy,&iy2,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 9: /* double binomial distribution */
	  iy=y[nm]; iy2=total[nm];
	  ddb(&iy,&iy2,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 10: /* normal distribution */
	  *like-=log(dnorm(y[nm],lmhu,lambda));
	  break;
	case 11: /* logistic distribution */
	  *like-=log(dlogis(y[nm],lmhu,lambda));
	  break;
	case 12: /* Cauchy distribution */
	  *like-=log(dcauchy(y[nm],lmhu,lambda));
	  break;
	case 13: /* Weibull distribution */
	  *like-=log(dweibull(y[nm],lambda,lmhu));
	  break;
	case 14: /* gamma distribution */
	  *like-=log(dgamma(y[nm],lambda,lmhu/lambda));
	  break;
	case 15: /* Laplace distribution */
	  *like+=fabs(y[nm]-lmhu)/lambda+log(2.*lambda);
	  break;
	case 16: /* inverse Gauss distribution */
	  *like+=(log(lambda)+pow((y[nm]-lmhu),2)/(y[nm]*lambda*lmhu*lmhu)+
		  log(6.283185*y[nm]*y[nm]*y[nm]))/2;
	  break;
	case 17: /* Pareto distribution */
	  tmp=1/(lmhu*(lambda-1.0));
	  *like-=log(lambda*tmp)-(lambda+1.0)*log(1.0+y[nm]*tmp);
	  break;
	case 18: /* Levy distribution */
	  *like-=0.5*log(lambda/(2.0*M_PI*pow(y[nm]-lmhu,3)))-
	    lambda/(2.0*(y[nm]-lmhu));
	  break;
	case 19: /* generalized gamma distribution */
	  *like-=log(lambda2)+lambda*log(lambda)+(lambda*lambda2-1)*
	    log(y[nm])-lambda*lambda2*log(lmhu)-lgamma(lambda)-lambda*
	    pow(y[nm]/lmhu,lambda2);
	  break;
	case 20: /* generalized logistic distribution */
	  tmp=(y[nm]-lmhu)/lambda;
	  *like-=log(lambda2)-tmp-log(lambda)-(lambda2+1)*log(1+exp(-tmp));
	  break;
	case 21: /* Hjorth distribution */
	  *like-=-lambda2*log(1+lambda*y[nm])/lambda-pow(y[nm]/lmhu,2)/2+
	    log(y[nm]/(lmhu*lmhu)+lambda2/(1+lambda*y[nm]));
	  break;
	case 22: /* Burr distribution */
	  tmp=y[nm]/lmhu;
	  *like-=log(lambda/lmhu)+(lambda-1)*log(tmp)-
	    (lambda2+1)*log(1+pow(tmp,lambda)/lambda2);
	  break;
	case 23: /* generalized Weibull distribution (Mudholkar et al, 1995) */
	  *like-=log(dweibull(y[nm],lambda,lmhu))+log(lambda2)+
	    (lambda2-1)*log(1-exp(-pow(y[nm]/lmhu,lambda)));
	  break;
	case 24: /* generalized extreme value distribution */
	  tmp=pow(y[nm],lambda2)/lambda2;
	  *like-=log(lambda)+lambda*(tmp-log(lmhu))-
	    pow(exp(tmp)/lmhu,lambda)+(lambda2-1)*log(y[nm])-((lambda2>0)?
	    -pow(lmhu,-lambda):log(1.0-exp(-pow(lmhu,-lambda))));
	  break;
	case 25: /* generalized inverse Gauss distribution */
	  *like-=(lambda2-1.0)*log(y[nm])-(1.0/y[nm]+y[nm]/(lmhu*lmhu))
	    /(2.0*lambda)-lambda2*log(lmhu)
	    -log(2*bessel_k(1/(lambda*lmhu),fabs(lambda2),1.0));
	  break;
	case 26: /* power exponential distribution */
	  *like+=pow(fabs(y[nm]-lmhu)/lambda,2.0*lambda2)/2.0
	    +log(lambda*pow(2.0,tmp))+lgamma(tmp);
	  break;}}
      else {
	switch(*model){
	case 3: /* exponential distribution */
	  lik=pexp(y[nm],lmhu);
	  break;
	case 10: /* normal distribution */
	  lik=pnorm(y[nm],lmhu,lambda);
	  break;
	case 11: /* logistic distribution */
	  lik=plogis(y[nm],lmhu,lambda);
	  break;
	case 12: /* Cauchy distribution */
	  lik=pcauchy(y[nm],lmhu,lambda);
	  break;
	case 13: /* Weibull distribution */
	  lik=pweibull(y[nm],lambda,lmhu);
	  break;
	case 14: /* gamma distribution */
	  lik=pgamma(y[nm],lambda,lmhu/lambda);
	  break;
	case 15: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-lmhu)/lambda)/2;
	  lik=y[nm]<lmhu?tmp:1-tmp;
	  break;
	case 16: /* inverse Gauss distribution */
	  tmp=y[nm]/lmhu;
	  tmp2=sqrt(y[nm]*lambda);
	  lik=pnorm((tmp-1)/tmp2,0,1)+exp(2/(lmhu*lambda))*
	    pnorm(-(tmp+1)/tmp2,0,1);
	  break;
	case 17: /* Pareto distribution */
	  lik=pow(1.0-(1.0+y[nm]/(lmhu*(lambda-1.0))),-lambda);
	  break;
	case 18: /* Levy distribution */
	  plevy(&y[nm],&lmhu,&lambda,&tmp,&nn,&eps,&pts,&max,&iy,&lik);
	  break;
	case 19: /* generalized gamma distribution */
	  lik=pgamma(pow(y[nm],lambda2),lambda,pow(lmhu/lambda,lambda2));
	  break;
	case 20: /* generalized logistic distribution */
	  lik=pow((1+exp(-(y[nm]-lmhu)/lambda)),-lambda2);
	  break;
	case 21: /* Hjorth distribution */
	  lik=1-pow((1+lambda*y[nm]),(-lambda2/lambda))*
	    exp(-pow((y[nm]/lmhu),2)/2);
	  break;
	case 22: /* Burr distribution */
	  lik=1-pow((1+pow(y[nm]/lmhu,lambda2)/lambda2),-lambda);
	  break;
	case 23: /* generalized Weibull distribution (Mudholkar et al, 1995) */
	  lik=1-pow((1-exp(-pow((y[nm]/lmhu),lambda))),lambda2);
	  break;
	case 24: /* generalized extreme value distribution */
	  lik=pweibull(exp(pow(y[nm],lambda2)/lambda2),lambda,lmhu)/
	    ((lambda2>0)?exp(-pow(lmhu,-lambda)):
	    (1.0-exp(-pow(lmhu,-lambda))));
	  break;
	case 25: /* generalized inverse Gauss distribution */
	  pginvgauss(&y[nm],&lmhu,&lambda,&lambda2,&nn,&eps,&pts,&max,&iy,&lik);
	  break;
	case 26: /* power exponential distribution */
	  ppowexp(&y[nm],&lmhu,&lambda,&lambda2,&nn,&eps,&pts,&max,&iy,&lik);
	  break;}
	if(censor[nm]==0)*like-=lik<1.?log(1.-lik):0;
	else *like-=lik>0.?log(lik):-35.;}
      nm++;}}}
