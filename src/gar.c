/*
 *  repeated : A Library of Repeated Measurements Models
 *  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
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
 *  void gar_c(double y[], double total[], int *my, int nobs[], int *nind,
 *	 double times[], int censor[], int *cens, double eta[], int *tm,
 *	 double theta[], double dep[], int *model, int *thp, double shape[],
 *	 int *sh, int *link, int *ar, int *order, double *parch,
 *	 int *arch, double pred[], double rpred[], double volatility[],
 *       double *like)
 *
 *  DESCRIPTION
 *
 *    Function to compute the likelihood function for generalized nonlinear
 *  autoregression models with various distributions.
 *
 */

#include <math.h>
#include <stddef.h>
#include "dist.h"
#include "R.h"
#include "Rmath.h"

void gar_c(double y[], double total[], int *my, int nobs[], int *nind,
	double times[], int censor[], int *cens, double eta[], int *tm,
	double theta[], double dep[], int *model, int *thp, double shape[],
	int *sh, int *link, int *ar, int *order, double *parch, int *arch,
	double pred[], double rpred[], double volatility[], double *like){
  int nm,ii,jj,pos,bin,nn,iy,iy2,pts,max;
  double cres,cresp,cw,cwp,lmhu,lmhu2,diff,diff2,tmp,tmp2,lik,lambda,lambda2,
    lambda3,th1,th2,th3,cres2,cresp2,wt,eps;

  nn=1;
  wt=1.0;
  pts=5;
  max=16;
  eps=1.0e-6;
  lambda3=th1=th2=th3=0.0;
  pos=(*model!=11)&&(*model!=12)&&(*model!=13)&&(*model!=16)&&(*model!=19)
    &&(*model!=24)&&(*model!=30)&&(*model!=33);
  bin=(*model==1)||(*model==8)||(*model==9)||(*model==10);
  if(*thp){
    th1=theta[0];
    th2=theta[1];}
  else if(!*tm)th2=theta[0];
  if(*order==2)th3=theta[1];
  if(*ar>2){
    th2=exp(-th2);
    if(*order==2)th3=exp(-th3);}
  if(*model>3&&!*sh)lambda=exp(theta[*thp+*order]);
  if(*model>22&&*model!=31)
    lambda2=*sh?theta[*thp+*order]:exp(theta[*thp+*order+1]);
  else if(*model==31)lambda2=theta[*thp+*order+1];
  if(*model==30)tmp=1.0+1.0/(2.0*lambda2);
  if(*arch>0&&!*sh)lambda3=lambda;
  *like=0.;
  nm=0;
  for(ii=0;ii<*nind;ii++){
    cres=cres2=cw=lmhu=lmhu2=0.;
    if(*arch>0)lambda=lambda3;
    for(jj=0;jj<nobs[ii];jj++){
      if(*sh)lambda=shape[nm];
      if(jj>0){
	diff=times[nm]-times[nm-1];
	cresp=diff>0.0?cres:0.0;
	cwp=cw;
	if(*thp){
	  cw=exp(-th1*diff);
	  cres=cw*cresp;
	  cw=cw*cwp;}
	else cres=cw=0.;
	if(*arch==1){
	  if(*sh)lambda+=exp(*parch*diff)*pow(cresp,2);
	  else lambda=lambda3+exp(*parch*diff)*pow(cresp,2);}
	else if(*arch==2){
	  if(*sh)lambda+=exp(*parch*diff)*fabs(cresp);
	  else lambda=lambda3+exp(*parch*diff)*fabs(cresp);}
	else if(*arch==3){
	  if(*sh)lambda*=exp(exp(*parch*diff)*pow(cresp,2));
	  else lambda=lambda3*exp(exp(*parch*diff)*pow(cresp,2));}
	if(*tm)th2=-dep[nm];
	if(*order==0&&!*tm)lmhu=0.0;
	else {
	  switch(*ar){
	  case 1: lmhu=exp(-th2*diff)*cresp/cwp; break;
	  case 2: lmhu=exp(-th2*diff*diff)*cresp/cwp; break;
	  case 3: lmhu=cresp/cwp/(1+th2*diff*diff); break;
	  case 4: lmhu=diff<=1/th2?(pow(diff*th2,3)-3*th2*diff+2)*
		    cresp/cwp/2:0;
	    break;
	  case 5: lmhu=(2*th2*times[nm-1]+exp(-th2*times[nm])
		    +exp(-th2*times[nm-1])-1-exp(-th2*diff))/(2*pow(th2,3))
		    *cresp/cwp; break;}}
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
      if(*arch>0||*sh)volatility[nm]=sqrt(lambda);
      cw++;
      switch(*link){
      case 1: pred[nm]=eta[nm]; break;
      case 2: pred[nm]=log(eta[nm]); break;
      case 3: pred[nm]=sqrt(eta[nm]); break;
      case 4: pred[nm]=eta[nm]*eta[nm]; break;
      case 5: pred[nm]=exp(eta[nm]); break;
      case 6: pred[nm]=1/(1+exp(-eta[nm])); break;
      case 7: pred[nm]=1-exp(-exp(eta[nm])); break;
      case 8: pred[nm]=exp(-exp(eta[nm])); break;}
      lmhu+=pred[nm]+lmhu2;
      if(pos&&lmhu<=0.0)lmhu=0.01;
      if((*model==20||*model==21||*model==22)&&lmhu>=1.)lmhu=0.99;
      if(*model==19&&lmhu>=y[nm])lmhu=y[nm]-0.01;
      if(bin){
	if(lmhu>=1.)lmhu=0.99;
	rpred[nm]=lmhu*total[nm];
	cres+=y[nm]-pred[nm]*total[nm];}
      else {
	rpred[nm]=lmhu;
	cres+=y[nm]-pred[nm];}
      if(!*cens||censor[nm]==1){
	switch(*model){
	case 1: /* binomial distribution */
	  if(total[nm]==1)*like-=(int)y[nm]?log(lmhu):log(1-lmhu);
	  else *like-=dbinom(y[nm],total[nm],lmhu,1);
	  break;
	case 2: /* Poisson distribution */
	  *like-=dpois(y[nm],lmhu,1);
	  break;
	case 3: /* exponential distribution */
	  *like-=dexp(y[nm],lmhu,1);
	  break;
	case 4: /* negative binomial distribution */
	  *like-=dnbinom(y[nm],lambda,lambda/(lambda+lmhu),1);
	  break;
	case 5: /* multiplicative Poisson distribution */
	  iy=y[nm];
	  dmp_c(&iy,my,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 6: /* double Poisson distribution */
	  iy=y[nm];
	  ddp_c(&iy,my,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 7: /* Consul generalized Poisson distribution */
	  *like-=log(lmhu)-(lmhu+y[nm]*(lambda-1))/lambda+(y[nm]-1)*
	    log(lmhu+y[nm]*(lambda-1))-y[nm]*log(lambda)-lgammafn(y[nm]+1);
	  break;
	case 8: /* beta binomial distribution */
	  *like-=lbeta(y[nm]+lambda*lmhu,total[nm]-y[nm]+lambda*
	    (1-lmhu))-lbeta(lambda*lmhu,lambda*(1-lmhu))+
	    lchoose(total[nm],y[nm]);
	  break;
	case 9: /* multiplicative binomial distribution */
	  iy=y[nm]; iy2=total[nm];
	  dmb_c(&iy,&iy2,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 10: /* double binomial distribution */
	  iy=y[nm]; iy2=total[nm];
	  ddb_c(&iy,&iy2,&lmhu,&lambda,&nn,&wt,&tmp);
	  *like-=tmp;
	  break;
	case 11: /* normal distribution */
	  *like-=dnorm(y[nm],lmhu,sqrt(lambda),1);
	  break;
	case 12: /* logistic distribution */
	  *like-=dlogis(y[nm],lmhu,lambda,1);
	  break;
	case 13: /* Cauchy distribution */
	  *like-=dcauchy(y[nm],lmhu,lambda,1);
	  break;
	case 14: /* Weibull distribution */
	  *like-=dweibull(y[nm],lambda,lmhu,1);
	  break;
	case 15: /* gamma distribution */
	  *like-=dgamma(y[nm],lambda,lmhu/lambda,1);
	  break;
	case 16: /* Laplace distribution */
	  *like+=fabs(y[nm]-lmhu)/lambda+log(2.*lambda);
	  break;
	case 17: /* inverse Gauss distribution */
	  *like+=(log(lambda)+pow((y[nm]-lmhu),2)/(y[nm]*lambda*lmhu*lmhu)+
		  log(6.283185*y[nm]*y[nm]*y[nm]))/2;
	  break;
	case 18: /* Pareto distribution */
	  tmp=1/(lmhu*(lambda-1.0));
	  *like-=log(lambda*tmp)-(lambda+1.0)*log(1.0+y[nm]*tmp);
	  break;
	case 19: /* Levy distribution */
	  *like-=0.5*log(lambda/(2.0*M_PI*pow(y[nm]-lmhu,3)))-
	    lambda/(2.0*(y[nm]-lmhu));
	  break;
	case 20: /* beta distribution */
	  *like-=dbeta(y[nm],lmhu*lambda/(1-lmhu),lambda,1);
	  break;
	case 21: /* simplex distribution */
	  *like+=(pow((y[nm]-lmhu)/(lmhu*(1-lmhu)),2)/(y[nm]*(1-y[nm])*lambda)
	    +log(2*M_PI*lambda*pow(y[nm]*(1-y[nm]),3)))/2;
	  break;
	case 22: /* two-sided power distribution */
	  *like-=log(lambda)+(y[nm]<lmhu?(lambda-1)*log(y[nm]/lmhu):
	    (lambda-1)*log((1-y[nm])/(1-lmhu)));
	  break;
	case 23: /* generalized gamma distribution */
	  *like-=log(lambda2)+lambda*log(lambda)+(lambda*lambda2-1)*
	    log(y[nm])-lambda*lambda2*log(lmhu)-lgammafn(lambda)-lambda*
	    pow(y[nm]/lmhu,lambda2);
	  break;
	case 24: /* generalized logistic distribution */
	  tmp=(y[nm]-lmhu)/lambda;
	  *like-=log(lambda2)-tmp-log(lambda)-(lambda2+1)*log(1+exp(-tmp));
	  break;
	case 25: /* Hjorth distribution */
	  *like-=-lambda2*log(1+lambda*y[nm])/lambda-pow(y[nm]/lmhu,2)/2+
	    log(y[nm]/(lmhu*lmhu)+lambda2/(1+lambda*y[nm]));
	  break;
	case 26: /* Burr distribution */
	  tmp=y[nm]/lmhu;
	  *like-=log(lambda2*lambda/lmhu)+(lambda-1)*log(tmp)-
	    (lambda2+1)*log(1+pow(tmp,lambda));
	  break;
	case 27: /* generalized Weibull distribution (Mudholkar et al, 1995) */
	  *like-=dweibull(y[nm],lambda,lmhu,1)+log(lambda2)+
	    (lambda2-1)*log(1-exp(-pow(y[nm]/lmhu,lambda)));
	  break;
	case 28: /* generalized extreme value distribution */
	  tmp=pow(y[nm],lambda2)/lambda2;
	  *like-=log(lambda)+lambda*(tmp-log(lmhu))-
	    pow(exp(tmp)/lmhu,lambda)+(lambda2-1)*log(y[nm])-((lambda2>0)?
	    -pow(lmhu,-lambda):log(1.0-exp(-pow(lmhu,-lambda))));
	  break;
	case 29: /* generalized inverse Gauss distribution */
	  *like-=(lambda2-1.0)*log(y[nm])-(1.0/y[nm]+y[nm]/(lmhu*lmhu))
	    /(2.0*lambda)-lambda2*log(lmhu)
	    -log(2*bessel_k(1/(lambda*lmhu),fabs(lambda2),1.0));
	  break;
	case 30: /* power exponential distribution */
	  *like+=pow(fabs(y[nm]-lmhu)/sqrt(lambda),2.0*lambda2)/2.0
	    +log(sqrt(lambda)*pow(2.0,tmp))+lgammafn(tmp);
	  break;
	case 31: /* power variance function Poisson distribution */
	  iy=y[nm];
	  tmp=lambda2/lmhu;
	  dpvfp_c(&iy,&lmhu,&lambda,&tmp,&nn,&wt,&lik);
	  *like-=log(lik);
	  break;
	case 32: /* skew Laplace distribution */
	  *like-=log(lambda2)+(y[nm]>lmhu?-lambda2*(y[nm]-lmhu):(y[nm]-lmhu)/lambda2)/lambda-log(1+lambda2*lambda2)-log(lambda);
	  break;
	case 33: /* Student t */
	  *like-=dt((y[nm]-lmhu)/lambda,lambda2,1)-log(lambda);
	  break;}}
      else {
	switch(*model){
	case 3: /* exponential distribution */
	  lik=pexp(y[nm],lmhu,1,0);
	  break;
	case 11: /* normal distribution */
	  lik=pnorm(y[nm],lmhu,sqrt(lambda),1,0);
	  break;
	case 12: /* logistic distribution */
	  lik=plogis(y[nm],lmhu,lambda,1,0);
	  break;
	case 13: /* Cauchy distribution */
	  lik=pcauchy(y[nm],lmhu,lambda,1,0);
	  break;
	case 14: /* Weibull distribution */
	  lik=pweibull(y[nm],lambda,lmhu,1,0);
	  break;
	case 15: /* gamma distribution */
	  lik=pgamma(y[nm],lambda,lmhu/lambda,1,0);
	  break;
	case 16: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-lmhu)/lambda)/2;
	  lik=y[nm]<lmhu?tmp:1-tmp;
	  break;
	case 17: /* inverse Gauss distribution */
	  tmp=y[nm]/lmhu;
	  tmp2=sqrt(y[nm]*lambda);
	  lik=pnorm((tmp-1)/tmp2,0,1,1,0)+exp(2/(lmhu*lambda))*
	    pnorm(-(tmp+1)/tmp2,0,1,1,0);
	  break;
	case 18: /* Pareto distribution */
	  lik=1.0-pow(1.0+y[nm]/(lmhu*(lambda-1.0)),-lambda);
	  break;
	case 19: /* Levy distribution */
	  lik=2*(1-pnorm(1/sqrt((y[nm]-lmhu)/lambda),0,1,1,0));
	  break;
	case 20: /* beta distribution */
	  lik=pbeta(y[nm],lmhu*lambda/(1-lmhu),lambda,1,0);
	  break;
	case 21: /* simplex distribution */
	  psimplex_c(&y[nm],&lmhu,&lambda,&tmp,&nn,&eps,&pts,&max,&iy,&lik);
	  break;
	case 22: /* two-sided power distribution */
	  lik=y[nm]<lmhu?log(lmhu)+lambda*log(y[nm]/lmhu):
	    log(1-(1-lmhu)*pow((1-y[nm])/(1-lmhu),lambda));
	  break;
	case 23: /* generalized gamma distribution */
	  lik=pgamma(pow(y[nm],lambda2),lambda,pow(lmhu/lambda,lambda2),1,0);
	  break;
	case 24: /* generalized logistic distribution */
	  lik=pow((1+exp(-(y[nm]-lmhu)/lambda)),-lambda2);
	  break;
	case 25: /* Hjorth distribution */
	  lik=1-pow((1+lambda*y[nm]),(-lambda2/lambda))*
	    exp(-pow((y[nm]/lmhu),2)/2);
	  break;
	case 26: /* Burr distribution */
	  lik=1-pow((1+pow(y[nm]/lmhu,lambda)),-lambda2);
	  break;
	case 27: /* generalized Weibull distribution (Mudholkar et al, 1995) */
	  lik=1-pow((1-exp(-pow((y[nm]/lmhu),lambda))),lambda2);
	  break;
	case 28: /* generalized extreme value distribution */
	  lik=pweibull(exp(pow(y[nm],lambda2)/lambda2),lambda,lmhu,1,0)/
	    ((lambda2>0)?exp(-pow(lmhu,-lambda)):
	    (1.0-exp(-pow(lmhu,-lambda))));
	  break;
	case 29: /* generalized inverse Gauss distribution */
	  pginvgauss_c(&y[nm],&lmhu,&lambda,&lambda2,&nn,&eps,&pts,&max,&iy,&lik);
	  break;
	case 30: /* power exponential distribution */
	  ppowexp_c(&y[nm],&lmhu,&lambda,&lambda2,&nn,&eps,&pts,&max,&iy,&lik);
	  lik=lik-lmhu>0?0.5+lik:0.5-lik;
	  break;
	case 32: /* skew Laplace distribution */
	  tmp=(y[nm]-lmhu)/lambda;
	  lik=tmp>0?1-exp(-lambda2*fabs(tmp))/(1+lambda2*lambda2):lambda2*lambda2*exp(-fabs(tmp)/lambda2)/(1+lambda2*lambda2);
	  break;
	case 33: /* Student t */
	  lik=pt((y[nm]-lmhu)/lambda,lambda2,1,0);
	  break;}
	if(censor[nm]==0)*like-=lik<1.?log(1.-lik):0;
	else *like-=lik>0.?log(lik):-35.;}
      nm++;}}}
