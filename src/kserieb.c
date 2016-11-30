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
 *  void kserie(double p[],double y[],double t[],double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *npv,int *model,int *link,
 *	int *density,int *pfamily,int *dep,int *torder,int inter[],int *tvc,
 *	double tvcov[],int *fit,double pred[],double rpred[],int *rf,
 *	double bb[],int *sf, double vv[], double *like)
 *  void krand(double p[],double y[],double t[],double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *npv,int *model,int *link,
 *	int *density,int *torder,int inter[],int *tvc,double tvcov[],
 *	int *fit,double pred[],double rpred[],int *rf,
 *	double bb[],int *sf, double vv[],int *frser,double *like)
 *
 *  DESCRIPTION
 *
 *    Functions to compute the likelihood function for various distributions
 * inserted in a beta distribution with serial or frailty dependence using
 * Kalman-type update for continuous longitudinal data.
 *
 */

#include <math.h>
#include <stddef.h>
#include "R.h"
#include "Rmath.h"

void kserie(double p[],double y[],double t[],double x[],int *nind,int nobs[],
	    int *nbs,int *nccov,int *npv,int *model,int *link,int *density,
	    int *pfamily,int *dep,int *torder,int inter[],int *tvc,
	    double tvcov[],int *fit,double pred[],double rpred[],int *rf,
	    double bb[],int *sf, double vv[], double *like){
  int i,j,k,k1,k2,nm;
  double a,a1,b,b1,delta,lambda,omega,om,beta,bet,bt,h,tmp,ly,
    plap,intercept,v,family;

  *like=0;
  nm=0;
  delta=exp(-p[*nccov+*npv+*tvc+1]);
  if(*dep>0)
    omega=exp(p[*nccov+*npv+*tvc+2])/(1+exp(p[*nccov+*npv+*tvc+2]));
  if(*pfamily)family=p[*nccov+*npv+*tvc+2+(*dep>0)];
  if(*model==4)intercept=exp(p[*nccov+*npv+*tvc+2+(*dep>0)+*pfamily]);
  if(*model>1&&!*sf){
    if(*model!=5&&*model!=9)
      lambda=exp(p[*nccov+*npv+*tvc+2+(*dep>0)+*pfamily+(*model==4)]);
    else lambda=exp(p[*nccov+*npv+*tvc+2+(*dep>0)+*pfamily+(*model==4)]/2);}
  for(i=0;i<*nind;i++){
    a=b=delta;
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*tvc==0&&*torder==0)switch(*link){
      case 1: bt=beta; break;
      case 2: bt=log(beta); break;
      case 3: bt=sqrt(beta); break;
      case 4: bt=beta*beta; break;
      case 5: bt=exp(beta); break;}
      else bet=beta;}
    else if(*tvc==0)bt=bb[i];
    for(j=0;j<nobs[i];j++){
      if(*model>=9)ly=log(y[nm]);
      if(*model>1&&*sf)lambda=vv[nm];
      a1=a+1;
      b1=b;
      /* add in time-varying covariates */
      if(!*rf){
	if(*torder){
	  beta=bet;
	  tmp=1;
	  k1=k2=0;
	  for(k=0;k<*npv;k++){
	    if(k<*torder)tmp*=t[nm];
	    else {
	      if(k2>inter[k1]){
		k1++;
		k2=0;}
	      if(k2==0){
		tmp=x[i+k1**nind]*t[nm];
		k2++;}
	      else {
		tmp*=t[nm];
		k2++;}}
	    beta+=p[*nccov+k+1]*tmp;}}
	if(*tvc>0){
	  beta=*torder?beta:bet;
	  for(k=0;k<*tvc;k++)beta+=p[*nccov+*npv+k+1]*tvcov[nm+*nbs*k];}
	if(*torder||*tvc){
	  switch(*link){
	  case 1: bt=beta; break;
	  case 2: bt=log(beta); break;
	  case 3: bt=sqrt(beta); break;
	  case 4: bt=beta*beta; break;
	  case 5: bt=exp(beta); break;}}}
      else if(*tvc>0)bt=bb[nm];
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1: /* exponential distribution */
	  b1+=y[nm]/bt;
	  h=-log(bt);
	  break;
	case 2: /* Weibull distribution */
	  b1+=pow(y[nm]/bt,lambda);
	  h=log(lambda/bt)+(lambda-1)*log(y[nm]/bt);
	  break;
	case 3: /* gamma distribution */
	  b1-=pgamma(y[nm],lambda,bt,0,1);
	  h=dgamma(y[nm],lambda,bt,1)-pgamma(y[nm],lambda,bt,0,1);
	  break;
	case 4: /* generalized logistic distribution */
	  b1+=(y[nm]+log(lambda+intercept*exp(-bt*y[nm]))/bt)/lambda;
	  h=-log(lambda+intercept*exp(-bt*y[nm]));
	  break;
	case 5: /* normal distribution */
	  b1-=pnorm(y[nm],bt,lambda,0,1);
	  h=dnorm(y[nm],bt,lambda,1)-pnorm(y[nm],bt,lambda,0,1);
	  break;
	case 6: /* logistic distribution */
	  b1-=plogis(y[nm],bt,lambda,0,1);
	  h=dlogis(y[nm],bt,lambda,1)-plogis(y[nm],bt,lambda,0,1);
	  break;
	case 7: /* Cauchy distribution */
	  b1-=pcauchy(y[nm],bt,lambda,0,1);
	  h=dcauchy(y[nm],bt,lambda,1)-pcauchy(y[nm],bt,lambda,0,1);
	  break;
	case 8: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  plap=y[nm]<bt?tmp:1-tmp;
	  b1-=log(1-plap);
	  h=log(tmp/(lambda*(1-plap)));
	  break;
	case 9: /* log normal distribution */
	  b1-=pnorm(ly,bt,lambda,0,1);
	  h=dnorm(ly,bt,lambda,1)-ly-pnorm(ly,bt,lambda,0,1);
	  break;
	case 10: /* log logistic distribution */
	  b1-=plogis(ly,bt,lambda,0,1);
	  h=dlogis(ly,bt,lambda,1)-ly-plogis(ly,bt,lambda,0,1);
	  break;
	case 11: /* log Cauchy distribution */
	  b1-=pcauchy(ly,bt,lambda,0,1);
	  h=dcauchy(ly,bt,lambda,1)-ly-pcauchy(ly,bt,lambda,0,1);
	  break;
	case 12: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  b1-=log(1-plap);
	  h=log(tmp/(lambda*y[nm]*(1-plap)));
	  break;
	case 13: /* inverse Gauss distribution */
	  tmp=y[nm]/bt;
	  v=sqrt(y[nm]*lambda);
	  h=1-pnorm((tmp-1)/v,0,1,1,0)-exp(2/(bt*lambda))*pnorm(-(tmp+1)/v,0,1,1,0);
	  b1-=log(h);
	  h=-(pow(y[nm]/bt-1,2)/(2*y[nm]*lambda))-log(2*PI*lambda*pow(y[nm],3))-log(h);
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1: /* exponential distribution */
	  b1+=pexp(y[nm],bt,1,0);
	  h=dexp(y[nm],bt,1);
	  break;
	case 2: /* Weibull distribution */
	  b1+=pweibull(y[nm],lambda,bt,1,0);
	  h=dweibull(y[nm],lambda,bt,1);
	  break;
	case 3: /* gamma distribution */
	  b1+=pgamma(y[nm],lambda,bt,1,0);
	  h=dgamma(y[nm],lambda,bt,1);
	  break;
	case 4: /* generalized logistic distribution */
	  b1+=exp(-y[nm]/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])),1/(lambda*bt));
	  h=-y[nm]/lambda+(1/(lambda*bt)+1)*log((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])));
	  break;
	case 5: /* normal distribution */
	  b1+=pnorm(y[nm],bt,lambda,1,0);
	  h=dnorm(y[nm],bt,lambda,1);
	  break;
	case 6: /* logistic distribution */
	  b1+=plogis(y[nm],bt,lambda,1,0);
	  h=dlogis(y[nm],bt,lambda,1);
	  break;
	case 7: /* Cauchy distribution */
	  b1+=pcauchy(y[nm],bt,lambda,1,0);
	  h=dcauchy(y[nm],bt,lambda,1);
	  break;
	case 8: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  b1+=y[nm]<bt?tmp:1-tmp;
	  h=log(tmp/lambda);
	  break;
	case 9: /* log normal distribution */
	  b1+=pnorm(ly,bt,lambda,1,0);
	  h=dnorm(ly,bt,lambda,1)-ly;
	  break;
	case 10: /* log logistic distribution */
	  b1+=plogis(ly,bt,lambda,1,0);
	  h=dlogis(ly,bt,lambda,1)-ly;
	  break;
	case 11: /* log Cauchy distribution */
	  b1+=pcauchy(ly,bt,lambda,1,0);
	  h=dcauchy(ly,bt,lambda,1)-ly;
	  break;
	case 12: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  b1+=ly<bt?tmp:1-tmp;
	  h=log(tmp/lambda/y[nm]);
	  break;
	case 13: /* inverse Gauss distribution */
	  tmp=y[nm]/bt;
	  v=sqrt(y[nm]*lambda);
	  b1+=pnorm((tmp-1)/v,0,1,1,0)+exp(2/(bt*lambda))*pnorm(-(tmp+1)/v,0,1,1,0);
	  h=-(pow(y[nm]/bt-1,2)/(y[nm]*lambda))-log(2*PI*lambda*pow(y[nm],3))/2;
	  break;}}
      /* calculate likelihood */
      *like-=h+log(a);
      if(*pfamily)*like-=(family-1.)*log(b1)-a*(pow(b1,family)-pow(b,family))/family;
      else *like-=a*log(b)-a1*log(b1);
      /* calculate fitted values */
      if(*fit){
	pred[nm]=bt;
	tmp=b/a;
	if(!*density){
	  switch(*model){
	  case 1: rpred[nm]=bt*tmp; break;
	  case 2: rpred[nm]=bt*pow(tmp,1/lambda); break;
	  case 3: rpred[nm]=qgamma(1-exp(-tmp),lambda,bt,1,0); break;
	  case 5: rpred[nm]=qnorm(1-exp(-tmp),bt,lambda,1,0); break;
	  case 6: rpred[nm]=qlogis(1-exp(-tmp),bt,lambda,1,0); break;
	  case 7: rpred[nm]=qcauchy(1-exp(-tmp),bt,lambda,1,0); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?exp(-tmp):1-exp(-tmp)));
	    break;
	  case 9: rpred[nm]=exp(qnorm(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 10: rpred[nm]=exp(qlogis(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 11: rpred[nm]=exp(qcauchy(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	else{
	  switch(*model){
	  case 1: rpred[nm]=qexp(tmp,bt,1,0); break;
	  case 2: rpred[nm]=qweibull(tmp,lambda,bt,1,0); break;
	  case 3: rpred[nm]=qgamma(tmp,lambda,bt,1,0); break;
	  case 5: rpred[nm]=qnorm(tmp,bt,lambda,1,0); break;
	  case 6: rpred[nm]=qlogis(tmp,bt,lambda,1,0); break;
	  case 7: rpred[nm]=qcauchy(tmp,bt,lambda,1,0); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?tmp:1-tmp)); break;
	  case 9: rpred[nm]=exp(qnorm(tmp,bt,lambda,1,0)); break;
	  case 10: rpred[nm]=exp(qlogis(tmp,bt,lambda,1,0)); break;
	  case 11: rpred[nm]=exp(qcauchy(tmp,bt,lambda,1,0)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
      /* update parameters */
      if(*dep){
	om=nm&&j&&t[nm]>t[nm-1]?pow(omega,t[nm]-t[nm-1]):1;
	a=om*a1+(1-om)*delta;
	if(*dep==1)b=delta+om*(b1-b);
	else if(*dep==2)b=om*b1+(1-om)*delta;}
      nm++;}}
  return;}

void krand(double p[],double y[],double t[],double x[],int *nind,int nobs[],
	   int *nbs,int *nccov,int *npv,int *model,int *link,int *density,
	   int *torder,int inter[],int *tvc,double tvcov[],
	   int *fit,double pred[],double rpred[],int *rf,
	   double bb[],int *sf, double vv[], int *frser, double *like){
  int i,j,k,k1,k2,nm,nn,pos;
  double b1,delta,lambda,beta,bet,bt,l1,tmp,ly,plap,intercept,
    H,btp,res,v;

  *like=0;
  nm=0;
  delta=exp(p[*nccov+*npv+*tvc+1]);
  pos=*model<=3||*model>=9;
  if(*model>1&&!*sf){
    if(*model!=5&&*model!=9)lambda=exp(p[*nccov+*npv+*tvc+*frser+2]);
    else lambda=exp(p[*nccov+*npv+*tvc+*frser+2]/2);}
  if(*model==4)intercept=exp(p[*nccov+*npv+*tvc+*frser+3]);
  for(nn=i=0;i<*nind;i++)nn+=nobs[i];
  for(i=0;i<*nind;i++){
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*tvc==0&&*torder==0)switch(*link){
      case 1: bt=beta; break;
      case 2: bt=log(beta); break;
      case 3: bt=sqrt(beta); break;
      case 4: bt=beta*beta; break;
      case 5: bt=exp(beta); break;}
      else bet=beta;}
    else if(!*tvc)bt=bb[i];
    res=b1=0;
    for(j=0;j<nobs[i];j++){
      l1=log(1+delta*j);
      if(*model>=9)ly=log(y[nm]);
      if(*model>1&&*sf)lambda=vv[nm];
      /* add in time-varying covariates */
      if(!*rf){
	if(*torder){
	  beta=bet;
	  tmp=1;
	  k1=k2=0;
	  for(k=0;k<*npv;k++){
	    if(k<*torder)tmp*=t[nm];
	    else {
	      if(k2>inter[k1]){
		k1++;
		k2=0;}
	      if(k2==0){
		tmp=x[i+k1**nind]*t[nm];
		k2++;}
	      else {
		tmp*=t[nm];
		k2++;}}
	    beta+=p[*nccov+k+1]*tmp;}}
	if(*tvc>0){
	  beta=bet;
	  for(k=0;k<*tvc;k++)beta+=p[*nccov+*npv+k+1]*tvcov[nm+*nbs*k];}
	if(*torder||*tvc){
	  switch(*link){
	  case 1: bt=beta; break;
	  case 2: bt=log(beta); break;
	  case 3: bt=sqrt(beta); break;
	  case 4: bt=beta*beta; break;
	  case 5: bt=exp(beta); break;}}}
      else if(*tvc>0)bt=bb[nm];
      /* if AR, add discounted previous residual */
      btp=bt;
      if(*frser&&j>0){
	bt+=exp(p[*nccov+*npv+*tvc+2]*(t[nm]-t[nm-1]))*res;
	if(pos&&bt<=0.0)bt=0.01;}
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1: /* exponential distribution */
	  H=y[nm]/bt;
	  l1+=-log(bt);
	  break;
	case 2: /* Weibull distribution */
	  H=pow(y[nm]/bt,lambda);
	  l1+=log(lambda/bt)+(lambda-1)*log(y[nm]/bt);
	  break;
	case 3: /* gamma distribution */
	  H=-pgamma(y[nm],lambda,bt,0,1);
	  l1+=dgamma(y[nm],lambda,bt,1)-pgamma(y[nm],lambda,bt,0,1);
	  break;
	case 4: /* generalized logistic distribution */
	  H=(y[nm]+log(lambda+intercept*exp(-bt*y[nm]))/bt)/lambda;
	  l1+=-log(lambda+intercept*exp(-bt*y[nm]));
	  break;
	case 5: /* normal distribution */
	  H=-pnorm(y[nm],bt,lambda,0,1);
	  l1+=dnorm(y[nm],bt,lambda,1)-pnorm(y[nm],bt,lambda,0,1);
	  break;
	case 6: /* logistic distribution */
	  H=-plogis(y[nm],bt,lambda,0,1);
	  l1+=dlogis(y[nm],bt,lambda,1)-plogis(y[nm],bt,lambda,0,1);
	  break;
	case 7: /* Cauchy distribution */
	  H=-pcauchy(y[nm],bt,lambda,0,1);
	  l1+=dcauchy(y[nm],bt,lambda,1)-pcauchy(y[nm],bt,lambda,0,1);
	  break;
	case 8: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  plap=y[nm]<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  l1+=log(tmp/(lambda*(1-plap)));
	  break;
	case 9: /* log normal distribution */
	  H=-pnorm(ly,bt,lambda,0,1);
	  l1+=dnorm(ly,bt,lambda,1)-ly-pnorm(ly,bt,lambda,0,1);
	  break;
	case 10: /* log logistic distribution */
	  H=-plogis(ly,bt,lambda,0,1);
	  l1+=dlogis(ly,bt,lambda,1)-ly-plogis(ly,bt,lambda,0,1);
	  break;
	case 11: /* log Cauchy distribution */
	  H=-pcauchy(ly,bt,lambda,0,1);
	  l1+=dcauchy(ly,bt,lambda,1)-ly-pcauchy(ly,bt,lambda,0,1);
	  break;
	case 12: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  l1+=log(tmp/(lambda*y[nm]*(1-plap)));
	  break;
	case 13: /* inverse Gauss distribution */
	  tmp=y[nm]/bt;
	  v=sqrt(y[nm]*lambda);
	  H=1-pnorm((tmp-1)/v,0,1,1,0)+exp(2/(bt*lambda))*pnorm(-(tmp+1)/v,0,1,1,0);
	  l1+=-(pow(y[nm]/bt-1,2)/(y[nm]*lambda))-log(H*sqrt(2*PI*lambda*pow(y[nm],3)));
	  H=-log(H);
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1: /* exponential distribution */
	  H=pexp(y[nm],bt,1,0);
	  l1+=dexp(y[nm],bt,1);
	  break;
	case 2: /* Weibull distribution */
	  H=pweibull(y[nm],lambda,bt,1,0);
	  l1+=dweibull(y[nm],lambda,bt,1);
	  break;
	case 3: /* gamma distribution */
	  H=pgamma(y[nm],lambda,bt,1,0);
	  l1+=dgamma(y[nm],lambda,bt,1);
	  break;
	case 4: /* generalized logistic distribution */
	  H=exp(-y[nm]/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])),1/(lambda*bt));
	  l1+=-y[nm]/lambda+log((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])))/((lambda*bt)+1);
	  break;
	case 5: /* normal distribution */
	  H=pnorm(y[nm],bt,lambda,1,0);
	  l1+=dnorm(y[nm],bt,lambda,1);
	  break;
	case 6: /* logistic distribution */
	  H=plogis(y[nm],bt,lambda,1,0);
	  l1+=dlogis(y[nm],bt,lambda,1);
	  break;
	case 7: /* Cauchy distribution */
	  H=pcauchy(y[nm],bt,lambda,1,0);
	  l1+=dcauchy(y[nm],bt,lambda,1);
	  break;
	case 8: /* Laplace distribution */
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  H=y[nm]<bt?tmp:1-tmp;
	  l1+=log(tmp/lambda);
	  break;
	case 9: /* log normal distribution */
	  H=pnorm(ly,bt,lambda,1,0);
	  l1+=dnorm(ly,bt,lambda,1)-ly;
	  break;
	case 10: /* log logistic distribution */
	  H=plogis(ly,bt,lambda,1,0);
	  l1+=dlogis(ly,bt,lambda,1)-ly;
	  break;
	case 11: /* log Cauchy distribution */
	  H=pcauchy(ly,bt,lambda,1,0);
	  l1+=dcauchy(ly,bt,lambda,1)-ly;
	  break;
	case 12: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  H=ly<bt?tmp:1-tmp;
	  l1+=log(tmp/lambda/y[nm]);
	  break;
	case 13: /* inverse Gauss distribution */
	  tmp=y[nm]/bt;
	  v=sqrt(y[nm]*lambda);
	  H=pnorm((tmp-1)/v,0,1,1,0)+exp(2/(bt*lambda))*pnorm(-(tmp+1)/v,0,1,1,0);
	  l1+=-(pow(y[nm]/bt-1,2)/(y[nm]*lambda))-log(sqrt(2*PI*lambda*pow(y[nm],3)));
	  break;}}
      /* calculate likelihood and residual */
      *like-=l1;
      if(*frser)res=(*model>=9?ly:y[nm])-btp;
      /* calculate fitted values */
      if(*fit){
	pred[nm]=btp;
	tmp=(b1+1/(nn*delta))/(1/(nn*delta)+j+1);
	if(!*density){
	  switch(*model){
	  case 1: rpred[nm]=bt*tmp; break;
	  case 2: rpred[nm]=bt*pow(tmp,1/lambda); break;
	  case 3: rpred[nm]=qgamma(1-exp(-tmp),lambda,bt,1,0); break;
	  case 5: rpred[nm]=qnorm(1-exp(-tmp),bt,lambda,1,0); break;
	  case 6: rpred[nm]=qlogis(1-exp(-tmp),bt,lambda,1,0); break;
	  case 7: rpred[nm]=qcauchy(1-exp(-tmp),bt,lambda,1,0); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?exp(-tmp):1-exp(-tmp)));
	    break;
	  case 9: rpred[nm]=exp(qnorm(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 10: rpred[nm]=exp(qlogis(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 11: rpred[nm]=exp(qcauchy(1-exp(-tmp),bt,lambda,1,0)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	else{
	  switch(*model){
	  case 1: rpred[nm]=qexp(tmp,bt,1,0); break;
	  case 2: rpred[nm]=qweibull(tmp,lambda,bt,1,0); break;
	  case 3: rpred[nm]=qgamma(tmp,lambda,bt,1,0); break;
	  case 5: rpred[nm]=qnorm(tmp,bt,lambda,1,0); break;
	  case 6: rpred[nm]=qlogis(tmp,bt,lambda,1,0); break;
	  case 7: rpred[nm]=qcauchy(tmp,bt,lambda,1,0); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?tmp:1-tmp)); break;
	  case 9: rpred[nm]=exp(qnorm(tmp,bt,lambda,1,0)); break;
	  case 10: rpred[nm]=exp(qlogis(tmp,bt,lambda,1,0)); break;
	  case 11: rpred[nm]=exp(qcauchy(tmp,bt,lambda,1,0)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
      b1+=H;
      nm++;}
    *like+=(nobs[i]+1/delta)*log(1+delta*b1);}
  return;}
