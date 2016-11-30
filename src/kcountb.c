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
 * void kcountb(double p[],double y[],double *origin,int c[],double x[],
 *	int *nind,int nobs[],int *nbs,int *nccov,int *model,
 *	int *density,int *pfamily,int *dep,int *birth,int *tvc,double tvcov[],
 *	int *fit,double pred[],double rpred[],int *rf, double bbb[],
 *	int *sf, double vv[], double *like)
 * void countfb(double p[],double y[],int c[], double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *model,int *density,
 *	int *tvc,double tvcov[],int *fit,double pred[],double rpred[],int *rf,
 *	double bbb[],int *sf,double vv[],int *frser,double *like)
 *
 *  DESCRIPTION
 *
 *    Functions to compute the likelihood function for various distributions
 * inserted in a beta distribution with serial or frailty dependence using
 * Kalman-type update for longitudinal count data.
 *
 */

#include <math.h>
#include <stddef.h>
#include "R.h"
#include "Rmath.h"

static double dpvfp2(int y, double d, double s, double s1, double f);

void kcountb(double p[],double y[],double *origin,int c[],double x[],
	int *nind,int nobs[],int *nbs,int *nccov,int *model,
	int *density,int *pfamily,int *dep,int *birth,int *tvc,double tvcov[],
	int *fit,double pred[],double rpred[],int *rf, double bbb[],
	int *sf, double vv[], double *like){
  int i,j,j0,k,nm;
  double a,a1,b,b1,bb,bb0,sc,delta,lambda,omega,om,beta,bt,H,yy,yy0,
    tmp,ly,ly0,plap,intercept,family;
  
  *like=0;
  nm=0;
  delta=exp(-p[*nccov+*birth+*tvc+1]);
  if(*dep>0)omega=exp(p[*nccov+*birth+*tvc+2])/(1+exp(p[*nccov+*birth+*tvc+2]));
  if(*pfamily)family=p[*nccov+*birth+*tvc+2+(*dep>0)];
  if(*model==4)intercept=exp(p[*nccov+*birth+*tvc+2+(*dep>0)+*pfamily]);
  if(*model>1&&!*sf){
    if(*model<5)
      lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)+*pfamily+(*model==4)]);
    else lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)+*pfamily+(*model==4)]/2);}
  for(i=0;i<*nind;i++){
    a=b=delta;
    sc=bb0=ly0=0;
    yy0=*origin;
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*model<4){
	if(beta>40) beta=40;
	if(beta<-40)beta=-40;
	beta=exp(beta);}}
    else if(!*tvc)bt=bbb[i];
    j0=origin==0||*tvc?0:-1;
    for(j=j0;j<nobs[i];j++){
      yy=*origin+(j>-1?y[nm]:0);
      if(*model>=5)ly=log(yy);
      if(*model>1&&*sf)lambda=vv[nm];
      a1=a+(j>-1?c[nm]:0);
      b1=b;
      /* add in birth and time-varying covariates */
      if(!*rf){
	if(*tvc){
	  bt=0;
	  for(k=0;k<*tvc;k++)bt+=p[*nccov+*birth+k+1]*tvcov[nm+*nbs*k];
	  if(*model<4){
	    if(bt>40) bt=40;
	    if(bt<-40)bt=-40;
	    bt=exp(bt)*beta;}
	  else bt+=beta;}
	else bt=beta;
	if(j>-1&&*birth){
	  sc+=c[nm];
	  if(*model<5)bt*=pow(sc,p[*nccov+1]);
	  else bt+=p[*nccov+1]*log(sc);}}
      else if(*tvc)bt=bbb[nm];
      if(!*tvc){
	if(!*density){
	  /* intensity models */
	  switch(*model){
	  case 1: /* exponential distribution */
	    H=yy/bt; break;
	  case 2: /* Weibull distribution */
	    H=pow(yy/bt,lambda); break;
	  case 3: /* gamma distribution */
	    H=-pgamma(yy,lambda,bt,0,1); break;
	  case 4: /* generalized logistic distribution */
	    H=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda; break;
	  case 5: /* log normal distribution */
	    H=-pnorm(ly,bt,lambda,0,1); break;
	  case 6: /* log logistic distribution */
	    H=-plogis(ly,bt,lambda,0,1); break;
	  case 7: /* log Cauchy distribution */
	    H=-pcauchy(ly,bt,lambda,0,1); break;
	  case 8: /* log Laplace distribution */
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    H=-log(1-plap);
	    break;}}
	else{
	  /* density models */
	  switch(*model){
	  case 1: /* exponential distribution */
	    H=pexp(yy,bt,1,0); break;
	  case 2: /* Weibull distribution */
	    H=pweibull(yy,lambda,bt,1,0); break;
	  case 3: /* gamma distribution */
	    H=pgamma(yy,lambda,bt,1,0); break;
	  case 4: /* generalized logistic distribution */
	    H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	    break;
	  case 5: /* log normal distribution */
	    H=pnorm(ly,bt,lambda,1,0); break;
	  case 6: /* log logistic distribution */
	    H=plogis(ly,bt,lambda,1,0); break;
	  case 7: /* log Cauchy distribution */
	    H=pcauchy(ly,bt,lambda,1,0); break;
	  case 8: /* log Laplace distribution */
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    H=ly<bt?tmp:1-tmp;
	    break;}}
	  bb=H;
	  H-=bb0;
	  bb0=bb;}
      else {
	/* if there are time-varying covariates, finesse the problem
	   by integrating over the interval, but with covariate values
	   from the end of the interval */
	if(!*density){
	  /* intensity models */
	  switch(*model){
	  case 1: /* exponential distribution */
	    H=(yy-yy0)/bt; break;
	  case 2: /* Weibull distribution */
	    H=pow(yy/bt,lambda)-pow(yy0/bt,lambda); break;
	  case 3: /* gamma distribution */
	    H=-pgamma(yy,lambda,bt,0,1)+pgamma(yy0,lambda,bt,0,1); break;
	  case 4: /* generalized logistic distribution */
	    H=(yy-yy0+(log(lambda+intercept*exp(-bt*yy))-log(lambda+intercept*exp(-bt*yy0)))/bt)/lambda;
	    break;
	  case 5: /* log normal distribution */
	    H=-pnorm(ly,bt,lambda,0,1)+(yy0>0?pnorm(ly0,bt,lambda,0,1):0);
	    break;
	  case 6: /* log logistic distribution */
	    H=-plogis(ly,bt,lambda,0,1)+(yy0>0?plogis(ly0,bt,lambda,0,1):0);
	    break;
	  case 7: /* log Cauchy distribution */
	    H=-pcauchy(ly,bt,lambda,0,1)+(yy0>0?pcauchy(ly0,bt,lambda,0,1):0);
	    break;
	  case 8: /* log Laplace distribution */
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    H=-log(1-plap);
	    if(yy0>0){
	      tmp=exp(-fabs(ly0-bt)/lambda)/2;
	      plap=ly0<bt?tmp:1-tmp;
	      H+=log(1-plap);}
	    break;}}
	else{
	  /* density models */
	  switch(*model){
	  case 1: /* exponential distribution */
	    H=pexp(yy-yy0,bt,1,0); break;
	  case 2: /* Weibull distribution */
	    H=pweibull(yy,lambda,bt,1,0)-pweibull(yy0,lambda,bt,1,0); break;
	  case 3: /* gamma distribution */
	    H=pgamma(yy,lambda,bt,1,0)-pgamma(yy0,lambda,bt,1,0); break;
	  case 4: /* generalized logistic distribution */
	    H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt))-exp(-yy0/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy0)),1/(lambda*bt));
	    break;
	  case 5: /* log normal distribution */
	    H=pnorm(ly,bt,lambda,1,0)-(yy0>0?pnorm(ly0,bt,lambda,1,0):0); break;
	  case 6: /* log logistic distribution */
	    H=plogis(ly,bt,lambda,1,0)-(yy0>0?plogis(ly0,bt,lambda,1,0):0);
	    break;
	  case 7: /* log Cauchy distribution */
	    H=pcauchy(ly,bt,lambda,1,0)-(yy0>0?pcauchy(ly0,bt,lambda,1,0):0);
	    break;
	  case 8: /* log Laplace distribution */
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    H=ly<bt?tmp:1-tmp;
	    if(yy0>0)H-=exp(-fabs(ly0-bt)/lambda)/2;
	    break;}}
	ly0=ly;}
      if(j>-1){
      b1+=H;
      /* calculate likelihood */
      if(*pfamily)*like-=dpvfp2(c[nm],a,b,b1,family);
      else *like-=lgammafn(a1)-lgammafn(a)+a*log(b)-a1*log(b1);
      if(c[nm]>0)*like-=c[nm]*log(H);
      if(c[nm]>1)*like+=lgammafn(c[nm]+1);
      /* calculate fitted values */
      if(*fit){
	if(!*density){
	  switch(*model){
	  case 1: pred[nm]=1/bt; break;
	  case 2: pred[nm]=lambda*pow(yy/bt,lambda-1)/bt; break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt,0)/(1-pgamma(yy,lambda,bt,1,0)); break;
	  case 4: pred[nm]=1/(lambda+intercept*exp(-bt*yy)); break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda,0)/(y[nm]*(1-pnorm(ly,bt,lambda,1,0)));
	    break;
	  case 6:
	    pred[nm]=dlogis(ly,bt,lambda,0)/(y[nm]*(1-plogis(ly,bt,lambda,1,0)));
	    break;
	  case 7:
	    pred[nm]=dcauchy(ly,bt,lambda,0)/(y[nm]*(1-pcauchy(ly,bt,lambda,1,0)));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    pred[nm]=tmp/(lambda*y[nm]*(1-plap));
	    break;}}
	else{
	  switch(*model){
	  case 1: pred[nm]=dexp(yy,bt,0); break;
	  case 2: pred[nm]=dweibull(yy,lambda,bt,0); break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt,0); break;
	  case 4:
	    pred[nm]=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt)+1);
	    break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda,0)/y[nm]; break;
	  case 6: pred[nm]=dlogis(ly,bt,lambda,0)/y[nm]; break;
	  case 7: pred[nm]=dcauchy(ly,bt,lambda,0)/y[nm]; break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    pred[nm]=tmp/(lambda*y[nm]);
	    break;}}
	rpred[nm]=j==0?pred[nm]:a*H/b;}
      /* update parameters */
      switch(*dep){
      case 1: a=omega*a1+(1-omega)*delta; break;
      case 2:
      case 6: 
	om=pow(omega,yy-yy0);
	a=om*a1+(1-om)*delta;
	break;
      case 3: a=a1; break;
      case 4:
      case 5: a=omega*a1; break;
      default: break;}
      switch(*dep){
      case 1: b=omega*b1+(1-omega)*delta; break;
      case 2: b=om*b1+(1-om)*delta; break;
      case 3:
      case 5: b=omega*b1; break;
      case 6: b=om*(b1-b)+delta; break;
      default: break;}
    nm++;}
    yy0=yy;}}
  return;}

void countfb(double p[],double y[],int c[], double x[],int *nind,
	int nobs[],int *nbs,int *nccov,int *model,int *density,
	int *tvc,double tvcov[],int *fit,double pred[],double rpred[],int *rf,
	double bbb[], int *sf, double vv[], int *frser, double *like){
  int i,j,k,nm;
  double a,b,bb,bb0,delta,lambda,beta,bt,H,yy,tmp,ly,
    plap,intercept,res;
  
  *like=0;
  nm=0;
  delta=exp(p[*nccov+*tvc+1]);
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*tvc+*frser+2]);
    else lambda=exp(p[*nccov+*tvc+*frser+2]/2);}
  if(*model==4)intercept=exp(p[*nccov+*tvc+*frser+3]);
  for(i=0;i<*nind;i++){
    a=delta;
    b=bb0=0;
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*model<4){
	if(beta>40) beta=40;
	if(beta<-40)beta=-40;
	beta=exp(beta);}}
    else if(!*tvc)bt=bbb[i];
    res=0;
    for(j=0;j<nobs[i];j++){
      yy=y[nm];
      if(*model>=5)ly=log(yy);
      if(*model>1&&*sf)lambda=vv[nm];
      a+=c[nm];
      /* add in time-varying covariates */
      if(!*rf){
	if(*tvc){
	  bt=0;
	  for(k=0;k<*tvc;k++)bt+=p[*nccov+k+1]*tvcov[nm+*nbs*k];
	  if(*model<4)bt=exp(bt)*beta;
	  else bt+=beta;}
	else bt=beta;}
      else if(*tvc)bt=bbb[nm];
      /* if AR, add discounted previous residual */
      if(*frser&&j>0){
	if(*model<4)bt=exp(log(bt)+exp(p[*nccov+*tvc+2]*(y[nm]-y[nm-1]))*res);
	else bt+=exp(p[*nccov+*tvc+2]*(y[nm]-y[nm-1]))*res;}
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1: /* exponential distribution */
	  H=yy/bt; break;
	case 2: /* Weibull distribution */
	  H=pow(yy/bt,lambda); break;
	case 3: /* gamma distribution */
	  H=-pgamma(yy,lambda,bt,0,1); break;
	case 4: /* generalized logistic distribution */
	  H=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda; break;
	case 5: /* log normal distribution */
	  H=-pnorm(ly,bt,lambda,0,1); break;
	case 6: /* log logistic distribution */
	  H=-plogis(ly,bt,lambda,0,1); break;
	case 7: /* log Cauchy distribution */
	  H=-pcauchy(ly,bt,lambda,0,1); break;
	case 8: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1: /* exponential distribution */
	  H=pexp(yy,bt,1,0); break;
	case 2: /* Weibull distribution */
	  H=pweibull(yy,lambda,bt,1,0); break;
	case 3: /* gamma distribution */
	  H=pgamma(yy,lambda,bt,1,0); break;
	case 4: /* generalized logistic distribution */
	  H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	  break;
	case 5: /* log normal distribution */
	  H=pnorm(ly,bt,lambda,1,0); break;
	case 6: /* log logistic distribution */
	  H=plogis(ly,bt,lambda,1,0); break;
	case 7: /* log Cauchy distribution */
	  H=pcauchy(ly,bt,lambda,1,0); break;
	case 8: /* log Laplace distribution */
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  H=ly<bt?tmp:1-tmp;
	  break;}}
      bb=H;
      H-=bb0;
      bb0=bb;
      b+=H;
      /* if(*tvc)b+=H;
      else if(j==nobs[i]-1)b=bb; */
      /* calculate likelihood */
      if(c[nm]>0)*like-=c[nm]*log(H)-lgammafn(c[nm]+1);
      /* calculate fitted values */
      if(*fit||*frser){
	if(!*density){
	  switch(*model){
	  case 1: pred[nm]=1/bt; break;
	  case 2: pred[nm]=lambda*pow(yy/bt,lambda-1)/bt; break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt,0)/(1-pgamma(yy,lambda,bt,1,0)); break;
	  case 4: pred[nm]=1/(lambda+intercept*exp(-bt*yy)); break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda,0)/(y[nm]*(1-pnorm(ly,bt,lambda,1,0)));
	    break;
	  case 6:
	    pred[nm]=dlogis(ly,bt,lambda,0)/(y[nm]*(1-plogis(ly,bt,lambda,1,0)));
	    break;
	  case 7:
	    pred[nm]=dcauchy(ly,bt,lambda,0)/(y[nm]*(1-pcauchy(ly,bt,lambda,1,0)));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    pred[nm]=tmp/(lambda*y[nm]*(1-plap));
	    break;}}
	else{
	  switch(*model){
	  case 1: pred[nm]=dexp(yy,bt,0); break;
	  case 2: pred[nm]=dweibull(yy,lambda,bt,0); break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt,0); break;
	  case 4:
	    pred[nm]=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt)+1);
	    break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda,0)/y[nm]; break;
	  case 6: pred[nm]=dlogis(ly,bt,lambda,0)/y[nm]; break;
	  case 7: pred[nm]=dcauchy(ly,bt,lambda,0)/y[nm]; break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    pred[nm]=tmp/(lambda*y[nm]);
	    break;}}
	rpred[nm]=a*H/(delta+b);
	if(*frser)res=(*model>=4?ly:y[nm])-pred[nm];}
      nm++;}
    *like-=lgammafn(a)-a*log(delta+b);}
  *like-=*nind*(-lgammafn(delta)+delta*log(delta));
  return;}

/* power variance function Poisson */

static double pvfc2(int y, double d, double s1, double f){
  int i,j;
  double r,*c,tmp1,tmp2,tmp3;
  c=(double*)R_alloc((size_t)(y*y),sizeof(double));
  tmp1=gammafn(1.-f);
  tmp2=log(d);
  tmp3=log(s1);
  for(i=0;i<y;i++){
    c[i*(y+1)]=1;
    if(i>0){
      c[i*y]=gammafn(i+1-f)/tmp1;
      if(i>1)
	for(j=1;j<i;j++)c[y*i+j]=c[y*(i-1)+j-1]+c[y*(i-1)+j]*(i-(j+1)*f);}}
  r=0.;
  for(i=1;i<=y;i++){
    r+=c[y*(y-1)+i-1]*exp(i*tmp2+(i*f-y)*tmp3);}
  return(log(r));}

static double dpvfp2(int y, double d, double s, double s1, double f){
  double res;
  if(f==0.)res=dnbinom(y,d,s/s1,1);
  else {
    res=-d*(pow(s1,f)-pow(s,f))/f;
    if(y>0)res+=pvfc2(y,d,s1,f);}
  return(res);}
