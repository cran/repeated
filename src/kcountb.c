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
 * void kcountb(double p[],double y[],double *origin,int c[],double x[],
 *	int *nind,int nobs[],int *nbs,int *nccov,int *model,
 *	int *density,int *dep,int *birth,int *tvc,double tvcov[],
 *	int *fit,double pred[],double rpred[],int *rf, double bbb[],
 *	int *sf, double vv[], double *like)
 * void countfb(double p[],double y[],int c[], double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *model,int *density,
 *	int *tvc,double tvcov[],int *fit,double pred[],double rpred[],int *rf,
 *	double bbb[], int *sf, double vv[], double *like)
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

extern double lgamma(double x);
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
extern double ihgamma(double x, double shape, double scale);
extern double ihlogis(double x, double location, double scale);

void kcountb(double p[],double y[],double *origin,int c[],double x[],
	int *nind,int nobs[],int *nbs,int *nccov,int *model,
	int *density,int *dep,int *birth,int *tvc,double tvcov[],
	int *fit,double pred[],double rpred[],int *rf, double bbb[],
	int *sf, double vv[], double *like){
  int i,j,j0,k,nm;
  double a,a1,b,b1,bb,bb0,sc,delta,lambda,omega,om,beta,bt,H,yy,yy0,
    kk,tmp,ly,ly0,plap,intercept;
  
  *like=0;
  nm=0;
  delta=exp(-p[*nccov+*birth+*tvc+1]);
  if(*dep>0)omega=exp(p[*nccov+*birth+*tvc+2])/(1+exp(p[*nccov+*birth+*tvc+2]));
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)]);
    else lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)]/2);}
  if(*model==4)intercept=exp(p[*nccov+*birth+*tvc+3+(*dep>0)]);
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
	  case 1: H=yy/bt; break;
	  case 2: H=pow(yy/bt,lambda); break;
	  case 3: H=ihgamma(yy,lambda,bt); break;
	  case 4: H=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda; break;
	  case 5: H=-log(1-pnorm(ly,bt,lambda)); break;
	  case 6: H=ihlogis(ly,bt,lambda); break;
	  case 7: H=-log(1-pcauchy(ly,bt,lambda)); break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    H=-log(1-plap);
	    break;}}
	else{
	  /* density models */
	  switch(*model){
	  case 1: H=pexp(yy,bt); break;
	  case 2: H=pweibull(yy,lambda,bt); break;
	  case 3: H=pgamma(yy,lambda,bt); break;
	  case 4:
	    H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	    break;
	  case 5: H=pnorm(ly,bt,lambda); break;
	  case 6: H=plogis(ly,bt,lambda); break; H=pcauchy(ly,bt,lambda);
	    break;
	  case 8:
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
	  case 1: H=(yy-yy0)/bt; break;
	  case 2: H=pow(yy/bt,lambda)-pow(yy0/bt,lambda); break;
	  case 3: H=ihgamma(yy,lambda,bt)-ihgamma(yy0,lambda,bt); break;
	  case 4:
	    H=(yy-yy0+(log(lambda+intercept*exp(-bt*yy))-log(lambda+intercept*exp(-bt*yy0)))/bt)/lambda;
	    break;
	  case 5:
	    H=-log(1-pnorm(ly,bt,lambda))+(yy0>0?log(1-pnorm(ly0,bt,lambda)):0);
	    break;
	  case 6:
	    H=ihlogis(ly,bt,lambda)-(yy0>0?ihlogis(ly0,bt,lambda):0);
	    break;
	  case 7:
	    H=-log(1-pcauchy(ly,bt,lambda))+(yy0>0?log(1-pcauchy(ly0,bt,lambda)):0);
	    break;
	  case 8:
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
	  case 1: H=pexp(yy-yy0,bt); break;
	  case 2: H=pweibull(yy,lambda,bt)-pweibull(yy0,lambda,bt); break;
	  case 3: H=pgamma(yy,lambda,bt)-pgamma(yy0,lambda,bt); break;
	  case 4:
	    H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt))-exp(-yy0/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy0)),1/(lambda*bt));
	    break;
	  case 5: H=pnorm(ly,bt,lambda)-(yy0>0?pnorm(ly0,bt,lambda):0); break;
	  case 6:
	    H=plogis(ly,bt,lambda)-(yy0>0?plogis(ly0,bt,lambda):0);
	    break;
	  case 7:
	    H=pcauchy(ly,bt,lambda)-(yy0>0?pcauchy(ly0,bt,lambda):0);
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    H=ly<bt?tmp:1-tmp;
	    if(yy0>0)H-=exp(-fabs(ly0-bt)/lambda)/2;
	    break;}}
	ly0=ly;}
      if(j>-1){
      b1+=H;
      /* calculate likelihood */
      *like-=lgamma(a1)-lgamma(a)+a*log(b)-a1*log(b1);
      if(c[nm]>0){
	*like-=c[nm]*log(H);
	if(c[nm]>1)*like+=lgamma(c[nm]+1);}
      /* calculate fitted values */
      if(*fit){
	if(!*density){
	  switch(*model){
	  case 1: pred[nm]=1/bt; break;
	  case 2: pred[nm]=lambda*pow(yy/bt,lambda-1)/bt; break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt)/(1-dgamma(yy,lambda,bt)); break;
	  case 4: pred[nm]=1/(lambda+intercept*exp(-bt*yy)); break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda)/(y[nm]*(1-pnorm(ly,bt,lambda)));
	    break;
	  case 6:
	    pred[nm]=dlogis(ly,bt,lambda)/(y[nm]*(1-plogis(ly,bt,lambda)));
	    break;
	  case 7:
	    pred[nm]=dcauchy(ly,bt,lambda)/(y[nm]*(1-pcauchy(ly,bt,lambda)));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    pred[nm]=tmp/(lambda*y[nm]*(1-plap));
	    break;}}
	else{
	  switch(*model){
	  case 1: pred[nm]=dexp(yy,bt); break;
	  case 2: pred[nm]=dweibull(yy,lambda,bt); break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt); break;
	  case 4:
	    pred[nm]=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt)+1);
	    break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda)/y[nm]; break;
	  case 6: pred[nm]=dlogis(ly,bt,lambda)/y[nm]; break;
	  case 7: pred[nm]=dcauchy(ly,bt,lambda)/y[nm]; break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    pred[nm]=tmp/(lambda*y[nm]);
	    break;}}
	rpred[nm]=a*H/b;}
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
      default:}
      switch(*dep){
      case 1: b=omega*b1+(1-omega)*delta; break;
      case 2: b=om*b1+(1-om)*delta; break;
      case 3:
      case 5: b=omega*b1; break;
      case 6: b=om*(b1-b)+delta; break;
      default:}
    nm++;}
    yy0=yy;}}
  return;}

void countfb(double p[],double y[],int c[], double x[],int *nind,
	int nobs[],int *nbs,int *nccov,int *model,int *density,
	int *tvc,double tvcov[],int *fit,double pred[],double rpred[],int *rf,
	double bbb[], int *sf, double vv[], double *like){
  int i,j,k,nm;
  double a,b,bb,bb0,delta,lambda,omega,om,beta,bt,H,yy,kk,tmp,ly,
    plap,intercept;
  
  *like=0;
  nm=0;
  delta=exp(p[*nccov+*tvc+1]);
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*tvc+2]);
    else lambda=exp(p[*nccov+*tvc+2]/2);}
  if(*model==4)intercept=exp(p[*nccov+*tvc+3]);
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
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1: H=yy/bt; break;
	case 2: H=pow(yy/bt,lambda); break;
	case 3: H=ihgamma(yy,lambda,bt); break;
	case 4: H=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda; break;
	case 5: H=-log(1-pnorm(ly,bt,lambda)); break;
	case 6: H=ihlogis(ly,bt,lambda); break;
	case 7: H=-log(1-pcauchy(ly,bt,lambda)); break;
	case 8:
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1: H=pexp(yy,bt); break;
	case 2: H=pweibull(yy,lambda,bt); break;
	case 3: H=pgamma(yy,lambda,bt); break;
	case 4:
	  H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	  break;
	case 5: H=pnorm(ly,bt,lambda); break;
	case 6: H=plogis(ly,bt,lambda); break;
	case 7: H=pcauchy(ly,bt,lambda); break;
	case 8:
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
      if(c[nm]>0)*like-=c[nm]*log(H)-lgamma(c[nm]+1);
      /* calculate fitted values */
      if(*fit){
	if(!*density){
	  switch(*model){
	  case 1: pred[nm]=1/bt; break;
	  case 2: pred[nm]=lambda*pow(yy/bt,lambda-1)/bt; break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt)/(1-dgamma(yy,lambda,bt)); break;
	  case 4: pred[nm]=1/(lambda+intercept*exp(-bt*yy)); break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda)/(y[nm]*(1-pnorm(ly,bt,lambda)));
	    break;
	  case 6:
	    pred[nm]=dlogis(ly,bt,lambda)/(y[nm]*(1-plogis(ly,bt,lambda)));
	    break;
	  case 7:
	    pred[nm]=dcauchy(ly,bt,lambda)/(y[nm]*(1-pcauchy(ly,bt,lambda)));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    pred[nm]=tmp/(lambda*y[nm]*(1-plap));
	    break;}}
	else{
	  switch(*model){
	  case 1: pred[nm]=dexp(yy,bt); break;
	  case 2: pred[nm]=dweibull(yy,lambda,bt); break;
	  case 3: pred[nm]=dgamma(yy,lambda,bt); break;
	  case 4:
	    pred[nm]=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt)+1);
	    break;
	  case 5: pred[nm]=dnorm(ly,bt,lambda)/y[nm]; break;
	  case 6: pred[nm]=dlogis(ly,bt,lambda)/y[nm]; break;
	  case 7: pred[nm]=dcauchy(ly,bt,lambda)/y[nm]; break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    pred[nm]=tmp/(lambda*y[nm]);
	    break;}}
	rpred[nm]=a*H/(delta+b);}
      nm++;}
    *like-=lgamma(a)-a*log(delta+b);}
  *like-=*nind*(-lgamma(delta)+delta*log(delta));
  return;}
