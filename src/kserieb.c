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
 *  void kserie(double p[],double y[],double t[],double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *npv,int *model,int *link,
 *	int *density,int *dep,int *torder,int inter[],int *tvc,
 *	double tvcov[],int *fit,double pred[],double rpred[],int *rf,
 *	double bb[],int *sf, double vv[], double *like)
 *  void krand(double p[],double y[],double t[],double x[],int *nind,
 *	int nobs[],int *nbs,int *nccov,int *npv,int *model,int *link,
 *	int *density,int *torder,int inter[],int *tvc,double tvcov[],
 *	int *fit,double pred[],double rpred[],int *rf,
 *	double bb[],int *sf, double vv[], double *like)
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

extern double lgamma(double x);
extern double dexp(double x, double scale);
extern double pexp(double x, double scale);
extern double qexp(double x, double scale);
extern double dweibull(double x, double shape, double scale);
extern double pweibull(double x, double shape, double scale);
extern double qweibull(double x, double shape, double scale);
extern double dgamma(double x, double shape, double scale);
extern double pgamma(double x, double shape, double scale);
extern double qgamma(double p, double shape, double scale);
extern double dnorm(double x, double mean, double sd);
extern double pnorm(double x, double mean, double sd);
extern double qnorm(double p, double mu, double sigma);
extern double dlogis(double x, double location, double scale);
extern double plogis(double x, double location, double scale);
extern double qlogis(double x, double location, double scale);
extern double dcauchy(double x, double location, double scale);
extern double pcauchy(double x, double location, double scale);
extern double qcauchy(double x, double location, double scale);
extern double ihgamma(double x, double shape, double scale);
extern double ihlogis(double x, double location, double scale);

void kserie(double p[],double y[],double t[],double x[],int *nind,int nobs[],
	    int *nbs,int *nccov,int *npv,int *model,int *link,int *density,
	    int *dep,int *torder,int inter[],int *tvc,double tvcov[],
	    int *fit,double pred[],double rpred[],int *rf,
	    double bb[],int *sf, double vv[], double *like){
  int i,j,k,k1,k2,nm;
  double a,a1,b,b1,delta,lambda,omega,om,beta,bet,bt,h,kk,c,tmp,ly,
    plap,intercept;

  *like=0;
  nm=0;
  delta=exp(-p[*nccov+*npv+*tvc+1]);
  if(*dep>0){
    omega=exp(p[*nccov+*npv+*tvc+2])/(1+exp(p[*nccov+*npv+*tvc+2]));}
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*npv+*tvc+2+(*dep>0)]);
    else lambda=exp(p[*nccov+*npv+*tvc+2+(*dep>0)]/2);}
  if(*model==4)intercept=exp(p[*nccov+*npv+*tvc+3+(*dep>0)]);
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
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1:
	  b1+=y[nm]/bt;
	  h=1/bt;
	  break;
	case 2:
	  b1+=pow(y[nm]/bt,lambda);
	  h=lambda*pow(y[nm]/bt,lambda-1)/bt;
	  break;
	case 3:
	  b1+=ihgamma(y[nm],lambda,bt);
	  h=dgamma(y[nm],lambda,bt)/(1-pgamma(y[nm],lambda,bt));
	  break;
	case 4:
	  b1+=(y[nm]+log(lambda+intercept*exp(-bt*y[nm]))/bt)/lambda;
	  h=1/(lambda+intercept*exp(-bt*y[nm]));
	  break;
	case 5:
	  b1-=log(1-pnorm(y[nm],bt,lambda));
	  h=dnorm(y[nm],bt,lambda)/(1-pnorm(y[nm],bt,lambda));
	  break;
	case 6:
	  b1+=ihlogis(y[nm],bt,lambda);
	  h=dlogis(y[nm],bt,lambda)/(1-plogis(y[nm],bt,lambda));
	  break;
	case 7:
	  b1-=log(1-pcauchy(y[nm],bt,lambda));
	  h=dcauchy(y[nm],bt,lambda)/(1-pcauchy(y[nm],bt,lambda));
	  break;
	case 8:
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  plap=y[nm]<bt?tmp:1-tmp;
	  b1-=log(1-plap);
	  h=tmp/(lambda*(1-plap));
	  break;
	case 9:
	  b1-=log(1-pnorm(ly,bt,lambda));
	  h=dnorm(ly,bt,lambda)/y[nm]/(1-pnorm(ly,bt,lambda));
	  break;
	case 10:
	  b1+=ihlogis(ly,bt,lambda);
	  h=dlogis(ly,bt,lambda)/y[nm]/(1-plogis(ly,bt,lambda));
	  break;
	case 11:
	  b1-=log(1-pcauchy(ly,bt,lambda));
	  h=dcauchy(ly,bt,lambda)/y[nm]/(1-pcauchy(ly,bt,lambda));
	  break;
	case 12:
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  b1-=log(1-plap);
	  h=tmp/(lambda*y[nm]*(1-plap));
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1:
	  b1+=pexp(y[nm],bt);
	  h=dexp(y[nm],bt);
	  break;
	case 2:
	  b1+=pweibull(y[nm],lambda,bt);
	  h=dweibull(y[nm],lambda,bt);
	  break;
	case 3:
	  b1+=pgamma(y[nm],lambda,bt);
	  h=dgamma(y[nm],lambda,bt);
	  break;
	case 4:
	  b1+=exp(-y[nm]/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])),1/(lambda*bt));
	  h=exp(-y[nm]/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])),1/(lambda*bt)+1);
	  break;
	case 5:
	  b1+=pnorm(y[nm],bt,lambda);
	  h=dnorm(y[nm],bt,lambda);
	  break;
	case 6:
	  b1+=plogis(y[nm],bt,lambda);
	  h=dlogis(y[nm],bt,lambda);
	  break;
	case 7:
	  b1+=pcauchy(y[nm],bt,lambda);
	  h=dcauchy(y[nm],bt,lambda);
	  break;
	case 8:
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  b1+=y[nm]<bt?tmp:1-tmp;
	  h=tmp/lambda;
	  break;
	case 9:
	  b1+=pnorm(ly,bt,lambda);
	  h=dnorm(ly,bt,lambda)/y[nm];
	  break;
	case 10:
	  b1+=plogis(ly,bt,lambda);
	  h=dlogis(ly,bt,lambda)/y[nm];
	  break;
	case 11:
	  b1+=pcauchy(ly,bt,lambda);
	  h=dcauchy(ly,bt,lambda)/y[nm];
	  break;
	case 12:
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  b1+=ly<bt?tmp:1-tmp;
	  h=tmp/lambda/y[nm];
	  break;}}
      /* calculate likelihood */
      *like-=log(h)+log(a)+a*log(b)-a1*log(b1);
      /* calculate fitted values */
      if(*fit){
	pred[nm]=bt;
	tmp=b/a;
	if(!*density){
	  switch(*model){
	  case 1: rpred[nm]=bt*tmp; break;
	  case 2: rpred[nm]=bt*pow(tmp,1/lambda); break;
	  case 3: rpred[nm]=qgamma(1-exp(-tmp),lambda,bt); break;
	  case 5: rpred[nm]=qnorm(1-exp(-tmp),bt,lambda); break;
	  case 6: rpred[nm]=qlogis(1-exp(-tmp),bt,lambda); break;
	  case 7: rpred[nm]=qcauchy(1-exp(-tmp),bt,lambda); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?exp(-tmp):1-exp(-tmp)));
	    break;
	  case 9: rpred[nm]=exp(qnorm(1-exp(-tmp),bt,lambda)); break;
	  case 10: rpred[nm]=exp(qlogis(1-exp(-tmp),bt,lambda)); break;
	  case 11: rpred[nm]=exp(qcauchy(1-exp(-tmp),bt,lambda)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	else{
	  switch(*model){
	  case 1: rpred[nm]=qexp(tmp,bt); break;
	  case 2: rpred[nm]=qweibull(tmp,lambda,bt); break;
	  case 3: rpred[nm]=qgamma(tmp,lambda,bt); break;
	  case 5: rpred[nm]=qnorm(tmp,bt,lambda); break;
	  case 6: rpred[nm]=qlogis(tmp,bt,lambda); break;
	  case 7: rpred[nm]=qcauchy(tmp,bt,lambda); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?tmp:1-tmp)); break;
	  case 9: rpred[nm]=exp(qnorm(tmp,bt,lambda)); break;
	  case 10: rpred[nm]=exp(qlogis(tmp,bt,lambda)); break;
	  case 11: rpred[nm]=exp(qcauchy(tmp,bt,lambda)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
      /* update parameters */
      if(*dep){
	om=j?pow(omega,t[nm]-t[nm-1]):1;
	a=om*a1+(1-om)*delta;
	if(*dep==1)b=delta+om*(b1-b);
	else if(*dep==2)b=om*b1+(1-om)*delta;}
      nm++;}}
  return;}

void krand(double p[],double y[],double t[],double x[],int *nind,int nobs[],
	   int *nbs,int *nccov,int *npv,int *model,int *link,int *density,
	   int *torder,int inter[],int *tvc,double tvcov[],
	   int *fit,double pred[],double rpred[],int *rf,
	   double bb[],int *sf, double vv[], double *like){
  int i,j,k,k1,k2,nm,nn;
  double b1,delta,lambda,beta,bet,bt,l1,kk,tmp,ly,plap,intercept,H;

  *like=0;
  nm=0;
  delta=exp(p[*nccov+*npv+*tvc+1]);
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*npv+*tvc+2]);
    else lambda=exp(p[*nccov+*npv+*tvc+2]/2);}
  if(*model==4)intercept=exp(p[*nccov+*npv+*tvc+3]);
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
    b1=0;
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
      if(!*density){
	/* intensity models */
	switch(*model){
	case 1:
	  H=y[nm]/bt;
	  l1+=-log(bt);
	  break;
	case 2:
	  H=pow(y[nm]/bt,lambda);
	  l1+=log(lambda/bt)+(lambda-1)*log(y[nm]/bt);
	  break;
	case 3:
	  H=ihgamma(y[nm],lambda,bt);
	  l1+=log(dgamma(y[nm],lambda,bt)/(1-pgamma(y[nm],lambda,bt)));
	  break;
	case 4:
	  H=(y[nm]+log(lambda+intercept*exp(-bt*y[nm]))/bt)/lambda;
	  l1+=-log(lambda+intercept*exp(-bt*y[nm]));
	  break;
	case 5:
	  H=-log(1-pnorm(y[nm],bt,lambda));
	  l1+=log(dnorm(y[nm],bt,lambda)/(1-pnorm(y[nm],bt,lambda)));
	  break;
	case 6:
	  H=ihlogis(y[nm],bt,lambda);
	  l1+=log(dlogis(y[nm],bt,lambda)/(1-plogis(y[nm],bt,lambda)));
	  break;
	case 7:
	  H=-log(1-pcauchy(y[nm],bt,lambda));
	  l1+=log(dcauchy(y[nm],bt,lambda)/(1-pcauchy(y[nm],bt,lambda)));
	  break;
	case 8:
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  plap=y[nm]<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  l1+=log(tmp/(lambda*(1-plap)));
	  break;
	case 9:
	  H=-log(1-pnorm(ly,bt,lambda));
	  l1+=log(dnorm(ly,bt,lambda)/y[nm]/(1-pnorm(ly,bt,lambda)));
	  break;
	case 10:
	  H=ihlogis(ly,bt,lambda);
	  l1+=log(dlogis(ly,bt,lambda)/y[nm]/(1-plogis(ly,bt,lambda)));
	  break;
	case 11:
	  H=-log(1-pcauchy(ly,bt,lambda));
	  l1+=log(dcauchy(ly,bt,lambda)/y[nm]/(1-pcauchy(ly,bt,lambda)));
	  break;
	case 12:
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  plap=ly<bt?tmp:1-tmp;
	  H=-log(1-plap);
	  l1+=log(tmp/(lambda*y[nm]*(1-plap)));
	  break;}}
      else{
	/* density models */
	switch(*model){
	case 1:
	  H=pexp(y[nm],bt);
	  l1+=log(dexp(y[nm],bt));
	  break;
	case 2:
	  H=pweibull(y[nm],lambda,bt);
	  l1+=log(dweibull(y[nm],lambda,bt));
	  break;
	case 3:
	  H=pgamma(y[nm],lambda,bt);
	  l1+=log(dgamma(y[nm],lambda,bt));
	  break;
	case 4:
	  H=exp(-y[nm]/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])),1/(lambda*bt));
	  l1+=-y[nm]/lambda+log((lambda+intercept)/(lambda+intercept*exp(-bt*y[nm])))/((lambda*bt)+1);
	  break;
	case 5:
	  H=pnorm(y[nm],bt,lambda);
	  l1+=log(dnorm(y[nm],bt,lambda));
	  break;
	case 6:
	  H=plogis(y[nm],bt,lambda);
	  l1+=log(dlogis(y[nm],bt,lambda));
	  break;
	case 7:
	  H=pcauchy(y[nm],bt,lambda);
	  l1+=log(dcauchy(y[nm],bt,lambda));
	  break;
	case 8:
	  tmp=exp(-fabs(y[nm]-bt)/lambda)/2;
	  H=y[nm]<bt?tmp:1-tmp;
	  l1+=log(tmp/lambda);
	  break;
	case 9:
	  H=pnorm(ly,bt,lambda);
	  l1+=log(dnorm(ly,bt,lambda)/y[nm]);
	  break;
	case 10:
	  H=plogis(ly,bt,lambda);
	  l1+=log(dlogis(ly,bt,lambda)/y[nm]);
	  break;
	case 11:
	  H=pcauchy(ly,bt,lambda);
	  l1+=log(dcauchy(ly,bt,lambda)/y[nm]);
	  break;
	case 12:
	  tmp=exp(-fabs(ly-bt)/lambda)/2;
	  H=ly<bt?tmp:1-tmp;
	  l1+=log(tmp/lambda/y[nm]);
	  break;}}
      /* calculate likelihood */
      *like-=l1;
      /* calculate fitted values */
      if(*fit){
	pred[nm]=bt;
	tmp=(b1+1/(nn*delta))/(1/(nn*delta)+j+1);
	if(!*density){
	  switch(*model){
	  case 1: rpred[nm]=bt*tmp; break;
	  case 2: rpred[nm]=bt*pow(tmp,1/lambda); break;
	  case 3: rpred[nm]=qgamma(1-exp(-tmp),lambda,bt); break;
	  case 5: rpred[nm]=qnorm(1-exp(-tmp),bt,lambda); break;
	  case 6: rpred[nm]=qlogis(1-exp(-tmp),bt,lambda); break;
	  case 7: rpred[nm]=qcauchy(1-exp(-tmp),bt,lambda); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?exp(-tmp):1-exp(-tmp)));
	    break;
	  case 9: rpred[nm]=exp(qnorm(1-exp(-tmp),bt,lambda)); break;
	  case 10: rpred[nm]=exp(qlogis(1-exp(-tmp),bt,lambda)); break;
	  case 11: rpred[nm]=exp(qcauchy(1-exp(-tmp),bt,lambda)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	else{
	  switch(*model){
	  case 1: rpred[nm]=qexp(tmp,bt); break;
	  case 2: rpred[nm]=qweibull(tmp,lambda,bt); break;
	  case 3: rpred[nm]=qgamma(tmp,lambda,bt); break;
	  case 5: rpred[nm]=qnorm(tmp,bt,lambda); break;
	  case 6: rpred[nm]=qlogis(tmp,bt,lambda); break;
	  case 7: rpred[nm]=qcauchy(tmp,bt,lambda); break;
	  case 8: rpred[nm]=bt+lambda*log(2*(y[nm]<bt?tmp:1-tmp)); break;
	  case 9: rpred[nm]=exp(qnorm(tmp,bt,lambda)); break;
	  case 10: rpred[nm]=exp(qlogis(tmp,bt,lambda)); break;
	  case 11: rpred[nm]=exp(qcauchy(tmp,bt,lambda)); break;
	  case 12: rpred[nm]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
      b1+=H;
      nm++;}
    *like+=(nobs[i]+1/delta)*log(1+delta*b1);}
  return;}
