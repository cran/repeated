c
c  repeated : A Library of Repeated Measurements Models
c  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
c
c  This program is free software; you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation; either version 2 of the License, or
c  (at your option) any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c
c  You should have received a copy of the GNU General Public License
c  along with this program; if not, write to the Free Software
c  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c  SYNOPSIS
c
c    subroutine hidden_f(x,m,iq,nobs,mobs,s,n,l,pgamma,pos,gamma,model,
c   +     lgam,cmu,tvmu,pshape,pfam,ppar,par,delta,nn,filter,cf,a,b,c,gmod,
c   +     rhs,pivot,qraux,work,work2,like)
c
c  DESCRIPTION
c
c    Function to compute the likelihood of a hidden Markov chain model
c  with various response types in discrete time
c
      subroutine hidden_f(x,m,iq,nobs,mobs,s,n,l,pgamma,pos,gamma,model,
     +     lgam,ismu,mu,cmu,tvmu,pshape,pfam,ppar,par,delta,nn,filter,
     +     cf,a,b,c,gmod,rhs,pivot,qraux,work,work2,like)
*************************************************************************
*     Function hidden computes minus the log likelihood of a            *
*     multivariate hidden Markov model with m states and iq individuals.*
*                                                                       *
*     Feel free to use or improve this program, provided that the       *
*     origin is acknowledged.                                           * 
*                              Iain MacDonald and Walter Zucchini       *
*     Modified by J.K. Lindsey for R, March, 1998, November, 1998       *
*                                   December, 1999, January, 2000       *
*     and by P.J. Lindsey for ordinal responses, December, 1999         *
*************************************************************************
      implicit none
      integer n(*),m,iq,i,j,k,l,model,nobs(*),mobs,ii,nm,pos(m),nn,
     +     pivot(m)
      logical cf,ismu,ppar
      double precision like,s(*),pi,sflog,av,pshape(m),par,
     +     gmod(m,m),rhs(m),qraux(m),mu(nn,m,l),cmu(iq,m,l),
     +     tvmu(mobs,m,l),ll,x(*),gamma(m,m),delta(m),a(m),b(m,m),
     +     c(m),pgamma(m,m),filter(m,nn),pfam,tmp,tmp2,tmp3,lgam(*),
     +     work(2*m),work2(m)
      double precision bernpr,poispr,multpr,cmultpr,contpr,proppr,binpr,
     +     exppr,bbinpr,nbinpr,normpr,invgpr,logispr,cauchpr,laplpr,
     +     levypr,paretpr,gammpr,weibpr,ggampr,glogpr,hjorpr,burrpr,
     +     gweipr,gextpr,ginvgpr,powexpr,slaplpr,studpr

      call fromx(x,m,gamma,pgamma,pos)

c calculate stationary distribution

      call deltas(gamma,delta,m,gmod,rhs,pivot,qraux,work)

c take logs of probabilities

      do 10 i=1,m
         delta(i)=dlog(delta(i))
         do 11 j=1,m
            gamma(i,j)=dlog(gamma(i,j))
 11      continue
 10   continue

c initial conditions

      like=0.
      nm=0
      do 20 i = 1, iq
         if(ppar)then
            do 21 j=1,m
               work2(j)=0
 21         continue
         endif
         nm=nm+1
         do 30 j = 1, m
            a(j)=delta(j)
            if(model.lt.3.or.model.gt.6)then
               if(ismu)then
                  tmp=mu(nm,j,1)
               else
                  tmp=cmu(i,j,1)+tvmu(1,j,1)
               endif
            endif
            if(ppar)work2(j)=s(nm)-tmp
            goto(201,202,203,204,205,206,207,208,209,210,211,212,213,
     +           214,215,216,217,218,219,220,221,222,223,224,225,226,
     +           227,228,229),model
 201        pi = bernpr(s(nm),tmp)
            goto 250
 202        pi = poispr(s(nm),tmp)
            goto 250
 203        pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
            goto 250
 204        pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn,
     +           lgam)
            goto 250
 205        pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
            goto 250
 206        pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
            goto 250
 207        pi = binpr(s(nm),n(nm),tmp)
            goto 250
 208        pi = exppr(s(nm),tmp)
            goto 250
 209        pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
            goto 250
 210        pi = nbinpr(s(nm),tmp,pshape(j))
            goto 250
 211        pi = normpr(s(nm),tmp,pshape(j))
            goto 250
 212        pi = invgpr(s(nm),tmp,pshape(j))
            goto 250
 213        pi = logispr(s(nm),tmp,pshape(j))
            goto 250
 214        pi = cauchpr(s(nm),tmp,pshape(j))
            goto 250
 215        pi = laplpr(s(nm),tmp,pshape(j))
            goto 250
 216        pi = levypr(s(nm),tmp,pshape(j))
            goto 250
 217        pi = paretpr(s(nm),tmp,pshape(j))
            goto 250
 218        pi = gammpr(s(nm),tmp,pshape(j))
            goto 250
 219        pi = weibpr(s(nm),tmp,pshape(j))
            goto 250
 220        pi = ggampr(s(nm),tmp,pshape(j),pfam)
            goto 250
 221        pi = glogpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 222        pi = hjorpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 223        pi = burrpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 224        pi = gweipr(s(nm),tmp,pshape(j),pfam)
            goto 250
 225        pi = gextpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 226        pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 227        pi = powexpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 228        pi = slaplpr(s(nm),tmp,pshape(j),pfam)
            goto 250
 229        pi = studpr(s(nm),tmp,pshape(j),pfam)
 250        a(j) = a(j) + pi
 30      continue

c filtered conditional probabilities of states

         if(cf)then
            ll = 0.
            do 31 j = 1, m
               filter(j,nm)=dexp(a(j))
               ll = ll + filter(j,nm)
 31         continue
            do 32 j = 1, m
               filter(j,nm)=filter(j,nm)/ll
 32         continue
         endif

c update likelihood at each subsequent time point

         sflog = 0.
         do 110 k = 2, nobs(i)
            nm=nm+1
            do 70 j = 1, m
               if(model.lt.3.or.model.gt.6)then
                  if(ismu)then
                     tmp=mu(nm,j,1)
                  else
                     tmp=cmu(i,j,1)+tvmu(k,j,1)
                  endif
               endif
               if(ppar)then
                  tmp2=tmp
                  tmp3=par*work2(j)
                  if(model.ne.11.and.model.ne.13.and.model.ne.14.and.
     +                 model.ne.15.and.model.ne.16.and.model.ne.21.and.
     +                 model.ne.27.and.model.ne.28.and.model.ne.29.and.
     +                 tmp+tmp3.le.0.0)tmp3=0.0
                  if(model.eq.1.and.tmp+tmp3.ge.1.0)tmp3=0.0
                  tmp=tmp+tmp3
                  work2(j)=s(nm)-tmp2
               endif
               goto(301,302,303,304,305,306,307,308,309,310,311,312,313,
     +              314,315,316,317,318,319,320,321,322,323,324,325,326,
     +              327,328,329),model
 301           pi = bernpr(s(nm),tmp)
               goto 350
 302           pi = poispr(s(nm),tmp)
               goto 350
 303           pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
     +              nn)
               goto 350
 304           pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
     +              nn,lgam)
               goto 350
 305           pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
     +              nn)
               goto 350
 306           pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
     +              nn)
               goto 350
 307           pi = binpr(s(nm),n(nm),tmp)
               goto 350
 308           pi = exppr(s(nm),tmp)
               goto 350
 309           pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
               goto 350
 310           pi = nbinpr(s(nm),tmp,pshape(j))
               goto 350
 311           pi = normpr(s(nm),tmp,pshape(j))
               goto 350
 312           pi = invgpr(s(nm),tmp,pshape(j))
               goto 350
 313           pi = logispr(s(nm),tmp,pshape(j))
               goto 350
 314           pi = cauchpr(s(nm),tmp,pshape(j))
               goto 350
 315           pi = laplpr(s(nm),tmp,pshape(j))
               goto 350
 316           pi = levypr(s(nm),tmp,pshape(j))
               goto 350
 317           pi = paretpr(s(nm),tmp,pshape(j))
               goto 350
 318           pi = gammpr(s(nm),tmp,pshape(j))
               goto 350
 319           pi = weibpr(s(nm),tmp,pshape(j))
               goto 350
 320           pi = ggampr(s(nm),tmp,pshape(j),pfam)
               goto 350
 321           pi = glogpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 322           pi = hjorpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 323           pi = burrpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 324           pi = gweipr(s(nm),tmp,pshape(j),pfam)
               goto 350
 325           pi = gextpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 326           pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 327           pi = powexpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 328           pi = slaplpr(s(nm),tmp,pshape(j),pfam)
               goto 350
 329           pi = studpr(s(nm),tmp,pshape(j),pfam)
 350           do 60 ii = 1, m
                  b(ii,j) = gamma(ii,j)+pi
 60            continue
 70         continue

c normalize to prevent underflow

            av = 0.
            do 90 j = 1, m
               c(j) = 0.
               do 80 ii = 1, m
                  c(j) = c(j) + dexp(a(ii)+b(ii,j))
 80            continue
               av = av + c(j)
 90         continue
            av = dlog(av/dble(m))
            do 100 j = 1, m
               a(j) = dlog(c(j))-av
 100        continue

c correction factor for normalization

            sflog = sflog + av

c filtered conditional probabilities of states

            if(cf)then
               ll = 0.
               do 101 j = 1, m
                  filter(j,nm)=dexp(a(j))
                  ll = ll + filter(j,nm)
 101           continue
               do 102 j = 1, m
                  filter(j,nm)=filter(j,nm)/ll
 102           continue
            endif
 110     continue

c calculate likelihood including correction factor

         ll = 0.
         do 120 j = 1, m
            ll = ll + dexp(a(j))
 120     continue
         like = like-(dlog(ll)+sflog)
 20   continue

c transform back to original values

      if(cf)then
         do 130 i = 1, m
            delta(i)=dexp(delta(i))
            do 131 j=1,m
               if(pgamma(i,j).le.1.d-30)then
                  gamma(i,j)=0.
               else if(pgamma(i,j).eq.1.)then
                  gamma(i,j)=1.
               else
                  gamma(i,j)=dexp(gamma(i,j))
               endif
 131        continue
 130     continue
      endif

      return
      end

      subroutine fromx(x, m, gamma, pgamma, pos)
***************************************************************************
*    Converts the vector of parameters into transition probability matrix *
***************************************************************************
      implicit none
      integer m,i,ii,j,pos(m)
      double precision sum,gamma(m,m),x(*),pgamma(m,m)

      ii=0
      do 30 i = 1, m
         sum = 1.
         do 20 j = 1, m
            if(j.eq.pos(i))then
               gamma(i,j) = 1.
            else if(pgamma(i,j).lt.1.d-30.or.pgamma(i,j).eq.1.)then
               gamma(i,j)=pgamma(i,j)
            else
               ii=ii+1
               gamma(i,j) =dexp(x(ii))
               sum = sum + gamma(i,j)
            endif
 20      continue
         do 40 j=1,m
            if(pgamma(i,j).gt.1.d-30.and.pgamma(i,j).ne.1.)
     +           gamma(i,j)=gamma(i,j)/sum
 40      continue
 30   continue
      return
      end

      subroutine deltas(gamma,delta,m,gmod,rhs,pivot,qraux,work)
****************************************************
* Computes stationary distribution of Markov chain *
****************************************************
      implicit none
      integer m,i,j,rank,info,pivot(*)
      double precision qraux(*),work(*),dummy,
     +     gamma(m,m),delta(*),gmod(m,m),rhs(*)

      do 20 i = 2, m
         do 10 j = 1, m
            gmod(i,j) = gamma(j,i)
 10      continue
         gmod(i,i) = gmod(i,i) - 1.
         rhs(i) = 0.
 20   continue
*
      do 30 j = 1, m
         pivot(j)=j
         gmod(1,j) = 1.
 30   continue
      rhs(1) = 1.
      call dqrdc2(gmod,m,m,m,1d-07,rank,qraux,pivot,work)
c      call dqrcf(gmod,m,rank,qraux,rhs,1,delta,info)
      call dqrsl(gmod,m,m,rank,qraux,rhs,dummy,rhs,delta,
     +     dummy,dummy,100,info)
      return
      end

      double precision function bernpr(svec,pvec)
***************************
* Bernoulli probabilities *
***************************
      implicit none
      double precision pi,svec,pvec

      pi = 1./(1.+dexp(-pvec))
      if(svec.eq.1.and.pi.gt.0.) then
         bernpr=dlog(pi)
      else if(svec.eq.0.and.pi.lt.1.) then
         bernpr=dlog(1.-pi)
      else
         bernpr=-35.
      endif
      return
      end

      double precision function poispr(svec,pvec)
*************************
* Poisson probabilities *
*************************
      implicit none
      integer j
      double precision svec,pvec

      poispr = -pvec
      do 1 j = 1, int(svec)
         poispr = poispr+dlog(pvec/dble(j))
 1    continue
      return
      end

C      do 1 j = 1, 999999999
C         if(j > svec) exit
C         poispr = poispr+dlog(pvec/dble(j))
C 1    continue
C      return
C      end

      double precision function multpr(svec,ismu,mu,cmu,tvmu,i,j,k,
     +     iq,m,l,mobs,nm,nn)
*****************************
* Multinomial probabilities *
*****************************
      implicit none
      integer i,j,k,iq,m,l,mobs,n,nn,nm
      double precision svec,mu(nn,m,l),cmu(iq,m,l),tvmu(mobs,m,l),
     +     denom
      logical ismu

      denom=1.
      if(ismu)then
         do 1 n=1,l
            denom=denom+dexp(mu(nm,j,n))
 1       continue
         if(svec.gt.0)then
            multpr = mu(nm,j,int(svec))-dlog(denom)
         else
            multpr = -dlog(denom)
         endif
      else
         do 2 n=1,l
            denom=denom+dexp(cmu(i,j,n)+tvmu(k,j,n))
 2       continue
         if(svec.gt.0)then
            multpr = cmu(i,j,int(svec))+tvmu(k,j,int(svec))-dlog(denom)
         else
            multpr = -dlog(denom)
         endif
      endif
      return
      end

      double precision function cmultpr(svec,ismu,mu,cmu,tvmu,i,j,k,
     +     iq,m,l,mobs,nm,nn,lgam)
***********************************
* Multinomial count probabilities *
***********************************
      implicit none
      integer i,j,k,iq,m,l,mobs,n,nn,nm
      double precision mu(nn,m,l),cmu(iq,m,l),tvmu(mobs,m,l),
     +     svec(*),lgam(*),pi,denom,sum
      logical ismu

      denom=1.
      pi=lgam(nm)
      sum=svec(1+(nm-1)*(l+1))
      if(ismu)then
         do 1 n=1,l
            sum=sum+svec(n+1+(nm-1)*(l+1))
            pi=pi+svec(n+1+(nm-1)*(l+1))*mu(nm,j,n)
            denom=denom+dexp(mu(nm,j,n))
 1       continue
      else
         do 2 n=1,l
            sum=sum+svec(n+1+(nm-1)*(l+1))
            pi=pi+svec(n+1+(nm-1)*(l+1))*(cmu(i,j,n)+tvmu(k,j,n))
            denom=denom+dexp(cmu(i,j,n)+tvmu(k,j,n))
 2       continue
      endif
      cmultpr=pi-sum*dlog(denom)
      return
      end

      double precision function contpr(svec,ismu,mu,cmu,tvmu,i,j,k,
     +     iq,m,l,mobs,nm,nn)
**********************************************
* Continuation-ratio probabilities (upwards) *
**********************************************
      implicit none
      integer i,j,k,iq,m,l,mobs,n,nn,nm
      double precision svec,mu(nn,m,l),cmu(iq,m,l),tvmu(mobs,m,l),pi
      logical ismu

      if(ismu)then
         if(svec.eq.0)then
            pi=1.
         else
            pi=1./(1.+dexp(mu(nm,j,int(svec))))
         endif
         do 2 n=int(svec)+1,l
            pi=pi/(1.+dexp(-mu(nm,j,n)))
 2       continue
      else
         if(svec.eq.0)then
            pi=1.
         else
            pi=1./(1.+dexp(cmu(i,j,int(svec))+tvmu(k,j,int(svec))))
         endif
         do 4 n=int(svec)+1,l
            pi=pi/(1.+dexp(-cmu(i,j,n)-tvmu(k,j,n)))
 4       continue
      endif
      if(pi.gt.0.) then
         contpr=dlog(pi)
      else
         contpr=-35.
      endif
      return
      end

      double precision function proppr(svec,ismu,mu,cmu,tvmu,i,j,k,
     +     iq,m,l,mobs,nm,nn)
***********************************
* Proportional-odds probabilities *
***********************************
      implicit none
      integer i,j,k,iq,m,l,mobs,nn,nm
      double precision svec,mu(nn,m,l),cmu(iq,m,l),tvmu(mobs,m,l),pi
      logical ismu

      if(ismu)then
         if(svec.eq.l)then
            pi=1./(1.+dexp(mu(nm,j,l)))
         else if(svec.gt.0)then
            pi=1/(1.+dexp(-mu(nm,j,int(svec+1))))
     +           -1/(1.+dexp(-mu(nm,j,int(svec))))
         else
            pi=1/(1.+dexp(-mu(nm,j,1)))
         endif
      else
         if(svec.eq.l)then
            pi=1./(1.+dexp(cmu(i,j,l)+tvmu(k,j,l)))
         else if(svec.gt.0)then
            pi=1/(1.+dexp(-cmu(i,j,int(svec+1))-tvmu(k,j,int(svec+1))))
     +           -1/(1.+dexp(-cmu(i,j,int(svec))-tvmu(k,j,int(svec))))
         else
            pi=1/(1.+dexp(-cmu(i,j,1)-tvmu(k,j,1)))
         endif
      endif
      if(pi.gt.0.) then
         proppr=dlog(pi)
      else
         proppr=-35.
      endif
      return
      end

      double precision function binpr(svec,nvec,pvec)
**************************
* Binomial probabilities *
**************************
      implicit none
      integer nvec,j
      double precision svec,pi,pvec,tmp,tmp2

      pi = 1./(1.+dexp(-pvec))
      if(pi.ne.0..and.pi.ne.1.) then
         if(svec.eq.0) then
            binpr=nvec*dlog(1.-pi)
         else if(svec.eq.nvec) then
            binpr=svec*dlog(pi)
         else
            binpr=1.0
            tmp2=nvec+1.0
            if(svec.lt.nvec/2)then
               tmp=svec+1.0
               do 1 j = 1, int(svec)
                  binpr=binpr*(tmp2-j)/(tmp-j)
 1             continue
            else
               tmp=nvec-svec+1.0
               do 2 j = 1, int(nvec-svec)
                  binpr=binpr*(tmp2-j)/(tmp-j)
 2             continue
            endif
            binpr=dlog(binpr)+svec*dlog(pi)+(nvec-svec)*dlog(1.-pi)
         endif
      else
         binpr=-35.
      endif
      return
      end

      double precision function exppr(svec,pvec)
*****************************
* Exponential probabilities *
*****************************
      implicit none
      double precision svec,pvec

      exppr = dlog(pvec)-svec*pvec
      return
      end

      double precision function bbinpr(svec,nvec,pvec,pshape)
*******************************
* Beta-binomial probabilities *
*******************************
      implicit none
      integer nvec,j
      double precision svec,pi,pi2,pvec,tmp,tmp2,tmp3,pshape

      pi = 1./(1.+dexp(-pvec))
      pi2 = pshape*(1-pi)
      pi = pshape*pi
      if(svec+pi.gt.0..and.nvec-svec+pi2.gt.0.)then
         call flbeta(svec+pi,nvec-svec+pi2,bbinpr)
      else
         bbinpr=0.
      endif
      if(pi.gt.0..and.pi2.gt.0.)then
         call flbeta(pi,pi2,tmp)
      else
         tmp=-35.
      endif
      bbinpr=bbinpr-tmp
      tmp=1.0
      tmp2=nvec+1.0
      if(svec.lt.nvec/2)then
         tmp3=svec+1.0
         do 1 j = 1, int(svec)
            tmp=tmp*(tmp2-j)/(tmp3-j)
 1       continue
      else
         tmp3=nvec-svec+1.0
         do 2 j = 1, int(nvec-svec)
            tmp=tmp*(tmp2-j)/(tmp3-j)
 2       continue
      endif
      bbinpr=bbinpr+dlog(tmp)
      return
      end

      double precision function nbinpr(svec,pvec,pshape)
***********************************
* Negative binomial probabilities *
***********************************
      implicit none
      double precision svec,pvec,pshape,tmp

      call flgamma(svec+pshape,nbinpr)
      call flgamma(pshape,tmp)
      nbinpr=nbinpr-tmp
      call flgamma(svec+1.,tmp)
      nbinpr=nbinpr-tmp
      nbinpr=nbinpr+pshape*dlog(pshape)+svec*dlog(pvec)
     +     -(pshape+svec)*dlog(pshape+pvec)
      return
      end

      double precision function normpr(svec,pvec,pshape)
************************
* Normal probabilities *
************************
      implicit none
      double precision svec,pvec,pshape

      normpr = -0.5*(dlog(2*3.14159265358979323846264338327950288419
     +     7169399375*pshape)+(svec-pvec)**2/pshape)
      return
      end

      double precision function invgpr(svec,pvec,pshape)
*******************************
* Inverse Gauss probabilities *
*******************************
      implicit none
      double precision svec,pvec,pshape

      invgpr = -0.5*(dlog(2*3.14159265358979323846264338327950288419
     +     7169399375*pshape*svec**3)+(svec-pvec)**2/
     +     (svec*pshape*pvec**2))
      return
      end

      double precision function logispr(svec,pvec,pshape)
**************************
* Logistic probabilities *
**************************
      implicit none
      double precision svec,pvec,pshape,tmp

      tmp=pshape*dsqrt(3.d0)/3.14159265358979323846264338327950288419
     +     7169399375
      logispr=(svec-pvec)/tmp
      logispr = -logispr-dlog(tmp)-2.0*dlog(1+dexp(-logispr))
      return
      end

      double precision function cauchpr(svec,pvec,pshape)
************************
* Cauchy probabilities *
************************
      implicit none
      double precision svec,pvec,pshape

      cauchpr = -dlog(pshape*(1+((svec-pvec)/pshape)**2)*3.141
     +     592653589793238462643383279502884197169399375)
      return
      end

      double precision function laplpr(svec,pvec,pshape)
*************************
* Laplace probabilities *
*************************
      implicit none
      double precision svec,pvec,pshape

      laplpr = -dabs(svec-pvec)/pshape-dlog(pshape*2.0)
      return
      end

      double precision function levypr(svec,pvec,pshape)
*************************
* Levy probabilities *
*************************
      implicit none
      double precision svec,pvec,pshape

      levypr = 0.5*dlog(pshape/(2.0*3.14159265358979323846264
     +     3383279502884197169399375))-1.5*dlog(svec-pvec)
     +     -pshape/(2*(svec-pvec))-dlog(2.d0)
      return
      end

      double precision function paretpr(svec,pvec,pshape)
************************
* Pareto probabilities *
************************
      implicit none
      double precision svec,pvec,pshape

      paretpr = 1.0/(pvec*(pshape-1.0))
      paretpr = dlog(pshape*paretpr)-(pshape+1.0)*dlog(1+svec*paretpr)
      return
      end

      double precision function gammpr(svec,pvec,pshape)
***********************
* Gamma probabilities *
***********************
      implicit none
      double precision svec,pvec,pshape

      call flgamma(pshape,gammpr)
      gammpr = pshape*(dlog(pshape/pvec)-svec/pvec)+(pshape-1.0)
     +     *dlog(svec)-gammpr
      return
      end

      double precision function weibpr(svec,pvec,pshape)
*************************
* Weibull probabilities *
*************************
      implicit none
      double precision svec,pvec,pshape

      weibpr = dlog(pshape)+(pshape-1.0)*dlog(svec)-pshape*dlog(pvec)
     +     -(svec/pvec)**pshape
      return
      end

      double precision function ggampr(svec,pvec,pshape,pfam)
***********************************
* Generalized gamma probabilities *
***********************************
      implicit none
      double precision svec,pvec,pshape,pfam

      call flgamma(pshape,ggampr)
      ggampr=pshape*pfam*dlog(pshape/pvec)
     +     -(pshape*svec/pvec)**pfam+log(pfam)
     +     +(pshape*pfam-1)*dlog(svec)-ggampr
      return
      end

      double precision function glogpr(svec,pvec,pshape,pfam)
**************************************
* Generalized logistic probabilities *
**************************************
      implicit none
      double precision svec,pvec,pshape,pfam

      glogpr=(svec-pvec)/pshape
      glogpr=dlog(pfam)-glogpr-dlog(pshape)-(pfam+1)*
     +     dlog(1+dexp(-glogpr))
      return
      end

      double precision function hjorpr(svec,pvec,pshape,pfam)
************************
* Hjorth probabilities *
************************
      implicit none
      double precision svec,pvec,pshape,pfam

      hjorpr=-pfam*dlog(1.0+pshape*svec)/pshape-(svec/pvec)**2/2.0
     +	    +dlog(svec/(pvec*pvec)+pfam/(1+pshape*svec))
      return
      end

      double precision function burrpr(svec,pvec,pshape,pfam)
**********************
* Burr probabilities *
**********************
      implicit none
      double precision svec,pvec,pshape,pfam

      burrpr=svec/pvec
      burrpr=dlog(pfam*pshape/pvec)+(pshape-1.0)*dlog(burrpr)-
     +	    (pfam+1.0)*dlog(1.0+burrpr**pshape)
      return
      end

      double precision function gweipr(svec,pvec,pshape,pfam)
*************************************
* Generalized Weibull probabilities *
*************************************
      implicit none
      double precision svec,pvec,pshape,pfam,tmp

      tmp=(svec/pvec)**pshape
      gweipr=dlog(pshape*pfam)+(pshape-1.0)*dlog(svec)
     +     -pshape*dlog(pvec)+(pfam-1.0)*dlog(1.0-dexp(-tmp))-tmp
      return
      end

      double precision function gextpr(svec,pvec,pshape,pfam)
*******************************************
* Generalized extreme value probabilities *
*******************************************
      implicit none
      double precision svec,pvec,pshape,pfam,tmp

      if(pfam.gt.0)then
         tmp=-pvec**(-pshape)
      else
         tmp=dlog(1.0-dexp(-pvec**(-pshape)))
      endif
      gextpr=svec**pfam/pfam
      gextpr=dlog(pshape)+pshape*(gextpr-dlog(pvec))-tmp
     +	    -(dexp(gextpr)/pvec)**pshape+(pfam-1.0)*dlog(svec)
      return
      end

      double precision function ginvgpr(svec,pvec,pshape,pfam)
*******************************************
* Generalized inverse Gauss probabilities *
*******************************************
      implicit none
      double precision svec,pvec,pshape,pfam

      call fbesselk(1/(pshape*pvec),dabs(pfam),ginvgpr)
      ginvgpr=(pfam-1.0)*dlog(svec)-(1.0/svec+svec/(pvec*pvec))
     +	    /(2.0*pshape)-pfam*dlog(pvec)-dlog(2*ginvgpr)
      return
      end

      double precision function powexpr(svec,pvec,pshape,pfam)
***********************************
* Power exponential probabilities *
***********************************
      implicit none
      double precision svec,pvec,pshape,pfam,tmp1,tmp2

      tmp1=dsqrt(pshape)
      tmp2=1.0+1.0/(2.0*dabs(pfam))
      call flgamma(tmp2,powexpr)
      powexpr=-(dabs(svec-pvec)/tmp1)**(2.0*dabs(pfam))/2.0
     +     -dlog(tmp1*2.0**tmp2)-powexpr
      return
      end

      double precision function slaplpr(svec,pvec,pshape,pfam)
******************************
* skew Laplace probabilities *
******************************
      implicit none
      double precision svec,pvec,pshape,pfam,tmp

      if(svec.gt.pvec)then
         tmp=-pfam*(svec-pvec)
      else
         tmp=(svec-pvec)/pfam
      endif
      slaplpr = dlog(pfam)+tmp/pshape-dlog((1+pfam**2)*pshape)
      return
      end

      double precision function studpr(svec,pvec,pshape,pfam)
***************************
* Student t probabilities *
***************************
      implicit none
      double precision svec,pvec,pshape,pfam,tmp

      call flgamma(0.5*(pfam+1.0),studpr)
      call flgamma(0.5*pfam,tmp)
      studpr = studpr-tmp-0.5*(dlog(pfam*3.1415926535897932384626433
     +     83279502884197169399375)+(pfam+1.0)*dlog(1+((svec-pvec)
     +     /pshape)**2/pfam))-dlog(pshape)
      return
      end
