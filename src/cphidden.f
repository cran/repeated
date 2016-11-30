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
c    subroutine chidden(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,gamma2,
c   +     val,vec,invec,model,lgam,ismu,mu,cmu,tvmu,pshape,pfam,ppar,par,
c   +     delta,nn,filter,cf,a,b,c,gmod,rhs,pivot,qraux,work,work2,work3,
c   +     like)
c
c  DESCRIPTION
c
c    Function to compute the likelihood of a two-state hidden Markov chain
c  for finding a changepoint, with various response types in continuous time
c
      subroutine cphidden(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,
     +     gamma2,val,vec,invec,model,lgam,ismu,mu,cmu,tvmu,pshape,pfam,
     +     ppar,par,delta,nn,filter,cf,a,b,c,gmod,rhs,pivot,qraux,work,
     +     work2,work3,like)
*************************************************************************
*     Function chidden computes minus the log likelihood of a           *
*     multivariate hidden Markov model with m states and iq individuals *
*     in continuous time.                                               *
*                                                                       *
*     Feel free to use or improve this program, provided that the       *
*     origin is acknowledged.                                           * 
*                              Iain MacDonald and Walter Zucchini       *
*     Modified by J.K. Lindsey for R, March, November, December 1998    *
*                                      December, 1999, January, 2000    *
*                                      November, 2001, April, 2003      *
*************************************************************************
      implicit none
      integer n(*),m,iq,i,j,k,l,model,nobs(*),mobs,ii,nm,nn,pivot(2)
      logical cf,ismu,ppar
      double precision like,s(*),pi,sflog,av,tt,pshape(2),pgamma(2,2),
     +     gmod(2,2),rhs(2),qraux(2),mu(nn,2,l),cmu(iq,2,l),par,
     +     tvmu(mobs,2,l),ll,x(*),gamma(2,2),delta(2),a(2),b(2,2),
     +     c(2),val(2),vec(2,2),invec(2,2),times(*),filter(2,nn),pfam,
     +     tmp,tmp2,tmp3,lgam(*),work(4),work2(2),work3(2,2),
     +     gamma2(2,2)
      double precision bernpr,poispr,multpr,cmultpr,contpr,proppr,binpr,
     +     exppr,bbinpr,nbinpr,normpr,invgpr,logispr,cauchpr,laplpr,
     +     levypr,paretpr,gammpr,weibpr,ggampr,glogpr,hjorpr,burrpr,
     +     gweipr,gextpr,ginvgpr,powexpr,slaplpr,studpr

      gamma2(1,2)=dexp(x(1))
      gamma2(1,1)=-gamma2(1,2)
      gamma2(2,1)=0.0
      gamma2(2,2)=0.0

c find eigenvalues/vectors and calculate gamma for unit time

      call geigen(gamma2, val, vec, invec, a, c, gmod, pivot, qraux,
     +     work, work3, 2)
      call mexp(gamma, val, vec, invec, 1.0d0, m, .false.)

c initial distribution

      delta(1)=1.0
      delta(2)=0.0

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
 203        pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,2,l,mobs,nm,nn)
            goto 250
 204        pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,1,iq,2,l,mobs,nm,nn,
     +           lgam)
            goto 250
 205        pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,2,l,mobs,nm,nn)
            goto 250
 206        pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,2,l,mobs,nm,nn)
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
            tt=times(nm)-times(nm-1)
            if(tt.ne.1.0)then
               call mexp(gmod, val, vec, invec, tt, m, .true.)
            endif
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
                  tmp3=par**tt*work2(j)
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
 303           pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,2,l,mobs,nm,
     +              nn)
               goto 350
 304           pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,k,iq,2,l,mobs,nm,
     +              nn,lgam)
               goto 350
 305           pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,2,l,mobs,nm,
     +              nn)
               goto 350
 306           pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,2,l,mobs,nm,
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
 350           if(tt.eq.1.0)then
                  do 60 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 60               continue
               else
                  do 64 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 64               continue
               endif
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
            av = dlog(av/dble(2))
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
 130     continue
      endif
      return
      end
