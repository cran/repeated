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
c    subroutine chidden_f(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,gamma2,
c   +     val,vec,invec,model,lgam,ismu,mu,cmu,tvmu,pshape,pfam,ppar,par,
c   +     delta,nn,filter,cf,a,b,c,gmod,rhs,pivot,qraux,work,work2,work3,
c   +     like)
c
c  DESCRIPTION
c
c    Function to compute the likelihood of a hidden Markov chain model
c  with various response types in continuous time
c
      subroutine chidden_f(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,
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
*                                      November, 2001                   *
*************************************************************************
      implicit none
      integer n(*),m,iq,i,j,k,l,model,nobs(*),mobs,ii,nm,nn,pivot(m)
      logical cf,ismu,ppar
      double precision like,s(*),pi,sflog,av,tt,pshape(m),pgamma(m,m),
     +     gmod(m,m),rhs(m),qraux(m),mu(nn,m,l),cmu(iq,m,l),par,
     +     tvmu(mobs,m,l),ll,x(*),gamma(m,m),delta(m),a(m),b(m,m),
     +     c(m),val(m),vec(m,m),invec(m,m),times(*),filter(m,nn),pfam,
     +     tmp,tmp2,tmp3,lgam(*),work(2*m),work2(m),work3(m,m),
     +     gamma2(m,m)
      double precision bernpr,poispr,multpr,cmultpr,contpr,proppr,binpr,
     +     exppr,bbinpr,nbinpr,normpr,invgpr,logispr,cauchpr,laplpr,
     +     levypr,paretpr,gammpr,weibpr,ggampr,glogpr,hjorpr,burrpr,
     +     gweipr,gextpr,ginvgpr,powexpr,slaplpr,studpr

      call cfromx(x,m,gamma2,pgamma)

c find eigenvalues/vectors and calculate gamma for unit time

      call geigen(gamma2, val, vec, invec, a, c, gmod, pivot, qraux,
     +     work, work3, m)
      call mexp(gamma, val, vec, invec, 1.0d0, m, .false.)

c calculate stationary distribution

      call deltas(gamma,delta,m,gmod,rhs,pivot,qraux,work)

c special case: discretized Poisson process
c G: gamma2, b: G-L, gmod: exp(G), gamma: exp(G-L) or exp(G)-exp(G-L)

      if(model.eq.30)then
         do 3 i=1,m
            c(i)=gamma2(i,i)
            do 4 j=1,m
               b(i,j)=gamma2(i,j)
 4          continue
 3       continue
         like=0.
         nm=0
         do 1 i = 1, iq
            do 12 j=1,m
               a(j)=delta(j)
 12         continue
            do 2 k = 1, nobs(i)
               nm=nm+1
               if(k.eq.1)then
                  tt=times(nm)
               else
                  tt=times(nm)-times(nm-1)
               endif
               do 5 j=1,m
                  if(ismu)then
                     b(j,j)=c(j)-exp(mu(nm,j,1))
                  else
                     b(j,j)=c(j)-exp(cmu(i,j,1)+tvmu(k,j,1))
                  endif
 5             continue
               call geigen(b, val, vec, invec, rhs, work2, gmod, pivot,
     +              qraux, work, work3, m)
               call mexp(gamma, val, vec, invec, tt, m, .false.)
               if(s(nm).eq.1.0)then
                  call geigen(gamma2, val, vec, invec, rhs, work2, gmod,
     +                 pivot, qraux, work, work3, m)
                  call mexp(gmod, val, vec, invec, tt, m, .false.)
                  do 8 ii=1,m
                     do 9 j=1,m
                        gamma(ii,j)=gmod(ii,j)-gamma(ii,j)
 9                   continue
 8                continue
               endif
               do 6 ii=1,m
                  filter(ii,nm)=0.0
                  do 7 j=1,m
                     filter(ii,nm)=filter(ii,nm)+a(j)*gamma(j,ii)
 7                continue
 6             continue
               tmp=0.0
               do 13 j=1,m
                  a(j)=filter(j,nm)
                  tmp=tmp+a(j)
 13            continue
               if(cf)then
                  do 15 j=1,m
                     filter(j,nm)=filter(j,nm)/tmp
 15               continue
               endif
 2          continue
            tmp=0.0
            do 16 j=1,m
               tmp=tmp+a(j)
 16         continue
            like=like-dlog(tmp)
 1       continue
         return
      endif

c all other models

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
c            goto(201,202,203,204,205,206,207,208,209,210,211,212,213,
c     +           214,215,216,217,218,219,220,221,222,223,224,225,226,
c     +           227,228,229),model
c 201        pi = bernpr(s(nm),tmp)
c            goto 250
c 202        pi = poispr(s(nm),tmp)
c            goto 250
c 203        pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
c            goto 250
c 204        pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn,
c     +           lgam)
c            goto 250
c 205        pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
c            goto 250
c 206        pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,mobs,nm,nn)
c            goto 250
c 207        pi = binpr(s(nm),n(nm),tmp)
c            goto 250
c 208        pi = exppr(s(nm),tmp)
c            goto 250
c 209        pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
c            goto 250
c 210        pi = nbinpr(s(nm),tmp,pshape(j))
c            goto 250
c 211        pi = normpr(s(nm),tmp,pshape(j))
c            goto 250
c 212        pi = invgpr(s(nm),tmp,pshape(j))
c            goto 250
c 213        pi = logispr(s(nm),tmp,pshape(j))
c            goto 250
c 214        pi = cauchpr(s(nm),tmp,pshape(j))
c            goto 250
c 215        pi = laplpr(s(nm),tmp,pshape(j))
c            goto 250
c 216        pi = levypr(s(nm),tmp,pshape(j))
c            goto 250
c 217        pi = paretpr(s(nm),tmp,pshape(j))
c            goto 250
c 218        pi = gammpr(s(nm),tmp,pshape(j))
c            goto 250
c 219        pi = weibpr(s(nm),tmp,pshape(j))
c            goto 250
c 220        pi = ggampr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 221        pi = glogpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 222        pi = hjorpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 223        pi = burrpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 224        pi = gweipr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 225        pi = gextpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 226        pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 227        pi = powexpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 228        pi = slaplpr(s(nm),tmp,pshape(j),pfam)
c            goto 250
c 229        pi = studpr(s(nm),tmp,pshape(j),pfam)
c 250        a(j) = a(j) + pi
             select case(model)
             case(1)
                 pi = bernpr(s(nm),tmp)
                 a(j) = a(j) + pi
             case(2)        
                 pi = poispr(s(nm),tmp)
                 a(j) = a(j) + pi
             case(3)        
                 pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,
     +                       mobs,nm,nn)
                 a(j) = a(j) + pi
             case(4)        
                 pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,1,
     +                        iq,m,l,mobs,nm,nn,lgam)
                 a(j) = a(j) + pi
             case(5)        
                 pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,
     +                       iq,m,l,mobs,nm,nn)
                 a(j) = a(j) + pi
             case(6)        
                 pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,1,
     +                       iq,m,l,mobs,nm,nn)
                 a(j) = a(j) + pi
             case(7)        
                 pi = binpr(s(nm),n(nm),tmp)
                 a(j) = a(j) + pi
             case(8)        
                 pi = exppr(s(nm),tmp)
                 a(j) = a(j) + pi
             case(9)        
                 pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(10)        
                 pi = nbinpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(11)        
                 pi = normpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(12)        
                 pi = invgpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(13)        
                 pi = logispr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(14)        
                 pi = cauchpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(15)        
                 pi = laplpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(16)        
                 pi = levypr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(17)        
                 pi = paretpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(18)        
                 pi = gammpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(19)        
                 pi = weibpr(s(nm),tmp,pshape(j))
                 a(j) = a(j) + pi
             case(20)        
                 pi = ggampr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(21)        
                 pi = glogpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(22)        
                 pi = hjorpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(23)        
                 pi = burrpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(24)        
                 pi = gweipr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(25)        
                 pi = gextpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(26)        
                 pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(27)        
                 pi = powexpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(28)        
                 pi = slaplpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case(29)        
                 pi = studpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             case default
                 pi = studpr(s(nm),tmp,pshape(j),pfam)
                 a(j) = a(j) + pi
             end select
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
c               goto(301,302,303,304,305,306,307,308,309,310,311,312,313,
c     +              314,315,316,317,318,319,320,321,322,323,324,325,326,
c     +              327,328,329),model
c 301           pi = bernpr(s(nm),tmp)
c               goto 350
c 302           pi = poispr(s(nm),tmp)
c               goto 350
c 303           pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
c     +              nn)
c               goto 350
c 304           pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
c     +              nn,lgam)
c               goto 350
c 305           pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
c     +              nn)
c               goto 350
c 306           pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,k,iq,m,l,mobs,nm,
c     +              nn)
c               goto 350
c 307           pi = binpr(s(nm),n(nm),tmp)
c               goto 350
c 308           pi = exppr(s(nm),tmp)
c               goto 350
c 309           pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
c               goto 350
c 310           pi = nbinpr(s(nm),tmp,pshape(j))
c               goto 350
c 311           pi = normpr(s(nm),tmp,pshape(j))
c               goto 350
c 312           pi = invgpr(s(nm),tmp,pshape(j))
c               goto 350
c 313           pi = logispr(s(nm),tmp,pshape(j))
c               goto 350
c 314           pi = cauchpr(s(nm),tmp,pshape(j))
c               goto 350
c 315           pi = laplpr(s(nm),tmp,pshape(j))
c               goto 350
c 316           pi = levypr(s(nm),tmp,pshape(j))
c               goto 350
c 317           pi = paretpr(s(nm),tmp,pshape(j))
c               goto 350
c 318           pi = gammpr(s(nm),tmp,pshape(j))
c               goto 350
c 319           pi = weibpr(s(nm),tmp,pshape(j))
c               goto 350
c 320           pi = ggampr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 321           pi = glogpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 322           pi = hjorpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 323           pi = burrpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 324           pi = gweipr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 325           pi = gextpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 326           pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 327           pi = powexpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 328           pi = slaplpr(s(nm),tmp,pshape(j),pfam)
c               goto 350
c 329           pi = studpr(s(nm),tmp,pshape(j),pfam)
c 350           if(tt.eq.1.0)then
c                  do 60 ii = 1, m
c                     b(ii,j) = gamma(ii,j)+pi
c 60               continue
c               else
c                  do 64 ii = 1, m
c                     b(ii,j) = gmod(ii,j)+pi
c 64               continue
c               endif
             select case(model)
             case(1)
                pi = bernpr(s(nm),tmp)
                if(tt.eq.1.0)then
                  do 601 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 601               continue
                else
                  do 641 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 641               continue
               endif
             case(2)        
                 pi = poispr(s(nm),tmp)
                 if(tt.eq.1.0)then
                  do 602 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 602               continue
                 else
                  do 642 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 642               continue
               endif
             case(3)        
                 pi = multpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,iq,m,l,
     +                       mobs,nm,nn)
                 if(tt.eq.1.0)then
                  do 603 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 603               continue
                 else
                  do 643 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 643               continue
               endif
             case(4)        
                 pi = cmultpr(s,ismu,mu,cmu,tvmu,i,j,1,
     +                        iq,m,l,mobs,nm,nn,lgam)
                 if(tt.eq.1.0)then
                  do 604 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 604               continue
                 else
                  do 644 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 644               continue
               endif
             case(5)        
                 pi = contpr(s(nm),ismu,mu,cmu,tvmu,i,j,1,
     +                       iq,m,l,mobs,nm,nn)
                 if(tt.eq.1.0)then
                  do 605 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 605               continue
                 else
                  do 645 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 645               continue
               endif
             case(6)        
                 pi = proppr(s(nm),ismu,mu,cmu,tvmu,i,j,1,
     +                       iq,m,l,mobs,nm,nn)
                 if(tt.eq.1.0)then
                  do 606 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 606               continue
                 else
                  do 646 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 646               continue
               endif
             case(7)        
                 pi = binpr(s(nm),n(nm),tmp)
                 if(tt.eq.1.0)then
                  do 607 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 607               continue
                 else
                  do 647 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 647               continue
               endif
             case(8)        
                 pi = exppr(s(nm),tmp)
                 if(tt.eq.1.0)then
                  do 608 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 608               continue
                 else
                  do 648 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 648               continue
               endif
             case(9)        
                 pi = bbinpr(s(nm),n(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 609 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 609               continue
                 else
                  do 649 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 649               continue
               endif
             case(10)        
                 pi = nbinpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6010 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6010               continue
                 else
                  do 6410 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6410               continue
               endif
             case(11)        
                 pi = normpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6011 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6011               continue
                 else
                  do 6411 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6411               continue
               endif
             case(12)        
                 pi = invgpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6012 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6012               continue
                 else
                  do 6412 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6412               continue
               endif
             case(13)        
                 pi = logispr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6013 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6013               continue
                 else
                  do 6413 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6413               continue
               endif
             case(14)        
                 pi = cauchpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6014 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6014               continue
                 else
                  do 6414 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6414               continue
               endif
             case(15)        
                 pi = laplpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6015 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6015               continue
                 else
                  do 6415 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6415               continue
               endif
             case(16)        
                 pi = levypr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6016 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6016               continue
                 else
                  do 6416 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6416               continue
               endif
             case(17)        
                 pi = paretpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6017 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6017               continue
                 else
                  do 6417 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6417               continue
               endif
             case(18)        
                 pi = gammpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6018 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6018               continue
                 else
                  do 6418 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6418               continue
               endif
             case(19)        
                 pi = weibpr(s(nm),tmp,pshape(j))
                 if(tt.eq.1.0)then
                  do 6019 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6019              continue
                 else
                  do 6419 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6419               continue
               endif
             case(20)        
                 pi = ggampr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6020 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6020               continue
                 else
                  do 6420 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6420               continue
               endif
             case(21)        
                 pi = glogpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6021 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6021               continue
                 else
                  do 6421 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6421               continue
               endif
             case(22)        
                 pi = hjorpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6022 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6022               continue
                 else
                  do 6422 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6422               continue
               endif
             case(23)        
                 pi = burrpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6023 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6023               continue
                 else
                  do 6423 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6423               continue
               endif
             case(24)        
                 pi = gweipr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6024 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6024               continue
                 else
                  do 6424 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6424               continue
               endif
             case(25)        
                 pi = gextpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6025 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6025               continue
                 else
                  do 6425 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6425               continue
               endif
             case(26)        
                 pi = ginvgpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6026 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6026               continue
                 else
                  do 6426 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6426               continue
               endif
             case(27)        
                 pi = powexpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6027 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6027               continue
                 else
                  do 6427 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6427               continue
               endif
             case(28)        
                 pi = slaplpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6028 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6028               continue
                 else
                  do 6428 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6428               continue
               endif
             case(29)        
                 pi = studpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6029 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6029               continue
                 else
                  do 6429 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6429               continue
               endif
             case default
                 pi = studpr(s(nm),tmp,pshape(j),pfam)
                 if(tt.eq.1.0)then
                  do 6030 ii = 1, m
                     b(ii,j) = gamma(ii,j)+pi
 6030            continue
                 else
                  do 6430 ii = 1, m
                     b(ii,j) = gmod(ii,j)+pi
 6430            continue
               endif
             end select
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
 130     continue
      endif
      return
      end

      subroutine cfromx(x, m, gamma, pgamma)
*******************************************************************
*    Convert the vector of parameters into transition rate matrix *
*******************************************************************
      implicit none
      integer m,i,ii,j
      double precision sum,gamma(m,m),pgamma(m,m),x(*)

      ii=0
      do 30 i = 1, m
         sum = 0.
         do 20 j = 1, m
            if(j.ne.i.and.pgamma(i,j).ne.0.)then
               ii=ii+1
               gamma(i,j) = dexp(x(ii))
               sum = sum + gamma(i,j)
            endif
 20      continue
         gamma(i,i)=-sum
 30   continue
      return
      end

      subroutine geigen(gamma, val, vec, invec, a, c, gmod, pivot,
     +     qraux, work, work3, m)
***********************************************************
*    Obtain eigenvalues and vectors and inverse of former *
***********************************************************
      implicit none
      integer info,i,j, rank
c      double precision rank
      integer m,pivot(m)
      double precision val(m),vec(m,m),invec(m,m),gamma(m,m),a(m),c(m),
     +     gmod(m,m),qraux(m),work(2*m),work3(m,m)

c     get eigenvalues and vectors
      do 3 i=1,m
         do 4 j=1,m
            work3(i,j)=gamma(i,j)
 4       continue
 3    continue
      call rg(m,m,work3,val,a,1,vec,pivot,c,info)

c invert matrix of eigenvectors

      do 1 i=1,m
         do 2 j=1,m
            gmod(i,j)=vec(i,j)
            if(i.eq.j)then
               work3(i,i)=1.
            else
               work3(i,j)=0.
            endif
 2       continue
 1    continue
      call dqrdc2(gmod,m,m,m,1d-07,rank,qraux,pivot,work)
      call dqrcf(gmod,m,rank,qraux,work3,m,invec,info)
      return
      end

      subroutine mexp(gamma, val, vec, invec, time, m, tlog)
**********************************************************
*    Matrix exponentiation given eigenvalues and vectors *
**********************************************************
      implicit none
      integer i,j,k,m
      logical tlog
      double precision val(m),vec(m,m),invec(m,m),gamma(m,m),time

      do 3 i=1,m
         do 2 j=1,m
            gamma(i,j)=0.0
            do 1 k=1,m
               gamma(i,j)=gamma(i,j)+vec(i,k)*dexp(time*val(k))*
     +              invec(k,j)
 1          continue
            if(tlog)then
               gamma(i,j)=dlog(gamma(i,j))
            endif
 2       continue
 3    continue
      return
      end
