c
c  repeated : A Library of Repeated Measurements Models
c  Copyright (C) 1998 J.K. Lindsey
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
c    subroutine chidden(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,val,
c   +     vec,invec,model,cmu,tvmu,pshape,pfam,delta,nn,filter,cf,a,b,c,
c   +     gmod,rhs,pivot,qraux,work,like)
c
c  DESCRIPTION
c
c    Function to compute the likelihood of a hidden Markov chain model
c  with various response types in continuous time
c
      subroutine chidden(x,m,iq,nobs,mobs,s,n,times,l,pgamma,gamma,val,
     +     vec,invec,model,cmu,tvmu,pshape,pfam,delta,nn,filter,cf,a,b,
     +     c,gmod,rhs,pivot,qraux,work,like)
*************************************************************************
*     Function chidden computes minus the log likelihood of a           *
*     multivariate hidden Markov model with m states and iq individuals *
*     in continuous time.                                               *
*                                                                       *
*     Feel free to use or improve this program, provided that the       *
*     origin is acknowledged.                                           * 
*                              Iain MacDonald and Walter Zucchini       *
*     Modified by J.K. Lindsey for R, March, November, December 1998    *
*************************************************************************
      implicit none
      integer n(1),m,iq,i,j,k,l,model,nobs(1),mobs,ii,nm,i1,nn,
     +     rank,info,pivot(m)
      logical cf
      double precision like,s(1),pi,sflog,av,tt,pshape(m),pgamma(m,m),
     +     gmod(m,m),rhs(m),qraux(m),work(2*m),cmu(iq,m,l),
     +     tvmu(mobs,m,l),ll,x(1),gamma(m,m),delta(m),a(m),b(m,m),
     +     c(m),val(m),vec(m,m),invec(m,m),times(1),filter(m,nn),pfam
      double precision bernpr,poispr,multpr,binpr,exppr,bbinpr,nbinpr,
     +     normpr,invgpr,logispr,cauchpr,laplpr,levypr,paretpr,gammpr,
     +     weibpr,ggampr,glogpr,hjorpr,burrpr,gweipr,gextpr,ginvgpr,
     +     powexpr

      call cfromx(x,m,gamma,pgamma)

c     get eigenvalues and vectors
      call rg(m,m,gamma,val,a,1,vec,m,c,info)

c invert matrix of eigenvectors

      do 5 i=1,m
         do 4 j=1,m
            gmod(i,j)=vec(i,j)
            if(i.eq.j)then
               gamma(i,i)=1.
            else
               gamma(i,j)=0.
            endif
 4       continue
 5    continue
      call dqrdc2(gmod,m,m,m,1d-07,rank,qraux,pivot,work)
      call dqrcf(gmod,m,rank,qraux,gamma,m,invec,info,1)

c calculate gamma for unit time

      do 3 i=1,m
         do 2 j=1,m
            gamma(i,j)=0.0
            do 1 k=1,m
               gamma(i,j)=gamma(i,j)+vec(i,k)*dexp(val(k))*invec(k,j)
 1          continue
 2       continue
 3    continue

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
         nm=nm+1
         do 30 j = 1, m
            a(j)=delta(j)
            goto(201,202,203,204,205,206,207,208,209,210,211,212,213,
     +           214,215,216,217,218,219,220,221,222,223,224),model
 201        pi = bernpr(s(nm),cmu(i,j,1)+tvmu(1,j,1))
            goto 250
 202        pi = poispr(s(nm),cmu(i,j,1)+tvmu(1,j,1))
            goto 250
 203        pi = multpr(s(nm),cmu,tvmu,i,j,1,iq,m,l,mobs)
            goto 250
 204        pi = binpr(s(nm),n(nm),cmu(i,j,1)+tvmu(1,j,1))
            goto 250
 205        pi = exppr(s(nm),cmu(i,j,1)+tvmu(1,j,1))
            goto 250
 206        pi = bbinpr(s(nm),n(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 207        pi = nbinpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 208        pi = normpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 209        pi = invgpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 210        pi = logispr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 211        pi = cauchpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 212        pi = laplpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 213        pi = levypr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 214        pi = paretpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 215        pi = gammpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 216        pi = weibpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j))
            goto 250
 217        pi = ggampr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 218        pi = glogpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 219        pi = hjorpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 220        pi = burrpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 221        pi = gweipr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 222        pi = gextpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 223        pi = ginvgpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
            goto 250
 224        pi = powexpr(s(nm),cmu(i,j,1)+tvmu(1,j,1),pshape(j),pfam)
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
               do 63 j=1,m
                  do 62 ii = 1, m
                     gmod(ii,j)=0.0
                     do 61 i1=1,m
                        gmod(ii,j)=gmod(ii,j)+vec(ii,i1)*
     +                       dexp(tt*val(i1))*invec(i1,j)
 61                  continue
                     gmod(ii,j)=dlog(gmod(ii,j))
 62               continue
 63            continue
            endif
            do 70 j = 1, m
               goto(301,302,303,304,305,306,307,308,309,310,311,312,313,
     +              314,315,316,317,318,319,320,321,322,323,324),model
 301           pi = bernpr(s(nm),cmu(i,j,1)+tvmu(k,j,1))
               goto 350
 302           pi = poispr(s(nm),cmu(i,j,1)+tvmu(k,j,1))
               goto 350
 303           pi = multpr(s(nm),cmu,tvmu,i,j,k,iq,m,l,mobs)
               goto 350
 304           pi = binpr(s(nm),n(nm),cmu(i,j,1)+tvmu(k,j,1))
               goto 350
 305           pi = exppr(s(nm),cmu(i,j,1)+tvmu(k,j,1))
               goto 350
 306           pi = bbinpr(s(nm),n(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 307           pi = nbinpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 308           pi = normpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 309           pi = invgpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 310           pi = logispr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 311           pi = cauchpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 312           pi = laplpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 313           pi = levypr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 314           pi = paretpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 315           pi = gammpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 316           pi = weibpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j))
               goto 350
 317           pi = ggampr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 318           pi = glogpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 319           pi = hjorpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 320           pi = burrpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 321           pi = gweipr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 322           pi = gextpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 323           pi = ginvgpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
               goto 350
 324           pi = powexpr(s(nm),cmu(i,j,1)+tvmu(k,j,1),pshape(j),pfam)
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
********************************************************************
*    Converts the vector of parameters into transition rate matrix *
********************************************************************
      implicit none
      integer m,i,ii,j
      double precision sum,gamma(m,m),pgamma(m,m),x(1)

      ii=0
      do 30 i = 1, m
         sum = 0.
         do 20 j = 1, m
            if(j.ne.i.and.pgamma(i,j).ne.0.)then
               ii=ii+1
               gamma(i,j) = x(ii)
               sum = sum + gamma(i,j)
            endif
 20      continue
         gamma(i,i)=-sum
 30   continue
      return
      end
