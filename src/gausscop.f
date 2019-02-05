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
c     subroutine gcopula_f(theta,like,x,y,nobs,nest,lnest,nind,
c    +     nld,npre,npar,ar,v,tmp1,tmp2,tmp3,warn)
c
c  DESCRIPTION
c
c    Function to compute the likelihood function for the Gaussian
c copula with various autocorrelation functions and
c one or two levels of random effects
c
c
      subroutine gcopula_f(theta,like,x,y,nobs,nest,lnest,nind,
     +     nld,npre,npar,ar,v,tmp1,tmp2,tmp3,warn)
c
c routine for computing -log(l),
c including both autoregression and (nested) random effects
c
      implicit none
      integer nind,i,j,k,ar,lnest,nm,npre,npar,nld
      integer nobs(nind),nest(*)
      double precision x(*),y(*),like,theta(*),v(nld,nld),
     +     tausq(2),rho,ldet,det(2),tmp1(nld),tmp2(nld),tmp3(nld)
      logical warn
      warn=.false.
c
c     compute parameters of the correlation matrix
c
      if(npre.gt.0)then
         do 3 i=1,npre
            tausq(i)=theta(i)
 3       continue
         if(npre.eq.1)tausq(2)=0.
      else
         tausq(1)=0.
         tausq(2)=0.
      endif
      if(npar.gt.0)then
         rho=theta(npre+1)
         if(rho.eq.1.)rho=0.9999
      else
         rho=0.
      endif
c
c     compute the likelihood components
c
      like=0.d0
      nm=0
      do 49 i=1,nind
c
c     compute the inverse of correlation matrix and its log determinant
c
         call cmpcorr(v,ldet,det,tausq,rho,nind,i,nm,x,nobs(i),nest,
     +        lnest,nld,npre,npar,ar,warn,tmp1,tmp2,tmp3)
c
c     compute the quadratic form
c
         do 29 k=1,nobs(i)
            do 39 j=1,nobs(i)
               if(j.eq.k)v(j,j)=v(j,j)-1.
               ldet=ldet+y(nm+k)*v(k,j)*y(nm+j)
 39         continue
 29      continue
         like=like+ldet
         nm=nm+nobs(i)
 49   continue
c
c     calculate negative log likelihood for multivariate normal copula
c
      like=like/2.
      return
      end
c
c
      subroutine cmpcorr(v,ldet,det,tausq,rho,nind,i,nm,x,nobs,nest,
     +     lnest,nld,npre,npar,ar,warn,tmp1,tmp2,tmp3)
c
c compute the correlation matrix for the unit
c
      implicit none
      integer nobs,i,j,k,nm,nest(*),lnest,nn1,nn2,k1,j1,
     +     nind,nld,info,npre,npar,ar,ierr
      double precision x(*),v(nld,nld),det(2),tausq(2),
     +     rho,ldet,tmp,tmp1(nld),tmp2(nld),tmp3(nld)
      logical warn
c
      k1=0
      nn1=nest(nm+1)
      do 9 k=1,nobs
         if(lnest.gt.0)then
            if(nest(nm+k).ne.nn1)then
               k1=k1+1
               nn1=nest(nm+k)
            endif
         endif
         j1=0
         nn2=nest(nm+1)
         do 19 j=1,k
            if(lnest.gt.0)then
               if(nest(nm+j).ne.nn2)then
                  j1=j1+1
                  nn2=nest(nm+j)
               endif
            endif
            v(k,j)=tausq(1)
            if(k1.eq.j1)then
               v(k,j)=v(k,j)+tausq(2)
               if(k.eq.j)then
                  v(k,k)=1.0
               else if(rho.gt.0.0)then
                  if(ar.eq.1)then
                     tmp=rho**dabs(x(nm+k)-x(nm+j))
                  else if(ar.eq.2)then
                     tmp=rho**((x(nm+k)-x(nm+j))**2)
                  else if(ar.eq.3)then
                     tmp=1./(1+rho*(x(nm+k)-x(nm+j))**2)
                  else if(ar.eq.4)then
                     if(dabs(x(nm+k)-x(nm+j)).le.1/rho)then
                        tmp=((dabs(x(nm+k)-x(nm+j))*rho)**3-3*rho*
     +                    dabs(x(nm+k)-x(nm+j))+2)/2
                     else
                        tmp=0.
                     endif
                  endif
                  v(k,j)=v(k,j)+tmp
               endif
            endif
 19      continue
 9    continue
      do 1 k=2,nobs
         do 3 j=1,k-1
            v(j,k)=v(k,j)
 3       continue
 1    continue
c
c     check that correlation matrix is positive definite
c
      if(npre+npar.gt.1.and..not.warn)then
         call rs(nld,nobs,v,tmp1,0,0,tmp2,tmp3,ierr)
         tmp=1.
         do 2 j=1,nobs
            tmp=tmp*tmp1(j)
 2       continue
         warn=tmp.le.0.
      endif
c
c     factor v and compute the inverse and determinant
c
      call dpofa(v,nld,nobs,info)
      call dpodi(v,nld,nobs,det,11)
      ldet=dlog(det(1)*10.0**det(2))
      do 8 k=2,nobs
         do 18 j=1,k-1
            v(k,j)=v(j,k)
 18      continue
 8    continue
      return
      end
