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
c           subroutine binnest_f(F_value,res,x,g,hess,p,
c     1    num_case1,num_case2,num_subj,max_kid,max_mother,totcol,
c     2    total1,total2,total3,uph1_in,uph2_in,fcalls,dim_w,par,
c     3    Case1,Case2,subject,iout,hab,had,hac,ha,
c     4    v1,v2,h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk,mother,
c     5    rr,r,sn,z,uu1,uu2,w,
c     6    habb,habs1,habs2,has1s1,has1s2,has2s2,ebb,ebs1,ebs2,fs1s1,
c     7    fs1s2,gs2s2,e2bb,e2bs1,e2bs2,f2s1s1,f2s1s2,g2s2s2)
c
c
c  DESCRIPTION
c
c   Fortran code to fit binary random effects models with two levels
c  of nesting.
c
c	Program Logit_bin_nest
      subroutine binnest_f(F_value,res,x,g,hess,p,
     1    num_case1,num_case2,num_subj,max_kid,max_mother,totcol,
     2    total1,total2,total3,uph1_in,uph2_in,fcalls,dim_w,par,
     3    Case1,Case2,subject,iout,hab,had,hac,ha,
     4    v1,v2,h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk,mother,
     5    rr,r,sn,z,uu1,uu2,w,
     6    habb,habs1,habs2,has1s1,has1s2,has2s2,ebb,ebs1,ebs2,fs1s1,
     7    fs1s2,gs2s2,e2bb,e2bs1,e2bs2,f2s1s1,f2s1s2,g2s2s2)

!===============================================================================
!   Program Name: 	Logit_bin_nest.f90   
!   Project Name:	Mixed Model
!   Author:		T.R. Ten Have
!                       modified for R by J.K. Lindsey, June 1999
!   Date:		01/05/97
!   Language:		Fortran 90
!   Platform:		Unix
!   Description:	
!	
!   Method:
!	
!   Inputs:		This program is reading a user-supplied file "user_x"
!			(locating in the same directory as this program) which
!			contains information about the name of the input datafile, 
!			value of epsilon, maximum number of function calls, number 
!			betas, number of sigmas, and initial values for those estimates.
!
!			***** User should supply the file "user_x" *****
!
!			- With all the information on the same line and separated
!			  by one space:
!			- Name of file      - Up to 20 characters
!			- Epsilon
!			- Maximum number of function calls
!			- Number of betas   - up to 10
!			- Number of sigmas  - up to 10
!			- Initial estimates.
!
!			Example:
!
!			xfile.dat 0.00001 300, 2 1 -1.4 1.4 0.9
!
!
!   Outputs:		The program will display the Vector of estimates (X),
!			Vector of gradients (G), Hessian matrix, inverse
!			Hessian matrix, and standard errors.
!
!
!==============================================================================

        ! Begin main program
   
        Implicit none

        integer max_t1,max_t2,max_t3

        parameter(max_t1=10,max_t2=10,max_t3=10)

	Integer total1, total2, total3, total

        integer res(3),i,totcol
        double precision par(3),p(total1+total2+total3)

	Integer FCALLS,Flag,Iter_N,Fun_N 
	Integer num_subj 
        Integer Method, iout

        Integer Num_Case1, Num_Case2
	Integer Uph1_in, Uph2_in        

c	Integer, allocatable :: Case1(:)
c	Integer, allocatable :: Case2(:)

	Integer Case1(Num_Case1)
	Integer Case2(Num_Case2)

	Integer Max_kid, Max_mother, Dim_W

c for calcbn subroutine
	Integer Mother(Max_mother)
	Double Precision HAB(Max_mother,total1)
	Double Precision HAC(Max_mother,total2)
	Double Precision HAD(Max_mother,total3)
	Double Precision HA(Max_mother)
	Double Precision v1(uph1_in),v2(uph2_in)
	Double Precision h1choo(uph1_in),h2choo(uph2_in)
	Double Precision Hn(uph2_in),H1(uph1_in),H2(uph2_in)
	Double Precision betakk(Max_mother,Max_Kid)  
	Double Precision sig1kk(Max_mother,Max_Kid)
	Double Precision sig2kk(Max_mother,Max_Kid)
	Double Precision Z(num_case1,Max_Mother,Max_Kid,Total1) 
	Double Precision Uu1(num_case1,Max_Mother,Max_Kid,Total2) 
	Double Precision Uu2(num_case1,Max_Mother,Max_Kid,Total3)

c for calcbn2 subroutine
	double precision RR(num_case1,Max_Mother,Max_kid)
	Double Precision R(num_case1,Max_Mother,Max_Kid)
	Double Precision Sn(num_case1,Max_Mother,Max_Kid)
	Double Precision Hess(total1+total2+total3,total1+total2+total3)
	Double Precision HABB(Max_mother,total1,total1),
     1   HABS1(Max_mother,total1,total2)
	Double Precision HABS2(max_mother,total1,total3),
     1   HAS1S1(max_mother,total2,total2)
	Double Precision HAS1S2(max_mother,total2,total3),
     1   HAS2S2(max_mother,total3,total3)
	Double Precision EBB(max_mother,total1,total1),
     1   EBS1(max_mother,total1,total2)
	Double Precision EBS2(max_mother,total1,total3),
     1   FS1S1(max_mother,total2,total2)
	Double Precision FS1S2(max_mother,total2,total3),
     1   GS2S2(max_mother,total3,total3)
	Double Precision E2BB(max_mother,max_mother,total1,total1),
     +  E2BS1(max_mother,max_mother,total1,total2)
	Double Precision E2BS2(max_mother,max_mother,total1,total3),
     +  F2S1S1(max_mother,max_mother,total2,total2)
	Double Precision F2S1S2(max_mother,max_mother,total2,total3),
     +  G2S2S2(max_mother,max_mother,total3,total3)

c for conminbn subroutine
	DOUBLE PRECISION w(dim_w)

	Double Precision q1_in, q2_in
c	double precision Accuracy, Dim_W
	double precision Accuracy

	Double precision EPS, F_value 
	
c	Double Precision, allocatable :: x(:)
c	Double Precision, allocatable :: G(:)
c	Double Precision, allocatable :: subject(:,:) 

	Double Precision x(max_t1+max_t2+max_t3)
	Double Precision G(max_t1+max_t2+max_t3)
	Double Precision subject(Num_Subj,totcol)

	! Reading information from the user supplied file 'users_file3'

        eps=par(1)
        q1_in=par(2)
        q2_in=par(3)
        do i=1,total1+total2+total3
           x(i)=p(i)
        enddo

        Total = Total1 + Total2 + Total3
c           Allocate (X(Total))
	
c	Allocate(G(total))

	! Initialize variables for Subroutine CONMINBN
	Accuracy = 10.E-20
	Method = 1

	! Call Subroutine Conminbn
	Call Conminbn(total1,total2,total3,num_case1,num_case2,
     1    num_subj,Case1,Case2,Subject,X,G,F_value,Iter_N,Fun_N,Flag,
     2    Eps,Fcalls,Accuracy,Method,Dim_W,Max_Kid,
     3    Max_Mother,Uph1_in,Uph2_in,q1_in,q2_in,
     4    mother,hab,had,hac,ha,v1,v2,h1choo,h2choo,hn,
     5    h1,h2,betakk,sig1kk,sig2kk,r,sn,z,uu1,uu2,w,iout,totcol)


	   ! Call Subroutine Calcbn2 for the second derivative.
        if(flag.eq.0)then
	   Call Calcbn2(Total1,Total2,Total3,num_case1,num_case2,
     1	     Num_subj,Case1,Case2,X,Subject,Max_Kid,Max_Mother,Uph1_in,
     2	     Uph2_in,q1_in,q2_in,mother,rr,r,sn,z,uu1,uu2,
     3       hess,habb,habs1,habs2,has1s1,has1s2,
     4       has2s2,ebb,ebs1,ebs2,fs1s1,fs1s2,gs2s2,e2bb,e2bs1,e2bs2,
     5       f2s1s1,f2s1s2,g2s2s2,hab,hac,had,ha,v1,v2,
     6       h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk,totcol)
	ENDIF

        res(1)=Iter_N
        res(2)=Fun_N
        res(3)=flag

	END !PROGRAM LOGIT_BIN_NEST ! End main program
!==========================================================================


	double precision function sum(x,n)
	integer i,n
	double precision x(n)
	sum=0.0
	do i=1,n
	   sum=sum+x(i)
	enddo
	return
	end

c      Double precision FUNCTION FACTOR(n1)	         
c available in logitord

!==========================================================================
!			Begin Subroutine Calcbn
!==========================================================================
	Subroutine Calcbn(T1,T2,T3,Num_F,Num_M,Num_S,Case1,Case2,Temp_X,
     1	  Gradient,Sli,R,Sn,Z,Uu1,Uu2,Max_K,Max_M,Uph1_Temp,
     2	  Uph2_Temp,q1,q2,mother,hab,had,hac,ha,v1,v2,h1choo,h2choo,hn,
     3    h1,h2,betakk,sig1kk,sig2kk)

	implicit none

        integer max_t1,max_t2,max_t3

        parameter(max_t1=10,max_t2=10,max_t3=10)

	Integer Max_k, Max_M

	Integer Total, T1, T2, T3, Num_F, Num_M,Num_S
	Integer i		!  index of Father
	Integer j,jj		!  index of mother
	Integer k		!  index of children
	Integer kk		!  index of estimate
	Integer lowh1,lowh2
	Integer Uph1_Temp,Uph2_Temp
	Integer Uph1,Uph2

	Integer Case2_count1, Case2_count2
 
c    Integer ,Allocatable :: Mother(:)

	Integer Mother(max_m)

c    One Dimensional Array
	Integer Case1(Num_F) 
	Integer Case2(Num_M)  

	Double Precision Uh1Fac, Uh2Fac
	Double Precision q1, q2,Li, sLi
	Double Precision DSQRT, Expon, FACTOR
	Double Precision P,Q,N

	Double Precision A_Prod,A_tmp, J_Prod, HA_tmp
	Double Precision B_Sum(max_T1), C_Sum(max_T2), D_Sum(max_T3)
	Double Precision E(max_T1), F(max_T2), G(max_T3)
	Double Precision D1_beta(max_T1), D1_Sig1(max_T2),
     1    D1_Sig2(max_T3)
	Double Precision HAj

c    Double Precision, Allocatable :: HAB(:,:)
c    Double Precision, Allocatable :: HAD(:,:)
c    Double Precision, Allocatable :: HAC(:,:)
c    Double Precision, Allocatable :: HA(:)
c    Double Precision, Allocatable :: v1(:),v2(:)
c    Double Precision, Allocatable :: h1choo(:),h2choo(:)
c    Double Precision, Allocatable :: Hn(:),H1(:),H2(:)
	
	Double Precision HAB(max_m,t1)
	Double Precision HAD(max_m,t3)
	Double Precision HAC(max_m,t2)
	Double Precision HA(max_m)
	Double Precision v1(uph1_temp),v2(uph2_temp)
	Double Precision h1choo(uph1_temp),h2choo(uph2_temp)
	Double Precision Hn(uph2_temp),H1(uph1_temp),H2(uph2_temp)

	Double Precision beta(max_T1),sigma1(max_T2),sigma2(max_T3) 
	Double Precision Temp_X(max_T1+max_T2+max_T3),
     1    Gradient(max_T1+max_T2+max_T3) 

c    Double Precision, Allocatable :: betakk(:,:)  
c    Double Precision, Allocatable :: sig1kk(:,:)
c    Double Precision, Allocatable :: sig2kk(:,:)
    
	Double Precision betakk(max_m,Max_K)  
	Double Precision sig1kk(max_m,Max_K)
	Double Precision sig2kk(max_m,Max_K)

	Double Precision  R(Num_F,Max_M,Max_K)  
	Double Precision  Sn(Num_F,Max_M,Max_K) 
	Double Precision  Z(Num_F,Max_M,Max_K,T1) 
	Double Precision  Uu1(Num_F,Max_M,Max_K,T2) 
	Double Precision  Uu2(Num_F,Max_M,Max_K,T3) 

	double precision sum

	Total = T1 + T2 + T3	! Number of estimate

C     Initialize Upper H and q and this value can change
	UpH1 = Uph1_Temp
	UpH2 = Uph2_Temp

c	beta = Temp_X(1:T1)	! Read Beta values from the vector x
c	sigma1 = Temp_X(1+T1:T1+T2) ! Read Sigma1 values from the vector x
c	sigma2 = Temp_X(1+T1+T2:Total) ! Read Sigma2 values from the vector x
	do i=1, t1
	   beta(i)=temp_x(i)
	enddo
	do i=1,t2
	   sigma1(i)=temp_x(t1+i)
	enddo
	do i=1,t3
	   sigma2(i)=temp_x(t1+t2+i)
	enddo

C     Allocate space for V and Hchoo arrays.
c    Allocate (v1(uph1))
c    Allocate (h1choo(uph1))
c    Allocate (v2(uph2))
c    Allocate (h2choo(uph2))
c    Allocate (Hn(uph2))
c    Allocate (H1(uph1))
c    Allocate (H2(uph2))
    
C     Initialize V and Hchoo arrays
	do i=1, uph1
	   v1(i) = 0.0
	   h1choo(i) = 0.0
	enddo
	do i=1, uph1
	   v2(i) = 0.0
	   h2choo(i) = 0.0
	   Hn(i) = 0.0
	enddo
    
C     Calculate V and Hchoo arrays
	Uph1 = Uph1 - 1
	Uph2 = Uph2 - 1
	Uh1Fac = Factor(Uph1)
	Uh2Fac = Factor(UpH2)
	Do LowH1 = 0, UpH1 
	   v1(lowh1+1) = (lowh1 - uph1 * q1)/(DSQRT(uph1*q1*(1-q1)))
	   h1choo(lowh1+1) = Uh1Fac/(Factor(lowh1) * Factor(uph1-lowh1))
	End Do
	Do lowh2 = 0, uph2
	   v2(lowh2+1) = (lowh2 - uph2 * q1)/(DSQRT(uph2*q2*(1-q2)))
	   h2choo(lowh2+1) = Uh2Fac/(Factor(lowh2) * Factor(uph2-lowh2))
	End do
	Do lowh2 = 0,Uph2
c  changed JKL
c	   Hn(lowh2+1) = h2choo(lowh2+1)*q2**lowh2*(1-q2)**(Uph2-lowh2)
	   Hn(lowh2+1) = dexp(dlog(h2choo(lowh2+1))+lowh2*dlog(q2)
     1          +(Uph2-lowh2)*dlog(1-q2))
	End do

C     Initialize variables for the main loop.
	J_Prod = 1.0
	A_Prod = 1.0
	A_tmp = 1.0
	do i=1,t1
	   B_Sum(i) = 0.0
	   E(i) = 0.0
	   D1_beta(i) = 0.0
	enddo
	do i=1,t2
	   C_Sum(i) = 0.0
	   F(i) = 0.0
	   D1_Sig1(i) = 0.0
	enddo
	do i=1,t3
	   D_Sum(i) = 0.0
	   G(i) = 0.0
	   D1_Sig2(i) = 0.0
	enddo
	HA_tmp = 1.0
	do i=1,total
	   Gradient(i) = 0.0
	enddo
	HAj = 1.0
	Li = 0.0
	sLi = 0.0

	Case2_count1 = 1
	Case2_count2 = 1
   
C     Begin main loop
	Do i = 1, Num_F

C       Allocate beta and sigma and initialize to zero
c      Allocate(betakk(Case1(i),Max_K))
c      Allocate(Sig1kk(Case1(i),Max_K))
c      Allocate(Sig2kk(Case1(i),Max_K))
	   do j=1,case1(i)
	      Mother(j) = 0
	      do k=1,max_k
		 betakk(j,k) = 0.0
		 sig1kk(j,k) = 0.0
		 sig2kk(j,k) = 0.0
	      enddo
	   enddo

c      allocate(Mother(Case1(i)))
c	   Mother = 0

C       Calculate beta and sigma
	   Do j = 1, Case1(i)	! Number of mothers with in a father
	      Mother(j) = Case2(Case2_count1)    
	      Do k = 1, Mother(j) ! Number of kids with in a mather
		 Do kk = 1, t1
		    betakk(j,k) = betakk(j,k) + Z(i,j,k,kk) * beta(kk)
		 End do
		 Do kk = 1, t2
		    sig1kk(j,k) = sig1kk(j,k) + Uu1(i,j,k,kk) * sigma1(kk)
		 End do
		 Do kk = 1, t3
		    sig2kk(j,k) = sig2kk(j,k) + Uu2(i,j,k,kk) * sigma2(kk)
		 End do
	      End do
	      case2_count1 = case2_count1 + 1 
	   End do
C       Allocate H A B C and D arrays and initialize to sero

c      Allocate(HAC(Case1(i),t2))
c      Allocate(HAD(Case1(i),t3))
c      Allocate(HA(Case1(i)))
c      Allocate(HAB(Case1(i),t1))

	   do j=1,case1(i)
	      HA(j) = 0.0
	      do k=1,t1
		 HAB(j,k) = 0.0
	      enddo
	      do k=1,t2
		 HAC(j,k) = 0.0
	      enddo
	      do k=1,t3
		 HAD(j,k) = 0.0
	      enddo
	   enddo

	   Do  LowH1 = 0, UpH1	! Loop to do the Summation of lowh1 to UpH1
	      Do j = 1, Case1(i) ! Loop to do the Product of J     
		 Do lowh2 = 0, uph2 ! Loop to do the Summation of lowh2 to UpH2
		    Do k = 1, Mother(j)	! Loop to do the Product of K
c  lines added JKL
		       Expon = betakk(j,k) + (sig1kk(j,k)*v1(lowH1+1))
     1		        +(sig2kk(j,k)*v2(lowh2+1))
                       if(expon.gt.35)then
                     if(betakk(j,k).gt.15.)betakk(j,k)=betakk(j,k)/2
                     if(sig1kk(j,k)*v1(lowH1+1).gt.15.)
     +                    sig1kk(j,k)=sig1kk(j,k)/dabs(2*v1(lowH1+1))
                     if(sig2kk(j,k)*v2(lowh2+1).gt.15.)
     +                    sig2kk(j,k)=sig2kk(j,k)/dabs(2*v2(lowh2+1))
		       Expon = betakk(j,k) + (sig1kk(j,k)*v1(lowH1+1))
     1		        +(sig2kk(j,k)*v2(lowh2+1))
c                         write(*,*)'expon > 35',betakk(j,k),sig1kk(j,k),
c     1                      v1(lowH1+1),sig2kk(j,k),v2(lowh2+1)
                       endif
                       expon=dexp(expon)
c		       Expon = DEXP(betakk(j,k) + (sig1kk(j,k)*v1(lowH1+1))
c     1		        +(sig2kk(j,k)*v2(lowh2+1)))
		       P  = Expon / (1 + Expon)
		       Q = 1 - P
		       N = R(i,j,k) - (Sn(i,j,k) * P)
c  changed JKL
c		       A_tmp = A_tmp *(P**(R(i,j,k)) * Q**(Sn(i,j,k)
c     1		         - R(i,j,k)))
		       A_tmp = dexp(dlog(A_tmp)+R(i,j,k)*dlog(P)+
     1                      (Sn(i,j,k)- R(i,j,k))*dlog(Q))

		  ! Summation of Bi, Ci, and Di for First derivative.  
		       Do kk = 1, t1
			  B_Sum(kk) = B_Sum(kk) + (Z(i,j,k,kk) * N)
		       End do
		       Do kk = 1, t2
			  C_Sum(kk) = C_Sum(kk) + (V1(lowh1+1) *
     1		            Uu1(i,j,k,kk) *  N)	
		       End do
		       Do kk = 1, t3
			  D_Sum(kk) = D_Sum(kk) + (V2(lowh2+1)
     1		           * Uu2(i,j,k,kk) *  N)	
		       End do
		    End do	! End Loop Product K


	       ! Calculate the Sum of HAB, HAC, HAD, and HA
		    Do kk = 1, T1
		       HAB(j,kk) = HAB(j,kk) + (Hn(lowh2+1) *
     1		         A_tmp * B_Sum(kk))
		    End Do
		    Do kk = 1, T2
		       HAC(j,kk) = HAC(j,kk) + (Hn(lowh2+1) *
     1		         A_tmp * C_Sum(kk))
		    End Do
		    Do kk = 1, T3
		       HAD(j,kk) = HAD(j,kk) + (Hn(lowh2+1) *
     1		         A_tmp * D_Sum(kk))
		    End Do 
		    HA(j) = HA(j) + (Hn(lowh2+1) * A_tmp)

	       ! Calculate Summation of lowh2 to UpH2
c  changed JKL
c		    H2(lowh2+1) = h2choo(lowh2+1) * q2**lowh2 *
c     1		      (1-q2)**(UpH2 - lowh2) * A_tmp
		    H2(lowh2+1) = dexp(dlog(h2choo(lowh2+1))+lowh2*dlog(q2)
     1		      +(UpH2-lowh2)*dlog(1-q2)+dlog(A_tmp))

	       ! Reset A,B,C, and D
		    A_tmp = 1.0
		    do kk=1,t1
		       B_Sum(kk) = 0.0
		    enddo
		    do kk=1,t2
		       C_Sum(kk) = 0.0
		    enddo
		    do kk=1,t3
		       D_Sum(kk) = 0.0
		    enddo
		 End do		! End Loop Sumamtion of lowh2 to UpH2
		 J_Prod = J_Prod * Sum(H2,uph2) ! Calculate the product of J
		 HAj =  Haj * ha(j) ! Calculate the Product of HA

	      End do		! End Loop of Product J

	 ! Calculate E, F, and G for first Derivative.
	      Do jj = 1, Case1(i)
		 Do kk = 1, T1
		    E(kk) = E(kk) + ((HAj/HA(jj)) * HAB(jj,kk))
		 End do
		 Do kk = 1, T2
		    F(kk) = F(kk) + ((HAj/HA(jj)) * HAC(jj,kk))
		 End do
		 Do kk = 1, T3
		    G(kk) = G(kk) + ((HAj/HA(jj)) * HAD(jj,kk))
		 End do
	      End do

	 ! Reset HA,HAB,HAC,HAD,HAj
	      do jj=1,case1(i)
		 HA(jj) = 0.0
		 do kk=1,t1
		    HAB(jj,kk) = 0.0
		 enddo
		 do kk=1,t2
		    HAC(jj,kk) = 0.0
		 enddo
		 do kk=1,t3
		    HAD(jj,kk) = 0.0
		 enddo
	      enddo
	      HAj = 1.0 

	 ! Calculate the Summation of lowh1 to UpH1
c  changed JKL
c	      H1(lowh1+1) = h1choo(lowh1+1)*q1**lowh1*(1-q1)**(UpH1-lowh1)
c     1	        * J_Prod
	      H1(lowh1+1) = dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     1             (UpH1-lowh1)*dlog(1-q1))*J_Prod
	      J_Prod = 1.0	! Reset product of J

	 ! Calculate the Summation of beta and sigma
	      Do kk = 1, t1
c  changed JKL
c		 D1_beta(kk) = D1_beta(kk) + (h1choo(lowh1+1)*q1**
c     1	           lowh1*(1-q1)**(UpH1-lowh1) * E(kk))
		 D1_beta(kk) = D1_beta(kk) + dexp(dlog(h1choo(lowh1+1))
     1             +lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1))*E(kk)
	      End do
	      Do kk = 1, t2
c  changed JKL
c		 D1_Sig1(kk) = D1_Sig1(kk) + (h1choo(lowh1+1)*q1**
c     1	            lowh1*(1-q1)**(UpH1-lowh1) * F(kk))
		 D1_Sig1(kk) = D1_Sig1(kk) + dexp(dlog(h1choo(lowh1+1))
     1	           +lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1))*F(kk)
	      End do
	      Do kk = 1, t3
c  changed JKL
c		 D1_Sig2(kk) = D1_Sig2(kk) + (h1choo(lowh1+1)*q1**
c     1	            lowh1*(1-q1)**(UpH1-lowh1) * G(kk))
		 D1_Sig2(kk) = D1_Sig2(kk) + dexp(dlog(h1choo(lowh1+1))
     1	            +lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1))*G(kk)
	      End do	 
	 ! Reset vectors E,F, and G
	      do kk=1,t1
		 E(kk) = 0.0
	      enddo
	      do kk=1,t2
		 F(kk) = 0.0
	      enddo
	      do kk=1,t3
		 G(kk) = 0.0
	      enddo

	   End do

c      deallocate(Mother)   !Deallocate a 

C       Release Unused space
c      Deallocate(HA)
c      Deallocate(HAB)
c      Deallocate(HAC)
c      Deallocate(HAD)

c      Deallocate(betakk)
c      Deallocate(Sig1kk)
c      Deallocate(Sig2kk)


	   Li = Sum(H1,uph1)		! Li = Suma(h1-H1)*[Prod(J)*{Suma(h2-H2)*(Prod(K)}]
	   sLi = sLi + LOG(Li)	! Summation of SLI


C       Calculate First Derivative
	   Do kk = 1, t1
	      Gradient(kk) = Gradient(kk) + (D1_beta(kk) / Li)
	   End do
	   Do kk = 1, t2
	      Gradient(t1+kk) = Gradient(t1+kk) + (D1_Sig1(kk) / Li)
	   End do
	   Do kk = 1, t3
	      Gradient(t1+t2+kk)=Gradient(t1+t2+kk)+(D1_Sig2(kk)/Li)
	   End do
      
C       Reset beta and sigma
	   do kk=1,t1
	      D1_beta(kk) = 0.0
	   enddo
	   do kk=1,t2
	      D1_Sig1(kk) = 0.0
	   enddo
	   do kk=1,t3
	      D1_Sig2(kk) = 0.0
	   enddo
	End do			! End Main Loop

c    DeAllocate (v1)
c    DeAllocate (h1choo)
c    DeAllocate (v2)
c    DeAllocate (h2choo)
c    DeAllocate (Hn)
c    DeAllocate (H1)
c    DeAllocate (H2)


	Sli = - Sli		! Take the '-' to change from Max to Min
	do kk=1,total
	   Gradient(kk) = - Gradient(kk) ! Take the '-' to change from Max to Min
	enddo
	Return			! Return Sli and Gradient to Conminbn

	End !Subroutine Calcbn	! End Subroutine Calcbn
C==========================================================================



C**************************************************************************
C			Begin Subroutine Calcbn2
C==========================================================================
      Subroutine Calcbn2(T1,T2,T3,Num_F,Num_M,Num_S,Case1,Case2,Temp_X,
     1   a,Max_K,Max_M,Uph1_Temp,Uph2_Temp,q1,q2,mother,rr,r,sn,z,uu1,
     2   uu2,hess,habb,habs1,
     7   habs2,has1s1,has1s2,has2s2,ebb,ebs1,ebs2,fs1s1,fs1s2,gs2s2,
     8   e2bb,e2bs1,e2bs2,f2s1s1,f2s1s2,g2s2s2,hab,hac,had,ha,v1,v2,
     9   h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk,totcol)

	implicit none	

        integer max_t1,max_t2,max_t3

        parameter(max_t1=10,max_t2=10,max_t3=10)

C     Constant declaration
	Integer Max_k, Max_M, totcol

C     Integer declaration
	Integer Total, T1, T2, T3, Num_F, Num_M,Num_S
	Integer i		!  index of Father
	Integer j,jj,jjj	!  index of mother
	Integer k		!  index of children
	Integer kk,kkk		!  index of estimate
	Integer UpH1, UpH2 
	Integer lowh1,lowh2
	Integer count
	Integer Uph1_Temp,Uph2_Temp

	Integer Case2_count1
	Integer Case2_count2
	Integer Case2_count3  
    
C    One Dimensional Array
	Integer Case1(Num_F)
	Integer Case2(Num_M)  

c    Integer ,Allocatable :: Mother(:)

	Integer Mother(max_m)

C     Double Precision declaration
	Double Precision Uh1Fac, Uh2Fac
	Double Precision q1, q2,Li
	Double Precision DSQRT, DEXP, Expon, FACTOR
	Double Precision P,Q,N,O

	Double Precision  A(Num_S,totcol) 

	double precision RR(Num_F,Max_M,Max_k)

	Double Precision  R(Num_F,Max_M,Max_K)  
	Double Precision  Sn(Num_F,Max_M,Max_K) 
	Double Precision  Z(Num_F,Max_M,Max_K,T1) 
	Double Precision  Uu1(Num_F,Max_M,Max_K,T2) 
	Double Precision  Uu2(Num_F,Max_M,Max_K,T3) 

	Double Precision A_Prod,A_tmp, J_Prod, HA_tmp
	Double Precision B_Sum(MAX_T1), C_Sum(MAX_T2), D_Sum(MAX_T3)
	Double Precision E(MAX_T1), F(MAX_T2), G(MAX_T3)
	Double Precision D1_beta(MAX_T1),D1_Sig1(MAX_T2),D1_Sig2(MAX_T3)
	Double Precision HAj

	Double Precision D2_BB(MAX_T1,MAX_T1),D2_BS1(MAX_T1,max_t2),
     1    D2_BS2(MAX_T1,MAX_T3)
	Double Precision D2_S1S1(MAX_T2,max_t2),D2_S1S2(max_t2,max_t3),
     1    D2_S2S2(max_t3,max_t3)

	Double Precision EEBB(MAX_T1,MAX_T1),EEBS1(MAX_T1,MAX_T2),
     1    EEBS2(MAX_T1,MAX_T3)
	Double Precision FFS1S1(MAX_T2,MAX_T2),FFS1S2(MAX_T2,MAX_T3),
     1    GGS2S2(MAX_T3,MAX_T3)

	Double Precision BB(MAX_T1,MAX_T1),BS1(MAX_T1,MAX_T2),
     1    BS2(MAX_T1,MAX_T3)
	Double Precision S1S1(MAX_T2,MAX_T2), S1S2(MAX_T2,MAX_T3),
     1    S2S2(MAX_T3,MAX_T3)

	Double Precision Beta_beta(max_t1,max_t1),
     1    beta_sig1(max_t1,max_t2),beta_sig2(max_t1,max_t3)
	Double Precision sig1_sig1(max_t2,max_t2),
     1    sig1_sig2(max_t2,max_t3),sig2_sig2(max_t3,max_t3)

	Double Precision D2_beta(max_t1),D2_sig1(MAX_T2), D2_sig2(MAX_T3)

	Double Precision Hess(t1+t2+t3,t1+t2+t3)

C     Double Precision allocatable array
c    Double Precision, Allocatable :: HABB(:,:,:),HABS1(:,:,:)
c    Double Precision, Allocatable :: HABS2(:,:,:),HAS1S1(:,:,:)
c    Double Precision, Allocatable :: HAS1S2(:,:,:),HAS2S2(:,:,:)

c    Double Precision, Allocatable :: EBB(:,:,:),EBS1(:,:,:)
c    Double Precision, Allocatable :: EBS2(:,:,:),FS1S1(:,:,:)
c    Double Precision, Allocatable :: FS1S2(:,:,:),GS2S2(:,:,:)

c    Double Precision, Allocatable :: E2BB(:,:,:,:),E2BS1(:,:,:,:)
c    Double Precision, Allocatable :: E2BS2(:,:,:,:),F2S1S1(:,:,:,:)
c    Double Precision, Allocatable :: F2S1S2(:,:,:,:),G2S2S2(:,:,:,:)

c    Double Precision, Allocatable :: HAB(:,:),HAC(:,:),HAD(:,:)
c    Double Precision, Allocatable :: HA(:)
c    Double Precision, Allocatable :: v1(:),v2(:)
c    Double Precision, Allocatable :: h1choo(:),h2choo(:)
c    Double Precision, Allocatable :: Hn(:),H1(:),H2(:)

	Double Precision HABB(max_m,t1,t1),HABS1(max_m,t1,t2)
	Double Precision HABS2(max_m,t1,t3),HAS1S1(max_m,t2,t2)
	Double Precision HAS1S2(max_m,t2,t3),HAS2S2(max_m,t3,t3)

	Double Precision EBB(max_m,t1,t1),EBS1(max_m,t1,t2)
	Double Precision EBS2(max_m,t1,t3),FS1S1(max_m,t2,t2)
	Double Precision FS1S2(max_m,t2,t3),GS2S2(max_m,t3,t3)

	Double Precision E2BB(max_m,max_m,t1,t1),
     +  E2BS1(max_m,max_m,t1,t2)
	Double Precision E2BS2(max_m,max_m,t1,t3),
     +  F2S1S1(max_m,max_m,t2,t2)
	Double Precision F2S1S2(max_m,max_m,t2,t3),
     +  G2S2S2(max_m,max_m,t3,t3)

	Double Precision HAB(max_m,t1),HAC(max_m,t2),
     +  HAD(max_m,t3)
	Double Precision HA(max_m)
	Double Precision v1(uph1_temp),v2(uph2_temp)
	Double Precision h1choo(uph1_temp),h2choo(uph2_temp)
	Double Precision Hn(uph2_temp),H1(uph1_temp),H2(uph2_temp)

	Double Precision beta(MAX_T1),sigma1(MAX_T2),sigma2(MAX_T3)
	Double Precision Temp_X(MAX_T1+MAX_T2+MAX_T3)

c    Double Precision, Allocatable :: betakk(:,:)  
c    Double Precision, Allocatable :: sig1kk(:,:)
c    Double Precision, Allocatable :: sig2kk(:,:)

	Double Precision betakk(max_m,Max_k)
	Double Precision sig1kk(max_m,Max_k)
	Double Precision sig2kk(max_m,Max_k)

	double precision sum

	Total = T1 + T2 + T3	! Number of estimate

C     Initalize R, Sn, Z, Uu1, and Uu2
	count = 0.0
	Case2_count1 = 1
	Do i = 1, Num_F
	   Do j = 1, Case1(i)
	      Do k = 1, Case2(Case2_count1)
		 Count = count + 1
		 R(i,j,k) = A(count,1)
		 Sn(i,j,k) = A(count,2)
                 Z(i,j,k,1) = 1.0
		 Do kk = 2, t1
		    Z(i,j,k,kk) = A(count,kk+1)
		 End do
                 if(t2.gt.0)Uu1(i,j,k,1) = 1.0
		 Do kk = 2, t2
		    Uu1(i,j,k,kk) = A(count,t1+kk)
		 End do
                 if(t3.gt.0)Uu2(i,j,k,1) = 1.0
		 Do kk = 2, t3
		    Uu2(i,j,k,kk) = A(count,t1+t2+kk-1)
		 End do
	      End do
	      Case2_count1 = Case2_count1 + 1
	   End do
	End do
	Case2_count2 = 1
C     This loop use to elimilate a bug in the program
	Do i = 1, Num_F
	   Do j = 1, Case1(i)
	      Do k = 1, Case2(Case2_count2)
		 RR(i,j,k) = R(i,j,k)
	      End do
	      Case2_count2 = Case2_count2 + 1
	   End do
	End do

C     Initialize Upper H and q and this value can change
	UpH1 = Uph1_Temp
	UpH2 = Uph2_Temp

c	beta = Temp_X(1:T1)	! Read Beta values from the vector x
c	sigma1 = Temp_X(1+T1:T1+T2) ! Read Sigma1 values from the vector x
c	sigma2 = Temp_X(1+T1+T2:Total) ! Read Sigma2 values from the vector x

	do i=1,t1
	   beta(i)=temp_x(i)
	enddo
	do i=1,t2
	   sigma1(i)=temp_x(t1+i)
	enddo
	do i=1,t3
	   sigma2(i)=temp_x(t1+t2+i)
	enddo
C     Allocate space for V and Hchoo arrays.
c    Allocate (v1(uph1))
c    Allocate (h1choo(uph1))
c    Allocate (v2(uph2))
c    Allocate (h2choo(uph2))
c    Allocate (Hn(uph2))
c    Allocate (H1(uph1))
c    Allocate (H2(uph2))
    
C     Initialize V and Hchoo arrays
	do i=1,uph1
	   v1(i) = 0.0
	   h1choo(i) = 0.0
	enddo
	do i=1,uph2
	   v2(i) = 0.0
	   h2choo(i) = 0.0
	   Hn(i) = 0.0
	enddo
    
C     Calculate V and Hchoo arrays
	Uph1 = Uph1 - 1
	Uph2 = Uph2 - 1
	Uh1Fac = Factor(Uph1)
	Uh2Fac = Factor(UpH2)
	Do LowH1 = 0, UpH1 
	   v1(lowh1+1) = (lowh1 - uph1 * q1)/(DSQRT(uph1*q1*(1-q1)))
	   h1choo(lowh1+1) = Uh1Fac/(Factor(lowh1) * Factor(uph1-lowh1))
	End Do
	Do lowh2 = 0, uph2
	   v2(lowh2+1) = (lowh2 - uph2 * q1)/(DSQRT(uph2*q2*(1-q2)))
	   h2choo(lowh2+1) = Uh2Fac/(Factor(lowh2) * Factor(uph2-lowh2))
	End do
	Do lowh2 = 0,Uph2
c  changed JKL
c	   Hn(lowh2+1) = h2choo(lowh2+1)*q2**lowh2*(1-q2)**(Uph2-lowh2)
	   Hn(lowh2+1) = dexp(dlog(h2choo(lowh2+1))+lowh2*dlog(q2)+
     1          (Uph2-lowh2)*dlog(1-q2))
	End do

C     Initialize variables for the main loop.
	J_Prod = 1.0
	A_Prod = 1.0
	A_tmp = 1.0
	do i=1,t1
	   B_Sum(i) = 0.0
	enddo
	do i=1,t2
	   C_Sum(i) = 0.0
	enddo
	do i=1,t3
	   D_Sum(i) = 0.0
	enddo

	do i=1,t1
	   E(i) = 0.0
	   D1_beta(i) = 0.0
	   D2_beta(i) = 0.0
	   do j=1,t1
	      BB(i,j) = 0.0
	      EEBB(i,j) = 0.0
	      D2_BB(i,j) = 0.0
	      beta_beta(i,j) = 0.0
	   enddo
	   do j=1,t2
	      BS1(i,j) = 0.0
	      EEBS1(i,j) = 0.0
	      D2_BS1(i,j) = 0.0
	      beta_sig1(i,j) = 0.0
	   enddo
	   do j=1,t3
	      BS2(i,j) = 0.0
	      EEBS2(i,j) = 0.0
	      D2_BS2(i,j) = 0.0
	   enddo
	enddo
	do i=1,t2
	   F(i) = 0.0
	   D1_Sig1(i) = 0.0
	   D2_Sig1(i) = 0.0
	   do j=1,t2
	      S1S1(i,j) = 0.0
	      FFS1S1(i,j) = 0.0
	      D2_S1S1(i,j) = 0.0
	      sig1_sig1(i,j) = 0.0
	   enddo
	   do j=1,t3
	      S1S2(i,j) = 0.0
	      FFS1S2(i,j) = 0.0
	      D2_S1S2(i,j) = 0.0
	      beta_sig2(i,j) = 0.0
	      sig1_sig2(i,j) = 0.0
	   enddo
	enddo
	do i=1,t3
	   G(i) = 0.0
	   D1_Sig2(i) = 0.0
	   D2_Sig2(i) = 0.0
	   do j=1,t3
	      S2S2(i,j) = 0.0
	      GGS2S2(i,j) = 0.0
	      D2_S2S2(i,j) = 0.0
	      sig2_sig2(i,j) = 0.0
	   enddo
	enddo
	do i=1,total
	   do j=1,total
	      Hess(i,j) = 0.0
	   enddo
	enddo

	HA_tmp = 1.0
	HAj = 1.0
	Li = 0.0
 
	Case2_count3 = 1

C     Begin main loop
	Do i = 1, Num_F
C       Allocate beta and sigma and initialize to zero
c      Allocate(betakk(Case1(i),Max_k))
c      Allocate(Sig1kk(Case1(i),Max_k))
c      Allocate(Sig2kk(Case1(i),Max_k))
	   do j=1,case1(i)
	      do k=1,max_k
		 betakk(j,k) = 0.0
		 sig1kk(j,k) = 0.0
		 sig2kk(j,k) = 0.0
	      enddo
	   enddo
C       Calculate beta and sigma

c      Allocate(Mother(Case1(i)))
	   do j=1,case1(i)
	      Mother(j) = 0
	   enddo

	   Do j = 1, Case1(i)	! Number of mothers with in a father
	      Mother(j) = Case2(Case2_count3)
	      Do k = 1, Mother(j) ! Number of kids with in a mother
		 Do kk = 1, t1
		    betakk(j,k) = betakk(j,k) + Z(i,j,k,kk) * beta(kk)
		 End do
		 Do kk = 1, t2
		    sig1kk(j,k) = sig1kk(j,k) + Uu1(i,j,k,kk) * sigma1(kk)
		 End do
		 Do kk = 1, t3
		    sig2kk(j,k) = sig2kk(j,k) + Uu2(i,j,k,kk) * sigma2(kk)
		 End do
	      End do
	      Case2_count3 = Case2_count3 + 1
	   End do
C       Allocate H A B C and D arrays and initialize to sero
c      Allocate(HAB(Case1(i),t1))
c      Allocate(HAC(Case1(i),t2))
c      Allocate(HAD(Case1(i),t3))

c      Allocate(HA(Case1(i)))

c      Allocate(HABB(Case1(i),t1,t1))
c      Allocate(HABS1(Case1(i),t1,t2))
c      Allocate(HABS2(Case1(i),t1,t3))
c      Allocate(HAS1S1(Case1(i),t2,t2))
c      Allocate(HAS1S2(Case1(i),t2,t3))
c      Allocate(HAS2S2(Case1(i),t3,t3))

c      Allocate(EBB(Case1(i),t1,t1))
c      Allocate(EBS1(Case1(i),t1,t2))
c      Allocate(EBS2(Case1(i),t1,t3))
c      Allocate(FS1S1(Case1(i),t2,t2))
c      Allocate(FS1S2(Case1(i),t2,t3))
c      Allocate(GS2S2(Case1(i),t3,t3))

c      Allocate(E2BB(Case1(i),Case1(i),t1,t1))
c      Allocate(E2BS1(Case1(i),Case1(i),t1,t2))
c      Allocate(E2BS2(Case1(i),Case1(i),t1,t3))
c      Allocate(F2S1S1(Case1(i),Case1(i),t2,t2))
c      Allocate(F2S1S2(Case1(i),Case1(i),t2,t3))
c      Allocate(G2S2S2(Case1(i),Case1(i),t3,t3))

	   do j=1,case1(i)
	      HA(j) = 0.0
	      do k=1,t1
		 HAB(j,k) = 0.0
		 do kk=1,t1
		    HABB(j,k,kk) = 0.0
		    EBB(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       E2BB(j,jj,k,kk) = 0.0
		    enddo
		 enddo
		 do kk=1,t2
		    HABS1(j,k,kk) = 0.0
		    EBS1(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       E2BS1(j,jj,k,kk) = 0.0
		    enddo
		 enddo
		 do kk=1,t3
		    HABS2(j,k,kk) = 0.0
		    EBS2(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       E2BS2(j,jj,k,kk) = 0.0
		    enddo
		 enddo
	      enddo
	      do k=1,t2
		 HAC(j,k) = 0.0
		 do kk=1,t2
		    HAS1S1(j,k,kk) = 0.0
		    FS1S1(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       F2S1S1(j,jj,k,kk) = 0.0
		    enddo
		 enddo
		 do kk=1,t3
		    HAS1S2(j,k,kk) = 0.0
		    FS1S2(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       F2S1S2(j,jj,k,kk) = 0.0
		    enddo
		 enddo
	      enddo
	      do k=1,t3
		 HAD(j,k) = 0.0
		 do kk=1,t3
		    HAS2S2(j,k,kk) = 0.0
		    GS2S2(j,k,kk) = 0.0
		    do jj=1,case1(i)
		       G2S2S2(j,jj,k,kk) = 0.0
		    enddo
		 enddo
	      enddo
	   enddo

	   Do  LowH1 = 0, UpH1	! Loop to do the Summation of lowh1 to UpH1
	      Do j = 1, Case1(i) ! Loop to do the Product of J     
		 Do lowh2 = 0, uph2 ! Loop to do the Summation of lowh2 to UpH2
		    Do k = 1, Mother(j)	! Loop to do the Product of K
c  lines added JKL
		       Expon = betakk(j,k) + (sig1kk(j,k)*v1(lowH1+1))
     1		        +(sig2kk(j,k)*v2(lowh2+1))
                       if(expon.gt.35)then
                     if(betakk(j,k).gt.15.)betakk(j,k)=betakk(j,k)/2
                     if(sig1kk(j,k)*v1(lowH1+1).gt.15.)
     +                    sig1kk(j,k)=sig1kk(j,k)/dabs(2*v1(lowH1+1))
                     if(sig2kk(j,k)*v2(lowh2+1).gt.15.)
     +                    sig2kk(j,k)=sig2kk(j,k)/dabs(2*v2(lowh2+1))
		       Expon = betakk(j,k) + (sig1kk(j,k)*v1(lowH1+1))
     1		        +(sig2kk(j,k)*v2(lowh2+1))
c                         write(*,*)'expon > 35',betakk(j,k),sig1kk(j,k),
c     1                      v1(lowH1+1),sig2kk(j,k),v2(lowh2+1)
                       endif
                       expon=dexp(expon)
c		       Expon = DEXP(betakk(j,k) + (sig1kk(j,k) * v1(lowH1+1))
c     1		         +(sig2kk(j,k)*v2(lowh2+1)))
		       P  = Expon / (1 + Expon)
		       Q = 1 - P
		       N = RR(i,j,k) - (Sn(i,j,k) * P)
		       O = Sn(i,j,k) * P * Q
c  changed JKL
c		       A_tmp = A_tmp *(P**(RR(i,j,k)) * Q**(Sn(i,j,k)
c     1		          - RR(i,j,k)))
		       A_tmp = A_tmp *dexp(RR(i,j,k)*dlog(P)+
     1		          (Sn(i,j,k)- RR(i,j,k))*dlog(Q))
		       Do kk = 1, t1 ! Summation of Bi 
			  B_Sum(kk) = B_Sum(kk) + (Z(i,j,k,kk) * N)
			  Do kkk = 1, t1
			     BB(kk,kkk) = BB(kk,kkk)+(Z(i,j,k,kk)
     1			       *Z(i,j,k,kkk) * O)
			  End do
			  Do kkk = 1, t2
			     BS1(kk,kkk) = BS1(kk,kkk) + (Z(i,j,k,kk) *
     1			       V1(lowh1+1) * Uu1(i,j,k,kkk) * O)
			  End do
			  Do kkk = 1, t3
			     BS2(kk,kkk) = BS2(kk,kkk)+(Z(i,j,k,kk)*
     1			       V2(lowh2+1)*Uu2(i,j,k,kkk) * O)
			  End do
		       End do	! Summation of Ci
		       Do kk = 1, t2
			  C_Sum(kk) = C_Sum(kk) + (V1(lowh1+1) *
     1		            Uu1(i,j,k,kk) * N)
			  Do kkk = 1, t2
			     S1S1(kk,kkk) = S1S1(kk,kkk)+(V1(lowh1+1)*
     1			       V1(lowh1+1)
     2			       * Uu1(i,j,k,kk) * Uu1(i,j,k,kkk) * O)
			  End do
			  Do kkk = 1, t3
			     S1S2(kk,kkk) = S1S2(kk,kkk)+(V1(lowh1+1)*
     1			       V2(lowh2+1)
     2			       * Uu1(i,j,k,kk) * Uu2(i,j,k,kkk) * O)
			  End do  
		       End do
		       Do kk = 1, t3 ! Summation of Di
			  D_Sum(kk) = D_Sum(kk) + (V2(lowh2+1) *
     1		            Uu2(i,j,k,kk) *  N)	
			  Do kkk = 1, t3
			     S2S2(kk,kkk) = S2S2(kk,kkk)+(V2(lowh2+1)*
     1			       V2(lowh2+1)
     2			       * Uu2(i,j,k,kk) * Uu2(i,j,k,kkk) * O)	
			  End do
		       End do
		    End do	! End Loop Product K

		    Do kk = 1, T1 ! Sum HABi Loop
		       HAB(j,kk) = HAB(j,kk) + (Hn(lowh2+1) * A_tmp *
     1		          B_Sum(kk))
		       Do kkk = 1, T1		
			  HABB(j,kk,kkk) = HABB(j,kk,kkk) +
     1		            (Hn(lowh2+1)*A_tmp*
     2		            (B_Sum(kk) * B_Sum(kkk) - BB(kk,kkk)))
		       End do
		       Do kkk = 1, T2
			  HABS1(j,kk,kkk) = HABS1(j,kk,kkk)+(Hn(lowh2+1)*
     1		           A_tmp *
     2		           (B_Sum(kk) * C_Sum(kkk) - BS1(kk,kkk)))
		       End do
		       Do kkk = 1, T3
			  HABS2(j,kk,kkk) = HABS2(j,kk,kkk)+(Hn(lowh2+1)*
     1		             A_tmp *
     2		             (B_Sum(kk) * D_Sum(kkk) - BS2(kk,kkk)))
		       End do
		    End Do
		    Do kk = 1, T2 ! Sum HACi Loop
		       HAC(j,kk) = HAC(j,kk) + (Hn(lowh2+1) * A_tmp *
     1		         C_Sum(kk))
		       Do kkk = 1, T2
			  HAS1S1(j,kk,kkk) = HAS1S1(j,kk,kkk)+
     1		            (Hn(lowh2+1)*A_tmp*
     2		            (C_Sum(kk) * C_Sum(kkk) - S1S1(kk,kkk)))
		       End do
		       Do kkk = 1, T3
			  HAS1S2(j,kk,kkk) = HAS1S2(j,kk,kkk)+
     1		            (Hn(lowh2+1)*A_tmp*
     2		            (C_Sum(kk) * D_Sum(kkk) - S1S2(kk,kkk)))
		       End do 
		    End Do
		    Do kk = 1, T3 ! Sum HADi Loop
		       HAD(j,kk) = HAD(j,kk) + (Hn(lowh2+1) * A_tmp *
     1		         D_Sum(kk))
		       Do kkk = 1, T3
			  HAS2S2(j,kk,kkk) = HAS2S2(j,kk,kkk)+
     1		            (Hn(lowh2+1)*A_tmp*
     2		            (D_Sum(kk) * D_Sum(kkk)- S2S2(kk,kkk)))
		       End do
		    End Do 
		    HA(j) = HA(j) + (Hn(lowh2+1) * A_tmp)

	       ! Calculate Summation of lowh2 to UpH2
c  changed JKL
c		    H2(lowh2+1) = h2choo(lowh2+1) * q2**lowh2 *
c     1		      (1-q2)**(UpH2 - lowh2) * A_tmp
		    H2(lowh2+1)=dexp(dlog(h2choo(lowh2+1))+lowh2*dlog(q2)+
     1		      (UpH2 - lowh2)*dlog(1-q2)) * A_tmp
	       ! Reset A,B,C, and D
		    A_tmp = 1.0
		    do kk=1,t1
		       B_Sum(kk) = 0.0
		       do kkk=1,t1
			  BB(kk,kkk) = 0.0
		       enddo
		       do kkk=1,t2
			  BS1(kk,kkk) = 0.0
		       enddo
		       do kkk=1,t3
			  BS2(kk,kkk) = 0.0
		       enddo
		    enddo
		    do kk=1,t2
		       C_Sum(kk) = 0.0
		       do kkk=1,t2
			  S1S1(kk,kkk) = 0.0
		       enddo
		       do kkk=1,t3
			  S1S2(kk,kkk) = 0.0
		       enddo
		    enddo
		    do kk=1,t3
		       D_Sum(kk) = 0.0
		       do kkk=1,t3
			  S2S2(kk,kkk) = 0.0
		       enddo
		    enddo
		 End do		! End Loop Summation of lowh2 to UpH2
		 J_Prod = J_Prod * Sum(H2,uph2) ! Calculate the product of J
		 HAj =  Haj * ha(j) ! Calculate the Product of HA
	      End do		! End Loop of Product J

	 ! Calculate E, F, and G for first and second Derivative.
	      Do jj = 1, Case1(i)
		 Do kk = 1, T1
		    E(kk) = E(kk) + ((HAj/HA(jj)) * HAB(jj,kk))
		    Do kkk = 1, T1
		       EBB(jj,kk,kkk) = ((HAj/HA(jj))*HABB(jj,kk,kkk))
		    End do
		    Do kkk = 1, T2
		       EBS1(jj,kk,kkk) = ((HAj/HA(jj))*HABS1(jj,kk,kkk))
		    End do
		    Do kkk = 1, T3
		       EBS2(jj,kk,kkk) = ((HAj/HA(jj))*HABS2(jj,kk,kkk))
		    End do
		 End do
		 Do kk = 1, T2
		    F(kk) = F(kk) + ((HAj/HA(jj)) * HAC(jj,kk))
		    Do kkk = 1, T2
		       FS1S1(jj,kk,kkk)=((HAj/HA(jj))*HAS1S1(jj,kk,kkk)) 
		    End do
		    Do kkk = 1, T3
		       FS1S2(jj,kk,kkk) = ((HAj/HA(jj))*HAS1S2(jj,kk,kkk))
		    End do
		 End do
		 Do kk = 1, T3
		    G(kk) = G(kk) + ((HAj/HA(jj)) * HAD(jj,kk))
		    Do kkk = 1, T3
		       GS2S2(jj,kk,kkk)= ((HAj/HA(jj))*HAS2S2(jj,kk,kkk))
		    End do
		 End do
	      End do
C Bruce changes first index from `j` to `jj` given J loop ends line 1257
C	      do jj=1,case1(i)
C		 do kk=1,t1
C		    do kkk=1,t1
C		       HABB(jj,k,kk) = 0.0
C		    enddo
C		    do kkk=1,t2
C		       HABS1(jj,k,kk) = 0.0
C		    enddo
C		    do kkk=1,t3
C		       HABS2(jj,k,kk) = 0.0
C		    enddo
C		 enddo
C		 do kk=1,t2
C		    do kkk=1,t2
C		       HAS1S1(jj,k,kk) = 0.0
C		    enddo
C		    do kkk=1,t3
C		       HAS1S2(jj,k,kk) = 0.0
C		    enddo
C		 enddo
C		 do kk=1,t3
C		    do kkk=1,t3
C		       HAS2S2(jj,k,kk) = 0.0
C		    enddo
C		 enddo
C	      enddo

	      Do jj = 1, Case1(i) 
		 Do jjj = 1, Case1(i)
		    if (jj .ne. jjj) then
		       Do kk = 1, T1
			  Do kkk = 1, T1
			     E2BB(jj,jjj,kk,kkk) = E2BB(jj,jjj,kk,kkk)
     1			       +((HAj/(HA(jj)
     2			       * HA(jjj))) * HAB(jj,kk) * HAB(jjj,kkk))
			  End do
			  Do kkk = 1, T2
			     E2BS1(jj,jjj,kk,kkk)=E2BS1(jj,jjj,kk,kkk)
     1			       +((HAj/(HA(jj)
     2			       * HA(jjj))) * HAB(jj,kk) * HAC(jjj,kkk))
			  End do
			  Do kkk = 1, T3
			     E2BS2(jj,jjj,kk,kkk)=E2BS2(jj,jjj,kk,kkk)
     1			       +((HAj/(HA(jj)
     2			       * HA(jjj))) * HAB(jj,kk) * HAD(jjj,kkk))
			  End do
		       End do
		       Do kk = 1, T2
			  Do kkk = 1, T2
			     F2S1S1(jj,jjj,kk,kkk)=F2S1S1(jj,jjj,kk,kkk)
     1			       +((HAj/(HA(jj)
     2			       * HA(jjj))) * HAC(jj,kk) * HAC(jjj,kkk))
			  End do
			  Do kkk = 1, T3
			     F2S1S2(jj,jjj,kk,kkk)=F2S1S2(jj,jjj,kk,kkk)
     1			      +((HAj/(HA(jj)
     2			      * HA(jjj))) * HAC(jj,kk) * HAD(jjj,kkk))
			  End do
		       End do
		       Do kk = 1, T3
			  Do kkk = 1, T3
			     G2S2S2(jj,jjj,kk,kkk)=G2S2S2(jj,jjj,kk,kkk)
     1			       +((HAj/(HA(jj)
     2			       * HA(jjj))) * HAD(jj,kk) * HAD(jjj,kkk))
			  End do
		       End do
		    End if
		 End do
	      End do
		
	 ! Reset HA,HAB,HAC,HAD,HAj
	      do jj=1,case1(i)
		 HA(jj) = 0.0
		 do jjj=1,t1
		    HAB(jj,jjj) = 0.0
		 enddo
		 do jjj=1,t2
		    HAC(jj,jjj) = 0.0
		 enddo
		 do jjj=1,t3
		    HAD(jj,jjj) = 0.0
		 enddo
	      enddo
	      HAj = 1.0
     	       	  
	      Do jj = 1, Case1(i)
		 Do jjj = 1, Case1(i)
		    If (jj .ne. jjj) then	
		       Do kk = 1, T1
			  Do kkk = 1, T1
			     EEBB(kk,kkk)=EEBB(kk,kkk)+E2BB(jj,jjj,kk,kkk)
			  End do
			  Do kkk = 1, T2
			     EEBS1(kk,kkk)=EEBS1(kk,kkk)
     1			       + E2BS1(jj,jjj,kk,kkk)
			  End do
			  Do kkk = 1, T3
			     EEBS2(kk,kkk) = EEBS2(kk,kkk)
     1			       + E2BS2(jj,jjj,kk,kkk)
			  End do
		       End do
		       Do kk = 1, T2
			  Do kkk = 1, T2
			     FFS1S1(kk,kkk) = FFS1S1(kk,kkk)
     1			       + F2S1S1(jj,jjj,kk,kkk)
			  End do
			  Do kkk = 1, T3
			     FFS1S2(kk,kkk) = FFS1S2(kk,kkk)
     1			       + F2S1S2(jj,jjj,kk,kkk)
			  End do
		       End do
		       Do kk = 1, T3
			  Do kkk = 1, T3
			     GGS2S2(kk,kkk) = GGS2S2(kk,kkk)
     1			      + G2S2S2(jj,jjj,kk,kkk)
			  End do
		       End do
		    End if
		 End do
		 Do kk = 1, T1
		    Do kkk = 1, T1
		       EEBB(kk,kkk) = EEBB(kk,kkk) + EBB(jj,kk,kkk)
		    End do
		    Do kkk = 1, T2
		       EEBS1(kk,kkk) = EEBS1(kk,kkk) + EBS1(jj,kk,kkk)
		    End do
		    Do kkk = 1, T3
		       EEBS2(kk,kkk) = EEBS2(kk,kkk) + EBS2(jj,kk,kkk)
		    End do
		 End do
		 Do kk = 1, T2
		    Do kkk = 1, T2
		       FFS1S1(kk,kkk) = FFS1S1(kk,kkk) + FS1S1(jj,kk,kkk)
		    End do
		    Do kkk = 1, T3
		       FFS1S2(kk,kkk) = FFS1S2(kk,kkk) + FS1S2(jj,kk,kkk)
		    End do
		 End do
		 Do kk = 1, T3
		    Do kkk = 1, T3
		       GGS2S2(kk,kkk) = GGS2S2(kk,kkk) + GS2S2(jj,kk,kkk)
		    End do
		 End do
	      End do

	      do jj=1,case1(i)
		 do kk=1,t1
		    do kkk=1,t1
		       EBB(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  E2BB(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		    do kkk=1,t2
		       EBS1(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  E2BS1(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		    do kkk=1,t3
		       EBS2(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  E2BS2(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		 enddo
		 do kk=1,t2
		    do kkk=1,t2
		       FS1S1(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  F2S1S1(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		    do kkk=1,t3
		       FS1S2(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  F2S1S2(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		 enddo
		 do kk=1,t3
		    do kkk=1,t3
		       GS2S2(jj,kk,kkk) = 0.0
		       do jjj=1,case1(i)
			  G2S2S2(jj,jjj,kk,kkk) = 0.0
		       enddo
		    enddo
		 enddo
	      enddo

	 ! Calculate the Summation of lowh1 to UpH1
c  changed JKL
c	      H1(lowh1+1) = h1choo(lowh1+1)*q1**lowh1*(1-q1)**(UpH1-lowh1)
c     1	        * J_Prod
	      H1(lowh1+1) = dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     1	        (UpH1-lowh1)*dlog(1-q1))* J_Prod
	
	      J_Prod = 1.0	! Reset product of J

	 ! Calculate the Summation of beta and sigma
	      Do kk = 1, t1
c  changed JKL
c		 D1_beta(kk) = D1_beta(kk) + (h1choo(lowh1+1)*q1**
c     1	           lowh1*(1-q1)**(UpH1-lowh1) * E(kk))
		 D1_beta(kk) = D1_beta(kk) + dexp(dlog(h1choo(lowh1+1))+
     1	           lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1)) * E(kk)
		 Do kkk = 1, T1
c  changed JKL
c		    D2_BB(kk,kkk) = D2_BB(kk,kkk) + (h1choo(lowh1+1)*q1**
c     1		     lowh1*(1-q1)**(UpH1-lowh1)*EEBB(kk,kkk))
		    D2_BB(kk,kkk) = D2_BB(kk,kkk) + 
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2                   (UpH1-lowh1)*dlog(1-q1))*EEBB(kk,kkk)
		 End do
		 Do kkk = 1, T2
c  changed JKL
c		    D2_BS1(kk,kkk)=D2_BS1(kk,kkk)+ (h1choo(lowh1+1)*q1**
c     1		       lowh1*(1-q1)**(UpH1-lowh1)*EEBS1(kk,kkk))
		    D2_BS1(kk,kkk)=D2_BS1(kk,kkk)+
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2                   (UpH1-lowh1)*dlog(1-q1))*EEBS1(kk,kkk)
		 End do
		 Do kkk = 1, T3
c  changed JKL
c		    D2_BS2(kk,kkk)=D2_BS2(kk,kkk)+ (h1choo(lowh1+1)*q1**
c     1		      lowh1*(1-q1)**(UpH1-lowh1)*EEBS2(kk,kkk))
		    D2_BS2(kk,kkk)=D2_BS2(kk,kkk)+
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2                   (UpH1-lowh1)*dlog(1-q1))*EEBS2(kk,kkk)
		 End do
	      End do
	      Do kk = 1, t2
c  changed JKL
c		 D1_Sig1(kk) = D1_Sig1(kk) + (h1choo(lowh1+1)*q1**
c     1	           lowh1*(1-q1)**(UpH1-lowh1) * F(kk))
		 D1_Sig1(kk) = D1_Sig1(kk)+dexp(dlog(h1choo(lowh1+1))+
     1	           lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1)) * F(kk)
		 Do kkk = 1, T2
c  changed JKL
c		    D2_S1S1(kk,kkk) = D2_S1S1(kk,kkk)+(h1choo(lowh1+1)*q1**
c     1		      lowh1*(1-q1)**(UpH1-lowh1)*FFS1S1(kk,kkk))
		    D2_S1S1(kk,kkk) = D2_S1S1(kk,kkk)+
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2                   (UpH1-lowh1)*dlog(1-q1))*FFS1S1(kk,kkk)
		 End Do
		 Do kkk = 1, T3
c  changed JKL
c		    D2_S1S2(kk,kkk) = D2_S1S2(kk,kkk)+(h1choo(lowh1+1)*q1**
c     1		      lowh1*(1-q1)**(UpH1-lowh1)*FFS1S2(kk,kkk))
		    D2_S1S2(kk,kkk) = D2_S1S2(kk,kkk)+
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2                   (UpH1-lowh1)*dlog(1-q1))*FFS1S2(kk,kkk)
		 End do
	      End do
	      Do kk = 1, t3
c  changed JKL
c		 D1_Sig2(kk) = D1_Sig2(kk) + (h1choo(lowh1+1)*q1**
c     1	           lowh1*(1-q1)**(UpH1-lowh1) * G(kk))
		 D1_Sig2(kk) = D1_Sig2(kk) + dexp(dlog(h1choo(lowh1+1))+
     1	           lowh1*dlog(q1)+(UpH1-lowh1)*dlog(1-q1)) * G(kk)
		 Do kkk = 1, T3
c  changed JKL
c		    D2_S2S2(kk,kkk) = D2_S2S2(kk,kkk)+(h1choo(lowh1+1)*q1**
c     1  		 lowh1*(1-q1)**(UpH1-lowh1)*GGS2S2(kk,kkk))
		    D2_S2S2(kk,kkk) = D2_S2S2(kk,kkk)+
     1                   dexp(dlog(h1choo(lowh1+1))+lowh1*dlog(q1)+
     2  		 (UpH1-lowh1)*dlog(1-q1))*GGS2S2(kk,kkk)
		 End do
	      End do
	 ! Reset vectors E,F, and G
	      do kk=1,t1
		 E(kk) = 0.0
		 do kkk=1,t1
		    EEBB(kk,kkk) = 0.0
		 enddo
		 do kkk=1,t2
		    EEBS1(kk,kkk) = 0.0
		 enddo
		 do kkk=1,t3
		    EEBS2(kk,kkk) = 0.0
		 enddo
	      enddo
	      do kk=1,t2
		 F(kk) = 0.0
		 do kkk=1,t2
		    FFS1S1(kk,kkk) = 0.0
		 enddo
		 do kkk=1,t3
		    FFS1S2(kk,kkk) = 0.0
		 enddo
	      enddo
	      do kk=1,t3
		 G(kk) = 0.0
		 do kkk=1,t3
		    GGS2S2(kk,kkk) = 0.0
		 enddo
	      enddo

	   End do

c      Deallocate(Mother)

C       Release Unused space
c      Deallocate(HA)
c      Deallocate(HAB)
c      Deallocate(HAC)
c      Deallocate(HAD)

c      Deallocate(HABB)
c      Deallocate(HABS1)
c      Deallocate(HABS2)
c      Deallocate(HAS1S1)
c      Deallocate(HAS1S2)
c      Deallocate(HAS2S2)

c      Deallocate(betakk)
c      Deallocate(Sig1kk)
c      Deallocate(Sig2kk)

c      Deallocate(EBB)
c      Deallocate(EBS1)
c      Deallocate(EBS2)
c      Deallocate(FS1S1)
c      Deallocate(FS1S2)
c      Deallocate(GS2S2)
      
c      Deallocate(E2BB)
c      Deallocate(E2BS1)
c      Deallocate(E2BS2)
c      Deallocate(F2S1S1)
c      Deallocate(F2S1S2)
c      Deallocate(G2S2S2)

	   Li = Sum(H1,uph1)		! Li = Suma(h1-H1)*[Prod(J)*{Suma(h2-H2)*(Prod(K)}]

C       Calculate First Derivative
	   Do kk = 1, t1
	      D2_beta(kk) = (D1_beta(kk) / Li)
	   End do
	   Do kk = 1, t2
	      D2_Sig1(kk) = (D1_Sig1(kk) / Li)
	   End do
	   Do kk = 1, t3
	      D2_Sig2(kk) = (D1_Sig2(kk) / Li)
	   End do
      
	   Do kk = 1, T1
	      Do kkk = 1, T1
		 beta_beta(kk,kkk) = D2_beta(kk) * D2_beta(kkk)
	      End do
	      Do kkk = 1, T2
		 Beta_Sig1(kk,kkk) = D2_beta(kk) * D2_Sig1(kkk)
	      End do
	      Do kkk = 1, T3
		 Beta_Sig2(kk,kkk) = D2_beta(kk) * D2_Sig2(kkk)
	      End do
	   End do
	   Do kk = 1, T2
	      Do kkk = 1, T2
		 Sig1_Sig1(kk,kkk) = D2_Sig1(kk) * D2_Sig1(kkk)
	      End do
	      Do kkk = 1, T3
		 Sig1_Sig2(kk,kkk) = D2_Sig1(kk) * D2_Sig2(kkk)
	      End do
	   End do
	   Do kk = 1, T3
	      Do kkk = 1, T3
		 Sig2_Sig2(kk,kkk) = D2_Sig2(kk) * D2_Sig2(kkk)
	      End do
	   End do 
      
	   Do kk = 1, T1
	      Do kkk = 1, T1	! Derivatives with respect to beta_beta
		 Hess(kk,kkk) = Hess(kk,kkk) +
     1	          ((D2_BB(kk,kkk) / Li) - beta_beta(kk,kkk))
	      End do
	      Do kkk = 1, T2	! Derivatives with respect to beta_Sigma1
		 Hess(kk,kkk+T1) = Hess(kk,kkk+t1) +
     1	          ((D2_BS1(kk,kkk) / Li) - beta_Sig1(kk,kkk)) 
	      End do
	      Do kkk = 1, T3	! Derivatives with respect to beta_Sigma2
		 Hess(kk,kkk+T1+T2)= Hess(kk,kkk+t1+t2) +
     1	           ((D2_BS2(kk,kkk) / Li) - beta_Sig2(kk,kkk)) 
	      End do
	   End do
	   Do kk = 1, T2
	      Do kkk = 1, T2	! Derivatives with respect to Sigma1_Sigma1
		 Hess(kk+T1,kkk+T1) = Hess(kk+t1,kkk+T1) +
     1	           ((D2_S1S1(kk,kkk) / Li) - Sig1_Sig1(kk,kkk)) 
	      End do
	      Do kkk = 1, T3	! Derivatives with respect to Sigma1_Sigma2
		 Hess(kk+T1,kkk+T1+T2) = Hess(kk+t1,kkk+t1+T2) +
     1	           ((D2_S1S2(kk,kkk) / Li) - Sig1_Sig2(kk,kkk))
	      End do
	   End do
	   Do kk = 1, T3
	      Do kkk = 1, T3	! Derivatives with respect to sigma2_sigma2
		 Hess(kk+T1+T2,kkk+T1+T2) = Hess(kk+t1+t2,kkk+T1+T2) +
     1	           ((D2_S2S2(kk,kkk) / Li) - Sig2_Sig2(kk,kkk))
	      End do
	   End do

C       Reset beta and sigma
	   do kk=1,t1
	      D1_beta(kk) = 0.0
	      D2_beta(kk) = 0.0
	      do kkk=1,t1
		 D2_BB(kk,kkk) = 0.0
	      enddo
	      do kkk=1,t2
		 D2_BS1(kk,kkk) = 0.0
	      enddo
	      do kkk=1,t3
		 D2_BS2(kk,kkk) = 0.0
	      enddo
	   enddo
	   do kk=1,t2
	      D1_Sig1(kk) = 0.0
	      D2_Sig1(kk) = 0.0
	      do kkk=1,t2
		 D2_S1S1(kk,kkk) = 0.0
	      enddo
	      do kkk=1,t3
		 D2_S1S2(kk,kkk) = 0.0
	      enddo
	   enddo
	   do kk=1,t3
	      D1_Sig2(kk) = 0.0
	      D2_Sig2(kk) = 0.0
	      do kkk=1,t3
		 D2_S2S2(kk,kkk) = 0.0
	      enddo
	   enddo

	End do			! End Main Loop

	Do kk = 1, total
	   Do kkk = kk+1, total
		 Hess(kkk,kk) = Hess(kk,kkk)
	   End do
	End do

	Return			! Return Sli and Gradient to Conminbn
	End !Subroutine Calcbn2	! End Subroutine Calcbn

C*************************************************************************
		      ! Begin Subroutine CONMINBN	
C=========================================================================

	Subroutine Conminbn(total1,total2,total3,Father,mother,kid,
     1    Case1,Case2,A,X,G,F,iter,ifun,nflag,
     2	  Eps,MxFun,Acc,Nmeth,Mdim,Max_K1,
     3    Max_M1,Uph1_in,Uph2_in,q1_in,q2_in,mother1,hab,had,
     4    hac,ha,v1,v2,h1choo,h2choo,hn,h1,h2,betakk,
     5    sig1kk,sig2kk,r,sn,z,uu1,uu2,w,iout,totcol)

	implicit none	

        integer max_t1,max_t2,max_t3

        parameter(max_t1=10,max_t2=10,max_t3=10)

	Integer Max_K1, Max_M1

	LOGICAL RSW

	Integer total1,total2,total3,Father,Mother,kid
	Integer N,i,j,k,ii, ij,kk,count,totcol

	INTEGER ifun, iter
	INTEGER nflag, nmeth, nx, ng, ncons, nry, ncons1, ncons2
	INTEGER nrst, ncalls, nxpi, ngpi, nrdpi, nrypi, ngpj, nrd
	INTEGER mxfun, mdim
	Integer Case2_count1
	Integer Uph1_in,Uph2_in
        integer iout,ioutk,idev

	Integer Case1(Father)
	Integer Case2(Mother)

c for calcbn subroutine
	Integer Mother1(max_m1)
	Double Precision HAB(max_m1,total1)
	Double Precision HAD(max_m1,total3)
	Double Precision HAC(max_m1,total2)
	Double Precision HA(max_m1)
	Double Precision v1(uph1_in),v2(uph2_in)
	Double Precision h1choo(uph1_in),h2choo(uph2_in)
	Double Precision Hn(uph2_in),H1(uph1_in),H2(uph2_in)
	Double Precision betakk(max_m1,Max_K1)  
	Double Precision sig1kk(max_m1,Max_K1)
	Double Precision sig2kk(max_m1,Max_K1)

	Double precision X(max_t1+max_t2+max_t3)
	Double Precision G(max_t1+max_t2+max_t3)    
	Double Precision A(Kid,Totcol)
 	  
C     Double Precision Variables
	DOUBLE PRECISION F,FP,FMIN,ALPHA,AT,AP,GSQ,DG,DG1
	DOUBLE PRECISION DP,STEP,ACC,DAL,U1,U2,U3,U4,EPS
	DOUBLE PRECISION XSQ,RTST,DSQRT,DMIN1,DMAX1,DABS
c	DOUBLE PRECISION, DIMENSION (mdim) :: w
	DOUBLE PRECISION w(mdim)

	Double Precision  q1_in, q2_in
	Double Precision  R(Father,Max_M1,Max_K1)  
	Double Precision  Sn(Father,Max_M1,Max_K1) 
	Double Precision  Z(Father,Max_M1,Max_K1,Total1) 
	Double Precision  Uu1(Father,Max_M1,Max_K1,Total2) 
	Double Precision  Uu2(Father,Max_M1,Max_K1,Total3) 

C     End Declaration Section

        idev=6
	count = 0.0
	Case2_count1 = 1

	Do i = 1, Father
	   Do j = 1, Case1(i)
	      Do K = 1, Case2(Case2_count1)
		 Count = count + 1
		 R(i,j,k) = A(count,1)
		 Sn(i,j,k) = A(count,2)
                 Z(i,j,k,1) = 1.0
		 Do kk = 2, total1
		    Z(i,j,k,kk) = A(count,kk+1)
		 End do
                 if(total2.gt.0)Uu1(i,j,k,1) = 1.0
		 Do kk = 2, total2
		    Uu1(i,j,k,kk) = A(count,total1+kk)
		 End do
                 if(total3.gt.0)Uu2(i,j,k,1) = 1.0
		 Do kk = 2, total3
		    Uu2(i,j,k,kk) = A(count,total1+total2+kk-1)
		 End do
	      End do
	      Case2_count1 = Case2_count1 + 1
	   End do
	End do
 
	n=total1+total2+total3

C     Initialize ITER,IFUN,NFLAG, and IOUTK, which counts output iterations
	ITER=0
	IFUN=0
	IOUTK=0
	NFLAG=0
C     SET PARAMETERS TO EXTRACT VECTORS FROM W.
C     W(I) HOLDS THE SEARCH VECTOR,W(NX+I) HOLDS THE BEST CURRENT
C     ESTIMATE TO THE MINIMIZER,AND W(NG+I) HOLDS THE GRADIENT
C     AT THE BEST CURRENT ESTIMATE.

	NX=N
	NG=NX+N

C     TEST WHICH METHOD IS BEING USED.
C     IF NMETH=0, W(NRY+I) HOLDS THE RESTART Y VECTOR AND
C     W(NRD+I) HOLDS THE RESTART SEARCH VECTOR.

	IF (NMETH.EQ.1) THEN 
	   NCONS = 3 * N 	!If NMETH=1, W(NONS+I) holds the appr. inverse HESSIAN
	ELSE 
	   NRY=NG+N
	   NRD=NRY+N
	   NCONS=5*N
	   NCONS1=NCONS+1
	   NCONS2=NCONS+2
	END IF
 
C     CALCULATE THE FUNCTION AND GRADIENT AT THE INITIAL
C     POINT AND INITIALIZE NRST,WHICH IS USED TO DETERMINE
C     WHETHER A BEALE RESTART IS BEING DONE. NRST=N MEANS THAT THIS
C     ITERATION IS A RESTART ITERATION. INITIALIZE RSW,WHICH INDICATES
C     THAT THE CURRENT SEARCH DIRECTION IS A GRADIENT DIRECTION.

 20	Call  Calcbn(Total1,Total2,Total3,father,mother,kid,
     1    Case1,Case2,X,G,F,R,Sn,Z,Uu1,Uu2,Max_K1,
     2    Max_M1,Uph1_in,Uph2_in,q1_in,q2_in,mother1,hab,had,
     3    hac,ha,v1,v2,h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk)

	IFUN = IFUN+1
	NRST = N
	RSW = .TRUE.

C     CALCULATE THE INITIAL SEARCH DIRECTION , THE NORM OF X SQUARED,
C     AND THE NORM OF G SQUARED. DG1 IS THE CURRENT DIRECTIONAL
C     DERIVATIVE,WHILE XSQ AND GSQ ARE THE SQUARED NORMS.

	DG1 = 0.
	XSQ = 0.
	DO I = 1,N
	   W(I) = -G(I)
	   XSQ = XSQ + X(I) * X(I)
	   DG1 = DG1 - G(I) * G(I)
	END DO
	GSQ = -DG1

C     TEST IF THE INITIAL POINT IS THE MINIMIZER.
	IF (GSQ .le. EPS*EPS*DMAX1(1.0D0,XSQ)) THEN
	   RETURN
	END IF

C     BEGIN THE MAJOR ITERATION LOOP. NCALLS IS USED TO GUARANTEE THAT
C     AT LEAST TWO POINTS HAVE BEEN TRIED WHEN NMETH=0. FMIN IS THE
C     CURRENT FUNCTION VALUE.

 40	FMIN=F
	NCALLS=IFUN

C     IF OUTPUT IS DESIRED,TEST IF THIS IS THE CORRECT ITERATION
C     AND IF SO, WRITE OUTPUT.

	IF (IOUT .eq. 0) THEN
	   ALPHA = ALPHA * DG / DG1 ! Set ALPHA to nonrestart conjugate gadient
	ELSE IF (IOUTK .ne. 0) THEN
	   IOUTK = IOUTK + 1
	   IF (IOUTK .eq. IOUT) THEN
	      IOUTK = 0
	   END IF
	   ALPHA = ALPHA * DG / DG1
	ELSE
C	   WRITE(IDEV,50)ITER,IFUN,FMIN,GSQ
 50	   FORMAT(10H ITERATION,I5,20H      FUNCTION CALLS,I6/5H F = ,
     1 D15.8,13H G-SQUARED = ,D15.8/)  
C	   WRITE(IDEV,60)(X(I),I=1,N)
 60	   FORMAT(/8HINTER X./1H ,20D16.8)
	END IF

C     IF NMETH=1 OR A RESTART HAS BEEN PERFORMED, SET ALPHA=1.0.
	IF (NRST .eq. 1.OR.NMETH .eq. 1) THEN
	   ALPHA=1.0
	END IF

C     IF A GRADIENT DIRECTION IS USED, SET ALPHA=1.0/DSQRT(GSQ),
C     WHICH SCALES THE INITIAL SEARCH VECTOR TO UNITY.
	IF (RSW) THEN
	   ALPHA=1.0/DSQRT(GSQ)
	END IF

C     THE LINEAR SEARCH FITS A CUBIC TO F AND DAL, THE FUNCTION AND ITS
C     DERIVATIVE AT ALPHA, AND TO FP AND DP,THE FUNCTION
C     AND DERIVATIVE AT THE PREVIOUS TRIAL POINT AP.
C     INITIALIZE AP ,FP,AND DP.

	AP=0.
	FP=FMIN
	DP=DG1

C     SAVE THE CURRENT DERIVATIVE TO SCALE THE NEXT SEARCH VECTOR.
	DG=DG1

C     UPDATE THE ITERATION.
	ITER=ITER+1

C     CALCULATE THE CURRENT STEPLENGTH  AND STORE THE CURRENT X AND G.
	STEP=0.
	DO I=1,N
	   STEP=STEP+W(I)*W(I)
	   NXPI=NX+I
	   NGPI=NG+I
	   W(NXPI)=X(I)
	   W(NGPI)=G(I)
	END DO
	STEP=DSQRT(STEP)

C     BEGIN THE LINEAR SEARCH ITERATION.
C     TEST FOR FAILURE OF THE LINEAR SEARCH.

 80	IF (ALPHA*STEP .le. ACC) THEN
C     TEST IF DIRECTION IS A GRADIENT DIRECTION.
	   IF (.NOT.RSW) THEN
	      GO TO 20		! Call subroutine CALCBN           
	   ELSE
	      NFLAG=2
	      RETURN
	   END IF
	END IF

C     CALCULATE THE TRIAL POINT.
	DO I = 1,N
	   NXPI = NX + I
	   X(I) = W(NXPI) + ALPHA * W(I)
	END DO

C     EVALUATE THE FUNCTION AT THE TRIAL POINT.
C     Call CALCBN
	Call  Calcbn(Total1,Total2,Total3,father,mother,kid,
     1    Case1,Case2,X,G,F,R,Sn,Z,Uu1,Uu2,Max_K1,
     2    Max_M1,Uph1_in,Uph2_in,q1_in,q2_in,mother1,hab,had,
     3    hac,ha,v1,v2,h1choo,h2choo,hn,h1,h2,betakk,sig1kk,sig2kk)
 
C     TEST IF THE MAXIMUM NUMBER OF FUNCTION CALLS HAVE BEEN USED.
	IFUN=IFUN+1
	IF(IFUN.gt.MXFUN) THEN 
	   NFLAG=1
	   RETURN
	END IF

C     COMPUTE THE DERIVATIVE OF F AT ALPHA.
	DAL=0.0
	DO I=1,N
	   DAL=DAL+G(I)*W(I)
	END DO

C     TEST WHETHER THE NEW POINT HAS A NEGATIVE SLOPE BUT A HIGHER
C     FUNCTION VALUE THAN ALPHA=0. IF THIS IS THE CASE,THE SEARCH
C     HAS PASSED THROUGH A LOCAL MAX AND IS HEADING FOR A DISTANT LOCAL
C     MINIMUM.
	IF (F.gt.FMIN .AND. DAL.lt.0.) GO TO 160

C     IF NOT, TEST WHETHER THE STEPLENGTH CRITERIA HAVE BEEN MET.
	IF(F.gt.(FMIN+.0001*ALPHA*DG) .OR. DABS(DAL/DG)
     1   .gt.(.92)) GO TO 130

C     IF THEY HAVE BEEN MET, TEST IF TWO POINTS HAVE BEEN TRIED
C     IF NMETH=0 AND IF THE TRUE LINE MINIMUM HAS NOT BEEN FOUND.
	IF ((IFUN-NCALLS) .le. 1 .AND. DABS(DAL/DG).gt.EPS .AND.
     1    NMETH .eq. 0) THEN 
	   GO TO 130
	ELSE
	   GO TO 170
	END IF
C     A NEW POINT MUST BE TRIED. USE CUBIC INTERPOLATION TO FIND
C     THE TRIAL POINT AT.
 130	U1=DP+DAL-3.0*(FP-F)/(AP-ALPHA)
	U2=U1*U1-DP*DAL
	IF(U2.LT.0.)U2=0.
	U2=DSQRT(U2)
	AT=ALPHA-(ALPHA-AP)*(DAL+U2-U1)/(DAL-DP+2.*U2)

C     TEST WHETHER THE LINE MINIMUM HAS BEEN BRACKETED.
	IF((DAL/DP).GT.0.)GO TO 140

C     THE MINIMUM HAS BEEN BRACKETED. TEST WHETHER THE TRIAL POINT LIES
C     SUFFICIENTLY WITHIN THE BRACKETED INTERVAL.
C     IF IT DOES NOT, CHOOSE AT AS THE MIDPOINT OF THE INTERVAL.

	IF(AT.LT.(1.01*DMIN1(ALPHA,AP)).OR.AT.GT.(.99*DMAX1
     1    (ALPHA,AP)))AT=(ALPHA+AP)/2.0
	GO TO 150

C     THE MINIMUM HAS NOT BEEN BRACKETED. TEST IF BOTH POINTS ARE
C     GREATER THAN THE MINIMUM AND THE TRIAL POINT IS SUFFICIENTLY
C     SMALLER THAN EITHER.

 140	IF (DAL .GT.0.0.AND.0.0.LT.AT.AND.AT.LT.
     1  (.99*DMIN1(AP,ALPHA))) GO TO 150

C     TEST IF BOTH POINTS ARE LESS THAN THE MINIMUM AND THE TRIAL POINT
C     IS SUFFICIENTLY LARGE.
	IF(DAL.LE.0.0.AND.AT.GT.(1.01*DMAX1(AP,ALPHA)))GO TO 150

C     IF THE TRIAL POINT IS TOO SMALL,DOUBLE THE LARGEST PRIOR POINT.
	IF(DAL.LE.0.)AT=2.0*DMAX1(AP,ALPHA)

C     IF THE TRIAL POINT IS TOO LARGE, HALVE THE SMALLEST PRIOR POINT.
	IF(DAL.GT.0.)AT=DMIN1(AP,ALPHA)/2.0

C     SET AP=ALPHA, ALPHA=AT,AND CONTINUE SEARCH.
 150	AP=ALPHA
	FP=F
	DP=DAL
	ALPHA=AT
	GO TO 80

C     A RELATIVE MAX HAS BEEN PASSED.REDUCE ALPHA AND RESTART THE SEARCH.
 160	ALPHA=ALPHA/3.
	AP=0.
	FP=FMIN
	DP=DG
	GO TO 80

C     THE LINE SEARCH HAS CONVERGED. TEST FOR CONVERGENCE OF THE ALGORITHM.
 170	GSQ=0.0
	XSQ=0.0
	DO I=1,N
	   GSQ=GSQ+G(I)*G(I)
	   XSQ=XSQ+X(I)*X(I)
	END DO

c        write(*,*)'test',eps,xsq,gsq,EPS*EPS*DMAX1(1.0D0,XSQ)
	IF (GSQ .le. EPS*EPS*DMAX1(1.0D0,XSQ)) THEN
	   RETURN
	END IF

C     SEARCH CONTINUES. SET W(I)=ALPHA*W(I),THE FULL STEP VECTOR.
	DO I=1,N
	   W(I)=ALPHA*W(I)
	END DO

C     COMPUTE THE NEW SEARCH VECTOR. FIRST TEST WHETHER A
C     CONJUGATE GRADIENT OR A VARIABLE METRIC VECTOR IS USED.
	IF (NMETH .ne. 1) THEN	! Begin if nmeth /= 1
C        CONJUGATE GRADIENT UPDATE SECTION.
C        TEST IF A POWELL RESTART IS INDICATED.
	   RTST=0.
	   DO I=1,N
	      NGPI=NG+I
	      RTST=RTST+G(I)*W(NGPI)
	   END DO
	   IF (DABS(RTST/GSQ).gt.0.2) THEN
	      NRST=N
	   END IF
C        IF A RESTART IS INDICATED, SAVE THE CURRENT D AND Y
C        AS THE BEALE RESTART VECTORS AND SAVE D'Y AND Y'Y
C        IN W(NCONS+1) AND W(NCONS+2).
	   IF (NRST .eq. N) THEN 
	      W(NCONS1)=0.
	      W(NCONS2)=0.
	      DO I=1,N
		 NRDPI=NRD+I
		 NRYPI=NRY+I
		 NGPI=NG+I
		 W(NRYPI)=G(I)-W(NGPI)
		 W(NRDPI)=W(I)
		 W(NCONS1)=W(NCONS1)+W(NRYPI)*W(NRYPI)
		 W(NCONS2)=W(NCONS2)+W(I)*W(NRYPI)
	      END DO
	   END IF
C        CALCULATE  THE RESTART HESSIAN TIMES THE CURRENT GRADIENT.
	   U1=0.0
	   U2=0.0
	   DO I=1,N
	      NRDPI=NRD+I
	      NRYPI=NRY+I
	      U1=U1-W(NRDPI)*G(I)/W(NCONS1)
	      U2=U2+W(NRDPI)*G(I)*2./W(NCONS2)-
     1          W(NRYPI)*G(I)/W(NCONS1)
	   END DO
	   U3 = W(NCONS2)/W(NCONS1)
	   DO I=1,N
	      NXPI=NX+I
	      NRDPI=NRD+I
	      NRYPI=NRY+I
	      W(NXPI)=-U3*G(I)-U1*W(NRYPI)-U2*W(NRDPI)
	   END DO
C        IF THIS IS A RESTART ITERATION,W(NX+I) CONTAINS THE NEW SEARCH
C        VECTOR.
	   IF (NRST .ne. N) THEN	! begin if nrst /= n
	  ! NOT A RESTART ITERATION. CALCULATE THE RESTART HESSIAN
	  ! TIMES THE CURRENT Y.
	      U1=0.
	      U2=0.
	      U3=0.
	      U4=0.
	      DO I=1,N
		 NGPI=NG+I
		 NRDPI=NRD+I
		 NRYPI=NRY+I
		 U1=U1-(G(I)-W(NGPI))*W(NRDPI)/W(NCONS1)
		 U2=U2-(G(I)-W(NGPI))*W(NRYPI)/W(NCONS1)
     1             +2.0*W(NRDPI)*(G(I)-W(NGPI))/W(NCONS2)
		 U3=U3+W(I)*(G(I)-W(NGPI))
	      END DO
	      STEP=0.
	      DO I=1,N
		 NGPI=NG+I
		 NRDPI=NRD+I
		 NRYPI=NRY+I
		 STEP=(W(NCONS2)/W(NCONS1))*(G(I)-W(NGPI))
     1              +U1*W(NRYPI)+U2*W(NRDPI)
		 U4=U4+STEP*(G(I)-W(NGPI))
		 W(NGPI)=STEP
	      END DO

	  ! CALCULATE THE DOUBLY UPDATED HESSIAN TIMES THE CURRENT
	  ! GRADIENT TO OBTAIN THE SEARCH VECTOR.
	      U1=0.0
	      U2=0.0
	      DO I=1,N
		 U1=U1-W(I)*G(I)/U3
		 NGPI=NG+I
		 U2=U2+(1.0+U4/U3)*W(I)*G(I)/U3-W(NGPI)*G(I)/U3
	      END DO
	      DO I=1,N
		 NGPI=NG+I
		 NXPI=NX+I
		 W(NXPI)=W(NXPI)-U1*W(NGPI)-U2*W(I)
	      END DO
C           CALCULATE THE DERIVATIVE ALONG THE NEW SEARCH VECTOR.
	   END IF		! End if nrst /= n

	   DG1=0.
	   DO I=1,N
	      NXPI=NX+I
	      W(I)=W(NXPI)
	      DG1=DG1+W(I)*G(I)
	   END DO
C        IF THE NEW DIRECTION IS NOT A DESCENT DIRECTION,STOP.
	   IF (DG1.gt.0.) THEN	!    GO TO 320
	      NFLAG = 3
	      RETURN
	   END IF 

C        UPDATE NRST TO ASSURE AT LEAST ONE RESTART EVERY N ITERATIONS.
	   IF (NRST .eq. N) NRST=0
	   NRST=NRST+1
	   RSW=.FALSE.
	   GO TO 40
C        A VARIABLE METRIC ALGORITM IS BEING USED. CALCULATE Y AND D'Y.
	END IF			! End if nmeth /= 1
      
	U1=0.0
	DO I=1,N
	   NGPI=NG+I
	   W(NGPI)=G(I)-W(NGPI)
	   U1=U1+W(I)*W(NGPI)
	END DO

C     IF RSW=.TRUE.,SET UP THE INITIAL SCALED APPROXIMATE HESSIAN.
	IF (RSW) THEN
C        CALCULATE Y'Y.
	   U2=0.
	   DO I=1,N
	      NGPI=NG+I
	      U2=U2+W(NGPI)*W(NGPI)
	   END DO
C        CALCULATE THE INITIAL HESSIAN AS H=(P'Y/Y'Y)*I
C        AND THE INITIAL U2=Y'HY AND W(NX+I)=HY.
	   IJ=1
	   U3=U1/U2
	   DO I=1,N
	      DO J=I,N
		 NCONS1=NCONS+IJ
		 W(NCONS1)=0.0
		 IF(I.EQ.J)W(NCONS1)=U3
		 IJ=IJ+1
	      END DO
	      NXPI=NX+I
	      NGPI=NG+I
	      W(NXPI)=U3*W(NGPI)
	   END DO
	   U2=U3*U2
	ELSE
C        CALCULATE W(NX+I)=HY AND U2=Y'HY.
	   U2=0.0
	   DO I=1,N
	      U3=0.0
	      IJ=I
	      IF (I .ne. 1) THEN
		 II=I-1
		 DO J=1,II
		    NGPJ=NG+J
		    NCONS1=NCONS+IJ
		    U3=U3+W(NCONS1)*W(NGPJ)
		    IJ=IJ+N-J
		 END DO 
	      END IF
	      DO J=I,N
		 NCONS1=NCONS+IJ
		 NGPJ=NG+J
		 U3=U3+W(NCONS1)*W(NGPJ)
		 IJ=IJ+1
	      END DO 
	      NGPI=NG+I
	      U2=U2+U3*W(NGPI)
	      NXPI=NX+I
	      W(NXPI)=U3
	   END DO
	END IF

C     CALCULATE THE UPDATED APPROXIMATE HESSIAN.
	U4=1.0+U2/U1
	DO I=1,N
	   NXPI=NX+I
	   NGPI=NG+I
	   W(NGPI)=U4*W(I)-W(NXPI)
	END DO
	IJ=1
	DO I=1,N
	   NXPI=NX+I
	   U3=W(I)/U1
	   U4=W(NXPI)/U1
	   DO J=I,N
	      NCONS1=NCONS+IJ
	      NGPJ=NG+J
	      W(NCONS1)=W(NCONS1)+U3*W(NGPJ)-U4*W(J)
	      IJ = IJ + 1
	   END DO
	END DO

C     CALCULATE THE NEW SEARCH DIRECTION W(I)=-HG AND ITS DERIVATIVE.
	DG1=0.0
	DO I=1,N
	   U3=0.0
	   IJ=I
	   IF (I .ne. 1) THEN
	      II=I-1
	      DO J=1,II
		 NCONS1=NCONS+IJ
		 U3=U3-W(NCONS1)*G(J)
		 IJ=IJ+N-J
	      END DO
	   END IF   	
	   DO J=I,N
	      NCONS1=NCONS+IJ
	      U3=U3-W(NCONS1)*G(J)
	      IJ=IJ+1
	   END DO
	   DG1=DG1+U3*G(I)
	   W(I)=U3
	END DO

C     TEST FOR A DOWNHILL DIRECTION.
	IF (DG1.gt.0.) THEN      
	   NFLAG = 3
	   RETURN
	ELSE
	   RSW=.FALSE.
	   GO TO 40
	END IF

	End !Subroutine CONMINBN	! End subroutine CONMINBN

C ========================================================================

C!! The End
