c     
c     repeated : A Library of Repeated Measurements Models
c     Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
c     
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or
c     (at your option) any later version.
c     
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c     
c     SYNOPSIS
c     
c     subroutine logitord_f(y,upk,EPS,FCALLS,iout,cg,total1,total2a,
c     &       total2b,nobs,x,ster,hess,hesinv,nflag,iter,ifun,f)
c     
c     DESCRIPTION
c     
c     Function to compute the likelihood of an longitudinal ordinal response
c     model with dropouts and common random effect for the two processes
c     
*     =========================================================================
*     name 			: Logitord.f  --- MULTIPLE SIGMAS
*     programmer 		: Euginia Zharichenko
*     date completed 		: July, 1995
*     Development History	: 
*     
*     
*     Modified by Ten Have for ordinal response, Oct, 1997
*     Modified by J.K. Lindsey for use with R, July, 1998
*     
*     
*     =========================================================================
*     This program is reading a user-suplied file 'users_x' 
*     (locating in the same directory as this program) which contains 
*     information about the name of the input datafile, value of epsilon, 
*     maximum number of funstion calls, number of betas, number of sigmas, 
*     and initial values for those estimates.
*     
*     !!!!         USER SHOULD SUPPLY THE FILE 'users_x'        !!!!
*     with all information on the same line and separated by one space: 
*     Name of file     - up to 20 characters,
*     Epsilon,
*     Maximum number of function calls,
*     Number of betas  - up to maxbet, 
*     Number of sigmas - up to maxsig, 
*     Initial estimates.   
*     
*     !!! Number of initial estimates suppose to be EQUAL to the
*     !!! "number of estimates".
*     
*     EXAMPLE :
*     xfile.dat 10.0 0.0001 300 2 1 -1.2 1.3 0.9
*     Do not include any other information in the "users_x" but that 
*     described previously.
*     =========================================================================


*     GLOBAL variables :
*     maxsub -- maximum number of subjects 
*     maxcas -- maximum number of cases inside each subject
*     maxbet -- maximum number of betas 
*     maxsig -- maximum number of sigmas 
*     maxest -- maximum number of estimates = maxbet+maxsig 
*     max_w  -- maximum size of vector w = maxest*(maxest+7)/2 
*     
*     OTHER variables : 
*     fname  -- name of the data file
*     total  -- total number of estimates = total1+total2
*     total1 -- total (real) number of betas 
*     total2a -- total (real) number of sigmas for random intercept
*     total2b -- total (real) number of sigmas for random slope
*     total2 -- total (real) number of sigmas = total2a + total2b
*     total3 -- total (real) number of subjects
*     total4 -- real dimensin of W
*     x      -- vector of estimates
*     g      -- vector of gradients
*     hess   -- hessian matrix
*     hesinv -- hessian invert matrix
*     ster   -- standard errors
*     EPS    -- epsilon - stopping criteria
*     FCALLS -- number of function calls
*     F      -- value of our function sli
*     id     -- subject's id number
*     numcas -- number of cases within one mother
*     
*     TEMPORARY variables :
*     aa,aaa, i,j,j1,j2,jj
*     ========================================================================

      subroutine logitord_f(y,upk,EPS,FCALLS,iout,cg,total1,total2a,
     &     total2b,nobs,p,x,ster,hess,hesinv,nflag,iter,ifun,f)
c     Begin main program

      implicit none  
c     Declaration Section

c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxbet,maxest,max_w	   
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)
      PARAMETER (max_w=maxest*(maxest+7)/2)
      
c     INTEGER Declaration

c     INTEGER
      INTEGER total,total1,total2a,total2b,total3,total4
      INTEGER total1x,cg,nobs,iout
      INTEGER i,j,jj,ii
      INTEGER nmeth,idev,IFun,iter,nflag,upk,FCALLS

c     INTEGER Array
      INTEGER id(maxsub), numcas(maxsub)

c     Double Precision Declaration

c     Double Precision
      Double Precision EPS,F,ACC,aa

c     One Dim Array
      DOUBLE PRECISION x(total1+total2a+total2b)
      DOUBLE PRECISION p(total1+total2a+total2b)
      DOUBLE PRECISION w(max_w),g(maxest)
      double precision ster(total1+total2a+total2b)

c     Two Dim Array      	  
      DOUBLE PRECISION hess(total1+total2a+total2b,
     +     total1+total2a+total2b), hesinv(total1+total2a+total2b,
     +     total1+total2a+total2b)
      DOUBLE PRECISION ri(maxsub,maxcas)
      double precision y(nobs,3+total1-cg+total2a+total2b)
      
c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     +     uu2(maxsub,maxcas,maxsig)
      
c     *** End Declaration Section                      

c     total # of estimates
      total=total1+total2a+total2b
      total1x=total1-cg+1
      do jj=1,total
         x(jj)=p(jj)
      enddo
c     Calculate dimension of vector W
      total4=total*(total+7)/2

      i=1
      j=1
      
      id(i)=y(1,1)
      ri(i,j)=y(1,2)
      do jj=1,total1x
         z(i,j,jj)=y(1,2+jj)
      enddo
      do jj=1,total2a
         uu1(i,j,jj)=y(1,2+total1x+jj)
      enddo
      do jj=1,total2b
         uu2(i,j,jj)=y(1,2+total1x+total2a+jj)
      enddo
c     Begin loop
 1    i=i+1
      do ii=2,nobs
         id(i)=y(ii,1)
         IF (id(i) .eq. id(i-1)) then
            j=j+1
            ri(i-1,j) = y(ii,2)
            IF (total1x .gt. 0) then
               DO jj=1,total1x
                  z(i-1,j,jj) = y(ii,jj+2)
               END DO
            ENDIF
            IF (total2a .ge. 1) then 
c                                          If there is more than 1 sigma
               DO jj=1,total2a
                  uu1(i-1,j,jj)=y(ii,2+total1x+jj)
               END DO
            ENDIF
            IF (total2b .ge. 1) then
               DO jj=1,total2b
                  uu2(i-1,j,jj)=y(ii,2+total1x+total2a+jj)
               END DO
            ENDIF
         else                   
c     Number of case inside i-th mother.
            numcas(i-1)=j
            j=1
            ri(i,j) = y(ii,2)
            IF (total1x .gt. 0)then
               DO jj=1,total1x
                  z(i,j,jj)=y(ii,jj+2)
               END DO
            ENDIF
            IF (total2a .ge. 1) then
               DO jj=1,total2a
                  uu1(i,j,jj)=y(ii,2+jj+total1x)
               END DO
            ENDIF
            IF (total2b .ge. 1) then
               DO jj=1,total2b
                  uu2(i,j,jj)=y(ii,2+jj+total1x+total2a)
               END DO
            ENDIF
            i=i+1
         ENDIF
      enddo

c     Total number of subjects (mothers) :
      total3=i-1
      numcas(total3)=j
      
      ACC = 10.D - 20
      nmeth = 1
      idev = 6

c     Call procedure CONMIN
      CALL CONMIN(upk,x,f,g,hess,ifun,iter,EPS,nflag,FCALLS,w,
     &     iout,total4,idev,ACC,nmeth,ri,z,uu1,uu2,
     &     total1,cg,total2a,total2b,total3,numcas)

      If (nflag .eq. 0) Then
         CALL CALCFG2(upk,X,total1,cg,total2a,total2b,total3,z,
     &        uu1,uu2,ri,numcas,F,G,HESS)
      End If 

c     invert HESSIAN matrix -- Call Invert procedure
      CALL invert(hess,hesinv,total)

c     find standard errors (square root of the negative main diagonal)
      DO  i=1,total
         aa= -hesinv(i,i)
         STER(i)=dsqrt(aa)
      END DO

      END
c     PROGRAM LOGITORD -- End main program
c==========================================================================
c      Begin Subroutine CALCFG	

      SUBROUTINE CALCFG(upk_temp,x,total1,cg,total2a,total2b,total3,z,
     &     uu1,uu2,ri,numcas,sli,g,hess)

*     =========================================================================
*     PURPOSE : Subroutine CALCFG calculates logistic-binomial derivatives
*     =======   and the object function value.
*     
*     INPUT :
*     =====
*     TOTAL1     - the number of variables in the function to be
*     minimized
*     TOTAL2
*     TOTAL3
*     X          - the vector containing the initial estimates BETAs and
*     SIGMA (must be dimensioned TOTAL1). 
*     BETA - a vector of fixed effects in logistic-regresion
*     model
*     SIGMA - a vector of variance components for random 
*     effects parameters
*     TOTAL3     - number of subjects
*     Z          - a matrix of covariates in logistic-regresion
*     model for total3 subjects (must be dimensioned
*     TOTAL3 x TOTAL1-1)
*     U
*     Ri         - vector of numbers of successes within subject i for
*     total3 subjects (must be dimensioned TOTAL3)
*     for total3 subjects (must be dimensioned TOTAL3)
*     NUMCAS     - number of different cases for each subject-mother
*     
*     OUTPUT :
*     ======
*     SLi        - the function value
*     G          - an array for the gradient with all first derivatives
*     in the following order: w.r.t. beta j's, w.r.t. sigma j's
*     (must be dimensioned TOTAL1)
*     
*     Other variable names used in this program : 
*     -----
*     u          - a covariate for the random effects
*     upperK     - a number of binomial distributions
*     qi         - proportion of binomial distribution (equal to 0.5)
*     Li         - the likelihood contribution from the i-th set
*     dbeta      - a vector of first derivatives of sli with respect to beta(j)
*     dsigma     - a vector of first derivatives of sli with respect to sigma(j)
*     j,jj       - loop variable
*     
*     ========================================================================

      implicit none
c     implicit DOUBLE PRECISION (A-Z)
c     Declaration Section

c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxbet,maxest 
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)
c     Integer Declaration

c     Integer Variables
      INTEGER total,total1,total1x,total2a,total2b,total3
      INTEGER j, jj, j2
      INTEGER i, j1, upk,upk_temp,cg 

c     Integer Array
      INTEGER numcas(maxsub)
      
c     Double Precision Declaration

c     Double Precision Variables
      DOUBLE PRECISION DLOG, li, sli, sqi1, sqi2 

c     One Dim Array
      DOUBLE PRECISION x(total1+total2a+total2b), g(maxest)	
      DOUBLE PRECISION sigmaj1(maxcas), sigmaj2(maxcas),betaj(maxcas)
      DOUBLE PRECISION dbeta(maxbet), sbetaj(maxbet),beta(maxbet)
      DOUBLE PRECISION dsigma1(maxsig), ssigmn1(maxsig)
      DOUBLE PRECISION dsigma2(maxsig), ssigmn2(maxsig)       
      DOUBLE PRECISION ssigmn12(maxsig),sigma1(maxsig),sigma2(maxsig) 

c     Two Dim Array
      DOUBLE PRECISION hess(total1+total2a+total2b,
     +     total1+total2a+total2b)
      DOUBLE PRECISION s2bjbl(maxbet,maxbet)
      DOUBLE PRECISION ri(maxsub,maxcas)
      DOUBLE PRECISION s2bjsn1(maxbet,maxsig),s2bjsn2(maxbet,maxsig)
      DOUBLE PRECISION s2snso2(maxsig,maxsig),s2snso12(maxsig,maxsig)
      DOUBLE PRECISION s2snso1(maxsig,maxsig)
      
c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     &     uu2(maxsub,maxcas,maxsig)
      
c     End Declaration Section
c     Upper K = 10 & can be changed      
c     Small Qi = 0.5 & can be changed
      total=total1+total2a+total2b           
      total1x=total1-cg+1
      upk = upk_temp            
      sqi1 = 0.5d0		
      sqi2 = 0.5d0

c     read BETA values from the vector X
      DO j=1,total1x
         beta(j)=x(cg-1+j)
      END DO 

c     read SIGMA values from the vector X
      DO j=1,total2a
         sigma1(j)=x(total1+j)
      END DO 

      DO j=1,total2b
         sigma2(j)=x(total1+total2a+j)
      END DO
c     Initialization of the G-vector,HESS matrix and function sli
      DO j=1,total
         G(j)=0.0d0
         DO j2=1,total
            HESS(j,j2)=0.0d0
         END DO
      END DO

      sli=0.0d0
c     calculate everything for each subject and 
c     sum all information at the end 
c     Begin Main Loop                    
      Do jj=1,total3
c     product of one row from matrix z to the       
c     vector beta for the formula (5)      
         DO i=1,numcas(jj)	    
            betaj(i)=0.0d0      
            DO j=1,total1x
               betaj(i)=betaj(i)+z(jj,i,j)*beta(j)    
            END DO				      
c     product of one row from matrix u to the 
c vector sigma for the formula (5) 
            sigmaj1(i)=0.0d0    
            sigmaj2(i)=0.0d0
            DO j=1,total2a
               sigmaj1(i)=sigmaj1(i)+uu1(jj,i,j)*sigma1(j)
            END DO 
            DO j=1,total2b
               sigmaj2(i)=sigmaj2(i)+uu2(jj,i,j)*sigma2(j)
            END DO
         END DO 
c     Call subroutine FORMUL
         CALL FORMUL(betaj,sigmaj1,sigmaj2,upK,sqi1,sqi2,
     &        ri,jj,numcas,
     &        total1,total1x,cg,total2a,
     &        total2b,z,uu1,uu2,Li,sbetaj,ssigmn1,
     &        s2bjbl,s2bjsn1,s2snso1,ssigmn2,s2bjsn2,
     &        s2snso2,ssigmn12,s2snso12,x)
         
         if(Li.gt.0.)then
c        formula (1)         
            sli=sli+dlog(Li)
            
            DO j=1,total1
c     first derivatives w.r.t. beta j's            
               dbeta(j)=sbetaj(j)/Li
c     FINAL array for the gradient G(j)               
               G(j)=G(j)+dbeta(j) 
            END DO
            
            DO j=1,total2a 
c     First derivative w.r.t. sigma            
               dsigma1(j) = ssigmn1(j)/Li 
               G(total1+j) = G(total1+j)+dsigma1(j)
            END DO

            DO j=1,total2b 
               dsigma2(j)=ssigmn2(j)/Li
               G(total1+total2a+j)=G(total1+total2a+j)+
     &              dsigma2(j)
            END DO
         endif
      END DO                    
c     end main loop
c     take "-" in order to change Max to Min
      sli=-sli                  
      
c     FINAL Hessian Matrix for all second partial derivatives
      DO j1=1,total             
         g(j1) = - g(j1)        
c     take "-" in order to change Max to Min
      END DO
      RETURN	

      End
c     Subroutine CALCFG     End Subroutine CALCFG
c=========================================================================

c      Begin Subroutine CALCFG2

      SUBROUTINE CALCFG2(upk_temp,x,total1,cg,total2a,total2b,total3,z,
     &     uu1,uu2,ri,numcas,sli,g,hess)

*     =========================================================================
*     PURPOSE : Subroutine CALCFG2 calculates logistic-binomial derivatives
*     =======   and the object function value.
*     
*     INPUT :
*     =====
*     TOTAL1     - the number of variables in the function to be
*     minimized
*     TOTAL2
*     TOTAL3
*     X          - the vector containing the initial estimates BETAs and
*     SIGMA (must be dimensioned TOTAL1). 
*     BETA - a vector of fixed effects in logistic-regresion
*     model
*     SIGMA - a vector of variance components for random 
*     effects parameters
*     TOTAL3     - number of subjects
*     Z          - a matrix of covariates in logistic-regresion
*     model for total3 subjects (must be dimensioned
*     TOTAL3 x TOTAL1-1)
*     U
*     Ri         - vector of numbers of successes within subject i for
*     total3 subjects (must be dimensioned TOTAL3)
*     for total3 subjects (must be dimensioned TOTAL3)
*     NUMCAS     - number of different cases for each subject-mother
*     
*     
*     OUTPUT :
*     ======
*     SLi        - the function value
*     G          - an array for the gradient with all first derivatives
*     in the following order: w.r.t. beta j's, w.r.t. sigma j's
*     (must be dimensioned TOTAL1)
*     
*     Other variable names used in this program : 
*     -----
*     u          - a covariate for the random effects
*     upperK     - a number of binomial distributions
*     qi         - proportion of binomial distribution (equal to 0.5)
*     Li         - the likelihood contribution from the i-th set
*     dbeta      - a vector of first derivatives of sli with respect to beta(j)
*     dsigma     - a vector of first derivatives of sli with respect to sigma(j)
*     j,jj       - loop variable
*     
*     ========================================================================

      implicit none
c     implicit DOUBLE PRECISION (A-Z)
c     Declaration Section

c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxbet,maxest 
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)
      
c     Integer Declaration

c     Integer Variables
      INTEGER total,total1,total1x,total2a,total2b,total3
      INTEGER j, jj, j2
      INTEGER i, j1, upk,upk_temp,cg 

c     Integer Array
      INTEGER numcas(maxsub)
      
c     Double Precision Declaration

c     Double Precision Variables
      DOUBLE PRECISION DLOG, li, sli, sqi1, sqi2

c     One Dim Array
      DOUBLE PRECISION x(total1+total2a+total2b), g(maxest)	
      DOUBLE PRECISION sigmaj1(maxcas), sigmaj2(maxcas),betaj(maxcas)
      DOUBLE PRECISION dbeta(maxbet), sbetaj(maxbet),beta(maxbet)
      DOUBLE PRECISION dsigma1(maxsig), ssigmn1(maxsig)
      DOUBLE PRECISION dsigma2(maxsig), ssigmn2(maxsig)       
      DOUBLE PRECISION ssigmn12(maxsig),sigma1(maxsig),sigma2(maxsig)

c     Two Dim Array
      DOUBLE PRECISION hess(total1+total2a+total2b,
     +     total1+total2a+total2b)
      DOUBLE PRECISION s2bjbl(maxbet,maxbet)
      DOUBLE PRECISION ri(maxsub,maxcas)
      DOUBLE PRECISION s2bjsn1(maxbet,maxsig),s2bjsn2(maxbet,maxsig)
      DOUBLE PRECISION s2snso2(maxsig,maxsig),s2snso12(maxsig,maxsig)
      DOUBLE PRECISION s2snso1(maxsig,maxsig)
      
c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     &     uu2(maxsub,maxcas,maxsig)
      
c     End Declaration Section
c     Upper K = 10 & can be changed      
c     Small Qi = 0.5 & can be changed
      total=total1+total2a+total2b           
      total1x=total1-cg+1
      upk = upk_temp            
      sqi1 = 0.5d0              
      sqi2 = 0.5d0

c     read BETA values from the vector X
      DO j=1,total1x
         beta(j)=x(cg-1+j)
      END DO 

c     read SIGMA values from the vector X
      DO j=1,total2a
         sigma1(j)=x(total1+j)
      END DO 

      DO j=1,total2b
         sigma2(j)=x(total1+total2a+j)
      END DO
      
c     Initialization of the G-vector,HESS matrix and function sli
      DO j=1,total
         G(j)=0.0d0
         DO j2=1,total
            HESS(j,j2)=0.0d0
         END DO
      END DO

      sli=0.0d0
c     calculate everything for each subject and 
c     sum all information at the end 
c     Begin Main Loop              
c     product of one row from matrix z to the 
c     vector beta for the formula (5)
      Do jj=1,total3            
         DO i=1,numcas(jj)      
            betaj(i)=0.0d0     
            DO j=1,total1x
               betaj(i)=betaj(i)+z(jj,i,j)*beta(j)    
            END DO
c     product of one row from matrix u to the 
            sigmaj1(i)=0.0d0      
c     vector sigma for the formula (5) 
            sigmaj2(i)=0.0d0
            DO j=1,total2a
               sigmaj1(i)=sigmaj1(i)+uu1(jj,i,j)*sigma1(j)
            END DO 
            DO j=1,total2b
               sigmaj2(i)=sigmaj2(i)+uu2(jj,i,j)*sigma2(j)
            END DO
         END DO 
         
c     Call subroutine FORMUL
         CALL FORMUL2(betaj,sigmaj1,sigmaj2,upK,sqi1,sqi2,
     &        ri,jj,numcas,
     &        total1,total1x,cg,total2a,
     &        total2b,z,uu1,uu2,Li,sbetaj,ssigmn1,
     &        s2bjbl,s2bjsn1,s2snso1,ssigmn2,s2bjsn2,
     &        s2snso2,ssigmn12,s2snso12,x)
c     formula (1)         
         sli=sli+dlog(Li)       
         
         DO j=1,total1
c     first derivatives w.r.t. beta j's         
            dbeta(j)=sbetaj(j)/Li
c     FINAL array for the gradient G(j)            
            G(j)=G(j)+dbeta(j)  
         END DO
         
         DO j=1,total2a 
c     First derivative w.r.t. sigma         
            dsigma1(j) = ssigmn1(j)/Li 
            G(total1+j) = G(total1+j)+dsigma1(j)
         END DO

         DO j=1,total2b 
            dsigma2(j)=ssigmn2(j)/Li
            G(total1+total2a+j)=G(total1+total2a+j)+
     &           dsigma2(j)
         END DO	      

c     second order derivatives
*     ====================================================================
*     matrix HESS(j1,j2) for the second partial derivatives 
*     HESS has the following order:
*     beta(1) beta(2) ... beta(n) sigma(1) sigma(2) ... sigma(m) 
*     beta(1)    x       x     ...   x       x        x      ...   x      
*     beta(2)    x       x           x       x        x            x      
*     ...        .       .     ...   .       .        .      ...   .      
*     beta(n)    x       x           x       x        x      ...   x      
*     sigma(1)   x       x     ...   x       x        x      ...   x      
*     sigma(2)   x       x     ...   x       x        x      ...   x      
*     ...        .       .     ...   .       .        .      ...   .      
*     sigma(m)   x       x     ...   x       x        x      ...   x      
*     where n=total1, m=total2.  
*     
*     1.Derivatives with respect to beta(j1) and beta(j2)
*     ====================================================================

         DO j1=1,total1 
            DO j2=j1,total1
               HESS(j1,j2)=HESS(j1,j2)+s2bjbl(j1,j2)/Li - 
     &              dbeta(j1)*dbeta(j2)
            END DO 
c     2.Derivatives with respect to beta(j1) and sigma(j2)
            DO j2=1,total2a
               HESS(j1,total1+j2)=HESS(j1,total1+j2)+
     &              s2bjsn1(j1,j2)/Li - dbeta(j1)*dsigma1(j2)
            END DO              
            DO j2=1,total2b
               HESS(j1,total1+total2a+j2)=HESS(j1,total1+
     & 		    total2a+j2)+ s2bjsn2(j1,j2)/Li -
     &              dbeta(j1)*dsigma2(j2)
            END DO
         END DO 

c     4. Derivatives with respect to sigma(j1) and sigma(j2)
         DO j1=1,total2a
            DO j2=1,total2a
               HESS(total1+j1,total1+j2)=HESS(total1+j1,total1+j2)
     &              +s2snso1(j1,j2)/Li-dsigma1(j1)*dsigma1(j2)
            END DO 
            DO j2=1,total2b
c     5. Der. w/respect to sigma1(j1) & sigma2(j2)
               HESS(total1+j1,total1+total2a+j2)=
     &              HESS(total1+j1,total1+total2a+j2)+
     &              s2snso12(j1,j2)/Li-dsigma1(j1)*dsigma2(j2)
c     bug corrected JKL
            enddo
         enddo
         do j1=1,total2b
            DO j2=1,total2b
               HESS(total1+total2a+j1,total1+total2a+j2)=
     &              HESS(total1+total2a+j1,total1+total2a+j2)+
     &              s2snso2(j1,j2)/Li-dsigma2(j1)*dsigma2(j2)
            END DO
         END DO    
      End DO                    
c     End main loop      
      sli=-sli                  
c     take "-" in order to change Max to Min      
c     FINAL Hessian Matrix for all second partial derivatives
      DO j1=1,total
         DO j2=1,total
            IF (j1 .gt. j2) THEN
               HESS(j1,j2) = HESS(j2,j1)
            END IF
         END DO 
         g(j1) = - g(j1)        
c     take "-" in order to change Max to Min
      END DO
      
      RETURN	

      End                        
c     Subroutine CALCFG2    End Subroutine CALCFG2
c=========================================================================

*     =========================================================================
*     This subroutine "FORMUL" calculates sums (from k=0 to K) of different
*     functions in order to calculate derivatives in the main program.
*     Input : same parameters as in the main program plus variable t and 
*     vectors dbeta(j) and dsigma(j) that were calculated in the 
*     main program.
*     Output: summations for the main program 
*     Li     - variable Li equal to the sum of all M(ri,ni,k)      
*     Sbetaj - sum of Pbetaj=M(ri,ni,k)*N(ri,ni,k) for the first 
*     derivative of function sli with respect to betaj   
*     Ssigmn - sum of Psigmn=M(ri,ni,k)*N(ri,ni,k)*v(k) for the first
*     derivative of function sli with respect to sigman 
*     S2bjbl - sum of P2bjbl=M(ri,ni,k)*O(ri,ni,k) for the second 
*     derivative of sli with respect to betaj and betal
*     S2bjsn - sum of P2bjsn=M(ri,ni,k)*O(ri,ni,k)*v(k) for the
*     second order derivative of sli with respect to
*     betaj and sigman
*     S2snso - sum of P2snso for the second derivative of sli with
*     respect to sigman and sigmao
*     ========================================================================



c     Begin Subroutine FORMUL	

      SUBROUTINE FORMUL(betaj,sigmaj1,sigmaj2,upK,sqi1,sqi2,
     &     ri,jj,numcas,
     &     total1,total1x,cg,total2a,total2b,
     &     z,uu1,uu2,Li,sbetaj,ssigmn1,s2bjbl,
     &     s2bjsn1,s2snso1,ssigmn2,s2bjsn2,s2snso2,
     &     ssigmn12,s2snso12,x)


      implicit none

c     Declaration Section
      
c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxbet,maxest
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)

c     Integer Declaration 

c     Integer Variables
      INTEGER i,j,jj
      INTEGER total1, total1x, total2a, total2b
      INTEGER upperk, lowerk1, lowerk2, upk,cg
      INTEGER        y1

c     Integer Array
      INTEGER numcas(maxsub)
      
      
c     Double Precision Declaration

c     Double Precision Variables
      DOUBLE PRECISION Li,Ai,Pi,Mi
      DOUBLE PRECISION V1,V2, sqi1, sqi2
      DOUBLE PRECISION ukfact, lkfact1, lkfact2
      DOUBLE PRECISION kkq1, kkq2, kkfact1, kkfact2
      DOUBLE PRECISION kchoos2, kchoos1, FACTOR
      DOUBLE PRECISION expon1, expon2, gamma1, gamma2
      DOUBLE PRECISION d1, d2, dd1, dd2

      DOUBLE PRECISION sigmaj1(maxcas),sigmaj2(maxcas)
      DOUBLE PRECISION betaj(maxcas),Nij(maxcas)
c     ,Oij(maxcas)
      DOUBLE PRECISION Niz(maxbet),sbetaj(maxbet)
      DOUBLE PRECISION x(total1+total2a+total2b)
      DOUBLE PRECISION ssigmn1(maxsig),ssigmn2(maxsig),ssigmn12(maxsig)
      DOUBLE PRECISION Niuu1(maxsig), Niuu2(maxsig)
      DOUBLE PRECISION uu(maxsig)

c     Two Dim Array
      DOUBLE PRECISION s2bjbl(maxbet,maxbet)
c     , Ozz(maxbet,maxbet)
      DOUBLE PRECISION ri(maxsub,maxcas)
      DOUBLE PRECISION s2snso1(maxsig,maxsig),s2snso2(maxsig,maxsig)
      DOUBLE PRECISION s2snso12(maxsig,maxsig)
c     ,Ouuu12(maxsig,maxsig)
c     DOUBLE PRECISION Ouuu1(maxsig,maxsig),Ouuu2(maxsig,maxsig)
      DOUBLE PRECISION s2bjsn1(maxbet,maxsig),s2bjsn2(maxbet,maxsig) 
c     DOUBLE PRECISION Ozuu1(maxbet,maxsig),Ozuu2(maxbet,maxsig)

c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     &     uu2(maxsub,maxcas,maxsig)
      
c     End Declaration Section

      
      DO j=1,cg-1
         uu(j)=x(j)
      end do
      

c     init of variables for deri. wrt beta      
      DO j=1,total1
         Sbetaj(j)=0.0d0      
      END DO 
      
c     Inititalization of Var for derivatives wrt sigma
      DO j = 1, total2a
         ssigmn1(j) = 0
      END DO   

      DO j = 1, total2b
         ssigmn2(j) = 0
      END DO   
      
      
      Li=0.0d0 
      
      upperK=upK-1
c     dec upperK by 1 to compare results from EGRET w/same K
      ukfact=FACTOR(upperK)
c     loop for summations from k=0 to K 
      DO lowerk1=0,upperK
         kkq1=lowerk1-upperK*sqi1
         v1=kkq1/dsqrt(upperK*sqi1*(1-sqi1))
         lkfact1=FACTOR(lowerk1)
         Kkfact1=FACTOR(upperK-lowerk1)
         kchoos1=uKfact/(lkfact1*Kkfact1)
         DO lowerk2=0,upperK		 
c     formula (4)                
            kkq2=lowerk2-upperK*sqi2                 
            v2=kkq2/dsqrt(upperK*sqi2*(1-sqi2))
c     initialization of some variables and arrays            
            Ai=1.0d0            	     
            DO j=1,total1 
               Niz(j)=0.0d0
            END DO

            DO j=1,total2a
               Niuu1(j)=0
            END DO  

            DO j=1,total2b
               Niuu2(j)=0
            END DO 

            DO i=1,numcas(jj)                        
               y1=ri(jj,i)

               if (y1 .eq. cg) then
c     lines added JKL
                  expon1=uu(y1-1)+betaj(i)+
     *                 sigmaj1(i)*v1+sigmaj2(i)*v2
                  if(expon1.gt.25)then
                     if(betaj(i).gt.15.)betaj(i)=betaj(i)/2
                     if(sigmaj1(i)*v1.gt.15.)
     +                    sigmaj1(i)=sigmaj1(i)/dabs(2*v1)
                     if(sigmaj2(i)*v2.gt.15.)
     +                    sigmaj2(i)=sigmaj2(i)/dabs(2*v2)
                     expon1=uu(y1-1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2
c     write(*,*)'overflow',i,jj,expon1,v1,v2,
c     +                    uu(y1-1),betaj(i),sigmaj1(i),sigmaj2(i)
                  endif
c     expon1=dexp(uu(y1-1)+betaj(i)+
c     *                  sigmaj1(i)*v1+sigmaj2(i)*v2)
                  expon1=dexp(expon1)
c     formula (5)
                  gamma2=1
                  gamma1=expon1/(1+expon1)
                  d2=0
                  d1=gamma1*(1-gamma1)
                  dd2=0
                  dd1=(1-2*gamma1)*d1
               else
                  if (y1 .eq. 1) then
                     expon2=dexp(uu(y1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2) 
                     gamma2=expon2/(1+expon2)
                     gamma1=0
                     d2=gamma2*(1-gamma2)
                     d1=0
                     dd2=(1-2*gamma2)*d2
                     dd1=0
                  else
                     expon2=dexp(uu(y1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2) 
                     expon1=dexp(uu(y1-1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2) 
                     gamma2=expon2/(1+expon2)
                     gamma1=expon1/(1+expon1)
                     d2=gamma2*(1-gamma2)
                     d1=gamma1*(1-gamma1)
                     dd2=(1-2*gamma2)*d2
                     dd1=(1-2*gamma1)*d1
                  endif
               endif
               pi=gamma2-gamma1
               Nij(i)=(d2-d1)/pi
               
c     for the formula (3)                                        
               lkfact2=FACTOR(lowerk2)                   
               Kkfact2=FACTOR(upperK-lowerk2)                   
               kchoos2=uKfact/(lkfact2*Kkfact2)

c     for the formula (3)
               Ai=Ai*Pi
               
               if (y1 .gt. 1) then
                  Niz(y1-1)=Niz(y1-1)-d1/pi
               endif
               
               if (y1 .lt. cg ) then
                  Niz(y1)=Niz(y1)+d2/pi
               endif
               
               
               DO j=cg,cg-1+total1x
                  Niz(j)=Niz(j)+z(jj,i,j-cg+1)*Nij(i)
               END DO
               
c     for the formula (17)
               DO j=1,total2a
                  Niuu1(j)=Niuu1(j)+uu1(jj,i,j)*Nij(i)
c     for the formulas (11),(15), and (17)
               END DO

c     for the formula (17)
               DO j=1,total2b
                  Niuu2(j)=Niuu2(j)+uu2(jj,i,j)*Nij(i)
c     for the formulas (11),(15), and (17)
               END DO

c     formula (3)
               Mi=Ai*
     &              kchoos1*sqi1**lowerk1*(1-sqi1)**(upperK-lowerk1)* 
     &              kchoos2*sqi2**lowerk2*(1-sqi2)**(upperK-lowerk2) 
            END DO              

            
c     products for derivatives of 1st order ( formulas (10)-(11) )
            DO j=1,total1
               Sbetaj(j)=Sbetaj(j)+Mi*Niz(j)
            END DO

            DO j=1,total2a
               Ssigmn1(j)=Ssigmn1(j)+v1*Mi*Niuu1(j)
            END DO

            DO j=1,total2b
               Ssigmn2(j)=Ssigmn2(j)+v2*Mi*Niuu2(j)
            END DO

            Li=Li+Mi            
         END DO
      END DO                    

      RETURN
      End                       

c======================================================================

*     =========================================================================
*     This subroutine "FORMUL2" calculates sums (from k=0 to K) of different
*     functions in order to calculate derivatives in the main program.
*     Input : same parameters as in the main program plus variable t and 
*     vectors dbeta(j) and dsigma(j) that were calculated in the 
*     main program.
*     Output: summations for the main program 
*     Li     - variable Li equal to the sum of all M(ri,ni,k)      
*     Sbetaj - sum of Pbetaj=M(ri,ni,k)*N(ri,ni,k) for the first 
*     derivative of function sli with respect to betaj   
*     Ssigmn - sum of Psigmn=M(ri,ni,k)*N(ri,ni,k)*v(k) for the first
*     derivative of function sli with respect to sigman 
*     S2bjbl - sum of P2bjbl=M(ri,ni,k)*O(ri,ni,k) for the second 
*     derivative of sli with respect to betaj and betal
*     S2bjsn - sum of P2bjsn=M(ri,ni,k)*O(ri,ni,k)*v(k) for the
*     second order derivative of sli with respect to
*     betaj and sigman
*     S2snso - sum of P2snso for the second derivative of sli with
*     respect to sigman and sigmao
*     ========================================================================



c     Begin Subroutine FORMUL2

      SUBROUTINE FORMUL2(betaj,sigmaj1,sigmaj2,upK,sqi1,sqi2,
     &     ri,jj,numcas, 
     &     total1,total1x,cg,total2a,total2b,
     &     z,uu1,uu2,Li,sbetaj,ssigmn1,s2bjbl,
     &     s2bjsn1,s2snso1,ssigmn2,s2bjsn2,s2snso2,
     &     ssigmn12,s2snso12,x)



      implicit none

c     Declaration Section
      
c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxbet,maxest
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)

c     Integer Declaration 

c     Integer Variables
      INTEGER i,j,j2,jj
      INTEGER total1,total1x,total2a, total2b
      INTEGER upk, upperk, lowerk1, lowerk2,cg 
      INTEGER        y1

c     Integer Array
      INTEGER numcas(maxsub)
      
      
c     Double Precision Declaration

c     Double Precision Variables
      DOUBLE PRECISION Li,Ai,Pi,Mi
      DOUBLE PRECISION V1,V2, sqi1, sqi2
      DOUBLE PRECISION ukfact, lkfact1, lkfact2
      DOUBLE PRECISION kkq1, kkq2, kkfact1, kkfact2
      DOUBLE PRECISION kchoos2, kchoos1, FACTOR
      DOUBLE PRECISION expon1, expon2, gamma1, gamma2
      DOUBLE PRECISION d1, d2, dd1, dd2


c     One Dim Array
      DOUBLE PRECISION sigmaj1(maxcas),sigmaj2(maxcas)
      DOUBLE PRECISION betaj(maxcas),Nij(maxcas),Oij(maxcas)
      DOUBLE PRECISION Niz(maxbet),sbetaj(maxbet) 
      DOUBLE PRECISION x(total1+total2a+total2b)
      DOUBLE PRECISION ssigmn1(maxsig),ssigmn2(maxsig),ssigmn12(maxsig)
      DOUBLE PRECISION Niuu1(maxsig), Niuu2(maxsig)
      DOUBLE PRECISION uu(maxsig)

c     Two Dim Array
      DOUBLE PRECISION s2bjbl(maxbet,maxbet), Ozz(maxbet,maxbet)
      DOUBLE PRECISION ri(maxsub,maxcas)
      DOUBLE PRECISION s2snso1(maxsig,maxsig),s2snso2(maxsig,maxsig)
      DOUBLE PRECISION s2snso12(maxsig,maxsig),Ouuu12(maxsig,maxsig)
      DOUBLE PRECISION Ouuu1(maxsig,maxsig),Ouuu2(maxsig,maxsig)
      DOUBLE PRECISION s2bjsn1(maxbet,maxsig),s2bjsn2(maxbet,maxsig) 
      DOUBLE PRECISION Ozuu1(maxbet,maxsig),Ozuu2(maxbet,maxsig)

c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     &     uu2(maxsub,maxcas,maxsig)
      
c     End Declaration Section

      
      DO j=1,cg-1
         uu(j)=x(j)
      end do
      
      
      DO j=1,total1
         Sbetaj(j)=0.0d0        
         DO j2=1,total1
            S2bjbl(j,j2)=0
         END DO 
         DO j2=1,total2a
            S2bjsn1(j,j2)=0     
         END DO 	   
         DO j2=1,total2b
            S2bjsn2(j,j2)= 0
         END DO
      END DO 
      
c     Inititalization of Var for derivatives wrt sigma
      DO j = 1, total2a
         ssigmn1(j) = 0
         DO j2 = 1, total2a
            s2snso1(j,j2) = 0
         END DO
      END DO   

      DO j = 1, total2b
         ssigmn2(j) = 0
         DO j2 = 1, total2b
            s2snso2(j,j2) = 0
         END DO
      END DO   
      
c     Initialization of Var for derivatives wrt sigma1 & sigma2
      DO J = 1, total2a
         DO j2 = 1, total2b
            S2snso12(j,j2) = 0
         END DO
      END DO
      
      Li=0.0d0 
      
      upperK=upK-1
c     dec upperK by 1 to compare results from EGRET w/same K
      ukfact=FACTOR(upperK)
c     loop for summations from k=0 to K 
      DO lowerk1=0,upperK
         kkq1=lowerk1-upperK*sqi1
         v1=kkq1/dsqrt(upperK*sqi1*(1-sqi1))
         lkfact1=FACTOR(lowerk1)
         Kkfact1=FACTOR(upperK-lowerk1)
         kchoos1=uKfact/(lkfact1*Kkfact1)
         DO lowerk2=0,upperK		 
c     formula (4)                
            kkq2=lowerk2-upperK*sqi2                 
            v2=kkq2/dsqrt(upperK*sqi2*(1-sqi2))         
            Ai=1.0d0            
            DO j=1,total1 
               Niz(j)=0.0d0
               DO j2=1,total1
                  Ozz(j,j2)=0
               END DO
               DO i=1,total2a
                  Ozuu1(j,i)=0
               END DO
               DO i=1,total2b
                  Ozuu2(j,i)=0
               END DO
            END DO

            DO j=1,total2a
               Niuu1(j)=0
               DO j2=1,total2a
                  Ouuu1(j,j2)=0
               END DO
               DO j2=1,total2b
                  Ouuu12(j,j2)=0
               END DO
            END DO  

            DO j=1,total2b
               Niuu2(j)=0
               DO j2=1,total2b
                  Ouuu2(j,j2)=0
               END DO
            END DO 

            DO i=1,numcas(jj)                        
               y1=ri(jj,i)

               if (y1 .eq. cg) then
c     lines added JKL
                  expon1=uu(y1-1)+betaj(i)+
     *                 sigmaj1(i)*v1+sigmaj2(i)*v2
                  if(expon1.gt.25)then
                     if(betaj(i).gt.15.)betaj(i)=betaj(i)/2
                     if(sigmaj1(i)*v1.gt.15.)
     +                    sigmaj1(i)=sigmaj1(i)/dabs(2*v1)
                     if(sigmaj2(i)*v2.gt.15.)
     +                    sigmaj2(i)=sigmaj2(i)/dabs(2*v2)
                     expon1=uu(y1-1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2
c     write(*,*)'overflow',i,jj,expon1,v1,v2,
c     +                    uu(y1-1),betaj(i),sigmaj1(i),sigmaj2(i)
                  endif
c     expon1=dexp(uu(y1-1)+betaj(i)+
c     *                  sigmaj1(i)*v1+sigmaj2(i)*v2)
                  expon1=dexp(expon1)
c     formula (5)
                  gamma2=1
                  gamma1=expon1/(1+expon1)
                  d2=0
                  d1=gamma1*(1-gamma1)
                  dd2=0
                  dd1=(1-2*gamma1)*d1
               else
                  if (y1 .eq. 1) then
                     expon2=dexp(uu(y1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2)
c     formula (5)
                     gamma2=expon2/(1+expon2)
                     gamma1=0
                     d2=gamma2*(1-gamma2)
                     d1=0
                     dd2=(1-2*gamma2)*d2
                     dd1=0
                  else
                     expon2=dexp(uu(y1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2)
c     formula (5)
                     expon1=dexp(uu(y1-1)+betaj(i)+
     *                    sigmaj1(i)*v1+sigmaj2(i)*v2)
c     formula (5)
                     gamma2=expon2/(1+expon2)
                     gamma1=expon1/(1+expon1)
                     d2=gamma2*(1-gamma2)
                     d1=gamma1*(1-gamma1)
                     dd2=(1-2*gamma2)*d2
                     dd1=(1-2*gamma1)*d1
                  endif
               endif
               pi=gamma2-gamma1
               Nij(i)=(d2-d1)/pi
               Oij(i)=-Nij(i)**2 + (dd2-dd1)/pi
               
c     for the formula (3)                                        
               lkfact2=FACTOR(lowerk2)                   
               Kkfact2=FACTOR(upperK-lowerk2)                   
               kchoos2=uKfact/(lkfact2*Kkfact2)

c     for the formula (3)
               Ai=Ai*Pi
               
               if (y1 .gt. 1) then
                  Niz(y1-1)=Niz(y1-1)-d1/pi
                  Ozz(y1-1,y1-1)=Ozz(y1-1,y1-1)+(d1/pi)**2+dd1/pi
               endif
               
               if (y1 .lt. cg ) then
                  Niz(y1)=Niz(y1)+d2/pi
                  Ozz(y1,y1)=Ozz(y1,y1)+(d2/pi)**2-dd2/pi
               endif
               
               if (y1 .gt. 1 .and. y1 .lt. cg) then
                  Ozz(y1-1,y1)=Ozz(y1-1,y1)-d1*d2/pi**2
               endif
               
               
               DO j=cg,cg-1+total1x
                  Niz(j)=Niz(j)+z(jj,i,j-cg+1)*Nij(i)
                  if (y1 .lt. cg) then
                     Ozz(y1,j)=Ozz(y1,j)-z(jj,i,j-cg+1)*
     &                    (-d2*(d2-d1)/pi**2+dd2/pi)
                  endif
                  if (y1 .gt. 1) then
                     Ozz(y1-1,j)=Ozz(y1-1,j)-z(jj,i,j-cg+1)*
     &                    (d1*(d2-d1)/pi**2-dd1/pi)
                  endif
                  
                  DO j2=cg,cg-1+total1x
                     Ozz(j,j2)=Ozz(j,j2)-Oij(i)*z(jj,i,j-cg+1)*
     &                    z(jj,i,j2-cg+1)
                  END DO
                  
                  DO j2=1,total2a
                     Ozuu1(j,j2)=Ozuu1(j,j2)-Oij(i)*z(jj,i,j-cg+1)*
     &                    uu1(jj,i,j2)
                  END DO
                  DO j2=1,total2b
                     Ozuu2(j,j2)=Ozuu2(j,j2)-Oij(i)*z(jj,i,j-cg+1)*
     &                    uu2(jj,i,j2)
                  END DO
               END DO
               
c     for the formula (17)
               DO j=1,total2a
                  Niuu1(j)=Niuu1(j)+uu1(jj,i,j)*Nij(i)
c     for the formulas (11),(15), and (17)
                  if (y1 .lt. cg) then
                     Ozuu1(y1,j)=Ozuu1(y1,j)-uu1(jj,i,j)*
     &                    (-d2*(d2-d1)/pi**2+dd2/pi)
                  endif
                  
                  if (y1 .gt. 1) then
                     Ozuu1(y1-1,j)=Ozuu1(y1-1,j)-uu1(jj,i,j)*
     &                    (d1*(d2-d1)/pi**2-dd1/pi)
                  endif
                  DO j2=1,total2a
                     Ouuu1(j,j2)=Ouuu1(j,j2)-Oij(i)*
     &                    uu1(jj,i,j)*uu1(jj,i,j2)
                     
                  END DO
                  DO j2=1,total2b
                     Ouuu12(j,j2)=Ouuu12(j,j2)-Oij(i)*
     &                    uu1(jj,i,j)*uu2(jj,i,j2)
                     
                  END DO
               END DO

c     for the formula (17)
               DO j=1,total2b
                  Niuu2(j)=Niuu2(j)+uu2(jj,i,j)*Nij(i)
c     for the formulas (11),(15), and (17)
                  if (y1 .lt. cg) then
                     Ozuu2(y1,j)=Ozuu2(y1,j)-uu2(jj,i,j)*
     &                    (-d2*(d2-d1)/pi**2+dd2/pi)
                  endif
                  
                  if (y1 .gt. 1) then
                     Ozuu2(y1-1,j)=Ozuu2(y1-1,j)-uu2(jj,i,j)*
     &                    (d1*(d2-d1)/pi**2-dd1/pi)
                  endif
                  DO j2=1,total2b
                     Ouuu2(j,j2)=Ouuu2(j,j2)-Oij(i)*
     &                    uu2(jj,i,j)*uu2(jj,i,j2)
                     
                  END DO
               END DO

c     formula (3)
               Mi=Ai*
     &              kchoos1*sqi1**lowerk1*(1-sqi1)**(upperK-lowerk1)* 
     &              kchoos2*sqi2**lowerk2*(1-sqi2)**(upperK-lowerk2) 

            END DO              
            
c     products for derivatives of 1st order ( formulas (10)-(11) )
            DO j=1,total1
               Sbetaj(j)=Sbetaj(j)+Mi*Niz(j)
            END DO

            DO j=1,total2a
               Ssigmn1(j)=Ssigmn1(j)+v1*Mi*Niuu1(j)
            END DO

            DO j=1,total2b
               Ssigmn2(j)=Ssigmn2(j)+v2*Mi*Niuu2(j)
            END DO

c     summations for derivatives of second order :
            
            DO j=1,total1
               DO j2=1,total1
c     for the formula (14)
                  S2bjbl(j,j2)=S2bjbl(j,j2)+Mi*
     &                 (Niz(j)*Niz(j2)-Ozz(j,j2))
               END DO
               DO j2=1,total2a
c     for the formula (15)
                  S2bjsn1(j,j2)=S2bjsn1(j,j2)+v1*Mi*(Niz(j)*
     &                 Niuu1(j2)-Ozuu1(j,j2))
               END DO
               DO j2=1,total2b
                  S2bjsn2(j,j2)=S2bjsn2(j,j2)+v2*Mi*(Niz(j)*
     &                 Niuu2(j2)-Ozuu2(j,j2))
               END DO
            END DO  

c     for the formula (17)
            DO j=1,total2a
               DO j2=1,total2a
                  S2snso1(j,j2)=S2snso1(j,j2)+v1*v1*Mi*
     &		       (Niuu1(j)* Niuu1(j2)-Ouuu1(j,j2))
               END DO
               DO j2=1,total2b
                  S2snso12(j,j2)=S2snso12(j,j2)+v1*v2*
     &                 Mi*(Niuu1(j)* Niuu2(j2)-Ouuu12(j,j2)) 		     
               END DO
            END DO 

            DO j=1,total2b
               DO j2=1,total2b
                  S2snso2(j,j2)=S2snso2(j,j2)+v2*v2*Mi*
     &		       (Niuu2(j)* Niuu2(j2)-Ouuu2(j,j2))
               END DO
            END DO
            
            Li=Li+Mi            
         END DO
      END DO                    

      RETURN
      End                       

c ======================================================================

* ======================================================================
* Function "FACTOR" calculates factorial of the number that comes
* in as a PARAMETER.
* Input  : n       - non-negative number which factorial should be taken
* Output : FACTOR  - equal to n!
* ======================================================================

c     Begin Function FACTOR()

      DOUBLE PRECISION function FACTOR(n1)

      implicit none
	   INTEGER n1
           DOUBLE PRECISION nfact
           INTEGER i

           nfact=1           
           IF (n1 .gt. 0) then 
              DO i=1,n1
                 nfact=nfact*i
              END DO	        
	   END IF  
	   FACTOR = nfact 
           RETURN

      End 
c     Function FACTOR     End Function FACTOR

c     ======================================================================

c     ======================================================================
*     INVERT: This inverts subroutine will inverts a matrix with 
*     leading dimension ld and actual dimension N.
*     
*     This was written by Geoff Vining (UF)
c     ======================================================================
      
      
      SUBROUTINE INVERT(A,AINV,n)
c     Begin Subroutine INVERT

      implicit none
      INTEGER maxsig,maxbet
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      
      integer n,i,j,i2,n2,i3
      real(8) A(n,n),AINV(n,n)
      real(8) MAX,TEMP,fact
      
      DO I = 1,N
         DO J = 1,N
            AINV(I,J) = 0
         END DO
         AINV(I,I) = 1.0
      END DO

      DO I = 1,N-1
         MAX = DABS(A(I,I))
         DO I2 = I+1,N
            IF (DABS(A(I2,I)) .gt. MAX) THEN
               MAX = DABS(A(I2,I))
               DO J = 1,N
                  TEMP = A(I,J)
                  A(I,J) = A(I2,J)
                  A(I2,J) = TEMP
                  TEMP = AINV(I,J)
                  AINV(I,J) = AINV(I2,J)
                  AINV(I2,J) = TEMP
               END DO
            END IF
         END DO	                      
         MAX = A(I,I)
         DO I2 = I+1,N
            FACT = A(I2,I)/MAX
            DO J = 1,N
               A(I2,J) = A(I2,J) - FACT*A(I,J)
               AINV(I2,J) = AINV(I2,J) - FACT*AINV(I,J)
            END DO
         END DO
      END DO

      DO I2 = 1,N
         MAX = A(I2,I2)
         IF (MAX .eq. 0) THEN
C            PRINT *,'MATRIX IS SINGULAR !' 
            return
         END IF
         IF (DABS(MAX) .lt. 0.000000001) THEN 
C            PRINT *,'MATRIX IS NEAR SINGULAR !'
         END IF
         
         DO J = 1,N
            A(I2,J) = A(I2,J)/MAX
            AINV(I2,J) = AINV(I2,J)/MAX
         END DO
      END DO

      N2 = N+1
      DO I = 1,N-1
         I2 = N2-I
         AINV(I2,I2) = AINV(I2,I2)/A(I2,I2)
         A(I2,I2) = 1.0
         DO I3 = 1,I2-1
            DO J = 1,N
               AINV(I3,J)=AINV(I3,J)-A(I3,I2)*AINV(I2,J)
            END DO
         END DO
      END DO 

      AINV(1,1) = AINV(1,1)/A(1,1)
      RETURN
      End
c     Subroutine INVERT  End subroutine INVERT

c=========================================================================

c=========================================================================
*     PURPOSE:    SUBROUTINE CONMIN MINIMIZES AN UNCONSTRAINED NONLINEAR
*     SCALAR VALUED FUNCTION OF A VECTOR VARIABLE X
*     EITHER BY THE BFGS VARIABLE METRIC ALGORITHM OR BY A
*     BEALE RESTARTED CONJUGATE GRADIENT ALGORITHM.
*     
*     USAGE:      CALL CONMIN(N,X,F,G,IFUN,ITER,EPS,NFLAG,MXFUN,W,
*     IOUT,MDIM,IDEV,ACC,NMETH)
*     
*     PARAMETERS: N      THE NUMBER OF VARIABLES IN THE FUNCTION TO
*     BE MINIMIZED.
*     X      THE VECTOR CONTAINING THE CURRENT ESTIMATE TO
*     THE MINIMIZER. ON ENTRY TO CONMIN,X MUST CONTAIN
*     AN INITIAL ESTIMATE SUPPLIED BY THE USER.
*     ON EXITING,X WILL HOLD THE BEST ESTIMATE TO THE
*     MINIMIZER OBTAINED BY CONMIN. X MUST BE DOUBLE
*     PRECISIONED AND DIMENSIONED N.
*     F      ON EXITING FROM CONMIN,F WILL CONTAIN THE LOWEST
*     VALUE OF THE OBJECT FUNCTION OBTAINED.
*     F IS DOUBLE PRECISIONED.
*     G      ON EXITING FROM CONMIN,G WILL CONTAIN THE
*     ELEMENTS OF THE GRADIENT OF F EVALUATED AT THE
*     POINT CONTAINED IN X. G MUST BE DOUBLE
*     PRECISIONED AND DIMENSIONED N.
*     IFUN   UPON EXITING FROM CONMIN,IFUN CONTAINS THE
*     NUMBER OF TIMES THE FUNCTION AND GRADIENT
*     HAVE BEEN EVALUATED.
*     ITER   UPON EXITING FROM CONMIN,ITER CONTAINS THE
*     TOTAL NUMBER OF SEARCH DIRECTIONS CALCULATED
*     TO OBTAIN THE CURRENT ESTIMATE TO THE MINIZER.
*     EPS    EPS IS THE USER SUPPLIED CONVERGENCE PARAMETER.
*     CONVERGENCE OCCURS WHEN THE NORM OF THE GRADIENT
*     IS LESS THAN OR EQUAL TO EPS TIMES THE MAXIMUM
*     OF ONE AND THE NORM OF THE VECTOR X. EPS
*     MUST BE DOUBLE PRECISIONED.
*     NFLAG  UPON EXITING FROM CONMIN,NFLAG STATES WHICH
*     CONDITION CAUSED THE EXIT.
*     IF NFLAG=0, THE ALGORITHM HAS CONVERGED.
*     IF NFLAG=1, THE MAXIMUM NUMBER OF FUNCTION
*     EVALUATIONS HAVE BEEN USED.
*     IF NFLAG=2, THE LINEAR SEARCH HAS FAILED TO
*     IMPROVE THE FUNCTION VALUE. THIS IS THE
*     USUAL EXIT IF EITHER THE FUNCTION OR THE
*     GRADIENT IS INCORRECTLY CODED.
*     IF NFLAG=3, THE SEARCH VECTOR WAS NOT
*     A DESCENT DIRECTION. THIS CAN ONLY BE CAUSED
*     BY ROUNDOFF,AND MAY SUGGEST THAT THE
*     CONVERGENCE CRITERION IS TOO STRICT.
*     MXFUN  MXFUN IS THE USER SUPPLIED MAXIMUM NUMBER OF
*     FUNCTION AND GRADIENT CALLS THAT CONMIN WILL
*     BE ALLOWED TO MAKE.
*     W      W IS A VECTOR OF WORKING STORAGE.IF NMETH=0,
*     W MUST BE DIMENSIONED 5*N+2. IF NMETH=1,
*     W MUST BE DIMENSIONED N*(N+7)/2. IN BOTH CASES,
*     W MUST BE DOUBLE PRECISIONED.
*     IOUT   IOUT IS A USER  SUPPLIED OUTPUT PARAMETER.
*     IF IOUT = 0, THERE IS NO PRINTED OUTPUT FROM
*     CONMIN. IF IOUT > 0,THE VALUE OF F AND THE
*     NORM OF THE GRADIENT SQUARED,AS WELL AS ITER
*     AND IFUN,ARE WRITTEN EVERY IOUT ITERATIONS.
*     MDIM   MDIM IS THE USER SUPPLIED DIMENSION OF THE
*     VECTOR W. IF NMETH=0,MDIM=5*N+2. IF NMETH=1,
*     MDIM=N*(N+7)/2.
*     IDEV   IDEV IS THE USER SUPPLIED NUMBER OF THE OUTPUT
*     DEVICE ON WHICH OUTPUT IS TO BE WRITTEN WHEN
*     IOUT>0.
*     ACC    ACC IS A USER SUPPLIED ESTIMATE OF MACHINE
*     ACCURACY. A LINEAR SEARCH IS UNSUCCESSFULLY
*     TERMINATED WHEN THE NORM OF THE STEP SIZE
*     BECOMES SMALLER THAN ACC. IN PRACTICE,
*     ACC=10.D-20 HAS PROVED SATISFACTORY. ACC IS
*     DOUBLE PRECISIONED.
*     NMETH  NMETH IS THE USER SUPPLIED VARIABLE WHICH
*     CHOOSES THE METHOD OF OPTIMIZATION. IF
*     NMETH=0,A CONJUGATE GRADIENT METHOD IS
*     USED. IF NMETH=1, THE BFGS METHOD IS USED.
*     
*     REMARKS:    IN ADDITION TO THE SPECIFIED VALUES IN THE ABOVE
*     ARGUMENT LIST, THE USER MUST SUPPLY A SUBROUTINE
*     CALCFG WHICH CALCULATES THE FUNCTION AND GRADIENT AT
*     X AND PLACES THEM IN F AND G(1),...,G(N) RESPECTIVELY.
*     THE SUBROUTINE MUST HAVE THE FORM:
*     SUBROUTINE CALCFG(N,X,F,G)
*     DOUBLE PRECISION X(N),G(N),F
*     
*     AN EXAMPLE SUBROUTINE FOR THE ROSENBROCK FUNCTION IS:
*     
*     SUBROUTINE CALCFG(N,X,F,G)
*     DOUBLE PRECISION X(N),G(N),F,T1,T2
*     T1=X(2)-X(1)*X(1)
*     T2=1.0-X(1)
*     F=100.0*T1*T1+T2*T2
*     G(1)=-400.0*T1*X(1)-2.0*T2
*     G(2)=200.0*T1
*     RETURN
*     END
*     
c     =======================================================================

c     Begin Subroutine CONMIN	

      SUBROUTINE CONMIN(upk_in,X,F,G,hess,IFUN,ITER,EPS,NFLAG,MXFUN,W,
     &     IOUT,MDIM,IDEV,ACC,NMETH,ri,z,uu1,
     &     uu2,total1,cg,total2a,total2b,total3,numcas)

      implicit none

c     Declaration Section

c     Constant Declaration
      INTEGER maxsub,maxcas,maxsig,maxest,maxbet,max_w
      PARAMETER (maxsub=5200)
      PARAMETER (maxcas=10)
      PARAMETER (maxsig=10)
      PARAMETER (maxbet=25)
      PARAMETER (maxest=maxsig+maxbet)
      PARAMETER (max_w=maxest*(maxest+7)/2)
      
c     Logical Declaration
      LOGICAL RSW

c     Integer Declaration

c     Integer Variables
      INTEGER total1,total2a,total2b,total3
      INTEGER i, n, ii, j, ij
      INTEGER ifun, iter, iout, idev, ioutk
      INTEGER nflag, nmeth, nx, ng, ncons, nry, ncons1, ncons2
      INTEGER nrst, ncalls, nxpi, ngpi, nrdpi, nrypi, ngpj, nrd
      INTEGER mdim, upk_in,cg ,mxfun

c     Integer Array
      INTEGER numcas(maxsub)

c     Double Precision Declaration

c     Double Precision Variables	   
      DOUBLE PRECISION F,FP,FMIN,ALPHA,AT,AP,GSQ,DG,DG1
      DOUBLE PRECISION DP,STEP,ACC,DAL,U1,U2,U3,U4,EPS
      DOUBLE PRECISION XSQ,RTST,DSQRT,DMIN1,DMAX1,DABS
      
c     One Dim Array
      DOUBLE PRECISION x(total1+total2a+total2b) ,g(maxest)
      DOUBLE PRECISION w(max_w)

c     Two Dim Array
      DOUBLE PRECISION hess(total1+total2a+total2b,
     +     total1+total2a+total2b)
      DOUBLE PRECISION ri(maxsub,maxcas)

c     Three Dim Array
      DOUBLE PRECISION z(maxsub,maxcas,maxbet)
      DOUBLE PRECISION uu1(maxsub,maxcas,maxsig),
     +     uu2(maxsub,maxcas,maxsig)          
      
c     End Declaration Section

      
      n=total1+total2a+total2b

c     Initialize ITER,IFUN,NFLAG, and IOUTK, which counts output iterations 
      ITER=0		
      IFUN=0
      IOUTK=0
      NFLAG=0

c     SET PARAMETERS TO EXTRACT VECTORS FROM W.
c     W(I) HOLDS THE SEARCH VECTOR,W(NX+I) HOLDS THE BEST CURRENT
c     ESTIMATE TO THE MINIMIZER,AND W(NG+I) HOLDS THE GRADIENT
c     AT THE BEST CURRENT ESTIMATE.

      NX=N
      NG=NX+N

c     TEST WHICH METHOD IS BEING USED.
c     IF NMETH=0, W(NRY+I) HOLDS THE RESTART Y VECTOR AND
c     W(NRD+I) HOLDS THE RESTART SEARCH VECTOR.

      IF (NMETH.EQ.1) THEN 
         NCONS = 3 * N
c     If NMETH=1, W(NCONS+I) holds the appr. inverse HESSIAN
      ELSE 
         NRY=NG+N
         NRD=NRY+N
         NCONS=5*N
         NCONS1=NCONS+1
         NCONS2=NCONS+2
      END IF
      
c     CALCULATE THE FUNCTION AND GRADIENT AT THE INITIAL
c     POINT AND INITIALIZE NRST,WHICH IS USED TO DETERMINE
c     WHETHER A BEALE RESTART IS BEING DONE. NRST=N MEANS THAT THIS
c     ITERATION IS A RESTART ITERATION. INITIALIZE RSW,WHICH INDICATES
c     THAT THE CURRENT SEARCH DIRECTION IS A GRADIENT DIRECTION.

 20   CALL CALCFG(upk_in,X,total1,cg,total2a,total2b,total3,z,
     &     uu1,uu2,ri,numcas,F,G,HESS)
      IFUN = IFUN+1
      NRST = N
      RSW = .TRUE.

c     CALCULATE THE INITIAL SEARCH DIRECTION , THE NORM OF X SQUARED,
c     AND THE NORM OF G SQUARED. DG1 IS THE CURRENT DIRECTIONAL
c     DERIVATIVE,WHILE XSQ AND GSQ ARE THE SQUARED NORMS.

      DG1 = 0.
      XSQ = 0.
      DO I = 1,N
         W(I) = -G(I)
         XSQ = XSQ + X(I) * X(I)
         DG1 = DG1 - G(I) * G(I)
      END DO
      GSQ = -DG1

c     TEST IF THE INITIAL POINT IS THE MINIMIZER.
      IF (GSQ .le.  EPS*EPS*DMAX1(1.0D0,XSQ)) THEN
         RETURN
      END IF

c     BEGIN THE MAJOR ITERATION LOOP. NCALLS IS USED TO GUARANTEE THAT
c     AT LEAST TWO POINTS HAVE BEEN TRIED WHEN NMETH=0. FMIN IS THE
c     CURRENT FUNCTION VALUE.

 40   FMIN=F
      NCALLS=IFUN

c     IF OUTPUT IS DESIRED,TEST IF THIS IS THE CORRECT ITERATION
c     AND IF SO, WRITE OUTPUT.

      IF (IOUT .eq. 0) THEN
         ALPHA = ALPHA * DG / DG1
c     Set ALPHA to nonrestart conjugate gadient
      ELSE IF (IOUTK .ne. 0) THEN
         IOUTK = IOUTK + 1
         IF (IOUTK .eq. IOUT) THEN
            IOUTK = 0
         END IF
         ALPHA = ALPHA * DG / DG1   	        
      ELSE
C         WRITE(IDEV,50)ITER,IFUN,FMIN,GSQ
 50      FORMAT(10H ITERATION,I5,20H      FUNCTION CALLS,I6/5H F = 
     &        ,D15.8,13H G-SQUARED = ,D15.8/)  
C         WRITE(IDEV,60)(X(I),I=1,total1+total2a+total2b)
 60      FORMAT(/8HINTER X./1H ,20D16.8)
      END IF

c     IF NMETH=1 OR A RESTART HAS BEEN PERFORMED, SET ALPHA=1.0.
      IF (NRST .eq. 1.OR.NMETH .eq. 1) THEN
         ALPHA=1.0
      END IF

c     IF A GRADIENT DIRECTION IS USED, SET ALPHA=1.0/DSQRT(GSQ),
c     WHICH SCALES THE INITIAL SEARCH VECTOR TO UNITY.
      IF (RSW) THEN
         ALPHA=1.0/DSQRT(GSQ)
      END IF

c     THE LINEAR SEARCH FITS A CUBIC TO F AND DAL, THE FUNCTION AND ITS
c     DERIVATIVE AT ALPHA, AND TO FP AND DP,THE FUNCTION
c     AND DERIVATIVE AT THE PREVIOUS TRIAL POINT AP.
c     INITIALIZE AP ,FP,AND DP.

      AP=0.
      FP=FMIN
      DP=DG1

c     SAVE THE CURRENT DERIVATIVE TO SCALE THE NEXT SEARCH VECTOR.
      DG=DG1

c     UPDATE THE ITERATION.
      ITER=ITER+1

c     CALCULATE THE CURRENT STEPLENGTH  AND STORE THE CURRENT X AND G.
      STEP=0.
      DO I=1,N
         STEP=STEP+W(I)*W(I)
         NXPI=NX+I
         NGPI=NG+I
         W(NXPI)=X(I)
         W(NGPI)=G(I)
      END DO
      STEP=DSQRT(STEP)

c     BEGIN THE LINEAR SEARCH ITERATION.
c     TEST FOR FAILURE OF THE LINEAR SEARCH.

 80   IF (ALPHA*STEP .le.  ACC) THEN
c     TEST IF DIRECTION IS A GRADIENT DIRECTION.
         IF (.NOT.RSW) THEN
            GO TO 20
c     Call subroutine CALCFG           
         ELSE 	
            NFLAG=2
            RETURN
         END IF
      END IF

c     CALCULATE THE TRIAL POINT.
      DO I = 1,N
         NXPI = NX + I
         X(I) = W(NXPI) + ALPHA * W(I)
      END DO 

c     EVALUATE THE FUNCTION AT THE TRIAL POINT.
c     Call CALCFG
      CALL CALCFG(upk_in,X,total1,cg,total2a,total2b,total3,z,
     &     uu1,uu2,ri,numcas,F,G,HESS)

c     TEST IF THE MAXIMUM NUMBER OF FUNCTION CALLS HAVE BEEN USED.
      IFUN=IFUN+1
      IF(IFUN .gt. MXFUN) THEN 
         NFLAG=1
         RETURN
      END IF

c     COMPUTE THE DERIVATIVE OF F AT ALPHA.
      DAL=0.0
      DO I=1,N
         DAL=DAL+G(I)*W(I)
      END DO

c     TEST WHETHER THE NEW POINT HAS A NEGATIVE SLOPE BUT A HIGHER
c     FUNCTION VALUE THAN ALPHA=0. IF THIS IS THE CASE,THE SEARCH
c     HAS PASSED THROUGH A LOCAL MAX AND IS HEADING FOR A DISTANT LOCAL
c     MINIMUM.
      IF (F .gt. FMIN .AND. DAL .lt. 0.) GO TO 160

c     IF NOT, TEST WHETHER THE STEPLENGTH CRITERIA HAVE BEEN MET.
      IF(F .gt. (FMIN+.0001*ALPHA*DG) .OR. DABS(DAL/DG)
     &     .gt. (.9)) GO TO 130

c     IF THEY HAVE BEEN MET, TEST IF TWO POINTS HAVE BEEN TRIED
c     IF NMETH=0 AND IF THE TRUE LINE MINIMUM HAS NOT BEEN FOUND.
      IF ((IFUN-NCALLS).le. 1 .AND. DABS(DAL/DG) .gt.  EPS .AND.
     &     NMETH .eq. 0) THEN 
         GO TO 130
      ELSE
         GO TO 170
      END IF

c     A NEW POINT MUST BE TRIED. USE CUBIC INTERPOLATION TO FIND
c     THE TRIAL POINT AT.
 130  U1=DP+DAL-3.0*(FP-F)/(AP-ALPHA)
      U2=U1*U1-DP*DAL
      IF(U2.LT.0.)U2=0.
      U2=DSQRT(U2)
      AT=ALPHA-(ALPHA-AP)*(DAL+U2-U1)/(DAL-DP+2.*U2)


c     TEST WHETHER THE LINE MINIMUM HAS BEEN BRACKETED.
      IF((DAL/DP).GT.0.)GO TO 140

c     THE MINIMUM HAS BEEN BRACKETED. TEST WHETHER THE TRIAL POINT LIES
c     SUFFICIENTLY WITHIN THE BRACKETED INTERVAL.
c     IF IT DOES NOT, CHOOSE AT AS THE MIDPOINT OF THE INTERVAL.

      IF(AT.LT.(1.01*DMIN1(ALPHA,AP)).OR.AT.GT.(.99*DMAX1
     &     (ALPHA,AP)))AT=(ALPHA+AP)/2.0
      GO TO 150

c     THE MINIMUM HAS NOT BEEN BRACKETED. TEST IF BOTH POINTS ARE
c     GREATER THAN THE MINIMUM AND THE TRIAL POINT IS SUFFICIENTLY
c     SMALLER THAN EITHER.

 140  IF (DAL .GT.0.0.AND.0.0.LT.AT.AND.AT.LT.
     &     (.99*DMIN1(AP,ALPHA))) GO TO 150

c     TEST IF BOTH POINTS ARE LESS THAN THE MINIMUM AND THE TRIAL POINT
c     IS SUFFICIENTLY LARGE.
      IF(DAL.LE.0.0.AND.AT.GT.(1.01*DMAX1(AP,ALPHA)))GO TO 150

c     IF THE TRIAL POINT IS TOO SMALL,DOUBLE THE LARGEST PRIOR POINT.
      IF(DAL.LE.0.)AT=2.0*DMAX1(AP,ALPHA)

c     IF THE TRIAL POINT IS TOO LARGE, HALVE THE SMALLEST PRIOR POINT.
      IF(DAL.GT.0.)AT=DMIN1(AP,ALPHA)/2.0

c     SET AP=ALPHA, ALPHA=AT,AND CONTINUE SEARCH.
 150  AP=ALPHA
      FP=F
      DP=DAL
      ALPHA=AT
      GO TO 80

c     A RELATIVE MAX HAS BEEN PASSED.REDUCE ALPHA AND RESTART THE SEARCH.
 160  ALPHA=ALPHA/3.
      AP=0.
      FP=FMIN
      DP=DG
      GO TO 80

c     THE LINE SEARCH HAS CONVERGED. TEST FOR CONVERGENCE OF THE ALGORITHM.
 170  GSQ=0.0
      XSQ=0.0
      DO I=1,N
         GSQ=GSQ+G(I)*G(I)
         XSQ=XSQ+X(I)*X(I)
      END DO

      IF (GSQ  .le.  EPS*EPS*DMAX1(1.0D0,XSQ)) THEN
         RETURN
      END IF

c     SEARCH CONTINUES. SET W(I)=ALPHA*W(I),THE FULL STEP VECTOR.
      DO I=1,N
         W(I)=ALPHA*W(I)
      END DO

c     COMPUTE THE NEW SEARCH VECTOR. FIRST TEST WHETHER A
c     CONJUGATE GRADIENT OR A VARIABLE METRIC VECTOR IS USED.
      IF (NMETH .ne. 1) THEN
c     Begin if nmeth /= 1
c     CONJUGATE GRADIENT UPDATE SECTION.
c     TEST IF A POWELL RESTART IS INDICATED.
         RTST=0.
         DO I=1,N
            NGPI=NG+I
            RTST=RTST+G(I)*W(NGPI)
         END DO
         IF (DABS(RTST/GSQ) .gt. 0.2) THEN
            NRST=N
         END IF
c     IF A RESTART IS INDICATED, SAVE THE CURRENT D AND Y
c     AS THE BEALE RESTART VECTORS AND SAVE D'Y AND Y'Y
c     IN W(NCONS+1) AND W(NCONS+2).
         IF (NRST .eq. N) THEN 
            W(NCONS+1)=0.
            W(NCONS+2)=0.
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
c     CALCULATE  THE RESTART HESSIAN TIMES THE CURRENT GRADIENT.
         U1=0.0
         U2=0.0
         DO I=1,N
            NRDPI=NRD+I
            NRYPI=NRY+I
            U1=U1-W(NRDPI)*G(I)/W(NCONS1)
            U2=U2+W(NRDPI)*G(I)*2./W(NCONS2)-
     &           W(NRYPI)*G(I)/W(NCONS1)
         END DO
         U3 = W(NCONS2)/W(NCONS1)
         DO I=1,N
            NXPI=NX+I
            NRDPI=NRD+I
            NRYPI=NRY+I
            W(NXPI)=-U3*G(I)-U1*W(NRYPI)-U2*W(NRDPI)
         END DO
c     IF THIS IS A RESTART ITERATION,W(NX+I) CONTAINS THE NEW SEARCH
c     VECTOR.
         IF (NRST .ne. N) THEN
c     begin if nrst /= n
c     NOT A RESTART ITERATION. CALCULATE THE RESTART HESSIAN
c     TIMES THE CURRENT Y.
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
     &              +2.0*W(NRDPI)*(G(I)-W(NGPI))/W(NCONS2)
               U3=U3+W(I)*(G(I)-W(NGPI))
            END DO
            STEP=0.
            DO I=1,N
               NGPI=NG+I
               NRDPI=NRD+I
               NRYPI=NRY+I
               STEP=(W(NCONS2)/W(NCONS1))*(G(I)-W(NGPI))
     &              +U1*W(NRYPI)+U2*W(NRDPI)
               U4=U4+STEP*(G(I)-W(NGPI))
               W(NGPI)=STEP
            END DO

c     CALCULATE THE DOUBLY UPDATED HESSIAN TIMES THE CURRENT
c     GRADIENT TO OBTAIN THE SEARCH VECTOR.
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
c     CALCULATE THE DERIVATIVE ALONG THE NEW SEARCH VECTOR.
         END IF
c     End if nrst /= n
         DG1=0.
         DO I=1,N
            NXPI=NX+I
            W(I)=W(NXPI)
            DG1=DG1+W(I)*G(I)
         END DO
c     IF THE NEW DIRECTION IS NOT A DESCENT DIRECTION,STOP.
         IF (DG1 .gt. 0.) THEN
c     GO TO 320
            NFLAG = 3
            RETURN
         END IF 

c     UPDATE NRST TO ASSURE AT LEAST ONE RESTART EVERY N ITERATIONS.
         IF (NRST .eq. N) NRST=0
         NRST=NRST+1
         RSW=.FALSE.
         GO TO 40	      
c     A VARIABLE METRIC ALGORITM IS BEING USED. CALCULATE Y AND D'Y.
      END IF
c     End if nmeth /= 1

      U1=0.0
      DO I=1,N
         NGPI=NG+I
         W(NGPI)=G(I)-W(NGPI)
         U1=U1+W(I)*W(NGPI)
      END DO

c     IF RSW=.TRUE.,SET UP THE INITIAL SCALED APPROXIMATE HESSIAN.
      IF (RSW) THEN
c     CALCULATE Y'Y.
         U2=0.
         DO I=1,N
            NGPI=NG+I
            U2=U2+W(NGPI)*W(NGPI)
         END DO
c     CALCULATE THE INITIAL HESSIAN AS H=(P'Y/Y'Y)*I
c     AND THE INITIAL U2=Y'HY AND W(NX+I)=HY.
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
c     CALCULATE W(NX+I)=HY AND U2=Y'HY.
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

c     CALCULATE THE UPDATED APPROXIMATE HESSIAN.
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

c     CALCULATE THE NEW SEARCH DIRECTION W(I)=-HG AND ITS DERIVATIVE.
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

c     TEST FOR A DOWNHILL DIRECTION.
      IF (DG1 .gt. 0.) THEN      
         NFLAG = 3
         RETURN
      ELSE
         RSW=.FALSE.
         GO TO 40
      END IF    


      End                 
c     Subroutine CONMIN         End subroutine CONMIN
c ========================================================================
c     The End
