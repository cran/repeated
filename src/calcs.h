/*
 * calcs.h
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1996 by Richard J. Cook, University of Waterloo
 */

#ifndef __CALCS_H__
#define __CALCS_H__


/*
#define MORE_PRECISE
#define UNSIGNED_METHOD
*/


#define pow dPow

#define MAX_PARAMS  5
#define P_BETA0     0
#define P_BETA1     1
#define P_PHI       2
#define P_DELTA     3
#define P_THETA     4
#define MAX_HESSIAN_FCNS   ( MAX_PARAMS * ( MAX_PARAMS + 1 )) / 2
#define MAX_SCORE_FCNS      MAX_PARAMS
#define MAX_LIKELIHOOD_FCNS 1
#define MAX_FCNS            ( MAX_HESSIAN_FCNS + MAX_SCORE_FCNS + MAX_LIKELIHOOD_FCNS )

#define ERR_NONE                0
#define ERR_BAD_DATA            100
#define ERR_DOUBLE_PRECISION    200

typedef struct _RECURSE_PARAMS {
    double  dParams[MAX_PARAMS];
    long    lCovariateX;
    long    lSubjectID;
    double  dPrdA;
    double  dSumB;
    double  dSumC;
    long    lL;
    int     iNumEqns;
    double  ( * pfnEquation[MAX_FCNS] )( struct _RECURSE_PARAMS * );
} RECURSE_PARAMS;

typedef double ( * FN_EQUATION )( RECURSE_PARAMS * );

#if !defined( MORE_PRECISE )

typedef struct _BRANCH_SUM {
    double dPartialSum;
    double dSum;
} BRANCH_SUM;

#else

typedef struct _BRANCH_SUM {
    int iAllocSize;
    int iNumElements;
    double *padValues;
} BRANCH_SUM;

#endif

extern double dPow( double, double );
extern double dChoose( long, long );
extern void CalcWorkTimeFactors( void );
void CalcRecurse( RECURSE_PARAMS *, BRANCH_SUM * );
double FcnAsubL( long, long, long, long, long, long );
double FcnBsubL( long, long, long, long, long, long );
double FcnCsubL( long, long, long, long, long, long );

int compare( const void *, const void * );
void AddNode( BRANCH_SUM *, double );
double SumNodes( BRANCH_SUM * );

#endif
