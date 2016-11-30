/*
 * calcs1.c
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#include <memory.h>
#include <math.h>
#include "defs.h"
#include "gaps.h"
#include "calcs.h"
#include "calcs1.h"

double
L1( RECURSE_PARAMS *pR )
{
    return( pR->dPrdA * 
            exp( -pR->dSumB * exp( pR->dParams[P_BETA0] )) * 
            exp( -pR->dSumC * exp( pR->dParams[P_PHI] )));
}

void
LogLikelihood1( double *pdParameters, double *padReturn, int *piError )
{
    BRANCH_SUM aBranches[1];
    long lSubject;
    double dValue;
    RECURSE_PARAMS RParams;
    
    padReturn[0] = 0.0;

    *piError = ERR_NONE;

    if ( glNumSubjects == 0 ) {
        *piError = ERR_BAD_DATA;
        return;
    }

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 1;
    RParams.pfnEquation[0] = L1;
    RParams.lCovariateX = 0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;

        CalcRecurse( &RParams, aBranches );

        dValue = SumNodes( aBranches );
        if ( dValue > 0.0 )
            padReturn[0] += log( dValue );
        else {
            padReturn[0] = -Inf;
            break;
        }
    }
}

double
S1Beta( RECURSE_PARAMS *pR )
{
    return( -pR->dSumB * exp( pR->dParams[P_BETA0] ) *
            L1( pR ));
}

double
S1Phi( RECURSE_PARAMS *pR )
{
    return( -pR->dSumC * exp( pR->dParams[P_PHI] ) *
            L1( pR ));
}

void
ScoreVector1( double *pdParameters, double *padReturn )
{
    double adSum[3];
    BRANCH_SUM aBranches[3];
    double *padLogLik, *padScoreVector;
    long lSubject;
    RECURSE_PARAMS RParams;
    int i;
    
    if ( glNumSubjects == 0 )
        return;

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 3;
    padScoreVector = adSum;
    RParams.pfnEquation[0] = S1Beta;
    RParams.pfnEquation[1] = S1Phi;
    padLogLik = adSum + 2;
    RParams.pfnEquation[2] = L1;
    RParams.lCovariateX = 0;

    for ( i = 0; i < 2; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        for ( i = 0; i < 2; i++ )
            padReturn[i] += padScoreVector[i] / padLogLik[0];
    }
}

double
dL1_dBeta_dBeta( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumB * exp( pR->dParams[P_BETA0] );
    return( dTmp * ( dTmp - 1.0 ) * L1( pR ));
}

double
dL1_dBeta_dPhi( RECURSE_PARAMS *pR )
{
    return( pR->dSumB * pR->dSumC * 
            exp( pR->dParams[P_BETA0] + pR->dParams[P_PHI] ) *
            L1( pR ));
}

double
dL1_dPhi_dPhi( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumC * exp( pR->dParams[P_PHI] );
    return( dTmp * ( dTmp - 1.0 ) * L1( pR ));
}

void
Hessian1( double *pdParameters, double *padReturn )
{
    double adSum[6];
    BRANCH_SUM aBranches[6];
    double *padHessian, *padLogLik, *padScoreVector;
    long lSubject;
    RECURSE_PARAMS RParams;
    int i;
    
    if ( glNumSubjects == 0 )
        return;

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 6;
    padHessian = adSum;
    RParams.pfnEquation[0] = dL1_dBeta_dBeta;
    RParams.pfnEquation[1] = dL1_dBeta_dPhi;
    RParams.pfnEquation[2] = dL1_dPhi_dPhi;
    padScoreVector = adSum + 3;
    RParams.pfnEquation[3] = S1Beta;
    RParams.pfnEquation[4] = S1Phi;
    padLogLik = adSum + 5;
    RParams.pfnEquation[5] = L1;
    RParams.lCovariateX = 0;

    for ( i = 0; i < 4; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        /* dL1_dBeta_dBeta */
        padReturn[0] += ( padHessian[0] * padLogLik[0] - 
                          padScoreVector[0] * padScoreVector[0] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL1_dBeta_dPhi */
        padReturn[1] += ( padHessian[1] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL1_dPhi_dPhi */
        padReturn[3] += ( padHessian[2] * padLogLik[0] -
                          padScoreVector[1] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
    }
    /* dL1_dPhi_dBeta */
    padReturn[2] = padReturn[1];
}
