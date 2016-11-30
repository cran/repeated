/*
 * calcs2.c
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
#include "calcs2.h"

double
L2( RECURSE_PARAMS *pR )
{
    return( pR->dPrdA * 
            exp( -pR->dSumB * exp( pR->dParams[P_BETA0] +
                                   pR->dParams[P_BETA1] *
                                   pR->lCovariateX )) *
            exp( -pR->dSumC * exp( pR->dParams[P_PHI] )));
}

void
LogLikelihood2( double *pdParameters, double *padReturn, int *piError )
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
    RParams.pfnEquation[0] = L2;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;
        RParams.lCovariateX = gaSubjects[lSubject].lCovariateX[0];

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
S2Beta0( RECURSE_PARAMS *pR )
{
    return( -pR->dSumB * 
            exp( pR->dParams[P_BETA0] +
                 pR->dParams[P_BETA1] *
                 pR->lCovariateX ) *
            L2( pR ));
}

double
S2Beta1( RECURSE_PARAMS *pR )
{
    return( S2Beta0( pR ) *
            pR->lCovariateX );
}

double
S2Phi( RECURSE_PARAMS *pR )
{
    return( -pR->dSumC * exp( pR->dParams[P_PHI] ) *
            L2( pR ));
}

void
ScoreVector2( double *pdParameters, double *padReturn )
{
    double adSum[4];
    BRANCH_SUM aBranches[4];
    double *padLogLik, *padScoreVector;
    long lSubject;
    RECURSE_PARAMS RParams;
    int i;
    int iNumDerivatives;
    
    if ( glNumSubjects == 0 )
        return;

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 4;
    padScoreVector = adSum;
    iNumDerivatives = 3;
    RParams.pfnEquation[0] = S2Beta0;
    RParams.pfnEquation[1] = S2Beta1;
    RParams.pfnEquation[2] = S2Phi;
    padLogLik = adSum + iNumDerivatives;
    RParams.pfnEquation[3] = L2;

    for ( i = 0; i < iNumDerivatives; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;
        RParams.lCovariateX = gaSubjects[lSubject].lCovariateX[0];

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        for ( i = 0; i < iNumDerivatives; i++ )
            padReturn[i] += padScoreVector[i] / padLogLik[0];
    }
}

double
dL2_dBeta0_dBeta0( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumB * 
                  exp( pR->dParams[P_BETA0] +
                       pR->dParams[P_BETA1] * 
                       pR->lCovariateX );
    return( dTmp * ( dTmp - 1.0 ) * L2( pR ));
}

double
dL2_dBeta0_dBeta1( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumB * 
                  pR->lCovariateX *
                  exp( pR->dParams[P_BETA0] +
                       pR->dParams[P_BETA1] * 
                       pR->lCovariateX );
    return( dTmp * ( dTmp - 1.0 ) * L2( pR ));
}

double
dL2_dBeta0_dPhi( RECURSE_PARAMS *pR )
{
    return( pR->dSumB * pR->dSumC * 
            exp( pR->dParams[P_BETA0] + 
                 pR->dParams[P_BETA1] *
                 pR->lCovariateX +
                 pR->dParams[P_PHI] ) *
            L2( pR ));
}

double
dL2_dBeta1_dBeta1( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumB * 
                  pR->lCovariateX *
                  exp( pR->dParams[P_BETA0] +
                       pR->dParams[P_BETA1] * 
                       pR->lCovariateX );
    return( dTmp * ( dTmp - pR->lCovariateX ) * L2( pR ));
}

double
dL2_dBeta1_dPhi( RECURSE_PARAMS *pR )
{
    return( pR->dSumB * pR->dSumC * 
            pR->lCovariateX *
            exp( pR->dParams[P_BETA0] + 
                 pR->dParams[P_BETA1] *
                 pR->lCovariateX +
                 pR->dParams[P_PHI] ) *
            L2( pR ));
}

double
dL2_dPhi_dPhi( RECURSE_PARAMS *pR )
{
    double dTmp = pR->dSumC * exp( pR->dParams[P_PHI] );
    return( dTmp * ( dTmp - 1.0 ) * L2( pR ));
}

void
Hessian2( double *pdParameters, double *padReturn )
{
    double adSum[10];
    BRANCH_SUM aBranches[10];
    double *padHessian, *padLogLik, *padScoreVector;
    long lSubject;
    RECURSE_PARAMS RParams;
    int iNumDerivatives;
    int i;
    
    if ( glNumSubjects == 0 )
        return;

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 10;
    padHessian = adSum;
    iNumDerivatives = 6;
    RParams.pfnEquation[0] = dL2_dBeta0_dBeta0;
    RParams.pfnEquation[1] = dL2_dBeta0_dBeta1;
    RParams.pfnEquation[2] = dL2_dBeta0_dPhi;
    RParams.pfnEquation[3] = dL2_dBeta1_dBeta1;
    RParams.pfnEquation[4] = dL2_dBeta1_dPhi;
    RParams.pfnEquation[5] = dL2_dPhi_dPhi;
    padScoreVector = adSum + iNumDerivatives;
    RParams.pfnEquation[6] = S2Beta0;
    RParams.pfnEquation[7] = S2Beta1;
    RParams.pfnEquation[8] = S2Phi;
    padLogLik = adSum + iNumDerivatives + 3;
    RParams.pfnEquation[9] = L2;

    for ( i = 0; i < iNumDerivatives; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;
        RParams.lCovariateX = gaSubjects[lSubject].lCovariateX[0];

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        /* dL2_dBeta0_dBeta0 */
        padReturn[0] += ( padHessian[0] * padLogLik[0] - 
                          padScoreVector[0] * padScoreVector[0] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL2_dBeta0_dBeta1 */
        padReturn[1] += ( padHessian[1] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL2_dBeta0_dPhi */
        padReturn[2] += ( padHessian[2] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL2_dBeta1_dBeta1 */
        padReturn[4] += ( padHessian[3] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL2_dBeta1_dPhi */
        padReturn[5] += ( padHessian[4] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL2_dPhi_dPhi */
        padReturn[8] += ( padHessian[5] * padLogLik[0] - 
                          padScoreVector[2] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
    }
    /* dL2_dBeta1_dBeta0 */
    padReturn[3] = padReturn[1];
    /* dL2_dPhi_dBeta0 */
    padReturn[6] = padReturn[2];
    /* dL2_dPhi_dBeta1 */
    padReturn[7] = padReturn[5];
}
