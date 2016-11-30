/*
 * calcs3.c
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
#include "calcs3.h"

double 
f3( RECURSE_PARAMS *pR, int *paiFactors) 
{
    return( dPow( pR->dSumB, paiFactors[0] ) *
            exp( paiFactors[1] * pR->dParams[P_DELTA] +
                 paiFactors[2] * pR->dParams[P_BETA0] +
                 paiFactors[3] * pR->dParams[P_BETA1] * pR->lCovariateX +
                 paiFactors[4] * pR->dSumC * exp( pR->dParams[P_PHI] )));
}

double
L3( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   }   };
    return( pR->dPrdA * 
            dPow( 1 + f3( pR, aiFactors[0] ),
                 -exp( -pR->dParams[P_DELTA] )) *
            f3( pR, aiFactors[1] ));
}

void
LogLikelihood3( double *pdParameters, double *padReturn, int *piError )
{
    BRANCH_SUM aBranches[1];
    long lSubject;
    double dValue;
    RECURSE_PARAMS RParams;
    
    padReturn[0] = 0.0;

    *piError = ERR_NONE;

    if ( glNumSubjects == 0 )
        return;

    memset( aBranches, 0, sizeof( aBranches ));
    memcpy( RParams.dParams, pdParameters, sizeof( RParams.dParams ));
    RParams.dPrdA = 1.0;
    RParams.dSumB = 0.0;
    RParams.dSumC = 0.0;
    RParams.lL = 0;
    RParams.iNumEqns = 1;
    RParams.pfnEquation[0] = L3;

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
S3Beta0( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1, -1   }   };

    return( -dPow( 1 + f3( pR, aiFactors[0] ),
                  ( -exp( -pR->dParams[P_DELTA] ) - 1 )) *
            pR->dPrdA *
            f3( pR, aiFactors[1] ));
}

double
S3Beta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1, -1   }   };

    return( -dPow( 1 + f3( pR, aiFactors[0] ),
                  ( -exp( -pR->dParams[P_DELTA] ) - 1 )) *
            pR->dPrdA * pR->lCovariateX *
            f3( pR, aiFactors[1] ));
}

double
S3Phi( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   }   };

    return( -dPow( 1 + f3( pR, aiFactors[0] ),
                  ( -exp( -pR->dParams[P_DELTA] ))) *
            pR->dPrdA *
            pR->dSumC *
            exp( pR->dParams[P_PHI] ) *
            f3( pR, aiFactors[1] ));
}

double
S3Delta( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = f3( pR, aiFactors[1] );

    return( pR->dPrdA * 
            dPow( dEval1, ( -exp( -pR->dParams[P_DELTA] ) - 1 )) *
            ( exp( -pR->dParams[P_DELTA] ) *
              log( dEval1 ) + log( dEval1 ) * dEval2 - dEval2 ) *
            f3( pR, aiFactors[2] ));
}

void
ScoreVector3( double *pdParameters, double *padReturn, int *cov)
{
    double adSum[5];
    BRANCH_SUM aBranches[5];
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
    RParams.iNumEqns = *cov?5:4;
    padScoreVector = adSum;
    iNumDerivatives = *cov?4:3;
    RParams.pfnEquation[0] = S3Beta0;
    RParams.pfnEquation[1] = *cov?S3Beta1:S3Phi;
    RParams.pfnEquation[2] = *cov?S3Phi:S3Delta;
    RParams.pfnEquation[3] = *cov?S3Delta:L3;
    if(*cov)RParams.pfnEquation[4] = L3;
    padLogLik = adSum + iNumDerivatives;

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
dL3_dBeta0_dBeta0( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  2,  2, -1   },
                                {   0,  1,  2,  2, -1   }   };

    double dEval1 = dPow( 1 + f3( pR, aiFactors[0] ),
                         -exp( -pR->dParams[P_DELTA] ) - 2 ) *
                    pR->dPrdA * pR->dSumB * pR->dSumB;
    return( dEval1 * f3( pR, aiFactors[1] ) +
            S3Beta0( pR ) +
            dEval1 * f3( pR, aiFactors[2] ));
}

double
dL3_dBeta0_dBeta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  2,  2, -1   },
                                {   0,  1,  2,  2, -1   }   };
    double dEval1 = dPow( 1 + f3( pR, aiFactors[0] ),
                         -exp( -pR->dParams[P_DELTA] ) - 2 ) *
                    pR->dPrdA * pR->dSumB * pR->dSumB *
                    pR->lCovariateX;
    return( dEval1 * f3( pR, aiFactors[1] ) +
            S3Beta1( pR ) +
            dEval1 * f3( pR, aiFactors[2] ));
}

double
dL3_dBeta0_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1, -1   }   };
    return( dPow( 1 + f3( pR, aiFactors[0] ),
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dPrdA * pR->dSumC *
            f3( pR, aiFactors[1] ) *
            exp( pR->dParams[P_PHI] ));
}

double
dL3_dBeta0_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0, -1,  1,  1, -1   },
                                {   1,  0,  2,  2, -1   },
                                {   1,  1,  2,  2, -1   },
                                {   2,  2,  3,  3, -1   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, 
                         -exp( -pR->dParams[P_DELTA] ) - 1 );
    return( pR->dPrdA * pR->dSumB *
            ( -dEval2 * log( dEval1 ) *
              f3( pR, aiFactors[1] ) -
              dEval2 * log( dEval1 ) * f3( pR, aiFactors[2] ) +
              dEval2 * f3( pR, aiFactors[2] ) +
              dEval2 * f3( pR, aiFactors[3] ) +
              dEval2 * f3( pR, aiFactors[4] )) /
            dEval1 );
}

double
dL3_dBeta1_dBeta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   2,  0,  2,  2, -1   },
                                {   2,  1,  2,  2, -1   }   };
    double dEval1 = dPow( 1 + f3( pR, aiFactors[0] ), 
                         -exp( -pR->dParams[P_DELTA] ) - 2 );
    double dEval2 = dEval1 * pR->dPrdA * pR->lCovariateX * pR->lCovariateX;
    return( dEval2 * f3( pR, aiFactors[1] ) + 
            S3Beta1( pR ) +
            dEval2 * f3( pR, aiFactors[2] ));
}

double
dL3_dBeta1_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1, -1   }   };
    double dEval1 = dPow( 1 + f3( pR, aiFactors[0] ), 
                         -exp( -pR->dParams[P_DELTA] ) - 1 );
    return( dEval1 * pR->dPrdA * pR->lCovariateX * pR->dSumC *
            f3( pR, aiFactors[1] ) * exp( pR->dParams[P_PHI] ));
}

double
dL3_dBeta1_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0, -1,  1,  1, -1   },
                                {   1,  0,  2,  2, -1   },
                                {   1,  1,  2,  2, -1   },
                                {   2,  2,  3,  3, -1   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 );
    double dEval3 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 2 );
    return( pR->dPrdA * pR->dSumB * pR->lCovariateX * 
            ( -dEval2 * log( dEval1 ) * f3( pR, aiFactors[1] ) -
              dEval2 * log( dEval1 ) * f3( pR, aiFactors[2] ) +
              dEval2 * f3( pR, aiFactors[2] ) +
              dEval3 * f3( pR, aiFactors[3] ) +
              dEval3 * f3( pR, aiFactors[4] )) /
            dEval1 );
}

double
dL3_dPhi_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ));
    return( -dEval2 * pR->dPrdA * pR->dSumC * 
            f3( pR, aiFactors[1] ) * 
            exp( pR->dParams[P_PHI] ) +
            dEval2 * pR->dPrdA * pR->dSumC * pR->dSumC *
            f3( pR, aiFactors[1] ) *
            exp( 2 * pR->dParams[P_PHI] ));
}

double
dL3_dPhi_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   1,  0,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 );
    return( -dEval2 * pR->dPrdA *
            ( exp( -pR->dParams[P_DELTA] ) *
              log( dEval1 ) +
              log( dEval1 ) * f3( pR, aiFactors[1] ) -
              f3( pR, aiFactors[1] )) *
            pR->dSumC * exp( pR->dParams[P_PHI] ) *
            f3( pR, aiFactors[2] ));
}

double
dL3_dDelta_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][5] =    {   {   1,  1,  1,  1,  0   },
                                {   0,  0,  0,  0, -1   },
                                {   0, -2,  0,  0,  0   },
                                {   1, -1,  1,  1,  0   },
                                {   2,  0,  2,  2,  0   },
                                {   3,  1,  3,  3,  0   },
                                {   4,  2,  4,  4,  0   },
                                {   1,  0,  1,  1,  0   },
                                {   2,  1,  2,  2,  0   }   };
    double dEval1 = 1 + f3( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 2 );
    double dEval3 = log( dEval1 ) * log( dEval1 );
    double dEval4 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ));
    return( pR->dPrdA * f3( pR, aiFactors[1] ) *
            ( dEval2 * f3( pR, aiFactors[2] ) * dEval3 +
              4 * dEval2 * dEval3 * f3( pR, aiFactors[3] ) +
              6 * dEval2 * dEval3 * f3( pR, aiFactors[4] ) +
              4 * dEval2 * dEval3 * f3( pR, aiFactors[5] ) -
              2 * dEval2 * log( dEval1 ) * f3( pR, aiFactors[3] ) -
              6 * dEval2 * log( dEval1 ) * f3( pR, aiFactors[4] ) -
              6 * dEval2 * log( dEval1 ) * f3( pR, aiFactors[5] ) +
              dEval2 * dEval3 * f3( pR, aiFactors[6] ) -
              2 * dEval2 * log( dEval1 ) * f3( pR, aiFactors[6] ) +
              dEval2 * f3( pR, aiFactors[4] ) +
              2 * dEval2 * f3( pR, aiFactors[5] ) +
              dEval2 * f3( pR, aiFactors[6] ) -
              dEval4 * exp( -pR->dParams[P_DELTA] ) * log( dEval1 ) -
              2 * dEval4 * log( dEval1 ) * f3( pR, aiFactors[7] ) -
              dEval4 * log( dEval1 ) * f3( pR, aiFactors[8] ) +
              dEval4 * f3( pR, aiFactors[7] ) +
              2 * dEval4 * f3( pR, aiFactors[8] )) /
            dPow( dEval1, 2 ));
}

void
Hessian3( double *pdParameters, double *padReturn, int *cov)
{
    double adSum[15];
    BRANCH_SUM aBranches[15];
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
    RParams.iNumEqns = *cov?15:10;
    padHessian = adSum;
    iNumDerivatives = *cov?10:6;
    RParams.pfnEquation[0] = dL3_dBeta0_dBeta0;
    RParams.pfnEquation[1] = *cov?dL3_dBeta0_dBeta1:dL3_dBeta0_dPhi;
    RParams.pfnEquation[2] = *cov?dL3_dBeta0_dPhi:dL3_dBeta0_dDelta;
    RParams.pfnEquation[3] = *cov?dL3_dBeta0_dDelta:dL3_dPhi_dPhi;
    RParams.pfnEquation[4] = *cov?dL3_dBeta1_dBeta1:dL3_dPhi_dDelta;
    RParams.pfnEquation[5] = *cov?dL3_dBeta1_dPhi:dL3_dDelta_dDelta;
    RParams.pfnEquation[6] = *cov?dL3_dBeta1_dDelta:S3Beta0;
    RParams.pfnEquation[7] = *cov?dL3_dPhi_dPhi:S3Phi;
    RParams.pfnEquation[8] = *cov?dL3_dPhi_dDelta:S3Delta;
    RParams.pfnEquation[9] = *cov?dL3_dDelta_dDelta:L3;
    if(*cov){
      RParams.pfnEquation[10] = S3Beta0;
      RParams.pfnEquation[11] = S3Beta1;
      RParams.pfnEquation[12] = S3Phi;
      RParams.pfnEquation[13] = S3Delta;
      RParams.pfnEquation[14] = L3;}
    padScoreVector = adSum + iNumDerivatives;
    padLogLik = adSum + iNumDerivatives + 3 + *cov;

    for ( i = 0; i < iNumDerivatives; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;
        RParams.lCovariateX = gaSubjects[lSubject].lCovariateX[0];

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        /* dL3_dBeta0_dBeta0 */
        padReturn[0] += ( padHessian[0] * padLogLik[0] - 
                          padScoreVector[0] * padScoreVector[0] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dBeta0_dBeta1 */
        padReturn[1] += ( padHessian[1] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dBeta0_dPhi */
        padReturn[2] += ( padHessian[2] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dBeta0_dDelta */
        padReturn[*cov?3:4] += ( padHessian[3] * padLogLik[0] -
                          padScoreVector[*cov?0:1] * padScoreVector[*cov?3:1])/
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dBeta1_dBeta1 */
        padReturn[5] += ( padHessian[4] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[*cov?1:2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dBeta1_dPhi */
        padReturn[*cov?6:8] += ( padHessian[5] * padLogLik[0] - 
                          padScoreVector[*cov?1:2] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
	if(*cov){
        /* dL3_dBeta1_dDelta */
        padReturn[7] += ( padHessian[6] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[3] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL3_dPhi_dPhi */
        padReturn[10] += ( padHessian[7] * padLogLik[0] - 
                           padScoreVector[2] * padScoreVector[2] ) /
                         ( padLogLik[0] * padLogLik[0] );
        /* dL3_dPhi_dDelta */
        padReturn[11] += ( padHessian[8] * padLogLik[0] - 
                           padScoreVector[2] * padScoreVector[3] ) /
                         ( padLogLik[0] * padLogLik[0] );
        /* dL3_dDelta_dDelta */
        padReturn[15] += ( padHessian[9] * padLogLik[0] - 
                           padScoreVector[3] * padScoreVector[3] ) /
                         ( padLogLik[0] * padLogLik[0] );
	}}
    /* dL3_dBeta1_dBeta0 */
    padReturn[*cov?4:3] = padReturn[1];
    /* dL3_dPhi_dBeta0 */
    padReturn[*cov?8:6] = padReturn[2];
    /* dL3_dPhi_dBeta1 */
    padReturn[*cov?9:7] = padReturn[*cov?6:5];
    if(*cov){
      /* dL3_dDelta_dBeta0 */
      padReturn[12] = padReturn[3];
      /* dL3_dDelta_dBeta1 */
      padReturn[13] = padReturn[7];
      /* dL3_dDelta_dPhi */
      padReturn[14] = padReturn[11];
    }}
