/*
 * calcs4.c
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
#include "calcs4.h"

double 
f4( RECURSE_PARAMS *pR, int *paiFactors ) 
{
    return( dPow( pR->dSumB, paiFactors[0] ) *
            dPow( pR->dSumC, paiFactors[1] ) *
            exp( paiFactors[2] * pR->dParams[P_DELTA] +
                 paiFactors[3] * pR->dParams[P_BETA0] +
                 paiFactors[4] * pR->dParams[P_BETA1] * pR->lCovariateX +
                 paiFactors[5] * pR->dParams[P_PHI] +
                 paiFactors[6] * pR->dParams[P_THETA] ));
}

double
L4( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   }   };
    return( pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] )) *
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] )));
}

void
LogLikelihood4( double *pdParameters, double *padReturn, int *piError )
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
    RParams.pfnEquation[0] = L4;

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
S4Beta0( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  0,  0,  1,  1,  0,  0   }   };
    return( -pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB * 
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] )) *
            f4( pR, aiFactors[2] ));
}

double
S4Beta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  0,  0,  1,  1,  0,  0   }   };

    return( -pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB * 
            pR->lCovariateX *
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] )) *
            f4( pR, aiFactors[2] ));
}

double
S4Phi( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   }   };

    return( -pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] )) *
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] ) - 1 ) *
            f4( pR, aiFactors[2] ));
}

double
S4Delta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   1,  0,  0,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   }   };
    double dEval1 = 1 + f4( pR, aiFactors[0] );

    return( pR->dPrdA * 
            dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            ( exp( -pR->dParams[P_DELTA] ) * log( dEval1 ) +
              log( dEval1 ) * f4( pR, aiFactors[1] ) -
              f4( pR, aiFactors[1] )) *
            dPow( 1 + f4( pR, aiFactors[2] ),
                 -exp( -pR->dParams[P_THETA] )));
}

double
S4Theta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   }   };
    double dEval1 = 1 + f4( pR, aiFactors[1] );
    double dEval2 = f4( pR, aiFactors[2] );

    return( pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] )) *
            dPow( dEval1, -exp( -pR->dParams[P_THETA] ) - 1 ) *
            ( exp( -pR->dParams[P_THETA] ) * log( dEval1 ) +
              log( dEval1 ) * dEval2 - dEval2 ));
}

void
ScoreVector4( double *pdParameters, double *padReturn, int *cov)
{
    double adSum[6];
    BRANCH_SUM aBranches[6];
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
    RParams.iNumEqns = *cov?6:5;
    padScoreVector = adSum;
    iNumDerivatives = *cov?5:4;
    RParams.pfnEquation[0] = S4Beta0;
    RParams.pfnEquation[1] = *cov?S4Beta1:S4Phi;
    RParams.pfnEquation[2] = *cov?S4Phi:S4Delta;
    RParams.pfnEquation[3] = *cov?S4Delta:S4Theta;
    RParams.pfnEquation[4] = *cov?S4Theta:L4;
    if(*cov)RParams.pfnEquation[5] = L4;
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
dL4_dBeta0_dBeta0( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   2,  0,  0,  2,  2,  0,  0   },
                                {   2,  0,  1,  2,  2,  0,  0   }   };
    double dEval1 = pR->dPrdA * dPow( 1 + f4( pR, aiFactors[0] ),
                                     -exp( -pR->dParams[P_DELTA] ) - 2 );
    double dEval2 = dPow( 1 + f4( pR, aiFactors[1] ), 
                         -exp( -pR->dParams[P_THETA] ));
        
    return( dEval1 * dEval2 * f4( pR, aiFactors[2] ) +
            S4Beta0( pR ) +
            dEval1 * dEval2 * f4( pR, aiFactors[3] ));
}

double
dL4_dBeta0_dBeta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   2,  0,  0,  2,  2,  0,  0   },
                                {   2,  0,  1,  2,  2,  0,  0   }   };
    double dEval1 = pR->dPrdA * dPow( 1 + f4( pR, aiFactors[0] ),
                                     -exp( -pR->dParams[P_DELTA] ) - 2 );
    double dEval2 = dPow( 1 + f4( pR, aiFactors[1] ), 
                         -exp( -pR->dParams[P_THETA] ));
        
    return( dEval1 * dEval2 * pR->lCovariateX * f4( pR, aiFactors[2] ) +
            S4Beta1( pR ) +
            dEval1 * dEval2 * pR->lCovariateX * f4( pR, aiFactors[3] ));
}

double
dL4_dBeta0_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  1,  1,  1,  0   }   };
    return( pR->dPrdA * dPow( 1 + f4( pR, aiFactors[0] ),
                             -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB * dPow( 1 + f4( pR, aiFactors[1] ),
                             -exp( -pR->dParams[P_THETA] ) - 1 ) *
            f4( pR, aiFactors[2] ));
}

double
dL4_dBeta0_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  0, -1,  1,  1,  0,  0   },
                                {   1,  0,  0,  2,  2,  0,  0   },
                                {   1,  0,  1,  2,  2,  0,  0   },
                                {   2,  0,  2,  3,  3,  0,  0   }   };
    double dEval1 = 1 + f4( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 );
    double dEval3 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 2 );

    return( pR->dPrdA * pR->dSumB * 
            dPow( 1 + f4( pR, aiFactors[1] ), -exp( -pR->dParams[P_THETA] )) *
            ( -dEval2 * log( dEval1 ) * f4( pR, aiFactors[2] ) -
              dEval2 * log( dEval1 ) * f4( pR, aiFactors[3] ) +
              dEval2 * f4( pR, aiFactors[3] ) +
              dEval3 * f4( pR, aiFactors[4] ) +
              dEval3 * f4( pR, aiFactors[5] )) /
            dEval1 );
}

double
dL4_dBeta0_dTheta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   },
                                {   0,  0,  0,  1,  1,  0,  0   }   };

    return( -pR->dPrdA *
            dPow( 1 + f4( pR, aiFactors[0] ),
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB *
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] ) - 1 ) *
            ( exp( -pR->dParams[P_THETA] ) *
              log( 1 + f4( pR, aiFactors[1] )) +
              log( 1 + f4( pR, aiFactors[1] )) *
              f4( pR, aiFactors[2] ) -
              f4( pR, aiFactors[2] )) *
            f4( pR, aiFactors[3] ));
}

double
dL4_dBeta1_dBeta1( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   2,  0,  0,  2,  2,  0,  0   },
                                {   1,  0,  0,  1,  1,  0,  0   },
                                {   2,  0,  1,  2,  2,  0,  0   }   };
    double dEval1 = dPow( 1 + f4( pR, aiFactors[1] ),
                         -exp( -pR->dParams[P_THETA] ));
    double dEval2 = 1 + f4( pR, aiFactors[0] );
    double dEval3 = dPow( dEval2, -exp( -pR->dParams[P_DELTA] ) - 2 );
    double dEval4 = dPow( dEval2, -exp( -pR->dParams[P_DELTA] ) - 1 );
    
    return( pR->dPrdA * dEval3 * dEval1 * 
            pR->lCovariateX * pR->lCovariateX *
            f4( pR, aiFactors[2] ) -
            pR->dPrdA * dEval4 * dEval1 *
            pR->lCovariateX * pR->lCovariateX *
            f4( pR, aiFactors[3] ) +
            pR->dPrdA * dEval3 * dEval1 *
            pR->lCovariateX * pR->lCovariateX *
            f4( pR, aiFactors[4] ));
}

double
dL4_dBeta1_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  1,  1,  1,  0   }   };

    return( pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ), 
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB * pR->lCovariateX *
            dPow( 1 + f4( pR, aiFactors[1] ),
                 -exp( -pR->dParams[P_THETA] ) - 1 ) *
            f4( pR, aiFactors[2] ));
}

double
dL4_dBeta1_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  0, -1,  1,  1,  0,  0   },
                                {   1,  0,  0,  2,  2,  0,  0   },
                                {   1,  0,  1,  2,  2,  0,  0   },
                                {   2,  0,  2,  3,  3,  0,  0   }   };
    double dEval1 = 1 + f4( pR, aiFactors[0] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 );
    double dEval3 = dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 2 );

    return( pR->dPrdA * pR->dSumB * pR->lCovariateX *
            dPow( 1 + f4( pR, aiFactors[1] ), -exp( -pR->dParams[P_THETA] )) *
            ( -dEval2 * log( dEval1 ) * f4( pR, aiFactors[2] ) -
              dEval2 * log( dEval1 ) * f4( pR, aiFactors[3] ) +
              dEval2 * f4( pR, aiFactors[3] ) +
              dEval3 * f4( pR, aiFactors[4] ) +
              dEval3 * f4( pR, aiFactors[5] ))  /
            dEval1 );
}

double
dL4_dBeta1_dTheta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   },
                                {   0,  0,  0,  1,  1,  0,  0   }   };
    double dEval1 = 1 + f4( pR, aiFactors[1] );

    return( -pR->dPrdA * 
            dPow( 1 + f4( pR, aiFactors[0] ),
                 -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            pR->dSumB * pR->lCovariateX *
            dPow( dEval1, -exp( -pR->dParams[P_THETA] ) - 1 ) * 
            ( exp( -pR->dParams[P_THETA] ) * log( dEval1 ) +
              log( dEval1 ) * f4( pR, aiFactors[2] ) -
              f4( pR, aiFactors[2] )) *
            f4( pR, aiFactors[3] ));
}

double
dL4_dPhi_dPhi( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  2,  0,  0,  0,  2,  0   },
                                {   0,  1,  0,  0,  0,  1,  0   },
                                {   0,  2,  0,  0,  0,  2,  1   }   };
    double dEval1 = dPow( 1 + f4( pR, aiFactors[0] ), 
                         -exp( -pR->dParams[P_DELTA] ));
    double dEval2 = 1 + f4( pR, aiFactors[1] );
    double dEval3 = dPow( dEval2, -exp( -pR->dParams[P_THETA] ) - 1 );
    double dEval4 = dPow( dEval2, -exp( -pR->dParams[P_THETA] ) - 2 );

    return( pR->dPrdA * dEval1 * dEval4 * f4( pR, aiFactors[2] ) -
            pR->dPrdA * dEval1 * dEval3 * f4( pR, aiFactors[3] ) +
            pR->dPrdA * dEval1 * dEval4 * f4( pR, aiFactors[4] ));
}

double
dL4_dPhi_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   1,  0,  0,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   }   };
    double dEval1 = 1 + f4( pR, aiFactors[0] );

    return( -pR->dPrdA * 
            dPow( dEval1, -exp( -pR->dParams[P_DELTA] ) - 1 ) *
            ( exp( -pR->dParams[P_DELTA] ) * log( dEval1 ) +
              log( dEval1 ) * f4( pR, aiFactors[1] ) -
              f4( pR, aiFactors[1] )) *
            dPow( 1 + f4( pR, aiFactors[2] ), -exp( -pR->dParams[P_THETA] ) - 1 ) *
            f4( pR, aiFactors[3] ));
}

double
dL4_dPhi_dTheta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  0,  0,  0,  0,  1, -1   },
                                {   0,  1,  0,  0,  0,  2,  0   },
                                {   0,  1,  0,  0,  0,  2,  1   },
                                {   0,  2,  0,  0,  0,  3,  2   }   };
    double dEval1 = 1 + f4( pR, aiFactors[1] );
    double dEval2 = dPow( dEval1, -exp( -pR->dParams[P_THETA] ) - 1 );
    double dEval3 = dPow( dEval1, -exp( -pR->dParams[P_THETA] ) - 2 );

    return( pR->dPrdA *
            dPow( 1 + f4( pR, aiFactors[0] ), -exp( -pR->dParams[P_DELTA] )) *
            pR->dSumC *
            ( -dEval2 * log( dEval1 ) * f4( pR, aiFactors[2] ) -
              dEval2 * log( dEval1 ) * f4( pR, aiFactors[3] ) +
              dEval2 * f4( pR, aiFactors[3] ) +
              dEval3 * f4( pR, aiFactors[4] ) +
              dEval3 * f4( pR, aiFactors[5] )) /
            dEval1 );
}

double
dL4_dDelta_dDelta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   1,  0, -1,  1,  1,  0,  0   },
                                {   2,  0,  0,  2,  2,  0,  0   },
                                {   3,  0,  1,  3,  3,  0,  0   },
                                {   4,  0,  2,  4,  4,  0,  0   },
                                {   1,  0,  0,  1,  1,  0,  0   },
                                {   2,  0,  1,  2,  2,  0,  0   }   };
    double dCall[8];
    double dEval0;
    double dEval1;
    double dEval2;
    double dEval3;
    double dEval4;
    int i;

    for( i = 0; i < sizeof( dCall ) / sizeof( dCall[0] ); i++ )
        dCall[i] = f4( pR, aiFactors[i] );
    dEval0 = exp( -pR->dParams[P_DELTA] );
    dEval1 = 1 + dCall[0];
    dEval2 = dPow( dEval1, -dEval0 - 2 );
    dEval3 = dPow( dEval1, -dEval0 );
    dEval4 = log( dEval1 );

    return( pR->dPrdA * 
            dPow( 1 + dCall[1], -exp( -pR->dParams[P_THETA] )) *
            ( dEval2 * exp( -2 * pR->dParams[P_DELTA] ) * dPow( dEval4, 2 ) +
              4 * dEval2 * dPow( dEval4, 2 ) * dCall[2] +
              6 * dEval2 * dPow( dEval4, 2 ) * dCall[3] +
              4 * dEval2 * dPow( dEval4, 2 ) * dCall[4] -
              2 * dEval2 * dEval4 * dCall[2] -
              6 * dEval2 * dEval4 * dCall[3] -
              6 * dEval2 * dEval4 * dCall[4] +
              dEval2 * dPow( dEval4, 2 ) * dCall[5] -
              2 * dEval2 * dEval4 * dCall[5] +
              dEval2 * dCall[3] +
              2 * dEval2 * dCall[4] +
              dEval2 * dCall[5] -
              dEval3 * dEval0 * dEval4 -
              2 * dEval3 * dEval4 * dCall[6] -
              dEval3 * dEval4 * dCall[7] +
              dEval3 * dCall[6] +
              2 * dEval3 * dCall[7] ) /
            dPow( dEval1, 2 ));
}

double
dL4_dDelta_dTheta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   1,  0,  1,  1,  1,  0,  0   },
                                {   1,  0,  0,  1,  1,  0,  0   },
                                {   0,  1,  0,  0,  0,  1,  1   },
                                {   0,  1,  0,  0,  0,  1,  0   }   };
    double dCall[4];
    double dEval0;
    double dEval1;
    int i;

    for( i = 0; i < sizeof( dCall ) / sizeof( dCall[0] ); i++ )
        dCall[i] = f4( pR, aiFactors[i] );
    dEval0 = exp( -pR->dParams[P_DELTA] );
    dEval1 = exp( -pR->dParams[P_THETA] );

    return( pR->dPrdA * 
            dPow( 1 + dCall[0], -dEval0 - 1 ) *
            ( dEval0 * log( 1 + dCall[0] ) +
              log( 1 + dCall[0] ) * dCall[1] -
              dCall[1] ) *
            dPow( 1 + dCall[2], -dEval1 - 1 ) *
            ( dEval1 * log( 1 + dCall[2] ) +
              log( 1 + dCall[2] ) * dCall[3] -
              dCall[3] ));
}

double
dL4_dTheta_dTheta( RECURSE_PARAMS *pR )
{
    int aiFactors[][7] =    {   {   0,  1,  0,  0,  0,  1,  1   },
                                {   1,  0,  1,  1,  1,  0,  0   },
                                {   0,  0,  0,  0,  0,  0, -2   },
                                {   0,  1,  0,  0,  0,  1, -1   },
                                {   0,  2,  0,  0,  0,  2,  0   },
                                {   0,  3,  0,  0,  0,  3,  1   },
                                {   0,  4,  0,  0,  0,  4,  2   },
                                {   0,  1,  0,  0,  0,  1,  0   },
                                {   0,  2,  0,  0,  0,  2,  1   }   };
    double dCall[9];
    double dEval0;
    double dEval1;
    double dEval2;
    double dEval3;
    double dEval4;
    int i;

    for( i = 0; i < sizeof( dCall ) / sizeof( dCall[0] ); i++ )
        dCall[i] = f4( pR, aiFactors[i] );
    dEval0 = exp( -pR->dParams[P_THETA] );
    dEval1 = 1 + dCall[0];
    dEval2 = dPow( dEval1, -dEval0 - 2 );
    dEval3 = dPow( dEval1, -dEval0 );
    dEval4 = log( dEval1 );

    return( dPow( 1 + dCall[1], -exp( -pR->dParams[P_DELTA] )) * pR->dPrdA *
            ( dEval2 * dCall[2] * dPow( dEval4, 2 ) +
              4 * dEval2 * dPow( dEval4, 2 ) * dCall[3] +
              6 * dEval2 * dPow( dEval4, 2 ) * dCall[4] +
              4 * dEval2 * dPow( dEval4, 2 ) * dCall[5] -
              2 * dEval2 * dEval4 * dCall[3] -
              6 * dEval2 * dEval4 * dCall[4] -
              6 * dEval2 * dEval4 * dCall[5] +
              dEval2 * dPow( dEval4, 2 ) * dCall[6] -
              2 * dEval2 * dEval4 * dCall[6] +
              dEval2 * dCall[4] +
              2 * dEval2 * dCall[5] +
              dEval2 * dCall[6] -
              dEval3 * dEval0 * dEval4 -
              2 * dEval3 * dEval4 * dCall[7] -
              dEval3 * dEval4 * dCall[8] +
              dEval3 * dCall[7] +
              2 * dEval3 * dCall[8] ) /
            dPow( dEval1, 2 ));
}

void
Hessian4( double *pdParameters, double *padReturn, int *cov)
{
    double adSum[21];
    BRANCH_SUM aBranches[21];
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
    RParams.iNumEqns = *cov?21:15;
    padHessian = adSum;
    iNumDerivatives = *cov?15:10;
    RParams.pfnEquation[0] = dL4_dBeta0_dBeta0;
    RParams.pfnEquation[1] = *cov?dL4_dBeta0_dBeta1:dL4_dBeta0_dPhi;
    RParams.pfnEquation[2] = *cov?dL4_dBeta0_dPhi:dL4_dBeta0_dDelta;
    RParams.pfnEquation[3] = *cov?dL4_dBeta0_dDelta:dL4_dBeta0_dTheta;
    RParams.pfnEquation[4] = *cov?dL4_dBeta0_dTheta:dL4_dPhi_dPhi;
    RParams.pfnEquation[5] = *cov?dL4_dBeta1_dBeta1:dL4_dPhi_dDelta;
    RParams.pfnEquation[6] = *cov?dL4_dBeta1_dPhi:dL4_dPhi_dTheta;
    RParams.pfnEquation[7] = *cov?dL4_dBeta1_dDelta:dL4_dDelta_dDelta;
    RParams.pfnEquation[8] = *cov?dL4_dBeta1_dTheta:dL4_dDelta_dTheta;
    RParams.pfnEquation[9] = *cov?dL4_dPhi_dPhi:dL4_dTheta_dTheta;
    RParams.pfnEquation[10] = *cov?dL4_dPhi_dDelta:S4Beta0;
    RParams.pfnEquation[11] = *cov?dL4_dPhi_dTheta:S4Phi;
    RParams.pfnEquation[12] = *cov?dL4_dDelta_dDelta:S4Delta;
    RParams.pfnEquation[13] = *cov?dL4_dDelta_dTheta:S4Theta;
    RParams.pfnEquation[14] = *cov?dL4_dTheta_dTheta:L4;
    RParams.pfnEquation[15] = S4Beta0;
    RParams.pfnEquation[16] = S4Beta1;
    RParams.pfnEquation[17] = S4Phi;
    RParams.pfnEquation[18] = S4Delta;
    RParams.pfnEquation[19] = S4Theta;
    RParams.pfnEquation[20] = L4;
    padScoreVector = adSum + iNumDerivatives;
    padLogLik = adSum + iNumDerivatives + 4 + *cov;

    for ( i = 0; i < iNumDerivatives; i++ )
        padReturn[i] = 0.0;

    for ( lSubject = 0; lSubject < glNumSubjects; lSubject++ ) {
        RParams.lSubjectID = lSubject;
        RParams.lCovariateX = gaSubjects[lSubject].lCovariateX[0];

        CalcRecurse( &RParams, aBranches );

        for ( i = 0; i < RParams.iNumEqns; i++ )
            adSum[i] = SumNodes( aBranches + i );

        /* dL4_dBeta0_dBeta0 */
        padReturn[0] += ( padHessian[0] * padLogLik[0] - 
                          padScoreVector[0] * padScoreVector[0] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta0_dBeta1 */
        padReturn[1] += ( padHessian[1] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[1] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta0_dPhi */
        padReturn[2] += ( padHessian[2] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta0_dDelta */
        padReturn[3] += ( padHessian[3] * padLogLik[0] -
                          padScoreVector[0] * padScoreVector[3] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta0_dTheta */
        padReturn[*cov?4:5] += ( padHessian[4] * padLogLik[0] -
                          padScoreVector[*cov?0:1] * padScoreVector[*cov?4:1])/
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta1_dBeta1 */
        padReturn[6] += ( padHessian[5] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[*cov?1:2] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta1_dPhi */
        padReturn[7] += ( padHessian[6] * padLogLik[0] - 
                          padScoreVector[1] * padScoreVector[*cov?2:3] ) /
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta1_dDelta */
        padReturn[*cov?8:10] += ( padHessian[7] * padLogLik[0] - 
                          padScoreVector[*cov?1:2] * padScoreVector[*cov?3:2])/
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dBeta1_dTheta */
        padReturn[*cov?9:11] += ( padHessian[8] * padLogLik[0] - 
                          padScoreVector[*cov?1:2] * padScoreVector[*cov?4:3])/
                        ( padLogLik[0] * padLogLik[0] );
        /* dL4_dPhi_dPhi */
        padReturn[*cov?12:15] += ( padHessian[9] * padLogLik[0] - 
                          padScoreVector[*cov?2:3] * padScoreVector[*cov?2:3])/
                         ( padLogLik[0] * padLogLik[0] );
	if(*cov){
	  /* dL4_dPhi_dDelta */
	  padReturn[13] += ( padHessian[10] * padLogLik[0] - 
                           padScoreVector[2] * padScoreVector[3] ) /
                         ( padLogLik[0] * padLogLik[0] );
	  /* dL4_dPhi_dTheta */
	  padReturn[14] += ( padHessian[11] * padLogLik[0] - 
                           padScoreVector[2] * padScoreVector[4] ) /
                         ( padLogLik[0] * padLogLik[0] );
	  /* dL4_dDelta_dDelta */
	  padReturn[18] += ( padHessian[12] * padLogLik[0] - 
                           padScoreVector[3] * padScoreVector[3] ) /
                         ( padLogLik[0] * padLogLik[0] );
	  /* dL4_dDelta_dTheta */
	  padReturn[19] += ( padHessian[13] * padLogLik[0] - 
                           padScoreVector[3] * padScoreVector[4] ) /
                         ( padLogLik[0] * padLogLik[0] );
	  /* dL4_dTheta_dTheta */
	  padReturn[24] += ( padHessian[14] * padLogLik[0] - 
                           padScoreVector[4] * padScoreVector[4] ) /
                         ( padLogLik[0] * padLogLik[0] );
	}}
    /* dL4_dBeta1_dBeta0 */
    padReturn[*cov?5:4] = padReturn[1];
    /* dL4_dPhi_dBeta0 */
    padReturn[*cov?1:8] = padReturn[2];
    /* dL4_dPhi_dBeta1 */
    padReturn[*cov?11:9] = padReturn[*cov?7:6];
    /* dL4_dDelta_dBeta0 */
    padReturn[*cov?15:12] = padReturn[3];
    /* dL4_dDelta_dBeta1 */
    padReturn[*cov?16:13] = padReturn[*cov?8:7];
    /* dL4_dDelta_dPhi */
    padReturn[*cov?17:14] = padReturn[*cov?13:11];
    if(*cov){
      /* dL4_dTheta_dBeta0 */
      padReturn[20] = padReturn[4];
      /* dL4_dTheta_dBeta1 */
      padReturn[21] = padReturn[9];
      /* dL4_dTheta_dPhi */
      padReturn[22] = padReturn[14];
      /* dL4_dTheta_dDelta */
      padReturn[23] = padReturn[19];
    }}
