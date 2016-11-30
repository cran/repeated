/*
 * calcs.c
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "gaps.h"
#include "calcs.h"

double 
dPow( double dX, double dY )
{
    if (dY == 0.0 ) {
        return( 1.0 );
    }
    if ( dX == 0.0 ) {
            return( 0.0 );
    }
    if ( dX > 0.0 ) {
        return( exp( dY * log( dX )));
    }
    return( -exp( dY * log( -dX )));
}

/*
 * dChoose
 * Calculate lNum choose lDenom
 */
double
dChoose( long lNum, long lDenom )
{
    long lStart;
    /* double dNum, dDenom; */
    double dResult;
    long i, j;

    if ( lNum < lDenom )
        return 0.0;
    /*
     * calculate lNum! / ( lDenom! * ( lNum - lDenom )! )
     */
    if ( lNum < 0.0l || lDenom < 0.0l )
        /* should never occur as we're checking for bad file format */
        return 0.0;

    if ( lDenom > lNum - lDenom ) {
        lStart = lDenom;
        lDenom = lNum - lDenom;
    } else {
        lStart = lNum - lDenom;
    }

	dResult = 1.0;
	i = lNum;
	j = lDenom;
	while ( i > lStart || j > 1 ) {
		if ( i > lStart ) {
			if ( j > 1 )
				dResult *= ( (double)i-- / j-- );
			else
				dResult *= i--;
		} else if ( j > 1 )
			dResult /= j--;
	}

	/*
    dNum = 1.0;
    dDenom = 1.0;
    for ( i = lNum; i > lStart; i-- )
        dNum *= i;
    for ( i = lDenom; i > 1; i-- )
        dDenom *= i;
	*/

    /* no need to check for div by 0 */
    /* return( dNum / dDenom ); */
	return( dResult );
}

double
FcnAsubL( long lSubject, long lH, long lI, long lJ, long lK, long lL )
{
    struct _GAP_DATA *pGap = &( gaSubjects[lSubject].palGaps[lL] );
    double dValue = dChoose( pGap->aalR[0][0], lJ ) *
                    dChoose( pGap->aalR[1][1], lK ) *
                    dChoose( pGap->aalR[1][0] + lK, lH ) *
                    dChoose( pGap->aalR[0][1] + pGap->aalR[1][0] + lJ + lK, lI ) *
                    /* pow( -1.0, (double)( lH +lI + lJ + lK ));*/
                    ( 1 - 2 * (( lH + lI + lJ + lK ) % 2 ));
    return( dValue );
}

double
FcnBsubL( long lxSubject, long lH, long lI, long lJ, long lK, long lL )
{
    struct _GAP_DATA *pGap = &( gaSubjects[lxSubject].palGaps[lL] );

    return( (double)( pGap->aalR[0][1] + lJ + lH ));
}

double
FcnCsubL( long lxSubject, long lH, long lI, long lJ, long lK, long lL )
{
    return( (double)lI );
}

void
CalcRecurse( RECURSE_PARAMS *pR, BRANCH_SUM *pBranches )
{
    struct _GAP_DATA *pGap;
    long lH, lI, lJ, lK;
    RECURSE_PARAMS RParams;
    int i;
/*
    if ( glNumSubjects == 0 ) {
        printf( "Please load a dataset first\n" );
        return;
    }
*/
    if ( pR->lL == gaSubjects[pR->lSubjectID].lNumGaps ) {
        for ( i = 0; i < pR->iNumEqns; i++ ) {
            AddNode( pBranches + i, (*( pR->pfnEquation[i] ))( pR ));
        }
    } else {
        pGap = &( gaSubjects[pR->lSubjectID].palGaps[pR->lL] );

        for ( lK = 0; lK <= pGap->aalR[1][1]; lK++ ) {
            for ( lJ = 0; lJ <= pGap->aalR[0][0]; lJ++ ) {
                for ( lH = 0; lH <= pGap->aalR[1][0] + lK; lH++ ) {
                    for ( lI = 0; 
                          lI <= pGap->aalR[0][1] + pGap->aalR[1][0] + lJ + lK;
                          lI++ ) {

                        memcpy( &RParams, pR, sizeof( RECURSE_PARAMS ));

                        RParams.dPrdA *= FcnAsubL( RParams.lSubjectID, lH, lI, lJ, lK, pR->lL );
                        RParams.dSumB += FcnBsubL( RParams.lSubjectID, lH, lI, lJ, lK, pR->lL );
                        RParams.dSumC += FcnCsubL( RParams.lSubjectID, lH, lI, lJ, lK, pR->lL ) 
                                         * pGap->lTime;

                        RParams.lL++;
                        CalcRecurse( &RParams, pBranches );
                    }
                }
            }
        }
    }
}

#if defined( UNSIGNED_METHOD )

int 
compare( const void *arg1, const void *arg2 )
{
    double dVal1 = fabs( *(double *)arg1 );
    double dVal2 = fabs( *(double *)arg2 );

    if ( dVal1 == dVal2 )
        return( 0 );

    if ( dVal1 > dVal2 )
        return( 1 );

    return( -1 );
}

#else

int 
compare( const void *arg1, const void *arg2 )
{
    double dVal1 = *(double *)arg1;
    double dVal2 = *(double *)arg2;

    if ( dVal1 == dVal2 )
        return( 0 );

    if ( dVal1 > dVal2 )
        return( 1 );

    return( -1 );
}

#endif

#if !defined( MORE_PRECISE )
 
void
AddNode( BRANCH_SUM *pBranch, double dValue )
{
    if ( pBranch->dPartialSum == 0.0 )
        pBranch->dPartialSum = dValue;
    else {
        if (( pBranch->dPartialSum < 0 && dValue > 0 ) ||
            ( pBranch->dPartialSum > 0 && dValue < 0 )) {
            pBranch->dSum += ( pBranch->dPartialSum + dValue );
            pBranch->dPartialSum = 0.0;
        } else
            pBranch->dPartialSum += dValue;
    }
}

double
SumNodes( BRANCH_SUM *pBranch )
{
    double dSum = pBranch->dSum + pBranch->dPartialSum;

    memset( pBranch, 0, sizeof( BRANCH_SUM ));

    return( dSum );
}

#else

void
AddNode( BRANCH_SUM *pBranch, double dValue )
{
    if ( pBranch->iAllocSize == 0 ) {
        pBranch->iAllocSize = SIZE_BRANCH_BLOCKS;
        pBranch->padValues = calloc( pBranch->iAllocSize, sizeof( double ));
    } else if ( pBranch->iNumElements >= pBranch->iAllocSize ) {
        pBranch->iAllocSize += SIZE_BRANCH_BLOCKS;
        pBranch->padValues = realloc( pBranch->padValues, 
                                     sizeof( double ) * 
                                     pBranch->iAllocSize );
    }
    pBranch->padValues[pBranch->iNumElements++] = dValue;
}

#if defined( UNSIGNED_METHOD )

double
SumNodes( BRANCH_SUM *pBranch )
{
    int i;
    double dPartialSum = 0.0, dTempSum = 0.0;
    double dReturn = 0.0;

    qsort( pBranch->padValues, pBranch->iNumElements, 
           sizeof( pBranch->padValues[0] ), compare );

    for ( i = 0; i < pBranch->iNumElements; i++ ) {
        if ( dPartialSum == 0.0 )
            dPartialSum = pBranch->padValues[i];
        else {
            dTempSum = dPartialSum + pBranch->padValues[i];
            if (( dPartialSum > 0 && dTempSum < dPartialSum ) ||
                ( dPartialSum < 0 && dTempSum > dPartialSum )) {
                dReturn += dTempSum;
                dPartialSum = 0.0;
            } else
                dPartialSum = dTempSum;
        }
    }
    dReturn += dPartialSum;

    free( pBranch->padValues );
    pBranch->padValues = NULL;
    pBranch->iNumElements = 0;
    pBranch->iAllocSize = 0;

    return( dReturn );
}

#else

double
SumNodes( BRANCH_SUM *pBranch )
{
    int i, j;
    double dReturn = 0.0;

    qsort( pBranch->padValues, pBranch->iNumElements, 
           sizeof( pBranch->padValues[0] ), compare );

    i = 0;
    j = pBranch->iNumElements - 1;
    while( i <= j ) {
        if ( i == j )
            dReturn += pBranch->padValues[i++];
        else
            dReturn += ( pBranch->padValues[i++] +
                         pBranch->padValues[j--] );
    }

    free( pBranch->padValues );
    pBranch->padValues = NULL;
    pBranch->iNumElements = 0;
    pBranch->iAllocSize = 0;

    return( dReturn );
}

#endif
#endif
