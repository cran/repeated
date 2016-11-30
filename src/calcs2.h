/*
 * calcs2.h
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __CALCS2_H__
#define __CALCS2_H__

extern void Likelihood2( double *, double * );
extern void LogLikelihood2( double *, double *, int * );
double L2( RECURSE_PARAMS * );

extern void ScoreVector2( double *, double * );
double S2Beta0( RECURSE_PARAMS * );
double S2Beta1( RECURSE_PARAMS * );
double S2Phi( RECURSE_PARAMS * );

extern void Hessian2( double *, double * );
double dL2_dBeta0_dBeta0( RECURSE_PARAMS * );
double dL2_dBeta0_dBeta1( RECURSE_PARAMS * );
double dL2_dBeta0_dPhi( RECURSE_PARAMS * );
double dL2_dBeta1_dBeta1( RECURSE_PARAMS * );
double dL2_dBeta1_dPhi( RECURSE_PARAMS * );
double dL2_dPhi_dPhi( RECURSE_PARAMS * );

#endif
