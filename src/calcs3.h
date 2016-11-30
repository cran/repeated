/*
 * calcs3.h
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __CALCS3_H__
#define __CALCS3_H__

extern void Likelihood3( double *, double * );
extern void LogLikelihood3( double *, double *, int * );
double L3( RECURSE_PARAMS * );

extern void ScoreVector3( double *, double *, int * );
double S3Beta0( RECURSE_PARAMS * );
double S3Beta1( RECURSE_PARAMS * );
double S3Phi( RECURSE_PARAMS * );
double S3Delta( RECURSE_PARAMS * );

extern void Hessian3( double *, double *, int * );
double dL3_dBeta0_dBeta0( RECURSE_PARAMS * );
double dL3_dBeta0_dBeta1( RECURSE_PARAMS * );
double dL3_dBeta0_dPhi( RECURSE_PARAMS * );
double dL3_dBeta0_dDelta( RECURSE_PARAMS * );
double dL3_dBeta1_dBeta1( RECURSE_PARAMS * );
double dL3_dBeta1_dPhi( RECURSE_PARAMS * );
double dL3_dBeta1_dDelta( RECURSE_PARAMS * );
double dL3_dPhi_dPhi( RECURSE_PARAMS * );
double dL3_dPhi_dDelta( RECURSE_PARAMS * );
double dL3_dDelta_dDelta( RECURSE_PARAMS * );

double f3( RECURSE_PARAMS *, int * );

#endif
