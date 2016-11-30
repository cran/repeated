/*
 * calcs4.h
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __CALCS4_H__
#define __CALCS4_H__

extern void Likelihood4( double *, double * );
extern void LogLikelihood4( double *, double *, int * );
double L4( RECURSE_PARAMS * );

extern void ScoreVector4( double *, double *, int * );
double S4Beta0( RECURSE_PARAMS * );
double S4Beta1( RECURSE_PARAMS * );
double S4Phi( RECURSE_PARAMS * );
double S4Delta( RECURSE_PARAMS * );
double S4Theta( RECURSE_PARAMS * );

extern void Hessian4( double *, double *, int * );
double dL4_dBeta0_dBeta0( RECURSE_PARAMS * );
double dL4_dBeta0_dBeta1( RECURSE_PARAMS * );
double dL4_dBeta0_dPhi( RECURSE_PARAMS * );
double dL4_dBeta0_dDelta( RECURSE_PARAMS * );
double dL4_dBeta0_dTheta( RECURSE_PARAMS * );
double dL4_dBeta1_dBeta1( RECURSE_PARAMS * );
double dL4_dBeta1_dPhi( RECURSE_PARAMS * );
double dL4_dBeta1_dDelta( RECURSE_PARAMS * );
double dL4_dBeta1_dTheta( RECURSE_PARAMS * );
double dL4_dPhi_dPhi( RECURSE_PARAMS * );
double dL4_dPhi_dDelta( RECURSE_PARAMS * );
double dL4_dPhi_dTheta( RECURSE_PARAMS * );
double dL4_dDelta_dDelta( RECURSE_PARAMS * );
double dL4_dDelta_dTheta( RECURSE_PARAMS * );
double dL4_dTheta_dTheta( RECURSE_PARAMS * );

double f4( RECURSE_PARAMS *, int * );


#endif
