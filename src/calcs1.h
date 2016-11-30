/*
 * calcs1.h
 * Contains the definitions of the calculator functions required to perform
 * logistic regression with mixed likelihoods.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __CALCS1_H__
#define __CALCS1_H__

extern void Likelihood1( double *, double * );
extern void LogLikelihood1( double *, double *, int * );
double L1( RECURSE_PARAMS * );

extern void ScoreVector1( double *, double * );
double S1Beta( RECURSE_PARAMS * );
double S1Phi( RECURSE_PARAMS * );

extern void Hessian1( double *, double * );
double dL1_dBeta_dBeta( RECURSE_PARAMS * );
double dL1_dBeta_dPhi( RECURSE_PARAMS * );
double dL1_dPhi_dPhi( RECURSE_PARAMS * );

#endif
