/*
 * defs.h
 * Contains the global definitions.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __DEFS_H__
#define __DEFS_H__

/*
 * Boolean definitions
 */
#define BOOL                    int
#define FALSE                   0
#define TRUE                    (!FALSE)
#define Inf                     1e308

/*
 * File defintions
 * MAX_DATALINE = maximum characters per line
 */
#define MAX_DATALINE            256
#ifndef _MAX_PATH
#define _MAX_PATH       512
#endif

/*
 * Default settings
 * TOLERANCE = floating point tolerance
 * MAX_ITER = maximum number of iterations to perform
 */
#define TOLERANCE               10E-10
#define MAX_ITER                100


#define NUM_IDENTIFIERS         3   /* Subject ID, Gap ID, Gap Time */
#define MAX_STATES              2   /* maximum number of states */
#define MAX_COVARIATES          1   /* maximum number of covariates */
#define MIN_NUM_COLUMNS         NUM_IDENTIFIERS + MAX_STATES * MAX_STATES
#define MAX_NUM_COLUMNS         MIN_NUM_COLUMNS + MAX_COVARIATES
#define COL_SUBJECT_ID          0
#define COL_GAP_ID              1
#define COL_GAP_TIME            2
#define COL_RESPONSES           3
#define COL_COVARIATE           MIN_NUM_COLUMNS


#define NUM_FACTORS             3
#define FACTOR_A                0
#define FACTOR_B                1
#define FACTOR_C                2


#define SIZE_SUBJECT_BLOCKS     100
#define SIZE_BRANCH_BLOCKS      100

#define ID_NO_ERROR             0

#define ID_ERR_MEMORY           1000

#define ID_ERR_NOFILENAME       2000
#define ID_ERR_NOFILE           2001
#define ID_ERR_INCORRECTFORMAT  2100

#endif
