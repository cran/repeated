/*
 * gaps.h
 * Defines the main application.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#ifndef __GAPS_H__
#define __GAPS_H__


/*
 * Structure definitions
 * SUBJECT_DATA = data of a single subject
 */
typedef struct _SUBJECT_DATA {
    long    lSubjectID;
    long    lNumGaps;

    struct _GAP_DATA {
        long    lID;
        long    lTime;

        long    aalR[MAX_STATES][MAX_STATES];
    } *palGaps;

    long    lCovariateX[MAX_COVARIATES];
} SUBJECT_DATA;

/*
 * Function definitions
 */
extern void LoadData( double *, int *, int *, long *, long * );
extern void PurgeSubjectData( void );

long LocateSubject( long, SUBJECT_DATA ** );

/*
 * Global variables
 * ganSubjects = ORDERED global array of subject data
 * gnMaxSubjects = maximum number of subject data to read from file
 * gdParameters = array of double precision parameter estimates
 * gnMaxParameters = maximum number of parameter estimates
 */
extern SUBJECT_DATA *gaSubjects;
extern long glAllocSubjects;
extern long glNumSubjects;
extern double gdTolerance;

#endif
