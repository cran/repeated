/*
 * main.c
 * Defines the main application.
 *
 * Author:  RPS
 *
 * Created: Jan/96
 * Copyright (c) 1999 by Richard J. Cook, University of Waterloo
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "gaps.h"
#include "calcs.h"

/*
 * Initialize global variables
 */
SUBJECT_DATA *gaSubjects = NULL;
long glAllocSubjects = 0l;
long glNumSubjects = 0l;
double gdTolerance = TOLERANCE;

/*
 * LoadData
 * Load the treatment-control subject data file
 *
 * pszFileName = [I] data filename
 */
void
LoadData( double *response, int *nrow, int *nSize,
          long * plNumSubjects, 
          long *plReturn )
{
    double alVector[8];
    SUBJECT_DATA *pSubject = NULL;
    int iCov;
    int i, j, k;
    long lIX;
    BOOL bNewSubject;

    *plReturn = ID_NO_ERROR;

    PurgeSubjectData();

    /*
     * load subject data structure
     */
    glAllocSubjects = SIZE_SUBJECT_BLOCKS;
    gaSubjects = calloc( SIZE_SUBJECT_BLOCKS, sizeof( SUBJECT_DATA ));
    if ( gaSubjects == NULL ) {
        *plReturn = ID_ERR_MEMORY;
        goto OnError;
    }
    
    k=0;
    for ( glNumSubjects = 0; k<*nrow; ) {

	for(i=0;i<*nSize;i++)alVector[i]=response[i+k* *nSize];
	k++;

        /* ensure subject array is large enough to add another subject */
        if ( glNumSubjects >= glAllocSubjects ) {
            glAllocSubjects += SIZE_SUBJECT_BLOCKS;
            gaSubjects = realloc( gaSubjects, 
                                  glAllocSubjects * sizeof( SUBJECT_DATA ));
        }

            /* locate the subject with the matching subject ID since subject 
               data may continue on multiple lines */
            *plReturn = LocateSubject( alVector[COL_SUBJECT_ID], &pSubject );
            if ( *plReturn != ID_NO_ERROR ) {
                /* error occurred while searching for subject */
                goto OnError;
            }

            if ( pSubject == NULL ) {
                /* subject not found - use a new element in the subject array */
                pSubject = &( gaSubjects[glNumSubjects++] );
                bNewSubject = TRUE; } else bNewSubject = FALSE;
            /* assign the subject elements */
            pSubject->lSubjectID = alVector[COL_SUBJECT_ID];
            if ( bNewSubject ) {
                pSubject->lNumGaps = 1;
                pSubject->palGaps = calloc( pSubject->lNumGaps, 
                                            sizeof( struct _GAP_DATA ));
            } else {
                pSubject->lNumGaps++;
                pSubject->palGaps = realloc( pSubject->palGaps, 
                                             sizeof( struct _GAP_DATA ) * 
                                             pSubject->lNumGaps );
            }

            if ( pSubject->palGaps == NULL ) {
                *plReturn = ID_ERR_MEMORY;
                goto OnError;
            }

            pSubject->palGaps[pSubject->lNumGaps - 1]
                      .lID = alVector[COL_GAP_ID];
            pSubject->palGaps[pSubject->lNumGaps - 1]
                      .lTime = alVector[COL_GAP_TIME];
            for ( i = 0; i < MAX_STATES; i++ ) {
                for ( j = 0; j < MAX_STATES; j++ ) {
                    lIX = COL_RESPONSES + j + i * MAX_STATES;
                    pSubject->palGaps[pSubject->lNumGaps - 1]
                              .aalR[i][j] = alVector[lIX];
                }
            }

            if ( *nSize > MIN_NUM_COLUMNS ) {
                for ( iCov = 0; iCov < MAX_COVARIATES; iCov++ )
                    pSubject->lCovariateX[iCov] = alVector[COL_COVARIATE + iCov];
            } else {
                for ( iCov = 0; iCov < MAX_COVARIATES; iCov++ )
                    pSubject->lCovariateX[iCov] = 0l;
			}
        }


    /* clean up and exit */
OnError:
    if ( *plReturn != ID_NO_ERROR ) {
        /* clean-up if error occurred */
        PurgeSubjectData();
    }

    *plNumSubjects = glNumSubjects;
}

long
LocateSubject( long lSubjectID, SUBJECT_DATA **ppSubject )
{
    long l;

    *ppSubject = NULL;

    if ( gaSubjects == NULL )
        return( ID_ERR_MEMORY );

    for ( l = 0;( l < glNumSubjects ) && ( *ppSubject == NULL ); l++ ) {
        if ( gaSubjects[l].lSubjectID == lSubjectID )
            *ppSubject = &( gaSubjects[l] );
    }

    return( ID_NO_ERROR );
}

void
PurgeSubjectData( void )
{
    long l;

    if ( gaSubjects == NULL ) {
        glNumSubjects = 0;
        glAllocSubjects = 0;
        return;
    }

    for ( l = 0; l < glNumSubjects; l++ ) {
        if ( gaSubjects[l].palGaps )
            free( gaSubjects[l].palGaps );
    }
    free( gaSubjects );
    gaSubjects = NULL;
    glNumSubjects = 0;
    glAllocSubjects = 0;
}
