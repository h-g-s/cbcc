/********************************************************************************
 * Copyright (c) 2019 Haroldo Gambini Santos
 * 
 * rcpsp
 *
 * Example of using the C-API of the COIN-OR CBC MIP solver to solve the 
 * solve the Resource Constrained Project Scheduling Problem - RCPSP
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0
 *
 ********************************************************************************/

/**
 * @file rcpsp.c
 * @author Haroldo Gambini Santos
 * @date 4 Sep 2019
 *
 * Solves the Resource Constrained Project Scheduling Problem,
 * developed with the COIN-OR CBC MIP solver.
 */

#include <stdio.h>
#include <stdlib.h>
#include "Cbc_C_Interface.h"

// resources
const int r = 2;

// processing times
const int p[] =
    { 0, 3, 2, 5, 4, 2, 3, 4, 2, 4, 6, 0 };

// jobs
const int n = sizeof(p) / sizeof(int);

// resource usage
const int u[][2] = {
                    {0, 0},
                    {5, 1},
                    {0, 4},
                    {1, 4},
                    {1, 3},
                    {3, 2},
                    {3, 1},
                    {2, 4},
                    {4, 0},
                    {5, 2},
                    {2, 5},
                    {0, 0},
                   };

// precedence constraints
const int P[][2] = {
                     {0, 1},
                     {0, 2},
                     {0, 3},
                     {1, 4},
                     {1, 5},
                     {2, 10},
                     {2, 9},
                     {3, 8},
                     {4, 6},
                     {4, 7},
                     {5, 10},
                     {5, 9},
                     {6, 9},
                     {6, 8},
                     {7, 8},
                     {10, 11},
                     {9, 11},
                     {8, 11}
                   };

const int c[] = {6, 8};

const int nPrec = sizeof(P) / (sizeof(int)*2);

#define NEW_VECTOR(type, size) ((type *) xmalloc((sizeof(type))*(size)))
static void *xmalloc( const size_t size );


int main( int argc, const char **argv ) 
{
    // computing time horizon
    int maxT = 0;
    for ( int j=0 ; (j<n) ; ++j )
        maxT += p[j];

    int *est = NEW_VECTOR( int, n );
    int *lst = NEW_VECTOR( int, n );

    // EST
    for ( int j=0 ; (j<n) ; ++j )
        est[j] = 0;
    for ( int j=0 ; (j<n) ; ++j )
        for ( int i=0 ; (i<nPrec) ; ++i ) 
            est[P[i][1]] = est[P[i][0]] + p[P[i][0]];
    // LST
    for ( int j=0 ; (j<n) ; ++j )
        lst[j] = maxT;
    for ( int j=0 ; (j<n) ; ++j )
        for ( int i=0 ; (i<nPrec) ; ++i ) 
            lst[P[i][0]] = lst[P[i][1]] - p[P[i][1]];
    
    int **x = NEW_VECTOR( int *, n );
    x[0] = NEW_VECTOR( int, n*maxT );
    for ( int j=1 ; (j<n) ; ++j )
        x[j] = x[j-1] + maxT;
    for ( int i=0 ; (i<n*maxT) ; ++i )
        x[0][i] = -1;

    Cbc_Model *mip = Cbc_newModel();

    // creating x vars
    for ( int j=0 ; (j<n) ; ++j )
    {
        for ( int t=est[j] ; t<lst[j];  ++t )
        {
            x[j][t] = Cbc_getNumCols(mip);
            char colName[256];
            sprintf( colName, "x(%d,%d)", j, t );
            double obj = ( j == n-1 ) ? t : 0.0;
            Cbc_addCol( mip, colName, 0.0, 1.0, obj, 1, 0, NULL, NULL);
        }
    }

    int *idx = NEW_VECTOR( int, Cbc_getNumCols(mip));
    double *coef = NEW_VECTOR( double, Cbc_getNumCols(mip));

    // select job allocation time
    for ( int j=0 ; (j<n) ; ++j ) {
        int nz = 0;
        for ( int t=est[j] ; t<lst[j];  ++t )
        {
            idx[nz] = x[j][t];
            coef[nz] = 1.0;
            ++nz;
        }

        char rowName[256];
        sprintf( rowName, "allocate(%d)", j );
        Cbc_addRow( mip, rowName, nz, idx, coef, 'E', 1.0 );
    }

    // resource usage
    for ( int t=0 ; (t<maxT) ; t++ ) 
    {
        for ( int k=0 ; (k<r) ; ++k )
        {
            int nz = 0;
            // jobs
            for ( int j=0 ; j<n ; ++j )
            {
                for ( int tl=t ; tl>=t-p[j]+1; --tl ) {
                    if ( x[j][tl] < 0 || u[j][k] == 0 )
                        continue;

                    idx[nz] = x[j][tl];
                    coef[nz] = u[j][k];
                    ++nz;

                }
            }
            char rowName[256];
            sprintf( rowName, "cap(%d,%d)", k, t );

            Cbc_addRow( mip, rowName, nz, idx, coef, 'L', c[k]);
        }
    }

    // precedence
    for ( int ip=0 ; ip<nPrec ; ++ip ) 
    {
        int j1 = P[ip][0];
        int j2 = P[ip][1];

        int nz = 0;
        for ( int t=est[j2] ; t<lst[j2] ; ++t ) {
            idx[nz] = x[j2][t];
            coef[nz] = t;
            ++nz;
        }
        for ( int t=est[j1] ; t<lst[j1] ; ++t ) {
            idx[nz] = x[j1][t];
            coef[nz] = -t;
            ++nz;
        }
        char rowName[256];
        sprintf( rowName, "prec(%d,%d)", j1, j2 );
        Cbc_addRow( mip, rowName, nz, idx, coef, 'G', p[j1] );
    }

    Cbc_writeLp( mip, "a");
    
    Cbc_deleteModel( mip );
    
    free( idx );
    free( coef );
    free( est );
    free( lst );
    free( x[0] );
    free( x );

    return 0;
}

static void *xmalloc( const size_t size )
{
   void *result = malloc( size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      abort();
   }

   return result;
}

