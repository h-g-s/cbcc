/********************************************************************************
 * Copyright (c) 2019 Haroldo Gambini Santos
 * 
 * tsp-cuts
 *
 * Example of using the C-API of the COIN-OR CBC MIP solver to solve the 
 * traveling salesmam problem using a branch-and-cut algorithm that separates
 * subtour elimination constraints. 
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0
 *
 ********************************************************************************/

/**
 * @file tsp-cuts.c
 * @author Haroldo Gambini Santos
 * @date 15 Aug 2019
 *
 * Solves the traveling salesman problem using a branch-and-cut algorithm
 * developed with the COIN-OR CBC MIP solver.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>
#include <Cbc_C_Interface.h>
#include "tsp-instance.h"
#include "spaths.h"
#include "mincut.h"


struct DistArc {
    int n1;
    int n2;
    int dist;
};

#define NEW_VECTOR(type, size) ((type *) xmalloc((sizeof(type))*(size)))

static void *xmalloc( const size_t size );

static int compute_fartest_points( const TSPInstance *inst, struct DistArc **da );

static int compDistArcs( const void *p1, const void *p2 );

// checks if a variable is an arc variable in the format x(i,j) and 
// gets the nodes i and j indexes
static char arc_nodes( const char *varName, int *i, int *j );

struct CutAppData 
{
    TSPInstance *inst;
    int nPairs;              // node pairs sorted: most distant first
    struct DistArc *pairs;  
};

static int cutIt = 0;

static void cut_callback( void *osiSolver, void *osiCuts, void *appData ) {
    printf("Starting cut iteration %d\n", cutIt++ );
    fflush(stdout);
    struct CutAppData *caData = appData;
    TSPInstance *inst = caData->inst;

    int maxArcs = tspi_size(inst) * tspi_size(inst);

    int *tail = NEW_VECTOR( int, maxArcs );
    int *head = NEW_VECTOR( int, maxArcs );
    int *cap = NEW_VECTOR( int, maxArcs );
    int nArcs = 0;
    char colName[256];
    char *iv = NEW_VECTOR( char, tspi_size(inst) );

    const double *x = Osi_getColSolution( osiSolver );

    const int nCols = Osi_getNumCols(osiSolver);
    for ( int i=0 ; (i<nCols) ; ++i ) {
        int ai, aj;

        Osi_getColName( osiSolver, i, colName, 256 );

        if ( !arc_nodes( colName, &ai, &aj ) )
            continue;

        int c = (int) (x[i]*10000.0);

        tail[nArcs] = ai;
        head[nArcs] = aj;
        cap[nArcs] = c;
        nArcs++;
    }

    int *idx = NEW_VECTOR( int, nCols );
    double *coef = NEW_VECTOR( double, nCols );

    // checking first conectivity between distant nodes
    int iPair = caData->nPairs -1;
    MinCut *mc = NULL;
    for ( ; iPair >= 0 ; --iPair ) {
        int s = caData->pairs[iPair].n1;
        int t = caData->pairs[iPair].n2;
        assert( s!=t );

        printf("searching min cut %d -> %d\n", s, t); fflush(stdout); fflush(stderr);
    
        if (mc)
            minc_free( &mc );
        mc = minc_create( nArcs, tail, head, cap, s, t );

        int capCut = minc_optimize( mc );

        printf("cap cut: %d\n", capCut);

        if ( (!minc_in_s(mc, s)) ) 
            continue;
        if ( (minc_in_s(mc, t)) ) 
            continue;

        if ( minc_n_cut(mc)==0 ) 
            continue;

        if (capCut == 10000)
            continue;

        double rhs = 0.0;

        memset( iv, 0, sizeof(char)*(tspi_size(inst)) );

        // found cut 
        int nz = 0;
        for ( int i=0 ; (i<nCols) ; ++i ) {
            Osi_getColName( osiSolver, i, colName, 256 );

            int ai, aj;

            if (!arc_nodes(colName, &ai, &aj))
                continue;

            if ( (!minc_in_s(mc, ai)) || (!minc_in_s(mc, aj)) )
                continue;

            if (!iv[ai]) {
                rhs += 1.0;
                iv[ai] = 1;
            }
            if (!iv[aj]) {
                rhs += 1.0;
                iv[aj] = 1;
            }

            idx[nz] = i;
            coef[nz] = 1.0;
            ++nz;
        }

        if (!nz) 
            continue;

        rhs -= 1.0;

        // printing cut
        printf("Adding cut:\n");
        for ( int ic=0 ; (ic<nz) ; ++ic ) {
            Osi_getColName( osiSolver, idx[ic], colName, 256);
            printf("+%g %s ", coef[ic], colName );
        }
        printf("<= %g\n", rhs);
        fflush(stdout);



        OsiCuts_addRowCut( osiCuts, nz, idx, coef, 'L', rhs );

        break;
    }
    minc_free( &mc );
    
    free( iv );
    free( idx );
    free( coef );
    free( tail );
    free( head );
    free( cap );
}

int main(int argc, char **argv)
{
    if (argc<2) 
    {
        fprintf(stderr, "usage: tsp-cbc instanceName");
        exit(1);
    }

    TSPInstance *inst = tspi_create(argv[1]);

    int n = tspi_size( inst );

    Cbc_Model *mip = Cbc_newModel();

    /* creating x binary variables */
    int **x;  // references to variables indexes
    x = malloc(sizeof(int*)*n);
    assert(x);
    x[0] = malloc(sizeof(int)*n*n);
    for ( int i=1 ; (i<n) ; ++i )
        x[i] = x[i-1] + n;

    for ( int i=0 ; (i<n) ; ++i )
    {
        for ( int j=0 ; (j<n) ; ++j )
        {
            x[i][j] = INT_MAX;
            if (i==j)
                continue;
            const int dist = tspi_dist(inst, i, j);
            if (dist==INT_MAX)
                continue;
            char vname[64];
            x[i][j] = Cbc_getNumCols(mip);
            sprintf(vname, "x(%d,%d)", i, j);
            Cbc_addCol(mip, vname, 0.0, 1.0, dist, 1, 0, NULL, NULL);
        }
    }

    /* creating y continuous variables */
    int *y = malloc(sizeof(int)*n);
    for ( int i=0 ; (i<n) ; ++i )
    {
        y[i] = Cbc_getNumCols(mip);
        char vname[64];
        sprintf(vname, "y(%d)", i);
        Cbc_addCol(mip, vname, 0.0, DBL_MAX, 0.0, 0, 0, NULL, NULL);
    }

    /* creating constraints */
    /* variables to store constraint contents */
    int nz;
    int *idx = malloc(sizeof(int)*n);
    double *coef = malloc(sizeof(double)*n);

    /* constraint : enter each city coming from another city */
    for ( int i=0 ; (i<n) ; ++i )
    {
        nz = 0;
        for ( int j=0 ; (j<n) ; ++j )
        {
            if (i==j)
                continue;
            const int dist = tspi_dist(inst, j, i);
            if (dist==INT_MAX)
                continue;
            idx[nz] = x[j][i];
            coef[nz++] = 1.0;
        }

        char rname[64];
        sprintf(rname, "in(%d)", i);
        Cbc_addRow(mip, rname, nz, idx, coef, 'E', 1.0);
    }

    /* constraint : leave each city coming from another city */
    for ( int i=0 ; (i<n) ; ++i )
    {
        nz = 0;
        for ( int j=0 ; (j<n) ; ++j )
        {
            if (i==j)
                continue;
            const int dist = tspi_dist(inst, i, j);
            if (dist==INT_MAX)
                continue;
            idx[nz] = x[i][j];
            coef[nz++] = 1.0;
        }

        char rname[64];
        sprintf(rname, "out(%d)", i);
        Cbc_addRow(mip, rname, nz, idx, coef, 'E', 1.0);
    }
    
    /* subtour elimination (weak initial constraints) */
    for ( int i=1 ; (i<n) ; ++i )
    {
        for ( int j=1 ; (j<n) ; ++j )
        {
            if (i==j)
                continue;
            const int dist = tspi_dist(inst, i, j);
            if (dist==INT_MAX)
                continue;
            nz = 0;
            idx[nz] = y[i];
            coef[nz++] = 1.0;

            idx[nz] = x[i][j];
            coef[nz++] = -(n+1);

            idx[nz] = y[j];
            coef[nz++] = -1.0;

            char rname[64];
            sprintf(rname, "noSub(%d,%d)", i, j);
            Cbc_addRow(mip, rname, nz, idx, coef, 'G', -n );
        }
    }
    
    printf("initial model has %d variables, %d of which are integral and %d rows", 
            Cbc_getNumCols(mip), Cbc_getNumIntegers(mip), Cbc_getNumRows(mip));
    
    struct CutAppData caData;
    caData.inst = inst;
    caData.nPairs = compute_fartest_points( inst, &caData.pairs );

    Cbc_addCutCallback(mip, cut_callback, "Sub-tour elimination", &caData, 1, 0);
    Cbc_solve(mip);

    printf("best route found has length %g, best possible (obj bound) is %g\n", 
            Cbc_getObjValue(mip), Cbc_getBestPossibleObjValue(mip));


    free( x[0] );
    free( x );
    free( y );
    free( idx );
    free( coef );
    free( caData.pairs );

    /* free cbc model */
    Cbc_deleteModel(mip);

    tspi_free(inst);
}

static int compDistArcs( const void *p1, const void *p2 ) {
    const struct DistArc *da1 = p1;
    const struct DistArc *da2 = p2;

    return da1->dist - da2->dist;
}

static int compute_fartest_points( const TSPInstance *inst, struct DistArc **da ) {
    ShortestPathsFinder *spf = spf_create();
    int *start, *to, *dist;
    int maxArcs = tspi_size(inst)*tspi_size(inst);
    start = NEW_VECTOR( int, (tspi_size(inst)+1) );
    to = NEW_VECTOR( int, maxArcs );
    dist = NEW_VECTOR( int, maxArcs );
    start[0] = 0;

    *da = NEW_VECTOR( struct DistArc, maxArcs );
    int na = 0;
    for ( int i=0 ; (i<tspi_size(inst)) ; ++i ) {
        for ( int j=0 ; (j<tspi_size(inst)) ; ++j ) {
            if ( i==j )
                continue;
            if ( tspi_dist( inst, i, j ) == INT_MAX )
                continue;

            to[na] = j;
            dist[na] = tspi_dist( inst, i, j);
            na++;
        }
        start[i+1] = na;
    }

    spf_update_graph( spf, tspi_size( inst ), na, start, to, dist );

    na = 0;

    // computing pairs of distant nodes
    for ( int i=0 ; (i<tspi_size(inst)) ; ++i ) {
        spf_find( spf, i );
        for ( int j=i+1 ; j<tspi_size(inst) ; ++j ) {

            if ( spf_get_dist(spf,  j) != SP_INFTY_DIST ) {
                (*da)[na].n1 = i;
                (*da)[na].n2 = j;
                (*da)[na].dist = spf_get_dist(spf,  j);
                ++na;
            }
        }
    }

    spf_free(&spf);
    free(start);
    free(to);
    free(dist);

    qsort( *da, na, sizeof(struct DistArc), compDistArcs );

    return na;
}


static char arc_nodes( const char *varName, int *i, int *j ) {
    char *s = strstr(varName, "x(");

    if (!s)
        return 0;

    s += 2;
    assert( isdigit(*s) );

    char *s2 = strstr( s, "," );
    *s2 = '\0';

    s2++;
    assert( isdigit(*s2) );

    char *s3 = strstr( s2, ")" );
    assert(s3);
    *s3 = '\0';

    *i = atoi(s);
    *j = atoi(s2);

    return 1;
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

