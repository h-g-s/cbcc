#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "tsp-instance.h"
#include <Cbc_C_Interface.h>

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
    
    /* subtour elimination */
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

    //printf("model has %d variables, %d of which are integral and %d rows", 
            //Cbc_getNumCols(mip), Cbc_getNumIntegers(mip), Cbc_getNumRows(mip));

    Cbc_setMaximumSeconds(mip, 120);
    Cbc_writeLp(mip, "tsp");
    Cbc_solve(mip);

    printf("best route found has length %g, best possible (obj bound) is %g\n", 
            Cbc_getObjValue(mip), Cbc_getBestPossibleObjValue(mip));


    free(x[0]);
    free(x);
    free(idx);
    free(coef);

    /* free cbc model */
    Cbc_deleteModel(mip);

    tspi_free(inst);
}

