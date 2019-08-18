/********************************************************************************
 * Copyright (c) 2019 Haroldo Gambini Santos
 * 
 * queens
 * Example of using the C-API of the COIN-OR CBC MIP solver to solve the 
 * n-queens problem using lazy constraints (added only when an integer solution
 * is found using callbacks)
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0
 *
 ********************************************************************************/

/**
 * @file queens-lazy.c
 * @author Haroldo Gambini Santos
 * @date 15 Aug 2019
 *
 * Solves the n-queens problem using a compact MIP formulation. 
 * Communication with the COIN-OR CBC MIP solver is done with its C API.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <Cbc_C_Interface.h>

void queens_cut_callback(void *osiSolver, void *osiCuts, void *appdata);

int n;

int **x;


int *idx;
double *coef;

/* computes the lef-hand-side of a constraint, i.e. the value
 * obtained replacing the variables with the values of the current 
 * solution */
double lhs( const double *sol, int nz, int idx[], double coef[] );

int main( int argc, char **argv )
{
    const double *xs;
    int i, j, k;

    if (argc<2)
    {
        fprintf(stderr, "usage: queens numberOfQueens\n");
        exit(1);
    }

    n = atoi(argv[1]);

    x = (int **) malloc( sizeof(int*)*n );
    x[0] = (int *) malloc( sizeof(int)*n*n );
    for ( i=1 ; (i<n) ; ++i )
        x[i] = x[i-1] + n;

    Cbc_Model *model = Cbc_newModel();
    
    /* adding variables */
    k = 0;
    for ( i=0 ; (i<n) ; ++i )
    {
        for ( j=0 ; (j<n) ; ++j )
        {
            char name[256];
            x[i][j] = k++;
            sprintf(name, "x(%d,%d)", i, j);
            Cbc_addCol(model, name, 0.0, 1.0, 0.0, 1, 0, NULL, NULL);
        }
    }

    // area to store constraint contents
    idx = (int*) malloc(sizeof(int)*n);
    coef = (double *) malloc(sizeof(double)*n);

    /* constraint one per row */
    for ( i=0 ; (i<n) ; ++i )
    {
        char name[256];
        for ( j=0 ; j<n ; ++j )
        {
            idx[j] = x[i][j];
            coef[j] = 1.0;
        }
        sprintf(name, "row(%d)", i);
        Cbc_addRow(model, name, n, idx, coef, 'E', 1.0);
    }

    /* constraint one per column */
    for ( j=0 ; (j<n) ; ++j )
    {
        char name[256];
        for ( i=0 ; i<n ; ++i )
        {
            idx[i] = x[i][j];
            coef[i] = 1.0;
        }
        sprintf(name, "col(%d)", j);
        Cbc_addRow(model, name, n, idx, coef, 'E', 1.0);
    }

    Cbc_addCutCallback( model, queens_cut_callback, "diagonalConstraints", NULL, 1, 1);

    Cbc_setMaximumSeconds(model, 100);
    // to keep original variable indexes valid
    Cbc_setParameter(model, "preprocess", "off");
    Cbc_setParameter(model, "heur", "off");
    Cbc_solve(model);
    xs = Cbc_getColSolution(model);
    if (xs) {
        for ( i=0 ; (i<n) ; ++i ) {
            for ( j=0 ; (j<n) ; ++j ) {
                printf("%c ", (xs[x[i][j]]>=0.99 ? '*' : '.'));
            }
            printf("\n");
        }
    }

    free(idx);
    free(coef);
    free(x[0]);
    free(x);

    Cbc_deleteModel(model);
}

double lhs( const double *sol, int nz, int idxc[], double coefc[] ) {
    double res = 0.0;
    for ( int i=0 ; (i<nz) ; ++i )
        res += sol[idxc[i]]*coefc[i];

    return res;
}

void queens_cut_callback(void *osiSolver, void *osiCuts, void *appdata) 
{
    const double *sol = Osi_getColSolution(osiSolver);

    /* diagonal  */
    int i, j,  k;
    int p = 0;
    for ( k=2-n ; k<(n-1) ; ++k, ++p )
    {
        int nz = 0;
        for ( i=0 ; (i<n) ; ++i )
        {
            for ( j=0 ; (j<n) ; ++j )
            {
                if (i-j==k)
                {
                    idx[nz] = x[i][j];
                    coef[nz] = 1.0;
                    ++nz;
                }
            }
        }
        if (lhs(sol, nz, idx, coef) >= 1.0001) {
            OsiCuts_addRowCut( osiCuts, nz, idx, coef, 'L', 1.0 );
            printf("> added lazy constraint\n"); fflush(stdout);
        }
    }

    /* diagonal */
    p = 0;
    for ( k=3 ; k<(n+n) ; ++k, ++p )
    {
        int nz = 0;
        for ( i=0 ; (i<n) ; ++i )
        {
            for ( j=0 ; (j<n) ; ++j )
            {
                if (i+j==k)
                {
                    idx[nz] = x[i][j];
                    coef[nz] = 1.0;
                    ++nz;
                }
            }
        }
 
        if (lhs(sol, nz, idx, coef) >= 1.0001) {
            OsiCuts_addRowCut( osiCuts, nz, idx, coef, 'L', 1.0 );
            printf("> added lazy constraint\n"); fflush(stdout);
        }
    } 
}

