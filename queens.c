/********************************************************************************
 * Copyright (c) 2019 Haroldo Gambini Santos
 * 
 * queens
 * Example of using the C-API of the COIN-OR CBC MIP solver to solve the 
 * n-queens problem
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0
 *
 ********************************************************************************/

/**
 * @file queens.c
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

int main( int argc, char **argv )
{
    int *idx;
    double *coef;
    const double *xs;
    int i, j, k, p;

    if (argc<2)
    {
        fprintf(stderr, "usage: queens numberOfQueens\n");
        exit(1);
    }

    int n = atoi(argv[1]);

    int **x = malloc( sizeof(int*)*n );
    x[0] = malloc( sizeof(int)*n*n );
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
    idx = malloc(sizeof(int)*n);
    coef = malloc(sizeof(double)*n);

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

    /* diagonal  */
    p = 0;
    for ( k=2-n ; k<(n-1) ; ++k, ++p )
    {
        char name[256];
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
        sprintf(name, "diag1(%d)", k);
        char *s = name;
        while (*s != '\0') {
            if (*s == '-')
                *s = 'm';
            ++s;
        }
        Cbc_addRow(model, name, nz, idx, coef, 'L', 1.0);
    }

    /* diagonal */
    p = 0;
    for ( k=3 ; k<(n+n) ; ++k, ++p )
    {
        char name[256];
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
        sprintf(name, "diag2(%d)", k);
        char *s = name;
        while (*s != '\0') {
            if (*s == '-')
                *s = 'm';
            ++s;
        }
 
        Cbc_addRow(model, name, nz, idx, coef, 'L', 1.0);
    } 

    Cbc_setMaximumSeconds(model, 100);
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

