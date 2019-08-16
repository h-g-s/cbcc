#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "tsp-instance.h"

#define PI 3.141592
static const double RRR = 6378.388;

#define nint(x) ((int)((x) + 0.5))

static double rad(double x)
{    
    return x*PI/180.0;
}

struct _TSPInstance
{
    int size;

    char *name;

    int **d;
};

static void compute_distances(TSPInstance *inst, double **coord) 
{
    int **d = inst->d;
    for (int i=0 ; i<inst->size; ++i )
    {
        d[i][i] = 0;
        for ( int j=0 ; j<inst->size ; ++j )
        {
               double latitude_i = coord[i][0];
               double latitude_j = coord[j][0];
               double longitude_i = coord[i][1];
               double longitude_j = coord[j][1];
               double q1 = cos(longitude_i - longitude_j);
               double q2 = cos(latitude_i - latitude_j);
               double q3 = cos(latitude_i + latitude_j);
               d[i][j] = (int)( RRR * acos( 0.5 * ( ( 1.0 + q1 ) * q2 - (1.0 - q1) * q3 ) ) + 1.0 );
        }
    }
}

TSPInstance *tspi_create(const char fileName[])
{
    TSPInstance *inst = calloc( sizeof(TSPInstance), 1 );
#define LINE_SIZE 512
    char line[LINE_SIZE], *s;
    FILE *f = fopen(fileName, "r");
    assert(f);
    char readingCoord = 0;

    double **coord = NULL;
    while ( (s=fgets(line, LINE_SIZE, f)) ) 
    {
        if ( readingCoord ) 
        {
            if (strstr(s, "EOF"))
                break;

            int i; double x,y;
            int nr = sscanf(s, "%d %lf %lf", &i, &x, &y);

            assert(nr == 3);
            coord[i-1][0] = rad(x);
            coord[i-1][1] = rad(y);
        }
        else
        {
            if (strstr(s, "NAME:")) 
            {
                char *p, *t = NULL;
                p = strtok_r(s, " ", &t);
                assert( p != NULL );
                p = strtok_r(NULL, " ", &t);
                int l = strlen(p);
                if (p[l-1]=='\n') 
                {
                    p[l-1]='\0';
                    --l;
                }
                inst->name = malloc(l+1);
                strcpy(inst->name, p);
                printf("instance name: %s\n", inst->name);
            }
            if (strstr(s, "DIMENSION:")) 
            {
                char *p, *t = NULL;
                p = strtok_r(s, " ", &t);
                assert( p != NULL );
                p = strtok_r(NULL, " ", &t);
                assert( p != NULL );
                inst->size = atoi(p);
                printf("instance size: %d\n", inst->size);

                coord = malloc(sizeof(double*)*inst->size);
                coord[0] = malloc(sizeof(double)*inst->size*2);
                for ( int i=1 ; (i<inst->size) ; ++i )
                    coord[i] = coord[i-1]+2;

                inst->d = malloc(sizeof(int*)*inst->size);
                assert(inst->d);
                inst->d[0] = malloc(sizeof(int)*inst->size*inst->size);
                for ( int i=1 ; (i<inst->size) ; ++i )
                    inst->d[i] = inst->d[i-1] + inst->size;
            }
            if (strstr(s, "NODE_COORD_SECTION"))
            {
                readingCoord = 1;
            }
        }
    }

    compute_distances(inst, coord);

    /*
    printf("distance matrix:\n");
    for ( int i=0 ; (i<inst->size) ; ++i )
    {
        for ( int j=0 ; (j<inst->size) ; ++j )
        {
            printf("%d ", inst->d[i][j]);
        }
        printf("\n");
    } */
    
    if (coord)
    {
        free(coord[0]);
        free(coord);
    }

    fclose(f);
    return inst;
#undef LINE_SIZE
}

int tspi_size(const TSPInstance *tspi)
{
    return tspi->size;
}

int tspi_dist(const TSPInstance *tspi, int i, int j)
{
    return tspi->d[i][j];
}

void tspi_free(TSPInstance *tspi)
{
    if (tspi->name)
        free(tspi->name);
    free(tspi->d[0]);
    free(tspi->d);
    free(tspi);
}


