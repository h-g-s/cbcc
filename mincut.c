/********************************************************************************
 * Copyright (c) 2019 Haroldo Gambini Santos
 * 
 * MinCut
 * A self-contained, C99 compatible minimum cut implementation 
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0
 *
 ********************************************************************************/


/* auxiliary code */

#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "mincut.h"

/* macros */
#define INT_RANDOM( n ) \
   ((int) floor( ((double)(n)) * (((double)rand())/(((double)RAND_MAX)+((double)1.0))) ))

/* generates a random number between 0 and 1 */
#define DBL_RANDOM( ) \
   ((double) rand()) / (((double)RAND_MAX))

#define MAX( a, b ) ( (a)>(b) ? (a) : (b) )
#define MIN( a, b ) ( (a)<(b) ? (a) : (b) )

#define FILE_NAME_SIZE 1024
#define STR_SIZE        256
#define LSTR_SIZE       512
#define LINE_SIZE      2048

#define ALLOCATE( ptr, type ) {\
    ptr = (type*) malloc( sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_INI( ptr, type ) {\
    ptr = calloc( 1, sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

#define ALLOCATE_VECTOR( ptr, type, nElements ) {\
    ptr = (type*) malloc( sizeof(type)*(nElements) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }

/* allocate filling with zeros */
#define ALLOCATE_VECTOR_INI( ptr, type, nElements ) {\
    ptr = (type*) calloc( (nElements), sizeof(type) ); \
    if (!ptr) { \
       fprintf( stderr, "ERROR: no more memory available. at: %s:%d\n", __FILE__, __LINE__ ); \
       abort(); \
    }; }


#define OPEN_FILE( f, fileName, mode )  \
    f = fopen( fileName, mode ); \
    if (!f) { \
        fflush(stdout); \
        fprintf( stderr, "ERROR: could not open file %s with mode %s. At %s:%d\n", fileName, mode, __FILE__, __LINE__ ); \
        fflush(stderr);  abort(); exit(EXIT_FAILURE); }

#define True  1
#define False 0

/* fills from start until the last element before end */
#define FILL( vector, start, end, value ) { \
    int i; \
    for ( i=start ; (i<end) ; ++i ) vector[i] = value; \
} \

#define EPS 1e-5

#define DBL_EQUAL( v1, v2 ) ( fabs(v1-v2)<=EPS )

static void *xmalloc( const size_t size );

static void *xcalloc( const size_t elements, const size_t size );

static void *xmalloc( const size_t size )
{
   void *result = malloc( size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
      exit(1);
   }

   return result;
}

static void *xcalloc( const size_t elements, const size_t size )
{
   void *result = calloc( elements, size );
   if (!result)
   {
      fprintf(stderr, "No more memory available. Trying to callocte %zu bytes.", size);
      exit(1);
   }

   return result;
}

typedef struct
{
    int a;
    int b;
} IntPair;

#define NOT_FOUND INT_MAX


typedef struct { 
    char str[256]; 
    int value; 
    int keyPos; 
} Dict_Bucket_int; 

typedef struct _StrV StrV;
static StrV* strv_create(int strSize);
static void strv_push_back(StrV* strv, const char* str);
static int strv_size(const StrV* strv);
static void strv_free(StrV** strv);

struct _StrV {
    int capacity;
    int size;
    char** sv;
    char* s;
    int strSize;
};


StrV* strv_create(int strSize)
{
    StrV* res = (StrV*)xmalloc(sizeof(StrV));
    res->capacity = 64;
    res->size = 0;
    res->strSize = strSize;
    res->sv = (char**)xmalloc(sizeof(char*) * res->capacity);
    res->s = (char*)xcalloc(res->capacity * strSize, sizeof(char));
    res->sv[0] = res->s;
    int i;
    for(i = 1; (i < res->capacity); i++)
        res->sv[i] = res->sv[i - 1] + strSize;
    return res;
}


static int strv_size(const StrV* strv)
{
    return strv->size;
}

static void strv_increase_capacity_to(StrV* strv, int newCapacity)
{
    char** sv = (char**)xmalloc(sizeof(char*) * newCapacity);
    char* s = (char*)xcalloc(newCapacity * strv->strSize, sizeof(char));
    sv[0] = s;
    int i;
    for(i = 1; (i < newCapacity); i++)
        sv[i] = sv[i - 1] + strv->strSize;
    memcpy(s, strv->s, sizeof(char) * strv->size * strv->strSize);
    free(strv->sv);
    free(strv->s);
    strv->sv = sv;
    strv->s = s;
    strv->capacity = newCapacity;
}


static void strv_push_back(StrV* strv, const char* str)
{
    if(strv->size + 1 > strv->capacity)
        strv_increase_capacity_to(strv, strv->capacity * 2);
    strncpy(strv->sv[strv->size], str, strv->strSize);
    strv->size++;
}

static void strv_free(StrV** strv)
{
    free((*strv)->s);
    free((*strv)->sv);
    free((*strv));
    (*strv) = NULL;
}

typedef struct { 
    Dict_Bucket_int **cont; 
    int *rowSize; 
    int *rowCap; 
    int defValue; 
    unsigned int 
    hashSize; 
    StrV *keys; 
} Dict_int; 

static Dict_int *dict_int_create( unsigned int hashSize, int defaultValue ); 

static void dict_int_set( Dict_int *dict, const char *key, const int value ); 

static int dict_int_get( const Dict_int *dict, const char *str ); 

static void dict_int_free( Dict_int **dict );

Dict_int* dict_int_create(unsigned int hashSize, int defaultValue)
{
    Dict_int* res;
    {
        res = (Dict_int*)malloc(sizeof(Dict_int));
        if(!res) {
            fprintf(stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 49);
            abort();
        };
    };
    {
        res->cont = (Dict_Bucket_int**)calloc((hashSize), sizeof(Dict_Bucket_int*));
        if(!res->cont) {
            fprintf(stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 49);
            abort();
        };
    };
    {
        res->rowSize = (int*)calloc((hashSize), sizeof(int));
        if(!res->rowSize) {
            fprintf(stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 49);
            abort();
        };
    };
    {
        res->rowCap = (int*)calloc((hashSize), sizeof(int));
        if(!res->rowCap) {
            fprintf(stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 49);
            abort();
        };
    };
    res->defValue = defaultValue;
    res->keys = strv_create(256);
    res->hashSize = hashSize;
    return res;
}

static const unsigned int hashval[] = { 11, 269, 3, 7, 31, 37, 131, 13, 17, 647, 653, 89, 97, 101, 39, 149, 151, 157,
    821, 257, 263, 389, 397, 457, 461, 463, 331, 337, 347, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
    9, 53, 59 };
    
static const unsigned int nHashvalues = sizeof(hashval) / sizeof(int);

unsigned int str_hash(const char* str, const unsigned int hashSize)
{
    const unsigned int len = (unsigned int)((nHashvalues) < (strlen(str)) ? (nHashvalues) : (strlen(str)));
    unsigned int sum = 0, i;
    for(i = 0; (i < len); i++)
        sum += (unsigned int)str[i] * hashval[i];
    return sum % hashSize;
}

static void dict_int_set(Dict_int* dict, const char* key, const int value)
{
    unsigned int hashPos = str_hash(key, dict->hashSize);
    int i;
    char found = 0;
    for(i = 0; (i < dict->rowSize[hashPos]); i++) {
        if(strcmp(key, dict->cont[hashPos][i].str) == 0) {
            found = 1;
            break;
        }
    }
    if(found) {
        dict->cont[hashPos][i].value = value;
        return;
    }
    if(dict->rowSize[hashPos] + 1 > dict->rowCap[hashPos]) {
        if(dict->rowSize[hashPos] == 0) {
            {
                dict->cont[hashPos] = (Dict_Bucket_int*)malloc(sizeof(Dict_Bucket_int) * (64));
                if(!dict->cont[hashPos]) {
                    fprintf(stderr, "ERROR: no more memory available. at: %s:%d\n", "containers.c", 49);
                    abort();
                };
            };
            dict->rowCap[hashPos] = 64;
        } else {
            dict->rowCap[hashPos] *= 2;
            Dict_Bucket_int* bigger =
                (Dict_Bucket_int*)realloc(dict->cont[hashPos], sizeof(Dict_Bucket_int) * dict->rowCap[hashPos]);
            if(!bigger) {
                fprintf(stderr, "ERROR: no more memory available");
                abort();
                exit(1);
            }
            dict->cont[hashPos] = bigger;
        }
    }
    strncpy(dict->cont[hashPos][dict->rowSize[hashPos]].str, key, 256);
    dict->cont[hashPos][dict->rowSize[hashPos]].value = value;
    dict->cont[hashPos][dict->rowSize[hashPos]].keyPos = strv_size(dict->keys);
    dict->rowSize[hashPos]++;
    strv_push_back(dict->keys, key);
}

static int dict_int_get(const Dict_int* dict, const char* str)
{
    unsigned int hashPos = str_hash(str, dict->hashSize);
    int i;
    for(i = 0; (i < dict->rowSize[hashPos]); i++) {
        if(strcmp(str, dict->cont[hashPos][i].str) == 0) {
            return dict->cont[hashPos][i].value;
            break;
        }
    }
    return dict->defValue;
}

static void dict_int_free(Dict_int** dict)
{
    Dict_int* d = *dict;
    unsigned int i;
    for(i = 0; (i < d->hashSize); ++i)
        if(d->rowCap[i])
            free(d->cont[i]);
    free((*dict)->cont);
    free((*dict)->rowSize);
    free((*dict)->rowCap);
    strv_free(&(*dict)->keys);
    free(*dict);
}

/* min cut code */
struct MinCArc
{
    int v;
    int cap;
    int rpos; // position of reverse arc
    char original; // if arc exists in original graph
};

struct _MinCut
{
    int n;

    // start[u] indicates the starting position outbond arcs from u
    int *start;

    struct MinCArc *arcs;

    // maps to original node indexes,
    // just to output the solution
    int *orig;
    
    // maps original indexes
    // to new indexes
    int *newIdx;

    // incidence vector and list of visited nodes
    char *ivVisited;
    int nVisited;
    int *visited;
    
    // queue of unvisited nodes
    int *queue;

    // path from t to s
    int *parent;

    // to quickly search for an (u,v) arc
    Dict_int *arcDict;
    
    // source node
    int s;

    // target node
    int t;

    // minimum cut answer
    int nCut;
    int *cutU;
    int *cutV;
};

static char *arcName( char *str, int u, int v )
{
    sprintf( str, "(%d,%d)",u,v );
    return str;
}

static int arcPos( const MinCut *minc, int u, int v )
{
    char str[256];
    const Dict_int *arcDict = minc->arcDict;
    const char *aname = arcName( str, u, v );
    return dict_int_get( arcDict, aname );
}

MinCut *minc_create( int nArcs, const int _tail[], const int _head[], const int _cap[], int s, int t )
{
    assert( s!=t );

    int maxN = -1;
    for ( int i=0 ; (i<nArcs) ; ++i )
        maxN = MAX( maxN, _tail[i] );
    for ( int i=0 ; (i<nArcs) ; ++i )
        maxN = MAX( maxN, _head[i] );

    // maps original nodes to new pre-processed nodes
    int *ppnode;
    ALLOCATE_VECTOR( ppnode, int, (maxN+1)  );
    for ( int i=0 ; (i<maxN+1) ; ++i )
        ppnode[i] = -1;

    int n=0;

    for ( int i=0 ; (i<nArcs) ; ++i )
        if (ppnode[_tail[i]]==-1)
            ppnode[_tail[i]] = n++;
    for ( int i=0 ; (i<nArcs) ; ++i )
        if (ppnode[_head[i]]==-1)
            ppnode[_head[i]] = n++;

    // may be larger due to insertion of reverse arcs
    int *tail;
    int *head;
    int *cap;
    ALLOCATE_VECTOR( tail, int, (6*nArcs) );
    head = tail + 2*nArcs;
    cap = head + 2*nArcs;

    // inverse arcs must be added if needed
    Dict_int *arcDict = dict_int_create( nArcs*2, -1 );

    // storing positions of existing arcs
    char str[256];
    for ( int i=0 ; i<nArcs ; ++i )
    {
        tail[i] = ppnode[_tail[i]];
        head[i] = ppnode[_head[i]];
        if (tail[i] == head[i] )
        {
            fprintf( stderr, "minc ERROR: arc (%d,%d) specified (self arc).\n", _tail[i], _head[i] );
            abort();
        }
        cap[i] = _cap[i];
        const char *aname = arcName(str, tail[i], head[i] );
        if ( dict_int_get( arcDict, aname ) != -1 )
        {
            fprintf( stderr, "minc ERROR: arc (%d,%d) specified twice.\n", _tail[i], _head[i] );
            abort();
        }
        dict_int_set( arcDict, aname, i );
    }

    int *orig;
    ALLOCATE_VECTOR( orig, int, n );

    
    for ( int i=0 ; (i<maxN+1) ; ++i )
        if (ppnode[i]!=-1)
            orig[ppnode[i]] = i;

    MinCut *minc;
    ALLOCATE( minc, MinCut );

    minc->s = ppnode[s];
    minc->t = ppnode[t];
    minc->newIdx = NULL;

    // adding missing arcs
    int nOrigArcs = nArcs;
    for ( int i=0 ; i<nOrigArcs ; ++i )
    {
        int u = tail[i];
        int v = head[i];
        const char *rname = arcName(str, v, u );

        if ( dict_int_get( arcDict, rname ) == -1 )
        {
            dict_int_set( arcDict, rname, nArcs );
            tail[nArcs] = v;
            head[nArcs] = u;
            cap[nArcs] = 0;
            ++nArcs;
        }
    }

    int *start, *nNeigh;
    struct MinCArc *arcs;
    ALLOCATE_VECTOR( arcs, struct MinCArc, nArcs );

    ALLOCATE_VECTOR( start, int, (n+1) );
    ALLOCATE_VECTOR_INI( nNeigh, int, n );

    // counting neighbors per node
    for ( int i=0 ; (i<nArcs) ; ++i )
        nNeigh[tail[i]]++;

    // setting up start
    start[0] = 0;
    for ( int i=1 ; (i<n+1) ; ++i )
        start[i] = start[i-1] + nNeigh[i-1];

    memset( nNeigh,  0, sizeof(int)*n );

    // storing arcs in positions
    // organized by tail
    for ( int i=0 ; (i<nArcs) ; ++i )
    {
        const int ctail = tail[i];
        const int chead = head[i];
        const int pos = start[ctail]+nNeigh[ctail];
        arcs[pos].v = chead;
        arcs[pos].cap = cap[i];
        const char *aname = arcName( str, ctail, chead );
        dict_int_set( arcDict, aname, pos );

        ++(nNeigh[ctail]);
        
    }

    free( nNeigh );

    // filling reverse arcs positions
    for ( int u=0 ; (u<n) ; ++u )
    {
        for ( int j=start[u] ; j<start[u+1] ; ++j )
        {
            int v = arcs[j].v;
            const char *rname = arcName( str, v, u );
            int rpos = dict_int_get( arcDict, rname );
            assert( rpos >= 0 && rpos < nArcs && rpos != j && rpos >= start[v] );
            assert( arcs[rpos].v == u );
            arcs[j].rpos = rpos;
        } // all v
    } // all us 


    minc->n = n;
    minc->orig = orig;
    minc->arcs = arcs;
    minc->start = start;
    minc->newIdx = ppnode;

    ALLOCATE_VECTOR_INI( minc->ivVisited, char, minc->n );
    ALLOCATE_VECTOR( minc->visited, int, minc->n );
    minc->nVisited = 0;

    ALLOCATE_VECTOR( minc->queue, int, n );
    ALLOCATE_VECTOR( minc->parent, int, n );

    minc->arcDict = arcDict;

    free( tail );

    ALLOCATE_VECTOR( minc->cutU, int, 2*n );
    minc->cutV = minc->cutU + n;

    minc->nCut = 0;
    
    /* info about original arcs */
    for ( int i=0 ; i<nArcs ; ++i )
        arcs[i].original = arcs[i].cap>0;

    return minc;
}

static void addVisited( MinCut *minc, int node )
{
#ifdef DEBUG
    assert( minc->ivVisited[node]==False );
#endif
    minc->ivVisited[node] = True;
    minc->visited[minc->nVisited++] = node;
}

static void clearVisited( MinCut *minc )
{
    for ( int i=0 ; (i<minc->nVisited) ; ++i )
        minc->ivVisited[minc->visited[i]] = False;
    minc->nVisited = 0;
#ifdef DEBUG
    for ( int i=0 ; (i<minc->n) ; ++i )
    {
        assert( minc->ivVisited[i]==False );
    }
#endif
}

static char bfs( MinCut *minc )
{
    minc->nVisited = 0;

    const int s = minc->s;
    const int t = minc->t;

    int *queue = minc->queue;
    int *parent = minc->parent;
    const int *start = minc->start;
    char *ivVisited = minc->ivVisited;
    const struct MinCArc *arcs = minc->arcs;


    queue[0] = s;
    int nQueue = 1;
    addVisited( minc, s );
    parent[s] = -1;

    while ( nQueue>0 )
    {
        int u = queue[--nQueue];

        // exploring neighbors of u
        for ( int p=start[u] ; p<start[u+1] ; ++p )
        {
            int v = arcs[p].v;

            if (ivVisited[v]==False && arcs[p].cap>0 )
            {
                queue[nQueue++] = v;
                parent[v] = u;
                addVisited( minc, v );
            }
        }
    }

    char reachedT = ivVisited[t];

    // clear visited
    clearVisited( minc );

    return reachedT;
}

static void dfs( MinCut *minc, int s )
{
    const int *start = minc->start;
    char *ivVisited = minc->ivVisited;

    addVisited( minc, s );
    
    // checking neighbors
    for ( int j=start[s] ; (j<start[s+1]) ; ++j )
    {
        const struct MinCArc *arc = minc->arcs+j;
        if ( ivVisited[arc->v]==False && arc->cap )
            dfs( minc, arc->v );
    }
}

int minc_optimize( MinCut *minc )
{
    const int s = minc->s;
    const int t = minc->t;
    const int *parent = minc->parent;
    struct MinCArc *arcs = minc->arcs;
    const int *start = minc->start;

    int totalFlow = 0;
    while ( bfs( minc ) )
    {
        int flow = INT_MAX;
       
        for ( int v=t; (v!=s) ; v=parent[v] )
        {
            int u = parent[v];

            int apos = arcPos( minc, u, v );
            flow = MIN( flow, arcs[apos].cap );
        } // checking path capacity
        assert( flow > 0 );
        
        totalFlow += flow;

        // updating residual capacities
        for ( int v=t; (v!=s) ; v=parent[v] )
        {
            int u = parent[v];
            int apos = arcPos( minc, u, v );
            arcs[apos].cap -= flow;

            int rpos = arcPos( minc, v, u );
            arcs[rpos].cap += flow;
        }

    } // while found a path

    if (totalFlow)
    {
        dfs( minc, minc->s );

        const char *ivVisited = minc->ivVisited;

        // checking arc cuts
        for ( int u=0 ; (u<minc->n) ; ++u )
        {
            if (!ivVisited[u])
                continue;

            for ( int pos=start[u] ; (pos<start[u+1]) ; ++pos )
            {
                const struct MinCArc *arc = arcs+pos;
                if (ivVisited[arc->v] || !arc->original)
                    continue;

                minc->cutU[minc->nCut] = minc->orig[u];
                minc->cutV[minc->nCut] =  minc->orig[arc->v];

                ++minc->nCut;
            } // destination side
        } // source side 
    }

    return totalFlow;
}

int minc_n_cut( MinCut *minc )
{
    return minc->nCut;
}


int minc_cut_arc_source( MinCut *minc, int i )
{
    return minc->cutU[i];
}


int minc_cut_arc_destination( MinCut *minc, int i )
{
    return minc->cutV[i];
}

int minc_n( MinCut *minc ) 
{
    return minc->n;
}

char minc_in_s(MinCut *minc, int i)
{
    return minc->ivVisited[minc->newIdx[i]];
}

void minc_free( MinCut **_minc )
{
    MinCut *minc = *_minc;

    free( minc->start );
    free( minc->arcs );
    free( minc->orig );
    free( minc->ivVisited );
    free( minc->visited );
    free( minc->queue );
    free( minc->parent );
    dict_int_free( &minc->arcDict );
    free( minc->cutU );
    if (minc->newIdx)
        free(minc->newIdx);
    free( minc );

    *_minc = NULL;
}