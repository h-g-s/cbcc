/*
 * spaths.h
 * Shortest Path Algorithm - Dijkstra Algorithm
 * specialized for sparse graphs
 * Developed by Haroldo Gambini Santos
 * hsantos@ic.uff.br
 */

#ifndef SPATHS_H
#define SPATHS_H

#define SP_INFTY_DIST ((INT_MAX/2)-1)

#define NULL_NODE -1

typedef struct _ShortestPathsFinder ShortestPathsFinder;
typedef  ShortestPathsFinder * ShortestPathsFinderPtr;

/*
 * creates Shortest Path Finder
 */
ShortestPathsFinderPtr spf_create(  );

/* arcStart[i] inditates the position in vector
 * toNode and dist where arcs of node i start
 **/
void spf_update_graph( ShortestPathsFinder* spf, const int nodes, const int arcs, const int *arcStart, const int *toNode, const int *dist );

/* updates just one arc
 **/
void spf_update_arc( ShortestPathsFinder* spf, const int tail, const int head, const int cost );

/* returns the distance of one arc
 **/
int spf_get_arc( ShortestPathsFinder* spf, const int tail, const int head );

/* temporarily removes
 * an arc
 **/
void spf_temp_remove_arc( ShortestPathsFinder* spf, const int tail, const int head );

/* temporarily removes
 * an arc
 **/
void spf_restore_arc( ShortestPathsFinder* spf, const int tail, const int head );

/*
 * loads a gr file to memory, creating a new ShortestPathsFinder object
 */
ShortestPathsFinder *spf_load_gr( const char *fileName );

/*
 * queries number of nodes
 */
int spf_nodes( ShortestPathsFinder* spf );

/*
 * queries number of arcs
 */
int spf_arcs( ShortestPathsFinder* spf );

//int spf_query_arc();

/*
 * executes the shortest path finder
 * using the Dijkstra algorithm
 */
void spf_find( ShortestPathsFinder* spf, const int origin );

/*
 * solution query: returns distance to a node after executing spf_find
 */
int spf_get_dist( const ShortestPathsFinderPtr spf, const int node );

/*
 * solution query: returns previous node and allows one to build a path after executing spf_find
 */
int spf_get_previous( const ShortestPathsFinderPtr spf, const int node );

int *spf_previous( const ShortestPathsFinder *spf );

/* returns all previous nodes
 * which should be steped
 * to arrive at a given node (this node is not included)
 * returns how many nodes were filles in indexes
 */
int spf_get_path( const ShortestPathsFinder *spf, const int toNode, int indexes[] );

/* returns all previous nodes
 * which should be steped
 * to arrive at a given node (this node is not included)
 * returns how many nodes were filles in indexes
 */
int spf_get_path_fw( const ShortestPathsFinder *spf, const int fromNode, const int toNode, int indexes[] );

/* executes All-Pairs shortest path finder using the
 * Floyd Warshall algorithm
 */
void spf_fw_find( ShortestPathsFinder* spf );

/*
 * returns true if floyd warshall ran before
 */
int spf_fw_ran( ShortestPathsFinder* spf );


/* returns the shortest path distance
 * between i and j after executing floyd-warshall
 */
int spf_fw_get_dist( ShortestPathsFinder* spf, const int i, const int j );

/*
 * releases Shortest Path Finder object
 */
void spf_free( ShortestPathsFinderPtr *spf );

#endif /* ifndef SPATHS_H */
