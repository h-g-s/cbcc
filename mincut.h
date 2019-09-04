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

#ifndef MINCUT
#define MINCUT

/**
 * @file mincut.h
 * @author Haroldo Gambini Santos
 * @date 15 Aug 2019
 *
 * A self-contained, C99 compatible minimum cut implementation 
 * 
 * @see https://en.wikipedia.org/wiki/Minimum_cut
 */

typedef struct _MinCut MinCut;

/** @brief creates a min cut solver
 * @param nArcs number of arcs
 * @param tail vector with arc sources
 * @param head vector with arc destinations
 * @param _cap vector with arc capacities
 * @param s source
 * @param t destination
 * @return minimum cut solver
 **/
MinCut *minc_create( int nArcs, const int tail[], const int head[], const int _cap[], int s, int t );


/** @brief solves the min cut problem
 * @param minc mincut solver object
 * @return capacity of of min cut (max flow)
 */
int minc_optimize( MinCut *minc );

/** @brief number of nodes in graph
 * @param minc mincut solver object
 * @return number of nodes in the graph
 **/
int minc_n( MinCut *minc );


/** @brief number of arcs in minimum cut
 * @param minc mincut solver object
 * @return number of arcs in minimum cut
 **/
int minc_n_cut( MinCut *minc );


/** @brief returns the i-th arc source in the cut
 * @param minc mincut solver object
 * @param i index of the arc
 * @return source of the i-th arc in the cut
 **/
int minc_cut_arc_source( MinCut *minc, int i );


/** @brief returns the i-th arc destination in the cut
 * @param minc mincut solver object
 * @param i index of the destination arc
 * @return source of the i-th arc in the cut
 **/
int minc_cut_arc_destination( MinCut *minc, int i );

/** @brief checks if a node is in the subset of the source node or not
 * @param minc mincut solver object
 * @param i node to be checked
 * @return 1 is node i is in subset S, 0 otherwise
 **/
char minc_in_s(MinCut *minc, int i);


/** @brief frees memory of mincut solver
 **/
void minc_free( MinCut **_minc );

#endif

