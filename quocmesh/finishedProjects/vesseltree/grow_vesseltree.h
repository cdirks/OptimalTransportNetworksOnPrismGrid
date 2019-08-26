/*
 * Copyright (c) 2001-2014 AG Rumpf, INS, Universitaet Bonn                      *
 *                                                                               *
 * The contents of this file are subject to the terms of the Common Development  *
 * and Distribution License Version 1.0 (the "License"); you may not use       *
 * this file except in compliance with the License. You may obtain a copy of     *
 * the License at http://www.opensource.org/licenses/CDDL-1.0                    *
 *                                                                               *
 * Software distributed under the License is distributed on an "AS IS" basis,  *
 * WITHOUT WARRANTY OF ANY KIND, either expressed or implied.                    *
 */

/* *******************************************************************
 * geomtree and vesseltree classes for generating trees according
 * to the Schreiner-Neumann model (CCO).
 * These classes provide functions used for this generating procedure
 * only that are not useful elsewhere.
 * ******************************************************************* */

#ifndef __GROW_VESSELTREE_H
#define __GROW_VESSELTREE_H

#include "vesseltree.h"
#include "parameterParser.h"
#include "epswriter.h"

template <typename real_type>
void compute_intersect_lammu ( const real_type x1_i, const real_type y1_i, const real_type x1_t, const real_type y1_t, const real_type x2_i, const real_type y2_i, const real_type x2_t, const real_type y2_t, real_type &lam, real_type &mu ) {
  const real_type
    A = x1_t - x1_i, // testtermnode.x    - testinitnode.x   ,
    B = x2_i - x2_t, // getXCO(arc2.init) - getXCO(arc2.term),
    C = y1_t - y1_i, // testtermnode.y    - testinitnode.y   ,
    D = y2_i - y2_t, // getYCO(arc2.init) - getYCO(arc2.term),
    E = x2_i - x1_i, // getXCO(arc2.init) - testinitnode.x   ,
    F = y2_i - y1_i, // getYCO(arc2.init) - testinitnode.y   ,

    G = A * D - B * C; // determinant of the matrix

  lam = ( 1.0 / G ) * ( D * E - B * F );
  mu  = ( 1.0 / G ) * ( A * F - C * E );
}


template <typename real_type>
class grow_geomtree: public geomtree<real_type> {
public:
  //  void grow_testtree(real_type x0, real_type y0, real_type x1, real_type y1);

  grow_geomtree() {
    params = new aol::ParameterParser ( "par/growscne.par" );
  }

  ~grow_geomtree() {
    delete ( params );
  }

  bool is_intersection_free (      const node<real_type> &testinitnode, const node<real_type> &testtermnode );
  bool is_intersection_free_scne ( const node<real_type> &testinitnode, const node<real_type> &testtermnode );
  bool is_intersection_free_scne ( const node<real_type> &testinitnode, const node<real_type> &testtermnode, const int );

  bool does_not_get_itself_too_close ( int i, real_type theta );
  bool does_not_get_othertree_too_close ( grow_geomtree<real_type> &other_tree, int i, real_type theta );

  void scne_set_levels_and_radii();
  void scne_set_levels_recurse ( int arc_index, char in_level );
  void scne_set_radii_recurse ( int arc_index, real_type &out_rad );

  inline bool is_in_domain ( const node<real_type> &testnode ) {
    return ( ( testnode.x >= 0 ) && ( testnode.x <= 1 ) && ( testnode.y >= 0 ) && ( testnode.y <= 1 ) );
  }

  inline bool is_in_domain ( const node<real_type> &/*testinitnode*/, const node<real_type> &testtermnode ) {
    // determines whether terminal node of test arc is inside domain of computation.
    // we can assume the initial node to be inside because it was set and verified a recursion before.
    bool res =
      ( testtermnode.x >= 0 ) && ( testtermnode.x <= 1 ) &&
      ( testtermnode.y >= 0 ) && ( testtermnode.y <= 1 ) ;
#ifdef VERBOSE_GRWOING
    if ( !res ) {
      cerr << "-" ;
    }
#endif
    return ( res );
  }

  inline bool is_acceptable ( const node<real_type> &testinitnode, const node<real_type> &testtermnode ) {
    return ( is_in_domain ( testinitnode, testtermnode ) && is_intersection_free ( testinitnode, testtermnode ) );
  }

  bool has_enough_space ( const node<real_type> &testnode, int level );
  bool scne_has_enough_space ( const node<real_type> &testnode, int attempt );

  inline real_type scne_local_volume ( const node<real_type> &node1, const real_type rad1,
                                       const node<real_type> &node2, const real_type rad2,
                                       const node<real_type> &node3, const real_type rad3,
                                       const node<real_type> &nodec ) {
    const real_type
      x1 = node1.x,   y1 = node1.y,   x  = nodec.x,   y  = nodec.y,   x2 = node2.x,   y2 = node2.y,   x3 = node3.x,   y3 = node3.y,
      R1 = sqrt ( pow ( x - x1, 2 ) + pow ( y - y1, 2 ) ),
      R2 = sqrt ( pow ( x - x2, 2 ) + pow ( y - y2, 2 ) ),
      R3 = sqrt ( pow ( x - x3, 2 ) + pow ( y - y3, 2 ) );

    /*     cerr << R1 << " " << rad1 << ", " << R2 << " " << rad2 << ", " << R3 << " " << rad3 << endl; */
    /*     cerr << rad1 * rad1 * R1 << ", " << rad2 * rad2 * R2 << ", " << rad3 * rad3 * R3 << endl << endl; */
    /*     cerr << aol::color::blue <<  rad1 * rad1 * R1 + rad2 * rad2 * R2 + rad3 * rad3 * R3 << aol::color::black << endl; */

    return ( aol::NumberTrait<real_type>::pi * ( rad1 * rad1 * R1 + rad2 * rad2 * R2 + rad3 * rad3 * R3 ) );
  }

  inline const real_type scne_local_volume_increment ( const node<real_type> &node1, const real_type rad1,
                                                       const node<real_type> &node2, const real_type rad2,
                                                       const node<real_type> &node3, const real_type rad3,
                                                       const node<real_type> &nodec ) {
    const real_type
      x1 = node1.x,   y1 = node1.y,   x  = nodec.x,   y  = nodec.y,   x2 = node2.x,   y2 = node2.y,   x3 = node3.x,   y3 = node3.y,
      R1 = sqrt ( pow ( x - x1, 2 ) + pow ( y - y1, 2 ) ),
      R2 = sqrt ( pow ( x - x2, 2 ) + pow ( y - y2, 2 ) ),
      R3 = sqrt ( pow ( x - x3, 2 ) + pow ( y - y3, 2 ) ),
      R0 = sqrt ( pow ( x1 - x3, 2 ) + pow ( y1 - y3, 2 ) );

    return ( aol::NumberTrait<real_type>::pi * ( rad1 * rad1 * R1 + rad2 * rad2 * R2 + rad3 * rad3 * R3 ) - aol::NumberTrait<real_type>::pi * rad1 * rad1 * R0 );
  }


  inline void scne_compute_neg_gradient ( const node<real_type> &node1, const real_type rad1,
                                          const node<real_type> &node2, const real_type rad2,
                                          const node<real_type> &node3, const real_type rad3,
                                          const node<real_type> &nodec,
                                          real_type &gx, real_type &gy ) {
    const real_type
      x1 = node1.x,   y1 = node1.y,   x  = nodec.x,   y  = nodec.y,   x2 = node2.x,   y2 = node2.y,   x3 = node3.x,   y3 = node3.y,
      R1 = sqrt ( pow ( x - x1, 2 ) + pow ( y - y1, 2 ) ),
      R2 = sqrt ( pow ( x - x2, 2 ) + pow ( y - y2, 2 ) ),
      R3 = sqrt ( pow ( x - x3, 2 ) + pow ( y - y3, 2 ) );

    gx = ( -1.0 ) * aol::NumberTrait<real_type>::pi * ( rad1 * rad1 / R1 * ( x - x1 ) +  rad2 * rad2 / R2 * ( x - x2 ) +  rad3 * rad3 / R3 * ( x - x3 ) );
    gy = ( -1.0 ) * aol::NumberTrait<real_type>::pi * ( rad1 * rad1 / R1 * ( y - y1 ) +  rad2 * rad2 / R2 * ( y - y2 ) +  rad3 * rad3 / R3 * ( y - y3 ) );
  }

  // from geomtree
  bool intersect ( const node<real_type> &testinitnode, const node<real_type> &testtermnode, const arc<real_type> &arc2 );
  bool is_self_intersection_free();
  bool is_mutually_intersection_free ( grow_geomtree<real_type> &othertree ); // why not const?
  inline real_type length ( int level ) {
    // length of vessel on level, ensure strictly positive!
    return ( static_cast<real_type> ( nrnd ( params->getDouble ( "length_minrnd" ), 1.0 ) ) * params->getDouble ( "length_max" )
             * pow ( ( 1.0 - ( static_cast<real_type> ( level ) / ( params->getInt ( "maxlevel" ) + 1.0 ) ) ), params->getDouble ( "length_exp" ) ) );
  }


public:
  real_type eta, theta;

  aol::ParameterParser* params;

};

template <typename real_type>
class grow_vesseltree: public vesseltree<real_type> {
public:

  void grow_randomtrees_scne ( int seed );
  void grow_randomtrees_add_node_scne ( grow_geomtree<real_type> &tree, grow_geomtree<real_type> &othertree );
  void scne_optimize_locally ( grow_geomtree<real_type> &tree, const node<real_type> &node1, const real_type rad1, const node<real_type> &node2, const real_type rad2, const node<real_type> &node3, const real_type rad3, node<real_type> &nodec );
  bool scne_is_still_feasible ( grow_geomtree<real_type> &tree, grow_geomtree<real_type> &othertree, node<real_type> &new_midpoint, node<real_type> &trial_term, int splitarc_index );

#if 0
  void enforce_relative_distances_to_arc ( grow_geomtree<real_type> &this_tree, grow_geomtree<real_type> &other_tree, int i, real_type theta );
  void enforce_min_length_of_arc ( grow_geomtree<real_type> &tree, int i, real_type eta );
#endif

#if 0
  inline bool has_enough_space ( const node<real_type> &testnode, int level ) {
    return ( gvtree.has_enough_space ( testnode, level ) && gatree.has_enough_space ( testnode, level ) );
  }
#endif

  inline bool is_acceptable ( const node<real_type> &testinitnode, const node<real_type> &testtermnode ) {
    return ( gvtree.is_acceptable ( testinitnode, testtermnode ) && gatree.is_acceptable ( testinitnode, testtermnode ) );
  }

  grow_geomtree<real_type> gvtree; // veins    - tree for growing!
  grow_geomtree<real_type> gatree; // arteries - tree for growing!


  void scne_set_levels_and_radii_on_both() {
    gatree.scne_set_levels_and_radii();
    gvtree.scne_set_levels_and_radii();
  }

  void epsout_grow_trees_radii ( char outfilename[] );

  void save_grow_trees ( int seed );

  void epsout_grow_trees_radii ( int seed ) {
    char outfilename[256];
    sprintf ( outfilename, "output/trees_%d_%d.eps", seed, static_cast<int> ( gatree.nodes.size() ) );
    epsout_grow_trees_radii ( outfilename );
  }
};

#endif
