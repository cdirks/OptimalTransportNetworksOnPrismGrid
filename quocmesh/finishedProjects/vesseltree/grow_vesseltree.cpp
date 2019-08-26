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
 * Implementation of grow_vesseltree
 * ******************************************************************* */

#include "grow_vesseltree.h"

#define partial_saving

template <typename real_type>
bool grow_geomtree<real_type>::is_intersection_free ( const node<real_type> &testinitnode,
                                                      const node<real_type> &testtermnode ) {
  bool result = true;
  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    if ( testinitnode.is_different_from ( this->nodes[ this->arcs[i].term ] ) ) {
      if ( intersect ( testinitnode, testtermnode, this->arcs[i] ) ) {
        result = false;
#ifdef VERBOSE_GROWING
        cerr << "." ;
#endif
        return ( result );
      }
      //      result = result && !( intersect(testinitnode, testtermnode, arcs[i]) );

      //       if (intersect(testinitnode, testtermnode, arcs[i] ) ) {
      //         //         cerr << "Intersection: " << testinitnode.x << " "  << testinitnode.y
      //         //              << " " << testtermnode.x << " " << testtermnode.y
      //         //              << " with arc " << i << " " << nodes[arcs[i].init].x << " "  << nodes[arcs[i].init].y
      //         //              << " " << nodes[arcs[i].term].x << " " << nodes[arcs[i].term].y << endl;
      //       }
    }
    ; // else arcs[i] is the parent arc of testarc, and those obviously may intersect
  }

  return ( result );
}


template <typename real_type>
bool grow_geomtree<real_type>::is_intersection_free_scne ( const node<real_type> &testinitnode,
                                                           const node<real_type> &testtermnode ) {
  // test whether the arc given by the two nodes intersects any arc in this tree
  // to be called on the tree which the arc NOT going to be added to.

  bool result = true;
  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    if ( intersect ( testinitnode, testtermnode, this->arcs[i] ) ) {
      result = false;
      return ( result );
    }
    //    result = result && !( intersect(testinitnode, testtermnode, arcs[i]) );
  }
  return ( result );
}

template <typename real_type>
bool grow_geomtree<real_type>::is_intersection_free_scne ( const node<real_type> &testinitnode,
                                                           const node<real_type> &testtermnode, const int exception_i ) {
  // test whether the arc given by the two nodes intersects any arc in this tree except for exception_i
  // to be called on the tree which the arc is going to be added to.
  bool result = true;
  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    if ( static_cast<int> ( i ) != exception_i ) {
      if ( intersect ( testinitnode, testtermnode, this->arcs[i] ) ) {
        result = false;
        return ( result );
      }
      //    result = result && !( intersect(testinitnode, testtermnode, arcs[i]) );
    } // else arcs[i] is the parent arc of the testarc given by the two nodes.
  }
  return ( result );
}


template <typename real_type>
void grow_geomtree<real_type>::scne_set_levels_and_radii() {
  real_type dummy;
  //  scne_set_levels_and_radii_recurse( (int)0, (char)0, dummy);
  scne_set_levels_recurse ( static_cast<int> ( 1 ), static_cast<char> ( 0 ) );
  this->nodes[0].set_level ( static_cast<char> ( 0 ) );
  scne_set_radii_recurse ( static_cast<int> ( 0 ), dummy );
}


template <typename real_type>
void grow_geomtree<real_type>::scne_set_radii_recurse ( int arc_index, real_type &out_rad ) {
  const real_type gamma = 3.0;
  int dau_1 = -1, dau_2 = -1; // can use stored data here.
  unsigned int i = 0;
  while ( ( dau_1 == -1 ) && i < this->arcs.size() ) {
    if ( this->arcs[i].init == this->arcs[arc_index].term ) dau_1 = i;
    i++;
  }
  while ( ( dau_2 == -1 ) && i < this->arcs.size() ) {
    if ( this->arcs[i].init == this->arcs[arc_index].term ) dau_2 = i;
    i++;
  }

  if ( ( dau_1 != -1 ) && ( dau_2 != -1 ) ) {
    // arc to bifurcation node
    real_type rd1, rd2;
    scne_set_radii_recurse ( dau_1, rd1 );
    scne_set_radii_recurse ( dau_2, rd2 );
    out_rad = static_cast<real_type> ( pow ( static_cast<real_type> ( pow ( rd1, gamma ) ) + static_cast<real_type> ( pow ( rd2, gamma ) ), static_cast<real_type> ( ( 1.0 / gamma ) ) ) );
    this->arcs[arc_index].set_rad ( out_rad );
  } else if ( ( dau_1 != -1 ) && ( dau_2 == -1 ) ) {
    // arc to intermediate node
    cout << "intermediate node" << endl;
    real_type rd1;
    scne_set_radii_recurse ( dau_1, rd1 );
    out_rad = rd1;
    this->arcs[arc_index].set_rad ( out_rad );
  } else if ( ( dau_1 == -1 ) && ( dau_2 == -1 ) ) {
    // arc to terminal node
    // out_rad = nrnd(0.5, 1.5) * (real_type) params->getDouble("LeafRadEst");
    out_rad = static_cast<real_type> ( params->getDouble ( "LeafRadEst" ) );
    this->arcs[arc_index].set_rad ( out_rad );
  } else {
    cerr << "grow_geomtree<real_type>::scne_set_radii_recurse: Houston, we have a problem ... " << endl;
  }
}

template <typename real_type>
void grow_geomtree<real_type>::scne_set_levels_recurse ( int arc_index, char in_level ) {
  this->nodes[ this->arcs[arc_index].term ].set_level ( in_level );

  int dau_1 = -1, dau_2 = -1; // can use stored data here.
  unsigned int i = 0;
  while ( ( dau_1 == -1 ) && i < this->arcs.size() ) {
    if ( this->arcs[i].init == this->arcs[arc_index].term ) dau_1 = i;
    i++;
  }
  while ( ( dau_2 == -1 ) && i < this->arcs.size() ) {
    if ( this->arcs[i].init == this->arcs[arc_index].term ) dau_2 = i;
    i++;
  }

  if ( ( dau_1 != -1 ) && ( dau_2 != -1 ) ) {
    // arc to bifurcation point
    scne_set_levels_recurse ( dau_1, in_level + 1 );
    scne_set_levels_recurse ( dau_2, in_level + 1 );
  } else if ( ( dau_1 != -1 ) && ( dau_2 == -1 ) ) {
    // arc to intermediate node
    scne_set_levels_recurse ( dau_1, in_level + 1 );
  } else if ( ( dau_1 == -1 ) && ( dau_2 == -1 ) ) {
    // arc to terminal node
    // do nothing
  } else {
    cerr << "grow_geomtree<real_type>::scne_set_levels_recurse: Houston, we have a problem ... " << endl;
  }
}


template <typename real_type>
bool grow_geomtree<real_type>::does_not_get_itself_too_close ( int i, real_type Theta ) {
  bool result = true;

  real_type
  x1_i = this->getXCO ( this->arcs[i].init ),
         y1_i = this->getYCO ( this->arcs[i].init ),
                x1_t = this->getXCO ( this->arcs[i].term ),
                       y1_t = this->getYCO ( this->arcs[i].term );

  unsigned int j;
  for ( j = 0; j < this->arcs.size(); j++ ) {
    if ( ( this->arcs[i].term != this->arcs[j].init ) && ( this->arcs[i].init != this->arcs[j].term ) && ( this->arcs[i].init != this->arcs[j].init ) ) {
      // not daughter, parent or sibling of arc i
      real_type
      x2_i = this->getXCO ( this->arcs[j].init ),
             y2_i = this->getYCO ( this->arcs[j].init ),
                    x2_t = this->getXCO ( this->arcs[j].term ),
                           y2_t = this->getYCO ( this->arcs[j].term );
      real_type lam = 0.0, mu = 0.0;

      compute_intersect_lammu ( x1_i, y1_i, x1_t, y1_t, x2_i, y2_i, x2_t, y2_t, lam, mu );
      if ( ( lam >= 0 ) && ( lam <= ( 1.0 + Theta ) ) && ( mu >= 0 ) && ( mu <= ( 1.0 + Theta ) ) ) {
        result = false;
        return ( result );
      }
    }

    for ( j = 0; j < this->arcs.size(); j++ ) {
      if ( ( this->arcs[i].init != this->arcs[j].term ) || ( this->arcs[i].init != this->arcs[j].init ) ) {
        //      if( (int) j != i ){
        // see if not too close to arcs: parent, sibling
        const real_type
          v1x = this->nodes[ this->arcs[j].term ].x - this->nodes[ this->arcs[j].init].x,
          v1y = this->nodes[ this->arcs[j].term ].y - this->nodes[ this->arcs[j].init].y,
          v2x = this->nodes[ this->arcs[i].term ].x - this->nodes[ this->arcs[j].init].x,
          v2y = this->nodes[ this->arcs[i].term ].x - this->nodes[ this->arcs[j].init].y,
          phi = acos ( v1x * v2x + v1y * v2y / ( sqrt ( v1x * v1x + v1y * v1y ) * sqrt ( v2x * v2x + v2y * v2y ) ) ),
          d   = sqrt ( v2x * v2x + v2y * v2y ) * abs ( sin ( phi ) );

        if ( d < params->getDouble ( "eta" ) ) {
          result = false;
          return ( result );
        }
      }
    }
  }

  return ( result );
}

template <typename real_type>
bool grow_geomtree<real_type>::does_not_get_othertree_too_close ( grow_geomtree<real_type> &other_tree, int i, real_type Theta ) {
  bool result = true;

  real_type
    x1_i = this->getXCO ( this->arcs[i].init ),
    y1_i = this->getYCO ( this->arcs[i].init ),
    x1_t = this->getXCO ( this->arcs[i].term ),
    y1_t = this->getYCO ( this->arcs[i].term );

  for ( unsigned int j = 0; j < other_tree.arcs.size(); j++ ) {
    real_type
      x2_i = other_tree.getXCO ( other_tree.arcs[j].init ),
      y2_i = other_tree.getYCO ( other_tree.arcs[j].init ),
      x2_t = other_tree.getXCO ( other_tree.arcs[j].term ),
      y2_t = other_tree.getYCO ( other_tree.arcs[j].term );
    real_type lam = 0.0, mu = 0.0;

    compute_intersect_lammu ( x1_i, y1_i, x1_t, y1_t, x2_i, y2_i, x2_t, y2_t, lam, mu );
    if ( ( lam >= 0 ) && ( lam <= ( 1.0 + Theta ) ) && ( mu >= 0 ) && ( mu <= ( 1.0 + Theta ) ) ) {
      result = false;
      return ( result );
    }
  }

  return ( result );
}


// from geomtree:

template <typename real_type>
bool grow_geomtree<real_type>::intersect ( const node<real_type> &testinitnode,
                                           const node<real_type> &testtermnode,
                                           const arc<real_type>  &arc2 ) {
  // this needs to be extended to some special cases (parallel...) where the determinant is zero!
  // solve the system:
  // testinitnode + lambda (testtermnode - testinitnode) = arc2.init + mu (arc2.term - arc2.init)
  // <=>
  // A  B . l     E
  // C  D   m  =  F

  const real_type
    x1_i = testinitnode.x,
    y1_i = testinitnode.y,
    x1_t = testtermnode.x,
    y1_t = testtermnode.y,
    x2_i = this->getXCO ( arc2.init ),
    y2_i = this->getYCO ( arc2.init ),
    x2_t = this->getXCO ( arc2.term ),
    y2_t = this->getYCO ( arc2.term );

  real_type lam = 0.0, mu = 0.0;

  compute_intersect_lammu ( x1_i, y1_i, x1_t, y1_t, x2_i, y2_i, x2_t, y2_t, lam, mu );

  return ( ( lam >= 0.0 )  && ( lam <= 1.0 )  && ( mu >= 0.0 ) && ( mu <= 1.0 ) );
}


template <typename real_type>
bool grow_geomtree<real_type>::is_self_intersection_free() {
  bool does_not_intersect_itself = true;

  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    for ( unsigned int j = i + 1; j < this->arcs.size(); j++ ) {
      if ( ( this->arcs[i].init != this->arcs[j].term ) && ( this->arcs[i].term != this->arcs[j].init )  && ( this->arcs[i].init != this->arcs[j].init ) ) {
        does_not_intersect_itself =
          ( does_not_intersect_itself && ! ( intersect ( this->nodes[ this->arcs[i].init ], this->nodes[ this->arcs[i].term], this->arcs[j] ) ) );
      }
    }
  }

  return ( does_not_intersect_itself );
}


template <typename real_type>
bool grow_geomtree<real_type>::is_mutually_intersection_free ( grow_geomtree<real_type> &othertree ) {
  bool does_not_intersect_other_tree = true;
  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    for ( unsigned int j = 0; j < othertree.arcs.size(); j++ ) {
      does_not_intersect_other_tree =
        ( does_not_intersect_other_tree && !othertree.intersect ( this->nodes[ this->arcs[i].init ], this->nodes[ this->arcs[i].term], othertree.arcs[j] ) );
    }
  }

  return ( does_not_intersect_other_tree );
}

// end from geomtree


template <typename real_type>
void grow_vesseltree<real_type>::grow_randomtrees_scne ( int seed ) {
  // initially set first arcs
  gatree.eta = gatree.params->getDouble ( "eta" );
  gvtree.eta = gvtree.params->getDouble ( "eta" );
  gatree.theta = gatree.params->getDouble ( "theta" );
  gvtree.theta = gvtree.params->getDouble ( "theta" );
#define first_nodes_det
#ifndef first_nodes_det
  // code deprecated!
  {
    node<real_type> first_a_node ( 0.0, 0.5 + nrnd ( -0.05, 0.05 ), 0 );
    gatree.nodes.push_back ( first_a_node );
    bool a_node_acceptable = false;
    while ( !a_node_acceptable ) {
      node<real_type> second_a_node ( nrnd ( 0.0, 1.0 ), nrnd ( 0.0, 1.0 ), 1 );
      a_node_acceptable = ( gatree.is_intersection_free_scne ( first_a_node, second_a_node ) &&
                            gvtree.is_intersection_free ( first_a_node, second_a_node ) );
      gatree.nodes.push_back ( second_a_node );
      // && has_enough_space(second_v_node, level) ); // need level in each step. so make tree "reasonable" in each step!
    }
    arc<real_type> first_a_arc ( 0, 1, 0.0625, 1.0 );
    gatree.arcs.push_back ( first_a_arc );

    node<real_type> first_v_node ( 1.0, 0.5 + nrnd ( -0.05, 0.05 ), 0 );
    gvtree.nodes.push_back ( first_v_node );

    bool v_node_acceptable = false;
    while ( !v_node_acceptable ) {
      node<real_type> second_v_node ( nrnd ( 0.0, 1.0 ), nrnd ( 0.0, 1.0 ), 1 );
      v_node_acceptable = ( gatree.is_intersection_free ( first_v_node, second_v_node ) &&
                            gvtree.is_intersection_free_scne ( first_v_node, second_v_node ) );
      gvtree.nodes.push_back ( second_v_node );
      // && has_enough_space(second_v_node, level) ); // need level in each step. so make tree "reasonable" in each step!
    }
    arc<real_type> first_v_arc ( 0, 1, 0.0625, 1.0 );
    gvtree.arcs.push_back ( first_v_arc );
  }
#else
  // initially set first arcs
  std::vector< node<real_type> > firstanodes ( 8 ), firstvnodes ( 8 );
  firstanodes[0].replaceby ( 0.00 , 0.50  , -1 );
  firstanodes[1].replaceby ( 0.05 , 0.50  , -1 );
  firstanodes[2].replaceby ( 0.25 , 0.20  , -1 );
  firstanodes[3].replaceby ( 0.25 , 0.70  , -1 );
  firstanodes[4].replaceby ( 0.60 , 0.10  , -2 );
  firstanodes[5].replaceby ( 0.60 , 0.30  , -2 );
  firstanodes[6].replaceby ( 0.70 , 0.55  , -2 );
  firstanodes[7].replaceby ( 0.60 , 0.80  , -2 );

  for ( unsigned int i = 0; i < 8; i++ ) {
    gatree.nodes.push_back ( firstanodes[i] );
  }

  firstvnodes[0].replaceby ( 1.00 , 0.50  , -1 );
  firstvnodes[1].replaceby ( 0.95 , 0.50  , -1 );
  firstvnodes[2].replaceby ( 0.75 , 0.30  , -1 );
  firstvnodes[3].replaceby ( 0.75 , 0.80  , -1 );
  firstvnodes[4].replaceby ( 0.40 , 0.20  , -2 );
  firstvnodes[5].replaceby ( 0.30 , 0.45  , -2 );
  firstvnodes[6].replaceby ( 0.40 , 0.70  , -2 );
  firstvnodes[7].replaceby ( 0.40 , 0.90  , -2 );

  for ( unsigned int i = 0; i < 8; i++ ) {
    gvtree.nodes.push_back ( firstvnodes[i] );
  }

  // the radius values given here will be ignored later on.
  std::vector< arc<real_type> > firstarcs ( 7 );
  firstarcs[0].replaceby ( 0, 1, 1.0 );
  firstarcs[1].replaceby ( 1, 2, 1.0 );
  firstarcs[2].replaceby ( 1, 3, 1.0 );
  firstarcs[3].replaceby ( 2, 4, 1.0 );
  firstarcs[4].replaceby ( 2, 5, 1.0 );
  firstarcs[5].replaceby ( 3, 6, 1.0 );
  firstarcs[6].replaceby ( 3, 7, 1.0 );

  for ( unsigned int i = 0; i < 7; i++ ) {
    gatree.arcs.push_back ( firstarcs[i] );
    gvtree.arcs.push_back ( firstarcs[i] );
  }
#endif

  gatree.scne_set_levels_and_radii();
  gvtree.scne_set_levels_and_radii();

  const int save_every = gatree.params->getInt ( "save_every" ),
                         orig_num = gatree.nodes.size(),
                                    orig_num_term = ( orig_num ) / 2; // sic, want integer division

  // add points according to ScNe Model
  for ( int i = orig_num_term; i < gatree.params->getInt ( "MaxScNeNodes" ); i++ ) {
#ifdef partial_saving
    if ( ( i == orig_num_term ) || ! ( i % save_every ) ) {
      save_grow_trees ( seed );
      epsout_grow_trees_radii ( seed );
    }
#endif

    aol::StopWatch timer;
    timer.start();
    cerr << "adding node " << i << ": ";
    grow_randomtrees_add_node_scne ( gatree, gvtree );
    cerr << " ... ";
    grow_randomtrees_add_node_scne ( gvtree, gatree );
    timer.stop();
    cerr << "done. Took " << timer.elapsedCpuTime() << " seconds." << endl;
  }

  // finally set node levels to reasonable values
  // and set segment radii to reasonable values
}


template <typename real_type>
void grow_vesseltree<real_type>::grow_randomtrees_add_node_scne ( grow_geomtree<real_type> &tree, grow_geomtree<real_type> &othertree ) {
  real_type optimal_value_so_far;
  int optimal_arc_index = -2;
  node<real_type> new_midpoint ( 2.0, 2.0, static_cast<unsigned char> ( 255 ) );  // for gcc 4.1??
  node<real_type> orig_midpoint ( 2.0, 2.0, static_cast<unsigned char> ( 255 ) );
  node<real_type> optimal_midpoint_so_far ( 2.0, 2.0, static_cast<unsigned char> ( 255 ) );
  node<real_type> trial_term ( 2.0, 2.0, static_cast<unsigned char> ( 255 ) );

  real_type radius_existing_half = 1e20;
  real_type radius_new_half = 1e20;
  real_type radius_old_half = 1e20;

  //   real_type init_point_dist = (real_type) tree.params->getDouble("init_point_dist");

  bool node_found = false; int count = 0; const int maxcount = gatree.params->getInt ( "MaxNodeIter" );

  while ( ( !node_found ) && ( count < maxcount ) ) {
    //     real_type point_dist = init_point_dist / sqrt( (real_type) tree.nodes.size() )  * pow( (real_type)0.97, count );
    trial_term.replaceby ( nrnd ( 0.0, 1.0 ), nrnd ( 0.0, 1.0 ), static_cast<unsigned char> ( 255 ) );
    //    if( tree.scne_has_enough_space(trial_term, point_dist) && othertree.scne_has_enough_space(trial_term, point_dist) ){
    if ( tree.scne_has_enough_space ( trial_term, count ) && othertree.scne_has_enough_space ( trial_term, count ) ) {
      optimal_arc_index = -1;
      optimal_value_so_far = 1e20;
      for ( unsigned int a = 0; a < tree.arcs.size(); a++ ) {
        orig_midpoint = tree.get_midpoint_of_arc ( static_cast<int> ( a ) );
        new_midpoint = orig_midpoint;
        radius_new_half = static_cast<real_type> ( gatree.params->getDouble ( "LeafRadEst" ) );
        radius_existing_half = tree.arcs[a].rad;
        radius_old_half = pow ( static_cast<real_type> ( ( pow ( radius_existing_half, 3 ) +  pow ( radius_new_half, 3 ) ) ), static_cast<real_type> ( ( 1.0 / 3.0 ) ) );

        // find arc with minimum add'l volume to midpoint (optimized location)
        // see whether new node can be connected to midpoint
        if ( tree.is_intersection_free_scne     ( orig_midpoint, trial_term, a )    &&
             othertree.is_intersection_free_scne ( orig_midpoint, trial_term   )   // &&
             //             scne_is_still_feasible(tree, othertree, orig_midpoint, trial_term, a)
           ) {

          scne_optimize_locally ( tree,
                                  tree.nodes[ tree.arcs[a].init ], radius_old_half,
                                  trial_term,                      radius_new_half,
                                  tree.nodes[ tree.arcs[a].term ], radius_existing_half,
                                  new_midpoint );

//            if( !scne_is_still_feasible(tree, othertree, new_midpoint, trial_term, a) ){
//              // infeasible after optimization: use original midpoint:
//              new_midpoint = orig_midpoint;
//            }

          real_type
          current_local_volume = tree.scne_local_volume_increment ( tree.nodes[ tree.arcs[a].init ], radius_old_half,
                                                                    trial_term,                      radius_new_half,
                                                                    tree.nodes[ tree.arcs[a].term ], radius_existing_half,
                                                                    new_midpoint );
          // if better, store
          if ( scne_is_still_feasible ( tree, othertree, new_midpoint, trial_term, a ) )
            if ( current_local_volume < optimal_value_so_far ) {
              optimal_value_so_far = current_local_volume;
              optimal_arc_index = a;
              optimal_midpoint_so_far = new_midpoint;
            }
        }
      }

      if ( ( optimal_arc_index != -1 ) ) {
        node_found = true;
      }
    }

    count++;
  }


  if ( optimal_arc_index != -1 ) {
    node<real_type> optimal_midpoint ( optimal_midpoint_so_far );
    int a = optimal_arc_index;

    if ( ! ( tree.is_intersection_free_scne     ( optimal_midpoint, trial_term, a ) &&
             othertree.is_intersection_free_scne ( optimal_midpoint, trial_term   ) ) ) {
      cerr << aol::color::red << "We're in trouble here" << aol::color::black << endl;
    }

    // add node as new midpoint.
    const int interm_node_index = tree.nodes.size();
    const int new_node_index = interm_node_index + 1;

    tree.nodes.push_back ( optimal_midpoint );
    tree.nodes.push_back ( trial_term );

    arc<real_type> remaining_half_arc ( interm_node_index, tree.arcs[a].term, radius_existing_half );
    arc<real_type> new_arc           ( interm_node_index, new_node_index,    radius_new_half );

    tree.arcs.push_back ( remaining_half_arc );
    tree.arcs.push_back ( new_arc );
    tree.arcs[a].replaceby ( tree.arcs[a].init, interm_node_index, radius_old_half );

    tree.scne_set_levels_and_radii();
  } else {
    cerr << "No node found in this iteration" << endl;
  }
}


template <typename real_type>
bool grow_vesseltree<real_type>::scne_is_still_feasible ( grow_geomtree<real_type> &tree, grow_geomtree<real_type> &othertree, node<real_type> &new_midpoint, node<real_type> &trial_term, int splitarc_index ) {
  // temporatily only.  // attention: this modifies the trees!
  // this is poorly designed. currently need to make sure tree is in the same state as before.
  // could improve this by setting up a separate "small" tree of the trial bifurcation.

  const int interm_node_index = tree.nodes.size();
  const int new_node_index = interm_node_index + 1;
  const int a = splitarc_index;
  real_type eta = static_cast<real_type> ( tree.params->getDouble ( "eta" ) );

  // temporarily adding nodes and arcs just for structural analysis, so we'll add useless values.
  tree.nodes.push_back ( new_midpoint );
  tree.nodes.push_back ( trial_term );

  arc<real_type> remaining_half_arc ( interm_node_index, tree.arcs[a].term, 1e20 );
  arc<real_type> new_arc           ( interm_node_index, new_node_index, 1e20 );
  arc<real_type> backup_arc ( tree.arcs[a] );

  const int new_arc_1 = tree.arcs.size();
  const int new_arc_2 = new_arc_1 + 1;
  tree.arcs.push_back ( remaining_half_arc );
  tree.arcs.push_back ( new_arc );
  tree.arcs[a].replaceby ( tree.arcs[a].init, interm_node_index, 1e20 );

  bool
  all_nodes_in_domain = true,
                        arcs_length_ok = true,
                                         doesnt_get_itself_too_close = true,
                                                                       doesnt_get_othertree_too_close = true;

  const real_type maxl = gatree.params->getDouble ( "iota" ) / sqrt ( static_cast<double> ( tree.arcs.size() ) );
  // arcs not too short or too long
  arcs_length_ok &= ( tree.length_of_arc (         a ) > eta  &&
                      tree.length_of_arc ( new_arc_1 ) > eta  &&
                      tree.length_of_arc ( new_arc_2 ) > eta  &&
                      //   tree.length_of_arc(         a ) < maxl && // arc to bifurcation does not matter
                      //   tree.length_of_arc( new_arc_1 ) < maxl && // "old" arc from bifurcation does not matter
                      tree.length_of_arc ( new_arc_2 ) < maxl    );  // only "new" arc from bifurcation may not be too long.

  if ( arcs_length_ok ) {

    for ( unsigned int i = 0; i < tree.nodes.size(); i++ ) {
      all_nodes_in_domain &= ( tree.getXCO ( i ) >= 0. ) && ( tree.getXCO ( i ) <= 1. ) && ( tree.getYCO ( i ) >= 0. ) && ( tree.getYCO ( i ) <= 1. );
    }

    if ( all_nodes_in_domain ) {
      doesnt_get_itself_too_close = (    tree.does_not_get_itself_too_close ( a, tree.theta ) &&
                                         tree.does_not_get_itself_too_close ( new_arc_1, tree.theta ) &&
                                         tree.does_not_get_itself_too_close ( new_arc_2, tree.theta )                   );

      if ( doesnt_get_itself_too_close ) {
        doesnt_get_othertree_too_close = ( tree.does_not_get_othertree_too_close ( othertree, a, tree.theta ) &&
                                           tree.does_not_get_othertree_too_close ( othertree, new_arc_1, tree.theta ) &&
                                           tree.does_not_get_othertree_too_close ( othertree, new_arc_2, tree.theta )     );

      }
    }
  }

  // go back to old tree.
  tree.nodes.pop_back(); tree.nodes.pop_back();
  tree.arcs.pop_back() ; tree.arcs.pop_back() ;
  tree.arcs[a].replaceby ( backup_arc.init, backup_arc.term, backup_arc.rad );

  return ( doesnt_get_itself_too_close && doesnt_get_othertree_too_close && all_nodes_in_domain && arcs_length_ok );

}

template <typename real_type>
bool grow_geomtree<real_type>::scne_has_enough_space ( const node<real_type> &testnode, int attempt ) {
  bool result = true;
  for ( unsigned int i = 0; i < this->nodes.size(); i++ ) {
    // see if testnode is not too close to other nodes
    if ( ( sqr ( this->nodes[i].x - testnode.x ) + sqr ( this->nodes[i].y - testnode.y ) ) <
         ( params->getDouble ( "pointspacing" ) / sqrt ( static_cast<real_type> ( this->nodes.size() ) ) * pow ( static_cast<real_type> ( 0.97 ) , attempt ) ) ) {
      result = false;
      return ( result );
    }
    //     result = result && ( ( sqr( nodes[i].x - testnode.x ) + sqr( nodes[i].y - testnode.y ) ) >
    //                          ( params->getDouble("pointspacing") / sqrt( static_cast<real_type> ( nodes.size() ) ) * pow( (real_type)0.97, attempt ) );
  }
  for ( unsigned int i = 0; i < this->arcs.size(); i++ ) {
    // see if testnode is not too close to arcs // is this correct??
    const real_type
      v1x = this->nodes[ this->arcs[i].term ].x - this->nodes[ this->arcs[i].init].x,
      v1y = this->nodes[ this->arcs[i].term ].y - this->nodes[ this->arcs[i].init].y,
      v2x = testnode.x - this->nodes[ this->arcs[i].init].x,
      v2y = testnode.y - this->nodes[ this->arcs[i].init].y,
      phi = acos ( v1x * v2x + v1y * v2y / ( sqrt ( v1x * v1x + v1y * v1y ) * sqrt ( v2x * v2x + v2y * v2y ) ) ),
      d   = sqrt ( v2x * v2x + v2y * v2y ) * abs ( sin ( phi ) );

    if ( ( d < ( ( params->getDouble ( "linespacing" ) / sqrt ( static_cast<real_type> ( ( this->arcs.size() ) ) ) * pow ( static_cast<real_type> ( 0.97 ), attempt ) ) ) ) ||
         ( d < params->getDouble ( "eta" ) ) ) {
      result = false;
      return ( result );
    }
    // CHANGE 1
    // result = result && ( d > ( params->getDouble("linespacing") / sqrt( (real_type)arcs.size() ) ) * pow( (real_type)0.97, attempt ) );
    //     result = result && ( d > params->getDouble("eta") );
  }
  return ( result );
}


template <typename real_type>
void grow_vesseltree<real_type>::scne_optimize_locally ( grow_geomtree<real_type> &tree,
                                                         const node<real_type> &node1, const real_type rad1,
                                                         const node<real_type> &node2, const real_type rad2,
                                                         const node<real_type> &node3, const real_type rad3,
                                                         node<real_type> &nodec ) {
  const real_type sigma = 0.5,
                          beta = 0.9; // good choice??
  real_type desc_x, desc_y;

  real_type current_local_volume,
  local_volume_step_before = tree.scne_local_volume ( node1, rad1,
                                                      node2, rad2,
                                                      node3, rad3,
                                                      nodec );
  bool stop_optimizing = false; int count_o = 0;
  while ( !stop_optimizing && ( count_o < 100 ) ) {
    cerr << "~";
    tree.scne_compute_neg_gradient ( node1, rad1, node2, rad2, node3, rad3, nodec, desc_x, desc_y );
    const real_type  s = 0.1 / sqrt ( pow ( desc_x, 2 ) + pow ( desc_y, 2 ) );
    //    cerr << "s = " << s << endl;

    real_type orig_volume = tree.scne_local_volume ( node1, rad1, node2, rad2, node3, rad3, nodec );
    node<real_type> moved_node ( 2.0, 2.0, static_cast<char> ( 255 ) );

    bool Armijo_stop = false;
    int k = 0;
    while ( !Armijo_stop && ( k < 100 ) ) {
      // cerr << "Armijo iteration " << k << " :";
      // use Armijo rule to locally optimize position of new cental node to obtain minimum local volume of the tree
      real_type alpha = s * pow ( beta, k );
      moved_node.replaceby ( nodec.x + alpha * desc_x, nodec.y + alpha * desc_y, nodec.level );
      //      cerr << "change : " <<  alpha * desc_x << ", " << alpha * desc_y << endl;
      real_type new_volume = tree.scne_local_volume ( node1, rad1, node2, rad2, node3, rad3, moved_node );
      const real_type secant_slope = ( new_volume - orig_volume ) / alpha;
      const real_type tangnt_slope = ( -1.0 ) * ( pow ( desc_x, 2 ) + pow ( desc_y, 2 ) );
      //       cerr << orig_volume <<  aol::color::red << new_volume << aol::color::black <<  secant_slope / tangnt_slope << endl;
      Armijo_stop = ( ( secant_slope / tangnt_slope ) > sigma  );
      k++;
    }

    if ( k < 100 ) {
      nodec = moved_node;
    }
    current_local_volume = tree.scne_local_volume ( node1, rad1,
                                                    node2, rad2,
                                                    node3, rad3,
                                                    nodec );
    //     cerr.precision(15);
    //    cerr << aol::color::green << current_local_volume << " " << local_volume_step_before << aol::color::black ;
    stop_optimizing = ( ( current_local_volume / local_volume_step_before ) > 0.999 );
    local_volume_step_before = current_local_volume;
    count_o++;
  }

  // trying to enforce minimal arc length
  const real_type
  min_len = static_cast<real_type> ( tree.params->getDouble ( "eta" ) ),
            len1 = sqrt ( ( node1.x - nodec.x ) * ( node1.x - nodec.x ) + ( node1.y - nodec.y ) * ( node1.y - nodec.y ) ),
                   len2 = sqrt ( ( node2.x - nodec.x ) * ( node2.x - nodec.x ) + ( node2.y - nodec.y ) * ( node2.y - nodec.y ) ),
                          len3 = sqrt ( ( node3.x - nodec.x ) * ( node3.x - nodec.x ) + ( node3.y - nodec.y ) * ( node3.y - nodec.y ) );

  if ( ( len1 <  min_len ) && ( len2 >= min_len ) && ( len3 >= min_len ) ) {
    nodec.x += ( 1.0 - min_len / len1 ) * ( nodec.x - node1.x );
    nodec.y += ( 1.0 - min_len / len1 ) * ( nodec.y - node1.y );
  } else if ( ( len1 >= min_len ) && ( len2 <  min_len ) && ( len3 >= min_len ) ) {
    nodec.x += ( 1.0 - min_len / len2 ) * ( nodec.x - node2.x );
    nodec.y += ( 1.0 - min_len / len2 ) * ( nodec.y - node2.y );
  } else if ( ( len1 >= min_len ) && ( len2 >= min_len ) && ( len3 <  min_len ) ) {
    nodec.x += ( 1.0 - min_len / len3 ) * ( nodec.x - node3.x );
    nodec.y += ( 1.0 - min_len / len3 ) * ( nodec.y - node3.y );
  } else {
    // do nothing. this tree will not be accepted.
  }

  // const real_type
  //     len1n = sqrt( ( node1.x - nodec.x ) * ( node1.x - nodec.x ) + ( node1.y - nodec.y ) * ( node1.y - nodec.y ) ),
  //     len2n = sqrt( ( node2.x - nodec.x ) * ( node2.x - nodec.x ) + ( node2.y - nodec.y ) * ( node2.y - nodec.y ) ),
  //     len3n = sqrt( ( node3.x - nodec.x ) * ( node3.x - nodec.x ) + ( node3.y - nodec.y ) * ( node3.y - nodec.y ) );

  //   if( (len1n > min_len) || (len2n > min_len) || (len3n > min_len) ){
  //     // do nothing. this tree will not be accepted.
  //   }
}


template <typename real_type>
void grow_vesseltree<real_type>::save_grow_trees ( int seed ) {
  char aname[256], vname[256];
  sprintf ( aname, "output/a_%d_%d.tree", seed, static_cast<int> ( gatree.nodes.size() ) );
  sprintf ( vname, "output/v_%d_%d.tree", seed, static_cast<int> ( gvtree.nodes.size() ) );

  gatree.save_tree ( aname );
  gvtree.save_tree ( vname );

}

template <typename real_type>
void grow_vesseltree<real_type>::epsout_grow_trees_radii ( char outfilename[] ) {
  int x1, x2, y1, y2, rd;
  //  int maxth = gatree.params->getInt("LineThickness");
  int maxth = gatree.params->getInt ( "LineThickness" );

  cerr << "Plotting tree to eps file " << outfilename << " ... ";

  aol::EpsWriter saver ( outfilename );

  real_type maxrad = 0.0;
  for ( unsigned int i = 0; i < gatree.arcs.size(); i++ ) {
    if ( gatree.arcs[i].rad > maxrad )        maxrad = gatree.arcs[i].rad;
  }
  for ( unsigned int i = 0; i < gvtree.arcs.size(); i++ ) {
    if ( gvtree.arcs[i].rad > maxrad )        maxrad = gatree.arcs[i].rad;
  }

  real_type factor = ( 1.0 * maxth ) / ( maxrad * 1000000.0 );

  for ( unsigned int i = 0; i < gatree.arcs.size(); i++ ) {
    x1 = static_cast<int> ( floor ( 1000000.0 * gatree.nodes[ gatree.arcs[i].init ].x + 0.5 ) );
    x2 = static_cast<int> ( floor ( 1000000.0 * gatree.nodes[ gatree.arcs[i].term ].x + 0.5 ) );
    y1 = static_cast<int> ( floor ( 1000000.0 * gatree.nodes[ gatree.arcs[i].init ].y + 0.5 ) );
    y2 = static_cast<int> ( floor ( 1000000.0 * gatree.nodes[ gatree.arcs[i].term ].y + 0.5 ) );
    rd = static_cast<int> ( floor ( 1000000.0 * factor * gatree.arcs[i].rad + 0.5 ) );

    saver.writeLine ( x1, y1, x2, y2, rd, 4 ); // changed
  }

  for ( unsigned int i = 0; i < gvtree.arcs.size(); i++ ) {
    x1 = static_cast<int> ( floor ( 1000000.0 * gvtree.nodes[ gvtree.arcs[i].init ].x + 0.5 ) );
    x2 = static_cast<int> ( floor ( 1000000.0 * gvtree.nodes[ gvtree.arcs[i].term ].x + 0.5 ) );
    y1 = static_cast<int> ( floor ( 1000000.0 * gvtree.nodes[ gvtree.arcs[i].init ].y + 0.5 ) );
    y2 = static_cast<int> ( floor ( 1000000.0 * gvtree.nodes[ gvtree.arcs[i].term ].y + 0.5 ) );
    rd = static_cast<int> ( floor ( 1000000.0 * factor * gvtree.arcs[i].rad + 0.5 ) );

    saver.writeLine ( x1, y1, x2, y2, rd, 1 );
  }

  //  saver.close();
  cerr << endl;
}


#if 0
template <typename real_type>
void grow_geomtree<real_type>::grow_testtree ( real_type x0, real_type y0, real_type x1, real_type y1 ) {
  // this test tree is for testing whether determining the grid cells in which an arc lies works correctly

#define TESTNEU


#ifdef TESTNEU
  node<real_type>
  node0 ( x0, y0, 0 ),
  node1 ( x1, y1, 0 );

  nodes.push_back ( node0 );
  nodes.push_back ( node1 );

  arc<real_type> arc0 ( 0, 1, 0.001, 1.0 );
  arcs.push_back ( arc0 );
#endif

#ifdef TESTNN
  node<real_type>
  node0 ( 0.35, 0.55, 0 ),
  node1 ( 0.35, 0.45, 0 ),
  node2 ( 0.40, 0.40, 0 ),
  node3 ( 0.45, 0.35, 0 ),
  node4 ( 0.55, 0.35, 0 ),
  node5 ( 0.60, 0.40, 0 ),
  node6 ( 0.65, 0.45, 0 ),
  node7 ( 0.65, 0.55, 0 ),
  node8 ( 0.60, 0.60, 0 ),
  node9 ( 0.55, 0.65, 0 ),
  node10 ( 0.45, 0.65, 0 ),
  node11 ( 0.40, 0.60, 0 );

  nodes.push_back ( node0 );
  nodes.push_back ( node1 );
  nodes.push_back ( node2 );
  nodes.push_back ( node3 );
  nodes.push_back ( node4 );
  nodes.push_back ( node5 );
  nodes.push_back ( node6 );
  nodes.push_back ( node7 );
  nodes.push_back ( node8 );
  nodes.push_back ( node9 );
  nodes.push_back ( node10 );
  nodes.push_back ( node11 );
#endif

#ifdef TEST8
  arc<real_type> arc0 ( 0, 8, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 1, 8, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 2, 8, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 3, 8, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 4, 8, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 5, 8, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 6, 8, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 7, 8, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 8, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 8, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 8, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST7
  arc<real_type> arc0 ( 0, 7, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 1, 7, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 2, 7, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 3, 7, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 4, 7, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 5, 7, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 6, 7, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 8, 7, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 7, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 7, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 7, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST5
  arc<real_type> arc0 ( 0, 5, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 1, 5, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 2, 5, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 3, 5, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 4, 5, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 6, 5, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 7, 5, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 8, 5, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 5, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 5, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 5, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST4
  arc<real_type> arc0 ( 0, 2, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 1, 2, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 3, 2, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 4, 2, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 5, 2, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 6, 2, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 7, 2, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 8, 2, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 2, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 2, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 2, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST3
  arc<real_type> arc0 ( 0, 1, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 2, 1, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 3, 1, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 4, 1, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 5, 1, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 6, 1, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 7, 1, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 8, 1, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 1, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 1, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 1, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST2
  arc<real_type> arc0 ( 1, 0, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 2, 0, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 3, 0, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 4, 0, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 5, 0, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 6, 0, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 7, 0, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 8, 0, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 9, 0, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 10, 0, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 11, 0, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
#ifdef TEST1
  arc<real_type> arc0 ( 0, 1, 1.0, 1.0 );
  arcs.push_back ( arc0 );
  arc<real_type> arc1 ( 0, 2, 1.0, 1.0 );
  arcs.push_back ( arc1 );
  arc<real_type> arc2 ( 0, 3, 1.0, 1.0 );
  arcs.push_back ( arc2 );
  arc<real_type> arc3 ( 0, 4, 1.0, 1.0 );
  arcs.push_back ( arc3 );
  arc<real_type> arc4 ( 0, 5, 1.0, 1.0 );
  arcs.push_back ( arc4 );
  arc<real_type> arc5 ( 0, 6, 1.0, 1.0 );
  arcs.push_back ( arc5 );
  arc<real_type> arc6 ( 0, 7, 1.0, 1.0 );
  arcs.push_back ( arc6 );
  arc<real_type> arc7 ( 0, 8, 1.0, 1.0 );
  arcs.push_back ( arc7 );
  arc<real_type> arc8 ( 0, 9, 1.0, 1.0 );
  arcs.push_back ( arc8 );
  arc<real_type> arc9 ( 0, 10, 1.0, 1.0 );
  arcs.push_back ( arc9 );
  arc<real_type> arc10 ( 0, 11, 1.0, 1.0 );
  arcs.push_back ( arc10 );
#endif
}
#endif


template class grow_geomtree<float>;
// template class grow_geomtree<double>;

template class grow_vesseltree<float>;
// template class grow_vesseltree<double>;
