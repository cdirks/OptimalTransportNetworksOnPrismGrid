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
 * Geometric binary trees - implementation
 * ******************************************************************* */

#include "vesseltree.h"

template <typename real_type>
geomtree<real_type>::geomtree() {}


template <typename real_type>
geomtree<real_type>::~geomtree() {}

template <typename real_type>
void geomtree<real_type>::clear() {
  arcs.clear();
  nodes.clear();
}

template <typename real_type>
void geomtree<real_type>::read_tree ( char* filename ) {
  FILE * indat_stream;

  if ( ! ( indat_stream = fopen ( filename, "r" ) ) ) {
    cerr << "Error: could not open " << filename << " for reading tree." << endl;
  }

  int nof_nodes, nof_arcs;

  fscanf ( indat_stream, "%d\n", &nof_nodes );
  fscanf ( indat_stream, "%d\n", &nof_arcs );

  nodes.resize ( nof_nodes );
  arcs.resize ( nof_arcs );

  double xco, yco;
  for ( int i = 0; i < nof_nodes; i++ ) {
    fscanf ( indat_stream, "%lf %lf\n", &xco, &yco );
    nodes[i].replaceby ( static_cast<real_type> ( xco ), static_cast<real_type> ( yco ) );
  }

  int init, term;
  double rad, speed;

  for ( int i = 0; i < nof_arcs; i++ ) {
    fscanf ( indat_stream, "%d %d %lf %lf\n", &init, &term, &rad, &speed );
    arcs[i].replaceby ( init, term, static_cast<real_type> ( rad ), static_cast<real_type> ( speed )  );
  }

  fclose ( indat_stream );
}


template <typename real_type>
void geomtree<real_type>::set_par_dau_sib() {
  // this implementation is inefficient.
  for ( unsigned int arc = 0; arc < arcs.size(); arc++ ) {
    int parent = -1, dau_1 = -1, dau_2 = -1, sibling = -1;
    unsigned int i = 0;
    while ( ( dau_1 == -1 ) && i < arcs.size() ) {
      if ( arcs[i].init == arcs[arc].term ) dau_1 = i;
      i++;
    }
    // i not set to 0

    while ( ( dau_2 == -1 ) && i < arcs.size() ) {
      if ( arcs[i].init == arcs[arc].term ) dau_2 = i;
      i++;
    }

    i = 0;
    while ( ( parent == -1 ) && i < arcs.size() ) {
      if ( arcs[i].term == arcs[arc].init ) parent = i;
      i++;
    }

    i = 0;
    while ( ( sibling == -1 ) && i < arcs.size() ) {
      if ( ( arcs[i].init == arcs[arc].init ) && ( arcs[i].term != arcs[arc].term ) ) sibling = i;
      i++;
    }

    arcs[arc].parent = parent;
    arcs[arc].dau_1 = dau_1;
    arcs[arc].dau_2 = dau_2;
    arcs[arc].sibling = sibling;
  }
}


template <typename real_type>
void geomtree<real_type>::get_terminal_arc_indices ( aol::Vector<int> &term_arc_indices, aol::BitVector &is_terminal_arc ) {
  int s_n = nodes.size();
  int s_a = arcs.size();
  int n_term_arcs = 0;

  bool* is_terminal_node = new bool[ s_n ];
  for ( int i = 0; i < s_n; i++ ) is_terminal_node[i] = true;

  for ( int j = 0; j < s_a; j++ ) {
    is_terminal_node[ arcs[j].init ] = false;
  }

  is_terminal_arc.resize ( arcs.size() );
  is_terminal_arc.setAll( false );

  for ( int j = 0; j < s_a; j++ ) {
    if ( is_terminal_node[ arcs[j].term ] ) {
      n_term_arcs++;
      is_terminal_arc.set( j, true);
    }
  }

  term_arc_indices.resize ( n_term_arcs );
  int n = 0;
  for ( int j = 0; j < s_a; j++ ) {
    if ( is_terminal_arc[j] ) {
      term_arc_indices[n] = j;
      n++;
    }
  }

  delete[] ( is_terminal_node );
}


template <typename real_type>
void geomtree<real_type>::set_terminal_speeds ( aol::Vector<int> &term, aol::Vector<real_type> &s, int offset ) {
  for ( int j = 0; j < term.size(); j++ ) {
    int i = term[j] ;
    real_type v = length_of_arc ( i ) * s[ j+offset ] / arcs[i].get_cross_sect_area() ;
    arcs[i].set_speed ( v );
  }
}


template <typename real_type>
void geomtree<real_type>::set_speeds_from_term( ) {
  real_type dummy = set_speeds_from_term_recurse ( 0 );
  dummy *= 1.0; // to avoid "unused" warning.
}


template <typename real_type>
real_type geomtree<real_type>::set_speeds_from_term_recurse ( int arc_index ) {
  // use:  int dau_1 = arcs[arc_index].dau_1, dau_2 = arcs[arc_index].dau_2;
  int dau_1 = -1, dau_2 = -1; // can use stored data here.
  unsigned int i = 0;
  while ( ( dau_1 == -1 ) && i < arcs.size() ) {
    if ( arcs[i].init == arcs[arc_index].term ) dau_1 = i;
    i++;
  }
  while ( ( dau_2 == -1 ) && i < arcs.size() ) {
    if ( arcs[i].init == arcs[arc_index].term ) dau_2 = i;
    i++;
  }

  if ( ( dau_1 != -1 ) && ( dau_2 != -1 ) ) {
    // arc to bifurcation node
    real_type v_here = ( ( set_speeds_from_term_recurse ( dau_1 ) * arcs[dau_1].get_cross_sect_area() +
                           set_speeds_from_term_recurse ( dau_2 ) * arcs[dau_2].get_cross_sect_area()   )
                         / arcs[arc_index].get_cross_sect_area()                                       );
    arcs[arc_index].set_speed ( v_here );

    // set theta on dau_1, dau_2

    return ( v_here );

  } else if ( ( dau_1 != -1 ) && ( dau_2 == -1 ) ) {
    // arc to intermediate node
    real_type v_here = ( set_speeds_from_term_recurse ( dau_1 ) * arcs[dau_1].get_cross_sect_area()
                         / arcs[arc_index].get_cross_sect_area()                                 );
    arcs[arc_index].set_speed ( v_here );
    cerr << "untested case of intermediate node" << endl;

    // set theta on dau_1

    return ( v_here );

  } else if ( ( dau_1 == -1 ) && ( dau_2 == -1 ) ) {
    // arc to terminal node

    // no theta to be set here.

    return ( arcs[arc_index].get_speed() );

  } else {
    cerr << "geomtree<real_type>::set_speeds_from_term_recurse: Houston, we have a problem ... " << endl;
    return ( -1 );
  }

  return ( -2 );
}


template <typename real_type>
void geomtree<real_type>::remove_slow_arcs ( geomtree<real_type> &dest_tree, real_type threshold, int arc_index, int cp_init_index ) {
  cerr << arcs[arc_index].speed << endl;
  if ( arcs[arc_index].speed > threshold ) {
    cerr << "recurring on " << arc_index << " , cp_init_index = " << cp_init_index << endl;
    int dau_1 = arcs[arc_index].dau_1, dau_2 = arcs[arc_index].dau_2;
    if ( dau_1 != -1 ) {
      cerr << "dau_1 = " << dau_1 << endl;
      node<real_type> cp_dau1node ( nodes[ arcs[dau_1].term ] );
      int cp_d1n_index = dest_tree.nodes.size();
      cerr << "before push_back: " << dest_tree.nodes.size() << endl;
      dest_tree.nodes.push_back ( cp_dau1node );
      cerr << "after push_back: " << dest_tree.nodes.size() << endl;
      arc<real_type> cp_dau1arc ( arcs[dau_1] );
      cp_dau1arc.init = cp_init_index;
      cp_dau1arc.term = cp_d1n_index;
      cerr << "new arc: " << cp_init_index << " to " << cp_d1n_index << endl;
      dest_tree.arcs.push_back ( cp_dau1arc );
      remove_slow_arcs ( dest_tree, threshold, dau_1, cp_d1n_index );
    }
    if ( dau_2 != -1 ) {
      cerr << "dau_2 = " << dau_2 << endl;
      node<real_type> cp_dau2node ( nodes[ arcs[dau_2].term ] );
      int cp_d2n_index = dest_tree.nodes.size();
      dest_tree.nodes.push_back ( cp_dau2node );
      arc<real_type> cp_dau2arc ( arcs[dau_2] );
      cp_dau2arc.init = cp_init_index;
      cp_dau2arc.term = cp_d2n_index;
      dest_tree.arcs.push_back ( cp_dau2arc );
      remove_slow_arcs ( dest_tree, threshold, dau_2, cp_d2n_index );
    }
  } // else do nothing
}

template <typename real_type>
int geomtree<real_type>::get_max_level() {
  int max_level = 0;
  get_max_level_recurse ( 0, 0, max_level );
  return ( max_level );
}


template <typename real_type>
void geomtree<real_type>::get_max_level_recurse ( int cur_arc, int cur_level, int &max_level_so_far ) {
  int dau_1 = arcs[ cur_arc ].dau_1;
  int dau_2 = arcs[ cur_arc ].dau_2;
  if ( dau_1 != -1 ) {
    if ( max_level_so_far < cur_level ) max_level_so_far = cur_level;
    get_max_level_recurse ( dau_1, cur_level + 1, max_level_so_far );
  }

  if ( dau_2 != -1 ) {
    if ( max_level_so_far < cur_level ) max_level_so_far = cur_level;
    get_max_level_recurse ( dau_2, cur_level + 1, max_level_so_far );
  }
}


template <typename real_type>
void geomtree<real_type>::compute_rhs_fphi ( qc::GridDefinition &grid, qc::ScalarArray<real_type, qc::QC_2D> &rhs ) {
  // assumes segments_in_grid set correctly

  qc::GridDefinition::OldFullElementIterator elit;
  typename std::vector< segment_in_cell<real_type> >::iterator segit;

  for ( elit = grid.begin(); elit != grid.end(); ++elit ) {
    qc::Element &El = *elit;
    int ix = El.x();
    int jy = El.y();

    for ( segit = segments_in_grid.at_begin ( ix, jy ); segit != segments_in_grid.at_end ( ix, jy ); ++segit ) {
      compute_rhs_fphi_single_term ( ix, jy, segit->x_in, segit->y_in, segit->x_out, segit->y_out,
                                     arcs[segit->arcno].get_value() * arcs[segit->arcno].get_cross_sect_area(), grid, rhs );
    }
  }
}


template <typename real_type>
void geomtree<real_type>::compute_rhs_fphi_single_term ( int ix,
                                                         int jy,
                                                         real_type lxa,
                                                         real_type lya,
                                                         real_type lxb,
                                                         real_type lyb,
                                                         real_type value,
                                                         qc::GridDefinition &grid,
                                                         qc::ScalarArray<real_type, qc::QC_2D> &rhs ) {
  // This uses 3 point quadrature. Can be changed because quadrature is invisible to the outside.

  real_type
    h = grid.H(),
    val1 =     lxa  *     lya  + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * (    lya  +     lyb ) ) +     lxb *    lyb ,
    val2 = ( 1. - lxa ) *     lya  + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * (    lya  +     lyb ) ) + ( 1. - lxb ) *    lyb ,
    val3 =     lxa  * ( 1. - lya ) + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) +     lxb * ( 1. - lyb ),
    val4 = ( 1. - lxa ) * ( 1. - lya ) + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) + ( 1. - lxb ) * ( 1. - lyb );

  // cerr << ix << " " << jy  << " " << xa  << " " <<  ya  << " " <<  xb  << " " <<  yb
  //        << " " << lxa  << " " <<  lya  << " " <<  lxb  << " " <<  lyb  << endl;

  val1 *= value * h / 6.0 ;
  val2 *= value * h / 6.0 ;
  val3 *= value * h / 6.0 ;
  val4 *= value * h / 6.0 ;

  //   cerr << val1 << " " << val2 << " " << val3 << " " << val4 << " " << endl;

  rhs.add ( ix  , jy  , val1 );
  if ( ix + 1 < grid.getWidth() ) {
    rhs.add ( ix + 1, jy  , val2 );
  }
  if ( jy + 1 < grid.getWidth() ) {
    rhs.add ( ix  , jy + 1, val3 );
  }
  if ( ( ix + 1 < grid.getWidth() ) && ( jy + 1 < grid.getWidth() ) ) {
    rhs.add ( ix + 1, jy + 1, val4 );
  }
}


template <typename real_type>
void geomtree<real_type>::update_segments_in_grid ( const qc::GridDefinition &grid ) {
  segments_in_grid.resize ( grid.getWidth(), grid.getWidth() );
  const real_type h = static_cast<real_type> ( grid.H() );

  for ( unsigned int current_index = 0; current_index < arcs.size(); current_index++ ) {
    real_type xg_init, yg_init, xg_term, yg_term;
    int  ix_init, jy_init, ix_term, jy_term;

    {
      qc::Element c1, c2;
      aol::Vec3<real_type> veca ( nodes[ arcs[ current_index ].init ].x, nodes[ arcs[ current_index ].init ].y );
      aol::Vec3<real_type> vecb ( nodes[ arcs[ current_index ].term ].x, nodes[ arcs[ current_index ].term ].y );

      grid.getElementByPoint ( aol::Vec3<double>( veca ), c1 );
      grid.getElementByPoint ( aol::Vec3<double>( vecb ), c2 );

      ix_init =  c1.x(); jy_init =  c1.y(); ix_term =  c2.x(); jy_term =  c2.y();
      xg_init = veca[0]; yg_init = veca[1]; xg_term = vecb[0]; yg_term = vecb[1];
    }

    if ( jy_init == jy_term ) {
      if ( ix_init == ix_term ) {
        // same cell
        update_segments_in_grid_stay ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );

      } else if ( ix_init < ix_term ) {
        update_segments_in_grid_go_E ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );

      } else if ( ix_init > ix_term ) {
        update_segments_in_grid_go_W ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );

      } else {
        cerr << "Impossible - 1 -" << endl;
      }
    } else if ( jy_init < jy_term ) {
      if ( ix_init == ix_term ) {
        update_segments_in_grid_go_N ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );

      } else if ( ix_init < ix_term ) {
        // going NE
        int  jyz_old, jyz_new;
        real_type yzg_old, yzg_new;

        yzg_new = yg_init + ( yg_term - yg_init ) * ( ( ix_init + 1 ) * h - xg_init ) / ( xg_term - xg_init );
        jyz_new = static_cast<int> ( floor ( yzg_new / h ) );

        update_segments_in_grid_go_N ( current_index, h,
                                       ix_init, jy_init, ix_init, jyz_new,
                                       xg_init, yg_init, ( h * ( ix_init + 1.0 ) ), yzg_new );
        yzg_old = yzg_new; jyz_old = jyz_new;

        for ( int lix = ( ix_init + 1 ); lix < ix_term; lix++ ) {
          yzg_new = yg_init + ( yg_term - yg_init ) * ( ( lix + 1 ) * h - xg_init ) / ( xg_term - xg_init );
          jyz_new = static_cast<int> ( floor ( yzg_new / h ) );
          update_segments_in_grid_go_N ( current_index, h,
                                         lix, jyz_old, lix, jyz_new,
                                         ( h * lix ), yzg_old, ( h * ( lix + 1 ) ), yzg_new );
          yzg_old = yzg_new; jyz_old = jyz_new;
        }

        update_segments_in_grid_go_N ( current_index, h,
                                       ix_term, jyz_old, ix_term, jy_term,
                                       ( h * ix_term ), yzg_old, xg_term, yg_term );

      } else if ( ix_init > ix_term ) {
        // going NW
        int  jyz_old, jyz_new;
        real_type yzg_old, yzg_new;

        yzg_new = yg_init + ( yg_term - yg_init ) * ( ix_init * h - xg_init ) / ( xg_term - xg_init );
        jyz_new = static_cast<int> ( floor ( yzg_new / h ) );

        update_segments_in_grid_go_N ( current_index, h,
                                       ix_init, jy_init, ix_init, jyz_new,
                                       xg_init, yg_init, ( h * ix_init ), yzg_new );
        yzg_old = yzg_new; jyz_old = jyz_new;

        for ( int lix = ( ix_init - 1 ); lix > ix_term; lix-- ) {
          yzg_new = yg_init + ( yg_term - yg_init ) * ( lix * h - xg_init ) / ( xg_term - xg_init );
          jyz_new = static_cast<int> ( floor ( yzg_new / h ) );
          update_segments_in_grid_go_N ( current_index, h,
                                         lix, jyz_old, lix, jyz_new,
                                         ( h * ( lix + 1.0 ) ), yzg_old, ( h * lix ), yzg_new );
          yzg_old = yzg_new; jyz_old = jyz_new;
        }

        update_segments_in_grid_go_N ( current_index, h,
                                       ix_term, jyz_old, ix_term, jy_term,
                                       ( h * ( ix_term + 1.0 ) ), yzg_old, xg_term, yg_term );

      } else {
        cerr << "Impossible - 2 -" << endl;
      }
    } else if ( jy_init > jy_term ) {
      if ( ix_init == ix_term ) {
        update_segments_in_grid_go_S ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );

      } else if ( ix_init < ix_term ) {
        // going SE
        int  jyz_old, jyz_new;
        real_type yzg_old, yzg_new;

        yzg_new = yg_init + ( yg_term - yg_init ) * ( ( ix_init + 1 ) * h - xg_init ) / ( xg_term - xg_init );
        jyz_new = static_cast<int> ( floor ( yzg_new / h ) );

        update_segments_in_grid_go_S ( current_index, h,
                                       ix_init, jy_init, ix_init, jyz_new,
                                       xg_init, yg_init, ( h * ( ix_init + 1.0 ) ), yzg_new );
        yzg_old = yzg_new; jyz_old = jyz_new;

        for ( int lix = ( ix_init + 1 ); lix < ix_term; lix++ ) {
          yzg_new = yg_init + ( yg_term - yg_init ) * ( ( lix + 1 ) * h - xg_init ) / ( xg_term - xg_init );
          jyz_new = static_cast<int> ( floor ( yzg_new / h ) );
          update_segments_in_grid_go_S ( current_index, h,
                                         lix, jyz_old, lix, jyz_new,
                                         ( h * lix ), yzg_old, ( h * ( lix + 1 ) ), yzg_new );
          yzg_old = yzg_new; jyz_old = jyz_new;
        }

        update_segments_in_grid_go_S ( current_index, h,
                                       ix_term, jyz_old, ix_term, jy_term,
                                       ( h * ix_term ), yzg_old, xg_term, yg_term );

      } else if ( ix_init > ix_term ) {
        // going SW
        int  jyz_old, jyz_new;
        real_type yzg_old, yzg_new;

        yzg_new = yg_init + ( yg_term - yg_init ) * ( ix_init * h - xg_init ) / ( xg_term - xg_init );
        jyz_new = static_cast<int> ( floor ( yzg_new / h ) );

        update_segments_in_grid_go_S ( current_index, h,
                                       ix_init, jy_init, ix_init, jyz_new,
                                       xg_init, yg_init, ( h * ix_init ), yzg_new );
        yzg_old = yzg_new; jyz_old = jyz_new;

        for ( int lix = ( ix_init - 1 ); lix > ix_term; lix-- ) {
          yzg_new = yg_init + ( yg_term - yg_init ) * ( lix * h - xg_init ) / ( xg_term - xg_init );
          jyz_new = static_cast<int> ( floor ( yzg_new / h ) );
          update_segments_in_grid_go_S ( current_index, h,
                                         lix, jyz_old, lix, jyz_new,
                                         ( h * ( lix + 1.0 ) ), yzg_old, ( h * lix ), yzg_new );
          yzg_old = yzg_new; jyz_old = jyz_new;
        }

        update_segments_in_grid_go_S ( current_index, h,
                                       ix_term, jyz_old, ix_term, jy_term,
                                       ( h * ( ix_term + 1.0 ) ), yzg_old, xg_term, yg_term );

      } else {
        cerr << "Impossible - 3 -" << endl;
      }
    } else {
      cerr << "Impossible - 4 -" << endl;
    }
  }
}



template <typename real_type>
void geomtree<real_type>::update_segments_in_grid_go_E ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term ) {
  if ( ( ix_init < ix_term ) && ( jy_init == jy_term )  ) {

    real_type zg, zl;
    segment_in_cell<real_type> seg_new, seg_cur;

    if ( abs ( yg_init - yg_term ) < ( 0.01 * h ) ) {
      // go horizontally
      zl = yg_init / h - 1.0 * jy_init;
      // first cell:
      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), 1.0, zl );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      // intermediate cells:
      for ( int lix = ( ix_init + 1 ); lix < ix_term; lix++ ) {
        seg_cur.replaceby ( current_index, 0.0, zl, 1.0, zl );
        segments_in_grid.append ( lix, jy_init, seg_cur );
      }
      // last cell:
      seg_cur.replaceby ( current_index, 0.0, zl, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    } else {
      // need to compute y-coordinates
      zg = yg_init + ( yg_term - yg_init ) * ( ( ix_init + 1 ) * h - xg_init ) / ( xg_term - xg_init );
      zl = zg / h - 1.0 * jy_init;

      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), 1.0, zl );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      seg_new.replaceby ( current_index, 0.0, zl, -1.0, -1.0 );

      for ( int lix = ( ix_init + 1 ); lix < ix_term; lix++ ) {
        zg = yg_init + ( yg_term - yg_init ) * ( ( lix + 1 ) * h - xg_init ) / ( xg_term - xg_init );
        zl = zg / h - 1.0 * jy_init;

        seg_cur.replaceby ( current_index, seg_new.x_in, seg_new.y_in, 1.0, zl );
        segments_in_grid.append ( lix, jy_init, seg_cur );
        seg_new.replaceby ( current_index, 0.0, zl, -1.0, -1.0 );
      }

      seg_cur.replaceby ( current_index, 0.0, zl, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    }
  } else if ( ( ix_init == ix_term ) && ( jy_init == jy_term ) ) {
    update_segments_in_grid_stay ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );
  } else  {
    cerr << "go_E is not the correct function!" << endl;
  }
}


template <typename real_type>
void geomtree<real_type>::update_segments_in_grid_go_W ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term ) {
  if ( ( ix_init > ix_term ) && ( jy_init == jy_term ) ) {
    real_type zg, zl;
    segment_in_cell<real_type> seg_new, seg_cur;

    if ( abs ( yg_init - yg_term ) < ( 0.01 * h ) ) {
      // go horizontally
      zl = yg_init / h - 1.0 * jy_init;
      // first cell:
      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), 0.0, zl );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      // intermediate cells:

      for ( int lix = ( ix_init - 1 ); lix > ix_term; lix-- ) {
        seg_cur.replaceby ( current_index, 1.0, zl, 0.0, zl );
        segments_in_grid.append ( lix, jy_init, seg_cur );
      }
      // last cell:
      seg_cur.replaceby ( current_index, 1.0, zl, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    } else {
      // need to compute y-coordinates
      zg = yg_init + ( yg_term - yg_init ) * ( ix_init * h - xg_init ) / ( xg_term - xg_init );
      zl = zg / h - 1.0 * jy_init;

      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), 0.0, zl );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      seg_new.replaceby ( current_index, 1.0, zl, -1.0, -1.0 );

      for ( int lix = ( ix_init - 1 ); lix > ix_term; lix-- ) {
        zg = yg_init + ( yg_term - yg_init ) * ( lix * h - xg_init ) / ( xg_term - xg_init );
        zl = zg / h - 1.0 * jy_init;

        seg_cur.replaceby ( current_index, seg_new.x_in, seg_new.y_in, 0.0, zl );
        segments_in_grid.append ( lix, jy_init, seg_cur );
        seg_new.replaceby ( current_index, 1.0, zl, -1.0, -1.0 );
      }

      seg_cur.replaceby ( current_index, 1.0, zl, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    }
  } else if ( ( ix_init == ix_term ) && ( jy_init == jy_term ) ) {
    update_segments_in_grid_stay ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );
  } else  {
    cerr << "go_W is not the correct function!" << endl;
  }
}


template <typename real_type>
void geomtree<real_type>::update_segments_in_grid_go_N ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term ) {
  if ( ( ix_init == ix_term ) && ( jy_init < jy_term ) ) {
    real_type zg, zl;
    segment_in_cell<real_type> seg_new, seg_cur;

    if ( abs ( xg_init - xg_term ) < ( 0.01 * h ) ) {
      // go horizontally
      zl = xg_init / h - 1.0 * ix_init;
      // first cell:
      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), zl, 1.0 );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      // intermediate cells:
      for ( int ljy = ( jy_init + 1 ); ljy < jy_term; ljy++ ) {
        seg_cur.replaceby ( current_index, zl, 0.0, zl, 1.0 );
        segments_in_grid.append ( ix_init, ljy, seg_cur );
      }
      // last cell:
      seg_cur.replaceby ( current_index, zl, 0.0, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    } else {
      // need to compute x-coordinates
      zg = xg_init + ( xg_term - xg_init ) * ( ( jy_init + 1 ) * h - yg_init ) / ( yg_term - yg_init );
      zl = zg / h - 1.0 * ix_init;

      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), zl, 1.0 );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      seg_new.replaceby ( current_index, zl, 0.0, -1.0, -1.0 );

      for ( int ljy = ( jy_init + 1 ); ljy < jy_term; ljy++ ) {
        zg = xg_init + ( xg_term - xg_init ) * ( ( ljy + 1 ) * h - yg_init ) / ( yg_term - yg_init );
        zl = zg / h - 1.0 * ix_init;

        seg_cur.replaceby ( current_index, seg_new.x_in, seg_new.y_in, zl, 1.0 );
        segments_in_grid.append ( ix_init, ljy, seg_cur );
        seg_new.replaceby ( current_index, zl, 0.0, -1.0, -1.0 );
      }

      seg_cur.replaceby ( current_index, zl, 0.0, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    }
  } else if ( ( ix_init == ix_term ) && ( jy_init == jy_term ) ) {
    update_segments_in_grid_stay ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );
  } else  {
    cerr << "go_N is not the correct function!" << endl;
  }
}


template <typename real_type>
void geomtree<real_type>::update_segments_in_grid_go_S ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term ) {
  if ( ( ix_init == ix_term ) && ( jy_init > jy_term ) ) {
    real_type zg, zl;
    segment_in_cell<real_type> seg_new, seg_cur;

    if ( abs ( xg_init - xg_term ) < ( 0.01 * h ) ) {
      // go horizontally
      zl = xg_init / h - 1.0 * ix_init;
      // first cell:
      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), zl, 0.0 );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      // intermediate cells:
      for ( int ljy = ( jy_init - 1 ); ljy > jy_term; ljy-- ) {
        seg_cur.replaceby ( current_index, zl, 1.0, zl, 0.0 );
        segments_in_grid.append ( ix_init, ljy, seg_cur );
      }
      // last cell:
      seg_cur.replaceby ( current_index, zl, 1.0, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    } else {
      // need to compute x-coordinates
      zg = xg_init + ( xg_term - xg_init ) * ( jy_init * h - yg_init ) / ( yg_term - yg_init );
      zl = zg / h - 1.0 * ix_init;

      seg_cur.replaceby ( current_index, ( xg_init / h - 1.0 * ix_init ), ( yg_init / h - 1.0 * jy_init ), zl, 0.0 );
      segments_in_grid.append ( ix_init, jy_init, seg_cur );
      seg_new.replaceby ( current_index, zl, 1.0, -1.0, -1.0 );

      for ( int ljy = ( jy_init - 1 ); ljy > jy_term; ljy-- ) {
        zg = xg_init + ( xg_term - xg_init ) * ( ljy * h - yg_init ) / ( yg_term - yg_init );
        zl = zg / h - 1.0 * ix_init;

        seg_cur.replaceby ( current_index, seg_new.x_in, seg_new.y_in, zl, 0.0 );
        segments_in_grid.append ( ix_init, ljy, seg_cur );
        seg_new.replaceby ( current_index, zl, 1.0, -1.0, -1.0 );
      }

      seg_cur.replaceby ( current_index, zl, 1.0, ( xg_term / h - 1.0 * ix_term ), ( yg_term / h - 1.0 * jy_term ) );
      segments_in_grid.append ( ix_term, jy_term, seg_cur );
    }
  } else if ( ( ix_init == ix_term ) && ( jy_init == jy_term ) ) {
    update_segments_in_grid_stay ( current_index, h, ix_init, jy_init, ix_term, jy_term, xg_init, yg_init, xg_term, yg_term );
  } else  {
    cerr << "go_S is not the correct function!" << endl;
  }
}


template <typename real_type>
void geomtree<real_type>::update_segments_in_grid_stay ( int current_index, real_type h, int ix_init, int jy_init, int ix_term, int jy_term, real_type xg_init, real_type yg_init, real_type xg_term, real_type yg_term ) {
  // stay in same cell
  if ( ( ix_init == ix_term ) && ( jy_init == jy_term ) ) {
    segment_in_cell<real_type>
    seg ( current_index,
          ( xg_init / h - 1.0 * ix_init ),
          ( yg_init / h - 1.0 * jy_init ),
          ( xg_term / h - 1.0 * ix_term ),
          ( yg_term / h - 1.0 * jy_term ) );
    segments_in_grid.append ( ix_init, jy_init, seg );
  } else {
    cerr << "stay is not the correct function!" << endl;
  }
}


template <typename real_type>
void geomtree<real_type>::update_grid_in_segments ( const qc::GridDefinition &grid, bool is_atree ) {
  real_type
  lambda_in  = -1.0,
               lambda_out = -1.0;

  qc::GridDefinition::OldFullElementIterator elit;
  typename std::vector< segment_in_cell<real_type> >::iterator segit;

  for ( elit = grid.begin(); elit != grid.end(); ++elit ) {
    const qc::Element &El = *elit;
    const int
    ix = El.x(),
         jy = El.y();

    //    cerr << "working on element " << ix << " " << jy << ":" << endl;

    for ( segit = segments_in_grid.at_begin ( ix, jy ); segit != segments_in_grid.at_end ( ix, jy ); ++segit ) {

      //      cerr << " segment " << segit->arcno << endl;

      const real_type
        xl_in  = segit->x_in,
        yl_in  = segit->y_in,
        xl_out = segit->x_out,
        yl_out = segit->y_out,

        xg_in  = grid.H() * ( static_cast<real_type> ( ix ) + xl_in ),
        yg_in  = grid.H() * ( static_cast<real_type> ( jy ) + yl_in ),
        xg_out = grid.H() * ( static_cast<real_type> ( ix ) + xl_out ),
        yg_out = grid.H() * ( static_cast<real_type> ( jy ) + yl_out ),

        xg_init = nodes[ arcs[segit->arcno].init ].x,
        xg_term = nodes[ arcs[segit->arcno].term ].x;

      if ( abs ( xg_term - xg_init ) > 1e-6 ) {
        // nodes not almost same x coordinate
        lambda_in  = ( xg_in  - xg_init ) / ( xg_term - xg_init );
        lambda_out = ( xg_out - xg_init ) / ( xg_term - xg_init );
      } else {
        const real_type
        yg_init = nodes[ arcs[segit->arcno].init ].y,
                  yg_term = nodes[ arcs[segit->arcno].term ].y;

        if ( abs ( yg_term - yg_init ) > 1e-6 ) {
          lambda_in  = ( yg_in  - yg_init ) / ( yg_term - yg_init );
          lambda_out = ( yg_out - yg_init ) / ( yg_term - yg_init );

        } else {
          cerr << "Error: Terminal and initial node too close" << endl;
        }
      }

      if ( !is_atree ) {
        // parametrization is different way around for vtree
        real_type temp_in = lambda_in, temp_out = lambda_out;
        lambda_in  = 1.0 - temp_out;
        lambda_out = 1.0 - temp_in;
      }

      if ( lambda_in < 0.0 ) {
        cerr << "lambda_in < 0.0 detected. " << endl;
        lambda_in = 0.0;
      }
      if ( lambda_out > 1.0 ) {
        cerr << "lambda_out > 1.0 detected. " << endl;
        lambda_out = 1.0;
      }

      //      cerr << "lambda_in = " << lambda_in << " lambda_out = " << lambda_out << endl;

      segit->set_lambdas ( lambda_in, lambda_out );

    }
  }


  grid_under_segments.resize ( arcs.size() );

  for ( elit = grid.begin(); elit != grid.end(); ++elit ) {
    const qc::Element &El = *elit;
    const int ix = El.x(), jy = El.y();
    for ( segit = segments_in_grid.at_begin ( ix, jy ); segit != segments_in_grid.at_end ( ix, jy ); ++segit ) {
      //      segment_in_cell<real_type> &dummy = *segit;
      segment_in_cell<real_type> dummy;
      dummy.replaceby ( segit->arcno, segit->x_in, segit->y_in, segit->lambda_in, segit->x_out, segit->y_out, segit->lambda_out );
      grid_under_segments.append ( ix, jy, dummy );
    }
  }
}


template <typename real_type>
void geomtree<real_type>::save_tree ( char fname[] ) {
  // Save tree data to file.
  cerr << "Saving tree in " << fname << endl;
  ofstream output ( fname );
  output.precision ( 15 );

  if ( output.is_open() ) {
    output << nodes.size() << endl;
    output << arcs.size() << endl;

    typename std::vector< node<real_type> >::iterator nit;
    for ( nit = nodes_at_begin(); nit != nodes_at_end(); ++nit ) {
      output << nit->x << " " << nit->y << endl;
    }
    typename std::vector< arc<real_type> >::iterator ait;
    for ( ait = arcs_at_begin(); ait != arcs_at_end(); ++ait ) {
      output << ait->init << " " << ait->term << " " << ait->rad << " " << ait->speed << endl ;
    }
  }

  output.close();
}


template <typename real_type>
void geomtree<real_type>::set_all_speeds ( real_type value ) {
  typename std::vector< arc<real_type> >::iterator ait;
  for ( ait = arcs_at_begin(); ait != arcs_at_end(); ++ait ) {
    ait->speed = value;
  }
}


template <typename real_type>
// void geomtree<real_type>::plot_tree_speeds( aol::EpsWriter *saver, const real_type v_min, const real_type v_max ){
void geomtree<real_type>::plot_tree_speeds ( aol::EpsWriter &saver, const real_type v_min, const real_type v_max ) {
  const real_type lvmin = speedfunc ( v_min ), lvmax = speedfunc ( v_max );

  for ( unsigned int i = 0; i < arcs.size(); i++ ) {
    const real_type val = speedfct ( i ) ;

    const int
      x1 = static_cast<int> ( floor ( 1000000.0 * nodes[ arcs[i].init ].x + 0.5 ) ),
      x2 = static_cast<int> ( floor ( 1000000.0 * nodes[ arcs[i].term ].x + 0.5 ) ),
      y1 = static_cast<int> ( floor ( 1000000.0 * nodes[ arcs[i].init ].y + 0.5 ) ),
      y2 = static_cast<int> ( floor ( 1000000.0 * nodes[ arcs[i].term ].y + 0.5 ) );

    int co = static_cast<int> ( floor ( 255.0 * (  val - lvmin  ) / ( lvmax - lvmin ) ) );

    if ( co > 255 ) {
      cerr << "color value too big: " << co << endl;
      co = 255;
    }
    if ( co < 0 ) {
      cerr << "color value too small: " << co << endl;
      co = 0;
    }
    //    saver->writeLine(x1, y1, x2, y2, 2000, co);
    saver.writeLine ( x1, y1, x2, y2, 2000, co );
  }
}


template <typename real_type>
double geomtree<real_type>::adapt_speeds_up_recurse ( int current_arc ) {
  // do not call directly!
  // make sure par-dau-sib information is up to date!
  // make sure all non-terminal arcs have speed < -450 !

  if ( arcs[ current_arc ].speed > -450 ) {

    return ( arcs[ current_arc ].speed );

  } else {

    double
      ret = 1e20,
      flow_dau1 = 0.0,
      flow_dau2 = 0.0;

    int dau1 = arcs[current_arc].dau_1;
    int dau2 = arcs[current_arc].dau_2;

    if ( dau1 != -1 ) {
      flow_dau1 = adapt_speeds_up_recurse ( dau1 ) * arcs[dau1].get_cross_sect_area();
    }

    if ( dau2 != -1 ) {
      flow_dau2 = adapt_speeds_up_recurse ( dau2 ) * arcs[dau2].get_cross_sect_area();
    }

    ret = ( flow_dau1 + flow_dau2 ) / arcs[ current_arc ].get_cross_sect_area();
    arcs[ current_arc ].speed = ret;

    return ( ret );

  }
}

template class geomtree<float>;
template class geomtree<double>;


// other general stuff

void normalize_speed ( geomtree<double> &tree, double maxi ) {
  // this is a workaround to avoid rho > 1.
  double maximum = 0;
  for ( unsigned int i = 0; i < tree.arcs.size(); i++ )
    if ( tree.arcs[i].get_speed() > maximum ) maximum = tree.arcs[i].get_speed();

  for ( unsigned int i = 0; i < tree.arcs.size(); i++ )
    tree.arcs[i].set_speed ( tree.arcs[i].get_speed() * maxi / maximum );
}


void img_scaled_save_ctrscl ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double ctr, double scl, const aol::colorTrans colTrans ) {
  qc::ScalarArray<double, qc::QC_2D> tmpimg ( img, aol::DEEP_COPY );
  tmpimg.addToAll ( -ctr );

  if ( scl == 0.0 ) { // this means no explicit scaling. use full color width.
    const double imin = tmpimg.getMinValue(), imax = tmpimg.getMaxValue();
    scl = 128. / ( max ( imax, -imin ) );
#ifdef VERBOSE
    cerr << "saving image: min = " << imin << ", max = " << imax << ", scalefactor = " << scl << endl;
#endif
  } else {
    scl *= 128.0;
  }
  tmpimg *= scl;
  tmpimg.addToAll ( 128.0 );

  tmpimg.setQuietMode ( true );
  tmpimg.setOverflowHandling ( aol::CLIP, 0, 255 );
  tmpimg.save ( filename, qc::PGM_DOUBLE_BINARY );

  FILE * outdat;
  outdat = fopen ( filename_color, "w" );
  fprintf ( outdat, "P3\n%d %d\n255\n", tmpimg.getNumX(), tmpimg.getNumY() );
  aol::RGBColorMap<double> colorMap(0.0, 255.0, colTrans);
  for ( int j = 0; j < tmpimg.getNumX(); j++ ) {
    for ( int i = 0; i < tmpimg.getNumY(); i++ ) {
      int
      val  = static_cast<int> ( tmpimg.get ( i, j ) );
      if ( val <   0 ) { cerr << "."; val =   0; }
      if ( val > 255 ) { cerr << "."; val = 255; }
      const int
        rval = static_cast<int> ( colorMap.scalarToColor ( val )[0] ),
        gval = static_cast<int> ( colorMap.scalarToColor ( val )[1] ),
  bval = static_cast<int> ( colorMap.scalarToColor ( val )[2] );

      fprintf ( outdat, "%d %d %d\n", rval, gval, bval );
    }
  }
  cerr << endl;
  fclose ( outdat );
}


void img_scaled_save ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const aol::colorTrans colTrans ) {
  img_scaled_save_ctrscl ( img, filename, filename_color, 0.5 * ( img.getMinValue() + img.getMaxValue() ), 0.0, colTrans );
}


void img_scaled_save_scl ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double scl, const aol::colorTrans colTrans ) {
  img_scaled_save_ctrscl ( img, filename, filename_color, 0.5 * ( img.getMinValue() + img.getMaxValue() ), scl, colTrans );
}


void img_scaled_save_ctr ( const qc::ScalarArray<double, qc::QC_2D> &img, char* filename, char* filename_color, const double ctr, const aol::colorTrans colTrans ) {
  img_scaled_save_ctrscl ( img, filename, filename_color, ctr, 0.0, colTrans );
}


void project_ ( aol::Vector<double> &Vec ) {
  // projects Vec to the subspace (1, ..., 1) Vec = 0
  double s = Vec.sum();
  cerr << "project_: sum = " << s << endl;
  int n = Vec.size();
  for ( int i = 0; i < n; i++ ) {
    Vec.add ( i, - ( s / static_cast< double > ( n ) ) );
  }

  //  cerr << "                                                       sum = " << s << ", after projection: " << Vec.sum() << endl;

}
