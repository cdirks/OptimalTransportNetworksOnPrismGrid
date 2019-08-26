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
 * Implementation of the tissue_vessel_geometry class
 * ******************************************************************* */

#include "tissue_vessel_geometry.h"
// #include <set>

template<> const int tissue_vessel_geometry<double>::ACOL = 20;
template<> const int tissue_vessel_geometry<double>::VCOL = 240;
template<> const int tissue_vessel_geometry<double>::TOPVIEW_RD = 2000;

template <typename real_type>
void tissue_vessel_geometry<real_type>::discretize_vessels() {
  A_narcs = vessels.atree.arcs.size();
  V_narcs = vessels.vtree.arcs.size();
  A_lenv.resize ( A_narcs );
  A_lenv.setZero();
  V_lenv.resize ( V_narcs );
  V_lenv.setZero();

  const real_type vessel_h = 1.0 / ( static_cast<real_type> ( 1 << vessellevel ) + 1.0 ); // approximate grid spacing desired on vessels
  real_type smallest_h = 2.0;

#ifdef VERBOSE
  cerr << "Vector lengths on the atree segments are: " << endl;
#endif
  for ( unsigned int i = 0; i < A_narcs; i++ ) {
    double len = vessels.atree.length_of_arc ( i );
    A_lenv[i] = static_cast<int> ( 1.0 + floor ( len / vessel_h ) );
    if ( i == 0 ) { // root segment has additional node
      A_lenv[i]++;
    }

    if ( i > 0 ) {
      if ( A_lenv[i] < 4 ) {
        A_lenv[i] = 4;
        cerr << "Attention: A_lenv[" << i << "] was < 4" << endl;
      }
      if ( ( len / A_lenv[i] ) < smallest_h ) smallest_h = len / A_lenv[i];
    } else {
      if ( A_lenv[i] < 5 ) {
        A_lenv[i] = 5;
        cerr << "Attention: A_lenv[" << i << "] was < 5" << endl;
      }
      if ( ( len / ( A_lenv[i] - 1.0 ) ) < smallest_h ) smallest_h = len / ( A_lenv[i] - 1.0 );
    }

#ifdef VERBOSE
    cerr << A_lenv[i] << endl;
#endif
  }

  // may require different treatment, so do not use procedure for these two.
#ifdef VERBOSE
  cerr << "Vector lengths on the vtree segments are: " << endl;
#endif
  for ( unsigned int i = 0; i < V_narcs; i++ ) {
    double len = vessels.vtree.length_of_arc ( i );
    V_lenv[i] = static_cast<int> ( 1.0 + floor ( len / vessel_h ) );
    if ( i == 0 ) { // root segment has additional node
      V_lenv[i]++;
    }

    if ( i > 0 ) {
      if ( V_lenv[i] < 4 ) {
        V_lenv[i] = 4;
        cerr << "Attention: V_lenv[" << i << "] was < 4" << endl;
      }
      if ( ( len / V_lenv[i] ) < smallest_h ) smallest_h = len / V_lenv[i];
    } else {
      if ( V_lenv[i] < 5 ) {
        V_lenv[i] = 5;
        cerr << "Attention: V_lenv[" << i << "] was < 5" << endl;
      }
      if ( ( len / ( V_lenv[i] - 1.0 ) ) < smallest_h ) smallest_h = len / ( V_lenv[i] - 1.0 );
    }

#ifdef VERBOSE
    cerr << V_lenv[i] << endl;
#endif
  }

  smallest_vessel_h = smallest_h;

  atree_data.reallocate ( A_lenv );
  vtree_data.reallocate ( V_lenv );

  vessels.atree.get_terminal_arc_indices ( A_term, A_is_term );
  vessels.vtree.get_terminal_arc_indices ( V_term, V_is_term );
}


template <typename real_type>
real_type tissue_vessel_geometry<real_type>::evaluate_mvec_on_atree_at_lambda ( aol::MultiVector<real_type> &mvec, int seg, real_type lambda ) {
  // evaluate using piecewise linear interpolation
  real_type ret = 0.0;
  int gc_int = 0;
  real_type gc_frac = 0.0;

  if ( seg == 0 ) {
    // evaluate on root segment

    gc_int  = static_cast<int> ( floor ( lambda * ( A_lenv[seg] - 1.0 ) ) );
    gc_frac = lambda * ( A_lenv[seg] - 1.0 ) - gc_int ;

    //    cerr << "lambda = " << lambda << ", gc_int = " << gc_int << ", gc_frac = " << gc_frac << endl;

    if ( gc_int == A_lenv[seg] - 1 ) {
      // exactly at terminal node
      if ( gc_frac > 1e-8 ) {
        cerr << "tissue_vessel_geometry::evaluate_aseg_at_lambda: attention: too far from terminal node" << endl;
      }
      ret = mvec[seg][gc_int];
    } else {
      // inside segment
      ret = ( 1.0 - gc_frac ) * mvec[seg][gc_int] + gc_frac * mvec[seg][gc_int+1];
    }

  } else {
    // evaluate on nonroot segment

    gc_int  = static_cast<int> ( floor ( lambda * A_lenv[seg] ) ) - 1 ;
    gc_frac = lambda * ( A_lenv[seg] ) - ( gc_int + 1 ) ;

    //    cerr << "lambda = " << lambda << ", gc_int = " << gc_int << ", gc_frac = " << gc_frac << endl;

    if ( gc_int == -1 ) {
      // get value from parent segment
      ret = ( ( 1.0 - gc_frac ) * mvec[ vessels.atree.arcs[seg].parent ][ A_lenv[vessels.atree.arcs[seg].parent] - 1 ]
              + gc_frac * mvec[seg][gc_int+1] );
    } else if ( gc_int == A_lenv[seg] - 1 ) {
      // exactly at terminal node
      if ( gc_frac > 1e-8 ) {
        cerr << "tissue_vessel_geometry::evaluate_aseg_at_lambda: attention: too far from terminal node" << endl;
      }
      ret = mvec[seg][gc_int];
    } else {
      // inside segment
      ret = ( 1.0 - gc_frac ) * mvec[seg][gc_int] + gc_frac * mvec[seg][gc_int+1];
    }

  }

  return ( ret );
}


template <typename real_type>
real_type tissue_vessel_geometry<real_type>::evaluate_mvec_on_vtree_at_lambda ( aol::MultiVector<real_type> &mvec, int seg, real_type lambda ) {
  // evaluate using piecewise linear interpolation
  real_type ret = 0.0;
  int gc_int = 0;
  real_type gc_frac = 0.0;

  if ( seg == 0 ) {
    // evaluate on root segment

    gc_int  = static_cast<int> ( floor ( lambda * ( V_lenv[seg] - 1.0 ) ) );
    gc_frac = lambda * ( V_lenv[seg] - 1.0 ) - gc_int ;

    //    cerr << "lambda = " << lambda << ", gc_int = " << gc_int << ", gc_frac = " << gc_frac << endl;

    if ( gc_int == V_lenv[seg] - 1 ) {
      // exactly at terminal node
      if ( gc_frac > 1e-8 ) {
        cerr << "tissue_vessel_geometry::evaluate_vseg_at_lambda: attention: too far from terminal node" << endl;
      }
      ret = mvec[seg][gc_int];
    } else {
      // inside segment
      ret = ( 1.0 - gc_frac ) * mvec[seg][gc_int] + gc_frac * mvec[seg][gc_int+1];
    }

  } else {
    // evaluate on nonroot segment

    gc_int  = static_cast<int> ( floor ( lambda * V_lenv[seg] ) )  ;
    gc_frac = lambda * ( V_lenv[seg] ) -  gc_int ;

    //    cerr << "lambda = " << lambda << ", gc_int = " << gc_int << ", gc_frac = " << gc_frac << endl;

    if ( gc_int == V_lenv[seg] ) {
      // exactly at terminal node
      if ( gc_frac > 1e-8 ) {
        cerr << "tissue_vessel_geometry::evaluate_vseg_at_lambda: attention: too far from terminal node" << endl;
      }
      ret = mvec[ vessels.vtree.arcs[seg].parent ][ 0 ]; // value from parent segment
    } else if ( gc_int == V_lenv[seg] - 1 ) {
      // need value from parent segment
      ret = ( ( 1.0 - gc_frac ) * mvec[seg][gc_int] +
              gc_frac * mvec[ vessels.vtree.arcs[seg].parent ][ 0 ] );
    } else {
      // inside segment
      ret = ( 1.0 - gc_frac ) * mvec[seg][gc_int] + gc_frac * mvec[seg][gc_int+1];
    }

  }

  return ( ret );
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::add_mvec_phi_products_onetree ( aol::MultiVector<real_type> &mvec, qc::ScalarArray<real_type, qc::QC_2D> &rhs, bool on_atree ) {

  qc::GridDefinition::OldFullElementIterator elit;
  typename std::vector< segment_in_cell<real_type> >::iterator segit, it_begin, it_end;

  for ( elit = tissue.begin(); elit != tissue.end(); ++elit ) {
    qc::Element &El = *elit;
    const int
    ix = El.x(),
         jy = El.y();

    it_begin = ( on_atree ? vessels.atree.segments_in_grid.at_begin ( ix, jy ) : vessels.vtree.segments_in_grid.at_begin ( ix, jy ) );
    it_end   = ( on_atree ? vessels.atree.segments_in_grid.at_end ( ix, jy )   : vessels.vtree.segments_in_grid.at_end ( ix, jy )   );

    for ( segit = it_begin; segit != it_end; ++segit ) {
      // this uses 3pt Lobatto quadrature. should do this piecewise to get exact integration.

      // not sure whether this is correct.

      const int
      arcno  = segit->arcno;
      const real_type
        lx_in  = segit->x_in,
        ly_in  = segit->y_in,
        lx_out = segit->x_out,
        ly_out = segit->y_out,
        lx_ctr = 0.5 * ( lx_in + lx_out ),
        ly_ctr = 0.5 * ( ly_in + ly_out ),

        ll_in  = segit->lambda_in,
        ll_out = segit->lambda_out,
        ll_ctr = 0.5 * ( ll_in + ll_out ),

        sval_in  = ( on_atree ? evaluate_mvec_on_atree_at_lambda ( mvec, arcno, ll_in )
                     : evaluate_mvec_on_vtree_at_lambda ( mvec, arcno, ll_in ) ),
        sval_ctr = ( on_atree ? evaluate_mvec_on_atree_at_lambda ( mvec, arcno, ll_ctr )
                     : evaluate_mvec_on_vtree_at_lambda ( mvec, arcno, ll_ctr ) ),
        sval_out = ( on_atree ? evaluate_mvec_on_atree_at_lambda ( mvec, arcno, ll_out )
                     : evaluate_mvec_on_vtree_at_lambda ( mvec, arcno, ll_out ) ),

        val1 = ( tissue_h / 6.0 *
                 ( sval_in * lx_in * ly_in + 4.0 * sval_ctr * (        lx_ctr  *        ly_ctr  ) + sval_out * lx_out * ly_out ) ),
        val2 = ( tissue_h / 6.0 *
                 ( sval_in * lx_in * ly_in + 4.0 * sval_ctr * ( ( 1.0 - lx_ctr ) *        ly_ctr  ) + sval_out * lx_out * ly_out ) ),
        val3 = ( tissue_h / 6.0 *
                 ( sval_in * lx_in * ly_in + 4.0 * sval_ctr * (        lx_ctr  * ( 1.0 - ly_ctr ) ) + sval_out * lx_out * ly_out ) ),
        val4 = ( tissue_h / 6.0 *
                 ( sval_in * lx_in * ly_in + 4.0 * sval_ctr * ( ( 1.0 - lx_ctr ) * ( 1.0 - ly_ctr ) ) + sval_out * lx_out * ly_out ) );

      rhs.add ( ix  , jy  , val1 );

      if ( ix + 1 < tissue.getWidth() ) {
        rhs.add ( ix + 1, jy  , val2 );
      }
      if ( jy + 1 < tissue.getWidth() ) {
        rhs.add ( ix  , jy + 1, val3 );
      }
      if ( ( ix + 1 < tissue.getWidth() ) && ( jy + 1 < tissue.getWidth() ) ) {
        rhs.add ( ix + 1, jy + 1, val4 );
      }
    }
  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step_topview_recursive ( char* ofn, aol::colorTrans colTrans, const bool temperatures, const real_type min_speed, const real_type ctr, const real_type scl ) {
  cerr << "Saving recursively to " << ofn << " " ;
  aol::EpsWriter saver ( ofn, colTrans );

  aol::EpsWriter alpha_saver ( "output/alpha.eps" );
  // plots black lines for the trees to be converted to alpha mask for pnmcomp

  // treat root segment
  {
    const real_type
      gxi = vessels.atree.nodes[ vessels.atree.arcs[0].init ].x,
      gyi = vessels.atree.nodes[ vessels.atree.arcs[0].init ].y,
      gxt = vessels.atree.nodes[ vessels.atree.arcs[0].term ].x,
      gyt = vessels.atree.nodes[ vessels.atree.arcs[0].term ].y;

    alpha_saver.writeLine ( static_cast<int> ( 1000000 * gxi ), static_cast<int> ( 1000000 * gyi ), static_cast<int> ( 1000000 * gxt ), static_cast<int> ( 1000000 * gyt ), TOPVIEW_RD, 0 );

    aol::Vector<int> g1 ( 2 ), g2 ( 2 );
    const real_type h_segm = 1.0 / ( static_cast<real_type> ( A_lenv[0] - 1 ) );

    for ( int j = 1; j < A_lenv[0]; j++ ) {
      g1.setZero();   g2.setZero();
      g1[0] = static_cast<int> ( 1000000 * ( gxi + ( j - 1 ) * h_segm * ( gxt - gxi ) ) );
      g2[0] = static_cast<int> ( 1000000 * ( gxi +   j   * h_segm * ( gxt - gxi ) ) );
      g1[1] = static_cast<int> ( 1000000 * ( gyi + ( j - 1 ) * h_segm * ( gyt - gyi ) ) );
      g2[1] = static_cast<int> ( 1000000 * ( gyi +   j   * h_segm * ( gyt - gyi ) ) );

      // FIXME: use ctr, scl
      real_type v_avg = 0.5 * ( atree_data[0][j-1] + atree_data[0][j] );
      int cl = static_cast<int> ( floor ( 128.0 + scl * 128 * ( v_avg - ctr ) ) ) ;
      if ( cl < 0 ) cl = 0;
      if ( cl > 255 ) cl = 255;

      saver.writeLine ( g1[0], g1[1], g2[0], g2[1], TOPVIEW_RD, cl );

    }
  }

  // recursive call of two daughters
  if ( vessels.atree.arcs[0].dau_1 != -1 )
    save_advect_step_topview_atree_recurse ( saver, alpha_saver, vessels.atree.arcs[0].dau_1, temperatures, min_speed, ctr, scl );
  if ( vessels.atree.arcs[0].dau_2 != -1 )
    save_advect_step_topview_atree_recurse ( saver, alpha_saver, vessels.atree.arcs[0].dau_2, temperatures, min_speed, ctr, scl );

  cerr << " ";

  // vtree:
  // treat root segment
  {
    const real_type
      gxt = vessels.vtree.nodes[ vessels.vtree.arcs[0].init ].x,
      gyt = vessels.vtree.nodes[ vessels.vtree.arcs[0].init ].y,
      gxi = vessels.vtree.nodes[ vessels.vtree.arcs[0].term ].x,
      gyi = vessels.vtree.nodes[ vessels.vtree.arcs[0].term ].y; // opposite way around

    alpha_saver.writeLine ( static_cast<int> ( 1000000 * gxi ), static_cast<int> ( 1000000 * gyi ), static_cast<int> ( 1000000 * gxt ), static_cast<int> ( 1000000 * gyt ), TOPVIEW_RD, 0 );

    aol::Vector<int> g1 ( 2 ), g2 ( 2 );
    const real_type h_segm = 1.0 / ( static_cast<real_type> ( V_lenv[0] - 1 ) );

    for ( int j = 0; j < V_lenv[0] - 1 ; j++ ) {
      g1.setZero();   g2.setZero();
      g1[0] = static_cast<int> ( 1000000 * ( gxi +   j   * h_segm * ( gxt - gxi ) ) );
      g2[0] = static_cast<int> ( 1000000 * ( gxi + ( j + 1 ) * h_segm * ( gxt - gxi ) ) );
      g1[1] = static_cast<int> ( 1000000 * ( gyi +   j   * h_segm * ( gyt - gyi ) ) );
      g2[1] = static_cast<int> ( 1000000 * ( gyi + ( j + 1 ) * h_segm * ( gyt - gyi ) ) );

      const real_type v_avg = 0.5 * ( vtree_data[0][ j ] + vtree_data[0][j+1] );
      int cl = static_cast<int> ( floor ( 128.0 + scl * 128 * ( v_avg - ctr ) ) ) ;
      if ( cl < 0 ) cl = 0;
      if ( cl > 255 ) cl = 255;

      saver.writeLine ( g1[0], g1[1], g2[0], g2[1], TOPVIEW_RD, cl );

    }
  }

  // recursive call of two daughters
  if ( vessels.vtree.arcs[0].dau_1 != -1 )
    save_advect_step_topview_vtree_recurse ( saver, alpha_saver, vessels.vtree.arcs[0].dau_1, temperatures, min_speed, ctr, scl );
  if ( vessels.vtree.arcs[0].dau_2 != -1 )
    save_advect_step_topview_vtree_recurse ( saver, alpha_saver, vessels.vtree.arcs[0].dau_2, temperatures, min_speed, ctr, scl );

  cerr << endl;

  //  saver.close();
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step_topview_atree_recurse ( aol::EpsWriter &saver, aol::EpsWriter &alpha_saver, int cur_arc, const bool temperatures, const real_type min_speed, const real_type ctr, const real_type scl ) {
  // not to be called on root segment
  if ( vessels.atree.arcs[cur_arc].speed >= min_speed ) {
    const real_type
      gxi = vessels.atree.nodes[ vessels.atree.arcs[cur_arc].init ].x,
      gyi = vessels.atree.nodes[ vessels.atree.arcs[cur_arc].init ].y,
      gxt = vessels.atree.nodes[ vessels.atree.arcs[cur_arc].term ].x,
      gyt = vessels.atree.nodes[ vessels.atree.arcs[cur_arc].term ].y;

    alpha_saver.writeLine ( static_cast<int> ( 1000000 * gxi ), static_cast<int> ( 1000000 * gyi ), static_cast<int> ( 1000000 * gxt ), static_cast<int> ( 1000000 * gyt ), TOPVIEW_RD, 0 );

    aol::Vector<int> g1 ( 2 ), g2 ( 2 );
    const real_type h_segm = 1.0 / ( static_cast<real_type> ( A_lenv[cur_arc] ) );
    real_type v_avg;

    for ( int j = 0; j < A_lenv[cur_arc]; j++ ) {
      g1.setZero();   g2.setZero();
      g1[0] = static_cast<int> ( 1000000 * ( gxi +   j   * h_segm * ( gxt - gxi ) ) ); // grid point 0, vector entry 0 is at position 1
      g2[0] = static_cast<int> ( 1000000 * ( gxi + ( j + 1 ) * h_segm * ( gxt - gxi ) ) );
      g1[1] = static_cast<int> ( 1000000 * ( gyi +   j   * h_segm * ( gyt - gyi ) ) );
      g2[1] = static_cast<int> ( 1000000 * ( gyi + ( j + 1 ) * h_segm * ( gyt - gyi ) ) );

      if ( j > 0 ) {
        v_avg = 0.5 * ( atree_data[cur_arc][j-1] + atree_data[cur_arc][j] ); // if discrete u = temperature, this is okay.
      } else {
        v_avg = 0.5 * ( atree_data[vessels.atree.arcs[cur_arc].parent][ A_lenv[vessels.atree.arcs[cur_arc].parent] - 1 ]
                        + atree_data[cur_arc][j] );
      }

      if ( temperatures && A_is_term[cur_arc] ) {
        v_avg *= ( A_lenv[cur_arc] * 1.0 ) / ( A_lenv[cur_arc] - 0.5 - j ) ; // temperature at center.
      }

      // FIXME: use ctr, scl
      int cl = static_cast<int> ( floor ( 128.0 + scl * 128 * ( v_avg - ctr ) ) ) ;
      if ( cl < 0 ) cl = 0;
      if ( cl > 255 ) cl = 255;

      saver.writeLine ( g1[0], g1[1], g2[0], g2[1], TOPVIEW_RD, cl );
    }

    if ( vessels.atree.arcs[cur_arc].dau_1 != -1 )
      save_advect_step_topview_atree_recurse ( saver, alpha_saver, vessels.atree.arcs[cur_arc].dau_1, temperatures, min_speed, ctr, scl );
    if ( vessels.atree.arcs[cur_arc].dau_2 != -1 )
      save_advect_step_topview_atree_recurse ( saver, alpha_saver, vessels.atree.arcs[cur_arc].dau_2, temperatures, min_speed, ctr, scl );
  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step_topview_vtree_recurse ( aol::EpsWriter &saver, aol::EpsWriter &alpha_saver, int cur_arc, const bool temperatures, const real_type min_speed, const real_type ctr, const real_type scl ) {
  if ( vessels.vtree.arcs[cur_arc].speed >= min_speed ) {
    const real_type
      gxt = vessels.vtree.nodes[ vessels.vtree.arcs[cur_arc].init ].x,
      gyt = vessels.vtree.nodes[ vessels.vtree.arcs[cur_arc].init ].y,
      gxi = vessels.vtree.nodes[ vessels.vtree.arcs[cur_arc].term ].x,
      gyi = vessels.vtree.nodes[ vessels.vtree.arcs[cur_arc].term ].y; // opposite way around

    alpha_saver.writeLine ( static_cast<int> ( 1000000 * gxi ), static_cast<int> ( 1000000 * gyi ), static_cast<int> ( 1000000 * gxt ), static_cast<int> ( 1000000 * gyt ), TOPVIEW_RD, 0 );

    aol::Vector<int> g1 ( 2 ), g2 ( 2 );
    const real_type h_segm = 1.0 / ( static_cast<real_type> ( V_lenv[cur_arc] ) );
    real_type v_avg;

    for ( int j = 0; j < V_lenv[cur_arc] ; j++ ) {
      g1.setZero();   g2.setZero();
      g1[0] = static_cast<int> ( 1000000 * ( gxi +   j   * h_segm * ( gxt - gxi ) ) );
      g2[0] = static_cast<int> ( 1000000 * ( gxi + ( j + 1 ) * h_segm * ( gxt - gxi ) ) );
      g1[1] = static_cast<int> ( 1000000 * ( gyi +   j   * h_segm * ( gyt - gyi ) ) );
      g2[1] = static_cast<int> ( 1000000 * ( gyi + ( j + 1 ) * h_segm * ( gyt - gyi ) ) );

      if ( j < ( V_lenv[cur_arc] - 1 ) ) {
        v_avg = 0.5 * ( vtree_data[cur_arc][ j ] + vtree_data[cur_arc][j+1] );
      } else {
        v_avg = 0.5 * ( vtree_data[cur_arc][j]  + vtree_data[ vessels.vtree.arcs[cur_arc].parent ][ 0 ] );
      }

      if ( temperatures && V_is_term[cur_arc] ) {
        v_avg *= ( V_lenv[cur_arc] * 1.0 ) / ( j + 0.5 ) ; // temperature at center.
      }

      // FIXME: use ctr, scl
      int cl = static_cast<int> ( floor ( 128.0 + scl * 128 * ( v_avg - ctr ) ) ) ;
      if ( cl < 0 ) cl = 0;
      if ( cl > 255 ) cl = 255;

      saver.writeLine ( g1[0], g1[1], g2[0], g2[1], TOPVIEW_RD, cl );
    }

    if ( vessels.vtree.arcs[cur_arc].dau_1 != -1 )
      save_advect_step_topview_vtree_recurse ( saver, alpha_saver, vessels.vtree.arcs[cur_arc].dau_1, temperatures, min_speed, ctr, scl );
    if ( vessels.vtree.arcs[cur_arc].dau_2 != -1 )
      save_advect_step_topview_vtree_recurse ( saver, alpha_saver, vessels.vtree.arcs[cur_arc].dau_2, temperatures, min_speed, ctr, scl );

  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step ( char* ofn, aol::colorTrans colTrans, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperatures, real_type min_speed, bool do_save_atree, bool do_save_vtree ) {
  cerr << "Saving to " << ofn << " " ;
  aol::EpsWriter saver ( ofn, colTrans );

  if ( do_save_atree )
    save_advect_step_atree ( saver, zzoom, color_scale, rot, shift, temperatures, min_speed );
  cerr << " ";
  if ( do_save_vtree )
    save_advect_step_vtree ( saver, zzoom, color_scale, rot, shift, temperatures, min_speed );
  cerr << endl;

  //  saver.close();
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step_atree ( aol::EpsWriter &saver, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperatures, real_type min_speed ) {
  int rd = 1000;

  for ( unsigned int i = 0; i < A_narcs; i++ ) {
    if ( vessels.atree.arcs[i].speed >= min_speed ) {
      real_type gxi, gyi, gxt, gyt, v_avg;
      gxi = vessels.atree.nodes[ vessels.atree.arcs[i].init ].x;
      gyi = /*1.0 -*/ vessels.atree.nodes[ vessels.atree.arcs[i].init ].y;
      gxt = vessels.atree.nodes[ vessels.atree.arcs[i].term ].x;
      gyt = /*1.0 -*/ vessels.atree.nodes[ vessels.atree.arcs[i].term ].y;

      if ( i == 0 ) { // root segment
        aol::Vector<real_type>
        g1 ( 3 ), g2 ( 3 ), g3 ( 3 ), g4 ( 3 );
        const real_type h_segm = 1.0 / ( static_cast<real_type> ( A_lenv[i] - 1 ) );

        g1[0] = gxi; g1[1] = gyi; g1[2] = 0.0;
        g2[0] = gxt; g2[1] = gyt; g2[2] = 0.0;
        _sas_segment ( g1, g2, shift, rot, rd, ACOL, saver );

        for ( int j = 1; j < A_lenv[i]; j++ ) {
          g1.setZero();   g2.setZero();   g3.setZero();   g4.setZero();
          g1[0] = gxi + ( j - 1 ) * h_segm * ( gxt - gxi );
          g2[0] = gxi +   j   * h_segm * ( gxt - gxi );
          g1[1] = gyi + ( j - 1 ) * h_segm * ( gyt - gyi );
          g2[1] = gyi +   j   * h_segm * ( gyt - gyi );

          g3 = g1;
          g4 = g2;

          g1[2] = zzoom * atree_data[i][j-1];
          g2[2] = zzoom * atree_data[i][j];
          v_avg = 0.5 * ( atree_data[i][j-1] + atree_data[i][j] );

          _sas_trapezoid ( g1, g2, g3, g4, shift, rot, v_avg, color_scale, saver );
        }

      } else { // intermediate or terminal segment
        aol::Vector<real_type>
        g1 ( 3 ), g2 ( 3 ), g3 ( 3 ), g4 ( 3 );
        const real_type h_segm = 1.0 / ( static_cast<real_type> ( A_lenv[i] ) );

        g1[0] = gxi; g1[1] = gyi; g1[2] = 0.0;
        g2[0] = gxt; g2[1] = gyt; g2[2] = 0.0;
        _sas_segment ( g1, g2, shift, rot, rd, ACOL, saver );

        for ( int j = 0; j < A_lenv[i]; j++ ) {
          g1.setZero();   g2.setZero();   g3.setZero();   g4.setZero();
          g1[0] = gxi +   j   * h_segm * ( gxt - gxi ); // grid point 0, vector entry 0 is at position 1
          g2[0] = gxi + ( j + 1 ) * h_segm * ( gxt - gxi );
          g1[1] = gyi +   j   * h_segm * ( gyt - gyi );
          g2[1] = gyi + ( j + 1 ) * h_segm * ( gyt - gyi );

          g3 = g1;
          g4 = g2;

          g2[2] = zzoom * atree_data[i][j];

          if ( j > 0 ) {
            g1[2] = zzoom * atree_data[i][j-1];
            v_avg = 0.5 * ( atree_data[i][j-1] + atree_data[i][j] ); // if discrete u = temperature, this is okay.
          } else {
            g1[2] = zzoom * atree_data[vessels.atree.arcs[i].parent][ A_lenv[vessels.atree.arcs[i].parent] - 1 ];
            v_avg = 0.5 * ( atree_data[vessels.atree.arcs[i].parent][ A_lenv[vessels.atree.arcs[i].parent] - 1 ] + atree_data[i][j] );
          }

          if ( temperatures && A_is_term[i] ) {
            g2[2] *= ( A_lenv[i] * 1.0 ) / ( 1.0 * A_lenv[i] - j ) ;
            g1[2] *= ( j == ( A_lenv[i] - 1 ) ? 0.0 : ( A_lenv[i] * 1.0 ) / ( A_lenv[i] - 1.0 - j ) );
            v_avg *= ( A_lenv[i] * 1.0 ) / ( A_lenv[i] - 0.5 - j ) ; // temperature at center.
          }

          _sas_trapezoid ( g1, g2, g3, g4, shift, rot, v_avg, color_scale, saver );
        }

      }
    }
  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_advect_step_vtree ( aol::EpsWriter &saver, real_type zzoom, real_type color_scale, aol::FullMatrix<real_type> &rot, aol::Vector<real_type> &shift, bool temperatures, real_type min_speed ) {
  int rd = 1000;

  for ( unsigned int i = 0; i < V_narcs; i++ ) {
    if ( vessels.vtree.arcs[i].speed >= min_speed ) {
      real_type gxi, gyi, gxt, gyt, v_avg;
      gxt = vessels.vtree.nodes[ vessels.vtree.arcs[i].init ].x;
      gyt = /*1.0 -*/ vessels.vtree.nodes[ vessels.vtree.arcs[i].init ].y;
      gxi = vessels.vtree.nodes[ vessels.vtree.arcs[i].term ].x;
      gyi = /*1.0 -*/ vessels.vtree.nodes[ vessels.vtree.arcs[i].term ].y; // opposite way around


      if ( i == 0 ) { // root segment
        aol::Vector<real_type>
        g1 ( 3 ), g2 ( 3 ), g3 ( 3 ), g4 ( 3 );
        const real_type h_segm = 1.0 / ( static_cast<real_type> ( V_lenv[i] - 1 ) );

        g1[0] = gxi; g1[1] = gyi; g1[2] = 0.0;
        g2[0] = gxt; g2[1] = gyt; g2[2] = 0.0;
        _sas_segment ( g1, g2, shift, rot, rd, VCOL, saver );

        for ( int j = 0; j < V_lenv[i] - 1 ; j++ ) {
          g1.setZero();   g2.setZero();   g3.setZero();   g4.setZero();
          g1[0] = gxi +   j   * h_segm * ( gxt - gxi );
          g2[0] = gxi + ( j + 1 ) * h_segm * ( gxt - gxi );
          g1[1] = gyi +   j   * h_segm * ( gyt - gyi );
          g2[1] = gyi + ( j + 1 ) * h_segm * ( gyt - gyi );

          g3 = g1;
          g4 = g2;

          g1[2] = zzoom * vtree_data[i][ j ];
          g2[2] = zzoom * vtree_data[i][j+1];
          v_avg = 0.5 * ( vtree_data[i][ j ] + vtree_data[i][j+1] );

          _sas_trapezoid ( g1, g2, g3, g4, shift, rot, v_avg, color_scale, saver );
        }

      } else { // intermediate or terminal segment
        aol::Vector<real_type> g1 ( 3 ), g2 ( 3 ), g3 ( 3 ), g4 ( 3 );
        const real_type h_segm = 1.0 / ( static_cast<real_type> ( V_lenv[i] ) );

        g1[0] = gxi; g1[1] = gyi; g1[2] = 0.0;
        g2[0] = gxt; g2[1] = gyt; g2[2] = 0.0;
        _sas_segment ( g1, g2, shift, rot, rd, VCOL, saver );

        for ( int j = 0; j < V_lenv[i] ; j++ ) {
          g1.setZero();   g2.setZero();   g3.setZero();   g4.setZero();
          g1[0] = gxi +   j   * h_segm * ( gxt - gxi );
          g2[0] = gxi + ( j + 1 ) * h_segm * ( gxt - gxi );
          g1[1] = gyi +   j   * h_segm * ( gyt - gyi );
          g2[1] = gyi + ( j + 1 ) * h_segm * ( gyt - gyi );

          g3 = g1;
          g4 = g2;

          g1[2] = zzoom * vtree_data[i][ j ];

          if ( j < ( V_lenv[i] - 1 ) ) {
            g2[2] = zzoom * vtree_data[i][j+1];
            v_avg = 0.5 * ( vtree_data[i][ j ] + vtree_data[i][j+1] );
          } else {
            g2[2] = zzoom * vtree_data[ vessels.vtree.arcs[i].parent ][ 0 ];
            v_avg = 0.5 * ( vtree_data[i][j]  + vtree_data[ vessels.vtree.arcs[i].parent ][ 0 ] );
          }

          if ( temperatures && V_is_term[i] ) {
            g1[2] *= ( j == 0 ? 0.0 : ( V_lenv[i] * 1.0 ) / ( 1.0 * j ) );
            g2[2] *= ( V_lenv[i] * 1.0 ) / ( j + 1.0 );
            v_avg *= ( V_lenv[i] * 1.0 ) / ( j + 0.5 ); // temperature at center.
          }

          _sas_trapezoid ( g1, g2, g3, g4, shift, rot, v_avg, color_scale, saver );
        }
      }
    }
  }
}



template <typename real_type>
void tissue_vessel_geometry<real_type>::_sas_trapezoid ( aol::Vector<real_type> &g1, aol::Vector<real_type> &g2, aol::Vector<real_type> &g3, aol::Vector<real_type> &g4, aol::Vector<real_type> &shift, aol::FullMatrix<real_type> &rot, real_type v_avg, real_type color_scale, aol::EpsWriter &saver ) {
  // constness!
  // rotate
  aol::Vector<real_type> dummy ( 3 );
  dummy = g1;          dummy -= shift;          rot.apply ( dummy, g1 );          g1 += shift;
  dummy = g2;          dummy -= shift;          rot.apply ( dummy, g2 );          g2 += shift;
  dummy = g3;          dummy -= shift;          rot.apply ( dummy, g3 );          g3 += shift;
  dummy = g4;          dummy -= shift;          rot.apply ( dummy, g4 );          g4 += shift;

  // project to first two components + save
  {
    int
      x1 = static_cast<int> ( floor ( 1000000.0 * g1[0] + 0.5 ) ),
      x2 = static_cast<int> ( floor ( 1000000.0 * g2[0] + 0.5 ) ),
      x3 = static_cast<int> ( floor ( 1000000.0 * g3[0] + 0.5 ) ),
      x4 = static_cast<int> ( floor ( 1000000.0 * g4[0] + 0.5 ) ),
      y1 = static_cast<int> ( floor ( 1000000.0 * g1[1] + 0.5 ) ),
      y2 = static_cast<int> ( floor ( 1000000.0 * g2[1] + 0.5 ) ),
      y3 = static_cast<int> ( floor ( 1000000.0 * g3[1] + 0.5 ) ),
      y4 = static_cast<int> ( floor ( 1000000.0 * g4[1] + 0.5 ) ),
      cl = static_cast<int> ( floor ( color_scale * 254.0 *  v_avg ) ) + 1 ;

    if ( cl <   0 ) {
      cl   = 0;
      cerr << "+";
    }
    if ( cl > 255 ) {
      cl = 255;
      cerr << "-";
    }

    saver.writeFilledTrapezoid ( x3, y3, x4, y4, x2, y2, x1, y1, cl );
  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::_sas_segment ( aol::Vector<real_type> &g1, aol::Vector<real_type> &g2, aol::Vector<real_type> &shift, aol::FullMatrix<real_type> &rot, int rd, int acol, aol::EpsWriter &saver ) {
  aol::Vector<real_type> dummy ( 3 );
  dummy = g1;            dummy -= shift;      rot.apply ( dummy, g1 );     g1 += shift;
  dummy = g2;     dummy -= shift;      rot.apply ( dummy, g2 );     g2 += shift;
  int
    x1 = static_cast<int> ( floor ( 1000000.0 * g1[0] + 0.5 ) ),
    x2 = static_cast<int> ( floor ( 1000000.0 * g2[0] + 0.5 ) ),
    y1 = static_cast<int> ( floor ( 1000000.0 * g1[1] + 0.5 ) ),
    y2 = static_cast<int> ( floor ( 1000000.0 * g2[1] + 0.5 ) );
  saver.writeLine ( x1, y1, x2, y2, rd, acol );
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::add_atree_outflow_rhs ( aol::MultiVector<real_type> &rhs, aol::MultiVector<real_type> &mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type tau, const real_type /*C_tis*/, const real_type C_ves ) {

  // something might be wrong here. temperature is accumulated towards
  // the end of outflow segments. probably could improve by better
  // averageing.

  for ( int i = 0; i < A_term.size(); i++ ) {
    const int
    arc = A_term[i],
          n = rhs[arc].size();

    if ( arc == 0 ) cerr << "Setting outflow on root segment. This will not work correctly, regardless of output." << endl;

    aol::Vector<real_type> dwarf ( n ), mutdwarf ( n );

    const real_type
      l = vessels.atree.length_of_arc ( arc ),
      v = vessels.atree.arcs[arc].speed,
      h = l / n,
      rho = v * tau / h,
      A = vessels.atree.arcs[arc].get_cross_sect_area();

    aol::CompositeOp< aol::Vector<real_type> > M_scaled, Me_scaled;
    aol::DiagonalMatrix<real_type> diag_m ( n ), diag_me ( n );

    for ( int j = 0; j < n; ++j ) {
      diag_m.set ( j, j, ( j == ( n - 1 ) ? 0.0 : n / ( n - j - 1.0 ) ) );
      diag_me.set ( j, j, n / ( n - j - 1.0 + rho ) );
    }

    M_scaled.appendReference (  ( M_block.getReference ( arc, arc ) ) );
    M_scaled.appendReference (  diag_m );

    Me_scaled.appendReference ( ( Me_block.getReference ( arc, arc ) ) );
    Me_scaled.appendReference ( diag_me );

    M_scaled.applyAdd ( atree_data[arc], dwarf );
    Me_scaled.applyAdd ( atree_data[arc], dwarf );
    dwarf *= vessels.atree.arcs[arc].speed * tau / ( ( -2.0 ) * vessels.atree.length_of_arc ( arc ) );

    for ( int j = 0; j < n; j++ ) {
      const real_type
        lambda_j   = ( 1.0 + j ) / ( n * 1.0 ),
        x_j        = ( vessels.atree.nodes[ vessels.atree.arcs[arc].init ].x + lambda_j *
                       ( vessels.atree.nodes[ vessels.atree.arcs[arc].term ].x - vessels.atree.nodes[ vessels.atree.arcs[arc].init ].x ) ),
        y_j        = ( vessels.atree.nodes[ vessels.atree.arcs[arc].init ].y + lambda_j *
                       ( vessels.atree.nodes[ vessels.atree.arcs[arc].term ].y - vessels.atree.nodes[ vessels.atree.arcs[arc].init ].y ) ),
        v_ves      = atree_data[arc][j],
        t_ves      = ( j == ( n - 1 ) ? 0.0 : v_ves * ( n * 1.0 ) / ( n - 1.0 - j ) ), // attention: denominator -> infty
        v_tis      = tissue_data.interpolate_on01 ( x_j, y_j ),
        t_tis      = ( j == ( n - 1 ) ? 0.0 : v_tis ),
        t_diff     = t_ves - t_tis;

      mutdwarf[j] += ( j > ( n - 5 ) ? 0.0 : ( ( -1.0 ) * h * A * t_diff * v * tau / l ) );
      // workaround to avoid artificial heating at the end of outflow segments.
    }

#ifdef OLD_VERSION
    aol::Vector<real_type> dwarf ( n ), mutdwarf ( n );
    Me_block.getReference ( arc, arc ).applyAdd ( atree_data[arc], dwarf );
    M_block.getReference ( arc, arc ).applyAdd ( atree_data[arc], dwarf );

    const real_type
      v = vessels.atree.arcs[arc].speed,
      l = vessels.atree.length_of_arc ( arc ),
      numerator = ( -1.0 ) * ( n - 1 ) * v * tau / ( 2.0 * l ),
      rho = vessels.atree.arcs[arc].speed * tau / ( vessels.atree.length_of_arc ( arc ) / A_lenv[arc] );
    for ( int j = 0; j < n; j++ ) {
      dwarf[j]    *= numerator / ( n - 1.0 + rho - j );

      const real_type
        lambda_j   = ( 1.0 + j ) / ( n * 1.0 ),
        x_j        = ( vessels.atree.nodes[ vessels.atree.arcs[arc].init ].x + lambda_j *
                       ( vessels.atree.nodes[ vessels.atree.arcs[arc].term ].x - vessels.atree.nodes[ vessels.atree.arcs[arc].init ].x ) ),
        y_j        = ( vessels.atree.nodes[ vessels.atree.arcs[arc].init ].y + lambda_j *
                       ( vessels.atree.nodes[ vessels.atree.arcs[arc].term ].y - vessels.atree.nodes[ vessels.atree.arcs[arc].init ].y ) ),
        v_ves      = atree_data[arc][j],
        t_ves      = ( j == ( n - 1 ) ? 0.0 : v_ves * ( n * 1.0 ) / ( n - 0.5 - j ) ), // attention: denominator -> infty
        v_tis      = tissue_data.interpolate_on01 ( x_j, y_j ),
        t_tis      = ( j == ( n - 1 ) ? 0.0 : v_tis ),
        t_diff     = t_ves - t_tis,
        h_x_A      = ( vessels.atree.length_of_arc ( arc ) / A_lenv[arc] ) * vessels.atree.arcs[arc].get_cross_sect_area();

      mutdwarf[j] += ( j > ( n - 5 ) ? 0.0 : ( ( -1.0 ) * h_x_A * t_diff * v * tau / l ) );
      // workaround to avoid artificial heating at the end of outflow segments.
    }
#endif

    rhs[arc] += dwarf;

    mutdwarf *= C_ves;
    mutsrc[arc] += mutdwarf;
    // outflow ~ temperature difference  into energy balancing
  }
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::add_vtree_inflow_rhs ( aol::MultiVector<real_type> &rhs, aol::MultiVector<real_type> &/*mutsrc*/, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type tau, const real_type /*C_tis*/, const real_type /*C_ves*/ ) {
  for ( int i = 0; i < V_term.size(); i++ ) {
    const int
    arc = V_term[i],
          n = rhs[arc].size();

    if ( arc == 0 ) cerr << "Setting inflow on root segment. This will not work correctly, regardless of output." << endl;

    aol::Vector<real_type> dwarf ( n ), eff ( n );
    for ( int j = 0; j < n ; j++ ) {
      eff[j] = ( n + 1.0 ) * vessels.vtree.arcs[arc].speed / ( n * vessels.vtree.length_of_arc ( arc ) ) * temp_surr_vtree ( arc, j );
      // eff[j] = ( n + 1.0 ) * vessels.vtree.arcs[arc].speed / ( 1.0 * n ) * u_surr( arc, j );
    }

    M_block.getReference ( arc, arc ).applyAdd ( eff, dwarf );
    Me_block.getReference ( arc, arc ).applyAdd ( eff, dwarf );
    dwarf *= ( tau / 2.0 );

    rhs[arc] += dwarf;

    //  dwarf *= C_ves;
    //  mutsrc[arc] += dwarf;
    // inflow: not considered in  energy balance
  }

}

template <typename real_type>
real_type tissue_vessel_geometry<real_type>::temp_surr_vtree ( const int arc, const int gridpt ) {
  // this only makes sense for the vtree!
  const real_type
    x_ = ( vessels.vtree.nodes[ vessels.vtree.arcs[arc].term ].x + ( 1.0 * gridpt ) / ( 1.0 * V_lenv[arc] ) *
           ( vessels.vtree.nodes[ vessels.vtree.arcs[arc].init ].x - vessels.vtree.nodes[ vessels.vtree.arcs[arc].term ].x ) ),
    y_ = ( vessels.vtree.nodes[ vessels.vtree.arcs[arc].term ].y + ( 1.0 * gridpt ) / ( 1.0 * V_lenv[arc] ) *
           ( vessels.vtree.nodes[ vessels.vtree.arcs[arc].init ].y - vessels.vtree.nodes[ vessels.vtree.arcs[arc].term ].y ) );

  return ( tissue_data.interpolate_on01 ( x_ , y_ ) );
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::add_ves_tis_transfer_rhs_onetree ( aol::MultiVector<real_type> &brhsmv, aol::MultiVector<real_type> &mutsrc, const aol::BlockOp< real_type, aol::Matrix<real_type> > &M_block, const aol::BlockOp< real_type, aol::Matrix<real_type> > &Me_block, const real_type /*tau*/, const real_type C_tis, const real_type /*C_ves*/, const real_type perc, const bool on_atree ) {

  aol::MultiVector<real_type> f_mv ( ( on_atree ? A_lenv : V_lenv ) ), f_mut_mv ( ( on_atree ? A_lenv : V_lenv ) ), dwarfmv ( ( on_atree ? A_lenv : V_lenv ) );

  for ( int a = 0; a < brhsmv.numComponents(); ++a ) {
    if ( on_atree ) {
      if ( a == 0 ) {
        // root segment of atree

        for ( int i = 0; i < A_lenv[a]; ++i ) {

          const real_type
            lambda_i   = ( 1.0 * i ) / ( A_lenv[a] - 1.0 ),
            x_i        = ( vessels.atree.nodes[ vessels.atree.arcs[a].init ].x + lambda_i *
                           ( vessels.atree.nodes[ vessels.atree.arcs[a].term ].x - vessels.atree.nodes[ vessels.atree.arcs[a].init ].x ) ),
            y_i        = ( vessels.atree.nodes[ vessels.atree.arcs[a].init ].y + lambda_i *
                           ( vessels.atree.nodes[ vessels.atree.arcs[a].term ].y - vessels.atree.nodes[ vessels.atree.arcs[a].init ].y ) ),
            v_ves      = atree_data[a][i],
            t_ves      = v_ves,
            v_tis      = tissue_data.interpolate_on01 ( x_i, y_i ),
            t_tis      = v_tis,
            f_trans    = perc * ( t_tis - t_ves );

          f_mv[a][i]   += f_trans;
          f_mut_mv[a][i] += f_trans;

        }

      } else if ( !A_is_term[a] ) {
        // intermediate segment of atree

        for ( int i = 0; i < A_lenv[a] ; ++i ) {

          const real_type
            lambda_i   = ( 1.0 + i ) / ( A_lenv[a] * 1.0 ),
            x_i        = ( vessels.atree.nodes[ vessels.atree.arcs[a].init ].x + lambda_i *
                           ( vessels.atree.nodes[ vessels.atree.arcs[a].term ].x - vessels.atree.nodes[ vessels.atree.arcs[a].init ].x ) ),
            y_i        = ( vessels.atree.nodes[ vessels.atree.arcs[a].init ].y + lambda_i *
                           ( vessels.atree.nodes[ vessels.atree.arcs[a].term ].y - vessels.atree.nodes[ vessels.atree.arcs[a].init ].y ) ),
            v_ves      = atree_data[a][i],
            t_ves      = v_ves,
            v_tis      = tissue_data.interpolate_on01 ( x_i, y_i ),
            t_tis      = v_tis,
            f_trans    = perc * ( t_tis - t_ves );

          f_mv[a][i]   += f_trans;
          f_mut_mv[a][i] += f_trans;

        }

      } else {
        // terminal segment of atree: no heat diffusion
        // do nothing.
      }

    } else {
      if ( a == 0 ) {
        // root segment of vtree
        for ( int i = 0; i < V_lenv[a] - 1 ; ++i ) {
          const real_type
            lambda_i   = ( 1.0 * i ) / ( V_lenv[a] - 1.0 ),
            x_i        = ( vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].x + lambda_i *
                           ( vessels.vtree.nodes[ vessels.vtree.arcs[a].init ].x - vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].x ) ),
            y_i        = ( vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].y + lambda_i *
                           ( vessels.vtree.nodes[ vessels.vtree.arcs[a].init ].y - vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].y ) ),
            v_ves      = vtree_data[a][i] ,
            t_ves      = v_ves,
            v_tis      = tissue_data.interpolate_on01 ( x_i, y_i ),
            t_tis      = v_tis,
            f_trans    = perc * ( t_tis - t_ves );

          f_mv[a][i]   += f_trans;
          f_mut_mv[a][i] += f_trans;

        }
      } else if ( !V_is_term[a] ) {
        // intermediate segment of vtree
        for ( int i = 0; i < V_lenv[a] ; ++i ) {
          const real_type
            lambda_i   = ( 1.0 * i ) / ( V_lenv[a] ),
            x_i        = ( vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].x + lambda_i *
                           ( vessels.vtree.nodes[ vessels.vtree.arcs[a].init ].x - vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].x ) ),
            y_i        = ( vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].y + lambda_i *
                           ( vessels.vtree.nodes[ vessels.vtree.arcs[a].init ].y - vessels.vtree.nodes[ vessels.vtree.arcs[a].term ].y ) ),
            v_ves      = vtree_data[a][i],
            t_ves      = v_ves,
            v_tis      = tissue_data.interpolate_on01 ( x_i, y_i ),
            t_tis      = v_tis,
            f_trans    = perc * ( t_tis - t_ves );

          f_mv[a][i]   += f_trans;
          f_mut_mv[a][i] += f_trans;

        }
      } else {
        // terminal segment of vtree: no heat diffusion
        // do nothing
      }
    }
  }
  // before: used averageing over grid cells on trees.

  dwarfmv.setZero();
  M_block.applyAdd (  f_mv, dwarfmv );
  Me_block.applyAdd ( f_mv, dwarfmv );
  dwarfmv *= 0.5 / C_tis ;          // temperature equilibration within one time step.
  brhsmv += dwarfmv;                // is this correct??

  dwarfmv.setZero();
  M_block.applyAdd (  f_mut_mv, dwarfmv );
  Me_block.applyAdd ( f_mut_mv, dwarfmv );
  dwarfmv *= 0.5;                   // temperature equilibration within one time step.
  mutsrc += dwarfmv;                // is this correct??
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::save_all_data ( char* ofn, const int prec ) {
  // save tissue and vessel data
  ofstream output ( ofn );
  output.precision ( prec );

  if ( output.is_open() ) {

    const int nx = tissue.getWidth();
    output << nx << " " << nx << endl;
    for ( int i = 0; i < nx; i++ )
      for ( int j = 0; j < nx; j++ )
        output << tissue_data.get ( i, j ) << endl;

    output << A_narcs << endl;
    for ( unsigned int i = 0; i < A_narcs; i++ ) {
      output << A_lenv[i] << endl;
      for ( int j = 0; j < A_lenv[i]; j++ )
        output << atree_data[i][j] << endl;
    }

    output << V_narcs << endl;
    for ( unsigned int i = 0; i < V_narcs; i++ ) {
      output << V_lenv[i] << endl;
      for ( int j = 0; j < V_lenv[i]; j++ )
        output << vtree_data[i][j] << endl;
    }

  }
  output.close();
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::read_all_data ( char* ifn ) {
  // save tissue and vessel data

  FILE * indat;

  if ( ! ( indat = fopen ( ifn, "r" ) ) ) {
    cerr << "Error: could not open " << ifn << " for reading tree." << endl;
  } else {
    int nx;
    real_type dummy;

    fscanf ( indat, "%d\n", &nx );
    fscanf ( indat, "%d\n", &nx );
    if ( nx != tissue.getWidth() ) {
      cerr << "tissue_vessel_geometry::read_all_data: dimensions of tissue do not match. No tissue data loaded." << endl;
      for ( int i = 0; i < nx; i++ ) {
        for ( int j = 0; j < nx; j++ ) {
          fscanf ( indat, "%lf\n", &dummy );
        }
      }
    } else {
      for ( int i = 0; i < nx; i++ ) {
        for ( int j = 0; j < nx; j++ ) {
          fscanf ( indat, "%lf\n", &dummy );
          tissue_data.set ( i, j, dummy );
        }
      }
    }
    int read_a_narcs;
    fscanf ( indat, "%d\n", &read_a_narcs );
    if ( static_cast<unsigned int> ( read_a_narcs ) != A_narcs ) {
      cerr << "tissue_vessel_geometry::read_all_data: number of arcs in atree does not match." << endl;
    }
    for ( unsigned int i = 0; i < A_narcs; i++ ) {
      int a_lenvi;
      fscanf ( indat, "%d\n", &a_lenvi );
      if ( a_lenvi != A_lenv[i] ) {
        cerr << "tissue_vessel_geometry::read_all_data: dimensions of vessel discretization do not match. No atree data loaded." << endl;
      } else {
        for ( int j = 0; j < A_lenv[i]; j++ ) {
          fscanf ( indat, "%lf\n", &dummy );
          atree_data[i][j] = dummy;
        }
      }
    }
    int read_v_narcs;
    fscanf ( indat, "%d\n", &read_v_narcs );
    if ( static_cast<unsigned int> ( read_v_narcs ) != V_narcs ) {
      cerr << "tissue_vessel_geometry::read_all_data: number of arcs in vtree does not match." << endl;
    }
    for ( unsigned int i = 0; i < V_narcs; i++ ) {
      int v_lenvi;
      fscanf ( indat, "%d\n", &v_lenvi );
      if ( v_lenvi != V_lenv[i] ) {
        cerr << "tissue_vessel_geometry::read_all_data: dimensions of vessel discretization do not match. No vtree data loaded." << endl;
      } else {
        for ( int j = 0; j < V_lenv[i]; j++ ) {
          fscanf ( indat, "%lf\n", &dummy );
          vtree_data[i][j] = dummy;
        }
      }
    }

  }

  fclose ( indat );
}


template <typename real_type>
void tissue_vessel_geometry<real_type>::add_corresponding_tissue_sources ( aol::MultiVector<real_type> &A_mutsrc, aol::MultiVector<real_type> &V_mutsrc, qc::ScalarArray<real_type, qc::QC_2D> &tis_rhs, const real_type C_tis, const real_type /*C_ves*/ ) {

  // this can be optimized, no vector copies necessary.
  aol::MultiVector<real_type> negwA ( A_mutsrc, aol::DEEP_COPY ), negwV ( V_mutsrc, aol::DEEP_COPY ); // deep copies

  for ( unsigned int i = 0; i < A_narcs; i++ ) {
    negwA[i] *= ( -1.0 ) * ( 1.0 / C_tis ) / vessels.atree.arcs[i].get_cross_sect_area();
  }
  for ( unsigned int i = 0; i < V_narcs; i++ ) {
    negwV[i] *= ( -1.0 ) * ( 1.0 / C_tis  ) / vessels.vtree.arcs[i].get_cross_sect_area();
  }

  add_mvec_phi_products ( negwA, negwV, tis_rhs );

}

template class tissue_vessel_geometry<double>;
// template class tissue_vessel_geometry<float>;

