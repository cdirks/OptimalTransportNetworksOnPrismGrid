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
 * Determining flow splitting ratios of vessel trees generated by
 * CCO via solving an optimization problem
 * ******************************************************************* */

#include "vesseltree.h"
#include "epswriter.h"
#include "FEOpInterface.h"
#include "solver.h"
#include "parameterParser.h"
#include "tissue_vessel_geometry.h"
#include <configurators.h>

// #define VERBOSE_COMPUT 1
// #define V_VERBOSE_TIME 1
#define SPEED_TEST 1
#define VERBOSE_TIME 1


inline bool components_positive ( aol::Vector<double> &Vec ) {
  for ( int i = 0; i < Vec.size(); i++ ) {
    if ( Vec[i] < 0 ) {
      return ( false );
    }
  }
  return ( true );
}


void find_descent_direction ( aol::Vector<double> &Dir, aol::Vector<double> &Grad, aol::Vector<double> &Point,
                              aol::Vector<double> &l_a, int s_a, aol::Vector<double> &l_v, int s_v, double ce ) {

                                aol::Vector<double> Vertex ( Point.size() );
                                aol::Vector<double> tempdir ( Point.size() );

  int ia_min = -1, iv_min = -1;
  double val_min = 0.0, val_here = 0.0;
  for ( int ia = 0; ia < s_a; ia++ ) {
    for ( int iv = 0; iv < s_v; iv++ ) {

      Vertex.setZero();
      Vertex.set ( ia    , ce / l_a[    ia] );
      Vertex.set ( s_a + iv, ce / l_v[s_a+iv] );
      tempdir = Vertex;
      tempdir -= Point;
      val_here = ( Grad * tempdir );

      if ( val_here < val_min ) { //found better vertex
        val_min = val_here;
        ia_min = ia;
        iv_min = iv;
      }
    }
  }

  Dir.setZero();
  if ( val_min < 0.0 ) { // descent direction found
    Dir.set ( ia_min    , ce / l_a[    ia_min] );
    Dir.set ( s_a + iv_min, ce / l_v[s_a+iv_min] );
    Dir -= Point;
    //    cerr << ia_min << " " << iv_min << " " << (Grad * Dir) << endl;
  } else {
    cerr << "no descent direction found" << endl;
  }
}


void Armijo_project_step ( aol::Op< aol::Vector<double> > &atlas, aol::Vector<double> &s_k, aol::Vector<double> &s_kp1,
                           aol::Vector<double> &l_a, int s_a, aol::Vector<double> &l_v, int s_v, double ce,
                           double sigma, double beta, double s, int max_steps, int &n_armijo_steps ) {
  // todo: modify this such that old step size is used and step size can be increased
  // currently, stepsize is decreased, starting with s * beta ^ (number of steps before / 2)

  aol::Vector<double> g_k ( s_k.size() ), d_k ( s_k.size() ), dummy ( s_k.size() );
  atlas.apply ( s_k, g_k );
  double f_orig = s_k * g_k;
  double f_new = f_orig;
  g_k *= 2.0;

  find_descent_direction ( d_k, g_k, s_k, l_a, s_a, l_v, s_v, ce );

  double tangnt_slope =  g_k * d_k ; // ??
#ifdef VERBOSE_COMPUT
  cerr << "f_orig = " << f_orig << ", tangent slope: " << tangnt_slope << " " ;
#endif

  bool Armijo_stop = false;
  int step = n_armijo_steps / 2 ; // integer division!

  while ( !Armijo_stop && ( step < max_steps ) ) {
    double alpha = s * pow ( beta, step );
    s_kp1  = d_k;
    s_kp1 *= alpha;

    s_kp1 += s_k;

    if ( components_positive ( s_kp1 ) ) {
      atlas.apply ( s_kp1, dummy );
      f_new = s_kp1 * dummy;
      double secant_slope = ( f_new - f_orig ) / alpha;
      // #ifdef VERBOSE_COMPUT
      cerr << "."; // did one Armijo step
      // #endif
      Armijo_stop = ( ( secant_slope / tangnt_slope ) > sigma  );
    } else {
      cerr << ","; // outside of feasible set
    }
    step++;
  }
  cerr << endl;
  // if step = 1, 0th step was accepted. maybe can do better with bigger step size.

  if ( step == max_steps ) {
    // to avoid numerical difficulties, set step size to zero.
    cerr << "max number of steps reached. use zero descent" << endl;
    f_new = f_orig;
    s_kp1 = s_k;
  }

  n_armijo_steps = step;

#ifdef VERBOSE_COMPUT
  cerr << endl << "Armijo iteration complete in " << step << " steps. f_orig was " << f_orig << ", f_new is " << f_new
  << ", relative decrease: " << f_new / f_orig << endl;
#endif
}


void set_up_A_matrix ( aol::Matrix<double> &A, tissue_vessel_geometry<double> &tvgeom ) {
  double h = tvgeom.tissue.H();
  int wi = tvgeom.tissue.getWidth(); // quadratic grids only.
  int s_a = tvgeom.A_term.size(), s_v = tvgeom.V_term.size();

  qc::GridDefinition::OldFullElementIterator elit;
  std::vector< segment_in_cell<double> >::iterator segit;

  for ( elit = tvgeom.tissue.begin(); elit != tvgeom.tissue.end(); ++elit ) {
    qc::Element &El = *elit;
    int ix = El.x();
    int jy = El.y();

    for ( segit = tvgeom.vessels.atree.segments_in_grid.at_begin ( ix, jy );
          segit != tvgeom.vessels.atree.segments_in_grid.at_end ( ix, jy ); ++segit ) {
      if ( tvgeom.A_is_term[segit->arcno] ) {
        int k = -1;
        for ( int m = 0; m < s_a; m++ ) {
          if ( tvgeom.A_term[m] == segit->arcno ) {
            k = m;
          }
        }
        double
          lxa = segit->x_in,
          lya = segit->y_in,
          lxb = segit->x_out,
          lyb = segit->y_out,

          val1 =     lxa  *     lya  + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * (    lya  +     lyb ) ) +     lxb *    lyb ,
          val2 = ( 1. - lxa ) *     lya  + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * (    lya  +     lyb ) ) + ( 1. - lxb ) *    lyb ,
          val3 =     lxa  * ( 1. - lya ) + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) +     lxb * ( 1. - lyb ),
          val4 = ( 1. - lxa ) * ( 1. - lya ) + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) + ( 1. - lxb ) * ( 1. - lyb );

        val1 *= h / 6.0 ;
        val2 *= h / 6.0 ;
        val3 *= h / 6.0 ;
        val4 *= h / 6.0 ;

        A.add ( ix + jy * wi, k, val1 );
        if ( ix + 1 < tvgeom.tissue.getWidth() ) {
          A.add ( ix + 1 + jy * wi, k, val2 );
        }
        if ( jy + 1 < tvgeom.tissue.getWidth() ) {
          A.add ( ix + ( jy + 1 ) * wi, k, val3 );
        }
        if ( ( ix + 1 < tvgeom.tissue.getWidth() ) && ( jy + 1 < tvgeom.tissue.getWidth() ) ) {
          A.add ( ( ix + 1 ) + ( jy + 1 ) * wi, k, val4 );
        }
      }
    }

    for ( segit = tvgeom.vessels.vtree.segments_in_grid.at_begin ( ix, jy );
          segit != tvgeom.vessels.vtree.segments_in_grid.at_end ( ix, jy ); ++segit ) {
      if ( tvgeom.V_is_term[segit->arcno] ) {
        int k = -1;
        for ( int m = 0; m < s_v; m++ ) {
          if ( tvgeom.V_term[m] == segit->arcno ) {
            k = s_a + m;
          }
        }

        double
          lxa = segit->x_in,
          lya = segit->y_in,
          lxb = segit->x_out,
          lyb = segit->y_out,

          val1 =     lxa  *     lya  + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * (    lya  +     lyb ) ) +     lxb *    lyb ,
          val2 = ( 1. - lxa ) *     lya  + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * (    lya  +     lyb ) ) + ( 1. - lxb ) *    lyb ,
          val3 =     lxa  * ( 1. - lya ) + 4. * ( .5 * (    lxa  +     lxb ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) +     lxb * ( 1. - lyb ),
          val4 = ( 1. - lxa ) * ( 1. - lya ) + 4. * ( .5 * ( ( 1. - lxa ) + ( 1. - lxb ) ) ) * ( .5 * ( ( 1. - lya ) + ( 1. - lyb ) ) ) + ( 1. - lxb ) * ( 1. - lyb );

        val1 *= -h / 6.0 ;
        val2 *= -h / 6.0 ;
        val3 *= -h / 6.0 ;
        val4 *= -h / 6.0 ; // why h^2 ??

        A.add ( ix + jy * wi, k, val1 );
        if ( ix + 1 < tvgeom.tissue.getWidth() ) {
          A.add ( ix + 1 + jy * wi, k, val2 );
        }
        if ( jy + 1 < tvgeom.tissue.getWidth() ) {
          A.add ( ix + ( jy + 1 ) * wi, k, val3 );
        }
        if ( ( ix + 1 < tvgeom.tissue.getWidth() ) && ( jy + 1 < tvgeom.tissue.getWidth() ) ) {
          A.add ( ( ix + 1 ) + ( jy + 1 ) * wi, k, val4 );
        }

      }
    }
  }
}

// seems to work correctly:

#ifndef SPEED_TEST
void solve_pressure_problem ( tissue_vessel_geometry<double> &tvgeom, aol::ParameterParser &parser, bool continue_calc ) {
#else
void solve_pressure_problem ( tissue_vessel_geometry<double> &tvgeom, aol::ParameterParser &parser, bool continue_calc, double t_offset ) {
  // for conv timing
  aol::StopWatch conv_timer;
  conv_timer.start();
#endif

  int n  = tvgeom.tissue.getNumberOfNodes();

  const double
  root_speed = parser.getDouble ( "root_speed" ),
               ce = root_speed * min ( tvgeom.vessels.atree.arcs[0].get_cross_sect_area(), tvgeom.vessels.vtree.arcs[0].get_cross_sect_area() );

  int s_a = tvgeom.A_term.size(),
            s_v = tvgeom.V_term.size(),
                  s = s_a + s_v;

  aol::SparseMatrix<double> A ( n, s );
  set_up_A_matrix ( A, tvgeom );

  typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature< double, qc::QC_2D, 3 > > ConfigType;
  const aol::StiffOp<ConfigType>  Lop ( tvgeom.tissue, aol::ASSEMBLED );

  //  cerr << "computing transpose of A  ... "; // (may take a while as transposeTo is not implemented efficiently for SparseMatrices)
  aol::SparseMatrix<double> A_t ( s, n );
  A.transposeTo ( A_t );

  //  cerr << "done." << endl;

//   qc::UniformGridSparseMatrix<double> L( tvgeom.tissue ); // n x n
//   L.clear();
//   Lop.assembleAddMatrix( L );
  aol::LagrangeZeroProjectSolver< double > L_in ( Lop, n, parser.getDouble ( "solv_epsilon" ), parser.getInt ( "max_solv_steps" ), true );
  // why not pass Lop??

  aol::Vector<double> l_a ( s ), l_v ( s );
  for ( int i = 0; i < s_a; i++ ) l_a[i]     = tvgeom.vessels.atree.length_of_arc ( tvgeom.A_term[i] );
  for ( int i = 0; i < s_v; i++ ) l_v[s_a+i] = tvgeom.vessels.vtree.length_of_arc ( tvgeom.V_term[i] );

  aol::Vector<double> s_vec ( s );

  if ( !continue_calc ) {

    for ( int i = 0; i < s_a; i++ ) s_vec[i]     = ce / ( s_a * tvgeom.vessels.atree.length_of_arc ( tvgeom.A_term[i] ) );
    for ( int i = 0; i < s_v; i++ ) s_vec[s_a+i] = ce / ( s_v * tvgeom.vessels.vtree.length_of_arc ( tvgeom.V_term[i] ) );

  } else {
    // read data from tree
    // do a conversion from speeds on the outflow segments to outflow per length rates
    // we can assume that the sum is ce, can't we?
    for ( int i = 0; i < s_a; i++ )
      s_vec[i]
      =  tvgeom.vessels.atree.arcs[ tvgeom.A_term[i] ].speed
         * tvgeom.vessels.atree.arcs[ tvgeom.A_term[i] ].get_cross_sect_area() / tvgeom.vessels.atree.length_of_arc ( tvgeom.A_term[i] );
    // correct?

    for ( int i = 0; i < s_v; i++ )
      s_vec[s_a+i]
      =  tvgeom.vessels.vtree.arcs[ tvgeom.V_term[i] ].speed
         * tvgeom.vessels.vtree.arcs[ tvgeom.V_term[i] ].get_cross_sect_area() / tvgeom.vessels.vtree.length_of_arc ( tvgeom.V_term[i] );
    // correct?
  }


  // set up problem to solve:
  {
    double armijo_sigma        = parser.getDouble ( "armijo_sigma" );
    double armijo_beta         = parser.getDouble ( "armijo_beta" );
    double armijo_s            = parser.getDouble ( "armijo_s" );
    int armijo_max_steps       = parser.getInt ( "armijo_max_steps" );
    double armijo_min_decrease = parser.getDouble ( "armijo_min_decrease" );

    aol::Composite_differentDim_Op< aol::Vector<double> > atlas;
    atlas.appendReference ( A, n );
    atlas.appendReference ( L_in, n );
    atlas.appendReference ( A_t, s );

    aol::Vector<double> s_rhs ( s ), s_vec1 ( s );

    atlas.apply ( s_vec, s_rhs );

    //    double r = s_vec * s_vec;
    //    cerr << "norm of s = " << r << endl;
    double r = s_rhs * s_vec;
#ifdef VERBOSE_COMPUT
    cerr << "function value = " << r  << endl; // no division by n
#endif

    bool stop = false;
    double r1 = r;
    int step = 0;
    int n_armijo_steps = 0;

    while ( !stop ) {
      step++;
      aol::StopWatch inner_timer;
      inner_timer.start();
      Armijo_project_step ( atlas, s_vec, s_vec1, l_a, s_a, l_v, s_v, ce, armijo_sigma, armijo_beta, armijo_s, armijo_max_steps, n_armijo_steps );
#ifdef VERBOSE_COMPUT
      cerr << l_a * s_vec1 << " " << l_v * s_vec1 << " CHECK: " << ( l_a * s_vec1 ) - ( l_v * s_vec1 ) << endl;
#endif
      atlas.apply ( s_vec1, s_rhs );
      r1 = s_vec1 * s_rhs;
      stop = ( ( r1 / r ) > armijo_min_decrease );
      // stop = ( r1 <= parser.getDouble("abs_stop") );
      s_vec = s_vec1;
      r = r1;
      inner_timer.stop();
#ifdef V_VERBOSE_TIME
      cerr << "step " << step << " took " << inner_timer.elapsedCpuTime() << " seconds." << endl;
#endif

#ifdef SPEED_TEST
      conv_timer.stop();
      double ela = conv_timer.elapsedCpuTime() + t_offset;
      conv_timer.cont();
      cout << ela << " " << log ( r  ) << endl; // no division by n
#endif
    }

    cerr << "function value = " << r  << endl; // no division by n

  }

  { // write results to tree:
    // do a conversion from outflow per length at outflow segments to speeds there
    tvgeom.vessels.atree.set_all_speeds ( -1000.0 );
    for ( int i = 0; i < tvgeom.A_term.size(); i++ ) {
      tvgeom.vessels.atree.arcs[ tvgeom.A_term[i] ].speed
      = tvgeom.vessels.atree.length_of_arc ( tvgeom.A_term[i] )
        / tvgeom.vessels.atree.arcs[ tvgeom.A_term[i] ].get_cross_sect_area() * s_vec[i];
    }

    tvgeom.vessels.vtree.set_all_speeds ( -1000.0 );
    for ( int i = 0; i < tvgeom.V_term.size(); i++ ) {
      tvgeom.vessels.vtree.arcs[ tvgeom.V_term[i] ].speed
      = tvgeom.vessels.vtree.length_of_arc ( tvgeom.V_term[i] )
        / tvgeom.vessels.vtree.arcs[ tvgeom.V_term[i] ].get_cross_sect_area() * s_vec[ s_a + i];
    }

    // combine up from first given segments
    tvgeom.adapt_speeds_up_from_termsegm();
  }

  { // compute pressure and save
    qc::ScalarArray<double, qc::QC_2D> dummy ( tvgeom.tissue );
    A.apply ( s_vec, dummy );
    L_in.apply ( dummy, tvgeom.tissue_data );

    //     char ofn_p[1024], ofn_pc[1024];
    //     sprintf(ofn_p, "output/pressure_%d_opti.pgm", tvgeom.tissue.getGridDepth() );
    //     sprintf(ofn_pc, "output/pressure_%d_opti.ppm", tvgeom.tissue.getGridDepth() );
    //     img_scaled_save_ctrscl( pressure, ofn_p, ofn_pc, 0.0, 0.0 );
  }
}


int main ( int argc, char** argv ) {
  if ( argc == 2 ) {
    try {
      aol::StopWatch timer;
      timer.start();
      cerr.precision ( 10 );

      aol::ParameterParser parser ( argv[1] );
      char fn_atree[1024], fn_vtree[1024];
      parser.getString ( "atree_file", fn_atree );
      parser.getString ( "vtree_file", fn_vtree );

      int min_grid_level = parser.getInt ( "mingridlevel" );
      int max_grid_level = parser.getInt ( "maxgridlevel" );

      bool
      firstlevel = true,
                   continue_calc = static_cast<bool> ( parser.getInt ( "continue" ) );

      for ( int i = min_grid_level; i <= max_grid_level; i++ ) {

        aol::StopWatch level_timer;

        //      cerr << "================================= " << endl;
        cerr << "solving problem on grid level " << i << ":" <<  endl;
        //      cerr << "================================= " << endl;


        // on first level, load tree from input file. save to partial result, next step: load from partial result.
        if ( !firstlevel ) {
          sprintf ( fn_atree, "output/pressure_%d_opti_a.tree", i - 1 );
          sprintf ( fn_vtree, "output/pressure_%d_opti_v.tree", i - 1 );
        }

        cerr << "Please ignore the following warnings (if any) - we will not need the discretization of the vessels" << endl;
        tissue_vessel_geometry<double> tvgeom ( fn_atree, fn_vtree, i ); // i = current level.

        level_timer.start();

#ifndef SPEED_TEST
        solve_pressure_problem ( tvgeom, parser, continue_calc );
#else
        timer.stop();
        double offs = timer.elapsedCpuTime();
        timer.cont();
        solve_pressure_problem ( tvgeom, parser, continue_calc, offs );
        cout << endl;
#endif
        level_timer.stop();


        { // save partial result.
          char ofn_a[1024], ofn_v[1024], ofn_e[1024], ofn_d[1024];
          sprintf ( ofn_a, "output/pressure_%d_opti_a.tree", i );
          sprintf ( ofn_v, "output/pressure_%d_opti_v.tree", i );
          sprintf ( ofn_e, "output/pressure_%d_opti_speeds.eps", i );
          sprintf ( ofn_d, "output/pressure_%d_opti_data.dat", i );
          tvgeom.vessels.save_trees ( ofn_a, ofn_v );
          tvgeom.vessels.plot_tree_speeds ( ofn_e );
          tvgeom.save_all_data ( ofn_d, 20 );
        }

#ifdef VERBOSE_TIME
        cerr << "Optimization on grid level " << i << " took " << level_timer.elapsedCpuTime() << " seconds." << endl;
#endif

        firstlevel = false;
        continue_calc = true;

      }

      timer.stop();
      cerr << endl << "Total optimization took " << timer.elapsedCpuTime() << " seconds." << endl;

    } catch ( aol::Exception &ex ) {
      ex.dump();
    }
  } else {
    cerr << "usage: usevesseltree_det_flowsplit_mgrid <parameter file>" << endl;
  }

  return ( 0 );
}

