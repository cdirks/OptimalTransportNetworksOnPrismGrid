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
 * 1D ELLAM on vessel trees, advection only
 * ******************************************************************* */

#include "parameterParser.h"

#include "ellam_matrices.h"
#include "tissue_vessel_geometry.h"

// part of this code may need to be adapted.

// #define VERBOSE 1

int main() {
  try {
    aol::ParameterParser params ( "par/advect.par" );
    char atree_fn[256], vtree_fn[256];
    params.getString ( "atree_file", atree_fn );
    params.getString ( "vtree_file", vtree_fn );

    int level_vessels = params.getInt ( "level_vessels" );

    cerr << "Setting up tvgeom ... " << endl;
    tissue_vessel_geometry<double> tvgeom ( atree_fn, vtree_fn, 5, level_vessels );
    // 5 = tissue level, irrelevant here because only vessels considered in calculation.

    normalize_speed ( tvgeom.vessels.atree, 0.85 );
    normalize_speed ( tvgeom.vessels.vtree, 0.85 );

    double smallest_h = tvgeom.smallest_vessel_h;

    double tau = smallest_h;
    if ( pow ( 0.5, level_vessels + 1.0 ) < tau ) tau = pow ( 0.5, level_vessels + 1.0 );

    cerr << "tau = " << tau << endl;
    int timesteps = params.getInt ( "num_timesteps" );

    aol::BlockOp< double, aol::Matrix<double> > A_M_block ( tvgeom.A_narcs, tvgeom.A_narcs );
    aol::BlockOp< double, aol::Matrix<double> > A_Me_block ( tvgeom.A_narcs, tvgeom.A_narcs );
    aol::BlockOp< double, aol::Matrix<double> > V_M_block ( tvgeom.V_narcs, tvgeom.V_narcs );
    aol::BlockOp< double, aol::Matrix<double> > V_Me_block ( tvgeom.V_narcs, tvgeom.V_narcs );

    build_blockops_recurse ( tvgeom.vessels.atree, tvgeom.A_lenv, 0, tau, A_M_block, A_Me_block, true );  // on atree
    build_blockops_recurse ( tvgeom.vessels.vtree, tvgeom.V_lenv, 0, tau, V_M_block, V_Me_block, false ); // not on atree

    cerr << "Setting up Vectors ... " << endl;
    aol::MultiVector<double>
    A_brhsmv ( tvgeom.A_lenv ),
    A_dwarfmv ( tvgeom.A_lenv ),
    V_brhsmv ( tvgeom.V_lenv ),
    V_dwarfmv ( tvgeom.V_lenv );

#if 0
    // set some initial profile if first vector is long enough
    if ( A_umv_new[0].size() > 50 ) { // in tvgeom!
      // a triangular pulse
      for ( int i = 10; i < 50; i++ ) {
        A_umv_new[0][i] = 1.0 - 0.05 * abs ( 30 - i );
      }
    }
    if ( A_umv_new[0].size() > 110 ) {
      // a rectangular pulse (discontinuity!)
      for ( int i = 70; i < 100; i++ ) {
        A_umv_new[0][i] = 1.0;
      }
    }
#endif

#if 0
    // initial profile for the vtree.
    for ( int k = 0; k < V_term.size(); k++ ) { // in tvgeom!
      //      cerr << "k = " << k << endl;
      int j = V_term[k];
      if ( V_umv_new[j].size() > 50 ) {
        // a triangular pulse
        for ( int i = 10; i < 50; i++ ) {
          V_umv_new[j][i] = 1.0 - 0.05 * abs ( 30 - i );
        }
      }
      if ( V_umv_new[j].size() > 110 ) {
        // a rectangular pulse
        for ( int i = 70; i < 100; i++ ) {
          V_umv_new[j][i] = 1.0;
        }
      }
    }
#endif

    // compute inflow values
    aol::Vector<double> Rp_in ( timesteps + 1 );
    for ( int t = 0; t < timesteps + 1; t++ ) {
      Rp_in[t] = 0.5 * ( 1.0 - cos ( ( 1.0 * t ) / params.getDouble ( "cos_denom" ) ) ); // adapt denominator to level
      //      Rp_in[t] =  ( (t < 50) ? ( 0.02 * t ) : 1.0 );
      // 0.5 * (1.0 - cos( (2.0 * PI * t) / 100.0) );
    }

    tvgeom.tissue_data.setAll ( 1.0 );

    bool solver_quiet = static_cast< bool >( params.getInt ( "solverquiet" ) );
    aol::CGInverse< aol::MultiVector< double > > A_solver ( A_M_block, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );
    A_solver.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    A_solver.setQuietMode ( solver_quiet );

    aol::CGInverse< aol::MultiVector< double > > V_solver ( V_M_block, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );
    V_solver.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    V_solver.setQuietMode ( solver_quiet );

    char ofn_pattern[256];
    params.getString ( "outfilenames", ofn_pattern );
    char datafn_pattern[256];
    params.getString ( "datafilenames", datafn_pattern );
    int save_every = params.getInt ( "save_every" );

    aol::FullMatrix<double> rot ( 3, 3 );
    {
      const double xang_deg =  params.getDouble ( "xang_deg" ),
                               yang_deg =  params.getDouble ( "yang_deg" ),
                                           zang_deg =  params.getDouble ( "zang_deg" );
      aol::EpsWriter::setRotationMatrix ( rot, xang_deg, yang_deg, zang_deg );
    }

    aol::Vector<double> shift ( 3 );
    shift[0] = 0.5; shift[1] = 0.5, shift[2] = 0.0;

    double zzoom = params.getDouble ( "zzoom" );

    aol::colorTrans colTrans = aol::HSV_BLUE_TO_RED;
    aol::EpsWriter::writeScalePPM ( colTrans, "output/heightscale.ppm" );

    aol::StopWatch timer;
    timer.start();

    for ( int t = 1; t <= timesteps; t++ ) {
#ifdef VERBOSE
      cout << "Computing time step " << t << endl;
      cout.precision ( 20 );

#if 0
      for ( unsigned int k = 0; k < my_vesseltree.atree.arcs.size(); k++ ) { // in tvgeom!
        for ( int l = 0; l < A_umv_new[k].length(); l++ ) {
          cout << A_umv_new[k][l] << endl;
        }
        cout << endl;
      }
#endif
#endif
      A_dwarfmv.setZero();      A_brhsmv.setZero();
      V_dwarfmv.setZero();      V_brhsmv.setZero();

      // inflow values for atree.
      {
        double h   = tvgeom.vessels.atree.length_of_arc ( 0 ) / ( tvgeom.A_lenv[0] - 1.0 ) ;
        double rho = tvgeom.vessels.atree.arcs[0].get_speed() * tau / h;
        double A   = tvgeom.vessels.atree.arcs[0].get_cross_sect_area();

        A_brhsmv[0][0] += A * ( Rp_in[t-1] * h * ( rho / 2.0 - rho * rho / 3.0 ) +
                                Rp_in[t] * h * ( rho / 2.0 - rho * rho / 6.0 )    ) ;
        A_brhsmv[0][1] += A * ( Rp_in[t-1] * h * rho * rho / 3.0 + Rp_in[t] * h * rho * rho / 6.0 );

      }

      // source terms
      // pseudo source terms: outflow out of terminal segments of atree

      tvgeom.add_atree_outflow_rhs ( A_brhsmv, A_dwarfmv, A_M_block, A_Me_block, tau, 1.0, 1.0 ); // C_tis, C_ves

      A_dwarfmv.setZero(); // compatibility with function, data written to this multivector is deleted.

      A_Me_block.apply ( tvgeom.atree_data, A_dwarfmv );
      A_dwarfmv += A_brhsmv;
      tvgeom.atree_data.setZero();
      A_solver.apply ( A_dwarfmv, tvgeom.atree_data );

      // pseudo source terms: inflow into terminal segments of vtree

      tvgeom.add_vtree_inflow_rhs ( V_brhsmv, V_dwarfmv, V_M_block, V_Me_block, tau, 1.0, 1.0 ); // C_tis, C_ves
      V_dwarfmv.setZero(); // compatibility with function, data written to this multivector is deleted.

      V_Me_block.apply ( tvgeom.vtree_data, V_dwarfmv );
      V_dwarfmv += V_brhsmv;
      tvgeom.vtree_data.setZero();
      V_solver.apply ( V_dwarfmv, tvgeom.vtree_data );

      if ( ! ( t % save_every ) ) { // only save every save_every timesteps

        double color_scale = params.getDouble ( "color_scale" );
        char outfilename[256];
        sprintf ( outfilename, ofn_pattern, t );

        tvgeom.save_advect_step ( outfilename, colTrans, zzoom, color_scale, rot, shift, true, false ); // last 2: save atree, save vtree.
        // tvgeom.save_advect_step_topview_recursive(outfilename, colScale, false ); // plot energy content

        sprintf ( outfilename, datafn_pattern, t );
        tvgeom.save_all_data ( outfilename, 15 );

      }
    }

    timer.stop();
    cerr << "took " << timer.elapsedCpuTime() << " seconds. T_final = " << timesteps * tau  << ". Cleaning up ... " << endl;

    delete_blockop_entries ( A_M_block );    delete_blockop_entries ( A_Me_block );
    delete_blockop_entries ( V_M_block );    delete_blockop_entries ( V_Me_block );

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return ( 0 );
}
