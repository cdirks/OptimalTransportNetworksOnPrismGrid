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
 * coupled advection-diffusion problem with pulsed velocities
 * ******************************************************************* */


#include "parameterParser.h"

#include "tissue_vessel_geometry.h"
#include "ellam_matrices.h"
#include <configurators.h>

// #define VERBOSE
#define USE_PCG 1

void compute_aroot_inflow_values ( aol::Vector<double> &Rp_in, double inflow_val, const bool inflow_pulse = false ) {
  for ( int t = 0; t < Rp_in.size(); t++ )
    Rp_in[t] = inflow_val;  // 0.5;

  if ( inflow_pulse ) {
    if ( Rp_in.size() > 1000 )
      for ( int t = 0; t < 1000; ++t )
        Rp_in[t] += 0.5 * ( 1 - cos ( 2.0 * aol::NumberTrait<double>::pi * t / 1000.0 ) );

    if ( Rp_in.size() > 5000 )
      for ( int t = 4000; t < 5000; ++t )
        Rp_in[t] += 0.5 * ( 1 - cos ( 2.0 * aol::NumberTrait<double>::pi * t / 1000.0 ) );

    if ( Rp_in.size() > 9000 )
      for ( int t = 8000; t < 9000; ++t )
        Rp_in[t] += 0.5 * ( 1 - cos ( 2.0 * aol::NumberTrait<double>::pi * t / 1000.0 ) );

    if ( Rp_in.size() > 13000 )
      for ( int t = 12000; t < 13000; ++t )
        Rp_in[t] += 0.5 * ( 1 - cos ( 2.0 * aol::NumberTrait<double>::pi * t / 1000.0 ) );

    if ( Rp_in.size() > 17000 )
      for ( int t = 16000; t < 17000; ++t )
        Rp_in[t] += 0.5 * ( 1 - cos ( 2.0 * aol::NumberTrait<double>::pi * t / 1000.0 ) );
  }
}


#if 0
void set_probe_influence ( aol::ParameterParser &params, qc::ScalarArray<double, qc::QC_2D> &tis_data ) {
  // RF probe keeps area at constant temperature
  const int
    start_x = static_cast<int> ( params.getDouble ( "start_x" ) * tis_data.getNumX() ),
    stop_x  = static_cast<int> ( params.getDouble ( "stop_x" )  * tis_data.getNumX() ),
    start_y = static_cast<int> ( params.getDouble ( "start_y" ) * tis_data.getNumY() ),
    stop_y  = static_cast<int> ( params.getDouble ( "stop_y" )  * tis_data.getNumY() );
  const double value = params.getDouble ( "appl_source" );

  for ( int i = start_x; i < stop_x; ++i )
    for ( int j = start_y; j < stop_y; ++j )
      tis_data.set ( i, j, value );
}
#endif


void add_probe_source_terms ( aol::ParameterParser &params, qc::ScalarArray<double, qc::QC_2D> &tis_rhs, qc::ScalarArray<double, qc::QC_2D> &tis_data, double tau, double h ) {
  // heating by RF probe
  const double
    d_sp_x  = params.getDouble ( "spos_x" ),
    d_sp_y  = params.getDouble ( "spos_y" ),
    s_sigma = params.getDouble ( "s_sigma" ),
    s_value = params.getDouble ( "s_value" ),
    e       = exp ( 1.0 ),
    //    r2p     = sqrt( 2 * aol::NumberTrait<double>::pi ),
    th2     = tau * h * h;
  const int
  spos_x  = static_cast<int> ( d_sp_x  * tis_data.getNumX() ),
            spos_y  = static_cast<int> ( d_sp_y  * tis_data.getNumX() );


  int
    start_x = static_cast<int> ( ( d_sp_x - s_sigma ) * tis_data.getNumX() ),
    stop_x  = static_cast<int> ( ( d_sp_x + s_sigma ) * tis_data.getNumX() ) + 1,
    start_y = static_cast<int> ( ( d_sp_y - s_sigma ) * tis_data.getNumY() ),
    stop_y  = static_cast<int> ( ( d_sp_y + s_sigma ) * tis_data.getNumY() ) + 1;
  if ( start_x < 0 ) start_x = 0;
  if ( stop_x > tis_data.getNumX() ) stop_x = tis_data.getNumX();
  if ( start_y < 0 ) start_y = 0;
  if ( stop_y > tis_data.getNumY() ) stop_y = tis_data.getNumY();

  double sum = 0.0;
  int numpts = 0;

  for ( int i = start_x; i < stop_x; ++i )
    for ( int j = start_y; j < stop_y; ++j )
      if ( ( h * ( i - spos_x ) * h * ( i - spos_x ) + h * ( j - spos_y ) * h * ( j - spos_y ) ) < ( 0.16 * s_sigma * s_sigma ) ) {
        numpts++;
        sum = tis_data.get ( i, j );
      }

  if ( ( sum / numpts ) < s_sigma ) {
    for ( int i = start_x; i < stop_x; ++i ) {
      for ( int j = start_y; j < stop_y; ++j ) {
        const double dist2 = h * ( i - spos_x ) * h * ( i - spos_x ) + h * ( j - spos_y ) * h * ( j - spos_y );
        if ( dist2 < ( s_sigma * s_sigma ) ) {
          double intensity = e;
          intensity *= exp ( - ( 1.0 / ( 1 - dist2 / ( s_sigma * s_sigma ) ) ) );
          intensity *= s_value * th2;
          tis_rhs.add ( i, j, intensity );
        }
      }
    }
  }
}


int main ( int argc, char** argv ) {
  if ( argc == 2 ) {
    try {

      // load trees etc

      aol::ParameterParser params ( argv[1] );
      char atree_fn[512], vtree_fn[512], treefn_pattern[512], tissfn_pattern[512], tisscfn_pattern[512], datafn_pattern[512];
      params.getString ( "atree_file", atree_fn );
      params.getString ( "vtree_file", vtree_fn );
      params.getString ( "tree_files", treefn_pattern );
      params.getString ( "tiss_files", tissfn_pattern );
      params.getString ( "tiss_c_files", tisscfn_pattern );
      params.getString ( "data_files", datafn_pattern );

      const double
        //        scalefactor = params.getDouble("colorscalefactor"),
        C_ves       = params.getDouble ( "C_ves" ),
        C_tis       = params.getDouble ( "C_tis" );

      const int
        level_tissue = params.getInt ( "leveltissue" ),
        level_vessels = params.getInt ( "levelvessels" ),
        timesteps = params.getInt ( "num_timesteps" );

      cerr << "Setting up tvgeom ... " << endl;
      tissue_vessel_geometry<double> tvgeom ( atree_fn, vtree_fn, level_tissue, level_vessels );


      // for test only:
      // =============
      for ( unsigned int i = 0; i < tvgeom.vessels.atree.arcs.size(); i++ )
        tvgeom.vessels.atree.arcs[i].rad *= params.getDouble ( "rad_adaption" );
      for ( unsigned int i = 0; i < tvgeom.vessels.vtree.arcs.size(); i++ )
        tvgeom.vessels.vtree.arcs[i].rad *= params.getDouble ( "rad_adaption" );
      cerr << "a root radius = " << tvgeom.vessels.atree.arcs[0].rad << ", v radius = " << tvgeom.vessels.vtree.arcs[0].rad << endl;


      double tau = tvgeom.smallest_vessel_h;
      if ( pow ( 0.5, level_tissue + 1.0 ) < tau ) tau = pow ( 0.5, level_tissue + 1.0 );    // make sure tau is sufficiently small.
      cerr << "tau = " << tau << endl;


      // set up solvers for diffusion problem

      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature< double, qc::QC_2D, 3 > > ConfigType;
      const aol::StiffOp<ConfigType>  L ( tvgeom.tissue, aol::ASSEMBLED );
      const aol::MassOp<ConfigType>   M ( tvgeom.tissue, aol::ASSEMBLED );

      const double kappa = params.getDouble ( "kappa" );

      aol::LinCombOp< aol::Vector <double>, aol::Vector <double> > A;
      A.appendReference ( M ); A.appendReference ( L, kappa * tau );
#ifdef USE_PCG
      // here, pcg is significantly faster than cg.
      aol::DiagonalMatrix<double> Diag ( tvgeom.tissue.getNumberOfNodes() );
      {
        qc::UniformGridSparseMatrix<double> Lmat ( tvgeom.tissue ), Mmat ( tvgeom.tissue );
        L.assembleAddMatrix ( Lmat ); M.assembleAddMatrix ( Mmat );
        for ( int d = 0; d < Lmat.getNumRows(); d++ )
          Diag.set ( d, d, 1.0 / ( Mmat.get ( d, d ) + kappa * tau * Lmat.get ( d, d ) ) );
        // do not need Lmat, Mmat any more.
      }
      aol::PCGInverse< aol::Vector <double> > Diffusion_solver ( A, Diag, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );
#else
      aol::CGInverse< aol::Vector <double> > Diffusion_solver ( A, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );
#endif
      Diffusion_solver.setQuietMode ( static_cast<bool> ( params.getInt ( "solverquiet" ) ) );
      Diffusion_solver.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

      // initial data

      tvgeom.tissue_data.setAll ( params.getDouble ( "init_tissue_temp" ) );

      for ( unsigned int i = 0; i < tvgeom.A_narcs; ++i ) {
        if ( tvgeom.A_is_term[i] ) {
          for ( int j = 0; j < tvgeom.A_lenv[i]; ++j ) {
            tvgeom.atree_data[i][j] = ( tvgeom.A_lenv[i] - 1.0 - j ) * params.getDouble ( "init_vessel_temp" ) / ( tvgeom.A_lenv[i] - 1.0 );
            // ATTENTION:  might be wrong. need 1 - (j+1)/n
          }
        } else {
          tvgeom.atree_data[i].setAll ( params.getDouble ( "init_vessel_temp" ) );
        }
      }
      for ( unsigned int i = 0; i < tvgeom.V_narcs; ++i ) {
        if ( tvgeom.V_is_term[i] ) {
          for ( int j = 0; j < tvgeom.V_lenv[i]; ++j ) {
            tvgeom.vtree_data[i][j] = j * params.getDouble ( "init_vessel_temp" ) / ( tvgeom.V_lenv[i] - 1.0 );
            // ATTENTION:  might be wrong. need sth different
          }
        } else {
          tvgeom.vtree_data[i].setAll ( params.getDouble ( "init_vessel_temp" ) );
        }
      }

      // compute atree root inflow values

      aol::Vector<double> Rp_in ( timesteps + 1 );
      compute_aroot_inflow_values ( Rp_in, params.getDouble ( "inflow_temp" ), static_cast<bool> ( params.getInt ( "inflow_pulse" ) ) );


      // set up saving stuff

      int save_every = params.getInt ( "save_every" );

      aol::FullMatrix<double> rot ( 3, 3 );
      {
        const double
          xang_deg =  params.getDouble ( "xang_deg" ),
          yang_deg =  params.getDouble ( "yang_deg" ),
          zang_deg =  params.getDouble ( "zang_deg" );
        aol::EpsWriter::setRotationMatrix ( rot, xang_deg, yang_deg, zang_deg );
      }

      aol::Vector<double> shift ( 3 );
      shift[0] = 0.5; shift[1] = 0.5, shift[2] = 0.0;

      //      double zzoom = params.getDouble("zzoom");

      aol::EpsWriter::writeScalePPM ( aol::HSV_BLUE_TO_RED, "output/heightscale.ppm" );

      aol::StopWatch timer;
      timer.start();


      //   ====================
      //    compute time steps
      //   ====================

      aol::Vector<double> a_orig_speeds ( tvgeom.A_narcs ), v_orig_speeds ( tvgeom.V_narcs );
      for ( unsigned int i = 0 ; i < tvgeom.A_narcs; ++i ) {
        a_orig_speeds[i] = tvgeom.vessels.atree.arcs[i].speed;
      }
      for ( unsigned int i = 0 ; i < tvgeom.V_narcs; ++i ) {
        v_orig_speeds[i] = tvgeom.vessels.vtree.arcs[i].speed;
      }

      aol::MultiVector<double>
      A_brhsmv ( tvgeom.A_lenv ),
      V_brhsmv ( tvgeom.V_lenv ),
      A_dwarfmv ( tvgeom.A_lenv ),
      V_dwarfmv ( tvgeom.V_lenv ),
      A_mutsrc ( tvgeom.A_lenv ),
      V_mutsrc ( tvgeom.V_lenv );

      qc::ScalarArray<double, qc::QC_2D>
      tis_rhs ( tvgeom.tissue ),
      tis_dwarf ( tvgeom.tissue );

      for ( int t = 1; t <= timesteps; t++ ) {

        // vary velocities
        const double cur_fac = 0.6 + 0.4 * cos ( 2 * aol::NumberTrait<double>::pi * t / 1000 );
        for ( unsigned int i = 0 ; i < tvgeom.A_narcs; ++i ) {
          tvgeom.vessels.atree.arcs[i].speed = cur_fac * a_orig_speeds[i];
        }
        for ( unsigned int i = 0 ; i < tvgeom.V_narcs; ++i ) {
          tvgeom.vessels.vtree.arcs[i].speed = cur_fac * v_orig_speeds[i];
        }

        // set up solvers for advection problem

        aol::BlockOp< double, aol::Matrix<double> > A_M_block ( tvgeom.A_narcs, tvgeom.A_narcs );
        aol::BlockOp< double, aol::Matrix<double> > A_Me_block ( tvgeom.A_narcs, tvgeom.A_narcs );
        aol::BlockOp< double, aol::Matrix<double> > V_M_block ( tvgeom.V_narcs, tvgeom.V_narcs );
        aol::BlockOp< double, aol::Matrix<double> > V_Me_block ( tvgeom.V_narcs, tvgeom.V_narcs );

        build_blockops_recurse ( tvgeom.vessels.atree, tvgeom.A_lenv, 0, tau, A_M_block, A_Me_block, true );  // on atree
        build_blockops_recurse ( tvgeom.vessels.vtree, tvgeom.V_lenv, 0, tau, V_M_block, V_Me_block, false ); // not on atree

        bool solver_quiet = static_cast< bool >( params.getInt ( "solverquiet" ) );

        aol::CGInverse< aol::MultiVector<double> > A_solver ( A_M_block, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );
        aol::CGInverse< aol::MultiVector<double> > V_solver ( V_M_block, params.getDouble ( "cg_epsilon" ), params.getInt ( "max_cg_iter" ) );

        A_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );    A_solver.setQuietMode ( solver_quiet );
        V_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );    V_solver.setQuietMode ( solver_quiet );

#ifdef VERBOSE
        cerr << "Computing time step " << t << ":";
#endif

        tis_rhs.setZero();
        tis_dwarf.setZero();

        A_brhsmv.setZero();
        A_dwarfmv.setZero();
        V_brhsmv.setZero();
        V_dwarfmv.setZero();
        A_mutsrc.setZero();
        V_mutsrc.setZero();

        // set up source terms for the vesse trees and tissue

        // pseudo source terms: outflow out of terminal segments of atree
        tvgeom.add_atree_outflow_rhs ( A_brhsmv, A_mutsrc, A_M_block, A_Me_block, tau, C_tis, C_ves );
        tvgeom.add_vtree_inflow_rhs (  V_brhsmv, V_mutsrc, V_M_block, V_Me_block, tau, C_tis, C_ves );

        // transfer source terms ("through vessel walls")
        tvgeom.add_ves_tis_transfer_rhs ( A_brhsmv, A_mutsrc, A_M_block, A_Me_block, V_brhsmv, V_mutsrc, V_M_block, V_Me_block, tau, C_ves, C_tis, params.getDouble ( "perc" ) );

        // compute and set up corresponding source terms for the tissue
        tvgeom.add_corresponding_tissue_sources ( A_mutsrc, V_mutsrc, tis_rhs, C_tis, C_ves );

        // probe source (RF Applikator)
        if ( static_cast<bool> ( params.getInt ( "probe_present" ) ) )
          add_probe_source_terms ( params, tis_rhs, tvgeom.tissue_data, tau, tvgeom.tissue.H() );

        // inflow values (temperature) for atree through root. do not need corresponding tissue sources!
        {
          double h   = tvgeom.vessels.atree.length_of_arc ( 0 ) / ( tvgeom.A_lenv[0] - 1.0 ) ;
          double rho = tvgeom.vessels.atree.arcs[0].get_speed() * tau / h;
          double Aa  = tvgeom.vessels.atree.arcs[0].get_cross_sect_area();

          A_brhsmv[0][0] += Aa * ( Rp_in[t-1] * h * ( rho / 2.0 - rho * rho / 3.0 ) +
                                   Rp_in[t] * h * ( rho / 2.0 - rho * rho / 6.0 )    ) ;
          A_brhsmv[0][1] += Aa * ( Rp_in[t-1] * h * rho * rho / 3.0 + Rp_in[t] * h * rho * rho / 6.0 );
        }


        // diffusion time step
#ifdef VERBOSE
        cerr << " computing diffusion ... ";
#endif
        tis_dwarf.setZero();
        M.apply ( tvgeom.tissue_data, tis_dwarf );
        tis_dwarf += tis_rhs;

        tvgeom.tissue_data.setZero();
        Diffusion_solver.apply ( tis_dwarf, tvgeom.tissue_data );

        // advection time step
#ifdef VERBOSE
        cerr << "done, computing advection ... ";
#endif
        A_Me_block.apply ( tvgeom.atree_data, A_dwarfmv );
        A_dwarfmv += A_brhsmv;
        tvgeom.atree_data.setZero();
        A_solver.apply ( A_dwarfmv, tvgeom.atree_data );

        V_Me_block.apply ( tvgeom.vtree_data, V_dwarfmv );
        V_dwarfmv += V_brhsmv;
        tvgeom.vtree_data.setZero();
        V_solver.apply ( V_dwarfmv, tvgeom.vtree_data );

        // save time step

        if ( ! ( t % save_every ) ) { // only save every save_every timesteps
          cerr << "saving time step " << t << ": " ;
          char
          // tmp[1024],
          // tmpc[1024],
          datafn[1024];
          // profilefn[1024];

          sprintf ( datafn, datafn_pattern, t );
          tvgeom.save_all_data ( datafn, 10 );


#if 0
          sprintf ( tmp, tissfn_pattern, t );
          sprintf ( tmpc, tisscfn_pattern, t );
          img_scaled_save_ctrscl ( tvgeom.tissue_data, tmp, tmpc, params.getDouble ( "plot_ctr_temp" ), scalefactor, colTrans );
#endif
          cerr << "max = " << tvgeom.tissue_data.getMaxValue()
          << ", min = " << tvgeom.tissue_data.getMinValue()
          << ", avg = " << tvgeom.tissue_data.sum() / tvgeom.tissue.getNumberOfNodes() << endl;

#if 0
          double color_scale = params.getDouble ( "color_scale" );

          sprintf ( profilefn, treefn_pattern, t );

          if ( !static_cast<bool> ( params.getInt ( "topview" ) ) ) {
            tvgeom.save_advect_step ( profilefn, colTrans, zzoom, color_scale, rot, shift, true ); // true = plot temperatures.
          } else {
            zzoom *= 1.0;
            tvgeom.save_advect_step_topview_recursive ( profilefn, colTrans, color_scale, true );
          }
#endif

          //         tis_rhs.setOverflowHandling( aol::SCALE, 0, 255 );
          //         tis_rhs.setQuietMode( true );

          //         sprintf(outfilename, "output/rhs_%04d.pgm", t);
          //         tis_rhs.save( outfilename );
        }

#ifdef VERBOSE
        cerr << " done." << endl;
#endif

        delete_blockop_entries ( A_M_block );    delete_blockop_entries ( A_Me_block );
        delete_blockop_entries ( V_M_block );    delete_blockop_entries ( V_Me_block );

      }

      timer.stop();
      cerr << "took " << timer.elapsedCpuTime() << " seconds. T_final = " << timesteps * tau  << ". Cleaning up ... " << endl;

    } catch ( aol::Exception &ex ) {
      ex.dump();
    }

  } else {
    cerr << "usage: usevesseltree_advect_and_heat <parameter file>" << endl;
  }
  return ( 0 );
}

