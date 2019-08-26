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
 * Heat conduction in tissue with line sources given on vessel tree
 * ******************************************************************* */

#include "vesseltree.h"
#include "tissue_vessel_geometry.h"
#include "parameterParser.h"
#include <configurators.h>

// #define VERBOSE
#define CONV_ANALY


void set_appropriate_source_values ( vesseltree<double> &my_vesseltree ) {
  double
    a_volume = 0.0,
    v_volume = 0.0,
    a_val = 0.0,
    v_val = 0.0;

  for ( unsigned int i = 0; i < my_vesseltree.atree.arcs.size(); i++ ) {
    a_volume += my_vesseltree.atree.length_of_arc ( i ) * my_vesseltree.atree.arcs[i].get_cross_sect_area();
  }
  for ( unsigned int i = 0; i < my_vesseltree.vtree.arcs.size(); i++ ) {
    v_volume += my_vesseltree.vtree.length_of_arc ( i ) * my_vesseltree.vtree.arcs[i].get_cross_sect_area();
  }

  if ( a_volume > v_volume ) {
    v_val = -1.0; a_val = 1.0 * v_volume / a_volume;
  } else {
    a_val = 1.0; v_val = -1.0 * a_volume / v_volume;
  }
  cerr << "Volume of atree: " << a_volume << " value on atree: " << a_val << endl;
  cerr << "Volume of vtree: " << v_volume << " value on vtree: " << v_val << endl;

  for ( unsigned int i = 0; i < my_vesseltree.atree.arcs.size(); i++ ) {
    my_vesseltree.atree.arcs[i].set_value ( a_val );
  }
  for ( unsigned int i = 0; i < my_vesseltree.vtree.arcs.size(); i++ ) {
    my_vesseltree.vtree.arcs[i].set_value ( v_val );
  }
}


int main() {
  try {
    aol::ParameterParser params ( "par/heatcond.par" );
    //    double scalefactor =  params.getDouble("SCALEFACTOR");
    double cgepsilon =  params.getDouble ( "CGEPSILON" );
    int cgiter = params.getInt ( "CGITER" );

    vesseltree<double> my_vesseltree;
    char filename[1024];

    char atreefilename[1024], vtreefilename[1024];
    params.getString ( "ATREEFILENAME", atreefilename );
    params.getString ( "VTREEFILENAME", vtreefilename );
    my_vesseltree.read_trees ( atreefilename, vtreefilename );

    char datafilename[1024]; // picturefilename[1024], picturefilename_color[1024];
    //     params.getString( "PICTUREFILENAME", picturefilename );
    //     params.getString( "COLORPICFN", picturefilename_color );
    params.getString ( "DATAFILENAMES", datafilename );

    set_appropriate_source_values ( my_vesseltree );

    int lev = params.getInt ( "LEVEL" );

    qc::GridDefinition mygrid ( lev, qc::QC_2D );

    tissue_vessel_geometry<double>
    tvgeom ( atreefilename, vtreefilename, lev ),
    tvgeom_steady ( atreefilename, vtreefilename, lev );    // misuse of tvgeoms to allow saving.

    qc::ScalarArray< double, qc::QC_2D > rhs ( mygrid ),
    u_new ( tvgeom.tissue_data, aol::FLAT_COPY ), // want reference
    u_steady ( tvgeom_steady.tissue_data, aol::FLAT_COPY ),
    dwarf ( mygrid );

    //    rhs.setOverflowHandling(aol::CLIP_THEN_SCALE, 0, 255);

    cerr << "computing RHS ... ";
    my_vesseltree.update_segments_in_grid ( mygrid );
    my_vesseltree.compute_rhs_fphi ( mygrid, rhs );
    cerr << "max value: " << rhs.getMaxValue()  << ", min value: " << rhs.getMaxValue() <<  endl;
    project_ ( rhs ); // why that?

    rhs *= params.getDouble ( "rhs_scale" ); // scale intensity of source terms

    cerr << "done. max value: " << rhs.getMaxValue()  << ", min value: " << rhs.getMaxValue() <<  endl;

    typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature< double, qc::QC_2D, 3 > > ConfigType;
    const aol::StiffOp<ConfigType>  L ( mygrid, aol::ASSEMBLED );
    const aol::MassOp<ConfigType>   M ( mygrid, aol::ASSEMBLED );
    aol::LinCombOp< aol::Vector<double>, aol::Vector<double> > kL;
    kL.appendReference ( L, params.getDouble ( "kappa" ) );

    //    u_new.setZero();
    dwarf.setZero();

#ifdef CONV_ANALY
    cerr << "Computing steady state solution" << endl;
    u_steady.setZero();

    aol::LagrangeZeroProjectSolver< double > Lin ( kL, mygrid.getNumberOfNodes(), cgepsilon, cgiter, false );
    Lin.apply ( rhs, u_steady );

    //    my_scaled_save(u_steady, "output/steady_state.pgm", "output/steady_state.ppm", scalefactor);
    sprintf ( filename, "output/steady_state.dat" );
    tvgeom_steady.save_all_data ( filename, 10 );
    double init_error = u_steady.norm();
    cerr << "l2 norm of steady state solution " << init_error << " on grid level " << mygrid.getGridDepth() <<  endl;
#endif

    double tau = params.getDouble ( "TAUOVERH" ) * mygrid.H();
    unsigned int numtimesteps = static_cast< unsigned int > ( params.getInt ( "NUMTIMESTEPS" ) );

    u_new.setZero();

    aol::LinCombOp< aol::Vector <double>, aol::Vector <double> > A;
    A.appendReference ( M ); A.appendReference ( L, params.getDouble ( "kappa" ) * tau );
    aol::CGInverse< aol::Vector <double> > Ain ( A, cgepsilon, cgiter );  // can use PCG with diagonal preconditioning here.
    Ain.setStopping( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    Ain.setQuietMode ( true );

    qc::ScalarArray< double, qc::QC_2D > tau_x_rhs ( rhs, aol::DEEP_COPY ); // deep copy
    tau_x_rhs *= tau;

    cerr << "sum of tau_x_rhs: " <<  tau_x_rhs.sum() << endl;

    double error = 0;

    for ( unsigned int t = 0; t < numtimesteps; t++ ) {
      cerr << "Computing time step " << t << endl;
      dwarf.setZero();

      //      cerr << "Sum u_new " << u_new.sum() << endl;

      M.apply ( u_new, dwarf );

      //      cerr << "Sum dwarf " << dwarf.sum() << endl;

      dwarf += tau_x_rhs;

      //      cerr << "Sum tau_x_rhs " << tau_x_rhs.sum() << endl;

      //      cerr << "Sum dwarf " << dwarf.sum() << endl; // this still has zero sum.

      u_new.setZero();
      Ain.apply ( dwarf, u_new );

      //      cerr << "Sum u_new " << u_new.sum() << endl; // at this point, (1, ..., 1) u_new is no longer zero.

      sprintf ( filename, datafilename, t );
      tvgeom.save_all_data ( filename, 10 );
      //       sprintf(tmp, picturefilename, t);
      //       sprintf(tmpc, picturefilename_color, t);
      //       //      my_scaled_save(u_new, tmp, tmpc, scalefactor);

      dwarf = u_new;
      dwarf -= u_steady;
      error = dwarf.norm();

#ifdef VERBOSE
      double sum = u_new.sum();
      cerr << "sum, sum/h^2 of u_new: " <<  sum << ", " << sum * mygrid.H() * mygrid.H() << endl;
#ifdef CONV_ANALY

      // force to zero. this is cheated ...
      //      u_new += (- sum / u_new.size() );

      cerr << "l2 error between current timestep and steady state solution: " << error * mygrid.H() * mygrid.H()
      << ", relative to init error: " << error / init_error << endl; // h^2 = length of image vector.
#endif
#endif
#ifdef CONV_ANALY
      cout << log ( error / init_error ) << endl;
      //      sprintf(tmp, "output/difference%04d.pgm", t);
      //      sprintf(tmpc, "output/difference%04d.ppm", t);
      //      my_scaled_save(dwarf, tmp, tmpc, scalefactor);
#endif
    }
#ifdef CONV_ANALY
#ifdef VERBOSE
    cout << "lev = " << lev << ", tau = " << tau << ", final relative error = " << error / init_error << endl;
#endif
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return ( 0 );
}
