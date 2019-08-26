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
 * Plot data from advectheat calculation
 * tissue data to ppm file, vessel data to eps file (2D topview, 3D)
 * ******************************************************************* */

#include "tissue_vessel_geometry.h"
// #include "usevesseltree_ellamadvect.h" // for now
#include "parameterParser.h"

int main ( int argc, char** argv ) {
  if ( argc == 6 || argc == 5 ) {
    try {

      aol::ParameterParser params ( argv[1] );
      char atree_fn[512], vtree_fn[512];
      params.getString ( "atree_file", atree_fn );
      params.getString ( "vtree_file", vtree_fn );

      const double
      scalefactor = params.getDouble ( "scalefactor" ),
                    color_scale = params.getDouble ( "color_scale" );

      const int
      level_tissue = params.getInt ( "leveltissue" ),
                     level_vessels = params.getInt ( "levelvessels" );

      //      system("sleep 10");

      cerr << "Setting up tvgeom ... " << endl;
      tissue_vessel_geometry<double> tvgeom ( atree_fn, vtree_fn, level_tissue, level_vessels );
      tvgeom.read_all_data ( argv[2] );

      //      system("sleep 20");

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

      if ( argc == 6 ) {
        if ( !static_cast<bool> ( params.getInt ( "topview" ) ) ) {
          tvgeom.save_advect_step ( argv[3], colTrans, zzoom, color_scale, rot, shift, true, params.getDouble ( "min_plot_speed" ) );
        } else {
          zzoom *= 1.0;
          tvgeom.save_advect_step_topview_recursive ( argv[3], colTrans, true, params.getDouble ( "min_plot_speed" ), params.getDouble ( "ves_avg_val" ), params.getDouble ( "ves_scalefactor" ) ); // true: plot temperatures
        }
      }

      cerr << "max = " << tvgeom.tissue_data.getMaxValue()
      << ", min = " << tvgeom.tissue_data.getMinValue()
      << ", avg = " << tvgeom.tissue_data.sum() / tvgeom.tissue.getNumberOfNodes() << endl;

      if ( argc == 6 )
        img_scaled_save_ctrscl ( tvgeom.tissue_data, argv[4], argv[5], params.getDouble ( "avg_val" ), scalefactor, colTrans );
      if ( argc == 5 )
        img_scaled_save_ctrscl ( tvgeom.tissue_data, argv[3], argv[4], params.getDouble ( "avg_val" ), scalefactor, colTrans );

    } catch ( aol::Exception &ex ) {
      ex.dump();
    }

  } else {
    cerr << "usage: usevesseltree_advect_and_heat <parameter file> <data file> [ <tree temp profile file> ] <tissue grey file> <tissue color file>" << endl;
  }
  return ( 0 );
}

