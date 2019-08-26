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

#include "parameterParser.h"
#include <tpCFEUtils.h>

static const bool OUTPUTVIRTUALNODEBALLS = false;

typedef double RealType;

const RealType NONPENETRATION_OFFSET = 0.0;

int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: " << argv[0] << " filename_root" << endl;
    abort();
  }
  try {

    char infilename[1024];
    sprintf ( infilename, "%s.dat.bz2", argv[1] );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( infilename );

    typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_CD > GridType;

    typedef     tpcfe::CFEInterfaceTriangulationGenerator<GridType> ITGType;

    ITGType::VertexVectorType vertices;
    ITGType::TriangVectorType triangles;

    GridType grid ( qc::GridSize<qc::QC_3D>::createFrom ( levelset ) );
    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    {
      // attention: waste of memory, copy only used for convenience
      ITGType itg( grid ); // not necessary to load any data
      itg.determineInterfaceAndBoundaryTriangulation( aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( grid.getWidth() - 1, grid.getWidth() - 1, grid.getWidth() - 1 ) );
      itg.getVertexVector( vertices );
      itg.getTriangVector( triangles );
    }

    if ( ( vertices.size() == 0 ) || ( triangles.size() == 0 ) ) {
      cerr << vertices.size() << " " << triangles.size() << endl;
      abort();
    }

    char povfilename[1024], povinsidefilename[1024];
    sprintf ( povfilename, "%s.pov", argv[1] );
    ofstream povstream ( povfilename );

    sprintf ( povinsidefilename, "%s_inside.inc", argv[1] );
    ofstream povinsidestream ( povinsidefilename );

    povstream << "#include \"colors.inc\" " << endl
              << "#include \"textures.inc\" " << endl
              << "#include \"glass.inc\" " << endl
              << endl
              << "light_source {<0, 0, -1000> rgb 0.46}" << endl
              << "light_source {<150, 50, -200> rgb 0.46}" << endl
              << "light_source {<-0.6, 1.6, 3.7>*10000 rgb 1.0 }" << endl
              << "light_source {< 1.5, -1.5, 1.5 > rgb 0.8 spotlight point_at < 0, 1, 0 > }" << endl
              << endl
              << "camera { orthographic" << endl
              << "  angle       60" << endl
              << "  location    < 1.3, -2.0, 0.5 >" << endl
              << "  sky         z" << endl
              << "  up          z" << endl
              << "  look_at     < 0.5, 0.5, 0.0 >" << endl
              << "}" << endl
              << endl
              << "background {" << endl
              << "  color rgb < 1.0, 1.0, 1.0 >" << endl
              << "}" << endl
              << endl
              << "#declare myInside = mesh2 {" << endl
              << "#include \"" << povinsidefilename << "\"" << endl
              << "}" << endl
              << endl
              << "object { myInside" << endl
              << "  texture { pigment { color rgbf < 0.15, 0.15, 0.15, 0.0 > } finish { Dull } } no_shadow" << endl
              << "}" << endl
              << "box { < 0.0001, 0.0001, 0.0001 >, < "
              << 0.9999 * grid.H() * ( levelset.getNumX() - 1 ) << ", "
              << 0.9999 * grid.H() * ( levelset.getNumY() - 1 ) << ", "
              << 0.9999 * grid.H() * ( levelset.getNumZ() - 1 ) << " >" << endl
              << "  texture { pigment { color rgbf < 1.0, 0.9, 0.4, 0.85 > } finish { Dull } }" << endl
              << "}" << endl;


    povinsidestream << "vertex_vectors {" << endl
                    << vertices.size() << endl
                    << "< " << vertices[0][0] << ", " << vertices[0][1] << ", " << vertices[0][2] << ">";

    for ( unsigned int i = 1 /*sic!*/; i < vertices.size(); ++i ) {
      povinsidestream << "," << endl << "< " << vertices[i][0] << ", " << vertices[i][1] << ", " << vertices[i][2] << ">";
    }

    povinsidestream << endl
                    << "}" << endl
                    << endl
                    << "face_indices {" << endl
                    << triangles.size() << endl
                    << "< " << triangles[0][0] << ", " << triangles[0][1] << ", " << triangles[0][2] << ">";

    for ( unsigned int i = 1 /*sic!*/; i < triangles.size(); ++i ) {
      povinsidestream << "," << endl << "< " << triangles[i][0] << ", " << triangles[i][1] << ", " << triangles[i][2] << ">";
    }

    povinsidestream << endl
                    << "}" << endl;

    if ( OUTPUTVIRTUALNODEBALLS ) {
      // output nodes

      char povnodefilename[1024];
      sprintf ( povnodefilename, "%s_vnodes.inc", argv[1] );
      ofstream povnodestream ( povnodefilename );
      for ( GridType::VNMapTypeConstIt vnit = grid.virtualNodeMapRef().begin(); vnit != grid.virtualNodeMapRef().end(); ++vnit ) {
        const RealType xc = grid.H() * vnit->second->_coord[0], yc = grid.H() * vnit->second->_coord[1], zc = grid.H() * vnit->second->_coord[2];
        povnodestream << "sphere { < " << xc << ", " << yc << "," << zc <<" >, ballRad" << endl
                      << "  texture { pigment { color Red } } no_shadow" << endl
                      << "}" << endl;
      }

      povstream << "#declare ballRad = 0.006;" << endl
                << "#include \"" << povnodefilename << "\"" << endl;
    }

  } catch ( aol::Exception e ) {
    e.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}

