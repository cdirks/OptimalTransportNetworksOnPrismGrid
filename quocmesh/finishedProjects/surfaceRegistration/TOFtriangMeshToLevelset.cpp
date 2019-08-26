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

/** \file tools/qcsmLevelset/obj2udply.cpp
 *  \brief converts ply into quoc-levelset function
 *
 *  Converts a ply-file to a quoc-levelset function.
 *
 *  Usage: triangMeshToLevelset sourceFile.ply destinationFile.bz2
 *
 *  \author Wirth, Bauer
 */

#include <aol.h>
#include <gridBase.h>
#include <scalarArray.h>
#include <configurators.h>
#include <triangMeshToDistFunc.h>
#include <triangMesh.h>
#include <eikonalNA3d.h>
#include <queue>
#include <bitArray.h>
#include <deformations.h>

typedef double RealType;

typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 4 ) {
      cerr << "Reads in a surface from a quoc-.ply file and exports it as quoc-.bz2 file.\n"
           << aol::color::red
           << "usage: " << argv[0] << " source-file.ply destination-file.bz2 grid-level [closedFlag scaleFactor offset_x offset_y offset_z]\n"
           << "If closedFlag is set, non-closed surfaces are ignored.\n"
           << "In that case, a local signed distance function (just\n"
           << "defined on the nodes near the interface) is additionally\n"
           << "saved as ""interface_""<destination-filename>.\n"
           << aol::color::reset;
      return EXIT_FAILURE;
    }

    cerr << aol::color::blue << "Read in .ply file ...\n" << aol::color::reset;
    aol::TriangMesh<RealType> mesh;
    mesh.loadFromUDPLY ( argv[1] );
    mesh.makeOrientationConsistent();
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

    cerr << aol::color::blue << "Initialize levelset function.\n" << aol::color::reset;

    qc::GridDefinition grid ( atoi ( argv[3] ), qc::QC_3D );
    qc::GridDefinition enlargedGrid ( atoi ( argv[3] ) + 1, qc::QC_3D );
    ConfiguratorType::ArrayType levelsetFnc ( enlargedGrid ), levelsetFncSign ( enlargedGrid ), levelsetFncAbs ( enlargedGrid );

    cerr << aol::color::blue << "Prepare conversion of mesh into levelset function.\n" << aol::color::reset;
    vector<qc::CoordType> seedPoints;
    qcsm::TriangMeshToDistFunc<RealType> qcdist ( mesh, levelsetFnc, seedPoints, true );

    // shift mesh center to (0,0,0) and scale mesh (0.5x)
    aol::Vec3<RealType> minXYZ, maxXYZ, avXYZ, test;
    mesh.computeBoundingBox ( minXYZ, maxXYZ );

    maxXYZ += minXYZ;
    maxXYZ *= -.5;
    mesh.shiftByOffset ( maxXYZ );

    RealType scaleFactor = 0.5;
    minXYZ.setAll ( scaleFactor );
    mesh.scaleSizeByFactor ( minXYZ );

    maxXYZ[0] = 0.5 + ( -maxXYZ[0] - 0.5 ) * scaleFactor;
    maxXYZ[1] = 0.5 + ( -maxXYZ[1] - 0.5 ) * scaleFactor;
    maxXYZ[2] = 0.5 + ( -maxXYZ[2] - 0.5 ) * scaleFactor;

    mesh.shiftByOffset ( maxXYZ );

    cerr << aol::color::blue << "Create a vector of nodes (""seedpoints"") containing all nodes closest to the mesh ... " << aol::color::reset;
    qcdist.makeDistFkt ( );
    // now, the levelset function is infinite everywhere except at the seed points
    //levelsetFnc.save( "intermediateResult.pgm", qc::PGM_DOUBLE_BINARY );
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

    cerr << aol::color::blue << "Extend the levelset function from the seedpoints to the whole cube ... " << aol::color::reset;// via solution of eikonal equation
    eik::EikonalNA3D<RealType> eikonal ( enlargedGrid );
    int numberOfSeeds = seedPoints.size();
    for ( int i = 0; i < numberOfSeeds ; i++ )
      eikonal.setSeedPoint ( seedPoints[i], aol::Abs ( levelsetFnc.get ( seedPoints[i] ) ) );
    eikonal.march();
    // now, the result is a distance function, which is positive everywhere
    //(eikonal.getTimeField()).save( "intermediateResult.pgm", qc::PGM_DOUBLE_BINARY );
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

    cerr << aol::color::blue << "Do marching cubes to detect all grid nodes outside the shape ... " << aol::color::reset;
    int signOutside = 0;
    levelsetFncSign.setAll ( 1. );
    std::queue<aol::Vec3<int> > activePoints;
    aol::Vec3<int> pointToBePushed ( 0, 0, 0 ), currentPoint;
    activePoints.push ( pointToBePushed );
    int gridWidthMin1 = enlargedGrid.getNumX() - 1;
    while ( !activePoints.empty() ) {
      currentPoint = activePoints.front();
      activePoints.pop();
      for ( int x = aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[0] - 1 ) ); x <= aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[0] + 1 ) ); x++ )
        for ( int y = aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[1] - 1 ) ); y <= aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[1] + 1 ) ); y++ )
          for ( int z = aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[2] - 1 ) ); z <= aol::Max ( 0, aol::Min ( gridWidthMin1, currentPoint[2] + 1 ) ); z++ )
            if ( levelsetFncSign.get ( x, y, z ) > -1. ) {
              levelsetFncSign.set ( x, y, z, -1. );
              if ( aol::isInf ( levelsetFnc.get ( x, y, z ) ) ) {
                //if ( levelsetFnc.get( x, y, z ) > .6 * grid.H() ) {
                pointToBePushed[0] = x;
                pointToBePushed[1] = y;
                pointToBePushed[2] = z;
                activePoints.push ( pointToBePushed );
              } else
                signOutside += aol::signum ( levelsetFnc.get ( x, y, z ) );
            }
    }
    if ( signOutside > 0 )
      levelsetFnc *= -1.;
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

    string outname ( argv[2] ), basename ( outname, 0, outname.length() - 3 );
    basename.append ( "interface.bz2" );
    if ( argc > 4 && atoi ( argv[4] ) )
      levelsetFnc.save ( basename.c_str(), qc::PGM_DOUBLE_BINARY );

    cerr << aol::color::blue << "Make all levelset function values outside the shape and away from the interface negative ... " << aol::color::reset;
    levelsetFncAbs = eikonal.getTimeField();
    for ( int i = 0; i < levelsetFncAbs.size(); i++ )
      if ( aol::isInf ( levelsetFnc[i] ) )
        levelsetFncAbs[i] *= levelsetFncSign[i];
      else
        levelsetFncAbs[i] = levelsetFnc[i];
    levelsetFnc = levelsetFncAbs;
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

    if ( argc > 4 && atoi ( argv[4] ) ) {
      cerr << aol::color::blue << "Make levelset function values negative on both sides of a non-closed surface ... " << aol::color::reset;
      for ( int i = 0; i < 2; i++ ) // two runs
        for ( int x = 0; x < enlargedGrid.getNumX() - 1; x++ )
          for ( int y = 0; y < enlargedGrid.getNumX() - 1; y++ )
            for ( int z = 0; z < enlargedGrid.getNumX() - 1; z++ ) {
              qc::CoordType pos ( x, y, z );
              for ( int j = 0; j < qc::QC_3D; j++ ) {
                qc::CoordType posFront ( pos ), posBack ( pos );
                posFront[j] += 1;
                posBack[j] -= 1;
                if ( aol::Abs ( levelsetFnc.get ( pos ) - levelsetFnc.get ( posFront ) ) > enlargedGrid.H() ) {
                  // this can only happen right at an interface, where nodes on one side are positive and negative on the other side,
                  // and where due to non-closedness of the surface, the positive nodes are surrounded by negative ones;
                  // hence there must be a double sign change at pos or at posFront
                  if ( levelsetFnc.get ( posBack ) * levelsetFnc.get ( posFront ) > 0. )
                    levelsetFnc.set ( pos, - levelsetFnc.get ( pos ) );
                  else
                    levelsetFnc.set ( posFront, - levelsetFnc.get ( posFront ) );
                }
              }
            }
      cerr << aol::color::blue << "done.\n" << aol::color::reset;
    }

    cerr << aol::color::blue << "Save the result ... " << aol::color::reset;

    ConfiguratorType::ArrayType levelsetFncRescaled ( grid );
    levelsetFncRescaled.padFrom ( levelsetFnc );

    levelsetFncRescaled.save ( argv[2], qc::PGM_DOUBLE_BINARY );
    cerr << aol::color::blue << "done.\n" << aol::color::reset;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
