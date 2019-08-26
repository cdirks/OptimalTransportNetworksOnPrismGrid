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

/** \file
 *  \brief converts ply into quoc-levelset function
 *
 *  Converts a ply-file to a quoc-levelset function.
 *
 *  Usage: triangMeshToLevelset sourceFile.ply destinationFile.bz2
 *
 *  \author Wirth
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
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    aol::StopWatch watch;
    watch.start();

    cerr << "Read in .ply file ...\n";
    aol::TriangMesh<RealType> mesh;
    mesh.loadBasedOnSuffix( argv[1] );
    mesh.makeOrientationConsistent();
    cerr << "done.\n";

    cerr << "Initialize levelset function.\n";
    qc::GridDefinition grid ( atoi ( argv[3] ), qc::QC_3D );
    ConfiguratorType::ArrayType levelsetFnc ( grid ), levelsetFncSign ( grid ), levelsetFncAbs ( grid );

    cerr << "Prepare conversion of mesh into levelset function.\n";
    vector<qc::CoordType> seedPoints;
    qcsm::TriangMeshToDistFunc<RealType> qcdist( mesh, levelsetFnc, seedPoints, true );

    cerr << "Scale mesh " << ( argc < 6 ? "to fit into the cube [0,1]^d ... " : " " );
    if( argc >= 6 )
      cerr << "by " << atof( argv[5] ) << endl;    
    aol::Vec3<RealType> minXYZ, maxXYZ, avXYZ;
    // shift center of bounding box to origin
    RealType maxWidth = mesh.computeBoundingBox( minXYZ, maxXYZ );
    maxXYZ += minXYZ;
    maxXYZ *= -.5;
    mesh.shiftByOffset( maxXYZ );
    // scale by prescribed factor and check if scaling factor has valid size
    RealType scaleFactor = argc < 6 ? ( 1. - 10. * grid.H() ) / maxWidth : atof( argv[5] );
    if( scaleFactor*maxWidth > 1. ){
        cerr << argv[1] <<  " has bounding box [ " <<  minXYZ[0] << ", " << maxXYZ[0] << " ] x [ "  <<  minXYZ[1] << ", " << maxXYZ[1] << " ] x [ "<<  minXYZ[2] << ", " << maxXYZ[2] << " ] " << endl;
	throw aol::Exception ( "tools/triangManipConv/triangMeshToLevelset.cpp: scaled mesh will not be contained in [0,1]^3", __FILE__, __LINE__ );
    }      
    minXYZ.setAll( scaleFactor );
    mesh.scaleSizeByFactor( minXYZ );
    // shift by (0.5, 0.5, 0.5), ie mesh should be contained in [0,1]^3 now
    minXYZ.setAll( .5 );
    mesh.shiftByOffset( minXYZ );
    cerr << "done.\n";

    if (argc > 6) {
      aol::Vec3<RealType> offset;
      offset[0] = atof ( argv[6] );
      offset[1] = atof ( argv[7] );
      offset[2] = atof ( argv[8] );
      cerr << "Shift by offset [" << offset << "] ... ";
      mesh.shiftByOffset ( offset );
      cerr << "done.\n";
    }

    cerr << "Create a vector of nodes (""seedpoints"") containing all nodes closest to the mesh ... ";
    qcdist.makeDistFkt ( );
    // now, the levelset function is infinite everywhere except at the seed points
    //levelsetFnc.save( "intermediateResult.pgm", qc::PGM_DOUBLE_BINARY );
    cerr << "done.\n";

    cerr << "Extend the levelset function from the seedpoints to the whole cube ... ";// via solution of eikonal equation
    eik::EikonalNA3D<RealType> eikonal ( grid );
    int numberOfSeeds = seedPoints.size();
    for ( int i = 0; i < numberOfSeeds ; i++ )
      eikonal.setSeedPoint ( seedPoints[i], aol::Abs ( levelsetFnc.get ( seedPoints[i] ) ) );
    eikonal.march();
    // now, the result is a distance function, which is positive everywhere
    //(eikonal.getTimeField()).save( "intermediateResult.pgm", qc::PGM_DOUBLE_BINARY );
    cerr << "done.\n";

    cerr << "Do marching cubes to detect all grid nodes outside the shape ... ";
    int signOutside = 0;
    levelsetFncSign.setAll ( 1. );
    std::queue<aol::Vec3<int> > activePoints;
    aol::Vec3<int> pointToBePushed ( 0, 0, 0 ), currentPoint;
    activePoints.push ( pointToBePushed );
    int gridWidthMin1 = grid.getNumX() - 1;
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
    cerr << "done.\n";

    string outname( argv[2] ), basename( outname, 0, outname.length() - 3 );
    basename.append( "interface.bz2" );
    if ( argc > 4 && atoi( argv[4] ) )
      levelsetFnc.save( basename.c_str(), qc::PGM_DOUBLE_BINARY );

    cerr << "Make all levelset function values outside the shape and away from the interface negative ... ";
    levelsetFncAbs = eikonal.getTimeField();
    for ( int i = 0; i < levelsetFncAbs.size(); i++ )
      if ( aol::isInf ( levelsetFnc[i] ) )
        levelsetFncAbs[i] *= levelsetFncSign[i];
      else
        levelsetFncAbs[i] = levelsetFnc[i];
    levelsetFnc = levelsetFncAbs;
    cerr << "done.\n";

    if ( argc > 4 && atoi( argv[4] ) ) {
      cerr << "Make levelset function values negative on both sides of a non-closed surface ... ";
      for ( int i = 0; i < 2; i++ ) // two runs
        for ( int x = 0; x < grid.getNumX() - 1; x++ )
          for ( int y = 0; y < grid.getNumX() - 1; y++ )
            for ( int z = 0; z < grid.getNumX() - 1; z++ ) {
              qc::CoordType pos( x, y, z );
              for ( int j = 0; j < qc::QC_3D; j++ ) {
                qc::CoordType posFront( pos ), posBack( pos );
                posFront[j] += 1;
                posBack[j] -= 1;
                if ( aol::Abs( levelsetFnc.get( pos ) - levelsetFnc.get( posFront ) ) > grid.H() ) {
                  // this can only happen right at an interface, where nodes on one side are positive and negative on the other side,
                  // and where due to non-closedness of the surface, the positive nodes are surrounded by negative ones;
                  // hence there must be a double sign change at pos or at posFront
                  if ( levelsetFnc.get( posBack ) * levelsetFnc.get( posFront ) > 0. )
                    levelsetFnc.set( pos, - levelsetFnc.get( pos ) );
                  else
                    levelsetFnc.set( posFront, - levelsetFnc.get( posFront ) );
                }
              }
            }
      cerr << "done.\n";
    }

    cerr << "Save the result ... ";
    levelsetFnc.save ( argv[2], qc::PGM_DOUBLE_BINARY );
    cerr << "done.\n";

    watch.stop();
    watch.printReport( cerr );
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
}
