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

/* this program loads a levelset file for a tpcfe geometry, applies
   precomputed shifts (here called deformation) to it and saves the
   deformed geometry as a ply file
   Maybe go back two versions where it was used with parameter file?
   \author Schwen
*/


#include "parameterParser.h"
#include <tpCFEUtils.h>


typedef double RealType;

int main ( int argc, char** argv) {

  if ( argc != 2 )
    cerr << "Usage: " << argv[0] << " filename_root" << endl;

  try {

    char filename[1024];

    static const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;
    typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;

    const RealType
      EMinus  = 70.0,
      nuMinus = 0.35,
      EPlus   = 3.0,
      nuPlus  = 0.38;

    const aol::Vec3<float> scaleFactor ( 20, 20, 20 );

    const RealType initVNStart = 1.0e-13, initVNStep = 1.0e-13;

    const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );

    sprintf ( filename, "%s.dat.bz2", argv[1] );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( filename );
    GridType grid ( qc::GridSize<qc::QC_3D>::createFrom ( levelset ) );
    grid.addStructureFrom ( levelset );

    qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > coeff ( grid );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > deformations ( grid );
    sprintf ( filename, "%s_Def%%d.dat.bz2", argv[1] );
    deformations.load ( filename );

    aol::TriangMesh<float> tMesh;

    cerr << "Writing slices ... ";
    const int yPos = grid.getNumY() / 6;
    bool deformSlices = false;
    const RealType aleph = 1.0;

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, deformations, qc::QC_X );
      itg.determineSliceTriangulation ( qc::QC_Y, yPos );
      if ( deformSlices ) itg.deformVertices ( deformations, scaleFactor );
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sX.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, deformations, qc::QC_Y );
      itg.determineSliceTriangulation ( qc::QC_Y, yPos );
      if ( deformSlices ) itg.deformVertices ( deformations, scaleFactor );
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sY.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, deformations, qc::QC_Z );
      itg.determineSliceTriangulation ( qc::QC_Y, yPos );
      if ( deformSlices ) itg.deformVertices ( deformations, scaleFactor );
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sZ.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithVonMisesStressGenerator<GridType> itg ( grid, coeff, deformations, aleph );
      itg.determineSliceTriangulation ( qc::QC_Y, yPos );
      if ( deformSlices ) itg.deformVertices ( deformations, scaleFactor );
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sVMS.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    cerr << "interfaces ... ";

    {
      tpcfe::CFEInterfaceTriangulationWithVonMisesStressGenerator<GridType> itg ( grid, coeff, deformations, aleph );
      itg.determineInterfaceAndBoundaryTriangulation();
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sInterfaceU.pov", argv[1] );
      tMesh.saveAsPov ( filename );
      sprintf ( filename, "%sInterfaceU.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );

      itg.deformVertices ( deformations, scaleFactor );
      itg.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sInterfaceD.pov", argv[1] );
      tMesh.saveAsPov ( filename );
      sprintf ( filename, "%sInterfaceD.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    cerr << "boundaries ... ";

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itgCube ( grid, deformations, qc::QC_Z );
      itgCube.determineBoundaryTriangulation( +1 );
      itgCube.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sCubeU.pov", argv[1] );
      tMesh.saveAsPov ( filename );
      sprintf ( filename, "%sCubeU.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );

      itgCube.deformVertices ( deformations, scaleFactor );
      itgCube.writeToTriangMesh ( tMesh );
      sprintf ( filename, "%sCubeD.pov", argv[1] );
      tMesh.saveAsPov ( filename );
      sprintf ( filename, "%sCubeD.ply.bz2", argv[1] );
      tMesh.saveAsUDPLY ( filename );
    }

    cerr << "done." << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
