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

#include <multiArray.h>
#include <parameterParser.h>
#include <shapeLevelsetGenerator.h>
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include "computeForce.h"


typedef double RealType;
static const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;
typedef tpcfe::CFEGrid < RealType, CT, tpcfe::VoigtElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;



const bool simulate = true;
const bool visualize = true;

int main ( const int, const char**  ) {

  const int level = 3, nRods = 1;

  GridType grid ( level  );
  grid.setCoarsenedFlag ( true ); // this is a workaround for storing the full matrix (necessary due to non-constant coefficient)

  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
  qc::ShapeLevelsetGenerator<RealType>::generateXRotated3DRodsLevelset ( levelset, nRods, aol::Vec3<RealType>( 0.12 / nRods, 0.12 / nRods, 0.12 / nRods ), atan ( 1. / nRods ) );
  grid.addStructureFrom ( levelset );

  levelset.save ( "out/jceTransIsoXRotRods.dat.bz2", qc::PGM_DOUBLE_BINARY );

  qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );

  {
    aol::Mat<6,6,RealType> VoigtTensorMinus[3], VoigtTensorPlus;
    { // plus is the surrounding material
      RealType
        Eone = 3.0,
        nu = 0.38,
        Gmix = Eone / ( 2 * ( 1 + nu ) );

      aol::Mat<6,6,RealType> VoigtTensorInverse;

      VoigtTensorInverse[0][0] = 1.  / Eone;   VoigtTensorInverse[0][1] = -nu / Eone;   VoigtTensorInverse[0][2] = -nu / Eone;
      VoigtTensorInverse[1][0] = -nu / Eone;   VoigtTensorInverse[1][1] = 1   / Eone;   VoigtTensorInverse[1][2] = -nu / Eone;
      VoigtTensorInverse[2][0] = -nu / Eone;   VoigtTensorInverse[2][1] = -nu / Eone;   VoigtTensorInverse[2][2] = 1   / Eone;
      VoigtTensorInverse[3][3] = 1.  / Gmix;
      VoigtTensorInverse[4][4] = 1.  / Gmix;
      VoigtTensorInverse[5][5] = 1.  / Gmix;

      VoigtTensorPlus.makeInverse ( VoigtTensorInverse );

      cerr << "VoigtTensorPlus" << endl << VoigtTensorPlus << endl << endl;

    }

    { // minus are the trabeculae
      VoigtTensorMinus[0][0][0] = 19.7785;   VoigtTensorMinus[0][0][1] =  8.4190;   VoigtTensorMinus[0][0][2] =  8.4190;
      VoigtTensorMinus[0][1][0] =  8.4190;   VoigtTensorMinus[0][1][1] = 16.1824;   VoigtTensorMinus[0][1][2] =  7.6152;
      VoigtTensorMinus[0][2][0] =  8.4190;   VoigtTensorMinus[0][2][1] =  7.6152;   VoigtTensorMinus[0][2][2] = 16.1824;
      VoigtTensorMinus[0][3][3] =  8.5671;
      VoigtTensorMinus[0][4][4] =  9.4713;
      VoigtTensorMinus[0][5][5] =  9.4713;

      VoigtTensorMinus[1][0][0] = 16.1824;   VoigtTensorMinus[1][0][1] =  8.4190;   VoigtTensorMinus[1][0][2] =  7.6152;
      VoigtTensorMinus[1][1][0] =  8.4190;   VoigtTensorMinus[1][1][1] = 19.7785;   VoigtTensorMinus[1][1][2] =  8.4190;
      VoigtTensorMinus[1][2][0] =  7.6152;   VoigtTensorMinus[1][2][1] =  8.4190;   VoigtTensorMinus[1][2][2] = 16.1824;
      VoigtTensorMinus[1][3][3] =  9.4713;
      VoigtTensorMinus[1][4][4] =  8.5671;
      VoigtTensorMinus[1][5][5] =  9.4713;

      VoigtTensorMinus[2][0][0] = 16.1824;   VoigtTensorMinus[2][0][1] =  7.6152;   VoigtTensorMinus[2][0][2] =  8.4190;
      VoigtTensorMinus[2][1][0] =  7.6152;   VoigtTensorMinus[2][1][1] = 16.1824;   VoigtTensorMinus[2][1][2] =  8.4190;
      VoigtTensorMinus[2][2][0] =  8.4190;   VoigtTensorMinus[2][2][1] =  8.4190;   VoigtTensorMinus[2][2][2] = 19.7785;
      VoigtTensorMinus[2][3][3] =  9.4713;
      VoigtTensorMinus[2][4][4] =  9.4713;
      VoigtTensorMinus[2][5][5] =  8.5671;

      cerr << "VoigtTensorsMinus" << endl << VoigtTensorMinus[0] << endl << endl << VoigtTensorMinus[1] << endl << endl << VoigtTensorMinus[2] << endl << endl;

    }

    tpcfe::setVaryingTransverselyIsotropicTensorCoeffs<GridType> ( grid, levelset, nRods, VoigtTensorMinus, VoigtTensorPlus, coeff, atan ( 1. / nRods ) );
  }

  grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-12, 1.0e-12 );

  if ( simulate == true ) {

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    // tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, qc::QC_Z, qc::QC_X, 0.01 );
    tpcfe::setTorsion ( grid, dirichletBCs, DirichletMask, aol::NumberTrait<RealType>::pi / 180 );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    // tpcfe::CFEJCEMassOp< ConfiguratorType > massOp ( grid ); // not needed
    ElastOpType elastOp ( grid, coeff );

    dwarf -= dirichletBCs;

    elastOp.apply ( dwarf, rhs );

    elastOp.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    aol::BlockGaussSeidelPreconditioner< aol::MultiVector<RealType>, ElastOpType::BlockMatrixType > prec ( elastOp.getBlockMatrixRef() );
    prec.setSweepingModes ( aol::GAUSS_SEIDEL_ZEBRA2_FORWARD, aol::GAUSS_SEIDEL_ZEBRA2_BACKWARD );

    aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, 1.0e-16, 100000 );

    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    solver.apply ( rhs, soln );

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

    soln += dirichletBCs;

    soln.save ( "out/jceTransIsoXRotRods_Def%d.dat.bz2", qc::PGM_DOUBLE_BINARY );

  }

  if ( visualize == true ) {

    qc::MultiArray< RealType, qc::QC_3D > soln ( grid );
    soln.load ( "out/jceTransIsoXRotRods_Def%d.dat.bz2" );

    const RealType displayFactor = 20;

    aol::TriangMesh<float> tMesh;

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_X );
      itg.determineSliceTriangulation ( qc::QC_Y, 5 * grid.getNumY() / 6 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsSliceDefX.ply" );
    }
    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Y );
      itg.determineSliceTriangulation ( qc::QC_Y, 5 * grid.getNumY() / 6 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsSliceDefY.ply" );
    }
    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineSliceTriangulation ( qc::QC_Y, 5 * grid.getNumY() / 6 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsSliceDefZ.ply" );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineInterfaceTriangulation();
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoXRotRodsInterfaceU.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsInterfaceU.ply" );

      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoXRotRodsInterfaceD.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsInterfaceD.ply" );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itgCube ( grid, soln, qc::QC_Z );
      itgCube.determineBoundaryTriangulation( +1 );
      itgCube.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoXRotRodsCubeU.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsCubeU.ply" );

      itgCube.deformVertices ( soln, displayFactor );
      itgCube.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoXRotRodsCubeD.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoXRotRodsCubeD.ply" );
    }

  }

  return ( EXIT_SUCCESS );
}
