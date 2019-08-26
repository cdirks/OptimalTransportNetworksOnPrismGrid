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


void setVaryingTransverselyIsotropicTensorCoeffs ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &levelset, qc::AArray< GridType::NodalCoeffType, qc::QC_3D > &coeff ) {

  aol::Mat<6,6,RealType> VoigtTensorMinusIso, VoigtTensorMinusOrtho, VoigtTensorPlus;

  { // plus is the surrounding material
    RealType
      Eone = 1.0,
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
  }

  { // minus are the trabeculae
    RealType
      Eone = 10.0,
      nu = 0.33,
      Gmix = Eone / ( 2 * ( 1 + nu ) );

    aol::Mat<6,6,RealType> VoigtTensorInverse;

    VoigtTensorInverse[0][0] = 1.  / Eone;   VoigtTensorInverse[0][1] = -nu / Eone;   VoigtTensorInverse[0][2] = -nu / Eone;
    VoigtTensorInverse[1][0] = -nu / Eone;   VoigtTensorInverse[1][1] = 1   / Eone;   VoigtTensorInverse[1][2] = -nu / Eone;
    VoigtTensorInverse[2][0] = -nu / Eone;   VoigtTensorInverse[2][1] = -nu / Eone;   VoigtTensorInverse[2][2] = 1   / Eone;
    VoigtTensorInverse[3][3] = 1.  / Gmix;
    VoigtTensorInverse[4][4] = 1.  / Gmix;
    VoigtTensorInverse[5][5] = 1.  / Gmix;

    VoigtTensorMinusIso.makeInverse ( VoigtTensorInverse );
  }

  { // minus are the trabeculae
    RealType
      Eone = 10.0,
      Eother = 5.0,
      nu = 0.33,
      Gmix = Eone / ( 2 * ( 1 + nu ) ),
      Gother = Eother / ( 2 * ( 1 + nu ) );

    aol::Mat<6,6,RealType> VoigtTensorInverse;

    VoigtTensorInverse[0][0] = 1.  / Eother;   VoigtTensorInverse[0][1] = -nu / Eother;   VoigtTensorInverse[0][2] = -nu / Eone;
    VoigtTensorInverse[1][0] = -nu / Eother;   VoigtTensorInverse[1][1] = 1   / Eother;   VoigtTensorInverse[1][2] = -nu / Eone;
    VoigtTensorInverse[2][0] = -nu / Eother;   VoigtTensorInverse[2][1] = -nu / Eother;   VoigtTensorInverse[2][2] = 1   / Eone;
    VoigtTensorInverse[3][3] = 1.  / Gmix;
    VoigtTensorInverse[4][4] = 1.  / Gmix;
    VoigtTensorInverse[5][5] = 1.  / Gother;

    VoigtTensorMinusOrtho.makeInverse ( VoigtTensorInverse );
  }

  for ( short i = 0; i < 6; ++i ) {
    for ( short j = 0; j < 6; ++j ) {
      cout << aol::longScientificFormat ( VoigtTensorPlus.get ( i, j ) ) << " ";
    }
    cout << endl;
  }
  cout << endl;

  for ( short i = 0; i < 6; ++i ) {
    for ( short j = 0; j < 6; ++j ) {
      cout << aol::longScientificFormat ( VoigtTensorMinusIso.get ( i, j ) ) << " ";
    }
    cout << endl;
  }
  cout << endl;

  for ( short i = 0; i < 6; ++i ) {
    for ( short j = 0; j < 6; ++j ) {
      cout << aol::longScientificFormat ( VoigtTensorMinusOrtho.get ( i, j ) ) << " ";
    }
    cout << endl;
  }
  cout << endl;

  // set interpolated tensor values

  for ( GridType::FullNodeIterator fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
    if ( levelset.get ( *fnit ) < 0 ) {
      aol::Mat<6,6,RealType> LocalVoigtTensor, dwarf;

      dwarf = VoigtTensorMinusIso;
      dwarf *= 2 * fabs ( grid.H() * (*fnit)[2] - 0.5 );
      LocalVoigtTensor += dwarf;

      dwarf = VoigtTensorMinusOrtho;
      // dwarf = VoigtTensorMinusIso;
      dwarf *= 1 - 2 * fabs ( grid.H() * (*fnit)[2] - 0.5 );
      LocalVoigtTensor += dwarf;

      coeff.getRef ( *fnit ) = GridType::NodalCoeffType( LocalVoigtTensor );

    } else {
      coeff.getRef ( *fnit ) = GridType::NodalCoeffType( VoigtTensorPlus );
    }

  }
}


const bool simulate = true;
const bool visualize = true;

int main ( const int, const char**  ) {

  const int level = 6;

  GridType grid ( level  );
  grid.setCoarsenedFlag ( true ); // this is a workaround for storing the full matrix (necessary due to non-constant coefficient)

  qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
  qc::ShapeLevelsetGenerator<RealType>::generateColumnLevelset ( levelset, 0.12 );
  grid.addStructureFrom ( levelset );

  levelset.save ( "out/jceTransIsoTrab.dat.bz2", qc::PGM_DOUBLE_BINARY );

  qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );

  setVaryingTransverselyIsotropicTensorCoeffs ( grid, levelset, coeff );

  grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-14, 1.0e-14 );

  if ( simulate == true ) {

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

    tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, qc::QC_Z, qc::QC_Z, -0.2 );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    tpcfe::CFEJCEMassOp< ConfiguratorType > massOp ( grid );
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

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MB" << endl;

    soln += dirichletBCs;

    soln.save ( "out/jceTransIsoTrab_Def%d.dat.bz2", qc::PGM_DOUBLE_BINARY );

  }

  if ( visualize == true ) {

    qc::MultiArray< RealType, qc::QC_3D > soln ( grid );
    soln.load ( "out/jceTransIsoTrab_Def%d.dat.bz2" );

    const RealType displayFactor = 1;

    aol::TriangMesh<float> tMesh;

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_X );
      itg.determineSliceTriangulation ( qc::QC_Y, grid.getNumY() / 2 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabSliceDefX.ply" );
    }
    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Y );
      itg.determineSliceTriangulation ( qc::QC_Y, grid.getNumY() / 2 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabSliceDefY.ply" );
    }
    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineSliceTriangulation ( qc::QC_Y, grid.getNumY() / 2 );
      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabSliceDefZ.ply" );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineInterfaceTriangulation();
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoTrabInterfaceU.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabInterfaceU.ply" );

      itg.deformVertices ( soln, displayFactor );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoTrabInterfaceD.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabInterfaceD.ply" );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itgCube ( grid, soln, qc::QC_Z );
      itgCube.determineBoundaryTriangulation( +1 );
      itgCube.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoTrabCubeU.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabCubeU.ply" );

      itgCube.deformVertices ( soln, displayFactor );
      itgCube.writeToTriangMesh ( tMesh );
      tMesh.saveAsPov ( "out/jceTransIsoTrabCubeD.pov" );
      tMesh.saveAsUDPLY ( "out/jceTransIsoTrabCubeD.ply" );
    }

  }

  return ( EXIT_SUCCESS );
}
