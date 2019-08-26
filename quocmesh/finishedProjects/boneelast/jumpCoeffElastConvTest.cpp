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
typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;


int main ( const int argc, const char** argv ) {
  if ( argc != 3 ) {
    cerr << "Usage: " << argv[0] << " {S,C,F} level" << endl;
    abort();
  }

  const RealType
    EMinus  = 2.0,
    nuMinus = 0.1,
    EPlus   = 1.0,
    nuPlus  = 0.3;

  cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;

  const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );

  const int GStdLevel = 8;
  const char* sfnmask = "out/JCEDeformation%%d_Lev%d.dat.bz2";
  const RealType radius = 1./3.;

  const RealType initVNStart = 1.0e-14, initVNStep = 1.0e-14;

  char solnfilenamemask[1024];

  aol::StopWatch timer;

  if ( argv[1][0] == 'S' ) {

    const int level = atoi ( argv[2] );

    GridType grid ( level );

    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset, radius );

    grid.addStructureFrom ( levelset );


    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( grid );

    tpcfe::setGeneralShearing ( grid, dirichletBCs, DirichletMask, 2, 2, 0.01 );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    tpcfe::CFEJCEMassOp< ConfiguratorType > massOp ( grid );
    ElastOpType elastOp ( grid, coeff );

    dwarf -= dirichletBCs;

    elastOp.apply ( dwarf, rhs );

    elastOp.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    timer.start();

#if 0
    aol::DiagonalPreconditioner < aol::MultiVector<RealType> > prec ( elastOp.getBlockMatrixRef() );
    // aol::BlockDiagonalPreconditioner < RealType, qc::QC_3D > prec ( elastOp.getBlockMatrixRef() ); // convergence observed to be worse than DiagonalPreconditioner
#else
    tpcfe::CFEBlockMultigridPreconditioner< ElastOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > prec ( grid, elastOp.getBlockMatrixRef(), 1, 3, 3 );
#endif
    aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, 1.0e-16, 10000 );

    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    // mg::boomerAMGBlockSolver< tpcfe::CFEJCElastOp< ConfiguratorType >::BlockMatrixType > solver ( elastOp.getBlockMatrixRef() ); // poor convergence, hardly surprising ...

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MB" << endl;

    solver.apply ( rhs, soln );

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MB" << endl;

    timer.stop();
    cerr << "Solver (including setup) took " << timer.elapsedWallClockTime() << " seconds" << endl;

    soln += dirichletBCs;

    sprintf ( solnfilenamemask, sfnmask, level );
    soln.save ( solnfilenamemask, qc::PGM_DOUBLE_BINARY );

  }

  if ( argv[1][0] == 'C' ) {

    GridType GStdGrid ( GStdLevel );
    qc::ScalarArray<RealType, qc::QC_3D> GStdLevelset ( GStdGrid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( GStdLevelset, radius );
    GStdGrid.addStructureFrom ( GStdLevelset );

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > GStdCoeff ( GStdGrid );
    tpcfe::setCoeffForLevelset ( GStdCoeff, GStdLevelset, ENuMinus, ENuPlus );
    GStdGrid.relaxedDetectAndInitVirtualNodes ( GStdCoeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > GStdSoln ( GStdGrid );
    sprintf ( solnfilenamemask, sfnmask, GStdLevel );
    GStdSoln.load ( solnfilenamemask );

    const int compLevel = atoi ( argv[2] );

    GridType compGrid ( compLevel );
    qc::ScalarArray<RealType, qc::QC_3D> compLevelset ( compGrid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( compLevelset, radius );
    compGrid.addStructureFrom ( compLevelset );

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > compCoeff ( compGrid );
    tpcfe::setCoeffForLevelset ( compCoeff, compLevelset, ENuMinus, ENuPlus );
    compGrid.relaxedDetectAndInitVirtualNodes ( compCoeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > compSoln ( compGrid );
    sprintf ( solnfilenamemask, sfnmask, compLevel );
    compSoln.load ( solnfilenamemask );

    // start comparison

    RealType lInfDiff = 0, l1Diff = 0, l2DiffS = 0, l2Sum = 0;

    aol::ProgressBar<> pb ( "Computing Difference" );
    pb.start ( GStdGrid.getNumberOfElements() );
    const RealType h = GStdGrid.H();

    const qc::GridSize<qc::QC_3D> gridSize ( GStdGrid );
    for ( GridType::FullElementIterator it ( GStdGrid ); it.notAtEnd(); ++it ) {
      tpcfe::CFEElement<RealType> el ( *it, gridSize, GStdGrid.getElType ( *it ) );

      el.computeAssembleData ( GStdGrid );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
        const tpcfe::CFETetra<RealType> &tetra = *tit;

        aol::Vec3<RealType> barycenter;
        tetra.computeGlobalCoordinateOfBarycenter ( barycenter, el );
        barycenter *= h; // global -> world coordinates; barycenter of current tetrahedron     // use better quadrature??

        aol::Vec<4, RealType> bary;  bary.setAll( 0.25 );

        aol::Vec3<RealType>
          GStdValue = interpolateDataAtPositionBary ( GStdGrid, GStdSoln, el, tetra, bary ),
          compValue = interpolateDataAtPositionWC ( compGrid, compSoln, barycenter );

        const RealType
          vol = tetra.getVolume() * aol::Cub ( h ),
          diffNorm = ( GStdValue - compValue ).norm();

        lInfDiff = aol::Max ( lInfDiff, diffNorm );
        l1Diff += vol * diffNorm;
        l2DiffS += vol * aol::Sqr ( diffNorm );

        l2Sum += vol * GStdValue.norm();

      }
    }
    pb.finish();

    const RealType
      L2Norm = sqrt ( l2Sum ),
      relLInfDiff = lInfDiff / L2Norm,
      relL1Diff = l1Diff / L2Norm,
      relL2Diff = sqrt ( l2DiffS ) / L2Norm;

    cerr << "l2Norm       = " << aol::longScientificFormat ( L2Norm      ) << endl
         << "rel lInfDiff = " << aol::longScientificFormat ( relLInfDiff ) << endl
         << "rel l1Diff   = " << aol::longScientificFormat ( relL1Diff   ) << endl
         << "rel l2Diff   = " << aol::longScientificFormat ( relL2Diff   ) << endl ;

  }

  if ( argv[1][0] == 'F' ) {

    const RealType aleph = 1;

    const int level = atoi ( argv[2] );

    GridType grid ( level );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset, radius );
    grid.addStructureFrom ( levelset );

    qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > coeff ( grid );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > soln ( grid );
    sprintf ( solnfilenamemask, sfnmask, level );
    soln.load ( solnfilenamemask );

    const unsigned int surfaceDirection = 2;
    const RealType tolerance = 1.0e-8;

    aol::Vec3<RealType> surfaceOuterNormal;
    surfaceOuterNormal [ surfaceDirection ] = 1.0;

    const unsigned int surfacePosition = -1;

    aol::Vector<double> dummy;

    tpcfe::computeForce ( grid, soln, surfaceDirection, surfaceOuterNormal, tolerance, surfacePosition, aleph, coeff, cerr, dummy );

  }

  if ( argv[1][0] == 'V' ) { // visualize

    const int level = atoi ( argv[2] );

    GridType grid ( level  );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset, radius );
    grid.addStructureFrom ( levelset );

    qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > coeff ( grid );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, initVNStart, initVNStep );

    qc::MultiArray< RealType, qc::QC_3D > soln ( grid );
    sprintf ( solnfilenamemask, sfnmask, level );
    soln.load ( solnfilenamemask );

    tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
    itg.determineSliceTriangulation ( qc::QC_X, grid.getNumX() / 2 );
    itg.deformVertices ( soln, -20 );
    itg.saveToPLYFile ( "out/test.ply" );

  }

  return ( EXIT_SUCCESS );
}
