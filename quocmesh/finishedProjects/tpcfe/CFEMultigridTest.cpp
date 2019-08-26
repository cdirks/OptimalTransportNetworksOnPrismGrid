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

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include <scalarArray.h>
#include <shapeLevelsetGenerator.h>

typedef double RealType;


void testFailure ( const bool failed, const char *errorText ) {
  if ( !failed ) {
    cerr << aol::color::green << "OK." << aol::color::reset << endl;
  } else {
    cerr << aol::color::red << errorText << endl << aol::color::reset;
  }
}


bool testComplicatedDomain ( const RealType prolongTolerance, const RealType massMatSumTolerance, const RealType coarsenedCompareTolerance, const RealType solutionOrigTolerance, const bool useDirichlet ) {
  bool cdOK = true;
  const RealType thres = 2.0 / 3.0;

  typedef tpcfe::CFEGrid<RealType, tpcfe::CFE_CD>                  GridType;
  typedef tpcfe::CFEHybridMatrix<GridType>                       MatrixType;
  typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;

  GridType fineGrid ( 4 );

  qc::ScalarArray<RealType, qc::QC_3D> fineLevelset ( fineGrid );  qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( fineLevelset, thres );
  fineGrid.setDomainFrom ( fineLevelset );                         fineGrid.detectAndInitVirtualNodes();

  qc::BitArray<qc::QC_3D> fineDirichlet ( qc::GridSize<qc::QC_3D>::createFrom ( fineGrid ) );
  if ( useDirichlet ) {
    for ( int i = 0; i < fineGrid.getNumX(); ++i ) {
      for ( int j = 0; j < fineGrid.getNumY(); ++j ) {
        fineDirichlet.set ( i, j, 0, true );
        fineDirichlet.set ( i, j, fineGrid.getNumZ() - 1, true );
      }
    }
  } else {
    fineDirichlet.setAll ( false );
  }
  fineGrid.setDirichletMask ( fineDirichlet );
  fineGrid.setDOFMaskFromDirichletAndDomainNodeMask();

  GridType coarseGrid ( fineGrid.getGridDepth() - 1 );

  qc::ScalarArray<RealType, qc::QC_3D> coarseLevelset ( coarseGrid );  qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( coarseLevelset, thres );
  coarseGrid.setDomainFrom ( coarseLevelset );                         coarseGrid.detectAndInitVirtualNodes();

  tpcfe::coarsenDirichletMask ( fineGrid, coarseGrid ); // does the right thing in non-Dirichlet case
  coarseGrid.setDOFMaskFromDirichletAndDomainNodeMask();

  {
    cerr << "--- Testing CFE (complicated domain) prolongation ... ";
    bool failed = false;

    tpcfe::CFEProlongOp< GridType > prOp ( coarseGrid, fineGrid, aol::ONTHEFLY );
    tpcfe::CFEProlongOp< GridType > prOpA ( coarseGrid, fineGrid, aol::ASSEMBLED );
    failed |= !compareOps ( prOp, prOpA, coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() );

    qc::ScalarArray<RealType, qc::QC_3D> fArr ( fineGrid ), cArr ( coarseGrid );
    cArr.setAll ( 1.0 ); coarseGrid.restrictToDofs ( cArr );
    prOp.apply ( cArr, fArr );
    const RealType expectedResult = ( useDirichlet ? 2856 : 3468 );
    //cerr << endl << fabs ( 0 - fArr.getMinValue() ) + fabs ( 1 - fArr.getMaxValue() ) + fabs ( fArr.sum() - expectedResult ) << endl;
    failed |= ( fabs ( 0 - fArr.getMinValue() ) + fabs ( 1 - fArr.getMaxValue() ) + fabs ( fArr.sum() - expectedResult ) > prolongTolerance );

    testFailure ( failed, "<----- failed!" );
    cdOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (complicated domain) restriction ... ";
    bool failed = false;

    tpcfe::CFERestrictOp< GridType > reOp ( coarseGrid, fineGrid, aol::ONTHEFLY );
    tpcfe::CFERestrictOp< GridType > reOpA ( coarseGrid, fineGrid, aol::ASSEMBLED );
    failed |= !compareOps ( reOp, reOpA, fineGrid.getNumberOfNodes(), coarseGrid.getNumberOfNodes() );

    qc::ScalarArray<RealType, qc::QC_3D> fArr ( fineGrid ), fMOne ( fineGrid ), cArr ( coarseGrid ), cOne ( coarseGrid );
    fArr.setAll ( 1.0 ); fineGrid.restrictToDofs ( fArr );
    cOne.setAll ( 1.0 ); coarseGrid.restrictToDofs ( cOne );
    MatrixType massMat ( fineGrid );
    {
      tpcfe::CFEMassOp<ConfiguratorType> massOp ( fineGrid );
      massOp._quietMode = true;
      massOp.assembleAddMatrix ( massMat );
    }
    tpcfe::restrictNonDOFEntries ( fineGrid, massMat, 1.0 );
    massMat.apply ( fArr, fMOne );
    reOp.apply ( fMOne, cArr );
    const RealType expectedResult = ( useDirichlet ? 5763. / 10000. + 8. / 90000. : 2.0 / 3.0 );
    //cerr << endl << fabs ( cArr * cOne - expectedResult ) << endl;
    failed |= ( fabs ( cArr * cOne - expectedResult ) > massMatSumTolerance );

    testFailure ( failed, "<----- failed!" );
    cdOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (complicated domain) operator coarsening ... ";
    bool failed = false;

    tpcfe::CFEMassOp<ConfiguratorType> fineMassOp ( fineGrid ), coarseMassOp ( coarseGrid );
    fineMassOp._quietMode = coarseMassOp._quietMode = true;
    MatrixType fineMassMat ( fineGrid ), coarseMassMat ( coarseGrid ), coarsenedMassMatMM ( coarseGrid ), coarsenedMassMatOTF ( coarseGrid );
    fineMassOp.assembleAddMatrix ( fineMassMat );
    coarseMassOp.assembleAddMatrix ( coarseMassMat );

    tpcfe::restrictNonDOFEntries ( fineGrid, fineMassMat, 1.0 );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarseMassMat, 1.0 );

    tpcfe::coarsenCFEMatrix<ConfiguratorType> ( fineMassMat, coarsenedMassMatOTF, fineGrid, coarseGrid, 0.0 );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarsenedMassMatOTF, 1.0 );
    //cerr << endl << aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarseMassMat, coarsenedMassMatOTF, coarseGrid.getNumberOfNodes() ) << endl;
    failed |= ( aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarseMassMat, coarsenedMassMatOTF, coarseGrid.getNumberOfNodes() ) > coarsenedCompareTolerance );

    tpcfe::CFEProlongOp< GridType > prOpA ( coarseGrid, fineGrid, aol::ASSEMBLED );
    tpcfe::CFERestrictOp< GridType > reOpA ( coarseGrid, fineGrid, aol::ASSEMBLED );
    aol::SparseMatrix<RealType> prMat ( fineGrid.getNumberOfNodes(), coarseGrid.getNumberOfNodes() ), reMat ( coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() ), tmpMat ( coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() );
    prOpA.assembleAddMatrix ( prMat );
    reOpA.assembleAddMatrix ( reMat );

    aol::SparseMatrix<RealType> prMatT ( coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() );
    prMat.transposeTo ( prMatT );
    //cerr << endl << aol::rowwiseOpLinfDifference< aol::SparseMatrix<RealType>, aol::SparseMatrix<RealType>, RealType > ( reMat, prMatT, coarseGrid.getNumberOfNodes() ) << endl;
    failed |= ( aol::rowwiseOpLinfDifference< aol::SparseMatrix<RealType>, aol::SparseMatrix<RealType>, RealType > ( reMat, prMatT, coarseGrid.getNumberOfNodes() ) > coarsenedCompareTolerance );

    tmpMat.addMatrixProduct ( reMat, fineMassMat );
    coarsenedMassMatMM.addMatrixProduct ( tmpMat, prMat );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarsenedMassMatMM, 1.0 );
    //cerr << endl << aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarsenedMassMatOTF, coarsenedMassMatMM, coarseGrid.getNumberOfNodes() ) << endl;
    failed |= ( aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarsenedMassMatOTF, coarsenedMassMatMM, coarseGrid.getNumberOfNodes() ) > coarsenedCompareTolerance );

    testFailure ( failed, "<----- failed!" );
    cdOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (complicated domain) multigrid solver ... ";
    bool failed = false;

    tpcfe::CFEMassOp<ConfiguratorType> fineMassOp ( fineGrid );
    fineMassOp._quietMode = true;
    tpcfe::CFEStiffOp<ConfiguratorType> fineStiffOp ( fineGrid );
    fineStiffOp._quietMode = true;
    ConfiguratorType::MatrixType sysMat ( fineGrid );
    fineStiffOp.assembleAddMatrix ( sysMat );
    sysMat *= fineGrid.H();
    fineMassOp.assembleAddMatrix ( sysMat );
    tpcfe::restrictNonDOFEntries ( fineGrid, sysMat, 1.0 );
    aol::NoiseOperator<RealType> noiseOp;
    qc::ScalarArray<RealType, qc::QC_3D> orig ( fineGrid ), rhs ( fineGrid ), soln ( fineGrid ), difference ( fineGrid );
    noiseOp.apply ( rhs, orig );
    fineGrid.restrictToDofs ( orig );
    sysMat.apply ( orig, rhs );

    tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ONTHEFLY, mg::MATRIXMULT_COARSENING > mgsolver ( fineGrid, sysMat, 0 /*explicite level*/ );
    mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mgsolver.setVerboseMode ( 1 );
    mgsolver.apply ( rhs, soln );

    difference = soln;
    difference -= orig;

    //cerr << endl << difference.norm() / orig.norm() << endl;
    failed |= ( difference.norm() / orig.norm() ) > solutionOrigTolerance;

    testFailure ( failed, "<----- failed!" );
    cdOK &= !failed;
  }

  return ( cdOK );
}


template< tpcfe::ConstraintType AT >
bool testJumpingCoefficients ( const RealType prolongTolerance, const RealType restrictTolerance, const RealType coarsenedCompareTolerance, const RealType coarsenedTolerance, const RealType solutionOrigTolerance, const bool useDirichlet ) {

  typedef tpcfe::CFEGrid < RealType, AT >                   GridType;
  typedef tpcfe::CFEHybridMatrix< GridType >                MatrixType;
  typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
  typedef tpcfe::CFEMassOp  < ConfiguratorType >            MassOpType;

  bool jcOK = true;

  const RealType DC_PLUS = 1.0, DC_MINUS = 2.0, thres = 1.0 / 3.0;

  const int fineLevel = 4, coarseLevel = fineLevel - 1;

  GridType fineGrid ( fineLevel );

  qc::ScalarArray<RealType, qc::QC_3D> fineLevelset ( fineGrid );  qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( fineLevelset, thres );
  fineGrid.addStructureFrom ( fineLevelset );

  qc::AArray<RealType, qc::QC_3D> fineCoeff ( fineGrid );
  for ( int i = 0; i < fineLevelset.size(); ++i )
    fineCoeff[i] = ( fineLevelset[i] < 0 ? DC_MINUS : DC_PLUS );
  fineGrid.detectAndInitVirtualNodes ( fineCoeff );

  qc::BitArray<qc::QC_3D> fineDirichlet ( qc::GridSize<qc::QC_3D>::createFrom ( fineGrid ) );
  if ( useDirichlet ) {
    for ( int i = 0; i < fineGrid.getNumX(); ++i ) {
      for ( int j = 0; j < fineGrid.getNumY(); ++j ) {
        fineDirichlet.set ( i, j, 0, true );
        fineDirichlet.set ( i, j, fineGrid.getNumZ() - 1, true );
      }
    }
  } else {
    fineDirichlet.setAll ( false );
  }
  fineGrid.setDirichletMask ( fineDirichlet );
  fineGrid.setDOFMaskFromDirichletAndDomainNodeMask();

  GridType coarseGrid ( coarseLevel );

  qc::ScalarArray<RealType, qc::QC_3D> coarseLevelset ( coarseGrid ); qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( coarseLevelset, thres );
  coarseGrid.addStructureFrom ( coarseLevelset );

  qc::AArray<RealType, qc::QC_3D> coarseCoeff ( coarseGrid );
  for ( int i = 0; i < coarseLevelset.size(); ++i )
    coarseCoeff[i] = ( coarseLevelset[i] < 0 ? DC_MINUS : DC_PLUS );
  coarseGrid.detectAndInitVirtualNodes ( coarseCoeff );

  tpcfe::coarsenDirichletMask ( fineGrid, coarseGrid ); // does the right thing in non-Dirichlet case
  coarseGrid.setDOFMaskFromDirichletAndDomainNodeMask();

  tpcfe::SlopeInterface<GridType> fineSlI ( fineGrid ); fineSlI.determineSlopesAndCoarseningWeights ( fineGrid );

  {
    cerr << "--- Testing CFE (jumping coefficient) prolongation ... ";
    bool failed = false;

    qc::ScalarArray<RealType, qc::QC_3D> fArr ( fineGrid ), cArr ( coarseGrid );
    for ( qc::RectangularIterator<qc::QC_3D> bit ( fArr ); bit.notAtEnd(); ++bit ) {
      const RealType y = fineGrid.H() * ( *bit ) [1];
      fArr.set ( *bit, ( y < thres ? y : thres + ( y - thres ) / DC_MINUS ) + fineGrid.H() * ( *bit ) [0] + 2 * fineGrid.H() * ( *bit ) [2] );
    }
    fineGrid.restrictToDofs ( fArr );
    for ( qc::RectangularIterator<qc::QC_3D> bit ( cArr ); bit.notAtEnd(); ++bit ) {
      const RealType y = coarseGrid.H() * ( *bit ) [1];
      cArr.set ( *bit, ( y < thres ? y : thres + ( y - thres ) / DC_MINUS ) + coarseGrid.H() * ( *bit ) [0] + 2 * coarseGrid.H() * ( *bit ) [2] );
    }
    coarseGrid.restrictToDofs ( cArr );

    tpcfe::CFEProlongOp< GridType > prOp ( coarseGrid, fineGrid, fineSlI, aol::ONTHEFLY );
    tpcfe::CFEProlongOp< GridType > prOpA ( coarseGrid, fineGrid, fineSlI, aol::ASSEMBLED );
    failed |= ( !compareOps ( prOp, prOpA, coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() ) );

    if ( !useDirichlet ) {
      qc::ScalarArray<RealType, qc::QC_3D> refinedArr ( fineGrid ), difference ( fineGrid );
      prOp.apply ( cArr, refinedArr );

      difference = refinedArr;
      difference -= fArr;

      // cerr << endl << difference.getMaxAbsValue() << endl;
      failed |= ( difference.getMaxAbsValue() > prolongTolerance );
    } // else do not compare; function used here incompatible with zero Dirichlet BC

    testFailure ( failed, "<----- failed!" );
    jcOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (jumping coefficient) restriction ... ";
    bool failed = false;

    tpcfe::CFERestrictOp< GridType > reOp ( coarseGrid, fineGrid, fineSlI, aol::ONTHEFLY );
    tpcfe::CFERestrictOp< GridType > reOpA ( coarseGrid, fineGrid, fineSlI, aol::ASSEMBLED );

    failed |= !compareOps ( reOp, reOpA, fineGrid.getNumberOfNodes(), coarseGrid.getNumberOfNodes() );

    qc::ScalarArray<RealType, qc::QC_3D> fArr ( fineGrid ), fMOne ( fineGrid ), cArr ( coarseGrid ), cOne ( coarseGrid );
    fArr.setAll ( 1.0 ); fineGrid.restrictToDofs ( fArr );
    cOne.setAll ( 1.0 ); coarseGrid.restrictToDofs ( cOne );
    MatrixType massMat ( fineGrid );
    {
      tpcfe::CFEMassOp<ConfiguratorType> massOp ( fineGrid );
      massOp._quietMode = true;
      massOp.assembleAddMatrix ( massMat );
    }
    tpcfe::restrictNonDOFEntries ( fineGrid, massMat, 1.0 );
    massMat.apply ( fArr, fMOne );
    reOp.apply ( fMOne, cArr );

    // CFE_LIEHR: 0.864617151331018241
    // CFE_TPOS:  0.864584275385668843
    const RealType DirichletEValue = ( AT == tpcfe::CFE_LIEHR ? 0.864617151331018241 : 0.864584275385668843 );
    const RealType expectedValue = ( useDirichlet ? DirichletEValue  : 1.0 );
    // cerr << endl << fabs ( cArr * cOne - expectedValue ) << endl;
    failed |= ( fabs ( cArr * cOne - expectedValue ) > restrictTolerance );

    testFailure ( failed, "<----- failed!" );
    jcOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (jumping coefficient) operator coarsening ... ";
    bool failed = false;

    tpcfe::CFEMassOp<ConfiguratorType> fineOp ( fineGrid ), coarseOp ( coarseGrid );

    fineOp._quietMode = coarseOp._quietMode = true;
    MatrixType fineMat ( fineGrid ), coarseMat ( coarseGrid ), coarsenedMatMM ( coarseGrid ), coarsenedMatOTF ( coarseGrid );
    fineOp.assembleAddMatrix ( fineMat );
    coarseOp.assembleAddMatrix ( coarseMat );

    tpcfe::restrictNonDOFEntries ( fineGrid, fineMat, 1.0 );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarseMat, 1.0 );

    tpcfe::coarsenCFEMatrix<ConfiguratorType> ( fineMat, coarsenedMatOTF, fineGrid, coarseGrid, 0.0, 0.0, tpcfe::MGWeightProvider< GridType > ( fineSlI ) );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarsenedMatOTF, 1.0 );

    // coarsened is not the same as coarse, but Frobenius norm of the difference is an OK measure
    MatrixType diffMat ( coarseGrid );
    diffMat += coarseMat;
    diffMat -= coarsenedMatOTF;
    // cerr << "difference vs coarse Frobenius " << endl << ( diffMat.getFrobeniusNormSqr() / coarseMat.getFrobeniusNormSqr() ) << endl;

    failed |= ( ( diffMat.getFrobeniusNormSqr() / coarseMat.getFrobeniusNormSqr() ) > coarsenedCompareTolerance );

    tpcfe::CFEProlongOp< GridType > prOpA ( coarseGrid, fineGrid, fineSlI, aol::ASSEMBLED );
    tpcfe::CFERestrictOp< GridType > reOpA ( coarseGrid, fineGrid, fineSlI, aol::ASSEMBLED );
    aol::SparseMatrix<RealType> prMat ( fineGrid.getNumberOfNodes(), coarseGrid.getNumberOfNodes() ), reMat ( coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() ), tmpMat ( fineGrid.getNumberOfNodes(), coarseGrid.getNumberOfNodes() );
    prOpA.assembleAddMatrix ( prMat );
    reOpA.assembleAddMatrix ( reMat );

    aol::SparseMatrix<RealType> prMatT ( coarseGrid.getNumberOfNodes(), fineGrid.getNumberOfNodes() );
    prMat.transposeTo ( prMatT );
    // cerr << endl << aol::rowwiseOpLinfDifference< aol::SparseMatrix<RealType>, aol::SparseMatrix<RealType>, RealType > ( reMat, prMatT, coarseGrid.getNumberOfNodes() ) << endl;
    failed |= ( aol::rowwiseOpLinfDifference< aol::SparseMatrix<RealType>, aol::SparseMatrix<RealType>, RealType > ( reMat, prMatT, coarseGrid.getNumberOfNodes() ) > coarsenedTolerance );

    tmpMat.addMatrixProduct ( fineMat, prMat );
    coarsenedMatMM.addMatrixProduct ( reMat, tmpMat );
    tpcfe::restrictNonDOFEntries ( coarseGrid, coarsenedMatMM, 1.0 );
    // cerr << "coarsenedMM vs coarsenedOTF " << endl << aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarsenedMatMM, coarsenedMatOTF, coarseGrid.getNumberOfNodes() ) << endl;
    failed |= ( aol::rowwiseOpLinfDifference< MatrixType, MatrixType, RealType > ( coarsenedMatMM, coarsenedMatOTF, coarseGrid.getNumberOfNodes() ) > coarsenedTolerance );

    testFailure ( failed, "<----- failed!" );
    jcOK &= !failed;
  }

  {
    cerr << "--- Testing CFE (jumping coefficient) multigrid solver ... ";
    bool failed = false;

    tpcfe::CFEMassOp<ConfiguratorType> fineMassOp ( fineGrid );
    tpcfe::CFEStiffOpWI<ConfiguratorType> fineStiffOp ( fineCoeff, fineGrid );
    fineMassOp._quietMode = true;
    fineStiffOp._quietMode = true;
    MatrixType sysMat ( fineGrid );
    fineStiffOp.assembleAddMatrix ( sysMat );
    sysMat *= fineGrid.H();
    fineMassOp.assembleAddMatrix ( sysMat );
    aol::NoiseOperator<RealType> noiseOp;
    qc::ScalarArray<RealType, qc::QC_3D> orig ( fineGrid ), rhs ( fineGrid ), soln ( fineGrid ), difference ( fineGrid );
    noiseOp.apply ( rhs, orig );
    sysMat.apply ( orig, rhs );

    tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( fineGrid, sysMat, 1 /*explicite level*/ ); // onthefly coarsening not implemented for CFE_TPOS
    mgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mgsolver.setVerboseMode ( 1 );
    mgsolver.apply ( rhs, soln );

    difference = soln;
    difference -= orig;

    // cerr << endl << difference.norm() / orig.norm() << endl;
    failed |= ( difference.norm() / orig.norm() > solutionOrigTolerance );
    testFailure ( failed, "<----- failed!" );
    jcOK &= !failed;
  }

  return ( jcOK );
}


int main ( int, char** ) {

  aol::StopWatch timer;
  timer.start();

  bool allOk = true;
  cerr << "Testing CFE_CD" << endl;
  allOk &= testComplicatedDomain ( 1.0e-17, 2.0e-15, 1.0e-17, 5.0e-7, false );
  allOk &= testComplicatedDomain ( 1.0e-17, 2.0e-15, 1.0e-17, 5.0e-7, true  );

  cerr << endl << "Testing CFE_LIEHR" << endl;
  allOk &= testJumpingCoefficients< tpcfe::CFE_LIEHR > ( 7.0e-02, 1.0e-15, 6.0e-05, 1.0e-15, 5.0e-07, false ); // prolongation does not give fine result
  allOk &= testJumpingCoefficients< tpcfe::CFE_LIEHR > ( 7.0e-02, 1.0e-15, 1.0e-10, 1.0e-15, 7.0e-07, true  );

  cerr << endl << "Testing CFE_TPOS" << endl;
  allOk &= testJumpingCoefficients< tpcfe::CFE_TPOS  > ( 5.0e-02, 1.0e-15, 5.0e-04, 1.0e-15, 6.0e-8, false );
  allOk &= testJumpingCoefficients< tpcfe::CFE_TPOS  > ( 5.0e-02, 1.0e-15, 1.0e-06, 1.0e-15, 8.0e-7, true  );

  timer.stop();

  cerr << "Comparison took " << timer.elapsedCpuTime() << " seconds (cputime): " << 3898.58 / timer.elapsedCpuTime() << " nupsi." << endl;

  if ( allOk ) {
    cerr << aol::color::ok
         << " ============== " << endl
         << " ===   OK   === " << endl
         << " ============== " << endl
         << aol::color::reset;
    return ( EXIT_SUCCESS );
  } else {
    cerr << aol::color::error
         << " ============== " << endl
         << " === FAILED === " << endl
         << " ============== " << endl
         << aol::color::reset;
    return ( EXIT_FAILURE );
  }


}
