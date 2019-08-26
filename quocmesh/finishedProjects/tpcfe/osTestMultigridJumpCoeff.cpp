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

#include <scalarArray.h>
#include <indexMapper.h>
#include <parameterParser.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include <boomerAMG.h>

#include "affine.h"
#include "jumpCoeffMG.h"
#include "tpcfe_utils.h"
#include "osTestfunction.h"


typedef double RealType;

static const tpcfe::ConstraintType                        AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT >                   GridType;
typedef tpcfe::CFEHybridMatrix<GridType>                  MatrixType;

typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEMassOp  < ConfiguratorType >            MassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >          WStiffOpType;

int main ( int argc, char** argv ) {
  try {

    aol::ParameterParser params ( argc < 2 ? "par/osTestMultigridJumpCoeff.par" : argv[1] );
    cerr << "Reading from " << ( argc < 2 ? "par/osTestMultigridJumpCoeff.par" : argv[1] ) << endl;

    params.dump();

    const RealType eps = params.getDouble ( "epsilon" );

    aol::Vector<RealType> dcMinusValues;
    params.getRealVec ( "dcMinusValues", dcMinusValues );

    const RealType DC_PLUS = 1.0;
    const int level = params.getInt ( "level" ), explicitLevel = params.getInt ( "explicitLevel" ), numSmooth = params.getInt ( "numSmooth" ), dcMinusSteps = dcMinusValues.size();


    aol::Vector<int>    MGSteps ( dcMinusSteps ), CGSteps ( dcMinusSteps ), PCGSteps ( dcMinusSteps ), AMGSteps ( dcMinusSteps );
    aol::Vector<double> MGTimes ( dcMinusSteps ), CGTimes ( dcMinusSteps ), PCGTimes ( dcMinusSteps ), AMGTimes ( dcMinusSteps ), kappaValues ( dcMinusSteps );


    for ( short s = 0; s < dcMinusSteps ; ++s ) {
      const RealType DC_MINUS = dcMinusValues[s];
      kappaValues[s] = DC_MINUS / DC_PLUS;

      cerr << "DC_MINUS = " << DC_MINUS << endl;

      GridType grid ( level );

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

      {
        tpcfe::InterfaceTestFunction<RealType> *testFct = tpcfe::selectITFByType<RealType> ( params.getInt ( "interfaceType" ) );
        testFct->setWidth ( grid.getNumXYZ() );

        if ( params.getInt ( "interfaceType" ) == 13 ) {
          dynamic_cast< tpcfe::SwissCheeseInterfaceTestFunction<RealType>* > ( testFct )->setFreq ( params.getInt ( "numHoles" ) );
        }

        testFct->createInterface ( levelset );
        delete ( testFct );
      }

      grid.addStructureFrom ( levelset );

      qc::AArray<RealType, qc::QC_3D> coeff ( grid );
      for ( int i = 0; i < levelset.size(); ++i )
        coeff[i] = ( levelset[i] < 0 ? DC_MINUS : DC_PLUS );

      grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), solnMG ( grid ), solnCG ( grid ), solnPCG ( grid ), solnAMG ( grid );

      const int problemType = params.getInt ( "problemType" );
      // preparations before creating matrix
      if ( problemType == -1 ) { // elliptic BVP

        qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
        for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
          if ( ( *bnit ) [1] == 0 )
            DirichletMask.set ( *bnit, true );

          if ( ( *bnit ) [1] == grid.getNumY() - 1 )
            DirichletMask.set ( *bnit, true );
        }

        grid.setDirichletMask ( DirichletMask );

      } else {
        // do nothing
      }

      MatrixType sysMat ( grid );

      if ( problemType == -1 ) { // elliptic BVP

        qc::ScalarArray<RealType, qc::QC_3D> bcond ( grid ), Lbcond ( grid );
        for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
          if ( ( *bnit ) [1] == 0 )
            bcond.set ( *bnit, 0.0 );

          if ( ( *bnit ) [1] == grid.getNumY() - 1 )
            bcond.set ( *bnit, 1.0 );
        }

        cerr << "Setting up stiffOp" << endl;
        WStiffOpType stiffOp ( coeff, grid, aol::ONTHEFLY );
        cerr << "Assembling matrix" << endl;
        stiffOp.assembleAddMatrix ( sysMat );

        sysMat.apply ( bcond, Lbcond );
        rhs -= Lbcond;

        // modify matrix for boundary conditions
        tpcfe::enforceZeroDirichletBCs ( grid, sysMat, rhs );

      } else { // heat conduction time step

        const RealType tau = ( problemType == 1 ? grid.H() : 1 );

        MassOpType fineMassOp ( grid );
        WStiffOpType fineStiffOp ( coeff, grid, aol::ONTHEFLY );
        fineMassOp._quietMode = true; fineStiffOp._quietMode = true;

        fineStiffOp.assembleAddMatrix ( sysMat );
        sysMat *= tau;
        fineMassOp.assembleAddMatrix ( sysMat );

        // start with noise
        aol::RandomGenerator rg;
        for ( int i = 0; i < rhs.size(); ++i )
          rhs[i] = rg.rReal<RealType>();

      }

      aol::StopWatch timer;


      if ( params.getInt ( "doCG" ) == 1 ) {
        timer.start();
        aol::CGInverse< aol::Vector<RealType> > cgsolver ( sysMat, eps, 1000 );
        cgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
        cgsolver.apply ( rhs, solnCG );
        timer.stop();
        CGTimes[s] = timer.elapsedCpuTime();
        CGSteps[s] = cgsolver.getCount();
        cerr << "TIME ConjugateGradient " << aol::intFormat ( level ) << " " << aol::shortFormat ( kappaValues[s] ) << " "  << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( cgsolver.getCount() ) << " iterations" <<  endl;
      }


      if ( params.getInt ( "doPCG" ) == 1 ) {
        timer.start();
        aol::StopWatch timer2;
        timer2.start();
        aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( sysMat );
        aol::PCGInverse< aol::Vector<RealType> > pcgsolver ( sysMat, prec, eps, 10000 );
        timer2.stop();
        cerr << "Diagonal Preconditioner Setup took " << aol::shortFormat ( timer2.elapsedCpuTime() ) << " seconds" << endl;

        timer2.start();
        pcgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
        pcgsolver.apply ( rhs, solnPCG );
        timer2.stop();
        cerr << "Diagonally Preconditioned CG Solve took " << aol::shortFormat ( timer2.elapsedCpuTime() ) << " seconds" << endl;
        timer.stop();
        PCGTimes[s] = timer.elapsedCpuTime();
        PCGSteps[s] = pcgsolver.getCount();
        cerr << "TIME DiagonallyPreconditionedCG " << aol::intFormat ( level ) << " " << aol::shortFormat ( kappaValues[s] ) << " "  << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( pcgsolver.getCount() ) << " iterations" <<  endl;
      }


      if ( params.getInt ( "doAMG" ) == 1 ) {
#ifdef USE_EXTERNAL_HYPRE
        timer.start();
        mg::BoomerAMGSolver< MatrixType, RealType > AMGSolver ( sysMat, eps, 5000 );
        AMGSolver.apply ( rhs, solnAMG );
        timer.stop();
        AMGTimes[s] = timer.elapsedCpuTime();
        AMGSteps[s] = AMGSolver.getCount();
        cerr << "TIME BoomerAMG "  << aol::intFormat ( level ) << " " << aol::shortFormat ( kappaValues[s] ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( AMGSolver.getCount() ) << " iterations" <<  endl;
#endif
      }

      if ( params.getInt ( "doCFEMG" ) == 1 ) {
        timer.start();
        aol::StopWatch timer2;
        timer2.start();
        const RealType relax = 1.0;
        tpcfe::CFEMultigrid< tpcfe::CFEMassOp<ConfiguratorType>, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > mgsolver ( grid, sysMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, relax,
            mg::MULTIGRID_V_CYCLE, eps, 5000, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
        timer2.stop();
        cerr << "Multigrid Setup took " << aol::shortFormat ( timer2.elapsedCpuTime() ) << " seconds" << endl;

        mgsolver.setVerboseMode ( 1 );
        timer2.start();
        mgsolver.apply ( rhs, solnMG );
        timer2.stop();
        cerr << "Multigrid Solve took " << aol::shortFormat ( timer2.elapsedCpuTime() ) << " seconds" << endl;
        timer.stop();
        MGTimes[s] = timer.elapsedCpuTime();
        MGSteps[s] = mgsolver.getCount();
        cerr << "TIME Multigrid " << aol::intFormat ( level ) << " " << aol::shortFormat ( kappaValues[s] ) << " " << aol::shortFormat ( timer.elapsedCpuTime() ) << " seconds " << aol::intFormat ( mgsolver.getCount() ) << " iterations" <<  endl;
      }

    }

    cerr << "Statistics:" << endl;
    cerr << "kappa    & CG Steps   & PCG Steps   & AMG Steps   & MG Steps   &" << endl;
    for ( int s = 0; s < dcMinusSteps; ++s )
      cerr << aol::detailedFormat ( kappaValues[s] ) << " & "
      << aol::intFormat (  CGSteps[s] ) << " & "
      << aol::intFormat ( PCGSteps[s] ) << " & "
      << aol::intFormat ( AMGSteps[s] ) << " & "
      << aol::intFormat (  MGSteps[s] ) << " &" << endl;

    cerr << "CG Times   & PCG Times   & AMG Times   & MG Times   &" << endl;
    for ( int s = 0; s < dcMinusSteps; ++s )
      cerr << aol::shortFormat (  CGTimes[s] ) << " & "
      << aol::shortFormat ( PCGTimes[s] ) << " & "
      << aol::shortFormat ( AMGTimes[s] ) << " & "
      << aol::shortFormat (  MGTimes[s] ) << " &" << endl;

  } catch ( aol::Exception &ex ) {

    ex.dump();
    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );

}
