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

// this is a test program by OS solving an elliptic problem for
// analytically computed right hand side to test the tpCFE operators

// goal: observe 2nd order convergence in L2 norm

#include <parameterParser.h>
#include <preconditioner.h>

#include <tpCFEGrid.h>
#include <tpCFEMatrices.h>
#include <tpCFEStandardOp.h>

#include "tpcfe_utils.h"

#include "osTestfunction.h"
#include "osApproximator.h"

typedef double RealType;

#define USE_PCG_SOLVER 1

static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT > GridType;

typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;

typedef tpcfe::CFEMassOp  < ConfiguratorType >       MassOpType;
typedef tpcfe::CFEStiffOpWI < ConfiguratorType >   WStiffOpType;


int main ( int argc, char** argv ) {
  try {

    aol::ParameterParser params ( argc < 2 ? "par/osTestElliptic.par" : argv[1] );
    cerr << "Reading from " << ( argc < 2 ? "par/osEestElliptic.par" : argv[1] ) << endl;

    aol::Vector<int> types;
    params.getIntVec( "types", types );

    for ( int typeno = 0; typeno < types.size(); ++typeno ) {
      const int TYPE = types[typeno];

      RealType
        DIFF_COEFF_PLUS  = 1.0,
        DIFF_COEFF_MINUS = params.getDouble ( "dc_minus_start" );

      for ( int dc_mn = 0; dc_mn < params.getInt ( "dc_minus_steps" ); ++dc_mn ) {

        RealType kappa = DIFF_COEFF_MINUS / DIFF_COEFF_PLUS;

        const int size = params.getInt( "stop_depth" );
        aol::Vector<double>
          pcg_l2_errors( size ), pcg_linf_errors( size ), pcg_h1_errors( size ), pcg_wh1_errors ( size );

        for ( int depth = params.getInt ( "start_depth" ); depth < params.getInt ( "stop_depth" ); ++depth ) {

          GridType grid ( depth );
          grid.setCoarsenedFlag ( true ); // so that the hybrid matrix can be used with non-constant coefficients
          grid.setAdjustLevelset ( 0.0 );
          grid.setAdjustVirtualNodes ( 0.0 );

          qc::ScalarArray<RealType, qc::QC_3D>
            interfc ( grid ),             // interface
            orig ( grid ),                // analytical solution U
            bcond ( grid ),               // analytical boundary values
            eff ( grid ),                 // right hand side f: analytical laplacian times corresponding coefficient
            L_bc ( grid ),                // L times boundary values
            rhs ( grid ),                 // right hand side vector for solver
            dummy ( grid ),               // dummy vector
            solution ( grid ),            // solution obtained by pcg solver
            difference ( grid );          // difference between pcg and analytical solution

          qc::AArray<RealType, qc::QC_3D> diff_coeff ( grid );  // diffusion coefficient at grid points

          qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
          DirichletMask.setAll ( false );

          const int width = grid.getWidth();

          tpcfe::InterfaceTestFunction<RealType> *testFct = tpcfe::selectITFByType<RealType> ( TYPE );

          testFct->setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS );
          testFct->setWidth ( width );
          testFct->createInterfaceFunctionCoeffs ( interfc, orig, diff_coeff );

          grid.addStructureFrom ( interfc );
          grid.relaxedDetectAndInitVirtualNodes ( diff_coeff, 1.0e-16, 1.0e-16 );

          for ( qc::RectangularIterator<qc::QC_3D> fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
            aol::Vec3<RealType> pos ( ( *fnit ) [0], ( *fnit ) [1], ( *fnit ) [2] );
            eff.set ( ( *fnit ), - ( testFct->coefficientValueAtGlobal ( pos ) ) * ( testFct->Laplace_valueAt_global ( pos ) ) );

            if ( fnit->x() == 0 || fnit->y() == 0 || fnit->z() == 0 || fnit->x() == ( width - 1 ) || fnit->y() == ( width - 1 ) || fnit->z() == ( width - 1 ) ) {
              bcond.set ( *fnit,  orig.get( *fnit ) ); // may comment out if BC = 0
              DirichletMask.set ( *fnit, true );
            }
          }

          grid.setDirichletMask ( DirichletMask );
          grid.setDOFMaskFromDirichletAndDomainNodeMask(); // probably not necessary

          {
            MassOpType massOp ( grid, aol::ASSEMBLED );
            massOp.apply ( eff, rhs );
          }

          MatrixType stiffMat ( grid );
          {
            WStiffOpType stiffOp ( diff_coeff, grid, aol::ONTHEFLY );
            stiffOp.assembleAddMatrix ( stiffMat );
          }

          stiffMat.apply ( bcond, L_bc ); // may comment out if BC = 0
          rhs -= L_bc; // may comment out if BC = 0 => L_bc = 0

          // modify matrix for boundary conditions
          tpcfe::enforceZeroDirichletBCs ( grid, stiffMat, rhs );


#ifdef USE_PCG_SOLVER
          // aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
          aol::SSORPreconditioner< aol::Vector<RealType>, MatrixType > prec ( stiffMat );
          aol::PCGInverse< aol::Vector<RealType> > solver ( stiffMat, prec, params.getDouble ( "solver_eps" ), params.getInt ( "solver_steps" ) );
          solver.setQuietMode ( false );
          solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
#endif

#ifdef USE_MULTIGRID_SOLVER
          const int explicitLevel = params.getInt ( "explicitLevel" ), numSmooth = params.getInt ( "numSmooth" );
          tpcfe::CFEMultigrid< WStiffOpType, aol::ASSEMBLED, mg::MATRIXMULT_COARSENING > solver ( grid, stiffMat, explicitLevel, numSmooth, numSmooth, aol::GAUSS_SEIDEL_FORWARD, 1.0,
                                                                                                  mg::MULTIGRID_V_CYCLE, params.getDouble ( "solver_eps" ), params.getInt ( "solver_steps" ), aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

#endif
          solver.apply ( rhs, solution );

          cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

          solution += bcond;

          difference  = solution;
          difference -=    orig;

          cerr << "difference.norm() = " << difference.norm() << ", solution.norm() =  " << solution.norm() << endl;

          RealType pcg_linf_difference = 0, pcg_l2_difference = 0, pcg_h1_difference = 0, pcg_wh1_difference = 0, pcg_linf_norm = 0, pcg_l2_norm = 0, pcg_h1_norm = 0, pcg_wh1_norm = 0;

          tpcfe::CFEApproximator< GridType > MyApproximator ( grid, solution );

#if 1
          tpcfe::ComputeApproximationError( grid, *testFct, MyApproximator,
                                            pcg_linf_difference, pcg_l2_difference, pcg_h1_difference, pcg_wh1_difference, pcg_linf_norm, pcg_l2_norm, pcg_h1_norm, pcg_wh1_norm );
#else
          // debugging output
          qc::ScalarArray<RealType, qc::QC_3D> errors ( grid );
          tpcfe::ComputeApproximationError( grid, *testFct, MyApproximator,
                                            pcg_linf_difference, pcg_l2_difference, pcg_h1_difference, pcg_wh1_difference,
                                            pcg_linf_norm, pcg_l2_norm, pcg_h1_norm, pcg_wh1_norm,
                                            &errors, 1, 0.0 );

          cerr << "ERRORS IN RANGE " << errors.getMinValue() << " " << errors.getMaxValue() << endl;
          errors.saveSlices( "out/errors%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, errors.getMinValue(), errors.getMaxValue() );
          errors.saveRaw ( "out/errors.raw", qc::PGM_UNSIGNED_CHAR_BINARY, errors.getMaxValue() );
          interfc.save ( "out/interface.dat.bz2", qc::PGM_DOUBLE_BINARY );
          errors.save ( "out/errors.dat.bz2", qc::PGM_DOUBLE_BINARY );
#endif
          cerr << pcg_l2_difference << " " << pcg_h1_difference << " " << pcg_linf_difference << endl;
          cout << "CFE " << aol::intFormat ( TYPE ) << "E " << aol::detailedFormat ( kappa ) << " " << depth
               << aol::detailedFormat ( pcg_linf_difference ) << " " << aol::detailedFormat ( pcg_l2_difference ) << " " << aol::detailedFormat ( pcg_h1_difference ) << aol::detailedFormat ( pcg_wh1_difference ) << " "
               << aol::detailedFormat ( pcg_linf_difference / pcg_l2_norm ) << " " << aol::detailedFormat ( pcg_l2_difference / pcg_l2_norm ) << " " << aol::detailedFormat ( pcg_h1_difference / pcg_l2_norm ) << aol::detailedFormat ( pcg_wh1_difference / pcg_l2_norm ) << endl;

          pcg_linf_errors[depth] = pcg_linf_difference;
          pcg_l2_errors[depth] = pcg_l2_difference;
          pcg_h1_errors[depth] = pcg_h1_difference;
          pcg_wh1_errors[depth] = pcg_wh1_difference;


        } // depth loop

        cout << "CFE LinfCR " << aol::intFormat ( TYPE ) << "E " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_linf_errors[i-1] / pcg_linf_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFE L2CR   " << aol::intFormat ( TYPE ) << "E " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_l2_errors[i-1] / pcg_l2_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFE H1CR   " << aol::intFormat ( TYPE ) << "E " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_h1_errors[i-1] / pcg_h1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFE wH1CR  " << aol::intFormat ( TYPE ) << "E " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_wh1_errors[i-1] / pcg_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        DIFF_COEFF_MINUS *= params.getDouble ( "dc_minus_stepfactor" );
      } // dc_plus loop
      cout << endl;
    } // type loop
    cout << endl;


  } catch ( aol::Exception &ex ) {

    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );

}
