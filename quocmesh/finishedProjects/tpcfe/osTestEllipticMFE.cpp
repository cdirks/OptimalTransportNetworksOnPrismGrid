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

#include <FEOpInterface.h>
#include <preconditioner.h>
#include <quocMatrices.h>
#include <configurators.h>

#include <parameterParser.h>
#include "tpcfe_utils.h"

#include "osTestfunction.h"
#include "osApproximator.h"


typedef double RealType;

typedef tpcfe::QuocCompatGridDefinition< RealType > QuocGridType;
typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_TPOS > CFEGridType;

typedef qc::MultilinFEBandMatrix<RealType,qc::QC_3D> MatrixType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 23> > ConfiguratorType;

typedef aol::MassOp  < ConfiguratorType >                 MassOpType;
typedef aol::MyMFEWeightedStiffOp < ConfiguratorType >  WStiffOpType;


int main ( int argc, char** argv ) {
  try {
    aol::ParameterParser params ( argc < 2 ? "par/os_test_elliptic.par" : argv[1] );
    cerr << "Reading from " << ( argc < 2 ? "par/os_test_elliptic.par" : argv[1] ) << endl;

    aol::Vector<int> types;
    params.getIntVec( "types", types );

    for ( int typeno = 0; typeno < types.size(); ++typeno ) {
      const int TYPE = types[typeno];

      RealType
        DIFF_COEFF_PLUS  = params.getDouble ( "dc_plus_start" ),
        DIFF_COEFF_MINUS = 1.0;

      for ( int dc_pl = 0; dc_pl < params.getInt ( "dc_plus_steps" ); ++dc_pl ) {

        const int size = params.getInt( "stop_depth" );
        aol::Vector<double>
          pcg_l2_errors( size ), pcg_linf_errors( size ), pcg_h1_errors( size );

        for ( int depth = params.getInt ( "start_depth" ); depth < params.getInt ( "stop_depth" ); ++depth ) {

          QuocGridType quocGrid ( depth );

          qc::ScalarArray<RealType, qc::QC_3D>
            interfc ( quocGrid ),             // interface
            orig ( quocGrid ),                // analytical solution U
            bcond ( quocGrid ),               // analytical boundary values
            eff ( quocGrid ),                 // right hand side f: analytical laplacian times corresponding coefficient
            L_bc ( quocGrid ),                // L times boundary values
            rhs ( quocGrid ),                 // right hand side vector for solver
            dummy ( quocGrid ),               // dummy vector
            pcg_soln ( quocGrid ),            // solution obtained by pcg solver
            pcg_diff ( quocGrid );            // difference between pcg and analytical solution

          qc::AArray<RealType, qc::QC_3D>
            diff_coeff ( quocGrid );          // diffusion coefficient at grid points

          qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( quocGrid ) );

          const int width = quocGrid.getWidth();

          tpcfe::InterfaceTestFunction<RealType> *testFct = tpcfe::selectITFByType<RealType> ( TYPE );

          testFct->setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS );
          testFct->setWidth ( width );
          testFct->createInterfaceFunctionCoeffs ( interfc, orig, diff_coeff );

          for ( qc::RectangularIterator<qc::QC_3D> fnit ( quocGrid ); fnit.notAtEnd(); ++fnit ) {
            aol::Vec3<RealType> pos ( ( *fnit ) [0], ( *fnit ) [1], ( *fnit ) [2] );
            eff.set ( ( *fnit ), - ( testFct->coefficientValueAtGlobal ( pos ) ) * ( testFct->Laplace_valueAt_global ( pos ) ) ); //

            if ( fnit->x() == 0 || fnit->y() == 0 || fnit->z() == 0 || fnit->x() == ( width - 1 ) || fnit->y() == ( width - 1 ) || fnit->z() == ( width - 1 ) ) {
              bcond.set ( *fnit,  orig.get( *fnit ) );
              DirichletMask.set ( *fnit, true );
            }
          }

          MatrixType massMat  ( qc::GridSize<qc::QC_3D>::createFrom ( quocGrid ) ),
                     stiffMat ( qc::GridSize<qc::QC_3D>::createFrom ( quocGrid ) );
          {
            MassOpType massOp ( quocGrid, aol::ONTHEFLY );
            massOp.assembleAddMatrix ( massMat );

            WStiffOpType stiffOp ( quocGrid, *testFct, aol::ONTHEFLY );
            stiffOp.assembleAddMatrix ( stiffMat );
          }
          massMat.apply ( eff, rhs );

          stiffMat.apply ( bcond, L_bc );
          rhs -= L_bc;

          // modify matrix for boundary conditions
          tpcfe::enforceZeroDirichletBCsWholeCube ( quocGrid, stiffMat, rhs );

          aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );

          aol::PCGInverse< aol::Vector<RealType> > pcg_solver ( stiffMat, prec, params.getDouble ( "solver_eps" ), params.getInt ( "solver_steps" ) );
          pcg_solver.setQuietMode ( false );
          pcg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

          pcg_solver.apply ( rhs, pcg_soln );

          pcg_soln += bcond;

          pcg_diff  = pcg_soln;
          pcg_diff -=    orig;


          cerr << "pcg_diff.norm() = " << pcg_diff.norm() << ", pcg_soln.norm() =  " << pcg_soln.norm() << endl;

          RealType pcg_linf_difference = 0, pcg_l2_difference = 0, pcg_h1_difference = 0, pcg_wh1_difference = 0, pcg_linf_norm = 0, pcg_l2_norm = 0, pcg_h1_norm = 0, pcg_wh1_norm = 0;

          CFEGridType cfeGrid ( depth );
          //           qc::ScalarArray<RealType, qc::QC_3D> dummyInterface ( cfeGrid ), dummyDiffCoeff ( cfeGrid );
          //           dummyInterface.setAll ( -1 );
          //           dummyDiffCoeff.setAll ( -1 );
          cfeGrid.addStructureFrom ( interfc );
          //           tpcfe::CFEWeightInterface<RealType> dummyWeightInterface ( dummyDiffCoeff );
          cfeGrid.detectAndInitVirtualNodes ( diff_coeff );

          tpcfe::MultilinFEApproximator< CFEGridType, RealType > MyApproximator ( cfeGrid, pcg_soln );

          tpcfe::ComputeApproximationError( cfeGrid, *testFct, MyApproximator, pcg_linf_difference, pcg_l2_difference, pcg_h1_difference, pcg_wh1_difference, pcg_linf_norm, pcg_l2_norm, pcg_h1_norm, pcg_wh1_norm );

          cerr << pcg_l2_difference << " " << pcg_h1_difference << " " << pcg_linf_difference << endl;
          cout << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( DIFF_COEFF_PLUS ) << " " << depth
               << aol::detailedFormat ( pcg_l2_difference ) << " " << aol::detailedFormat ( pcg_h1_difference ) << " " << aol::detailedFormat ( pcg_linf_difference ) << endl;

          pcg_l2_errors[depth] = pcg_l2_difference;
          pcg_h1_errors[depth] = pcg_h1_difference;
          pcg_linf_errors[depth] = pcg_linf_difference;


        } // depth loop

        char cfet[1024];
        sprintf ( cfet, "%s", "MFE" );

        cout << cfet << " L2CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( DIFF_COEFF_PLUS / DIFF_COEFF_MINUS );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_l2_errors[i-1] / pcg_l2_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << cfet << " H1CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( DIFF_COEFF_PLUS / DIFF_COEFF_MINUS );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_h1_errors[i-1] / pcg_h1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << cfet << " LinfCR " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( DIFF_COEFF_PLUS / DIFF_COEFF_MINUS );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( pcg_linf_errors[i-1] / pcg_linf_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        DIFF_COEFF_PLUS *= params.getDouble ( "dc_plus_stepfactor" );
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
