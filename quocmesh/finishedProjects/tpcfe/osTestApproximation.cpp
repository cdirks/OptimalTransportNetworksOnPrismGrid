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

// determine the error between a continuous function satisfying a given kink condition and the CFE representation of that function.

// the tpCFE approximation is compared to an approximation with pcw affine FE on the same cubic grid and pcw trilinear FE.

#include <parameterParser.h>

#include <tpCFEGrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include "osTestfunction.h"
#include "osApproximator.h"


static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;

#define DO_CFE 1
#define DO_AFE 1
#define DO_MFE 1

#define USE_AFFINE_VNRECONSTRUCTION 1

typedef double RealType;

int main ( int argc, char** argv) {
  try {

    aol::ParameterParser params ( argc < 2 ? "par/osTestApproximation.par" : argv[1] );
    cerr << "Reading from " << ( argc < 2 ? "par/osTestApproximation.par" : argv[1] ) << endl;

    params.dump();

    aol::Vector<int> types;
    params.getIntVec( "types", types );

    for ( int typeno = 0; typeno < types.size(); ++typeno ) {
      const int TYPE = types[typeno];

      const RealType DIFF_COEFF_PLUS  = 1.0;
      RealType DIFF_COEFF_MINUS = params.getDouble ( "dc_minus_start" );

      for ( int dc_mn = 0; dc_mn < params.getInt ( "dc_minus_steps" ); ++dc_mn ) {

        const RealType kappa = DIFF_COEFF_MINUS / DIFF_COEFF_PLUS;

        const int size = params.getInt ( "stop_depth" );
        aol::Vector<RealType>
          AFE_linf_errors( size ), CFE_linf_errors( size ), MFE_linf_errors( size ),
          AFE_l2_errors( size ), CFE_l2_errors( size ), MFE_l2_errors( size ),
          AFE_h1_errors( size ), CFE_h1_errors( size ), MFE_h1_errors( size ),
          AFE_wh1_errors( size ), CFE_wh1_errors( size ), MFE_wh1_errors( size );

        for ( int depth = params.getInt ( "start_depth" ); depth < params.getInt ( "stop_depth" ); ++depth ) {

          tpcfe::CFEGrid < RealType, AT > grid ( depth );

          grid.setVerboseMode(true);
          grid.setAdjustLevelset( 0.0 );

          qc::ScalarArray<RealType, qc::QC_3D>
            interfc ( grid ),             // levelset function defining the interface
            function ( grid );            // values of the "affine" function

          qc::AArray<RealType, qc::QC_3D>
            coeff ( grid );               // diffusion coefficient at grid points

          tpcfe::InterfaceTestFunction<RealType> *testFct = tpcfe::selectITFByType<RealType> ( TYPE );

          testFct->setCoeffs ( DIFF_COEFF_PLUS, DIFF_COEFF_MINUS );
          testFct->setWidth ( grid.getWidth() );
          testFct->createInterfaceFunctionCoeffs ( interfc, function, coeff );

#if   defined USE_AFFINE_VNRECONSTRUCTION

          grid.addStructureFrom ( interfc );

#elif defined USE_TRILINEAR_VNRECONSTRUCTION

          tpcfe::CFEStructureTrilin<RealType> TrilinLevelset;
          TrilinLevelset.setImageFrom ( interfc );
          grid.addStructureRef ( TrilinLevelset );

#elif defined USE_ANALYTIC_VNRECONSTRUCTION

          interfc.setAll ( aol::NumberTrait<RealType>::NaN ); // make sure values are not used any longer
          tpcfe::CFEStructureAnalytic<RealType> AnalyticLevelset;
          AnalyticLevelset.setITFPointer ( testFct, grid.getNumX() );
          grid.addStructureRef ( AnalyticLevelset );

#else
          abort();
#endif

          grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

          RealType linf_error = 0.0, l2_error = 0.0, h1_error = 0.0, wh1_error = 0, linf_norm = 0.0, l2_norm = 0.0, h1_norm = 0.0, wh1_norm = 0.0;
          cerr.precision ( 10 );

#ifdef DO_CFE
          tpcfe::CFEApproximator< tpcfe::CFEGrid < RealType, AT > > CFEapproximator ( grid, function );
          tpcfe::ComputeApproximationError( grid, *testFct, CFEapproximator, linf_error, l2_error, h1_error, wh1_error, linf_norm, l2_norm, h1_norm, wh1_norm );

          cerr << "L2 norm of function = " << l2_norm << " H1 norm of function = " << h1_norm << " Linf norm of function = " << linf_norm << " weighted H1 norm of function = " << wh1_norm <<  endl
               << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << ", level " << depth << " Linf/ L2 / H1 / wh1      TPOS_CFE approximation error: "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << endl;


          cout << "CFE " << aol::intFormat ( TYPE ) << " " << aol::detailedFormat ( kappa ) << " " << depth << " "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << " "
               << aol::detailedFormat( linf_error / l2_norm ) << " " << aol::detailedFormat( l2_error / l2_norm ) << " " << aol::detailedFormat( h1_error / l2_norm ) << " " << aol::detailedFormat ( wh1_error / l2_norm )
               << endl;

          CFE_linf_errors[depth] = linf_error;
          CFE_l2_errors[depth] = l2_error;
          CFE_h1_errors[depth] = h1_error;
          CFE_wh1_errors[depth] = wh1_error;
#endif

#ifdef DO_AFE
          tpcfe::AffineFEApproximator< tpcfe::CFEGrid < RealType, AT > > AFEapproximator ( grid, function );
          tpcfe::ComputeApproximationError( grid, *testFct, AFEapproximator, linf_error, l2_error, h1_error, wh1_error, linf_norm, l2_norm, h1_norm, wh1_norm );
          cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << ", level " << depth << " Linf/ L2 / H1 / wh1      AffineFE approximation error: "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << endl;

          cout << "AFE " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa ) << " " << depth << " "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << " "
               << aol::detailedFormat( linf_error / l2_norm ) << " " << aol::detailedFormat( l2_error / l2_norm ) << " " << aol::detailedFormat( h1_error / l2_norm ) << " " << aol::detailedFormat ( wh1_error / l2_norm )
               << endl;

          AFE_linf_errors[depth] = linf_error;
          AFE_l2_errors[depth] = l2_error;
          AFE_h1_errors[depth] = h1_error;
          AFE_wh1_errors[depth] = wh1_error;
#endif

#ifdef DO_MFE
          tpcfe::MultilinFEApproximator< tpcfe::CFEGrid < RealType, AT >, RealType > MFEapproximator ( grid, function );
          tpcfe::ComputeApproximationError( grid, *testFct, MFEapproximator, linf_error, l2_error, h1_error, wh1_error, linf_norm, l2_norm, h1_norm, wh1_norm );
          cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << ", level " << depth << " Linf/ L2 / H1 / wh1    MultilinFE approximation error: "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << endl;

          cout << "MFE " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa ) << " " << depth << " "
               << aol::detailedFormat( linf_error ) << " " << aol::detailedFormat( l2_error ) << " " << aol::detailedFormat( h1_error ) << " " << aol::detailedFormat ( wh1_error ) << " "
               << aol::detailedFormat( linf_error / l2_norm ) << " " << aol::detailedFormat( l2_error / l2_norm ) << " " << aol::detailedFormat( h1_error / l2_norm ) << " " << aol::detailedFormat ( wh1_error / l2_norm )
               << endl;

          MFE_linf_errors[depth] = linf_error;
          MFE_l2_errors[depth] = l2_error;
          MFE_h1_errors[depth] = h1_error;
          MFE_wh1_errors[depth] = wh1_error;
#endif

          delete ( testFct );

        } // end depth loop

#ifdef DO_CFE
        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " CFE Linf  convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( CFE_linf_errors[i-1] / CFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " CFE L2    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( CFE_l2_errors[i-1] / CFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " CFE H1    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( CFE_h1_errors[i-1] / CFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " CFE wH1   convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( CFE_wh1_errors[i-1] / CFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;
#endif

#ifdef DO_AFE
        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " AFE Linf  convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( AFE_linf_errors[i-1] / AFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " AFE L2    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( AFE_l2_errors[i-1] / AFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " AFE H1    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( AFE_h1_errors[i-1] / AFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " AFE wH1   convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( AFE_wh1_errors[i-1] / AFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;
#endif

#ifdef DO_MFE
        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " MFE Linf  convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( MFE_linf_errors[i-1] / MFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " MFE L2    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( MFE_l2_errors[i-1] / MFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " MFE H1    convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( MFE_h1_errors[i-1] / MFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;

        cerr << "Type " << aol::intFormat ( TYPE ) << " dc- / dc+ = " << aol::shortFormat ( kappa ) << " MFE wH1   convRates:";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cerr << " " << aol::detailedFormat( log ( MFE_wh1_errors[i-1] / MFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cerr << endl;
#endif

        cout << "Approx    Type DC ratio      ";
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << aol::intFormat ( i-1 ) << " to " << aol::intFormat ( i ) << "          ";
        }
        cout << endl;

#ifdef DO_CFE
        cout << "CFELinfCR " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( CFE_linf_errors[i-1] / CFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFEL2CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( CFE_l2_errors[i-1] / CFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFEH1CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( CFE_h1_errors[i-1] / CFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "CFEwH1CR  " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( CFE_wh1_errors[i-1] / CFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;
#endif

#ifdef DO_AFE
        cout << "AFELinfCR " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( AFE_linf_errors[i-1] / AFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "AFEL2CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( AFE_l2_errors[i-1] / AFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "AFEH1CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( AFE_h1_errors[i-1] / AFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "AFEwH1CR  " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( AFE_wh1_errors[i-1] / AFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;
#endif

#ifdef DO_MFE
        cout << "MFELinfCR " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( MFE_linf_errors[i-1] / MFE_linf_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "MFEL2CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( MFE_l2_errors[i-1] / MFE_l2_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "MFEH1CR   " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( MFE_h1_errors[i-1] / MFE_h1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;

        cout << "MFEwH1CR  " << aol::intFormat ( TYPE ) << " " << aol::shortFormat ( kappa );
        for ( int i= params.getInt ( "start_depth" ) + 1 ; i < params.getInt ( "stop_depth" ); ++i ) {
          cout << " " << aol::detailedFormat( log ( MFE_wh1_errors[i-1] / MFE_wh1_errors[i] ) / log ( 2.0 ) );
        }
        cout << endl;
#endif

        DIFF_COEFF_MINUS *= params.getDouble ( "dc_minus_stepfactor" );
      } // dc_plus loop

    } // type loop

  } catch ( aol::Exception &exc ) {

    exc.dump();

    return ( EXIT_FAILURE );

  }

  return ( EXIT_SUCCESS );
}
