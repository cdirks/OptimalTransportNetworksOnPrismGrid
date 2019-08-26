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

// include all narrowband header files
#include <narrowBandConfigurators.h>
#include <narrowBandBoundaryIntegration.h>
#include <narrowBandGrid.h>
#include <narrow.h>
#include <subGridSparseMatrix.h>

#include <FEOpInterface.h>
#include <gridBase.h>
#include <quoc.h>

void testSuccess( const bool success, const char *errorText ) {
  if ( success ) {
    cerr << "OK." << endl;
  } else {
    cerr << aol::color::red << errorText << endl << aol::color::reset;
  }
}

int main( int, char** ){

  try {

    bool success = true;

    {
      cerr << "--- Testing nb::SubGridSparseMatrix and the computation only on a subgrid ... " ;

      qc::GridDefinition grid( 3, qc::QC_2D );

      typedef nb::NarrowBandGrid<qc::GridDefinition, qc::QC_2D> NarrowBandGridType;
      typedef nb::NarrowBandConfiguratorTraitMultiLin<double, qc::QC_2D,
                  aol::GaussQuadrature<double,qc::QC_2D,7>, NarrowBandGridType > NBConfType;

      nb::SubGridSparseMatrix<double,NBConfType> matrix( grid.getNumberOfNodes() );
      NarrowBandGridType testNB( grid );

      // define a config file for the two narrowbands (inner and outer)
      NBConfType nbTestConfig( testNB );

      // now fill the hash-table with elements
      qc::GridDefinition::OldFullElementIterator it;
      for ( it = grid.begin(); it != grid.end(); ++it )
      {
        if ( (*it)[0] > 2 && (*it)[0] < 6 && (*it)[1] > 2 && (*it)[1] < 6 )
          testNB.insert( *it );
//      cerr << (*it)[0] << ", " << (*it)[1] << endl;
      }

      // now rebuild the subgridmatrices according to their hash
      cerr << "\nRebuilding matrices ...";
      matrix.rebuild( nbTestConfig );
      // define an operator and assemble it only on the narrow band
      aol::LumpedMassOp<NBConfType> MassOp( testNB, aol::DO_NOT_INVERT );
      matrix.setZero();
      MassOp.assembleAddMatrix( matrix );

      qc::ScalarArray<double, qc::QC_2D> arg( grid );
      qc::ScalarArray<double, qc::QC_2D> dest( grid );
      arg.setAll( 1. );
      dest.setAll( 0. );

      MassOp.apply( arg, dest );

      // count the nodes, where a computation took place
      int num = 0;
      for (int i=0; i<dest.size(); i++)
        if ( dest[i] != 0 ) num++;

      if ( num != 16 ) success = false;


      if(success)
        cerr << "OK" << endl;
    }

    if(success) {
      aol::printSelfTestSuccessMessage ( "--                 NARROWBAND Self Test Successful                            --" );
      aol::callSystemPauseIfNecessaryOnPlatform();
      return 0;
    }
  }
  catch(std::exception &ex){
    cerr << "\n\nstd::exception caught:\n";
    cerr << ex.what () << endl;
  }
  catch(aol::Exception &ex){
    cerr << "\n\naol::Exception caught:\n";
    ex.dump ();
  }
  catch (...){
    cerr << "\n\nUnknown exception caught.\n";
  }
  aol::printSelfTestFailureMessage ( "!!                 NARROWBAND Self Test FAILED                                !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;

}
