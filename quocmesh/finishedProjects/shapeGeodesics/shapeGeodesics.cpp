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

#include "shapeGeodesics.h"

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype, finite element type, used quadrature rules
typedef double RealType;
const qc::Dimension Dim = qc::QC_2D;
//typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::GaussQuadrature<RealType, Dim, 3> > ConfiguratorType;
typedef qc::simplex::ConfiguratorTraitLinear<RealType, Dim, qc::simplex::MidpointQuadrature<RealType, Dim>, qc::simplex::GridStructure<qc::GridDefinition, Dim> > ConfiguratorType;

/**
 * \brief This main function coordinates the geodesic computation between level set function objects.
 * The parameter data is read in and the geodesic algorithm is executed.
 * If started from the command line, argc is the number of passed arguments,
 * and argv is a pointer on the passed arguments (which are of type char[]).
 * The name of the parameter file should be passed as the only argument.
 *
 * \author Wirth
 */
int main( int argc, char *argv[] ) {

  try {
    char parameterFileName[1024];

    // read in file names of parameter files
    if ( argc > 50 ) {             // too many input arguments specified; explain correct syntax
      cerr << "Too many input files. Use n<=49: averagePhaseField <parameter filename 1> ... <parameter filename n>" << endl;
      return EXIT_FAILURE;
    }
    else {
      for ( int i=1; i<argc; i++ ) {
        sprintf( parameterFileName, "%s",  argv[i] );

        // read in parameters from specified parameter file
        cerr << endl << "Reading parameters from '" << parameterFileName << "'." << endl;
        aol::ParameterParser parameterParser( parameterFileName );

        // execute the algorithm
        switch ( parameterParser.getInt( "programType" ) ) {
        // computation of geodesics in shape space
        case 6: ( ShapeGeodesicMultiGridOp<ConfiguratorType>( parameterParser ) ).execute(); break;
        case 7: ( ShapeGeodesicOp<ConfiguratorType>( parameterParser ) ).execute(); break;
        case 8: ( ShapeGeodesicOp<ConfiguratorType>( parameterParser ) ).execute(); break;
        case 9: ( ShapeGeodesicOp<ConfiguratorType>( parameterParser ) ).execute(); break;
        default:
          throw aol::Exception ( "missing default case", __FILE__, __LINE__ );
        }
      }
    }
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
