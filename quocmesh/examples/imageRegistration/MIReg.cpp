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

/**
 * \file
 *
 * \brief Registers two images using mutual information as similarity measure and a cascadic multiscale minimization algorithm.
 *
 * \author Berkels
 */

#include <mutualInformation.h>

template <typename ConfiguratorType, template<class> class RegistrationConfiguratorType, template<class> class RegularizationConfiguratorType>
void doMLRegistration ( const aol::ParameterParser &Parser ) {
  qc::StandardRegistrationMultilevelDescent<ConfiguratorType, RegistrationConfiguratorType<ConfiguratorType>, RegularizationConfiguratorType<ConfiguratorType> > mld( Parser );
  mld.solve();
}

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
    char parameterfilename[1024];

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf( parameterfilename, "%s",  argv[1] );
    }
    if ( argc == 1 ) {
      if( ConfType::Dim == 3)
        sprintf( parameterfilename, "ParameterFiles/MIReg_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/MIReg_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    aol::StopWatch watch;
    watch.start();

    const int mode = parser.hasVariable ( "mode" ) ? parser.getInt ( "mode" ) : 0;
    switch ( mode ) {
      case 0:
        doMLRegistration<ConfType, qc::MIRegistrationConfigurator, qc::DirichletRegularizationConfigurator> ( parser );
        break;
      case 1:
        doMLRegistration<ConfType, qc::NCCRegistrationConfigurator, qc::DirichletRegularizationConfigurator> ( parser );
        break;
      case 2:
        doMLRegistration<ConfType, qc::NCCRegistrationConfigurator, qc::ROFRegularizationConfigurator> ( parser );
        break;
      case 3:
        doMLRegistration<ConfType, qc::SSDRegistrationConfigurator, qc::ROFRegularizationConfigurator> ( parser );
        break;
      case 4:
        doMLRegistration<ConfType, qc::SSDRegistrationConfigurator, qc::DirichletRegularizationConfigurator> ( parser );
        break;
      case 5:
        doMLRegistration<ConfType, qc::SSDRegistrationConfigurator, qc::HyperelasticEnergyConfigurator> ( parser );
        break;
      case 6:
        doMLRegistration<ConfType, qc::RGBSSDRegistrationConfigurator, qc::DirichletRegularizationConfigurator> ( parser );
        break;
      case 7:
        doMLRegistration<ConfType, qc::NCCRegistrationConfigurator, qc::HyperelasticEnergyConfigurator> ( parser );
        break;
      default:
        throw aol::Exception ( "Unknown mode", __FILE__, __LINE__ );
    } 


    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;

    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}
