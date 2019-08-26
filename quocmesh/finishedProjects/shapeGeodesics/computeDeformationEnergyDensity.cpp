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

/*### ********************************************************************************************
 computes the elastic deformation energy density for a given displacement and saves it as ASCII

 Author: Benedikt Wirth
 ******************************************************************************************** ###*/

#include <scalarArray.h>
#include <hyperelastic.h>
#include <configurators.h>

typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
typedef qc::HyperelasticEnergyDensityDefault<ConfiguratorType> MaterialLawType;

int main ( int argc, char **argv ) {

  try {
    if ( argc != 7 && argc != 8 ) {
      cerr << aol::color::red << "USAGE: " << argv[0] << " dispX dispY elastParam1 elastParam2 elastParam3 [weightingMatrix] output" << endl << aol::color::reset;
      return EXIT_FAILURE;
    }

    // create a grid
    ConfiguratorType::ArrayType aux( argv[1] );
    int d = qc::logBaseTwo( aux.getNumX() );
    ConfiguratorType::InitType grid( d, qc::QC_2D );

    // load displacement
    aol::MultiVector<RealType> displacement( ConfiguratorType::Dim, grid.getNumberOfNodes() );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      ConfiguratorType::ArrayType( displacement[i], grid, aol::FLAT_COPY ).load( argv[1+i] );

    // load elastic parameters
    RealType param1 = strtod( argv[3], NULL ), param2 = strtod( argv[4], NULL ), param3 = strtod( argv[5], NULL );
    ConfiguratorType::ArrayType weight( grid );
    if ( argc == 8 )
      weight.load( argv[6] );
    else
      weight.setAll( 1. );
    MaterialLawType elasticEnergyDensity( param1, param2, param3 );

    // compute elastic deformation energy density
    ConfiguratorType::ArrayType energy( grid );
    typedef qc::WeightedHyperelasticEnergyDensity<ConfiguratorType,MaterialLawType> WeightedEnergyDensityType;
    WeightedEnergyDensityType weightedDensity( elasticEnergyDensity, grid, weight );
    qc::HyperelasticEnergy<ConfiguratorType,WeightedEnergyDensityType>( grid, weightedDensity ).applyAddIntegrand( displacement, energy );

    // save result
    energy.save( argv[argc-1], qc::PGM_FLOAT_ASCII );
  } catch ( aol::Exception &el ) {
    el.dump();
  }
}
