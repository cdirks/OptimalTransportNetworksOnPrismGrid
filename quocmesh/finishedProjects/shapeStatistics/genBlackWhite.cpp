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

#include <aol.h>
#include <configurators.h>
#include <deformations.h>

typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;

int main( int argc, char *argv[] ) {

  try {

    // read in file names of parameter files
    if ( argc != 4 ) {              // wrong syntax
      cerr << "Wrong syntax. Use: genBlackWhite <grid depth> <filename of levelset function> <filename of image>" << endl;
      return EXIT_FAILURE;
    } else {
      qc::GridDefinition grid( atoi(argv[1]), ConfiguratorType::Dim );
      ConfiguratorType::ArrayType levelset( grid );
      levelset.load( argv[2] );
      levelset /= - grid.H() * 2;
      levelset.addToAll( .5 );
      for ( int i = 0; i < levelset.size(); i++ )
        if ( levelset[i]<0 )
          levelset[i] = 0;
        else if ( levelset[i]>1 )
          levelset[i] = 1;
      levelset *= 255;

      // if needed, deform the levelset
      // generate a deformation
      qc::DataGenerator< ConfiguratorType > deformationGenerator ( grid );
      aol::MultiVector<RealType> disp( ConfiguratorType::Dim, levelset.size() );
      // shear deformation
      //deformationGenerator.generateShearDeformation ( .5, disp );
      /*// mirror image
      deformationGenerator.generateLineLevelset ( disp[0], 0 );
      disp *= -2.;
      // inversely deform the levelset
      qc::TransformFunction<RealType,ConfiguratorType::Dim> inverseDeformation( grid );
      inverseDeformation.setDeformation( disp );
      aol::MultiVector<RealType> levelsetMultiVector;
      levelsetMultiVector.appendReference( levelset );
      aol::MultiVector<RealType> aux( levelsetMultiVector );
      inverseDeformation.apply( aux, levelsetMultiVector );
      */
      levelset.save( argv[3], qc::PGM_UNSIGNED_CHAR_BINARY );
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
