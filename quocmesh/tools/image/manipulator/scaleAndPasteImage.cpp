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

/** \file
 *  \brief pastes a (3D) image into another one
 *
 * Scales an input image and pastes it into new image of same size with the origin lying at x0,y0,z0.
 * Undefined areas are filled up with fillValue.
 *
 *  \author Wirth
 */

#include <aol.h>
#include <scalarArray.h>
#include <generator.h>
#include <deformations.h>
#include <configurators.h>

typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;

int main ( int argc, char **argv ) {

  try {
    if ( argc != ( ConfiguratorType::Dim == qc::QC_2D ? 8 : 10 ) ) {
      cerr << aol::color::red << "USAGE: " << argv[0] << " sourceImage destinationImage x0 y0";
      cerr << ( ConfiguratorType::Dim == qc::QC_2D ? "" : " z0" ) << " scaleFactorX scaleFactorY";
      cerr << ( ConfiguratorType::Dim == qc::QC_2D ? "" : " scaleFactorZ" ) << " fillValue" << endl << aol::color::reset;
      return EXIT_FAILURE;
    }

    ConfiguratorType::ArrayType source( argv[1] ), dest( source, aol::STRUCT_COPY ), extendImage( source, aol::STRUCT_COPY );

    ConfiguratorType::InitType grid( source.getSize() );
    aol::MultiVector<RealType> displacement( grid );
    ConfiguratorType::MatType scale;
    ConfiguratorType::VecType shift;
    scale.setIdentity();
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      scale[i][i] /= strtod( argv[argc-1-ConfiguratorType::Dim+i], NULL );
      shift[i] = - strtod( argv[3+i], NULL ) * scale[i][i];
    }
    qc::DataGenerator<ConfiguratorType>( grid ).generateAffineDisplacement( scale, shift, displacement );
    extendImage.setAll( strtod( argv[argc-1], NULL ) );

    qc::DeformImage<ConfiguratorType>( source, grid, dest, displacement, extendImage );

    dest.save( argv[2], qc::PGM_DOUBLE_BINARY );

  } catch ( aol::Exception &el ) {
    el.dump();
  }
}
