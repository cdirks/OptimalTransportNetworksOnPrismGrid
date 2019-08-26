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

#include "atreg3d.h"
#include <configurators.h>

typedef float RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
/*
    typedef double TestRealType;
    aol::Vector<TestRealType> vec(1001);

    vec.setAll(1.);
    aol::StopWatch watch;
    watch.start();
    TestRealType temp = 0.;
    for( int i = 0; i < 1000000; i++ ){
      temp += vec.normSqr();
    }
    watch.stop();
    cerr << temp << endl;
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
*/
/*
    qc::GridDefinition Grid( 8, ConfType::Dim );
    ConfType::ArrayType image( Grid );
    image.setAll(1.);


    const int n = Grid.getWidth();
    aol::StopWatch watch;
    watch.start();

    float temp = 0.;
    for( int j = 0; j < n; j ++ ){
      for( int i = 0; i < n; i ++ ){
        for( int k = 0; k < n; k ++ ){
          temp += image.interpolate( i, j, k);
        }
      }
    }
    watch.stop();
    cerr << temp;
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
*/
/*
    //---------------------------------------------------------------
    qc::GridDefinition Grid( 4, ConfType::Dim );
    aol::MultiVector<RType> phi( ConfType::Dim, Grid.getNumberOfNodes() );
    for( int i = 0; i < ConfType::Dim; i++){
      phi[i].setAll( Grid.H() );
    }
    qc::TransformFunction<RType, ConfType::Dim> transform( Grid );
    transform.setDeformation( phi );
    ConfType::ArrayType image( Grid );

    const RType h= Grid.H();

    qc::GridDefinition::OldFullNodeIterator fnit;
    qc::FastILexMapper<ConfType::Dim> mapper(Grid);
    for(fnit = Grid.begin(); fnit != Grid.end(); fnit++){
      bool check = true;
      for( int i = 0; i < ConfType::Dim; i++){
        if (!( ( (*fnit)[i] * h > 0.3 ) && ( (*fnit)[i] * h < 0.7 )  ))
          check = false;
      }
      if ( check )
        image[ mapper.getGlobalIndex(*fnit) ] = 10+(*fnit)[0]+(*fnit)[1];
      else
        image[ mapper.getGlobalIndex(*fnit) ] = -(*fnit)[0]-(*fnit)[1];
    };
    image -= image.getMinValue();
    image /= image.getMaxValue();
    qc::writeImage<RType>( Grid, image, "test/image", 1 );

    ConfType::ArrayType tmp( Grid );
    transform.apply( image, tmp );
    qc::writeImage<RType>( Grid, tmp, "test/test", 1 );
    aol::callSystemPauseIfNecessaryOnPlatform();
    //---------------------------------------------------------------
*/
/*
    qc::ScalarArray<RType, qc::QC_3D> a(256, 256, 176);
    ifstream in( "f:/bPlath.img", ios::binary);
    a.loadRaw( in, 10, 256, 256, 176);
    a.save( "a.raw", 5);
*/
/*
    qc::ScalarArray<RType, qc::QC_3D> a("results/fine_checkbox_04_000.raw");
    a.save( "a.raw", 5);
*/

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
        sprintf( parameterfilename, "ParameterFiles/Reg_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/Reg_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    aol::StopWatch watch;
    watch.start();

    //qc::MIRegistrationMultilevelDescent<ConfType> mld( parser );
    atRegMultilevelDescent<ConfType> mld( parser );
    mld.solve();

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
