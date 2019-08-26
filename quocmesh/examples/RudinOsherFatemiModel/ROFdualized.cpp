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

#include <finiteDifferences.h>
#include <parameterParser.h>
#include <configurators.h>
#include <generator.h>
#include <timestepSaver.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main( int argc, char **argv ) {

  try {
    aol::ParameterParser parser( argc, argv, "ROFdualized.par" );

    aol::StopWatch watch;
    watch.start();

    typedef RType RealType;
    typedef ConfType ConfiguratorType;

    const RealType gamma = parser.getDouble( "gamma" );

    string saveDirectory = parser.getString( "saveDirectory" );
    aol::makeDirectory( saveDirectory.c_str() );
    parser.dumpToFile( "parameter-dump.txt", saveDirectory.c_str() );

    ConfiguratorType::ArrayType u0 ( parser.getString( "input-image" ) );
    ConfiguratorType::InitType grid ( qc::GridSize2d::createFrom ( u0 ) );
    u0.load( parser.getString( "input-image" ).c_str() );
    u0 /= u0.getMaxValue();

    qc::MultiArray<RealType, 2> pOld ( grid );
    qc::MultiArray<RealType, 2> pNew ( grid );
    qc::ScalarArray<RealType, qc::QC_2D> temp ( grid );
    qc::MultiArray<RealType, 2> gradient ( grid );
    const RealType tau = 0.25;
    u0 /= gamma;

    aol::TimestepSaver<RealType> timestepSaver( parser.getInt( "saveOffset"), parser.getInt( "numSaveFirst"), "test", true );
    timestepSaver.setSaveDirectory ( saveDirectory.c_str() );

    for ( int fixpointIterations = 0; fixpointIterations < 1000; fixpointIterations++ ) {
      qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( pOld, temp );
      temp -= u0;
      qc::calculateForwardFDGradient<RealType, ConfiguratorType::Dim> ( temp, gradient );

      for ( int j = 0; j < temp.size(); ++j ) {
        aol::Vec2<RealType> gradientVec ( gradient[0][j], gradient[1][j] );
        const RealType gradientVecNorm = gradientVec.norm();
        for ( int i = 0; i < 2; ++i ) {
          pNew[i][j] = ( pOld[i][j] + tau * gradientVec[i] ) / ( 1 + tau * gradientVecNorm );
        }
      }
      pOld -= pNew;

      // If the change is small enough, we consider the gradient descent to have converged.
      if ( pOld.norm() < 0.01 ) {
        cerr << "Stopping after " << fixpointIterations + 1 << " iterations.\n";
        break;
      }

      pOld = pNew;
    }

    qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( pOld, temp );
    temp *= gamma;
    u0 *= gamma;
    u0 -= temp;
    timestepSaver.savePNG ( u0, grid );

    watch.stop();
    watch.printReport( cerr );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
