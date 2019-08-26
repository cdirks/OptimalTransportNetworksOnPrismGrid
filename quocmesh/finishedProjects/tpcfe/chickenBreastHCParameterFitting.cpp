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

#include <vec.h>
#include <randomGenerator.h>
#include <scalarArray.h>

static const float thetaShiftMin   = -900.0;
static const float thetaShiftRange = 1200.0;

template< typename RealType >
RealType scaledTimeInterpolation ( const aol::Vector<RealType> &vec, const RealType tau, const RealType t ) {
  if ( t <= 0 )
    return ( vec[0] );
  else if ( ( t / tau ) >= vec.size() - 1 )
    return ( vec[ vec.size() - 1 ] );
  else
    return ( vec.interpolate ( t / tau ) );
}

template< typename RealType >
RealType mismatch ( const aol::Vector<RealType> &measured, const aol::Vector<RealType> &simulated, const RealType tauMeasured, const RealType tauSimulated, const RealType theta, const RealType shift, const RealType tauMinForMismatch,  const RealType tauMaxForMismatch ) {
  RealType ret = 0;
  for ( int k = 0; k < measured.size(); ++k ) {
    const RealType time = tauMeasured * k;
    if ( ( time >= tauMinForMismatch ) && ( time <= tauMaxForMismatch ) ) {
      ret += aol::Sqr ( scaledTimeInterpolation ( measured, tauMeasured, time - shift ) - scaledTimeInterpolation ( simulated, tauSimulated, theta * time ) );
    }
  }

  return ( ret );
}


template< typename RealType >
void thetaSearch ( const aol::Vector<RealType> &measured, const aol::Vector<RealType> &simulated, const RealType tauMeasured, const RealType tauSimulated,
                   const RealType tauMinForMismatch,  const RealType tauMaxForMismatch,
                   const RealType minTheta, const RealType maxTheta, const int steTheta, RealType &optTheta, RealType &optTimeShift  ) {

  qc::ScalarArray<RealType, qc::QC_2D> myValues ( steTheta );
  cerr << "Theta search in range " << aol::longScientificFormat ( minTheta ) << ", " << aol::longScientificFormat ( maxTheta ) << ", discretized in " << steTheta << "." << endl;

  for ( int i = 0; i < steTheta; ++i ) {
    for ( int j = 0; j < steTheta; ++j ) {

      const RealType theta = minTheta + i * ( maxTheta - minTheta ) / ( steTheta - 1. );
      const RealType shift = thetaShiftMin + j * thetaShiftRange / ( steTheta - 1.0 );
      const RealType value = mismatch ( measured, simulated, tauMeasured, tauSimulated, theta, shift, tauMinForMismatch, tauMaxForMismatch );
      myValues.set ( i, j, value );
    }
  }

  aol::Vec2<int> minPos ( -1, -1 );
  RealType minVal = + aol::NumberTrait<RealType>::Inf;
  for ( int i = 0; i < steTheta; ++i ) {
    for ( int j = 0; j < steTheta; ++j ) {
      if ( myValues.get(i,j) < minVal ) {
        minPos.set ( i, j );
        minVal = myValues.get ( i, j );
      }
    }
  }

  cerr << minPos << endl;

  if ( !( minPos[1] > 0 ) || !( minPos[1] < steTheta - 1 ) )
    throw aol::Exception ( "minimum at boundary of admissible set", __FILE__, __LINE__ );

  optTheta = minTheta + minPos[0] * ( maxTheta - minTheta ) / ( steTheta - 1. );
  optTimeShift = thetaShiftMin + minPos[1] * thetaShiftRange / ( steTheta - 1.0 );

  cerr << "Minimum " << aol::longScientificFormat ( minVal ) << " found at " << minPos << ", corresponding to " << aol::longScientificFormat ( optTheta ) << ", " << aol::detailedFormat ( optTimeShift ) << endl;

}

template< typename RealType >
void iteratedThetaSearch ( const aol::Vector<RealType> &measured, const aol::Vector<RealType> &simulated, const RealType tauMeasured, const RealType tauSimulated,
                           const RealType tauMinForMismatch,  const RealType tauMaxForMismatch,
                           const RealType initMinTheta, const RealType initMaxTheta, const int steTheta, const int numIter, RealType &optTheta,
                           RealType &optTimeShift ) {

  RealType currentRange = initMaxTheta - initMinTheta;
  thetaSearch ( measured, simulated, tauMeasured, tauSimulated, tauMinForMismatch, tauMaxForMismatch, initMinTheta, initMaxTheta, steTheta, optTheta, optTimeShift );

  currentRange /= ( steTheta - 1 );

  for ( int i = 1; i < numIter; ++i ) {
    thetaSearch ( measured, simulated, tauMeasured, tauSimulated, tauMinForMismatch, tauMaxForMismatch, optTheta - 0.5 * currentRange, optTheta + 0.5 * currentRange, steTheta, optTheta, optTimeShift );
    currentRange /= ( steTheta - 1 );
  }
}


typedef double RealType;


int main ( int argc, char** argv ) {
  if ( argc != 7 ) {
    cerr << "Usage: " << argv[0] << " [measured data file] [computed data file] [measurement time step] [simulated time step] [start time for match] [stop time for match]" << endl;
    return ( EXIT_FAILURE );
  }

  try {

    const char
      *MeasFilename = argv[1],
      *CompFilename = argv[2];

    const RealType
      tauMeas = atof ( argv[3] ),
      tauSimu = atof ( argv[4] ),
      tauMinForMismatch = atof ( argv[5] ),
      tauMaxForMismatch = atof ( argv[6] ),
      theta = 1.0;

    const RealType
      thetaLower = 0.0,
      thetaUpper = 10.0;

    const int
      thetaIters = 10,
      thetaSteps = 1001;

    aol::Vector< RealType > measVec, compVec;

    {
      ifstream measData ( MeasFilename );
      int counter = 0;
      do {
        measVec.resize ( measVec.size() + 1 );
        measData >> measVec[counter];
        ++counter;
      } while ( !measData.eof() );
    }
    {
      ifstream compData ( CompFilename );
      int counter = 0;
      do {
        compVec.resize ( compVec.size() + 1 );
        compData >> compVec[counter];
        ++counter;
      } while ( !compData.eof() );
    }

    if ( !( measVec[ measVec.size() - 1 ] > 1.0e-7 ) || !( compVec[ compVec.size() - 1 ] > 1.0e-7 ) )
      throw aol::Exception ( "Strange value found. Please check input file (must not end with newline)", __FILE__, __LINE__ );

    cerr << measVec << endl << compVec << endl;

    cerr << "Initial mismatch: " << aol::detailedFormat ( mismatch ( measVec, compVec, tauMeas, tauSimu, theta, 0., tauMinForMismatch, tauMaxForMismatch ) ) << endl;

    RealType optTheta = + aol::NumberTrait<RealType>::Inf;
    RealType optTimeShift = -1000;
    iteratedThetaSearch ( measVec, compVec, tauMeas, tauSimu, tauMinForMismatch, tauMaxForMismatch, thetaLower, thetaUpper, thetaSteps, thetaIters, optTheta, optTimeShift );

    cerr << "optimal theta = " << aol::detailedFormat ( optTheta ) << endl;
    cerr << "optimal time shift = " << aol::detailedFormat ( optTimeShift ) << endl;

    // match output, can be piped to file
    for ( int k = 0; k < measVec.size(); ++k )
      cout << aol::detailedFormat ( k * tauMeas ) << " "
           << aol::longScientificFormat ( scaledTimeInterpolation ( measVec, tauMeas, (k*tauMeas) - optTimeShift ) ) << " "
           << aol::longScientificFormat ( scaledTimeInterpolation ( compVec, tauSimu, (optTheta*k*tauMeas) ) ) << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}


#if 0
./chickenBreastHCParameterFitting ~/CBData/20080819V1_measured.dat ~/CBData/20080819V1_ellipsoid_simulated.dat 30 1 > ~/CBData/20080819V1_ellipsoid_match.dat
./chickenBreastHCParameterFitting ~/CBData/20080819V2_measured.dat ~/CBData/20080819V2_ellipsoid_simulated.dat 30 1 > ~/CBData/20080819V2_ellipsoid_match.dat
./chickenBreastHCParameterFitting ~/CBData/20080901V1_measured.dat ~/CBData/20080901V1_ellipsoid_simulated.dat 30 1 > ~/CBData/20080901V1_ellipsoid_match.dat
./chickenBreastHCParameterFitting ~/CBData/20080901V2_measured.dat ~/CBData/20080901V2_ellipsoid_simulated.dat 30 1 > ~/CBData/20080901V2_ellipsoid_match.dat
#endif
