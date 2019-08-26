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

#include <smallMat.h>
#include <vec.h>
#include <matrix.h>
#include <scalarArray.h>
#include <iterators.h>
#include <progressBar.h>
#include "VoigtTensor.h"

typedef double RealType;


RealType computeOrthotropyViolation ( const VoigtElastOp<RealType> &tensor ) {

  // RiOdKa96 approach: minimize sum of squared undesired entries divided by sum of squared desired entries
  // here weighted by the number of entries in the non-Voigtized tensor they represent.
  RealType numer = 0., denom = 0.;

  const int myONE = 1, myTWO = 2, myFOUR = 4;

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      numer +=              myTWO  * aol::Sqr ( tensor.get ( i+0, j+3 ) );
      numer +=              myTWO  * aol::Sqr ( tensor.get ( i+3, j+0 ) );
      numer += ( i != j ) * myFOUR * aol::Sqr ( tensor.get ( i+3, j+3 ) );

      denom +=              myONE  * aol::Sqr ( tensor.get ( i+0, j+0 ) );
      denom += ( i == j ) * myFOUR * aol::Sqr ( tensor.get ( i+3, j+3 ) );
    }
  }

  return ( numer / denom );

}

RealType orthotropyViolation ( const VoigtElastOp<RealType> &tensor, const aol::Matrix33<RealType> &rotMat ) {
  VoigtElastOp<RealType> rotatedTensor;
  applyRotation ( rotMat, tensor, rotatedTensor );
  return ( computeOrthotropyViolation ( rotatedTensor ) );
}

// objective function penalizing deviation of 6x6 tensor from orthotropic one; angles in rad
RealType orthotropyViolation ( const VoigtElastOp<RealType> &tensor, const aol::Vec3<RealType> &angles ) {
  VoigtElastOp<RealType> rotatedTensor;
  applyBackwardRotation ( angles, tensor, rotatedTensor );
  return ( computeOrthotropyViolation ( rotatedTensor ) );
}

// angles in rad
void angleSearch ( const VoigtElastOp<RealType> &tensor, const aol::Vec3<RealType> &minAngles, const aol::Vec3<RealType> &maxAngles, const int N, aol::Vec3<RealType> &optAngles ) {
  qc::ScalarArray<RealType, qc::QC_3D> myValues ( N, N, N );

  cerr << "Angle search in range (" << minAngles << "), (" << maxAngles << "), discretized in " << N << "." << endl;

  aol::ProgressBar<> pb( "Trying angles" );
  pb.start ( myValues.size() );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( aol::Vec3<short>( 0, 0, 0 ), aol::Vec3<short>( N, N, N ) ); bit.notAtEnd(); ++bit, pb++ ) {
    aol::Vec3<RealType> angles ( minAngles[0] + (*bit)[0] / (N-1.) * ( maxAngles[0] - minAngles[0] ), minAngles[1] + (*bit)[1] / (N-1.) * ( maxAngles[1] - minAngles[1] ), minAngles[2] + (*bit)[2] / (N-1.) * ( maxAngles[2] - minAngles[2] ) );
    const RealType value = orthotropyViolation ( tensor, angles );
    myValues.set ( *bit, value );
  }
  cerr << endl;

  qc::CoordType minPos;
  RealType minVal = aol::NumberTrait<RealType>::Inf;
  for ( qc::RectangularIterator<qc::QC_3D> bit ( myValues ); bit.notAtEnd(); ++bit ) {
    if ( myValues.get( *bit ) < minVal ) {
      minVal = myValues.get ( *bit );
      minPos = *bit;
    }
  }

  optAngles[0] =  minAngles[0] + minPos[0] / (N-1.) * ( maxAngles[0] - minAngles[0] );
  optAngles[1] =  minAngles[1] + minPos[1] / (N-1.) * ( maxAngles[1] - minAngles[1] );
  optAngles[2] =  minAngles[2] + minPos[2] / (N-1.) * ( maxAngles[2] - minAngles[2] );

  cerr << "Minimum " << minVal << " found at " << minPos << ", corresponding to angles " << optAngles << " (in rad)" << endl;
}

// search on nested intervals
void iteratedAngleSearch ( const VoigtElastOp<RealType> &tensor, const aol::Vec3<RealType> &initMinAngles, const aol::Vec3<RealType> &initMaxAngles, const int initN, const int laterN, const int numIter, aol::Vec3<RealType> &optAngles ) {
  aol::Vec3<RealType> currentRange = initMaxAngles - initMinAngles;

  angleSearch ( tensor, initMinAngles, initMaxAngles, initN, optAngles );
  currentRange /= initN;

  for ( int i = 1; i < numIter; ++i ) {
    angleSearch ( tensor, optAngles - 0.5 * currentRange, optAngles + 0.5 * currentRange, laterN, optAngles );
    currentRange /= laterN;
  }
}


void findAndPrintOrthoParams ( const VoigtElastOp<RealType> &tensor ) {

  VoigtElastOp<RealType> compliance;
  compliance.makeInverse ( tensor );
  compliance.printTensor ( cout );

  const RealType
    Exx = 1. / compliance.get ( 0, 0 ),
    Eyy = 1. / compliance.get ( 1, 1 ),
    Ezz = 1. / compliance.get ( 2, 2 ),
    Gyz = 1. / compliance.get ( 3, 3 ),
    Gzx = 1. / compliance.get ( 4, 4 ),
    Gxy = 1. / compliance.get ( 5, 5 ),
    nuxy = - compliance.get ( 1, 0 ) / compliance.get ( 0, 0 ),
    nuyx = - compliance.get ( 0, 1 ) / compliance.get ( 1, 1 ),
    nuxz = - compliance.get ( 2, 0 ) / compliance.get ( 0, 0 ),
    nuzx = - compliance.get ( 0, 2 ) / compliance.get ( 2, 2 ),
    nuyz = - compliance.get ( 2, 1 ) / compliance.get ( 1, 1 ),
    nuzy = - compliance.get ( 1, 2 ) / compliance.get ( 2, 2 );

    cout << "Exx = " << aol::detailedFormat ( Exx ) << endl
         << "Eyy = " << aol::detailedFormat ( Eyy ) << endl
         << "Ezz = " << aol::detailedFormat ( Ezz ) << endl
         << "Gyz = " << aol::detailedFormat ( Gyz ) << endl
         << "Gzx = " << aol::detailedFormat ( Gzx ) << endl
         << "Gxy = " << aol::detailedFormat ( Gxy ) << endl
         << "nuxy = " << aol::detailedFormat ( nuxy ) << endl
         << "nuyx = " << aol::detailedFormat ( nuyx ) << endl
         << "nuxz = " << aol::detailedFormat ( nuxz ) << endl
         << "nuzx = " << aol::detailedFormat ( nuzx ) << endl
         << "nuyz = " << aol::detailedFormat ( nuyz ) << endl
         << "nuzy = " << aol::detailedFormat ( nuzy ) << endl << endl;

    cout << "Anisotropy x/z = " << aol::detailedFormat ( Exx / Ezz ) << endl
         << "Anisotropy y/z = " << aol::detailedFormat ( Eyy / Ezz ) << endl
         << "Anisotropy x/y = " << aol::detailedFormat ( Exx / Eyy ) << endl;
}

int main ( int argc, char** argv ) {
  if ( argc != 3 ) {
    cerr << "usage: findBestRotationToOrthotropicTensor infile rangeAngle" << endl;
    return ( EXIT_FAILURE );
  }

  try {

    VoigtElastOp<RealType> inTensor( argv[1] );
    const int rangeAngle = atoi ( argv[2] );

    cout << "Loaded Tensor:" << endl;
    inTensor.print ( cout );

    cout << "Loaded Tensor (once again):" << endl;
    inTensor.printLatexTensor ( cout );


    cerr << "Initial orthotropy violation: " << orthotropyViolation( inTensor, aol::Vec3<RealType>() ) << endl;

    if ( rangeAngle > 0 ) {

      VoigtElastOp<RealType> backRotatedTensor;

#if 0
      aol::Vec3<RealType> myAngles (0, 0, 0 );
      myAngles *= aol::NumberTrait<RealType>::pi / 180;

      applyBackwardRotation ( myAngles, inTensor, backRotatedTensor );
      cerr << "orthotropy violation: " << orthotropyViolation( backRotatedTensor, aol::Vec3<RealType>() ) << endl;
#else


      aol::Vec3<RealType> optAngles, minAngles ( -rangeAngle, -rangeAngle, -rangeAngle ), maxAngles ( rangeAngle, rangeAngle, rangeAngle );
      minAngles *= aol::NumberTrait<RealType>::pi / 180;
      maxAngles *= aol::NumberTrait<RealType>::pi / 180;

      iteratedAngleSearch ( inTensor, minAngles, maxAngles, aol::Max ( 201, 2 * rangeAngle + 1 ), 81, 8, optAngles );

      applyBackwardRotation( optAngles, inTensor, backRotatedTensor );

      backRotatedTensor.print( cerr );

      aol::Matrix33<RealType> rotMat;
      getRotationMatrix( optAngles, rotMat );

      cerr << "Final orthotropy violation: " << orthotropyViolation( inTensor, optAngles ) << endl;

      optAngles *= 180 / aol::NumberTrait<RealType>::pi;

      cerr << "optAngles = ( "
           << aol::longScientificFormat ( optAngles[0] ) << ", "
           << aol::longScientificFormat ( optAngles[1] ) << ", "
           << aol::longScientificFormat ( optAngles[2] ) << ") (degrees)." << endl;

      cerr << "rotation Matrix = " << endl << rotMat << endl;
#endif

      cout << "BackRotatedTensor:" << endl;
      backRotatedTensor.printTensor ( cout );

      cout << "BackRotatedTensor (once again):" << endl;
      backRotatedTensor.printLatexTensor ( cout );

    } else {

      // inTensor.restrictToOrtho();
      cout << "E, G, nu for tensor restricted to orthotropy entries:" << endl;
      findAndPrintOrthoParams ( inTensor );

    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}

/* Input file should be of the form:
  9.0  5.0  2.0  0.0  0.0  0.0
  5.0 11.0  3.0  0.0  0.0  0.0
  2.0  3.0 12.0  0.0  0.0  0.0
  0.0  0.0  0.0  3.5  0.0  0.0
  0.0  0.0  0.0  0.0  2.5  0.0
  0.0  0.0  0.0  0.0  0.0  1.8
*/
