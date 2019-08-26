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

#include <segment3d.h>
#include <linearSmoothOp.h>

namespace eik {

Segment3d::Segment3d ( qc::ScalarArray<float, qc::QC_3D> *Original )
    : Eikonal3d(),
    greysigma ( 100. ),
    gradsigma ( 200. ),
    MAX_SEEDPOINTS ( 1024 ),
    NUM_IDF_POINTS ( 1025 ) // 2^k + 1 bc boundaries used
{
  // createLookUpTables( ); called in constructor of base class
  nop = Original->getNumX();
  image = new qc::ScalarArray<float, qc::QC_3D> ( nop, nop, nop );
  ( *image ) = ( *Original ); // we want a complete copy here. can copy constructor be used?

  ( *image ) /= image->getMaxValue(); // to scale grey values to [0; 1]

  gvAtSeedPoints = new float[MAX_SEEDPOINTS];
  gvAtStopSeed = new float[MAX_SEEDPOINTS];

  numOfSeedPoints = 0;
  numOfStopSeed = 0;

  speed = new aol::Vector<float> ( NUM_IDF_POINTS );

  for ( int i = 0; i < MAX_SEEDPOINTS; i++ ) {
    gvAtSeedPoints[i] = -1.0;
    gvAtStopSeed[i] = -1.0;
  }

  updateGreyspeed();
}

Segment3d::~Segment3d() {
  // cerr << "starting destructor Segment3d ... ";
  if ( image ) delete ( image ); // delete is not called if pointers are == NULL
  if ( speed ) delete ( speed );
  if ( gvAtSeedPoints ) delete[] gvAtSeedPoints;
  if ( gvAtStopSeed ) delete[] gvAtStopSeed;
  //  cerr << "done." << endl;
}

void Segment3d::setGreySigma ( float Greysigma ) {
  greysigma = Greysigma;
}

void Segment3d::setGradSigma ( float Gradsigma ) {
  gradsigma = Gradsigma;
}

void Segment3d::newSeedPoint ( int X, int Y, int Z ) {
  addSeedGV ( image->get ( X, Y, Z ) );
  addSeedPoint ( X, Y, Z );
}

void Segment3d::addDrySeedPoint ( int X, int Y, int Z ) {
  addSeedPoint ( X, Y, Z );
}

void Segment3d::addStopSeed ( int X, int Y, int Z ) {
  addStopGV ( image->get ( X, Y, Z ) );
}

void Segment3d::addSeedGV ( float Gv ) {
  // store grey value where speed is to be set 1
  bool isnew = true;
  if ( Gv < 0. ) Gv = 0.; if ( Gv > 1. ) Gv = 1.;

  for ( int i = 0; i < numOfSeedPoints; i++ )
    if ( gvAtSeedPoints[i] == Gv ) isnew = false;
  if ( numOfSeedPoints == MAX_SEEDPOINTS ) isnew = false;

  if ( isnew ) {
    gvAtSeedPoints[numOfSeedPoints] = Gv;
    numOfSeedPoints++;
    cerr << "Added new seedPoint with " << numOfSeedPoints << ". grey value " << Gv << endl;
  } else {
    cerr << "Grey value " << Gv << " already stored as seed point. " << endl;
  }

  updateGreyspeed();
}

void Segment3d::addStopGV ( float Gv ) {
  // store grey value where speed is to be set 0
  bool isnew = true;
  if ( Gv < 0. ) Gv = 0.; if ( Gv > 1. ) Gv = 1.;

  for ( int i = 0; i < numOfStopSeed; i++ )
    if ( gvAtStopSeed[i] == Gv ) isnew = false;
  if ( numOfStopSeed == MAX_SEEDPOINTS ) isnew = false;

  if ( isnew ) {
    gvAtStopSeed[numOfStopSeed] = Gv;
    numOfStopSeed++;
    cerr << "Added new StopSeed with " << numOfStopSeed << ". grey value " << Gv << endl;
  } else {
    cerr << "Grey value " << Gv << " already stored as stopSeed. " << endl;
  }

  updateGreyspeed();
}

void Segment3d::resetSeedGV( ) {
  numOfSeedPoints = 0;
  for ( int i = 0; i < MAX_SEEDPOINTS; i++ ) {
    gvAtSeedPoints[i] = -1.;
  }
  cerr << "gv@SeedPoints reset" << endl;
  updateGreyspeed();

}

void Segment3d::resetStopGV( ) {
  numOfStopSeed = 0;
  for ( int i = 0; i < MAX_SEEDPOINTS; i++ ) {
    gvAtStopSeed[i] = -1.;
  }
  cerr << "gv@StopSeed reset" << endl;
  updateGreyspeed();
}

void Segment3d::updateGreyspeed( ) {
  // update discrete values for speed,
  // i. e. evaluate grey value dependent factor of speed at NUM_IDF_POINTS  equidistant values in [0; 1]

  cerr << "updating greyspeed with greysigma = " << greysigma  << " and "
  << numOfSeedPoints << " seedPoint(s), " << numOfStopSeed << " stopSeed value(s) ..." ;

  for ( int g = 0; g < NUM_IDF_POINTS; g++ ) {
    float gv_h = static_cast<float> ( g ) / ( NUM_IDF_POINTS - 1 ); // _h here, i. e. for this grey value
    float speed_here = 1.;

    if ( numOfSeedPoints > 0 ) {
      float temp = 1.;
      float speed_m =  1.f / ( 1.f + greysigma * ( gv_h - gvAtSeedPoints[0] ) * ( gv_h - gvAtSeedPoints[0] ) );

      // first find maximum of functions of type 1 / (1 + (\delta greyvalue)^2)
      for ( int i = 1; i < numOfSeedPoints; i++ ) {
        float gv = gvAtSeedPoints[i];
        temp =  1.f / ( 1.f + greysigma * ( gv_h - gv ) * ( gv_h - gv ) );
        if ( temp > speed_m ) // found greater value
          speed_m = temp;
      }

      speed_here = speed_m;
    }

    if ( numOfStopSeed > 0 ) {
      float temp = 1.;
      float speed_m = ( greysigma * ( gv_h - gvAtStopSeed[0] ) * ( gv_h - gvAtStopSeed[0] ) ) / ( 1 +  greysigma * ( gv_h - gvAtStopSeed[0] ) * ( gv_h - gvAtStopSeed[0] ) );
      for ( int i = 1; i < numOfStopSeed; i++ ) {
        float gv = gvAtStopSeed[i];
        temp =  ( greysigma * ( gv_h - gv ) * ( gv_h - gv ) ) / ( 1.f +  greysigma * ( gv_h - gv ) * ( gv_h - gv ) );
        if ( temp < speed_m ) speed_m = temp;
      }

      speed_here *= speed_m;
    }

    ( *speed ) [g] = speed_here;
  }
  cerr << "done (4info: gradSigma is " << gradsigma << ")" << endl;
}

float Segment3d::getSpeed ( int X, int Y, int Z, int /*Filter*/ ) {
  // Filter is for backwards compatibility only
  // returns speed at position.

  // find a different way to update parameters.
  //   static float old_greysigma = -1.0;

  float s = 0.0;
  //   if (greysigma != old_greysigma){
  //     // if greysigma has been changed since last call of function greyspeed needs to be updated.
  //     updateGreyspeed();
  //     old_greysigma = greysigma;
  //   }

  // grey value dependent factor comes from greyspeed class (linear piecewise interpolation)
  // using median of 1-environment
  s = speed->interpolateInRange ( image->getMedianFilterValue ( X, Y, Z ), 0.0f, 1.0f );


  // gradient dependent factor is done here:
  // s *= ( 1.f / ( 1.f + gradsigma * image->gradSqr ( X, Y, Z ) ) ); // gradSqr actually was two times the gradient
  aol::Vec3<float> curGrad;
  image->gradientFD( X, Y, Z, curGrad );
  s *= ( 1.f / ( 1.f + gradsigma * ( 4 * ( curGrad * curGrad ) ) ) );
  return ( s );
}

void Segment3d::makeTrialNodeActive ( TrialNode<float> &Node, int /*SweepMode*/ ) {
  // for backwards compatiblity
  Eikonal3d::makeTrialNodeActive ( Node );
}


void Segment3d::saveFinalTimeField() {
  try {
    qc::ScalarArray<float, qc::QC_3D> saveField ( 129, 129, 129 );
    qc::ScalarArray<float, qc::QC_3D> tempField ( 129, 129, 129 );
    qc::GridDefinition grid129 ( 7, qc::QC_3D );

    qc::LinearSmoothOp<float> smoothen;

    // if (size > 129): use qc::RestrictOp to calculate a 129^3 timefield
    if ( GRID_DEPTH > 7 ) {
      cout << "dim > 2^7. Using qc::RestrictOp to shrink array ... " << endl;
      qc::GridDefinition origGrid ( GRID_DEPTH, qc::QC_3D );

      qc::RestrictOp<float, qc::STD_QUOC_RESTRICT> shrink ( grid129, origGrid );
      shrink.apply ( ( *finalTimeField ), tempField );
    } else {
      tempField = ( *finalTimeField );
    }

    cout << "Using qc::LinearSmoothOp to smoothen image ... " << endl;
    smoothen.setCurrentGrid ( grid129 );
    smoothen.setSigma ( 1.f / ( 1 << GRID_DEPTH ) );
    smoothen.apply ( tempField, saveField );

    saveField /= saveField.getMaxValue();

    char filename[ 1024 ];
    cout << "filename: ";
    cin >> filename;
    saveField.save ( filename, qc::PGM_FLOAT_BINARY );
    cout << "file saved successfully" << endl;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}

}   // end namespace eik
