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

#ifndef __SEGMENT3D_H
#define __SEGMENT3D_H

#include <eikonal3d.h>
#include <math.h>
#include <linearSmoothOp.h>

namespace eik {

class Segment3d: public Eikonal3d {
protected:

  float greysigma, gradsigma;
  const int MAX_SEEDPOINTS;
  const int NUM_IDF_POINTS;
  int numOfSeedPoints, numOfStopSeed; // not exactly, but number of Points with different grey values
  float *gvAtSeedPoints, *gvAtStopSeed;

  aol::Vector<float> *speed;
  int nop; // necessary?
  qc::ScalarArray<float, qc::QC_3D> *image;

public:
  Segment3d ( qc::ScalarArray<float, qc::QC_3D> *Original );

  ~Segment3d( );

  void resetSeedGV( );
  void resetStopGV( );

  void setGreySigma ( float Greysigma );
  void setGradSigma ( float Gradsigma );

  void addDrySeedPoint ( int X, int Y, int Z );
  void addStopSeed ( int X, int Y, int Z );
  void addSeedGV ( float Gv );
  void addStopGV ( float Gv );
  void updateGreyspeed( );

  void newSeedPoint ( int X, int Y, int Z );

  float getSpeed ( int X, int Y, int Z , int Filter = 0 );

  void makeTrialNodeActive ( TrialNode<float> &Node, int SweepMode );

  void saveFinalTimeField();
};

}   // end namespace eik

#endif
