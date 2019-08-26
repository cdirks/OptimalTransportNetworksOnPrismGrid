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

#ifndef __TPCFETOVTK_H
#define __TPCFETOVTK_H

#include <parameterParser.h>
#include <tpCFEGrid.h>

#include <vtkIncludes.h>

#include "polyDataGroup.h"

/* class PolyDataGroup; */
class vtkPolyData;

namespace tpcfe {

/** This class is used for loading data from a tpcfe elasticity computation in 3D and preparing it in such a way that it can be displayed using vtkfox
 *  TODO: find a better name ...
 *  \author Schwen
 */
template< typename GridType >
class TpcfeToVTK {
public:
  typedef typename GridType::RealType RealType;

  TpcfeToVTK ( char* filename, PolyDataGroup *pdgr ) : parser ( filename ), _polyDataGroup ( pdgr ) {
    cerr << "TpcfeToVTK: Loading parameter file " << filename << endl;
  }

public:
  vtkPolyData* makePolyData ( );

protected:
  aol::ParameterParser parser;
  PolyDataGroup* _polyDataGroup;

private:
  TpcfeToVTK ( ) : parser ( "" ), _polyDataGroup ( NULL ) {
    cerr << "No parameter file specified" << endl;
  }

};


template< typename RealType >
class LevelsetToVTK {

public:
  LevelsetToVTK ( qc::ScalarArray<RealType, qc::QC_3D> &levelset ) : _levelset ( levelset ) {}

public:
  vtkPolyData* makePolyData ( );
  vtkPolyData* makePolyData ( const aol::Vec3<int> &clip_lower, const aol::Vec3<int> &clip_upper );

protected:
  qc::ScalarArray<RealType, qc::QC_3D> &_levelset;
  // maybe allow other value, then need deep copy for modifying...

};

}

#endif
