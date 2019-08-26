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

#ifndef __GRAPHTOTRIANGMESH_H
#define __GRAPHTOTRIANGMESH_H

#include <scalarArray.h>
#include <triangMesh.h>


namespace qcsm {

template <typename RealType >
class GraphToTriangMesh : public aol::Op< qc::ScalarArray<RealType, qc::QC_2D>, aol::TriangMesh<RealType> > {
  RealType _scalingFactor;
public:
  GraphToTriangMesh( ) : aol::Op< qc::ScalarArray<RealType, qc::QC_2D>, aol::TriangMesh<RealType> > () {
    _scalingFactor = 1.;
  }

  virtual ~GraphToTriangMesh() {}

  void apply ( const qc::ScalarArray<RealType, qc::QC_2D>& Arg, aol::TriangMesh<RealType> &Dest ) const;

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_2D>&, aol::TriangMesh<RealType>& ) const {
    throw aol::UnimplementedCodeException ("qcsm::GraphToTriangMesh::applyAdd does not make sense", __FILE__, __LINE__ );
  }

  void setScalingFactor ( const RealType scalingFactor ) { _scalingFactor = scalingFactor; }
};

}
#endif
