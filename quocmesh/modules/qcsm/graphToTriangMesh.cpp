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

#include <graphToTriangMesh.h>
#include <scalarArray.h>
#include <configurators.h>
#include <progressBar.h>

namespace qcsm {

template <typename RealType>
void GraphToTriangMesh<RealType>::apply ( const qc::ScalarArray<RealType, qc::QC_2D>& Arg, aol::TriangMesh<RealType> &Dest ) const {
  Dest.clear();

  if ( Arg.getNumX() != Arg.getNumY() )
    throw aol::Exception ( "Sorry, only for quadratic arrays.", __FILE__, __LINE__ );

  int N = Arg.getNumX ();
  int d = qc::logBaseTwo ( N );
  qc::GridDefinition grid ( d, qc::QC_2D );

  // Configurator
  qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 7> > config ( grid );

  // reserving enough space to avoid reallocation
  Dest.reserve ( grid.getNumberOfNodes(), 2*grid.getNumberOfElements() );

  // get the gridDefinition and define some variables
  qc::GridDefinition::OldFullElementIterator it;
  RealType n = static_cast<RealType> ( N );
  aol::Vec3<int> intCoords;
  aol::Vec3<RealType> realCoords;
  aol::Vec3<int> ind;


  // traverse all nodes and store all vertices
  cerr << "\nThe graph will be scaled with factor " << aol::color::red << _scalingFactor << aol::color::reset << endl;
  aol::ProgressBar<> pbV ( "Extracting Vertices: " );
  pbV.start ( grid.getNumberOfNodes() );
  for ( int i = 0; i < grid.getNumberOfNodes(); i++ ) {
    grid.indexToCoordGlobal ( i, intCoords );                                 // get the grid-coordinates
    for ( int j = 0; j < 2; j++ ) realCoords[j] = static_cast<RealType> ( intCoords[j] ) / n;
    realCoords[2] = _scalingFactor * Arg[i];                                  // get the height of the graph
    Dest.pushBackVertex ( realCoords );                                       // store the vertex
    pbV++;
  }
  pbV.finish();

  // now traverse all elements and split each element into two triangles
  aol::ProgressBar<> pbT ( "Extracting Triangles:" );
  pbT.start ( grid.getNumberOfElements() );
  int indList1[] = { 0, 1, 3 };       // lower right triangle of the element
  int indList2[] = { 0, 3, 2 };       // upper left one

  for ( it = grid.begin(); it != grid.end(); ++it ) {
    // first triangle (right lower one)
    for ( int i = 0; i < 3; i++ ) ind[i] = config.localToGlobal ( *it, indList1[i] );
    Dest.pushBackTriang ( ind );                                                   // store the facet itself

    // second triangle (right lower one)
    for ( int i = 0; i < 3; i++ ) ind[i] = config.localToGlobal ( *it, indList2[i] );
    Dest.pushBackTriang ( ind );                                                   // store the facet itself
    pbT++;
  }
  pbT.finish();

}

template class GraphToTriangMesh<float>;
template class GraphToTriangMesh<double>;

}
