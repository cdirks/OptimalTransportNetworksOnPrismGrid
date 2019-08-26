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

#include <levelsetToTriangMesh.h>
#include <tpCFEGrid.h>
#include <tpCFEUtils.h>

namespace qcsm {

template <typename RealType>
void LevelsetToTriangMesh<RealType>::apply ( const qc::ScalarArray<RealType, qc::QC_3D>& Arg, aol::TriangMesh<RealType> &Dest ) const {
  Dest.clear();

  // create tpCFEGrid
  typedef tpcfe::CFEGrid< RealType, tpcfe::CFE_NONE > GridType;
  GridType grid ( Arg.getSize() );
  grid.setDomainFrom ( Arg );
  grid.detectVirtualNodes();

  // create interface triangulation
  tpcfe::CFEInterfaceTriangulationGenerator< GridType, RealType > itg ( grid );
  itg.determineInterfaceTriangulation();

  typedef typename tpcfe::CFEInterfaceTriangulationGenerator<GridType, RealType>::VertexVectorType VertexVectorType;
  typedef typename tpcfe::CFEInterfaceTriangulationGenerator<GridType, RealType>::TriangVectorType TriangVectorType;
  const VertexVectorType& vertices  = itg.getInterfaceVertexVectorRef();
  const TriangVectorType& triangles = itg.getInterfaceTriangVectorRef();

  Dest.reserve ( itg.getNumberOfVertices(), itg.getNumberOfTriangles() );

  for ( typename VertexVectorType::const_iterator it = vertices.begin(); it != vertices.end(); ++it )
    Dest.pushBackVertex ( ( *it ) );

  for ( typename TriangVectorType::const_iterator it = triangles.begin(); it != triangles.end(); ++it )
    Dest.pushBackTriang ( ( *it ) );

}

template <typename RealType>
void LevelsetToTriangMesh<RealType>::addVertexData
( aol::TriangMesh<RealType> &mesh, const qc::ScalarArray<RealType, qc::QC_3D>& vertexData, string vertexDataDescr ) const {

  int newIndex = mesh.getNumVertexDataVectors();
  if (vertexDataDescr == "")
    mesh.addVertexData ();
  else
    mesh.addVertexData ( vertexDataDescr );

  for (int i = 0; i < mesh.getNumVertices(); ++i)
    mesh.getVertexData(newIndex)[i] = vertexData.interpolate_on01 ( mesh.getVertex(i) );
}

template class LevelsetToTriangMesh<float>;
template class LevelsetToTriangMesh<double>;
template class LevelsetToTriangMesh<long double>;

}
