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

#include "triangMeshToVTK.h"

namespace aol {

template<class DataType>
vtkPolyData* TriangMeshToVTK<DataType>::makePolyData ( const TriangMesh<DataType> &mesh, DataType &minC, DataType &maxC, bool copyVertexData ) {

  cerr << mesh.getNumTriangs() << " " << mesh.getNumVertices() << " " << mesh.hasVertexData() << endl;

  vtkPolyData*   pd = vtkPolyData::New();  // make new vtkPolyData to return
  vtkPoints*     p  = vtkPoints::New();    // make the points-obj for the poly-data

  p->SetNumberOfPoints ( mesh.getNumVertices() );
  for ( int i = 0; i < mesh.getNumVertices(); i++ )   // copy the mesh coords into the vtkPoints
  {
    aol::Vec3<DataType> vertex = mesh.getVertex(i);
    p->SetPoint ( i, static_cast<float> ( vertex[0] ), static_cast<float> ( vertex[1] ), static_cast<float> ( vertex[2] ) );
  }

  pd->SetPoints ( p );
  p->Delete();  // pass the vtkPoints to the poly-data

  pd->Allocate ( mesh.getNumTriangs() );  // copy the triangles
  for ( int i = 0; i < mesh.getNumTriangs(); ++i ) {
    aol::Vec3<int> triang = mesh.getTriang(i);
    vtkIdType tmpvert[3] = { triang[0], triang[1], triang[2] };
    pd->InsertNextCell ( VTK_TRIANGLE, 3, tmpvert );
  }

  if ( copyVertexData && mesh.hasVertexData() ) {                        // set scalar data, if any passed
    setPolyDataVertexScalars ( pd, mesh.getVertexData(), minC, maxC );   // (data is deep-copied)
  }

  if ( copyVertexData && mesh.hasTriangData() ) {
    setPolyDataTriangScalars ( pd, mesh.getTriangData(), minC, maxC );
  }

  pd->BuildCells();
  pd->BuildLinks();
  pd->Update();

  return pd;

}

//! copies content of \a v into \a pd's scalar field
//!
//! \warning vtkPolyData has only one scalar field. When you call
//! this function, the previously stored scalar field will get lost.
template<class DataType>
vtkPolyData* TriangMeshToVTK<DataType>::setPolyDataVertexScalars ( vtkPolyData* pd, const aol::Vector<DataType> &v, DataType &minC, DataType &maxC ) {
  int i;
  int N = pd->GetNumberOfPoints();

  vtkFloatArray* f = vtkFloatArray::New();

  minC = aol::NumberTrait<DataType>::Inf;
  maxC = -aol::NumberTrait<DataType>::Inf;

  f->SetNumberOfValues ( N );
  for ( i = 0; i < N; ++i ) {
    const DataType val = v[i] ;
    if ( val > maxC )
      maxC = val;
    if ( val < minC )
      minC = val;

    f->SetValue ( i, static_cast<float> ( val )  );
  }
  pd->GetPointData()->SetScalars ( f );

  f->Delete();

  return pd;
}

//! copies content of \a v into \a pd's scalar field
//!
//! \warning vtkPolyData has only one scalar field. When you call
//! this function, the previously stored scalar field will get lost.
template<class DataType>
vtkPolyData* TriangMeshToVTK<DataType>::setPolyDataTriangScalars ( vtkPolyData* pd, const aol::Vector<DataType> &v, DataType &minC, DataType &maxC ) {
  int i;
  int N = pd->GetNumberOfCells();

  vtkFloatArray* f = vtkFloatArray::New();

  minC = aol::NumberTrait<DataType>::Inf;
  maxC = -aol::NumberTrait<DataType>::Inf;

  f->SetNumberOfValues ( N );
  for ( i = 0; i < N; ++i ) {
    const DataType val = v[i] ;
    if ( val > maxC )
      maxC = val;
    if ( val < minC )
      minC = val;

    f->SetValue ( i, static_cast<float> ( val )  );
  }
  pd->Update();
  pd->GetCellData()->SetScalars ( f );
  pd->Update();

  f->Delete();

  return pd;
}

//! set texture coordinates
//!
//! This does not affect the scalar field which you can set via TriangMeshToVTK::setPolyDataScalars().
template<class DataType>
vtkPolyData* TriangMeshToVTK<DataType>::setPolyDataTCoords ( vtkPolyData* pd, const aol::Vector<DataType> & tCoordX, const aol::Vector<DataType> & tCoordY ) {

  vtkFloatArray * tCoords = vtkFloatArray::New();
  tCoords->SetNumberOfComponents ( 2 );
  tCoords->SetNumberOfTuples ( pd->GetNumberOfPoints() );
  for ( int i = 0; i < pd->GetNumberOfPoints(); i++ ) {
    float tCoord[2] = { static_cast<float> ( tCoordX[i] ), static_cast<float> ( tCoordY[i] ) };
    tCoords->SetTuple ( i, tCoord );
  }

  // set this array as texture coordinates
  pd->GetPointData()->SetTCoords ( tCoords );

  pd->BuildCells();
  pd->BuildLinks();
  pd->Update();

  return pd;
}

template class TriangMeshToVTK<float>;

template class TriangMeshToVTK<double>;

}
