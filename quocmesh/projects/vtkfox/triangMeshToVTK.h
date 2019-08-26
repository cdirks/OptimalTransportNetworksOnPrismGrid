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

#ifndef __TRIANGMESHTOVTK_H
#define __TRIANGMESHTOVTK_H

#include <vtkIncludes.h>

#include <aol.h>
#include <triangMesh.h>

using namespace std;

namespace aol {

// forward declaration
template< typename DataType, typename TriangleType > class TriangMesh;

//! a class to convert TriangMeshes for use in VTK
template<class DataType>
class TriangMeshToVTK {
public:
  TriangMeshToVTK() {};

  ~TriangMeshToVTK() {};

  static vtkPolyData* makePolyData ( const TriangMesh<DataType> &mesh, DataType &minC, DataType &maxC, bool copyVertexData = true );

  static vtkPolyData* setPolyDataVertexScalars ( vtkPolyData* pd, const aol::Vector<DataType> &v, DataType &minC, DataType &maxC );
  static vtkPolyData* setPolyDataTriangScalars ( vtkPolyData* pd, const aol::Vector<DataType> &v, DataType &minC, DataType &maxC );

  static vtkPolyData* setPolyDataTCoords ( vtkPolyData* pd, const aol::Vector<DataType> & tCoordX, const aol::Vector<DataType> & tCoordY );

};

}



#endif

