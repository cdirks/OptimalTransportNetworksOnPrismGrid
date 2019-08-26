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

#include <tpCFEElement.h>
#include <tpCFEGrid.h>

namespace tpcfe {

template < typename RealType >
void CFEElement<RealType>::getTetraVertices ( CFETetra<RealType> &t ) const {
  aol::Vec3<RealType> vertex[4];
  for ( int l = 0; l < 4; l++ ) {
    const short a = t ( l, 0 );
    const short b = t ( l, 1 );
    if ( b == NODE_NOT_VIRTUAL ) { // Is not a virtual node
      vertex[l][0] = CFELookup<RealType>::_hexNbOffs[a][0];
      vertex[l][1] = CFELookup<RealType>::_hexNbOffs[a][1];
      vertex[l][2] = CFELookup<RealType>::_hexNbOffs[a][2];
    } else {       // Is a virtual node and thus an interpolation between two real nodes
      RealType w1 = _cutRelation[a][b];
      RealType w2 = _cutRelation[b][a];

      vertex[l][0] = w2 * CFELookup<RealType>::_hexNbOffs[a][0] + w1 * CFELookup<RealType>::_hexNbOffs[b][0];
      vertex[l][1] = w2 * CFELookup<RealType>::_hexNbOffs[a][1] + w1 * CFELookup<RealType>::_hexNbOffs[b][1];
      vertex[l][2] = w2 * CFELookup<RealType>::_hexNbOffs[a][2] + w1 * CFELookup<RealType>::_hexNbOffs[b][2];
    }
  }
  t._edge[0] = vertex[1] - vertex[0];
  t._edge[1] = vertex[2] - vertex[0];
  t._edge[2] = vertex[3] - vertex[0];
}

template < typename RealType >
aol::Vec3<RealType> CFEElement<RealType>::interfaceNormalOnRegTetra ( const int i ) const {
  aol::Vec3<RealType> normal;

  for ( int j = 0; j < 4; ++j ) {
    const int vertexIndex = CFELookup<RealType>::_stdTetraVertex[i][j];
    const RealType value = this->_structureValue[vertexIndex];
    const aol::Vec3<RealType> localNormal = CFELookup<RealType>::_stdTetraBasisGradients[i][j];
    normal += value * localNormal;
  }
  normal.normalize();

  return ( normal );
}


template class CFEElement< float >;
template class CFEElement< double >;
template class CFEElement< long double >;

}

