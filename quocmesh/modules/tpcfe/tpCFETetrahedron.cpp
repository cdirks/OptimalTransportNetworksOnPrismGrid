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

#include <tpCFETetrahedron.h>
#include <tpCFEElement.h>

namespace tpcfe {

const char CFETopoTetra::volumeTypeName[23][8] = {  "P", "PPPN", "PPPN1", "PPPN2", "PPPN3", "NNNP",
                                                    "NNNP1", "NNNP2", "NNNP3", "NNNPPP", "NNNPPP1",
                                                    "NNNPPP2", "NNNPPP3", "NNNPPP4", "NNNPPP5", "PPPNNN",
                                                    "PPPNNN1", "PPPNNN2", "PPPNNN3", "PPPNNN4", "PPPNNN5", "N", "NO_TYPE" };

template <typename RealType>
void CFETopoTetra::computeLocalCoordinate ( aol::Vec3<RealType> &tn0x, const CFEElement<RealType> &el, const int index ) const {
  aol::Vec3<RealType> tn1x;
  const short n0 = _node[index][0], n1 = _node[index][1];

  if ( n1 == NODE_NOT_VIRTUAL ) {
    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = CFELookup<RealType>::_hexNbOffs[n0][i];
    }
  } else {
    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = CFELookup<RealType>::_hexNbOffs[n0][i];
      tn1x[i] = CFELookup<RealType>::_hexNbOffs[n1][i];
    }

    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = el.getCutRelation ( n0, n1 ) * tn1x[i] + el.getCutRelation ( n1, n0 ) * tn0x[i];
    }
  }
}

template <typename RealType>
void CFETopoTetra::computeLocalCoordinate ( aol::Vec3<RealType> &tn0x, RealType _cutRelation[8][8], const int index ) const {
  aol::Vec3<RealType> tn1x;
  const short n0 = _node[index][0], n1 = _node[index][1];

  if ( n1 == NODE_NOT_VIRTUAL ) {
    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = CFELookup<RealType>::_hexNbOffs[n0][i];
    }
  } else {
    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = CFELookup<RealType>::_hexNbOffs[n0][i];
      tn1x[i] = CFELookup<RealType>::_hexNbOffs[n1][i];
    }

    for ( short i = 0; i < 3; ++i ) {
      tn0x[i] = _cutRelation[n0][n1] * tn1x[i] + _cutRelation[n1][n0] * tn0x[i];
    }
  }
}


template <typename RealType>
void CFETopoTetra::computeGlobalCoordinate ( aol::Vec3<RealType> &tn0x, const CFEElement<RealType> &el, const int index ) const {
  computeLocalCoordinate ( tn0x, el, index );
  for ( short i = 0; i < 3; ++i ) {
    tn0x[i] += static_cast<RealType> ( el[i] ) ;
  }
}

template void CFETopoTetra::computeLocalCoordinate<float> ( aol::Vec3<float> &tn0x, float _cutRelation[8][8], const int index ) const;
template void CFETopoTetra::computeLocalCoordinate<double> ( aol::Vec3<double> &tn0x, double _cutRelation[8][8], const int index ) const;
template void CFETopoTetra::computeLocalCoordinate<long double> ( aol::Vec3<long double> &tn0x, long double _cutRelation[8][8], const int index ) const;

template void CFETopoTetra::computeLocalCoordinate<float> ( aol::Vec3<float> &tn0x, const CFEElement<float> &el, const int index ) const;
template void CFETopoTetra::computeLocalCoordinate<double> ( aol::Vec3<double> &tn0x, const CFEElement<double> &el, const int index ) const;
template void CFETopoTetra::computeLocalCoordinate<long double> ( aol::Vec3<long double> &tn0x, const CFEElement<long double> &el, const int index ) const;

template void CFETopoTetra::computeGlobalCoordinate<float> ( aol::Vec3<float> &tn0x, const CFEElement<float> &el, const int index ) const;
template void CFETopoTetra::computeGlobalCoordinate<double> ( aol::Vec3<double> &tn0x, const CFEElement<double> &el, const int index ) const;
template void CFETopoTetra::computeGlobalCoordinate<long double> ( aol::Vec3<long double> &tn0x, const CFEElement<long double> &el, const int index ) const ;




template <typename RealType>
RealType CFETetra<RealType>::determinant() const {
  const aol::Vec3<RealType> &a = _edge[0];
  const aol::Vec3<RealType> &b = _edge[1];
  const aol::Vec3<RealType> &c = _edge[2];

  RealType det = 0.0;

  det += a.x() * ( b.y() * c.z() - c.y() * b.z() );
  det -= a.y() * ( b.x() * c.z() - c.x() * b.z() );
  det += a.z() * ( b.x() * c.y() - b.y() * c.x() );

  return ( det );
}


template <typename RealType>
void CFETetra<RealType>::computeInverseTransformation() const {
  const aol::Vec3<RealType> &a = _edge[0];
  const aol::Vec3<RealType> &b = _edge[1];
  const aol::Vec3<RealType> &c = _edge[2];

  _barycentricCoordinates.setRow ( 0, b.y() * c.z() - b.z() * c.y(),
                                   c.x() * b.z() - b.x() * c.z(),
                                   b.x() * c.y() - c.x() * b.y() );

  _barycentricCoordinates.setRow ( 1, c.y() * a.z() - a.y() * c.z(),
                                   a.x() * c.z() - c.x() * a.z(),
                                   c.x() * a.y() - a.x() * c.y() );

  _barycentricCoordinates.setRow ( 2, a.y() * b.z() - b.y() * a.z(),
                                   b.x() * a.z() - a.x() * b.z(),
                                   a.x() * b.y() - b.x() * a.y() );
}


template <typename RealType>
void CFETetra<RealType>::computeStiffnessMatrix() const {
  RealType value = 0.0, all = 0.0;

  for ( short i = 0; i < 3; ++i ) {
    for ( short j = 0; j < 3; ++j ) {
      value = 0.0;
      for ( short k = 0; k < 3; ++k ) {
        value += _barycentricCoordinates[i][k] * _barycentricCoordinates[j][k];
      }
      _lsm.set ( i + 1, j + 1, value );
      _lsm.set ( j + 1, i + 1, value );
    }
  }
  for ( short i = 1; i < 4; ++i ) {
    value = 0.0;
    for ( short j = 1; j < 4; ++j ) {
      value -= _lsm.get ( i, j );
    }
    _lsm.set ( i, 0, value );
    _lsm.set ( 0, i, value );
    all += value;
  }
  _lsm.set ( 0, 0, -all );
}


template <typename RealType>
void CFETetra<RealType>::dump( ostream& out /* = cout */ ) const {
  out << "Dumping CFETetra: " << endl
      << "Edges: " << _edge[0] << " " << _edge[1] << " " << _edge[2] << endl
      << "local stiffness matrix: " << _lsm << endl
      << "barycentric coordinates: " << _barycentricCoordinates << endl
      << "Volume type " << _volType << " = " << volumeTypeName[_volType ] << endl
      << "Volume = " << _volume << endl
      << "Number of virtual nodes = " << _virtualNodeNum << endl
      << "sign = " << static_cast<short> ( _sign ) << endl
      << "parent = " << _parent << endl
      << "nodes: " << _node[0][0] << " " << _node[0][1] << ",  " << _node[1][0] << " " << _node[1][1] << ",  " << _node[2][0] << " " << _node[2][1] << ",  " << _node[3][0] << " " << _node[3][1] << endl << endl;
}


template class CFETetra < float >;
template class CFETetra < double >;
template class CFETetra < long double >;

}
