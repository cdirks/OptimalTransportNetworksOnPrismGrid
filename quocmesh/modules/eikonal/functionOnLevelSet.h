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

#ifndef __FUNCTIONONLEVELSET_H
#define __FUNCTIONONLEVELSET_H

#include <quoc.h>
#include <isolineIterator2d.h>
#include <signedDistanceOp.h>
#include <geom.h>
#include <hashSet.h>

namespace qc {

template <typename RealType>
class EdgeNode {
public:
  EdgeNode ( const qc::CoordType &node1, const qc::CoordType &node2, RealType locCoord )
      : _node1 ( node1 ), _node2 ( node2 ), _locCoord ( locCoord ), _v ( 0. ), _w ( 0. ) {}

  qc::CoordType _node1, _node2; // local Nodes of edges
  RealType _locCoord;     // where the node is on [x[loc1],x[loc2]]

  RealType _v;
  RealType _w;

  bool operator== ( const EdgeNode<RealType> &el ) const {
    return ( ( _node1 == el._node1 ) && ( _node2 == el._node2 ) );
  }

};

#if 0
// this does not work. can possibly repaired like in aol/hashSet
//! GNU C++ Extensions
namespace __gnu_cxx {

template <typename RealType> struct hash< qc::EdgeNode<RealType> > {
  size_t operator() ( const qc::EdgeNode<RealType>& t ) const {
    return t._node1.x() + ( t._node1.z() * 769 + t._node1.y() ) * 1543;
  }
};

}
#endif


template <typename RealType, typename Imp>
class ComputeFunctionOnLevelSetInterface {
public:
#ifdef  _MSC_VER
  typedef stdext::hash_set<EdgeNode<RealType>, stdext::_Hash<EdgeNode<RealType> > > hashType;
#elif defined ( _LIBCPP_VERSION )
  typedef __1::unordered_set< EdgeNode<RealType> > hashType;
#else
  typedef tr1::unordered_set< EdgeNode<RealType> > hashType;
#endif

  typedef typename hashType::iterator hashIteratorType;

protected:
  hashType _hash;
  const GridDefinition &_grid;
  const ScalarArray<RealType, qc::QC_2D> &_ls;
public:
  ComputeFunctionOnLevelSetInterface ( const GridDefinition &grid, ScalarArray<RealType, qc::QC_2D> &ls )
      : _grid ( grid ), _ls ( ls ) {}

  // this is the only interface function, which has to be defined in derived classes
  RealType evaluateFunction ( const Element &el, const aol::Vec2<RealType> locCoords ) const {
    return asImp().evaluateFunction ( el, locCoords );
  }

  //! call this function to initialize this class. It traverses all edges of the level-set
  //! and integrates the function along these edges
  void computeFunction( ) {
    IsoLineManager2d<RealType> isoManager ( _grid, _ls );
    const RealType alpha[2] = { 0.21132487, 0.78867513 };

    for ( IsoLineIterator2d<RealType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {

      RealType len;
      // RealType w[2] = { isoIt->weights[0], isoIt->weights[1] };
      aol::Vec3<RealType> pt ( isoIt->el.x(), isoIt->el.y(), 0. );
      aol::Vec3<RealType> qpts[2];

      // calculate length of segment
      aol::Vec3<RealType> tmp;
      tmp =   isoIt->points[0];
      tmp -=  isoIt->points[1];
      len = tmp.norm();

      for ( int j = 0; j < 2; j++ ) {
        // compute quadrature points on line
        qpts[j]  = isoIt->points[0];
        qpts[j] *= alpha[j] / ( 1. - alpha[j] );
        qpts[j] += isoIt->points[1];
        qpts[j] *= 1. - alpha[j];

        qpts[j] -= pt; // make reference coords
      }

      // quadrature
      RealType a = 0.;
      for ( int j = 0; j < 2; j++ ) {
        aol::Vec2<RealType> locCoords ( qpts[j][0], qpts[j][1] );
        a += evaluateFunction ( isoIt->el, locCoords ) * ( 1. - alpha[j] );
      }
      // cerr << a << " ";
      a *= len * _grid.H();

      EdgeNode<RealType> en1 ( isoIt->intersect[0][0], isoIt->intersect[0][1], isoIt->weights[0] );
      EdgeNode<RealType> en2 ( isoIt->intersect[1][0], isoIt->intersect[1][1], isoIt->weights[1] );

      if ( _hash.find ( en1 ) != _hash.end() ) {
        hashIteratorType it = _hash.find ( en1 );
        en1 = *it;
        en1._v += a;
        en1._w += len * _grid.H();
        _hash.erase ( it );
        _hash.insert ( en1 );
      } else {
        en1._v += a;
        en1._w += len * _grid.H();
        _hash.insert ( en1 );
      }

      if ( _hash.find ( en2 ) != _hash.end() ) {
        hashIteratorType it = _hash.find ( en1 );
        en2 = *it;
        en2._v += a;
        en2._w += len * _grid.H();
        _hash.erase ( it );
        _hash.insert ( en2 );
      } else {
        en2._v += a;
        en2._w += len * _grid.H();
        _hash.insert ( en2 );
      }
    }
  }

  //! call this function to output the extension of the function on the level-set
  //! to a ScalarArray<QC_2D>
  void extend ( ScalarArray<RealType, qc::QC_2D> &extension ) {
    IsoLineManager2d<RealType> isoManager ( _grid, _ls );
    ScalarArray<RealType, qc::QC_2D> dist ( _grid );
    dist.setAll ( 10. ); // just must be larger than diam of reference element
    for ( IsoLineIterator2d<RealType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {

      EdgeNode<RealType> en1 ( isoIt->intersect[0][0], isoIt->intersect[0][1], isoIt->weights[0] );
      EdgeNode<RealType> en2 ( isoIt->intersect[1][0], isoIt->intersect[1][1], isoIt->weights[1] );
      RealType vel[2];
      if ( _hash.find ( en1 ) != _hash.end() ) {
        hashIteratorType it = _hash.find ( en1 );
        vel[0] = it->_v / it->_w;
      } else {
        throw aol::Exception ( "edgenode not found", __FILE__, __LINE__ );
      }
      if ( _hash.find ( en2 ) != _hash.end() ) {
        hashIteratorType it = _hash.find ( en2 );
        vel[1] = it->_v / it->_w;
      } else {
        throw aol::Exception ( "edgenode not found", __FILE__, __LINE__ );
      }

      for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 2; j++ ) {
          short x = isoIt->intersect[i][j][0];
          short y = isoIt->intersect[i][j][1];
          aol::Vec3<RealType> pt ( x, y, 0. );
          aol::Vec3<RealType> a, b;
          RealType d, v;
          a = isoIt->points[0];
          a -= pt;
          b = isoIt->points[1];
          b -= isoIt->points[0];
          RealType lambda = - ( a * b ) / ( b * b );
          if ( lambda >= 0. && lambda <= 1. ) {
            b *= lambda;
            a += b;
            d =  a.norm();
            v = ( 1.0 - lambda ) * vel[0] + lambda * vel[1];
          } else {
            RealType d1, d2;
            a = isoIt->points[1];
            a -= pt;
            d1 = a.norm();
            a = isoIt->points[1];
            a -= pt;
            d2 = a.norm();
            if ( d1 <= d2 ) {
              v = vel[0];
              d = d1;
            } else {
              v = vel[1];
              d = d2;
            }
          }

          if ( dist.get ( x, y ) > d ) {
            dist.set ( x, y, d );
            extension.set ( x, y, v );
          }
        }
      }
    }
    SignedDistAndExtVelocityOp<RealType> signedDistOp ( _grid );
    signedDistOp.setExtVelocityField (  extension );
    signedDistOp.apply ( _ls, dist ); // throw away dist anyways
  }

protected:
// barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

template <typename RealType, typename DiscFuncType>
class ExtendFromLevelSet {
protected:
  const GridDefinition &_grid;
  const ScalarArray<RealType, qc::QC_2D> &_ls;
public:
  ExtendFromLevelSet ( const GridDefinition &grid, const ScalarArray<RealType, qc::QC_2D> &ls )
      : _grid ( grid ), _ls ( ls ) {
  }


  void extend ( const DiscFuncType &interfaceFunction, ScalarArray<RealType, qc::QC_2D> &extension ) {
    IsoLineManager2d<RealType> isoManager ( _grid, _ls );
    ScalarArray<RealType, qc::QC_2D> dist ( _grid );
    dist.setAll ( 1000. ); // just must be larger than diam of reference element

    for ( IsoLineIterator2d<RealType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {

      const qc::Element &el = isoIt->el;

      aol::Vec2<RealType> endPoints[2] = {
                                           aol::Vec2<RealType> ( isoIt->points[0][0], isoIt->points[0][1] ),
                                           aol::Vec2<RealType> ( isoIt->points[1][0], isoIt->points[1][1] )
                                         };

      aol::LineSegment<RealType, qc::QC_2D> seg ( endPoints[0], endPoints[1] );

      for ( qc::ElementNodeIterator<qc::QC_2D> enit = qc::GridDefinition::elementNodeBegin<qc::QC_2D> ( el );
            enit != qc::GridDefinition::elementNodeEnd<qc::QC_2D> ( el ); ++enit ) {


        aol::Vec2<RealType> pt ( enit->x(), enit->y() ), projPt;

        seg.projectTo ( pt, projPt );
        RealType d = euclidianDist ( pt, projPt );

        if ( dist.get ( *enit ) > d ) {
          dist.set ( *enit, d );

          aol::Vec2<RealType> localCoords ( projPt[0] - el.x(), projPt[1] - el.y() );

          extension.set ( *enit, interfaceFunction.evaluate ( el, localCoords ) );
        }
      }
    }
    SignedDistAndExtVelocityOp<RealType> signedDistOp ( _grid );
    signedDistOp.setExtVelocityField (  extension );
    signedDistOp.apply ( _ls, dist ); // throw away dist anyways
  }
};

}

#endif
