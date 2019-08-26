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

#ifndef __SIGNEDDISTANCEOP_H
#define __SIGNEDDISTANCEOP_H

#include <gridBase.h>
#include <op.h>
#include <eikonalNA.h>
#include <eikonalNA3d.h>
#include <isolineIterator2d.h>
#include <qmException.h>
#include <triangle.h>
#include <configurators.h>

#include <tpCFEGrid.h>
#include <tpCFEUtils.h>

namespace qc {

/** Base class of classes computing the signed distance function.
 * Provides initialization along level set and works for general
 * types of eikonal solvers.
 * \author Droske
 */
template <typename ConfiguratorType, class EikonalType>
class SignedDistanceBase {
protected:
  typedef typename ConfiguratorType::RealType RealType;
public:
  void setIsoValue ( RealType Value ) {
    _isoValue = Value;
  }
protected:
  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  mutable EikonalType _eik;
  RealType _isoValue;

  SignedDistanceBase ( const typename ConfiguratorType::InitType &Grid )
      : _config ( Grid ), _grid ( Grid ), _eik ( Grid ), _isoValue ( 0. ) {}

  void initializeBoundary ( const Array<RealType> &Func ) const {
    typename ConfiguratorType::ElementIteratorType it;
    for ( it = _config.begin(); it != _config.end(); ++it ) {
      int numpos = 0;
      short x = it->x(), y = it->y();
      if ( Func.get ( x, y ) > _isoValue ) numpos++;
      if ( Func.get ( x + 1, y ) > _isoValue ) numpos++;
      if ( Func.get ( x + 1, y + 1 ) > _isoValue ) numpos++;
      if ( Func.get ( x, y + 1 ) > _isoValue ) numpos++;

      if ( numpos == 1 || numpos == 3 ) {
        qc::CoordType corner;
        qc::CoordType nb1, nb2;
        if ( numpos == 1 ) {
          if ( Func.get ( x, y ) > _isoValue ) {
            corner.set ( x, y, 0 );
            nb1.set ( x + 1, y, 0 );
            nb2.set ( x, y + 1, 0 );
          } else if ( Func.get ( x + 1, y ) > _isoValue ) {
            corner.set ( x + 1, y, 0 );
            nb1.set ( x, y, 0 );
            nb2.set ( x + 1, y + 1, 0 );
          } else if ( Func.get ( x + 1, y + 1 ) > _isoValue ) {
            corner.set ( x + 1, y + 1, 0 );
            nb1.set ( x, y + 1, 0 );
            nb2.set ( x + 1, y, 0 );
          } else {
            corner.set ( x, y + 1, 0 );
            nb1.set ( x, y, 0 );
            nb2.set ( x + 1, y + 1, 0 );
          }
        } else {
          if ( Func.get ( x, y ) <= _isoValue ) {
            corner.set ( x, y, 0 );
            nb1.set ( x + 1, y, 0 );
            nb2.set ( x, y + 1, 0 );
          } else if ( Func.get ( x + 1, y ) <= _isoValue ) {
            corner.set ( x + 1, y, 0 );
            nb1.set ( x, y, 0 );
            nb2.set ( x + 1, y + 1, 0 );
          } else if ( Func.get ( x + 1, y + 1 ) <= _isoValue ) {
            corner.set ( x + 1, y + 1, 0 );
            nb1.set ( x, y + 1, 0 );
            nb2.set ( x + 1, y, 0 );
          } else {
            corner.set ( x, y + 1, 0 );
            nb1.set ( x, y, 0 );
            nb2.set ( x + 1, y + 1, 0 );
          }
        }
        initDistToCorner ( Func, corner, nb1, nb2 );
      } else if ( numpos == 2 ) {
        RealType v11 = Func.get ( x, y );
        RealType v21 = Func.get ( x + 1, y );
        RealType v22 = Func.get ( x + 1, y + 1 );
        //RealType v12 = Func.get( x, y+1 );
        if ( ( v11 > _isoValue ) == ( v22 > _isoValue ) ) {
          // the ambiguous case: values of opposite corners have same sign
          cerr << "\n\n\n ambiguous case!!!! \n\n\n";
          qc::CoordType coord ( x, y, 0 );
          initAmbiguousCase ( Func, coord );
        } else {
          qc::CoordType c1, c2, c3, c4;
          if ( ( v11 > _isoValue ) == ( v21 > _isoValue ) ) {
            c1.set ( x, y, 0 );
            c2.set ( x, y + 1, 0 );
            c3.set ( x + 1, y, 0 );
            c4.set ( x + 1, y + 1, 0 );
          } else {
            c1.set ( x, y, 0 );
            c2.set ( x + 1, y, 0 );
            c3.set ( x, y + 1, 0 );
            c4.set ( x + 1, y + 1, 0 );
          }
          initDistIntersectTwoEdges ( Func, c1, c2, c3, c4 );
        }
      }
    }
  }

  void initDistToCorner ( const Array<RealType> &Func,
                          const qc::CoordType &Corner,
                          const qc::CoordType &Nb1,
                          const qc::CoordType &Nb2 ) const {
    RealType vc = Func.get ( Corner );
    RealType vnb1 = Func.get ( Nb1 );
    RealType vnb2 = Func.get ( Nb2 );

    RealType d1 = ( _isoValue - vc ) / ( vnb1 - vc ) * _grid.H();
    RealType d2 = ( _isoValue - vc ) / ( vnb2 - vc ) * _grid.H();

    if ( vnb1 == vc ) {
      d1 = 0.5 * _grid.H();
    }
    if ( vnb2 == vc ) {
      d2 = 0.5 * _grid.H();
    }
    RealType d = sqrt ( d1 * d1 * d2 * d2 / ( d1 * d1 + d2 * d2 ) );
    if ( d1 == 0. || d2 == 0. ) {
      d = 0.;
    }

    _eik.setSeedPoint ( Corner, d );
    _eik.setSeedPoint ( Nb1, _grid.H() - d1 );
    _eik.setSeedPoint ( Nb2, _grid.H() - d2 );
  }

  void initDistIntersectTwoEdges ( const Array<RealType> &Func,
                                   const qc::CoordType &C1, const qc::CoordType &C2,
                                   const qc::CoordType &C3, const qc::CoordType &C4 ) const {
    const RealType v1 = Func.get ( C1 );
    const RealType v2 = Func.get ( C2 );
    const RealType v3 = Func.get ( C3 );
    const RealType v4 = Func.get ( C4 );

    const RealType w1 = ( _isoValue - v1 ) / ( v2 - v1 );
    const RealType w2 = ( _isoValue - v3 ) / ( v4 - v3 );

    aol::Vec3<RealType> pt1 ( static_cast< RealType >( C2.x() * w1 ) + static_cast< RealType >( C1.x() * ( 1. - w1 ) ), static_cast< RealType >( C2.y() * w1 ) + static_cast< RealType >( C1.y() * ( 1. - w1 ) ), 0. );
    aol::Vec3<RealType> pt2 ( static_cast< RealType >( C4.x() * w2 ) + static_cast< RealType >( C3.x() * ( 1. - w2 ) ), static_cast< RealType >( C4.y() * w2 ) + static_cast< RealType >( C3.y() * ( 1. - w2 ) ), 0. );


    aol::LineSegment<RealType,qc::QC_3D> segment ( pt1, pt2 );

    aol::Vec3<RealType> pt;

    pt.set ( C1.x(), C1.y(), 0. );
    _eik.setSeedPoint ( C1, segment.dist ( pt ) * _grid.H() );
    pt.set ( C2.x(), C2.y(), 0. );
    _eik.setSeedPoint ( C2, segment.dist ( pt ) * _grid.H() );
    pt.set ( C3.x(), C3.y(), 0. );
    _eik.setSeedPoint ( C3, segment.dist ( pt ) * _grid.H() );
    pt.set ( C4.x(), C4.y(), 0. );
    _eik.setSeedPoint ( C4, segment.dist ( pt ) * _grid.H() );
  }

  void initAmbiguousCase ( const Array<RealType> &Func,
                           const qc::CoordType &Coord ) const {
    const RealType v11 = Func.get ( Coord.x(), Coord.y() );
    const RealType v21 = Func.get ( Coord.x() + 1, Coord.y() );
    const RealType v22 = Func.get ( Coord.x() + 1, Coord.y() + 1 );
    const RealType v12 = Func.get ( Coord.x(), Coord.y() + 1 );

    const short x = Coord.x();
    const short y = Coord.y();


    RealType len1, len2;
    const RealType y1 = ( _isoValue - v11 ) / ( v12 - v11 );
    const RealType y2 = ( _isoValue - v21 ) / ( v22 - v21 );
    const RealType x1 = ( _isoValue - v11 ) / ( v21 - v11 );
    const RealType x2 = ( _isoValue - v12 ) / ( v22 - v12 );

    // case 1: [\\]
    aol::Vec3<RealType> pt1 ( 0., y1, 0. );
    aol::Vec3<RealType> pt2 ( x1, 0., 0. );
    aol::LineSegment<RealType,qc::QC_3D> segment ( pt1, pt2 );
    len1 = segment.length();

    pt1.set ( x2, 1., 0. );
    pt2.set ( 1., y2, 0. );
    len1 += segment.length();

    // case 2: [//]
    pt1.set ( x1, 0., 0. );
    pt2.set ( 1., y2, 0. );
    len2 = segment.length();

    pt1.set ( 0., y1, 0. );
    pt2.set ( x2, 1., 0. );
    len2 += segment.length();

    if ( len1 < len2 ) {
      qc::CoordType c;
      aol::Vec3<RealType> pt;
      pt1.set ( 0., y1, 0. );
      pt2.set ( x1, 0., 0. );

      c.set ( x, y, 0 );
      pt.set ( 0., 0., 0. );
      _eik.setSeedPoint ( c, segment.dist ( pt ) * _grid.H() );

      c.set ( x + 1, y, 0 );
      _eik.setSeedPoint ( c, aol::Min ( static_cast<RealType>(1.0) - x1, y2 ) * _grid.H() );

      pt1.set ( x2, 1., 0. );
      pt2.set ( 1., y2, 0. );
      c.set ( x + 1, y + 1, 0 );
      pt.set ( 1., 1., 0. );
      _eik.setSeedPoint ( c, segment.dist ( pt ) * _grid.H() );

      c.set ( x, y + 1, 0 );
      _eik.setSeedPoint ( c, aol::Min ( static_cast<RealType>(1.0) - y1, x2 ) * _grid.H() );
    } else {
      // [//]
      qc::CoordType c;
      aol::Vec3<RealType> pt;
      pt1.set ( x1, 0., 0. );
      pt2.set ( 1., y2, 0. );

      c.set ( x + 1, y, 0 );
      pt.set ( 1., 0., 0. );
      _eik.setSeedPoint ( c, segment.dist ( pt ) * _grid.H() );

      c.set ( x, y, 0 );
      _eik.setSeedPoint ( c, aol::Min ( x1, y1 ) * _grid.H() );

      c.set ( x + 1, y + 1, 0 );
      _eik.setSeedPoint ( c, aol::Min ( 1. - y2, 1. - x2 ) * _grid.H() );

      pt1.set ( 0., y1, 0. );
      pt2.set ( x2, 1., 0. );
      c.set ( x, y + 1, 0 );
      pt.set ( 0., 1., 0. );
      _eik.setSeedPoint ( c, segment.dist ( pt ) * _grid.H() );
    }
  }


};

/**
 * compute the signed distance function to the level lines of the function given by the argument
 * with specified isovalue and write the result in the destination vector.
 * 2d-only.
 * \author Droske
 */
template <typename ConfiguratorType>
class SignedDistanceOp : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > ,
      public SignedDistanceBase<ConfiguratorType, eik::EikonalNA<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType> > {
  typedef typename ConfiguratorType::RealType RealType;
public:
  SignedDistanceOp ( const typename ConfiguratorType::InitType &Grid ) :
      SignedDistanceBase<ConfiguratorType, eik::EikonalNA<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType> > ( Grid ) {
    if ( Grid.getDimOfWorld() != 2 ) {
      throw aol::Exception ( "\nERROR: qc::SignedDistanceOp works only in 2D! \n"
                             "For 3D use for example qc::SignedDistanceOp3D instead. \n",
                             __FILE__,  __LINE__ );
    }
  }

  virtual ~SignedDistanceOp() {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg_, aol::Vector<RealType> &Dest_ ) const {
    const Array<RealType> Arg ( Arg_, this->_grid );
    Array<RealType> Dest ( Dest_, this->_grid );

    this->_eik.reset();
    this->initializeBoundary ( Arg );
    this->_eik.march();

    typename ConfiguratorType::InitType::OldAllNodeIterator it;
    for ( it = this->_grid._nBeginIt; it != this->_grid._nEndIt; ++it ) {
      RealType v = this->_eik.getTimeField().get ( *it );
      if ( Arg.get ( *it ) <= this->_isoValue ) {
        Dest.set( *it, -v );
      } else {
        Dest.set( *it, v );
      }
    }
  }
};


/** This version of signed distance operator additionaly computes
 * extension velocities, which are given in the neighborhood of the
 * level set by solving \f$\nabla F^{ext}\cdot \phi^{ext} = 0\f$.
 * \author Droske
 */
template <typename ConfiguratorType>
class SignedDistAndExtVelocityOp : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > ,
      public SignedDistanceBase<ConfiguratorType, eik::EikonalNAWithExtensionVelocity<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
public:
  SignedDistAndExtVelocityOp ( const typename ConfiguratorType::InitType &Grid ) :
      SignedDistanceBase<ConfiguratorType, eik::EikonalNAWithExtensionVelocity<RealType> > ( Grid ) {}

  virtual ~SignedDistAndExtVelocityOp() {}

  void setExtVelocityField ( ScalarArray<RealType, qc::QC_2D> &extVelocity ) {
    this->_eik.setExtVelocityField ( extVelocity );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg_, aol::Vector<RealType> &Dest_ ) const {
    const Array<RealType> Arg ( Arg_, this->_grid );
    Array<RealType> Dest ( Dest_, this-> _grid );

    this->_eik.reset();
    this->_eik.checkTrialNodes();
    this->initializeBoundary ( Arg );
    this->_eik.march();

    GridDefinition::OldFullNodeIterator it;
    for ( it = this->_grid.begin(); it != this->_grid.end(); ++it ) {
      RealType v = this->_eik.getTimeField().get( *it );
      if ( Arg.get( *it ) <= this->_isoValue ) {
        Dest.set ( *it, -v );
      } else {
        Dest.set ( *it, v );
      }
    }
  }

};




template <typename ConfiguratorType, class EikonalType>
class SignedDistanceBase3D {
  typedef typename ConfiguratorType::RealType RealType;
public:
  void setIsoValue ( RealType Value ) { _isoValue = Value; }
protected:
  const typename ConfiguratorType::InitType &_grid;
  mutable EikonalType _eik;
  RealType _isoValue;

  RealType dist ( const aol::Vec3<RealType> &v1, const aol::Vec3<RealType> &v2 ) const {
    return sqrt ( aol::Sqr ( v1[0] - v2[0] ) + aol::Sqr ( v1[1] - v2[1] ) + aol::Sqr ( v1[2] - v2[2] ) );
  }

  SignedDistanceBase3D ( const typename ConfiguratorType::InitType &Grid )
      : _grid ( Grid ), _eik ( Grid ), _isoValue ( 0. ) { }

  void initializeBoundary ( const qc::ScalarArray<RealType, qc::QC_3D> &Func ) const {

    tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> cfeGrid ( qc::GridSize<qc::QC_3D>::createFrom ( Func ) );
    cfeGrid.setAdjustLevelset ( 1e-10 );
    cfeGrid.setDomainFrom ( Func );
    const RealType _h = _grid.H();

    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit( cfeGrid ); itit.notAtEnd(); ++itit ) {
      std::vector< aol::Vec3< RealType > > interfaceVertex = itit.getTriangRef().getGlobalCoordinates();
      for ( int i = 0; i < 3; ++i ) {
        interfaceVertex[i] *= _h;
      }
      const aol::Triangle<RealType> triangle ( interfaceVertex[0], interfaceVertex[1], interfaceVertex[2] );
      // interfaceVertex contains world coordinates of the vertices of the current interface triangle

      for ( int l = 0; l < 2 ; ++l ) { // attention: strange alphabetical ordering of loop indices
        for ( int j = 0 ; j < 2 ; ++j ) {
          for ( int k = 0 ; k < 2 ; ++k ) {
            const qc::CoordType elCoord ( itit.getTriangRef().getElement().x(), itit.getTriangRef().getElement().y(), itit.getTriangRef().getElement().z() );
            const qc::CoordType qcp = elCoord + qc::CoordType ( k, j, l );                          // one cube vertex in global coordinates
            const aol::Vec3<RealType> coord ( qcp[0] * _h, qcp[1] * _h, qcp[2] * _h );  // the same cube vertex in world coordinates

            RealType tmpdist = triangle.calcDist ( coord );
            if ( aol::isNaN ( tmpdist ) )
              throw aol::Exception( "tmpdist is NaN! Before 08-03-13 here a code was implemented which catched this case. If you are interested in what it did please use the history of this file ", __FILE__, __LINE__);

            _eik.setSeedPoint ( qcp, aol::Abs ( tmpdist ) + 1.0 ); // ON and OS have no idea why 1.0 is added here. But it is subtracted
                                                                   // again in the signedDistanceOp3d
          } // k
        } // j
      } // l

    }
  }
};


template <typename ConfiguratorType>
class SignedDistanceOp3DGeneral : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > ,
  public SignedDistanceBase3D<ConfiguratorType, eik::EikonalNA3D<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType> > {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
public:
  SignedDistanceOp3DGeneral ( const GridType &Grid ) :
      SignedDistanceBase3D<ConfiguratorType, eik::EikonalNA3D<RealType, GridType> > ( Grid ) { }

  virtual ~SignedDistanceOp3DGeneral( ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg_, aol::Vector<RealType> &Dest_ ) const {
#ifdef DEBUG
    if ( ( Arg_.getMinValue() >= 0 ) || ( Arg_.getMaxValue() <= 0 ) )
      throw aol::Exception ( "qc::SignedDistanceOp3DGeneral::apply: The argument doesn't have a zero level set!\n", __FILE__,  __LINE__ );
#endif
    const qc::Array<RealType> Arg ( Arg_, this->_grid );
    qc::Array<RealType> Dest ( Dest_, this->_grid );

    this->_eik.reset();
    const qc::ScalarArray<RealType, qc::QC_3D> arg3d ( Arg, aol::FLAT_COPY ); // why not interpret Arg_ as qc::ScalarArray<QC_3D> immediately ??
    this->initializeBoundary ( arg3d );
    this->_eik.march();

    typename ConfiguratorType::InitType::OldAllNodeIterator it;
    for ( it = this->_grid._nBeginIt; it != this->_grid._nEndIt; ++it ) {
      RealType v = this->_eik.getTimeField().get ( *it ) - 1.;    // here the mysterious 1 added in the call of setSeedPoint is subtracted again
      if ( Arg.get ( *it ) <= this->_isoValue ) {
        Dest.set( *it, -v );
      } else {
        Dest.set( *it, v );
      }
    }
  }
};

template <typename RealType>
class SignedDistanceOp3D : public SignedDistanceOp3DGeneral <qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> > > {
public:
  SignedDistanceOp3D ( const qc::GridDefinition &Grid ) :
    SignedDistanceOp3DGeneral <qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> > > ( Grid ) { }
};

template<typename RealType>
class SignedDistanceOp3DWithExtVel : public aol::Op<aol::Vector<RealType> >,
                                       public SignedDistanceBase3D<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> >, eik::EikonalNA3DWithExtVelocity<RealType> >  {
public:
  SignedDistanceOp3DWithExtVel( const qc::GridDefinition &Grid ):SignedDistanceBase3D<qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D, 3> >, eik::EikonalNA3DWithExtVelocity<RealType> > (Grid){
  }

  virtual ~SignedDistanceOp3DWithExtVel(){
  }

  void setExtVelocityField( qc::ScalarArray<RealType, qc::QC_3D> &extVelocity ){
    this->_eik.setExtVelocityField( extVelocity );
  }

 void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp( Dest.size() );
    apply( Arg, tmp );
    Dest += tmp;
  }

  void apply( const aol::Vector<RealType> &Arg_, aol::Vector<RealType> &Dest_ ) const {
    const qc::Array<RealType> Arg( Arg_, this->_grid );
    qc::Array<RealType> Dest( Dest_, this->_grid );

    qc::ScalarArray<RealType, qc::QC_3D> arg3d( Arg, this->_grid );

    this->_eik.reset();cerr<<" l1 ";
    this->_eik.checkTrialNodes();cerr<<" l1 ";
    initializeBoundary( arg3d );cerr<<" l1 ";
    this->_eik.march();cerr<<" l1 ";
    cerr<<" l1 ";
    //    qc::GridDefinition::OldFullNodeIterator it;
    //    for ( it = _grid.begin(); it != _grid.end(); ++it ) {
    //      RealType v = _eik.getTimeField()( *it ) - 1.;
    //      if ( Arg( *it ) <= _isoValue ) {
    //        Dest( *it ) = -v;
    //      } else {
    //      Dest( *it ) = v;
    //      }
    //    }
  }
};


/**
 * Trait to select the proper signed distance op based on the Configurator template.
 *
 *  \author Berkels
 */
template <typename ConfiguratorType, qc::Dimension Dim = ConfiguratorType::Dim>
class SignedDistanceOpTrait {};

template <typename ConfiguratorType>
class SignedDistanceOpTrait<ConfiguratorType, qc::QC_2D> {
public:
  typedef qc::SignedDistanceOp<ConfiguratorType> OpType;
};

template <typename ConfiguratorType>
class SignedDistanceOpTrait<ConfiguratorType, qc::QC_3D> {
public:
  typedef qc::SignedDistanceOp3DGeneral<ConfiguratorType> OpType;
};

} // end namespace qc

#endif
