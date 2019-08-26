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

#ifndef __NARROWBANDBOUNDARYINTEGRATION_H
#define __NARROWBANDBOUNDARYINTEGRATION_H

#include <aol.h>
#include <quoc.h>
#include <op.h>
#include <gridBase.h>
#include <discreteFunction.h>
#include <smallVec.h>
#include <baseFunctionSet.h>

namespace nb {


template <qc::Dimension>
class DecrementDim { };

template <>
class DecrementDim<qc::QC_3D> {
public:
  static const qc::Dimension d = qc::QC_2D;
};
template <>
class DecrementDim<qc::QC_2D> {
public:
  static const qc::Dimension d = qc::QC_1D;
};




//! Interface class for computing \f$ \int_{\partial \Omega} f \psi \, da \f$. The function \f$ f \f$
//! has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
template <typename ConfigType, typename GridType, typename VecType, qc::Dimension, typename QuadType, typename Imp>
class BoundaryIntegrationInterface : public aol::Op<VecType, VecType> { };

//! Interface class for computing \f$ \int_{\partial \Omega} f \psi \, da \f$. The function \f$ f \f$
//! has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
//! 3d-specification
template <typename ConfigType, typename GridType, typename VecType, typename QuadType, typename Imp>
class BoundaryIntegrationInterface<ConfigType, GridType, VecType, qc::QC_3D, QuadType, Imp> : public aol::Op<VecType, VecType> {
protected:
  const ConfigType &_config;
  const GridType &_grid;
public:
  typedef typename ConfigType::RealType RealType;
  typedef typename ConfigType::ElementType ElementType;
  typedef typename ConfigType::FaceType FaceType;

  BoundaryIntegrationInterface ( const ConfigType &config, GridType &grid )
      : _config ( config ), _grid ( grid ) {}

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec3<RealType> &normal, const ElementType &/*el*/, const aol::Vec3<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {
    // This lookUp-table represents the four corners of the the unit-square and is needed
    // for the computation of the boundary-vertex.
    static RealType lookUp[8][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };

    // eiterator is a face-iterator that traverses the boundary
    for ( typename GridType::eiterator eit = _grid.ebegin(); eit != _grid.eend(); ++eit ) {
      // _locNode stores, which nodes of the element belong to the boundary
      FaceType face = *eit;
      aol::Vec3<RealType> edgePt1 ( lookUp[face._locNode[0]][0], lookUp[face._locNode[0]][1], lookUp[face._locNode[0]][2] );
      aol::Vec3<RealType> edgePt2 ( lookUp[face._locNode[1]][0], lookUp[face._locNode[1]][1], lookUp[face._locNode[1]][2] );
      aol::Vec3<RealType> edgePt3 ( lookUp[face._locNode[2]][0], lookUp[face._locNode[2]][1], lookUp[face._locNode[2]][2] );
      aol::Vec3<RealType> edgePt4 ( lookUp[face._locNode[3]][0], lookUp[face._locNode[3]][1], lookUp[face._locNode[3]][2] );

      // computation of the normal
      aol::Vec3<RealType> edgePt ( edgePt1 );
      edgePt += edgePt2;
      edgePt += edgePt3;
      edgePt += edgePt4;
      edgePt /= 4.;

      aol::Vec3<RealType> normal ( edgePt[0] - 0.5, edgePt[1] - 0.5, edgePt[2] - 0.5 );
      normal *= 2.;

      QuadType quad;
      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate the 3d refcoords of the 2d quadpoint ***/
        const aol::Vec2<RealType> &rc = quad.getRefCoord ( q );
        for ( int i = 0; i < 3; i++ ) {
          edgePt[i] = edgePt1[i] * ( 1. - rc[0] ) * ( 1. - rc[1] ) +
                      edgePt2[i] * rc[0] * ( 1. - rc[1] ) +
                      edgePt3[i] * ( 1. - rc[0] ) * rc[1] +
                      edgePt4[i] * rc[0] * rc[1];
        }

        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, face._el, edgePt );

//         for ( int b=0; b<_config.getNumLocalDofs( face._el ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 4; nodeIndex++ ) {
          int b = face._locNode[nodeIndex];
          dest[ _config.localToGlobal ( face._el, b ) ] += quad.getWeight ( q ) * integrand * aol::Sqr ( _config.H ( face._el ) )
                                                           * _config.getBaseFunctionSet ( face._el ).evaluate ( b, edgePt );
        }
      }
    }
  }

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};

//! Interface class for computing \f$ \int_{\partial \Omega} f \psi \, da \f$. The function \f$ f \f$
//! has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
//! 2d-specification
template <typename ConfigType, typename GridType, typename VecType, typename QuadType, typename Imp>
class BoundaryIntegrationInterface<ConfigType, GridType, VecType, qc::QC_2D, QuadType, Imp> : public aol::Op<VecType, VecType> {
protected:
  const ConfigType &_config;
  const GridType &_grid;
public:
  typedef typename ConfigType::RealType RealType;
  typedef typename ConfigType::ElementType ElementType;

  BoundaryIntegrationInterface ( const ConfigType &config, GridType &grid )
      : _config ( config ), _grid ( grid ) {}

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec2<RealType> &normal, const ElementType &/*el*/, const aol::Vec2<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {
    // This lookUp-table represents the four corners of the the unit-square and is needed
    // for the computation of the boundary-vertex.
    static RealType lookUp[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };

    // eiterator is an face-iterator that traverses the boundary
    for ( typename GridType::eiterator eit = _grid.ebegin(); eit != _grid.eend(); ++eit ) {
      // _locNode stores, which nodes of the element belong to the boundary
      aol::Vec2<RealType> edgePt1 ( lookUp[eit->_locNode[0]][0], lookUp[eit->_locNode[0]][1] );
      aol::Vec2<RealType> edgePt2 ( lookUp[eit->_locNode[1]][0], lookUp[eit->_locNode[1]][1] );

      // computation of the normal
      aol::Vec2<RealType> edgePt ( edgePt1 );
      edgePt += edgePt2;
      edgePt /= 2.;

      aol::Vec2<RealType> normal ( edgePt[0] - 0.5, edgePt[0] - 0.5 );
      normal *= 2;

      QuadType quad;
      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate the 2d-refcoords from the 1d-quadpoint ***/
        const RealType rc = quad.getRefCoord ( q )[0];
        for ( int i = 0; i < 2; i++ ) {
          edgePt[i] = edgePt1[i] * ( 1. - rc )  +
                      edgePt2[i] * rc;
        }

        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, eit->_el, edgePt );

//         for ( int b=0; b<_config.getNumLocalDofs( eit->_el ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 2; nodeIndex++ ) {
          int b = eit->_locNode[nodeIndex];
          dest[ _config.localToGlobal ( eit->_el, b ) ] += quad.getWeight ( q ) * integrand * _config.H ( eit->_el )
                                                           * _config.getBaseFunctionSet ( eit->_el ).evaluate ( b, edgePt );
        }
      }
    }
  }

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};


/**
 * This class computes the following integral:
 * \f[ \int_{\partial\Omega} \tilde \gamma_z(\nabla \phi, \epsilon) \vec n \varphi da \f]
 *
 */
template <typename RealType, typename ConfigType, typename SubGridType, typename AnisoType, qc::Dimension Dim>
class IntegrateAnisotropyOverBoundary :
      public BoundaryIntegrationInterface < ConfigType, SubGridType, aol::Vector<RealType>,
      Dim, aol::GaussQuadrature<RealType, DecrementDim<Dim>::d, 7 > ,
      IntegrateAnisotropyOverBoundary < RealType, ConfigType, SubGridType,
      AnisoType, Dim > > {
  RealType _eps;
  const AnisoType &_anisotropy;
public:
  IntegrateAnisotropyOverBoundary ( const ConfigType &config, SubGridType &subGrid, RealType Epsilon,
                                    const AnisoType &anisotropy )
      : BoundaryIntegrationInterface < ConfigType, SubGridType, aol::Vector<RealType>,
      Dim, aol::GaussQuadrature<RealType, DecrementDim<Dim>::d, 7 > ,
      IntegrateAnisotropyOverBoundary < RealType, ConfigType, SubGridType,
      AnisoType, Dim > > ( config, subGrid ),
      _eps ( Epsilon ), _anisotropy ( anisotropy ) {}

  typedef typename ConfigType::ElementType ElementType;
  typedef typename ConfigType::VecType VecType;

  //! computes \f[ \tilde \gamma_z(\nabla \phi, \epsilon) \f]
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &arg, const VecType &normal,
                                  const ElementType &el, const VecType edgePoint ) const {
    // compute gradient of boundaryfunction @ edgePt
    aol::DiscreteFunctionDefault<ConfigType> discrU ( this->_grid, arg );
    VecType gradU;
    VecType gammaZ;

    discrU.evaluateGradient ( el, edgePoint, gradU );
    _anisotropy.gammaFirstDerivative ( gradU, gammaZ );   // gamma_z( grad U )

    return gammaZ * normal;
  }
};


/**
 * This class computes the following integral:
 * \f[ \int_{\partial\Omega} \frac{\nabla \Phi}{\| \nabla \Phi \|} \cdot n \vartheta da \f]
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename SubGridType>
class IntegrateMCMGradUOverBoundary :
      public BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateMCMGradUOverBoundary<ConfiguratorType, BoundaryQuadratureType, SubGridType> > {
  typename ConfiguratorType::RealType _epsSqr;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrPhiOld;

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;

  IntegrateMCMGradUOverBoundary ( const ConfiguratorType &config, SubGridType &subGrid, RealType Epsilon )
      : BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<RealType>,
      ConfiguratorType::Dim, BoundaryQuadratureType,
      IntegrateMCMGradUOverBoundary< ConfiguratorType, BoundaryQuadratureType, SubGridType> > ( config, subGrid ),
      _epsSqr ( Epsilon*Epsilon ), _discrPhiOld( NULL ) {}

  ~IntegrateMCMGradUOverBoundary ( ) {
    if ( _discrPhiOld ) delete _discrPhiOld;
  }

  void setImageReference ( const aol::Vector<RealType> &Phi ) {
    if ( static_cast< int >( Phi.size() ) != this->_config.getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Phi.size() << " _config.getNumGlobalDofs() = " << this->_config.getNumGlobalDofs()  << endl;
      throw aol::Exception ( "IntegrateProjectionUOverBoundary: Array Phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrPhiOld )
      delete _discrPhiOld;
    _discrPhiOld = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->_config.getInitializer(), Phi );
  }

  //! computes \f$ \frac{\nabla \Phi}{\| \nabla \Phi \|} \cdot n \f$
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &arg, const VecType &normal,
                                  const ElementType &el, const VecType edgePoint ) const {
    if ( !_discrPhiOld ) {
      throw aol::Exception ( "IntegrateMCMGradUOverBoundary::integrandAtQuadPoint: set image first!", __FILE__, __LINE__ );
    }

    // compute gradient of boundaryfunction @ edgePt
    aol::DiscreteFunctionDefault<ConfiguratorType> discrPhi ( this->_grid, arg );
    VecType gradPhi, gradPhiOld;
    discrPhi.evaluateGradient ( el, edgePoint, gradPhi );
    _discrPhiOld->evaluateGradient ( el, edgePoint, gradPhiOld );

    RealType norm = sqrt ( gradPhiOld.normSqr() + _epsSqr );       // sqrt( nabla PhiOld^2 + eps^2 )
    gradPhi /= norm;

    return gradPhi * normal;
  }
};



/**
 * This class computes the following integral:
 * \f[ \int_{\partial\Omega} P[\Phi] \nabla u \cdot n \vartheta da \f]
 *  \f$ P[\Phi]  \f$ is the projection on the tangent space defined by  \f$ \Phi. \f$
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename SubGridType>
class IntegrateProjectionUOverBoundary :
      public BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateProjectionUOverBoundary<ConfiguratorType, BoundaryQuadratureType, SubGridType> > {
  typename ConfiguratorType::RealType _epsSqr;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrPhi;

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  IntegrateProjectionUOverBoundary ( const ConfiguratorType &config, SubGridType &subGrid, RealType Epsilon )
      : BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<RealType>,
      ConfiguratorType::Dim, BoundaryQuadratureType,
      IntegrateProjectionUOverBoundary< ConfiguratorType, BoundaryQuadratureType, SubGridType> > ( config, subGrid ),
      _epsSqr ( Epsilon*Epsilon ), _discrPhi ( NULL ) {}

  ~IntegrateProjectionUOverBoundary ( ) {
    if ( _discrPhi ) delete _discrPhi;
  }

  void setImageReference ( const aol::Vector<RealType> &Phi ) {
    if ( static_cast<int> ( Phi.size() ) != this->_config.getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Phi.size() << " _config.getNumGlobalDofs() = " << this->_config.getNumGlobalDofs()  << endl;
      throw aol::Exception ( "IntegrateProjectionUOverBoundary: Array Phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrPhi )
      delete _discrPhi;
    _discrPhi = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->_config.getInitializer(), Phi );
  }


  //! computes \f[ \| \nabla \Phi \| P[\Phi] \nabla U \cdot \nu \f]
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &arg, const VecType &normal,
                                  const ElementType &el, const VecType edgePoint ) const {
    if ( !_discrPhi ) {
      throw aol::Exception ( "IntegrateProjectionUOverBoundary::integrandAtQuadPoint: set image first!", __FILE__, __LINE__ );
    }

    // compute gradient of boundaryfunction @ edgePt
    aol::DiscreteFunctionDefault<ConfiguratorType> discrU ( this->_grid, arg );
    VecType gradPhi, gradU, PGradU;
    discrU.evaluateGradient ( el, edgePoint, gradU );
    _discrPhi->evaluateGradient ( el, edgePoint, gradPhi );

    RealType normGradPhi = sqrt ( gradPhi.normSqr() + _epsSqr );       // sqrt( nabla Phi^2 + eps^2 )
    gradPhi /= normGradPhi;

    MatType mat;

    // compute the projection: Id - (nabla Phi)/||...|| \times (nabla Phi)/||...||
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      for ( int j = 0; j < ConfiguratorType::VecType::dim; j++ )
        mat[i][j] = -gradPhi[i] * gradPhi[j];
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      mat[i][i] += 1.;

    mat.mult ( gradU, PGradU );

    return normGradPhi * PGradU * normal;
  }
};


/**
 * This class computes the following integral:
 * \f[ \int_{\partial\Omega} x \cdot n \varphi da \f]
 * If you devide the summed result vector by the dimension you should
 * get the volume of the domain.
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename SubGridType>
class IntegrateIdOverBoundary :
      public BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateIdOverBoundary<ConfiguratorType, BoundaryQuadratureType, SubGridType> > {

  SubGridType _grid;                          // just for getting h

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;

  IntegrateIdOverBoundary ( const ConfiguratorType &config, SubGridType &subGrid )
      : BoundaryIntegrationInterface < ConfiguratorType, SubGridType, aol::Vector<RealType>,
      ConfiguratorType::Dim, BoundaryQuadratureType,
      IntegrateIdOverBoundary< ConfiguratorType, BoundaryQuadratureType, SubGridType> > ( config, subGrid ),
      _grid ( subGrid ) { }

  //! computes just \f[ x \cdot \nu \f]
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &/*arg*/, const VecType &normal,
                                  const ElementType &el, const VecType edgePoint ) const {
    VecType x ( el.x(), el.y(), el.z() );
    x += edgePoint;
    x *= _grid.H();

    return x * normal;
  }
};


} // end namespace nb


#endif
