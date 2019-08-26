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

#ifndef __DEFORMATIONENERGIES_H
#define __DEFORMATIONENERGIES_H

#include <aol.h>
#include <triangMeshConfigurators.h>
#include <triMesh.h>
#include "omp.h"

// basic headers
#include "geometry.h"
#include "auxiliaryOps.h"


//!===============================================================================================================================

using namespace aol;
using namespace qc;
using namespace om;

//!==========================================================================================================
//! MEMBRANE ENERGY (with hyperelastic, nonlinear density)
//!==========================================================================================================
//! hyperelastic energy density with log-term to penalize degeneration of triangles
template <typename ConfiguratorType>
 class HyperelasticEnergyDensityLog{
 protected:
   typedef typename ConfiguratorType::RealType RealType;
   // Lame constants and weight of term preventing complete volume compression (which in fact does not belong to St Venant-Kirchhoff)
   const RealType _lambdaQuarter, _muHalf, _dimMuHalf, _muHalfPlusLambdaQuarter;
   
 public:
   HyperelasticEnergyDensityLog( const RealType Mu, const RealType Lambda, const qc::Dimension dim = qc::QC_2D) :
     _lambdaQuarter( Lambda / 4. ),
     _muHalf( Mu / 2. ),
     _dimMuHalf( dim*_muHalf ),
     _muHalfPlusLambdaQuarter(_muHalf + _lambdaQuarter) {}
 
   inline RealType evaluate ( const RealType tr, 
                              const RealType det ) const {
        //cerr<<"old: tr = "<<tr<<", det = "<<det<<endl;
     return _muHalf *  tr + _lambdaQuarter * aol::Abs( det ) - _muHalfPlusLambdaQuarter * std::log( det ) - _dimMuHalf - _lambdaQuarter;
   }
 
   inline RealType evaluate ( const aol::Matrix22<RealType>& G ) const {
     return _muHalf *  G.tr() + _lambdaQuarter * aol::Abs( G.det() ) - _muHalfPlusLambdaQuarter * std::log( G.det() ) - _dimMuHalf - _lambdaQuarter;
   }

};

//! \brief Membrane energy between two thin shells given as triangular meshes.
//! \author Heeren
//!
//! Membrane energy between two shells $S$ and $\tilde S$ given by $E[S, \tilde S] = \int_S W(G[\phi]) da$, where 
//!   - $\phi: S \rightarrow \tilde S $ denotes the unique piecewise linear deformation
//!   - $G[\phi] = g^{-1}g_\phi \in R^{2,2}$ is Cauchy-Green tensor with g being the metric on $S$ and $g_\phi$ being the metric on $\tilde S = \phi(S)$
//!   - $W()$ denotes a hyperelastic energy density given e.g. by class HyperelasticEnergyDensityLog<> above
//!
//! Note that the energy might either be thought of as $S \mapsto E[S, \tilde S]$ (active shell is undeformed shell) or $\tilde S \mapsto E[S, \tilde S]$ (active shell is deformed shell).
//! The active shell is considered the argument whereas the inactive shell is given in the constructor.

template< typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneEnergy : public aol::Op< aol::MultiVector<typename TriMeshConfType::RealType>, aol::Scalar<typename TriMeshConfType::RealType> > {

public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _inactiveGeometry;
  const MaterialLawType _elasticEnergyDensity;
  const bool _activeShellIsDeformed;

public:

  MembraneEnergy( const MeshTopologySaver<MeshType>& topology,
                  const aol::ParameterParser& pparser,
                  const VectorType& InactiveGeometry,
		  const bool ActiveShellIsDeformed ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _elasticEnergyDensity( pparser.getDoubleOrDefault("lengthWeight", 1.0), pparser.getDoubleOrDefault("volWeight", 1.0) ), 
    _activeShellIsDeformed( ActiveShellIsDeformed) {}

  // energy evaluation
  void applyAdd( const VectorType& ActiveGeometry, aol::Scalar<RealType>& Dest ) const {

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    for ( int faceIdx = 0; faceIdx < this->_topology.getNumFaces(); ++faceIdx ){

      aol::RandomAccessContainer<Vec> P(3), Q(3);
      for ( int i = 0; i < 3; i++ ) {
        int idx = this->_topology.getNodeOfTriangle( faceIdx, i );
        defShellP->getTo( idx, P[i] );
        undefShellP->getTo( idx, Q[i] );
      }

      aol::Matrix22< RealType > G, temp;
      getMetric( Q, G );
      RealType area = std::sqrt( G.det() );
      G.invert();
      getMetric( P, temp );
      G *= temp;
    
      // NOTE here we already integrate (area of unit triangle is 0.5)!
      Dest[0] += _elasticEnergyDensity.evaluate( G.tr(), G.det() ) * 0.5 * area;
      
    }
  }

 };
 
//==========================================================================================================
//! \brief First derivative of membrane energy w.r.t. the deformed configuration (cf. class MembraneEnergy above).
//! \author Heeren
template< typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneGradientDef :
    public UnitTriangleCenterQuadFEVectorInterface< typename TriMeshConfType::InitType, MembraneGradientDef<TriMeshConfType> >{
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::InitType    TriMeshType;
  typedef typename TriMeshConfType::RealType    RealType;
  typedef aol::MultiVector<RealType>            VectorType;
  typedef aol::Vec3<RealType>                   Vec;
  
  const MeshTopologySaver<TriMeshType>& _topology;
  const aol::MultiVector<RealType>&  _undefShell;
  const RealType _mu, _lambdaHalf;

public:
  MembraneGradientDef( const MeshTopologySaver<TriMeshType>& topology,
                       const aol::ParameterParser& pparser,
                       const VectorType& undefShell ) :
      UnitTriangleCenterQuadFEVectorInterface< typename TriMeshConfType::InitType, MembraneGradientDef<TriMeshConfType> >( topology ),
      _topology( topology ),
      _undefShell( undefShell),
      _mu( pparser.getDoubleOrDefault("lengthWeight", 1.0) ),
      _lambdaHalf( pparser.getDoubleOrDefault("volWeight", 1.0)/2.) {}

  void getNonlinearity ( const VectorType &Arg,
                         const int &ElIdx,
                         aol::Mat<3, 2, RealType> &Matrix) const {

    aol::Matrix22< RealType >  gDef, gRef;
    
    aol::RandomAccessContainer<Vec> P(3), Q(3);
    for ( int i = 0; i < 3; i++ ) {
      int idx = this->_topology.getNodeOfTriangle( ElIdx, i );
      _undefShell.getTo( idx, Q[i] );
      Arg.getTo( idx, P[i] );
    }

    getMetric( P, gDef );
    getMetric( Q, gRef );
    getDx( P, Matrix );    
    Matrix *= -1. * std::sqrt( gRef.det() ); // why -1 ???

    gRef.invert();
    RealType factor = _lambdaHalf * gRef.det() * gDef.det() - _mu - _lambdaHalf;
    gDef.invert();
      
    gRef *= _mu;
    gRef.addMultiple( gDef, factor );
    Matrix *= gRef;
  }

};
 
//!
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneSubHessianDef :
  public UnitTriangleCenterQuadFEMatrixInterface<typename TriMeshConfType::InitType, MembraneSubHessianDef<TriMeshConfType,_MaterialLawType> > {
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefGeometry, _defGeometry;
  const RealType _mu, _lambdaHalf;
  const MaterialLawType _elasticEnergyDensity;
  const RealType _factor;
  const int _k, _l;

public:
   MembraneSubHessianDef(     const MeshTopologySaver<MeshType>& Topology,
                              const aol::ParameterParser& Parser,
                              const VectorType& UndefGeometry,
			      const VectorType& DefGeometry,
			      const RealType Factor,
                              const int K,
                              const int L ) :
    UnitTriangleCenterQuadFEMatrixInterface<MeshType, MembraneSubHessianDef<TriMeshConfType,MaterialLawType> > (  Topology ),
    _topology( Topology ),
    _undefGeometry( UndefGeometry ),
    _defGeometry( DefGeometry ),
    _mu( Parser.getDoubleOrDefault("lengthWeight", 1.0) ),
    _lambdaHalf( Parser.getDoubleOrDefault("volWeight", 1.0)/2.),
    _elasticEnergyDensity( _mu, 2.0*_lambdaHalf ),
    _factor( Factor ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const int& ElementIdx, aol::Mat<2, 2, RealType> &Matrix ) const
  {
    aol::Mat< 3, 2, RealType > DxDefgDefInv, DxDef;
    aol::Mat< 3, 3, RealType > auxMat33;
    aol::Matrix22< RealType >  gDefInv, gRef, tensorProduct;
    
    aol::RandomAccessContainer<Vec> P(3), Q(3);
    for ( int i = 0; i < 3; i++ ) {
      int idx = this->_topology.getNodeOfTriangle( ElementIdx, i );
      _undefGeometry.getTo( idx, Q[i] );
      _defGeometry.getTo( idx, P[i] );
    }
    
    getMetric( P, gDefInv );
    getMetric( Q, gRef );
    getDx( P, DxDef );
       
    RealType detgRef = gRef.det();
    RealType detG( gDefInv.det() / detgRef );
    RealType minAlpha = _mu + _lambdaHalf - _lambdaHalf * detG;
    gDefInv.invert();
    DxDefgDefInv.makeProduct( DxDef, gDefInv );

    Matrix.setZero();
    // \lambda * \det\G * \tr\( A D\phi) * \tr( A D\psi \) - \alpha \tr( A D\phi A D\psi ), where A = g_2^{-1} Dx_2^T
    tensorProduct.makeTensorProduct( DxDefgDefInv[_k], DxDefgDefInv[_l] );
    Matrix.addMultiple( tensorProduct, minAlpha );
    tensorProduct.transpose();
    Matrix.addMultiple( tensorProduct, 2.* _lambdaHalf * detG );
     
    //  - \alpha \tr( g_2^{-1}D\phi^T B D\psi ), where B = Dx_2 g_2^{-1} Dx_2^T
    auxMat33.makeProductABtransposed( DxDefgDefInv, DxDef );    
    Matrix.addMultiple( gDefInv,  minAlpha * auxMat33.get(_k,_l) );  
    
    // \tr( dW D\phi^T \id D\psi ), dW = \mu g_1^{-1} + alpha g_2^{-1}
    if( _k == _l ){
      gRef.invert();
      Matrix.addMultiple( gRef, _mu );
      Matrix.addMultiple( gDefInv, -1. * minAlpha );
    }

    // transpose (why??)
    Matrix.transpose();
    Matrix *= _factor * std::sqrt( detgRef );

  }
};

//==========================================================================================================
//! \brief Second derivative of membrane energy w.r.t. the deformed configuration (cf. class MembraneEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType>, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class MembraneHessianDef
      : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {
	      
public:
  typedef _MaterialLawType MaterialLawType;
    
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType    MeshType;
  typedef aol::MultiVector<RealType>            VectorType;
  typedef aol::Vec3<RealType>                   Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::ParameterParser& _pparser;
  const VectorType& _undefGeometry;
  const RealType _factor;

public:
  MembraneHessianDef ( const MeshTopologySaver<MeshType>& Topology,
                       const aol::ParameterParser& Parser,
                       const VectorType& UndefGeometry,
		       const RealType Factor = 1. ) :
    _topology( Topology ),
    _pparser( Parser ),
    _undefGeometry( UndefGeometry ),
    _factor( Factor ){}


  void applyAdd( const VectorType &Arg, BlockMatrixType &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < TriMeshConfType::Dim; k++ )
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int l = 0; l < TriMeshConfType::Dim; l++ )
        MembraneSubHessianDef<TriMeshConfType,MaterialLawType>( _topology, _pparser, _undefGeometry, Arg, _factor, k, l ).assembleAddMatrix( Dest.getReference( k, l ) );
  }

  void assemble( const VectorType &Arg, BlockMatrixType &Dest ) const {
    apply( Arg, Dest );
  }

};
//==========================================================================================================

//==========================================================================================================
//! \brief First derivative of membrane energy w.r.t. the undeformed configuration (cf. class MembraneEnergy above).
//! \author Heeren
template< typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneGradientUndef :
    public UnitTriangleCenterQuadFEVectorInterface< typename TriMeshConfType::InitType, MembraneGradientUndef<TriMeshConfType> >{
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::InitType    TriMeshType;
  typedef typename TriMeshConfType::RealType    RealType;
  typedef aol::MultiVector<RealType>            VectorType;
  typedef aol::Vec3<RealType>                   Vec;
  
  const MeshTopologySaver<TriMeshType>& _topology;
  const aol::MultiVector<RealType>&  _defShell;
  const RealType _mu, _lambdaHalf;
  const MaterialLawType _elasticEnergyDensity;

public:
  MembraneGradientUndef( const MeshTopologySaver<TriMeshType>& topology,
                         const aol::ParameterParser& pparser,
                         const VectorType& defShell ) :
      UnitTriangleCenterQuadFEVectorInterface< typename TriMeshConfType::InitType, MembraneGradientUndef<TriMeshConfType> >( topology ),
      _topology( topology ),
      _defShell( defShell),
      _mu( pparser.getDoubleOrDefault("lengthWeight", 1.0) ),
      _lambdaHalf( pparser.getDoubleOrDefault("volWeight", 1.0)/2.),
      _elasticEnergyDensity( _mu, 2.*_lambdaHalf ){}

  void getNonlinearity ( const VectorType &Arg,
                         const int &ElIdx,
                         aol::Mat<3, 2, RealType> &Matrix) const {

    aol::Matrix22< RealType >  gRefInv, G;
    
    aol::RandomAccessContainer<Vec> P(3), Q(3);
    for ( int i = 0; i < 3; i++ ) {
      int idx = this->_topology.getNodeOfTriangle( ElIdx, i );
      _defShell.getTo( idx, P[i] );
      Arg.getTo( idx, Q[i] );
    }
    
    getMetric( P, G );
    getMetric( Q, gRefInv );
    getDx( Q, Matrix );    
    Matrix *= std::sqrt( gRefInv.det() ); 
    gRefInv.invert();
  
    G.leftMult( gRefInv );    
    RealType detG( G.det() ), traceG( G.tr() );

    G *= _mu;
    G.addToDiagonal( _lambdaHalf * detG - _mu - _lambdaHalf - _elasticEnergyDensity.evaluate( traceG, detG ) );
    G *= gRefInv;
  
    Matrix *= G;    
  }

};

//!
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneSubHessianUndef :
  public UnitTriangleCenterQuadFEMatrixInterface<typename TriMeshConfType::InitType, MembraneSubHessianUndef<TriMeshConfType,_MaterialLawType> > {
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefGeometry, _defGeometry;
  const RealType _mu, _lambdaHalf;
  const MaterialLawType _elasticEnergyDensity;
  const RealType _factor;
  const int _k, _l;

public:
   MembraneSubHessianUndef( const MeshTopologySaver<MeshType>& Topology,
                            const aol::ParameterParser& Parser,
                            const VectorType& UndefGeometry,
			    const VectorType& DefGeometry,
			    const RealType Factor,
                            const int K,
                            const int L ) :
    UnitTriangleCenterQuadFEMatrixInterface< MeshType, MembraneSubHessianUndef<TriMeshConfType,MaterialLawType> > ( Topology ),
    _topology( Topology ),
    _undefGeometry( UndefGeometry ),
    _defGeometry( DefGeometry ),
    _mu( Parser.getDoubleOrDefault("lengthWeight", 1.0) ),
    _lambdaHalf( Parser.getDoubleOrDefault("volWeight", 1.0)/2.),
    _elasticEnergyDensity( _mu, 2.0*_lambdaHalf ),
    _factor( Factor ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const int& ElementIdx, aol::Mat< 2, 2, RealType > &Matrix ) const
  {
    aol::Mat< 3, 2, RealType > DxRef, auxMat1, auxMat2;
    aol::Matrix22< RealType >  G, gRefInv, auxMat22;
    aol::Mat< 3, 3, RealType > auxMat33;
    
    aol::RandomAccessContainer<Vec> P(3), Q(3);
    for ( int i = 0; i < 3; i++ ) {
      int idx = this->_topology.getNodeOfTriangle( ElementIdx, i );
      _undefGeometry.getTo( idx, Q[i] );
      _defGeometry.getTo( idx, P[i] );
    }
    
    getMetric( P, G );
    getMetric( Q, gRefInv );
    getDx( Q, DxRef );   
    RealType area = std::sqrt( gRefInv.det() ); 
    gRefInv.invert();
  
    // G
    G.leftMult( gRefInv );   
    RealType detG = G.det();
    // dW 
    aol::Matrix22< RealType >  dW( G );
    dW *= _mu;
    dW.addToDiagonal( _lambdaHalf * detG - _mu - _lambdaHalf - _elasticEnergyDensity.evaluate( G.tr(), detG ) );

    // A = g_1^{-1}; B = dW g_1^{-1}; C = g_1^{-1} g_2 g_1^{-1}
    // M_r = M Dx_1^T; M_l = Dx_1 M; M_{rl} = Dx_1 M Dx_1^T
    dW *= gRefInv;
    G *= gRefInv;

    // \tr( A_r D\phi B_l^T D\psi ) - \tr(A_r D\phi)\tr(B_r D\psi)
    auxMat1.makeProduct( DxRef, gRefInv );
    auxMat2.makeProduct( DxRef, dW );
    auxMat22.makeTensorProduct( auxMat1[_k], auxMat2[_l] );
    Matrix = auxMat22;
    auxMat22.transpose( );
    Matrix -= auxMat22;
    
    // tr( A D\phi^T B_{lr}^T D\psi ) 
    auxMat33.makeProductABtransposed( auxMat2, DxRef );
    Matrix.addMultiple( gRefInv, auxMat33.get(_l, _k) );
    
    // \alpha * \tr(A_r D\phi) \tr(A_r D\psi)
    auxMat22.makeTensorProduct( auxMat1[_l], auxMat1[_k] );
    Matrix.addMultiple( auxMat22, _lambdaHalf * detG + _mu + _lambdaHalf );
    
    // \mu * (  \tr(C_r D\phi A_r D\psi) - tr(C_r D\phi) \tr( A_r D\psi)  )
    auxMat2.makeProduct( DxRef, G );
    auxMat22.makeTensorProduct( auxMat2[_k], auxMat1[_l] );
    Matrix.addMultiple( auxMat22, _mu );
    auxMat22.transpose();
    Matrix.addMultiple( auxMat22, -1. * _mu );
    
    // \mu \tr(C D\phi^T A_{lr} D\psi)
    auxMat33.makeProductABtransposed( auxMat1, DxRef );
    Matrix.addMultiple( G, _mu * auxMat33.get(_k, _l) );
    
    // - \tr(B^T D\phi^T \id D\psi)
    if( _k == _l )
      Matrix -= dW.transposed();  

    // transpose (why??)
    Matrix.transpose();
    Matrix *= _factor * area;
  }
};

//! \brief Second derivative of membrane energy w.r.t. the undeformed configuration (cf. class MembraneEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType>, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class MembraneHessianUndef
      : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType    MeshType;
  typedef aol::MultiVector<RealType>            VectorType;
  typedef aol::Vec3<RealType>                   Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::ParameterParser& _pparser;
  const VectorType& _defGeometry;
  const RealType _factor;

public:
  MembraneHessianUndef ( const MeshTopologySaver<MeshType>& Topology,
                         const aol::ParameterParser& Parser,
                         const VectorType& DefGeometry,
			 const RealType Factor = 1.  ) :
    _topology( Topology ),
    _pparser( Parser ),
    _defGeometry( DefGeometry ),
    _factor( Factor ){}


  void applyAdd( const VectorType &Arg, BlockMatrixType &Dest ) const {
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < TriMeshConfType::Dim; k++ )
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int l = 0; l < TriMeshConfType::Dim; l++ )
        MembraneSubHessianUndef<TriMeshConfType,MaterialLawType>( _topology, _pparser, Arg, _defGeometry, _factor, k, l ).assembleAddMatrix( Dest.getReference( k, l ) );
  }

  void assemble( const VectorType &Arg, BlockMatrixType &Dest ) const {
    apply( Arg, Dest );
  }

};

//==========================================================================================================
//! wrapper in order to test mixed second derivatives
//! S_2 \mapsto \partial_1 W[S_1, S_2]
template <typename TriMeshConfType>
class MembraneGradientUndefWrapper : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::ParameterParser& _pparser;
  const VectorType& _undefShell;

public:
MembraneGradientUndefWrapper( const MeshTopologySaver<MeshType>& topology,
			      const aol::ParameterParser& pparser,
                              const VectorType& undefShell ) : _topology( topology), _pparser( pparser ), _undefShell(undefShell) {}



void applyAdd( const VectorType& defShell, VectorType& Dest ) const {
  MembraneGradientUndef<TriMeshConfType>( _topology, _pparser, defShell ).applyAdd( _undefShell, Dest );
}

};

//! wrapper in order to test mixed second derivatives
//! S_1 \mapsto \partial_2 W[S_1, S_2]
template <typename TriMeshConfType>
class MembraneGradientDefWrapper : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::ParameterParser& _pparser;
  const VectorType& _defShell;

public:
MembraneGradientDefWrapper( const MeshTopologySaver<MeshType>& topology,
			    const aol::ParameterParser& pparser,
                            const VectorType& defShell ) : _topology( topology), _pparser( pparser ), _defShell(defShell) {}


void applyAdd( const VectorType& UndefShell, VectorType& Dest ) const {
  MembraneGradientDef<TriMeshConfType>( _topology, _pparser, UndefShell ).applyAdd( _defShell, Dest );
}

};

//!
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneSubHessianMixed :
  public UnitTriangleCenterQuadFEMatrixInterface<typename TriMeshConfType::InitType, MembraneSubHessianMixed<TriMeshConfType,_MaterialLawType> > {
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefGeometry, _defGeometry;
  const bool _firstDerivWRTUndef;
  const RealType _mu, _lambdaHalf;
  const MaterialLawType _elasticEnergyDensity;
  const RealType _factor;
  const int _k, _l;

public:
   MembraneSubHessianMixed(   const MeshTopologySaver<MeshType>& Topology,
                              const aol::ParameterParser& Parser,
                              const VectorType& UndefGeometry,
			      const VectorType& DefGeometry,
			      const bool FirstDerivWRTUndef,
			      const RealType Factor,
                              const int K,
                              const int L ) :
    UnitTriangleCenterQuadFEMatrixInterface<MeshType, MembraneSubHessianMixed<TriMeshConfType,MaterialLawType> > ( Topology ),
    _topology( Topology ),
    _undefGeometry( UndefGeometry ),
    _defGeometry( DefGeometry ),
    _firstDerivWRTUndef( FirstDerivWRTUndef ),
    _mu( Parser.getDoubleOrDefault("lengthWeight", 1.0) ),
    _lambdaHalf( Parser.getDoubleOrDefault("volWeight", 1.0)/2.),
    _elasticEnergyDensity( _mu, 2.0*_lambdaHalf ),
    _factor( Factor ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const int& ElementIdx, aol::Mat<2, 2, RealType> &Matrix ) const
  {
    aol::Mat< 3, 2, RealType > DxRefgRefInv, DxDef, DxDefgRefInv;
    aol::Mat< 3, 3, RealType > DxRefgRefInvDxDefT;
    aol::Matrix22< RealType >  gDefInv, gRefInv, auxMat22;
    
    aol::RandomAccessContainer<Vec> P(3), Q(3);
    for ( int i = 0; i < 3; i++ ) {
      int idx = this->_topology.getNodeOfTriangle( ElementIdx, i );
      _undefGeometry.getTo( idx, Q[i] );
      _defGeometry.getTo( idx, P[i] );
    }
    
    getMetric( P, gDefInv );
    getMetric( Q, gRefInv );
    getDx( P, DxDef );
    getDx( Q, DxRefgRefInv );
       
    RealType detgRef = gRefInv.det();
    RealType detG( gDefInv.det() / detgRef );
    RealType beta = _lambdaHalf * detG + _mu + _lambdaHalf;
    
    gDefInv.invert();
    gRefInv.invert();    
    DxDefgRefInv.makeProduct( DxDef, gRefInv );
    DxRefgRefInv *= gRefInv;
    DxRefgRefInvDxDefT.makeProductABtransposed( DxRefgRefInv, DxDef );
    
    int k = _firstDerivWRTUndef ? _l : _k;
    int l = _firstDerivWRTUndef ? _k : _l;
    
    // \tr( dW Dx_2^T D\phi) * \tr( g_1^{-1} Dx_1^T  D\psi), where dW  = \mu g_1^{-1} + \alpha g_2^{-1}
    getWeightedMatrixSum( _mu, gRefInv, -1. * beta, gDefInv, auxMat22 );
    DxDef *= auxMat22;    
    Matrix.makeTensorProduct( DxRefgRefInv[l], DxDef[k] );
    
    // -\mu \tr(g_1^{-1} D\phi^T B D\psi ), where B = Dx_2 g_1^{-1} Dx_1^T (note B \neq B^T !)
    Matrix.addMultiple( gRefInv, -1. * _mu *  DxRefgRefInvDxDefT.get(l, k) );
    // -\mu \tr( g_1^{-1} Dx_2^T D\phi g_1^{-1} Dx_1^T D\psi ) 
    auxMat22.makeTensorProduct( DxDefgRefInv[k], DxRefgRefInv[l] );
    Matrix.addMultiple( auxMat22, -1. * _mu  );

    // transpose (why??)
    if( !_firstDerivWRTUndef )
      Matrix.transpose();
    Matrix *= _factor * std::sqrt( detgRef );

  }
};

//! \brief Second (mixed) derivative of membrane energy (cf. class MembraneEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename _MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType>, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class MembraneHessianMixed
      : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {
      
public:
  typedef _MaterialLawType MaterialLawType;
  
protected:
  typedef typename TriMeshConfType::RealType    RealType;
  typedef typename TriMeshConfType::InitType    MeshType;
  typedef aol::MultiVector<RealType>            VectorType;
  typedef aol::Vec3<RealType>                   Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::ParameterParser& _pparser;
  const VectorType& _inactiveGeometry;
  const bool _activeShellIsDeformed, _firstDerivWRTUndef;
  const RealType _factor;

public:
  MembraneHessianMixed ( const MeshTopologySaver<MeshType>& Topology,
                         const aol::ParameterParser& Parser,
                         const VectorType& InactiveGeometry,
			 const bool ActiveShellIsDeformed,
			 const bool FirstDerivWRTUndef,
			 const RealType Factor = 1.  ) :
    _topology( Topology ),
    _pparser( Parser ),
    _inactiveGeometry( InactiveGeometry ),
    _activeShellIsDeformed( ActiveShellIsDeformed ),
    _firstDerivWRTUndef( FirstDerivWRTUndef ),
    _factor( Factor ){}


  void applyAdd( const VectorType &ActiveGeometry, BlockMatrixType &Dest ) const {
    
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int k = 0; k < TriMeshConfType::Dim; k++ )
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int l = 0; l < TriMeshConfType::Dim; l++ )
        MembraneSubHessianMixed<TriMeshConfType,MaterialLawType>( _topology, _pparser, *undefShellP, *defShellP, _firstDerivWRTUndef, _factor, k, l ).assembleAddMatrix( Dest.getReference( k, l ) );
  }

  void assemble( const VectorType &Arg, BlockMatrixType &Dest ) const {
    apply( Arg, Dest );
  }

};
//==========================================================================================================



//!==========================================================================================================
//! DISCRETE SHELLS BENDING ENERGY (HINGE ENERGY)
//!==========================================================================================================

//! \brief Discrete bending energy between two thin shells given as triangular meshes.
//! \author Heeren
//!
//! Smooth bending energy between two shells $S$ and $\tilde S$ is given by $E[S, \tilde S] = \int_S ( h - h_\phi])^2 da$, 
//! where $h$ and $h_\phi$ denote the mean curvature on $S$ and $\tilde S = \phi(S)$ respectively.
//! Discrete bending energy is taken from "Discrete shells" paper by Grinspun et al., 2003.
//! Here $E[S, \tilde S] = \sum_e \frac{(\theta_e - \theta_{\phi(e)})^2 |e|^2}{A_e}$, 
//! where the sum is over all edges $e \in S$, where $\theta_e$ is the dihedral angle at edge $e$ and $A_e = |T_1| + |T_2|$, if $T_i$ are the triangles sharing edge $e$.
//!
//! Note that the energy might either be thought of as $S \mapsto E[S, \tilde S]$ (active shell is undeformed shell) or $\tilde S \mapsto E[S, \tilde S]$ (active shell is deformed shell).
//! The active shell is considered the argument whereas the inactive shell is given in the constructor.
template< typename TriMeshConfType >
class SimpleBendingEnergy : public aol::Op< aol::MultiVector<typename TriMeshConfType::RealType>, aol::Scalar<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _inactiveGeometry;
  const bool _activeShellIsDeformed;

public:

  SimpleBendingEnergy( const MeshTopologySaver<MeshType>& topology,
                       const VectorType& InactiveGeometry,
		       const bool ActiveShellIsDeformed ) 
  : _topology( topology), 
    _inactiveGeometry(InactiveGeometry), 
    _activeShellIsDeformed( ActiveShellIsDeformed) {}

  // energy evaluation
  void applyAdd( const VectorType& ActiveGeometry, aol::Scalar<RealType>& Dest ) const {

    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry : &_inactiveGeometry;
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    
    for ( int edgeIdx = 0; edgeIdx < this->_topology.getNumEdges(); ++edgeIdx ){

      int pi( this->_topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( this->_topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( this->_topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( this->_topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      aol::Vec3<RealType> Pi, Pj, Pk, Pl, temp;

      defShellP->getTo( pi, Pi ); defShellP->getTo( pj, Pj );
      defShellP->getTo( pk, Pk ); defShellP->getTo( pl, Pl );
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      undefShellP->getTo( pi, Pi ); undefShellP->getTo( pj, Pj );
      undefShellP->getTo( pk, Pk ); undefShellP->getTo( pl, Pl );

      // compute volume, length of edge and theta difference
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) );       
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      
      Dest[0] += delTheta * delTheta * elengthSqr / vol;

    }
  }

 };
   
//========================================================================================================== 
//! \brief First derivative of simple bending energy w.r.t. the undeformed configuration (cf. class SimpleBendingEnergy above).
//! \author Heeren
template< typename TriMeshConfType >
class SimpleBendingGradientUndef : public aol::Op< aol::MultiVector<typename TriMeshConfType::RealType> > {

  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _defShell;

public:
SimpleBendingGradientUndef( const MeshTopologySaver<MeshType>& topology,
                            const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}

SimpleBendingGradientUndef( const MeshType& /*dummy*/,
			    const MeshTopologySaver<MeshType>& topology,
                            const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}

void applyAdd( const VectorType& undefShell, VectorType& Dest ) const {

  assert( _defShell.getTotalSize() == undefShell.getTotalSize() );
  if( (Dest.numComponents() != 3) || ( Dest[0].size() != undefShell[0].size() ) )
    Dest.reallocate( 3, undefShell[0].size() );

  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
            pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      // first get defomed quantities
      Vec Pi, Pj, Pk, Pl, temp;
      _defShell.getTo( pi, Pi ); _defShell.getTo( pj, Pj );
      _defShell.getTo( pk, Pk ); _defShell.getTo( pl, Pl );
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );

      //!now get the undeformed values
      undefShell.getTo( pi, Pi ); undefShell.getTo( pj, Pj );
      undefShell.getTo( pk, Pk ); undefShell.getTo( pl, Pl );

      // compute Ak + Al, |e|^2 and diff. o dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) ); 
      // note signe here!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      

      // derivatives    
      Vec gradk, gradl, gradi, gradj, gradTheta, gradArea;
      RealType factorGradTheta = -2. * delTheta * elengthSqr / vol;    
      RealType factorGradArea = -1. * delTheta * delTheta * elengthSqr / (vol*vol);
      RealType factorGradEdgeLengthSqr = 2. * delTheta * delTheta / vol;   
         
      // d_k
      getThetaGradK( Pi, Pj, Pk, gradTheta );
      getAreaGradK( Pi, Pj, Pk, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradk );
      
      // d_l
      getThetaGradK( Pj, Pi, Pl, gradTheta );
      getAreaGradK( Pj, Pi, Pl, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradl );
      
      // d_i
      getThetaGradI( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pj, Pk, Pi, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradi );
      getAreaGradK( Pl, Pj, Pi, gradArea );
      gradi.addMultiple( gradArea, factorGradArea );
      gradi.addMultiple( Pi-Pj, factorGradEdgeLengthSqr );
      
      // d_j
      getThetaGradJ( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pk, Pi, Pj, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradj );
      getAreaGradK( Pi, Pl, Pj, gradArea );
      gradj.addMultiple( gradArea, factorGradArea );
      gradj.addMultiple( Pj-Pi, factorGradEdgeLengthSqr );      
      
      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i][pi] += gradi[i];
        Dest[i][pj] += gradj[i];
        Dest[i][pk] += gradk[i];
        Dest[i][pl] += gradl[i];
      }

  }
}

};

//! \brief Second derivative of simple bending energy w.r.t. the undeformed configuration (cf. class SimpleBendingEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class SimpleBendingHessianUndef : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _defShell;
  const RealType _factor;

public:
SimpleBendingHessianUndef( const MeshTopologySaver<MeshType>& topology,
                           const VectorType& defShell,
			   const RealType Factor = 1. ) : _topology( topology), _defShell(defShell), _factor( Factor ) {}


void applyAdd( const VectorType& undefShell, BlockMatrixType& Dest ) const {

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      Vec Pi, Pj, Pk, Pl, temp;

      // first get defomed quantities
      _defShell.getTo( pi, Pi ); _defShell.getTo( pj, Pj );
      _defShell.getTo( pk, Pk ); _defShell.getTo( pl, Pl );
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );

      
      // get undeformed vertex positions
      undefShell.getTo( pi, Pi ); undefShell.getTo( pj, Pj );
      undefShell.getTo( pk, Pk ); undefShell.getTo( pl, Pl );      
            
      // compute Ak + Al, |e|^2 and dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr = Vec(Pj-Pi).normSqr();  
      // note the sign!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );
      delTheta *= -1.;
      
      // compute first derivatives of dihedral angle
      Vec thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );
      
      // compute first derivatives of area
      Vec areak, areal, areai, areaj;
      getAreaGradK( Pi, Pj, Pk, areak );
      getAreaGradK( Pj, Pi, Pl, areal );
      getAreaGradK( Pj, Pk, Pi, areai );
      getAreaGradK( Pl, Pj, Pi, temp );
      areai += temp;
      getAreaGradK( Pk, Pi, Pj, areaj );
      getAreaGradK( Pi, Pl, Pj, temp );
      areaj += temp;

      // now compute second derivatives
      aol::Matrix33<RealType> H, auxMat;
      Vec auxVec, e(Pj-Pi);
      
      //*k
      getWeightedVectorSum( elengthSqr, thetak, -1. * delTheta * elengthSqr / vol, areak,  auxVec );
      getWeightedVectorSum( 2., thetak, -1. * delTheta / vol, areak,  temp );
      
      //kk      
      H.makeTensorProduct( thetak, auxVec );            
      getHessThetaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );     
      auxMat.makeTensorProduct( areak, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );  
      getHessAreaKK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );           
      H *= 2./vol;
      localToGlobal( Dest, pk, pk, H ); 
      
      //lk      
      H.makeTensorProduct( thetal, auxVec );            
      auxMat.makeTensorProduct( areal, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );               
      H *= 2./vol;
      localToGlobal( Dest, pl, pk, H ); 

      //ik
      H.makeTensorProduct( thetai, auxVec );            
      getHessThetaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );    
      getHessAreaIK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );     
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta );      
      H *= 2./vol;
      localToGlobal( Dest, pi, pk, H ); 

      //jk
      H.makeTensorProduct( thetaj, auxVec );            
      getHessThetaJK( Pi, Pj, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pj, Pi, Pk, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );    
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta );         
      H *= 2./vol;
      localToGlobal( Dest, pj, pk, H ); 
      
      //*l
      getWeightedVectorSum( elengthSqr, thetal, -1. * delTheta * elengthSqr / vol, areal,  auxVec );
      getWeightedVectorSum( 2., thetal, -1. * delTheta / vol, areal,  temp );
      
      //ll      
      H.makeTensorProduct( thetal, auxVec );            
      getHessThetaKK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );     
      auxMat.makeTensorProduct( areal, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol ); 
      getHessAreaKK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );           
      H *= 2./vol;
      localToGlobal( Dest, pl, pl, H ); 
      
      //il
      H.makeTensorProduct( thetai, auxVec );            
      getHessThetaJK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pi, Pj, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );     
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta );  
      H *= 2./vol;
      localToGlobal( Dest, pi, pl, H ); 
      
      //jl
      H.makeTensorProduct( thetaj, auxVec );            
      getHessThetaIK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );      
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );   
      getHessAreaIK( Pj, Pi, Pl, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );   
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta );  
      H *= 2./vol;
      localToGlobal( Dest, pj, pl, H ); 
      
      //*j
      getWeightedVectorSum( elengthSqr, thetaj, -1. * delTheta * elengthSqr / vol, areaj,  auxVec );
      auxVec.addMultiple( e, 2.*delTheta );// caution with factor 2!!!!!!!
      getWeightedVectorSum( 2., thetaj, -1. * delTheta / vol, areaj,  temp );
     
      //jj     
      H.makeTensorProduct( thetaj, auxVec );   
      getHessThetaII( Pj, Pi, Pl, Pk, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );   
      auxVec.addMultiple( e, -1.*delTheta );
      auxMat.makeTensorProduct( areaj, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaKK( Pk, Pi, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaKK( Pi, Pl, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol ); 
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, delTheta ); 
      H.addToDiagonal( delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( Dest, pj, pj, H );
      
      //ij     
      auxVec.addMultiple( e, delTheta );
      H.makeTensorProduct( thetai, auxVec );   
      getHessThetaJI( Pi, Pj, Pk, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr );
      auxVec.addMultiple( e, -1.*delTheta );
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaIK( Pi, Pk, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaIK( Pi, Pl, Pj, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );      
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta ); 
      H.addToDiagonal( -1. * delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( Dest, pi, pj, H );
      
      //*i
      getWeightedVectorSum( elengthSqr, thetai, -1. * delTheta * elengthSqr / vol, areai,  auxVec );
      auxVec.addMultiple( e, -2.*delTheta ); // caution with factor 2!!!!!!!
      getWeightedVectorSum( 2., thetai, -1. * delTheta / vol, areai,  temp );
      
      //ii     
      H.makeTensorProduct( thetai, auxVec );   
      getHessThetaII( Pi, Pj, Pk, Pl, auxMat );
      H.addMultiple( auxMat,  delTheta * elengthSqr ); 
      auxVec.addMultiple( e, 1.*delTheta );
      auxMat.makeTensorProduct( areai, auxVec );
      H.addMultiple( auxMat, -1. * delTheta / vol );       
      getHessAreaKK( Pl, Pj, Pi, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      getHessAreaKK( Pj, Pk, Pi, auxMat );
      H.addMultiple( auxMat, -0.5 * delTheta * delTheta * elengthSqr / vol );  
      auxMat.makeTensorProduct( e, temp );
      H.addMultiple( auxMat, -1.*delTheta ); 
      H.addToDiagonal( delTheta * delTheta );
      H *= 2./vol;
      localToGlobal( Dest, pi, pi, H );
     
    }
  }

protected:
  void localToGlobal( BlockMatrixType& BlockOp, int k, int l, const aol::Matrix33<RealType>& localMatrix ) const {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( k, l, _factor * localMatrix.get(i,j) );	
	
    if( k != l)
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( l, k, _factor * localMatrix.get(j,i) );	
  }
  
};

//==========================================================================================================
//! \brief First derivative of simple bending energy w.r.t. the deformed configuration (cf. class SimpleBendingEnergy above).
//! \author Heeren
template< typename TriMeshConfType >
class SimpleBendingGradientDef : public aol::Op< aol::MultiVector<typename TriMeshConfType::RealType> > {

  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   MeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef aol::Vec3<RealType>                  Vec;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _undefShell;

public:
SimpleBendingGradientDef( const MeshTopologySaver<MeshType>& topology,
                          const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}

SimpleBendingGradientDef( const MeshType& /*dummy*/,
			  const MeshTopologySaver<MeshType>& topology,
                          const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}

void applyAdd( const VectorType& defShell, VectorType& Dest ) const {

  assert( _undefShell.getTotalSize() == defShell.getTotalSize() );
  if( (Dest.numComponents() != 3) || ( Dest[0].size() != defShell[0].size() ) )
    Dest.reallocate( 3, defShell[0].size() );

  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
            pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      //! first get undefomed quantities
      Vec Pi, Pj, Pk, Pl, temp;
      _undefShell.getTo( pi, Pi ); _undefShell.getTo( pj, Pj );
      _undefShell.getTo( pk, Pk ); _undefShell.getTo( pl, Pl );
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      // compute Ak + Al
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) );     
      
      //! now get the deformed values
      defShell.getTo( pi, Pi ); defShell.getTo( pj, Pj );
      defShell.getTo( pk, Pk ); defShell.getTo( pl, Pl );

      // compute weighted differnce of dihedral angles
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );
      delTheta *= -2. * elengthSqr / vol ;

      // compute first derivatives of dihedral angle
      Vec thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      // assemble in global matrix
      for( int i = 0; i < 3; i++ ){
        Dest[i][pi] += delTheta * thetai[i];
        Dest[i][pj] += delTheta * thetaj[i];
        Dest[i][pk] += delTheta * thetak[i];
        Dest[i][pl] += delTheta * thetal[i];
      }

  }
}

};

//! \brief Second derivative of simple bending energy w.r.t. the deformed configuration (cf. class SimpleBendingEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class SimpleBendingHessianDef : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefShell;
  const RealType _factor;

public:
SimpleBendingHessianDef( const MeshTopologySaver<MeshType>& topology,
                         const VectorType& undefShell,
			 const RealType Factor = 1. ) : _topology( topology), _undefShell(undefShell), _factor( Factor ) {} 



void applyAdd( const VectorType& defShell, BlockMatrixType& Dest ) const {

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      Vec Pi, Pj, Pk, Pl, temp;

      // first get defomed quantities
      _undefShell.getTo( pi, Pi ); _undefShell.getTo( pj, Pj );
      _undefShell.getTo( pk, Pk ); _undefShell.getTo( pl, Pl );
      RealType delThetaDouble = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      // compute Ak + Al
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr = Vec(Pj-Pi).normSqr();  
      
      // get undeformed vertex positions
      defShell.getTo( pi, Pi ); defShell.getTo( pj, Pj );
      defShell.getTo( pk, Pk ); defShell.getTo( pl, Pl );
      
      // compute difference in dihedral angles
      delThetaDouble -= getDihedralAngle( Pi, Pj, Pk, Pl );
      delThetaDouble *= -2. * elengthSqr / vol;
      RealType factor = 2. * elengthSqr / vol;

      // compute first derivatives of dihedral angle
      Vec thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      // now compute second derivatives of dihedral angle
      aol::Matrix33<RealType> tensorProduct, H, aux;
      
      //kk
      getHessThetaKK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetak, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pk, pk, H ); 
            
      //ik & ki (Hki = Hik)
      getHessThetaIK( Pi, Pj, Pk, aux);
      tensorProduct.makeTensorProduct( thetai, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pi, pk, H );
      
      //jk & kj (Hkj = Hjk)
      getHessThetaJK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetaj, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pj, pk, H );       
      
      //ll
      getHessThetaKK( Pj, Pi, Pl, aux );
      tensorProduct.makeTensorProduct( thetal, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pl, pl, H );      
      
      //il & li (Hli = Hil)
      getHessThetaJK( Pj, Pi, Pl, aux);
      tensorProduct.makeTensorProduct( thetai, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pi, pl, H );
      
      //jl & lj (Hlj = Hjl)
      getHessThetaIK( Pj, Pi, Pl, aux);
      tensorProduct.makeTensorProduct( thetaj, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pj, pl, H );            
      
      //kl/lk: Hkl = 0 and Hlk = 0
      tensorProduct.makeTensorProduct( thetak, thetal );
      tensorProduct *= factor;
      localToGlobal( Dest, pk, pl, tensorProduct );
        
      //ii  
      getHessThetaII( Pi, Pj, Pk, Pl, aux );
      tensorProduct.makeTensorProduct( thetai, thetai );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pi, pi, H );              

      //jj
      getHessThetaII( Pj, Pi, Pl, Pk, aux );       
      tensorProduct.makeTensorProduct( thetaj, thetaj );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H);
      localToGlobal( Dest, pj, pj, H );

      //ij & ji (Hij = Hji)
      getHessThetaJI( Pi, Pj, Pk, Pl, H );     
      H *= delThetaDouble;
      tensorProduct.makeTensorProduct( thetai, thetaj );
      H.addMultiple( tensorProduct, factor );
      localToGlobal( Dest, pi, pj, H );  

    }
  }

protected:
  void localToGlobal( BlockMatrixType& BlockOp, int k, int l, const aol::Matrix33<RealType>& localMatrix ) const {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( k, l, _factor * localMatrix.get(i,j) );	
	
    if( k != l)
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( l, k, _factor * localMatrix.get(j,i) );	
  }
  
};
   
//==========================================================================================================
//!
template <typename TriMeshConfType>
class SimpleBendingGradientUndefWrapper : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _undefShell;

public:
SimpleBendingGradientUndefWrapper( const MeshTopologySaver<MeshType>& topology,
                                   const VectorType& undefShell ) : _topology( topology), _undefShell(undefShell) {}



void applyAdd( const VectorType& defShell, VectorType& Dest ) const {
  SimpleBendingGradientUndef<TriMeshConfType>( _topology, defShell ).applyAdd( _undefShell, Dest );
}

};

//!
template <typename TriMeshConfType>
class SimpleBendingGradientDefWrapper : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const VectorType& _defShell;

public:
SimpleBendingGradientDefWrapper( const MeshTopologySaver<MeshType>& topology,
                              const VectorType& defShell ) : _topology( topology), _defShell(defShell) {}


void applyAdd( const VectorType& UndefShell, VectorType& Dest ) const {
  SimpleBendingGradientDef<TriMeshConfType>( _topology, UndefShell ).applyAdd( _defShell, Dest );
}

};

//! \brief Second (mixed) derivative of simple bending energy (cf. class SimpleBendingEnergy above).
//! \author Heeren
template <typename TriMeshConfType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class SimpleBendingHessianMixed : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType MeshType;
  typedef aol::Vec3<RealType>                Vec;
  typedef aol::MultiVector<RealType>         VectorType;

  const MeshTopologySaver<MeshType>& _topology;
  const aol::MultiVector<RealType>&  _inactiveGeometry;
  const bool _activeShellIsDeformed, _firstDerivWRTUndef;
  const RealType _factor;

public:
SimpleBendingHessianMixed( const MeshTopologySaver<MeshType>& topology,
                           const VectorType& InactiveGeometry,
		           const bool ActiveShellIsDeformed,
			   const bool FirstDerivWRTUndef,
			   const RealType Factor = 1. ) : _topology( topology), _inactiveGeometry(InactiveGeometry), _activeShellIsDeformed(ActiveShellIsDeformed), _firstDerivWRTUndef( FirstDerivWRTUndef ), _factor( Factor ) {}

//
void applyAdd( const VectorType& ActiveGeometry, BlockMatrixType& Dest ) const {
  
    const VectorType* undefShellP = _activeShellIsDeformed ? &_inactiveGeometry : &ActiveGeometry;
    const VectorType* defShellP   = _activeShellIsDeformed ? &ActiveGeometry    : &_inactiveGeometry;
      
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( min( pl, pk) < 0 )
        continue;

      // set up vertices and edges
      Vec Pi, Pj, Pk, Pl;

      // first get defomed quantities
      defShellP->getTo( pi, Pi ); defShellP->getTo( pj, Pj );
      defShellP->getTo( pk, Pk ); defShellP->getTo( pl, Pl );
      RealType delTheta = getDihedralAngle( Pi, Pj, Pk, Pl );
      
      Vec defGradi, defGradj, defGradk, defGradl, temp;
      getThetaGradK( Pi, Pj, Pk, defGradk );
      getThetaGradK( Pj, Pi, Pl, defGradl );
      getThetaGradI( Pi, Pj, Pk, Pl, defGradi );
      getThetaGradJ( Pi, Pj, Pk, Pl, defGradj );
      
      // get undeformed vertex positions
      undefShellP->getTo( pi, Pi ); undefShellP->getTo( pj, Pj );
      undefShellP->getTo( pk, Pk ); undefShellP->getTo( pl, Pl );      
            
      // compute Ak + Al, |e|^2 and diff. o dihedral angles
      temp.makeCrossProduct( Pk-Pj, Pi-Pk );
      RealType vol = temp.norm() / 2.;
      temp.makeCrossProduct( Pl-Pi, Pj-Pl );
      vol += temp.norm() / 2.;
      RealType elengthSqr( dotProduct(Pj-Pi,Pj-Pi) ); 
      // note signe here!
      delTheta -= getDihedralAngle( Pi, Pj, Pk, Pl );      

      // derivatives    
      Vec gradk, gradl, gradi, gradj, gradTheta, gradArea;
      RealType factorGradTheta = -2. * elengthSqr / vol;    
      RealType factorGradArea = -2. * delTheta * elengthSqr / (vol*vol);
      RealType factorGradEdgeLengthSqr = 4. * delTheta / vol;       
      
      // d_k
      getThetaGradK( Pi, Pj, Pk, gradTheta );
      getAreaGradK( Pi, Pj, Pk, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradk );
      
      // d_l
      getThetaGradK( Pj, Pi, Pl, gradTheta );
      getAreaGradK( Pj, Pi, Pl, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradl );
      
      // d_i
      getThetaGradI( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pj, Pk, Pi, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradi );
      getAreaGradK( Pl, Pj, Pi, gradArea );
      gradi.addMultiple( gradArea, factorGradArea );
      gradi.addMultiple( Pi-Pj, factorGradEdgeLengthSqr );
      
      // d_j
      getThetaGradJ( Pi, Pj, Pk, Pl, gradTheta );
      getAreaGradK( Pk, Pi, Pj, gradArea );
      getWeightedVectorSum( factorGradTheta, gradTheta, factorGradArea, gradArea, gradj );
      getAreaGradK( Pi, Pl, Pj, gradArea );
      gradj.addMultiple( gradArea, factorGradArea );
      gradj.addMultiple( Pj-Pi, factorGradEdgeLengthSqr );
      
      // k*
      aol::Matrix33<RealType> matrix;
      matrix.makeTensorProduct( gradk, defGradk );
      localToGlobal( Dest, pk, pk, matrix );
      matrix.makeTensorProduct( gradk, defGradl );
      localToGlobal( Dest, pk, pl, matrix );
      matrix.makeTensorProduct( gradk, defGradi );
      localToGlobal( Dest, pk, pi, matrix );
      matrix.makeTensorProduct( gradk, defGradj );
      localToGlobal( Dest, pk, pj, matrix );
      
      // l*
      matrix.makeTensorProduct( gradl, defGradk );
      localToGlobal( Dest, pl, pk, matrix );
      matrix.makeTensorProduct( gradl, defGradl );
      localToGlobal( Dest, pl, pl, matrix );
      matrix.makeTensorProduct( gradl, defGradi );
      localToGlobal( Dest, pl, pi, matrix );
      matrix.makeTensorProduct( gradl, defGradj );
      localToGlobal( Dest, pl, pj, matrix ); 
      
      // i*
      matrix.makeTensorProduct( gradi, defGradk );
      localToGlobal( Dest, pi, pk, matrix );
      matrix.makeTensorProduct( gradi, defGradl );
      localToGlobal( Dest, pi, pl, matrix );
      matrix.makeTensorProduct( gradi, defGradi );
      localToGlobal( Dest, pi, pi, matrix );
      matrix.makeTensorProduct( gradi, defGradj );
      localToGlobal( Dest, pi, pj, matrix );
      
      // j*
      matrix.makeTensorProduct( gradj, defGradk );
      localToGlobal( Dest, pj, pk, matrix );
      matrix.makeTensorProduct( gradj, defGradl );
      localToGlobal( Dest, pj, pl, matrix );
      matrix.makeTensorProduct( gradj, defGradi );
      localToGlobal( Dest, pj, pi, matrix );
      matrix.makeTensorProduct( gradj, defGradj );
      localToGlobal( Dest, pj, pj, matrix );

    }
  }

protected:
  void localToGlobal( BlockMatrixType& BlockOp, int k, int l, const aol::Matrix33<RealType>& localMatrix ) const {    
    if( _firstDerivWRTUndef )   
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( k, l, _factor * localMatrix.get(i,j) );
    }
    else
    {
      for( int i = 0; i < 3; i++ )
        for( int j = 0; j < 3; j++ )
          BlockOp.getReference(i,j).add( l, k, _factor * localMatrix.get(j,i) );
    }
  }
  
};

//!==========================================================================================================
//! DEFORMATION ENERGIES
//!==========================================================================================================
//! \brief Abstract base class to prescribe structure of deformation energies.
//! \author Heeren
//!
//! Deformation energy $E$ is considered as some functional $(S_1, S_2) \mapsto E[S_1, S_2]$,
//! where $S_1$ and $S_2$ are refered to as undeformed and deformed configuration, respectively.
//!
//! Several functions are supposed to be provided by derived classes:
//! - applyEnergy()
//! - applyUndefGradient (), i.e. D_1 E[S_1, S_2]
//! - applyDefGradient (), i.e. D_2 E[S_1, S_2]
//! - assembleAddDefHessian (), i.e. D_2^2 E[S_1, S_2]
//! - assembleAddUndefHessian (), i.e. D_1^2 E[S_1, S_2]
//! - assembleAddMixedHessian (), i.e. D_1 D_2 E[S_1, S_2]
//! where D_i denotes the derivative w.r.t. the ith argument.
template <typename TriMeshConfType >
class DeformationBase {

public:
  typedef aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > BlockMatrixType;
  
protected:
  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   TriMeshType;
  typedef aol::MultiVector<RealType>           VectorType;  
  
  const MeshTopologySaver<TriMeshType>& _topology;

public:
  DeformationBase( const MeshTopologySaver<TriMeshType>& Topology ) : _topology( Topology ){}
  
  virtual ~DeformationBase () {}
    
  virtual void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest ) const = 0;
  
  void applyAddEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest, RealType factor = 1.0 ) const {
    aol::Scalar<RealType> Temp;
    applyEnergy( UndeformedGeom, DeformedGeom, Temp );
    Dest.addMultiple( Temp, factor );
  }
  
  virtual void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const = 0;
  
  void applyAddUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest, RealType factor = 1.0 ) const {
    VectorType Temp( Dest, aol::STRUCT_COPY );
    applyUndefGradient( UndeformedGeom, DeformedGeom, Temp );
    Dest.addMultiple( Temp, factor );
  }
  
  virtual void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const = 0;
  
  void applyAddDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest, RealType factor = 1.0 ) const {
    VectorType Temp( Dest, aol::STRUCT_COPY );
    applyDefGradient( UndeformedGeom, DeformedGeom, Temp );
    Dest.addMultiple( Temp, factor );
  }
   
  virtual void assembleAddDefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const = 0;
  
  virtual void assembleAddUndefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const = 0;
  
  virtual void assembleAddMixedHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTUndef, RealType factor = 1.0 ) const = 0;
  
  void assembleMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTUndef ) const {
    allocateMatrix( Dest );
    Dest.setZero();
    assembleAddMixedHessian( UndeformedGeom, DeformedGeom, Dest, FirstDerivWRTUndef );
  }
  
protected:
  void allocateMatrix( BlockMatrixType& Matrix ) const {
    for ( int i = 0; i < TriMeshConfType::Dim; i++ )
      for ( int j = 0; j < TriMeshConfType::Dim; j++ )
        if( !Matrix.getPointer(i, j) )
          Matrix.allocateMatrix ( i, j, _topology.getNumVertices(), _topology.getNumVertices() );
  }
};

//! \brief Deformation energy representing a membrane energy given by the class MembraneEnergy<>
//! \author Heeren
template <typename TriMeshConfType, typename MaterialLawType = HyperelasticEnergyDensityLog<TriMeshConfType> >
class MembraneDeformation : public DeformationBase<TriMeshConfType> {

protected:
  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   TriMeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef typename DeformationBase<TriMeshConfType>::BlockMatrixType BlockMatrixType;
  
  const aol::ParameterParser& _pparser;

public: 
  MembraneDeformation( const MeshTopologySaver<TriMeshType>& Topology, const aol::ParameterParser& Pparser ) : DeformationBase<TriMeshConfType>( Topology ), _pparser( Pparser) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest ) const {
    MembraneEnergy<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, UndeformedGeom, true ).apply( DeformedGeom, Dest );
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    MembraneGradientUndef<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, DeformedGeom ).apply( UndeformedGeom, Dest );
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    MembraneGradientDef<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, UndeformedGeom ).apply( DeformedGeom, Dest );
  }
  
  void assembleAddDefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    MembraneHessianDef<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, UndeformedGeom, factor ).applyAdd( DeformedGeom, Dest );
  }
  
  void assembleAddUndefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    MembraneHessianUndef<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, DeformedGeom, factor ).applyAdd( UndeformedGeom, Dest );
  }
  
  void assembleAddMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTUndef, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    MembraneHessianMixed<TriMeshConfType, MaterialLawType>( this->_topology, _pparser, UndeformedGeom, true, FirstDerivWRTUndef, factor ).applyAdd( DeformedGeom, Dest );
  }
};

//! \brief Deformation energy representing a bending energy given by the class SimpleBendingEnergy<>
//! \author Heeren
template <typename TriMeshConfType>
class SimpleBendingDeformation : public DeformationBase<TriMeshConfType> {

protected:
  typedef typename TriMeshConfType::RealType   RealType;
  typedef typename TriMeshConfType::InitType   TriMeshType;
  typedef aol::MultiVector<RealType>           VectorType;
  typedef typename DeformationBase<TriMeshConfType>::BlockMatrixType BlockMatrixType;

public:
  SimpleBendingDeformation( const MeshTopologySaver<TriMeshType>& Topology ) : DeformationBase<TriMeshConfType>( Topology ) {}
		            
  SimpleBendingDeformation( const MeshTopologySaver<TriMeshType>& Topology, const aol::ParameterParser& /*Pparser*/ ) : DeformationBase<TriMeshConfType>( Topology ) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, aol::Scalar<RealType>& Dest  ) const {
    SimpleBendingEnergy<TriMeshConfType>( this->_topology, UndeformedGeom, true ).apply( DeformedGeom, Dest );
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    SimpleBendingGradientUndef<TriMeshConfType>( this->_topology, DeformedGeom ).apply( UndeformedGeom, Dest );
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
    SimpleBendingGradientDef<TriMeshConfType>( this->_topology, UndeformedGeom ).apply( DeformedGeom, Dest );
  }
  
  void assembleAddDefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    SimpleBendingHessianDef<TriMeshConfType>( this->_topology, UndeformedGeom, factor ).applyAdd( DeformedGeom, Dest ); 
  }
  
  void assembleAddUndefHessian  ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    SimpleBendingHessianUndef<TriMeshConfType>( this->_topology, DeformedGeom, factor ).applyAdd( UndeformedGeom, Dest ); 
  }
  
  void assembleAddMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, BlockMatrixType& Dest, const bool FirstDerivWRTUndef, RealType factor = 1.0 ) const {
    this->allocateMatrix( Dest );
    bool ActiveShellIsDeformed = true;
    SimpleBendingHessianMixed<TriMeshConfType>( this->_topology, UndeformedGeom, ActiveShellIsDeformed, FirstDerivWRTUndef, factor ).applyAdd( DeformedGeom, Dest );
  }

};
#endif