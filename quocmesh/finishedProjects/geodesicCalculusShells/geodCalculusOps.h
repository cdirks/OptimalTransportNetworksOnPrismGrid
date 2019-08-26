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

#ifndef __GEODCALCULUSOPS_H
#define __GEODCALCULUSOPS_H

#include <aol.h>
#include "deformationEnergies.h"
#include "gradientDescent.h"

using namespace aol;
using namespace qc;
using namespace om;

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes discrete Exp2 operator.
//! \author Heeren
//!
//! Given two shells S_0 and S_1 the Exp2 operator is defined as F[S] = d_2 W[S_0, S_1] + d_1 W[S_1, S],
//! where W = \alpha W_{mem} + \beta W_{bend} is a deformation energy consisting of:
//!   - a MembraneDeformationType W_{mem} (realized as template argument) and a corresponding membrane weight \alpha (to be parsed in the parameter parser)
//!   - a BendingDeformationType W_{bend} (realized as template argument) and a corresponding bending weight \beta (to be parsed in the parameter parser)
//! Here d_i denotes the derivative w.r.t. the ith argument of W[.,.].
//!
//! Note that (S_0, S_1, S_2) is a 3-point-geodesic iff. F[S_2] = 0.

template <typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType >
class Exp2Energy : public aol::Op< aol::MultiVector< typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType TriMeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const aol::ParameterParser& _pparser;
  const MeshTopologySaver<TriMeshType>& _topology;

  const VectorType& _shell0, _shell1;
  RealType  _memWeight, _bendWeight;
  VectorType _constPart;
  
  aol::BitVector *_bdryMask;
  bool _fixBoundary;

public:
    Exp2Energy ( const MeshTopologySaver<TriMeshType>& topology,
                 const aol::ParameterParser& pparser,
                 const VectorType& shell0,
                 const VectorType& shell1 ) :
      _pparser( pparser ),
      _topology( topology ),
      _shell0(shell0),
      _shell1(shell1),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _bdryMask(NULL),
      _fixBoundary(false){
        // initialization
        calcConstPartOfEnergy();
      }
      
  ~Exp2Energy(){ 
    if( _bdryMask ) 
      delete _bdryMask; 
  }
  
  //! if there is a boundary to be fixed
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
  }

  //! The vertex positions of S are given as argument.
  void applyAdd( const VectorType& shell2, VectorType& Dest ) const {

    assert( shell2.numComponents()  == TriMeshConfType::Dim );
    assert( Dest.numComponents()  == TriMeshConfType::Dim );
    
    // add constant partial
    Dest += _constPart;
    
    // add elastic energy variation
    MembraneDeformationType( _topology, _pparser ).applyAddUndefGradient( _shell1, shell2, Dest, _memWeight );

    // add bending gradient
    BendingDeformationType( _topology, _pparser ).applyAddUndefGradient( _shell1, shell2, Dest, _bendWeight );
      
    // fix boundary?
    if( _fixBoundary )
      Dest.setAllMasked( 0., *_bdryMask ); 
  }

protected:
  // pre-cimpute constant part of energy W[S_0, S_1]
  void calcConstPartOfEnergy() {
    _constPart.resize( 3, _topology.getNumVertices() );
    _constPart.setZero();
        
    MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( _shell0, _shell1, _constPart, _memWeight );
    BendingDeformationType(  _topology, _pparser ).applyAddDefGradient( _shell0, _shell1, _constPart, _bendWeight );
  }

};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes derivative of discrete Exp2 operator.
//! \author Heeren
//!
//! Derivative of Exp2 operator F[S] (see documentation of Exp2Energy above), i.e. DF[S] = d_2 d_1 W[S_1,S].
//! Here d_i denotes the derivative w.r.t. the ith argument of W[.,.].
//!
//! As the model is invariant w.r.t. rigid bidy motions (rbm) the matrix representing the derivative is singular.
//! Hence one might either add constraints to prescribe rbm, or regularize the matrix by adding epsilon times the identity matrix.

template < typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType, typename MatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class Exp2Deriv : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, MatrixType> {

protected:
  typedef typename TriMeshConfType::RealType  RealType;
  typedef typename TriMeshConfType::InitType  TriMeshType;
  typedef typename aol::MultiVector<RealType> VectorType;

  const MeshTopologySaver<TriMeshType>& _topology;
  const aol::ParameterParser& _pparser;
  
  const VectorType& _shell1;
  RealType  _bendWeight, _memWeight;
  aol::BitVector *_bdryMask;
  bool _fixBoundary, _regularizeMatrix;
  RealType _gridSize;

public:
  Exp2Deriv( const MeshTopologySaver<TriMeshType>& topology,
             const aol::ParameterParser& pparser,
             const VectorType &shell1,
	     bool regularizeMatrix = false ) :
      _topology( topology ),
      _pparser( pparser ),
      _shell1( shell1 ),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bdryMask(NULL),
      _fixBoundary(false),
      _regularizeMatrix( regularizeMatrix ),
      _gridSize(0.){
	// determine grid size
	if( regularizeMatrix ){
	  TriMeshType topMesh; 
          createMeshFromGeometryAndTopology<TriMeshType>( _topology, _shell1, topMesh );
	  _gridSize = topMesh.H();
	}
      }
      
  ~Exp2Deriv(){ 
    if( _bdryMask ) 
      delete _bdryMask; 
  }  
  
  //! if there is a boundary to be fixed
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
  }

  //! The vertex positions of S are given as argument.
  void applyAdd( const VectorType &shell2, MatrixType &Dest ) const {
    
    // assemble the different components of the Derivative
    MembraneDeformationType( _topology, _pparser ).assembleAddMixedHessian( _shell1, shell2, Dest, true, _memWeight );
    BendingDeformationType( _topology, _pparser ).assembleAddMixedHessian( _shell1, shell2, Dest, true, _bendWeight );

     // add eps*Id to matrix
     if( _regularizeMatrix ){       
       RealType eps( _pparser.getDoubleOrDefault( "regEpsForHessian", 1e-4 ) );
       for ( int j = 0; j < 3; j++ )
         for( int i = 0; i < _topology.getNumVertices(); i++ )
           Dest.getReference( j, j ).add( i, i, eps * _gridSize );
     }     
     
     if( !_fixBoundary )
       return;
     
     // fix boundary ?
     for( int k = 0; k < _bdryMask->size(); k++ ){
      if( (*_bdryMask)[k] ){
        for ( int i = 0; i < 3; i++ ){
	  Dest.getReference( i, i ).setRowColToDiagonal( k );
          Dest.getReference( i, (i+1)%3 ).setRowColToZero( k );
	  Dest.getReference( i, (i+2)%3 ).setRowColToZero( k );
	}
      }
    }

  }
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes 3-point path energy.
//! \author Heeren
//!
//! Given two shells S_0 and S_2 the 3-point path energy is defined as E[S] = W[S_0, S] + W[S, S_2],
//! where W = \alpha W_{mem} + \beta W_{bend} is a deformation energy consisting of:
//!   - a MembraneDeformationType W_{mem} (realized as template argument) and a corresponding membrane weight \alpha (to be parsed in the parameter parser)
//!   - a BendingDeformationType W_{bend} (realized as template argument) and a corresponding bending weight \beta (to be parsed in the parameter parser)

template <typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType >
class Geod3Energy : public aol::Op< aol::MultiVector< typename TriMeshConfType::RealType>, aol::Scalar<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType TriMeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const aol::ParameterParser& _pparser;
  const MeshTopologySaver<TriMeshType>& _topology;

  const VectorType& _shell0, _shell2;
  RealType  _memWeight, _bendWeight;

public:
    Geod3Energy ( const MeshTopologySaver<TriMeshType>& topology,
                  const aol::ParameterParser& pparser,
                  const VectorType& shell0,
                  const VectorType& shell2 ) :
      _pparser( pparser ),
      _topology( topology ),
      _shell0(shell0),
      _shell2(shell2),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bendWeight( pparser.getDouble("bendWeight") ){}

  //! The vertex positions of \f$ \shell_1 \f$ are given as argument.
  void applyAdd( const VectorType& shell1, aol::Scalar<RealType>& Dest ) const {

    assert( shell1.numComponents()  == TriMeshConfType::Dim );    

    // add elastic energy variation
    MembraneDeformationType( _topology, _pparser ).applyAddEnergy ( _shell0, shell1, Dest, _memWeight );
    MembraneDeformationType( _topology, _pparser ).applyAddEnergy ( shell1, _shell2, Dest, _memWeight );
    
    // add bending gradient
    BendingDeformationType( _topology, _pparser ).applyAddEnergy ( _shell0, shell1, Dest, _bendWeight );
    BendingDeformationType( _topology, _pparser ).applyAddEnergy ( shell1, _shell2, Dest, _bendWeight );
  }

};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes first derivative of 3-point path energy.
//! \author Heeren
//!
//! Given two shells S_0 and S_2 the 3-point path energy is defined as E[S] = W[S_0, S] + W[S, S_2] (cf. documentation of class Geod3Energy above).
//! The first derivative w.r.t. S is defined as DE[S] = d_2 W[S_0, S] + d_1 W[S, S_2], where d_i denotes the derivative w.r.t. the ith argument of W[.,.].
template <typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType >
class Geod3Grad : public aol::Op< aol::MultiVector< typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType TriMeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const aol::ParameterParser& _pparser;
  const MeshTopologySaver<TriMeshType>& _topology;

  const VectorType& _shell0, _shell2;
  RealType  _memWeight, _bendWeight;
  aol::BitVector *_bdryMask;
  bool _fixBoundary;

public:
    Geod3Grad ( const MeshTopologySaver<TriMeshType>& topology,
                const aol::ParameterParser& pparser,
                const VectorType& shell0,
                const VectorType& shell2 ) :
      _pparser( pparser ),
      _topology( topology ),
      _shell0(shell0),
      _shell2(shell2),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _bdryMask(NULL),
      _fixBoundary(false){}
      
  ~Geod3Grad(){ 
    if( _bdryMask ) 
      delete _bdryMask; 
  }
  
  //! if there is a boundary to be fixed
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
  }

  //! The vertex positions of \f$ \shell_1 \f$ are given as argument.
  void applyAdd( const VectorType& shell1, VectorType& Dest ) const {

    assert( shell1.numComponents()  == TriMeshConfType::Dim );
    assert( Dest.numComponents()  == TriMeshConfType::Dim );

    // add elastic energy variation
    MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( _shell0, shell1, Dest, _memWeight );
    MembraneDeformationType( _topology, _pparser ).applyAddUndefGradient( shell1, _shell2, Dest, _memWeight );
    
    // add bending gradient
    BendingDeformationType( _topology, _pparser ).applyAddDefGradient( _shell0, shell1, Dest, _bendWeight );
    BendingDeformationType( _topology, _pparser ).applyAddUndefGradient( shell1, _shell2, Dest, _bendWeight );
    
    // fix boundary?
    if( _fixBoundary )
      Dest.setAllMasked( 0., *_bdryMask ); 

  }

};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes second derivative of 3-point path energy.
//! \author Heeren
//!
//! Given two shells S_0 and S_2 the 3-point path energy is defined as E[S] = W[S_0, S] + W[S, S_2] (cf. documentation of class Geod3Energy above).
//! The second derivative w.r.t. S is defined as D^2E[S] = d_2 d_2 W[S_0, S] + d_1 d_1 W[S, S_2], where d_i denotes the derivative w.r.t. the ith argument of W[.,.].
template < typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType, typename MatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class Geod3Hess : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, MatrixType> {

protected:
  typedef typename TriMeshConfType::RealType  RealType;
  typedef typename TriMeshConfType::InitType  TriMeshType;
  typedef typename aol::MultiVector<RealType> VectorType;

  const MeshTopologySaver<TriMeshType>& _topology;
  const aol::ParameterParser& _pparser;
  
  const VectorType& _shell0, _shell2;
  RealType  _bendWeight, _memWeight;
  
  aol::BitVector *_bdryMask;
  bool _fixBoundary, _regularizeMatrix;
  RealType _gridSize;

public:
  Geod3Hess( const MeshTopologySaver<TriMeshType>& topology,
             const aol::ParameterParser& pparser,
             const VectorType& shell0,
             const VectorType& shell2,
	     bool regularizeMatrix = false ) :
      _topology( topology ),
      _pparser( pparser ),
      _shell0(shell0),
      _shell2(shell2),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bdryMask(NULL),
      _fixBoundary(false),
      _regularizeMatrix(regularizeMatrix),
      _gridSize(0.){
	// determine grid size
	if( regularizeMatrix ){
	  TriMeshType topMesh; 
          createMeshFromGeometryAndTopology<TriMeshType>( _topology, _shell0, topMesh );
	  _gridSize = topMesh.H();
	}	
      }
      
  ~Geod3Hess(){ 
    if( _bdryMask ) 
      delete _bdryMask; 
  }
  
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
  }


  //! in Arg steckt das Argument \f$ \shell_1 \f$
  void applyAdd( const VectorType &shell1, MatrixType &Dest ) const {
    
    // assemble the different components of the Derivative
    MembraneDeformationType( _topology, _pparser ).assembleAddDefHessian( _shell0, shell1, Dest, _memWeight );
    MembraneDeformationType( _topology, _pparser ).assembleAddUndefHessian( shell1, _shell2, Dest, _memWeight );

    BendingDeformationType( _topology, _pparser ).assembleAddDefHessian( _shell0, shell1, Dest, _bendWeight );
    BendingDeformationType( _topology, _pparser ).assembleAddUndefHessian( shell1, _shell2, Dest, _bendWeight );

     // add eps*Id to matrix
     if( _regularizeMatrix ){
       RealType eps( _pparser.getDoubleOrDefault( "regEpsForHessian", 1e-4 ) );
       for ( int j = 0; j < 3; j++ )
         for( int i = 0; i < _topology.getNumVertices(); i++ )
           Dest.getReference( j, j ).add( i, i, eps * _gridSize );
     }
      
     if( !_fixBoundary )
       return;
      
     // fix boundary?
     for( int k = 0; k < _bdryMask->size(); k++ ){
       if( (*_bdryMask)[k] ){
         for ( int i = 0; i < 3; i++ ){
	   Dest.getReference( i, i ).setRowColToDiagonal( k );
           Dest.getReference( i, (i+1)%3 ).setRowColToZero( k );
	   Dest.getReference( i, (i+2)%3 ).setRowColToZero( k );
	 }
       }
     }
  }
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes (K+1)-point path energy.
//! \author Heeren
//!
//! Given a start shell S_0 = S_A and an end shell S_K = S_B, the (K+1)-point path energy is defined as E[S_1, \ldots, S_{K-1} ] = \sum_{k=1}^K W[S_{k-1}, S_{k}]
//! where W = \alpha W_{mem} + \beta W_{bend} is a deformation energy consisting of:
//!   - a MembraneDeformationType W_{mem} (realized as template argument) and a corresponding membrane weight \alpha (to be parsed in the parameter parser)
//!   - a BendingDeformationType W_{bend} (realized as template argument) and a corresponding bending weight \beta (to be parsed in the parameter parser)
//! The number of shells (K+1) has to be parsed in the constructor.
template <typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicEnergy : public aol::Op< aol::MultiVector< typename TriMeshConfType::RealType>, aol::Scalar<typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType TriMeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const aol::ParameterParser& _pparser;
  const MeshTopologySaver<TriMeshType>& _topology;

  const VectorType& _start, _end;
  RealType  _memWeight, _bendWeight;
  const int _numOfShells;

public:
    GeodesicEnergy ( const MeshTopologySaver<TriMeshType>& topology,
                  const aol::ParameterParser& pparser,
                  const VectorType& start,
                  const VectorType& end,
		  const int numOfShells ) :
      _pparser( pparser ),
      _topology( topology ),
      _start( start ),
      _end( end ),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _numOfShells( numOfShells ){}

  //! The vertex positions of \f$ \shell_1, \ldots, \shell_{K-1} \f$ are given as arguments.
  void applyAdd( const VectorType& Arg, aol::Scalar<RealType>& Dest ) const {
    
    const int numOfFreeShells = _numOfShells-2;
    const int dim = TriMeshConfType::Dim;
    assert( Arg.numComponents()  ==  dim * numOfFreeShells );      
    
    // bring varaiables into more convenient form
    aol::RandomAccessContainer< VectorType > arg( _numOfShells );
    arg[0].appendReference( _start );
    for( int i = 0; i < numOfFreeShells; i++ )
      for( int j = 0; j < dim; j++ )
	arg[i+1].appendReference( Arg[i*dim+j] );
    arg[_numOfShells-1].appendReference( _end );

    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int shellIdx = 1; shellIdx < _numOfShells; shellIdx++ ){

      aol::Scalar<RealType> aux;
      aux.setZero();

      // add energy terms
      MembraneDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _memWeight );
      BendingDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _bendWeight );

#ifdef _OPENMP
#pragma omp critical (GenericGeodesicEnergy_applyAdd)
#endif
      Dest += aux;
    }

  }
  
  void evaluateSingleEnergies( const VectorType& Arg, aol::Vector<RealType>& Dest, bool print = false ) const {
    const int numOfFreeShells = _numOfShells-2;
    const int dim = TriMeshConfType::Dim;
    assert( Arg.numComponents()  ==  dim * numOfFreeShells );     
    
    Dest.reallocate( numOfFreeShells + 1 );
    
    // bring varaiables into more convenient form
    aol::RandomAccessContainer< VectorType > arg( _numOfShells );
    arg[0].appendReference( _start );
    for( int i = 0; i < numOfFreeShells; i++ )
      for( int j = 0; j < dim; j++ )
	arg[i+1].appendReference( Arg[i*dim+j] );
    arg[_numOfShells-1].appendReference( _end );

    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int shellIdx = 1; shellIdx < _numOfShells; shellIdx++ ){

      aol::Scalar<RealType> aux;
      aux.setZero();

      // add energy terms
      MembraneDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _memWeight );
      BendingDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _bendWeight );
      
      Dest[shellIdx-1] = aux[0];
    }
    
    if( print ){
      RealType pathenergy = 0.;
      for ( int shellIdx = 0; shellIdx < _numOfShells - 1; shellIdx++ ){
	pathenergy += Dest[shellIdx];
	cerr << shellIdx << "th deformation energy is " << Dest[shellIdx] << endl;
      }
      cerr << "Path energy is " << (_numOfShells - 1) * pathenergy << endl;
      cerr << "Path length is " << sqrt( (_numOfShells - 1) * pathenergy ) << endl;
    }
  }
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes first derivative of (K+1)-point path energy.
//! \author Heeren
//!
//! Given a start shell S_0 = S_A and an end shell S_K = S_B, the (K+1)-point path energy is defined as E[S_1, \ldots, S_{K-1} ] = \sum_{k=1}^K W[S_{k-1}, S_{k}] (cf. documentation of class GeodesicEnergy above).
//! The first derivative w.r.t. S_i is defined as D_i E = d_2 W[S_{i-1}, S_i] + d_1 W[S_i, S_{i+1}], where d_i denotes the derivative w.r.t. the ith argument of W[.,.].
template <typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicGradient : public aol::Op< aol::MultiVector< typename TriMeshConfType::RealType> > {

protected:
  typedef typename TriMeshConfType::RealType RealType;
  typedef typename TriMeshConfType::InitType TriMeshType;
  typedef aol::MultiVector<RealType>         VectorType;

  const aol::ParameterParser& _pparser;
  const MeshTopologySaver<TriMeshType>& _topology;

  const VectorType& _start, _end;
  RealType  _memWeight, _bendWeight;
  const int _numOfShells;
  aol::BitVector *_bdryMask;
  bool _fixBoundary;

public:
    GeodesicGradient ( const MeshTopologySaver<TriMeshType>& topology,
                const aol::ParameterParser& pparser,
                const VectorType& start,
                const VectorType& end,
		const int numOfShells  ) :
      _pparser( pparser ),
      _topology( topology ),
      _start( start ),
      _end( end ),
      _memWeight( pparser.getDouble("tangWeight") ),
      _bendWeight( pparser.getDouble("bendWeight") ),
      _numOfShells( numOfShells ),
      _bdryMask(NULL),
      _fixBoundary( false ){}
      
     
  ~GeodesicGradient(){ 
    if( _fixBoundary ) 
      delete _bdryMask; 
  }
  
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
    if( _bdryMask->size() != _topology.getNumVertices() )
      throw aol::Exception ( "GeodesicGradient::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
  }

  //! The vertex positions of \f$ \shell_1, \ldots, \shell_{K-1} \f$ are given as arguments.
  void applyAdd( const VectorType& Arg, VectorType& Dest ) const {

    const int numOfFreeShells = _numOfShells-2;
    const int dim = TriMeshConfType::Dim;
    assert( Arg.numComponents()  ==  dim * numOfFreeShells );  
    assert( Dest.numComponents()  == Arg.numComponents() );
    
    // bring varaiables into more convenient form
    aol::RandomAccessContainer< VectorType > arg( _numOfShells ), dest( numOfFreeShells );
    arg[0].appendReference( _start );
    for( int i = 0; i < numOfFreeShells; i++ )
      for( int j = 0; j < dim; j++ ){
	arg[i+1].appendReference( Arg[i*dim+j] );
	dest[i].appendReference( Dest[i*dim+j] );
      }
    arg[_numOfShells-1].appendReference( _end );


#ifdef _OPENMP
#pragma omp parallel for
#endif
    // deformierte shells
    for ( int shellIdx = 0; shellIdx < _numOfShells-2; shellIdx++ ) {

      // add variations with respect to deformed geometry
      MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx], _memWeight );
      BendingDeformationType( _topology, _pparser ).applyAddDefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx], _bendWeight );

    }
    
#ifdef _OPENMP
#pragma omp parallel for
#endif    
    for ( int shellIdx = 1; shellIdx < _numOfShells-1; shellIdx++ ){
      
      // add variations with respect to undeformed geometry
      MembraneDeformationType( _topology, _pparser ).applyAddUndefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx-1], _memWeight );
      BendingDeformationType( _topology, _pparser ).applyAddUndefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx-1], _bendWeight );
      
    }
    
    if( !_fixBoundary )
      return;
    
    // fix boundary ?
    for ( int shellIdx = 0; shellIdx < _numOfShells-2; shellIdx++ )
      dest[shellIdx].setAllMasked( 0., *_bdryMask );

  }
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class that realizes second derivative of (K+1)-point path energy.
//! \author Heeren
//!
//! Given a start shell S_0 = S_A and an end shell S_K = S_B, the (K+1)-point path energy is defined as E[S_1, \ldots, S_{K-1} ] = \sum_{k=1}^K W[S_{k-1}, S_{k}] (cf. documentation of class GeodesicEnergy above).
//! The second derivative is defined as D_{i,i+1} E =  d_2 d_1 W[S_i, S_{i+1}] or D_{i,i} E =  d_1 d_1 W[S_i, S_{i+1}] + d_2 d_2 W[S_{i-1}, S_i], where d_i denotes the derivative w.r.t. the ith argument of W[.,.].
//! Note that D_{i,j}E = 0 if |i-j| > 1.
template < typename TriMeshConfType, typename MembraneDeformationType, typename BendingDeformationType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename TriMeshConfType::RealType> > >
class GeodesicHessian : public aol::Op<aol::MultiVector<typename TriMeshConfType::RealType>, BlockMatrixType> {

protected:
  typedef typename TriMeshConfType::RealType  RealType;
  typedef typename TriMeshConfType::InitType  TriMeshType;
  typedef typename aol::MultiVector<RealType> VectorType;

  const MeshTopologySaver<TriMeshType>& _topology;
  const aol::ParameterParser& _pparser;
  
  const VectorType& _start, _end;
  RealType  _bendWeight, _memWeight;
  const int _numOfShells;        
  aol::BitVector *_bdryMask;
  bool _fixBoundary, _regularizeMatrix;
  

public:
  GeodesicHessian( const MeshTopologySaver<TriMeshType>& topology,
                   const aol::ParameterParser& pparser,
                   const VectorType& start,
                   const VectorType& end,
		   const int numOfShells   ) :
      _topology( topology ),
      _pparser( pparser ),
      _start( start ),
      _end( end ),      
      _bendWeight( pparser.getDouble("bendWeight") ),
      _memWeight( pparser.getDouble("tangWeight") ),
      _numOfShells( numOfShells ),
      _bdryMask(NULL),
      _fixBoundary( false ),
      _regularizeMatrix( false ){}
     
  ~GeodesicHessian(){ 
    if( _fixBoundary ) 
      delete _bdryMask; 
  }
  
  void setBoundaryMask( const aol::BitVector& mask ){
    _fixBoundary = true; 
    if ( _bdryMask )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( mask );
    if( _bdryMask->size() != _topology.getNumVertices() )
      throw aol::Exception ( "GeodesicHessian::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
  }
  
  //! The vertex positions of \f$ \shell_1, \ldots, \shell_{K-1} \f$ are given as arguments.
  void applyAdd( const VectorType &Arg, BlockMatrixType &Dest ) const {
    
    const int numOfFreeShells = _numOfShells-2;
    const int dim = TriMeshConfType::Dim;
    assert( Arg.numComponents()  ==  dim * numOfFreeShells );  
    assert( Dest.getNumRows() == Dest.getNumCols());
    
    // bring varaiables into more convenient form
    aol::RandomAccessContainer< VectorType > arg( _numOfShells );
    arg[0].appendReference( _start );
    for( int i = 0; i < numOfFreeShells; i++ )
      for( int j = 0; j < dim; j++ )
	arg[i+1].appendReference( Arg[i*dim+j] );
    arg[_numOfShells-1].appendReference( _end );
    
    // deformierte shells
#ifdef _OPENMP
#pragma omp parallel for
#endif    
    for ( int shellIdx = 0; shellIdx < _numOfShells-2; shellIdx++ ) {
      
      // get reference on subblock of Dest
      std::vector< std::vector< std::vector< int > > > blockPositions( dim );       
      getBlockPositions( blockPositions, shellIdx );
      BlockMatrixType dest(dim, dim);
      Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

      // add second derivatives w.r.t. the deformed geometries
      MembraneDeformationType( _topology, _pparser ).assembleAddDefHessian( arg[shellIdx], arg[shellIdx+1], dest, _memWeight );
      BendingDeformationType( _topology, _pparser ).assembleAddDefHessian( arg[shellIdx], arg[shellIdx+1], dest, _bendWeight );

    }
    
#ifdef _OPENMP
#pragma omp parallel for
#endif    
    for ( int shellIdx = 1; shellIdx < _numOfShells-1; shellIdx++ ){
      
      // get reference on subblock of Dest
      std::vector< std::vector< std::vector< int > > > blockPositions( dim );       
      getBlockPositions( blockPositions, shellIdx-1 );
      BlockMatrixType dest(dim, dim);
      Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  
      
      // add second derivatives w.r.t. the undeformed geometries
      MembraneDeformationType( _topology, _pparser ).assembleAddUndefHessian( arg[shellIdx], arg[shellIdx+1],  dest, _memWeight );
      BendingDeformationType( _topology, _pparser ).assembleAddUndefHessian( arg[shellIdx], arg[shellIdx+1], dest, _bendWeight );
      
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif              
    for ( int shellIdx = 1; shellIdx < _numOfShells-2; shellIdx++ ){
      
      // get reference on subblock of Dest
      std::vector< std::vector< std::vector< int > > > blockPositions( dim ); 
      BlockMatrixType dest(dim, dim);
      
      getBlockPositions( blockPositions, shellIdx-1, 1 );      
      Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

      // assemble the different components of the Derivative
      MembraneDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, true, _memWeight );
      BendingDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, true, _bendWeight );           
     
      //! TODO use symmetry instead of computing again!
      getBlockPositions( blockPositions, shellIdx-1, -1 );      
      Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

      // assemble the different components of the Derivative
      MembraneDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, false, _memWeight );
      BendingDeformationType(  _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, false, _bendWeight ); 
    }
   
    
    if( !_fixBoundary )
     return;
      
     // fix boundary?
     for ( int shellIdx = 0; shellIdx < _numOfShells-2; shellIdx++ ) {
       std::vector< std::vector< std::vector< int > > > blockPositions( dim );   
       BlockMatrixType dest(dim, dim);
       
       // diagonal blocks
       getBlockPositions( blockPositions, shellIdx );       
       Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );
       applyBdryMaskToSubmatrix( dest );
       
       if( shellIdx == 0 )
	 continue;
       
       // off diagonal blocks
       getBlockPositions( blockPositions, shellIdx-1, 1 );       
       Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );
       applyBdryMaskToSubmatrix( dest );
       
       getBlockPositions( blockPositions, shellIdx-1, -1 );       
       Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );
       applyBdryMaskToSubmatrix( dest ); 
     }     

  }
  
protected:
  void getBlockPositions( std::vector< std::vector< std::vector< int > > >& blockPositions, int k, int off = 0 ) const {
    const int dim = TriMeshConfType::Dim;
    int colShift = off > 0 ? off : 0;
    int rowShift = off < 0 ? off : 0;
    for ( int i = 0; i < dim; ++i ){
      blockPositions[i].resize ( dim );
      for ( int j = 0; j < dim; ++j ) {
        blockPositions[i][j].resize ( 2 );
        blockPositions[i][j][0] = (k - rowShift) * dim + i;
        blockPositions[i][j][1] = (k + colShift) * dim + j;
      }
    }
  }
  
  void applyBdryMaskToSubmatrix( BlockMatrixType& Submatrix ) const {
    for( int k = 0; k < _bdryMask->size(); k++ )
      if( (*_bdryMask)[k] )
         for ( int i = 0; i < 3; i++ ){
	     Submatrix.getReference( i, i ).setRowColToDiagonal( k );
             Submatrix.getReference( i, (i+1)%3 ).setRowColToZero( k );
	     Submatrix.getReference( i, (i+2)%3 ).setRowColToZero( k );
	 }
 
  }
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class to compute discrete geodesics.
//! \author Heeren
//!
//! Given a start shell S_0 = S_A and an end shell S_K = S_B, the (K+1)-point path energy (cf. documentation of GeodesicEnergy above) is minimized.
//! This can be done:
//!   - by a straight forward optimization via computeGeodesic()
//!   - using progressive methods as a hierachical scheme in SPACE via computeGeodesicMultiLevel()
//!   - using a red-black-approach as a hierachical scheme in TIME via computeGeodesicRedBlack()
template<typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType>
class GeodesicOp{
  
   typedef typename ConfiguratorType::RealType RealType;
   typedef aol::MultiVector<RealType> VectorType;
   typedef typename ConfiguratorType::InitType MeshType;
   typedef aol::RandomAccessContainer< MeshType > MeshContainer;
    
   typedef SparseMatrix<RealType> SubMatrixType;
   typedef SparseBlockMatrix<SubMatrixType> BlockMatrixType;
   
   const aol::ParameterParser& _pparser;
   const MeshTopologySaver<MeshType>* _topology;
   mutable MeshType _topolMesh;    
   const VectorType *_startShape, *_endShape;
   mutable int _numOfShells, _numOfFreeShells;   
   int _levels, _alternatingSteps;
   RealType _theta;
   int _solverType;
   int _solvingMethod;
   int _numOfDescentSteps;
   RealType _stopCriterion;
   char _destDirectory[1024], _savenameStem[1024];
   bool _fixBoundary;
   mutable aol::BitVector *_bdryMask; 
   bool _regularizeMatrix, _lagrangeSetup;
   bool _quiet, _deleteTopology;
   mutable bool _saving;
   
   const static int _NumOfSingleConstraints = 6;
   
public:
  GeodesicOp( const aol::ParameterParser& pparser ) 
    : _pparser( pparser), 
    _topolMesh( _pparser.getString("endShape") ), 
    _numOfShells( pparser.getInt("numOfShells") ),
    _numOfFreeShells( _numOfShells - 2 ),
    _levels( _pparser.getIntOrDefault("multilevel", 1) ),
    _alternatingSteps( _pparser.getIntOrDefault("redBlackSteps", 0) ),
    _theta( _pparser.getDoubleOrDefault("theta", 0.) ),
    _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
    _solvingMethod( _pparser.getIntOrDefault( "solverMethod", GRADIENTDESCENT) ),
    _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
    _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
    _fixBoundary( _pparser.checkAndGetBool("fixBoundary") || _pparser.hasVariable("bdryMask") ),
    _bdryMask( NULL ),
    _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
    _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
    _quiet( false ),
    _deleteTopology( true ),
    _saving( true ){
      
      // dymanic allocating
      _topology = new MeshTopologySaver<MeshType>( _topolMesh );      

      VectorType temp( 3, _topolMesh.getNumVertices() );
      _topolMesh.toVector( temp);
      _endShape = new VectorType( temp );
      
      _topolMesh.loadFromPLY( _pparser.getString("startShape") );
      _topolMesh.toVector( temp );
      _startShape = new VectorType( temp );
      
      pparser.getString( "destDirectory", _destDirectory );
      pparser.getString( "savenameStem", _savenameStem );
      
      // boundary mask? 
      if( _fixBoundary ){
        if( pparser.hasVariable("bdryMask") ){
	  cerr << "Load boundary mask..." << endl;	
	  _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
        }
        else{
	  cerr << "Set boundary mask..." << endl;	
	   _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	   _topolMesh.fillBoundaryMask( *_bdryMask );
        }
      } 
      
      cerr << "Geodesic computation between: " << endl << _pparser.getString("startShape")<< endl << _pparser.getString("endShape") << endl << endl;
      
      printSettings();
    }
  
  GeodesicOp( const aol::ParameterParser& pparser, const MeshTopologySaver<MeshType>& Topology, bool Quiet = true ) 
    : _pparser( pparser), 
    _topology( &Topology ),
    _topolMesh( Topology.getGrid() ),
    _startShape( NULL ),
    _endShape( NULL ),
    _numOfShells( -1 ),
    _levels( _pparser.getIntOrDefault("multilevel", 1) ),
    _alternatingSteps( _pparser.getIntOrDefault("redBlackSteps", 1) ),
    _theta( _pparser.getDoubleOrDefault("theta", 0.) ),
    _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
    _solvingMethod( _pparser.getIntOrDefault( "solverMethod", GRADIENTDESCENT) ),
    _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
    _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
    _fixBoundary( _pparser.checkAndGetBool("fixBoundary") || _pparser.hasVariable("bdryMask") ),
    _bdryMask( NULL ),
    _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
    _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
    _quiet( Quiet ),
    _deleteTopology( false ),
    _saving( false ){ 
      // boundary mask? 
      if( _fixBoundary ){
        if( pparser.hasVariable("bdryMask") ){
	  if( !_quiet ) cerr << "Load boundary mask..." << endl;	
	  _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
        }
        else{
	  if( !_quiet ) cerr << "Set boundary mask..." << endl;	
	   _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	   _topolMesh.fillBoundaryMask( *_bdryMask );
        }
      } 
    }
  
  ~GeodesicOp() {
    if( _startShape )
      delete _startShape;
    if( _endShape )
      delete _endShape;
    if( _deleteTopology && _topology )
      delete _topology;
    if( _fixBoundary ) 
      delete _bdryMask;
  }
  
  void setBoundaryMask( const BitVector& Mask ){
    if( !_quiet ) cerr << "Set boundary mask..." << endl;	    
    if( _fixBoundary )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( Mask );
    _fixBoundary = true;
  }
  
  void printSettings() const {
    cerr << endl << "===========================" << endl;
    cerr << "Solving flags:" << endl;	
    cerr << " - max number of steps = " << _numOfDescentSteps << endl;
    cerr << " - stoping criterion = " << _stopCriterion << endl;
    if( _lagrangeSetup )
      cerr << " - using Lagrange setup" << endl;
    if( _regularizeMatrix )
      cerr << " - regularize matrix by adding eps*Id" << endl;
    cerr << "===========================" << endl << endl;
  }

  void setQuietMode( bool Quiet ){
    _quiet = Quiet;
  }
  
  //! Computes discrete geodesic according to optimization mode specified in parameter parser
  void execute( ) const { 
    execute( _pparser.getIntOrDefault("geodesicMode", 0) );
  }
  
  //! Computes discrete geodesic according to optimization mode:
  //!   (0) computeGeodesic()
  //!   (1) computeGeodesicMultiLevel()
  //!   (2) computeGeodesicRedBlack()
  //!   (3) computeGeodesicCascadic()  
  void execute( int optimizationMode ) const {
    // check setup
    if( !_startShape ) throw aol::Exception ( "GeodesicOp::execute(): start shape not set!", __FILE__, __LINE__ );
    if( !_endShape ) throw aol::Exception ( "GeodesicOp::execute(): end shape not set!", __FILE__, __LINE__ );    
    if( _topology->getNumVertices() == 0 ) throw aol::Exception ( "GeodesicOp::execute(): topology not set!", __FILE__, __LINE__ );
    if( _topolMesh.getNumVertices() == 0 ) throw aol::Exception ( "GeodesicOp::execute(): mesh not set!", __FILE__, __LINE__ );
    
    if( !_quiet ){
      cerr << "Start geodesic computation";
      if( optimizationMode > 0 ) 
	cerr << " with " << ( (optimizationMode == 1)  ? "multilevel" : ( (optimizationMode == 2) ? "red-black" : "cascadic" ) ) << " refinement..." << endl;
      else
	cerr << ".\n";
      //cerr << "Length of geodesic = " << _numOfShells << endl;
      cerr << "#nodes = " << _topology->getNumVertices() << " and #faces = " << _topology->getNumFaces() << endl << endl;
    }
    
    VectorType Path;    
    switch( optimizationMode ){
      case 0: computeGeodesic( *_startShape, *_endShape, _numOfShells, Path ); break;
      case 1: computeGeodesicMultiLevel( *_startShape, *_endShape, _numOfShells, Path ); break;
      case 2: computeGeodesicRedBlack( *_startShape, *_endShape, _numOfShells, Path ); break;
      case 3: computeGeodesicCascadic( *_startShape, *_endShape, _pparser.getIntOrDefault( "cascadicLevels", 1 ), Path ); break;      
      default: throw aol::Exception ( "GeodesicOp::execute(): unknown optimization mode!", __FILE__, __LINE__ );
    }    
  }
  
  //! Computes discrete geodesic via straight forward optmization of path energy E = K \sum_k W[S_{k-1}, S_k]; 
  //! Discrete geodesic is stored in Path, minimal energy value is returned.
  //! If initialize = false, Path contains already contains a suitable initialization when calling.
  RealType computeGeodesic( const VectorType& startShape, const VectorType& endShape, int numOfShells, VectorType& Path, bool initialize = true ) const {    
         
    if( !_quiet ) cerr <<"Start geodesic optimization with Newton method..." << endl;
    
    // update boundary mask if necessary
    if( _fixBoundary )
      checkAndUpdateBdryMask( _topology->getGrid() );
    
    _numOfShells = numOfShells;
    _numOfFreeShells = _numOfShells - 2;    
    
    if( Path.numComponents() != ConfiguratorType::Dim * _numOfFreeShells )
      Path.reallocate( ConfiguratorType::Dim * _numOfFreeShells, _topology->getNumVertices() );
    
    // initialization
    if( initialize ){
      if( !_quiet ) cerr << "Start initialization... "  << endl;
      initializePath( startShape, endShape, Path );
    }
      
    aol::Scalar<RealType> energy;
    GeodesicEnergy<ConfiguratorType, MembraneDeformationType, BendingDeformationType > E( *_topology, _pparser, startShape, endShape, _numOfShells  );
    GeodesicGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, startShape, endShape, _numOfShells  );
    GeodesicHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE2( *_topology, _pparser, startShape, endShape, _numOfShells  );
    
    // fix boundary?
    if( _fixBoundary ){
      dE.setBoundaryMask( *_bdryMask );
      dE2.setBoundaryMask( *_bdryMask );
    }

    if( !_quiet && _lagrangeSetup ) cerr << "Using Lagrange setup." << endl;
      
    // intial energy
    E.apply( Path, energy );
    if( !_quiet ) cerr << "Initial energy = " <<  energy[0] * (_numOfFreeShells+1) << endl;       
   
    // rigid body motion handler
    int NumOfSingleConstraints = _lagrangeSetup ? _NumOfSingleConstraints : 0;
    int NumOfConstraints = _numOfFreeShells * NumOfSingleConstraints;
    typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
    RBMHandlerType CHandler( *_topology, startShape, NumOfSingleConstraints, _numOfFreeShells );

    GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( dE, CHandler );
    GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( dE2, CHandler );
   
    // append Lagrange multipliers
    VectorType Destination;           
    VectorType lagrangeMult( NumOfConstraints, 1 );
    Destination.appendReference( Path );
    if( _lagrangeSetup ) Destination.appendReference( lagrangeMult );
    VectorType Argument( Destination );
    
    // initialize solver and solve
    RealType stopCriterion = max( _stopCriterion, 1e-8 * sqrt( _topology->getNumVertices() * _numOfFreeShells) );
    if( !_quiet ) cerr <<"Start solving with stopping criterion " << stopCriterion << " ... " << endl;
    NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( *_topology, dL, dL2, ConfiguratorType::Dim * _numOfFreeShells, NumOfConstraints, _solverType, _numOfDescentSteps, stopCriterion );
    zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
    zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
    zeroFinder.setQuietMode( _quiet );      
    zeroFinder.apply( Argument, Destination );   
     
    if( !_quiet && _lagrangeSetup ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Path ) << endl;
   
    // save solution    
    if( _saving ){
      if( !_quiet ) cerr <<"Save final solution... " << endl;
      
      VectorType fullPath;
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
        fullPath.appendReference( startShape[j] );
      for( int i = 0; i < _numOfFreeShells; i++ )
        for( int j = 0; j < ConfiguratorType::Dim; j++ )
          fullPath.appendReference( Path[ ConfiguratorType::Dim*i + j ] );
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
        fullPath.appendReference( endShape[j] );
      
      ostringstream savename;
      savename << _destDirectory << _savenameStem;       
      saveResults( fullPath, savename.str() );
    }
              
    // final energy
    E.apply( Path, energy );
    if( !_quiet ) cerr << "Final energy = " <<  energy[0] * (_numOfFreeShells+1) << endl;       
    return energy[0] * (_numOfFreeShells+1);
  }
  
  //! Computes discrete geodesic via multilevel optimization of path energy E = K \sum_k W[S_{k-1}, S_k]; 
  //! Here multilevel refers to a spatial hierarchy realized by progressive meshes.
  RealType computeGeodesicMultiLevel( const VectorType& startShape, const VectorType& endShape, int numOfShells, VectorType& Path ) const {    
         
    if( !_topology )
      throw aol::Exception ( "GeodesicOp::executeMultiLevel(): Topology not set!", __FILE__, __LINE__ );    
    
    if( (_levels > 1) && _fixBoundary && !_quiet )
      cerr << "WARNING: when doing multilevel and fixing boundaries loading boundary masks is not possible! " << endl << endl;
    
    _numOfShells = numOfShells;
    _numOfFreeShells = numOfShells - 2;   
    aol::Scalar<RealType> finalEnergy;
    
    int NumOfSingleConstraints = _lagrangeSetup ? _NumOfSingleConstraints : 0;
    int NumOfConstraints = _numOfFreeShells * NumOfSingleConstraints;
    
    // multilevel preparations    
    MeshContainer startShapes, endShapes;        
    _topolMesh.fromVector( endShape );
    endShapes.pushBack( _topolMesh );
    _topolMesh.fromVector( startShape );
    startShapes.pushBack( _topolMesh );
    
    // compute decimation
    aol::RandomAccessContainer< ProlongInfo<RealType> > prolongationInfos( 2 * (_levels-1) );
    if( _levels > 1 )
      computeDecimation( startShapes, endShapes, prolongationInfos );    
    
    // initialize arguments
    if( !_quiet ) cerr <<"Start initialization..." << endl;
    VectorType startShapeOnCurrentLevel( ConfiguratorType::Dim, startShapes[_levels-1].getNumVertices() );
    startShapes[_levels-1].toVector( startShapeOnCurrentLevel );
    VectorType endShapeOnCurrentLevel( ConfiguratorType::Dim, startShapes[_levels-1].getNumVertices() );
    endShapes[_levels-1].toVector( endShapeOnCurrentLevel );
    
    VectorType InitialGuessOnCurrentLevel( ConfiguratorType::Dim * _numOfFreeShells, startShapes[_levels-1].getNumVertices() );
    Path.reallocate( InitialGuessOnCurrentLevel );
    
    // initialize with start shape on coarsest level
    for( int i = 0; i < _numOfFreeShells; i++ )
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
        InitialGuessOnCurrentLevel[ ConfiguratorType::Dim*i + j ] = startShapeOnCurrentLevel[ j ];   
    // initialize linearly  
    if( !_pparser.checkAndGetBool("initializeTrivially") )
      initializePathLinearly( startShapeOnCurrentLevel, endShapeOnCurrentLevel, InitialGuessOnCurrentLevel );    
    
    // start multilevel
    for( int level = 1; level < _levels + 1; level++ ){
      
      if( !_quiet ) cerr << "=============================================================" << endl;
            
      // topology on this level
      int currentLevel = _levels-level;
      MeshTopologySaver<MeshType> topology( startShapes[currentLevel] );
      if( !_quiet ) cerr << "Start multilevel method with " << numOfShells << " shapes on level " << _levels-currentLevel << " of " << _levels  << " with #nodes = " << topology.getNumVertices() << endl;      
      
      // save initialization      
      if( _saving && _pparser.checkAndGetBool("saveInitialization") ){
	if( !_quiet ) cerr <<"Save initialization... " << endl;
        ostringstream savename;
	if( _levels > 1 )
          savename << _destDirectory << _savenameStem << "_level" << level << "_initial";    
	else
	  savename << _destDirectory << _savenameStem << "_initial";  
        saveResults( startShapes[currentLevel], endShapes[currentLevel], InitialGuessOnCurrentLevel, savename.str() );
      }      
      
      // geodesic energy and derivatives
      GeodesicEnergy<ConfiguratorType, MembraneDeformationType, BendingDeformationType >     E( topology, _pparser, startShapeOnCurrentLevel, endShapeOnCurrentLevel, numOfShells  );
      GeodesicGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( topology, _pparser, startShapeOnCurrentLevel, endShapeOnCurrentLevel, numOfShells  );
      GeodesicHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE2( topology, _pparser, startShapeOnCurrentLevel, endShapeOnCurrentLevel, numOfShells  );
      
      // fix boundary?
      if( _fixBoundary ){
	if( !_quiet ) cerr <<"Update boundary mask if necessary " << endl;
	checkAndUpdateBdryMask( startShapes[currentLevel] );
	dE.setBoundaryMask( *_bdryMask );
	dE2.setBoundaryMask( *_bdryMask );
      }
    
      bool solved = false;      
      
      // initialize with cascadic approach on coarsest level
      if( (level == 1) && (numOfShells > 3) && (numOfShells%2 == 1) ){
	if( !_quiet ) cerr <<"Start solving with red-black approach on coarsest level... " << endl;
	bool oldSaving = _saving;
	_saving = false;
	GeodesicOp<ConfiguratorType,MembraneDeformationType,BendingDeformationType>( _pparser, topology, _quiet ).computeGeodesicRedBlack( startShapeOnCurrentLevel, endShapeOnCurrentLevel, numOfShells, Path, true ); 	
	solved = true;
	_saving = oldSaving;
      }
                  
      // gradient descent
      if( !solved && (_solvingMethod == GRADIENTDESCENT) ){
	if( !_quiet ) cerr <<"Start solving with gradient descent... " << endl;
	//throw aol::Exception ( "GeodesicOp::executeMultiLevel(): not implemented!", __FILE__, __LINE__ ); 
	typedef aol::GridlessGradientDescent<RealType,VectorType,VectorType > GDType;
	GDType descent( E, dE, _numOfDescentSteps, 1., _stopCriterion);
        descent.setTimestepController ( GDType::POWELL_WOLFE );
        descent.setConfigurationFlags ( GDType::USE_NONLINEAR_CG|GDType::LOG_GRADIENT_NORM_AT_OLD_POSITION );
        descent.apply( InitialGuessOnCurrentLevel, Path );
	solved = true;
      }
      
      // quasi Newton
      if( !solved && (_solvingMethod == QUASINEWTON) ){
	if( !_quiet ) cerr <<"Start solving with Quasi Newton method... " << endl;
        int restart = _pparser.getIntOrDefault("restart", 1);
        aol::QuasiNewtonIteration< RealType, VectorType, VectorType  > quasiNewtonDescent( E, dE, _numOfDescentSteps, _stopCriterion, restart, false );
        quasiNewtonDescent.setTauMin( .1 ); // experience shows: either tau=1 or tau=2 or tau=0
	if( _pparser.hasVariable("sigma") )  quasiNewtonDescent.setSigma( _pparser.getDouble("sigma") );
        if( _pparser.hasVariable("beta") ) quasiNewtonDescent.setBeta( _pparser.getDouble("beta") );
        quasiNewtonDescent.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );
        quasiNewtonDescent.apply( InitialGuessOnCurrentLevel, Path );
	solved = true;
      }
      
      // Newton method
      if( !solved && (_solvingMethod == NEWTON) ){
	
	if( _lagrangeSetup ){
          // rigid body motion handler
          typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
          RBMHandlerType CHandler( topology, startShapeOnCurrentLevel, NumOfSingleConstraints, _numOfFreeShells );

          //GenericLagrangeFunction< ConfiguratorType, RBMHandlerType >  L( E, CHandler );
          GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( dE, CHandler );
          GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( dE2, CHandler );
      
          // append Lagrange multipliers
          VectorType Argument, Destination;           
          VectorType lagrangeMultArg( NumOfConstraints, 1 ), lagrangeMultDest( NumOfConstraints, 1 );
          Argument.appendReference( InitialGuessOnCurrentLevel );
          Argument.appendReference( lagrangeMultArg );
          Destination.appendReference( Path );
          Destination.appendReference( lagrangeMultDest );
    
          // initialize solver and solve	
          if( !_quiet ) cerr <<"Start solving with Newton method (using Lagrange setup)... " << endl;
          NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( topology, dL, dL2, ConfiguratorType::Dim * _numOfFreeShells, NumOfConstraints, _solverType, _numOfDescentSteps, _stopCriterion );
          zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
          zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
          zeroFinder.setQuietMode( _quiet );      
          zeroFinder.apply( Argument, Destination ); 	
	  if( !_quiet ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Path ) << endl;
	}
	else{
	  // initialize solver and solve	
          if( !_quiet ) cerr <<"Start solving with Newton method... " << endl;
          //NewtonMinimizer<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( topology, E, dE, dE2, ConfiguratorType::Dim * _numOfFreeShells, _solverType, _numOfDescentSteps, _stopCriterion );
	  NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( topology, dE, dE2, ConfiguratorType::Dim * _numOfFreeShells, 0, _solverType, _numOfDescentSteps, _stopCriterion );
          zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
          zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
          zeroFinder.setQuietMode( _quiet );      
          zeroFinder.apply( InitialGuessOnCurrentLevel, Path ); 
	}
	
	solved = true;
      }
    
      if( !solved )
	throw aol::Exception ( "GeodesicOp::executeMultiLevel(): unknown solver method!", __FILE__, __LINE__ ); 
    
      // save solution      
      if( _saving ){
	if( !_quiet ) cerr <<"Save solution... " << endl;
        ostringstream savename;
        if( _levels > 1 )
          savename << _destDirectory << _savenameStem << "_level" << level << "_";    
	else
	  savename << _destDirectory << _savenameStem;       
        saveResults( startShapes[currentLevel], endShapes[currentLevel], Path, savename.str() );
      }
            
      // no prolongation on last level
      if( level == _levels ){
	// final results
        E.apply( Path, finalEnergy );
        if( !_quiet ) cerr << "Final energy = " <<  finalEnergy[0] * (_numOfFreeShells+1) << endl;
        continue;
      }

      // compute prolongation
      if( !_quiet ) cerr <<"Prolongation... " << endl;
      computeProlongation( startShapes, endShapes, prolongationInfos, currentLevel, Path );
      
      InitialGuessOnCurrentLevel.reallocate( Path ); 
      InitialGuessOnCurrentLevel = Path;
      
      startShapeOnCurrentLevel.reallocate( ConfiguratorType::Dim, startShapes[currentLevel-1].getNumVertices() );
      startShapes[currentLevel-1].toVector( startShapeOnCurrentLevel );
      endShapeOnCurrentLevel.reallocate( ConfiguratorType::Dim, startShapes[currentLevel-1].getNumVertices() );
      endShapes[currentLevel-1].toVector( endShapeOnCurrentLevel );
    
    }// end multilevel
    
    return finalEnergy[0] * (_numOfFreeShells+1);
  }
  
  //! Computes discrete geodesic via an alternating red-black optimization of path energy E = K \sum_k W[S_{k-1}, S_k];
  //! Here we assume the K is an odd number.
  //! Alternating between
  //!  - fixing all even shapes S_0, S_2, ... S_K and compute 3-point geodesics (S_i, S_{i+1}, S_{i+1}) for i = 0, 2, ..., K-2; hence update all odd shapes S_1, S_3, ..., S_{K-1}
  //!  - fixing all odd shapes S_1, S_3, ..., S_{K-1} and compute 3-point geodesics (S_i, S_{i+1}, S_{i+1}) for i = 1, 3, ..., K-3; hence update all odd shapes S_2, S_4, ..., S_{K-2}
  //! Note that within each of the two alternating steps we can parallelize all computations of thees 3-point-geodesics which is not done so far.
  RealType computeGeodesicRedBlack( const VectorType& startShape, const VectorType& endShape, int numOfShells, VectorType& Path, bool initialize = true ) const {    
         
    if( !_quiet ) cerr <<"Start red black optimization with " << numOfShells <<" shapes (using Newton method)..." << endl;
      
    // update boundary mask if necessary
    if( _fixBoundary )
      checkAndUpdateBdryMask( _topology->getGrid() );
    
    _numOfShells = numOfShells;
    _numOfFreeShells = _numOfShells - 2;    
    int NumOfSingleConstraints = _lagrangeSetup ? _NumOfSingleConstraints : 0;
    int NumOfConstraints = _numOfFreeShells * NumOfSingleConstraints;
    
    // reallocation if necessary
    if( (Path.numComponents() != (ConfiguratorType::Dim * _numOfFreeShells)) || (Path[0].size() != _topology->getNumVertices()) )
      Path.reallocate( ConfiguratorType::Dim * _numOfFreeShells, _topology->getNumVertices() );
    
    // initialization
    if( initialize ){
      if( !_quiet ) cerr << "Start initialization... "  << endl;
      initializePath( startShape, endShape, Path );
    }
    
    // full path includes start and end shape
    VectorType fullPath;
    for( int j = 0; j < ConfiguratorType::Dim; j++ )
      fullPath.appendReference( startShape[j] );
    for( int i = 0; i < _numOfFreeShells; i++ )
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
        fullPath.appendReference( Path[ ConfiguratorType::Dim*i + j ] );
    for( int j = 0; j < ConfiguratorType::Dim; j++ )
      fullPath.appendReference( endShape[j] );
    
    // save solution    
    if( _saving && _pparser.checkAndGetBool("saveInitialization") ){
      if( !_quiet ) cerr <<"Save initialization... " << endl;
      ostringstream savename;
      savename << _destDirectory << _savenameStem << "_initial";       
      saveResults( fullPath, savename.str() );
    }  
    
    int numOfRedSteps = (_numOfShells - 1) / 2;
    int numOfBlackSteps = numOfRedSteps - 1;
    if( numOfRedSteps + numOfBlackSteps + 2 != _numOfShells )
      throw aol::Exception ( "GeodesicOp::executeRedBlack(): wrong number of shells!", __FILE__, __LINE__ );
    
    // set timestep controller
    //typename aol::NewtonInfo<RealType>::TIMESTEP_CONTROLLER timestepController = aol::NewtonInfo<RealType>::ARMIJO;
    typename aol::NewtonInfo<RealType>::TIMESTEP_CONTROLLER timestepController = aol::NewtonInfo<RealType>::NEWTON_OPTIMAL;
      
    // alternating steps
    for( int step = 0; step < _alternatingSteps; step++ ){  
      
      if( !_quiet ) cerr << "=============================================================" << endl;
      if( !_quiet ) cerr << "Start alternating step " << step+1 << " of " << _alternatingSteps << endl;
      
      //! red
      if( !_quiet ) cerr << "-------------------------------" << endl;
      if( !_quiet ) cerr << "Do " << numOfRedSteps << " red steps... "<< endl;
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
      for( int i = 0; i < numOfRedSteps; i++ ){

        VectorType leftShape, rightShape;
	for( int j = 0; j < ConfiguratorType::Dim; j++ ){
	  leftShape.appendReference( fullPath[(2*i) * ConfiguratorType::Dim + j] );
	  rightShape.appendReference( fullPath[(2*i+2) * ConfiguratorType::Dim + j] );
	}
	
	Geod3Grad<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, leftShape, rightShape );
        Geod3Hess<ConfiguratorType, MembraneDeformationType, BendingDeformationType  >  dE2( *_topology, _pparser, leftShape, rightShape, _regularizeMatrix );
	if( _fixBoundary ){
	  dE.setBoundaryMask( *_bdryMask );
	  dE2.setBoundaryMask( *_bdryMask );
	}
	
	// solving
	VectorType singleSolution;
	for( int j = 0; j < ConfiguratorType::Dim; j++ )
	  singleSolution.appendReference( fullPath[ (2*i+1) * ConfiguratorType::Dim + j] );
        VectorType InitialGuess( singleSolution );

        NewtonMethodSolver<ConfiguratorType>( _topology->getGrid(), dE, dE2, _solverType, _numOfDescentSteps, _stopCriterion, _quiet, timestepController ).solve( InitialGuess, singleSolution, _lagrangeSetup );   
      
      }
      
      //! black
      if( !_quiet ) cerr << "-------------------------------" << endl;
      if( !_quiet ) cerr << "Do " << numOfBlackSteps << " black steps... "<< endl;
      for( int i = 0; i < numOfBlackSteps; i++ ){

        VectorType leftShape, rightShape;
	for( int j = 0; j < ConfiguratorType::Dim; j++ ){
	  leftShape.appendReference( fullPath[(2*i+1) * ConfiguratorType::Dim + j] );
	  rightShape.appendReference( fullPath[(2*i+3) * ConfiguratorType::Dim + j] );
	}
	
	Geod3Grad<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, leftShape, rightShape );
        Geod3Hess<ConfiguratorType, MembraneDeformationType, BendingDeformationType  >  dE2( *_topology, _pparser, leftShape, rightShape, _regularizeMatrix );
	if( _fixBoundary ){
	  dE.setBoundaryMask( *_bdryMask );
	  dE2.setBoundaryMask( *_bdryMask );
	}

	// solving
	VectorType singleSolution;
	for( int j = 0; j < ConfiguratorType::Dim; j++ )
	  singleSolution.appendReference( fullPath[ (2*i+2) * ConfiguratorType::Dim + j] );
        VectorType InitialGuess( singleSolution );

        NewtonMethodSolver<ConfiguratorType>( _topology->getGrid(), dE, dE2, _solverType, _numOfDescentSteps, _stopCriterion, _quiet, timestepController ).solve( InitialGuess, singleSolution, _lagrangeSetup );   
      
      }
    
      // save solution      
      if( _saving && !_pparser.checkAndGetBool("saveOnlyFinalResult") ){
	if( !_quiet ) cerr <<"Save red black solution... " << endl;
        ostringstream savename;
        savename << _destDirectory << _savenameStem << "_rb" << step << "_";       
        saveResults( fullPath, savename.str() );
      }    
      
    } // end alternating
    
    
    //! final relaxation   
    aol::Scalar<RealType> finalEnergy;
    GeodesicEnergy<ConfiguratorType, MembraneDeformationType, BendingDeformationType >     E( *_topology, _pparser, startShape, endShape, _numOfShells  );
    
    if( numOfBlackSteps > 0 ){ 
      if( !_quiet ) cerr << "=============================================================" << endl;
      if( !_quiet ) cerr << "Final relaxation with " << _numOfShells << " shapes..." << endl;
      if( !_quiet && _lagrangeSetup ) cerr << "Using Lagrange setup." << endl;
      
      // save initialization      
      if( _saving && _pparser.checkAndGetBool("saveInitialization") ){
	if( !_quiet ) cerr <<"Save initialization ... " << endl;
        ostringstream savename;
        savename << _destDirectory << _savenameStem << "_unrelaxed";  
        saveResults( fullPath, savename.str() );
      }   
      
      E.apply( Path, finalEnergy );
      if( !_quiet ) cerr << "Initial energy = " <<  finalEnergy[0] * (_numOfFreeShells+1) << endl;   
      
      GeodesicGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, startShape, endShape, _numOfShells  );
      GeodesicHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE2( *_topology, _pparser, startShape, endShape, _numOfShells  );

      // fix boundary?
      if( _fixBoundary ){
        dE.setBoundaryMask( *_bdryMask );
        dE2.setBoundaryMask( *_bdryMask );
      }
    
      // rigid body motion handler
      typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
      RBMHandlerType CHandler( *_topology, startShape, NumOfSingleConstraints, _numOfFreeShells );

      GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( dE, CHandler );
      GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( dE2, CHandler );
   
      // append Lagrange multipliers
      VectorType Destination;           
      VectorType lagrangeMult( NumOfConstraints, 1 );
      Destination.appendReference( Path );
      if( _lagrangeSetup ) Destination.appendReference( lagrangeMult );
      VectorType Argument( Destination );
    
      // initialize solver and solve
      RealType stopCriterion = max( _stopCriterion, 1e-8 * sqrt( _topology->getNumVertices() * _numOfFreeShells) );
      if( !_quiet ) cerr <<"Start solving with stopping criterion " << stopCriterion << " ... " << endl;
      NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( *_topology, dL, dL2, ConfiguratorType::Dim * _numOfFreeShells, NumOfConstraints, _solverType, _numOfDescentSteps, stopCriterion );
      zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
      zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
      zeroFinder.setQuietMode( _quiet );      
      zeroFinder.apply( Argument, Destination ); 
    
      if( !_quiet && _lagrangeSetup ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Path ) << endl;    
    }
     
    // save solution    
    if( _saving ){
      if( !_quiet ) cerr <<"Save final solution... " << endl;
      ostringstream savename;
       savename << _destDirectory << _savenameStem;       
      saveResults( fullPath, savename.str() );
    }
    
    // final results
    E.apply( Path, finalEnergy );
    if( !_quiet ) cerr << "Final energy = " <<  finalEnergy[0] * (_numOfFreeShells+1) << endl;   
    return finalEnergy[0] * (_numOfFreeShells+1);
  }
   
  //! Computes discrete geodesic via a cascadic optimization of path energy E = K \sum_k W[S_{k-1}, S_k];
  //! We have K = 2^L + 1 where L is the number of cascadic levels.
  //! The cascadic algorithm starts calculating one (= 2^0) 3-point geodesic (S_A, S^0_0, S_B) by optimizing S^0_0
  //! Then we compute two (=2^1) 3-point geodesics (S_A, S^1_0, S^0_1) and (S^0_1, S^1_1, S_B) by optimizing S^1_0 and S^1_1
  //! Hence on the l-th level we compute 2^l 3-point geodesics optimizing 2^l shapes S^l_0, S^l_1, ...
  //! After the L-th level we optimize the resulting K-point geodesic.
  RealType computeGeodesicCascadic( const VectorType& startShape, const VectorType& endShape, int cascadicLevels, VectorType& Path, bool initialize = true ) const {    
         
    if( !_quiet ) cerr <<"Start cascadic optimization..." << endl;
    
    if( cascadicLevels > 3 )
      cerr << "WARNING: very high cascadic level chosen!" << endl << endl;   
    
    // update boundary mask if necessary
    if( _fixBoundary )
      checkAndUpdateBdryMask( _topology->getGrid() );
    
    _numOfShells = pow(2, cascadicLevels ) + 1;
    _numOfFreeShells = _numOfShells - 2;    
    int NumOfSingleConstraints = _lagrangeSetup ? _NumOfSingleConstraints : 0;
    int NumOfConstraints = _numOfFreeShells * NumOfSingleConstraints;
    
    Path.reallocate( ConfiguratorType::Dim * _numOfFreeShells, _topology->getNumVertices() );
    
    // initialization
    if( initialize ){
      if( !_quiet ) cerr << "Start initialization... "  << endl;
      initializePath( startShape, endShape, Path );
    }
    
    VectorType fullPath;
    for( int j = 0; j < ConfiguratorType::Dim; j++ )
      fullPath.appendReference( startShape[j] );
    for( int i = 0; i < _numOfFreeShells; i++ )
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
        fullPath.appendReference( Path[ ConfiguratorType::Dim*i + j ] );
    for( int j = 0; j < ConfiguratorType::Dim; j++ )
      fullPath.appendReference( endShape[j] );
      
    // alternating steps
    for( int level = 0; level < cascadicLevels; level++ ){  
      
      if( !_quiet ) cerr << "=============================================================" << endl;
      if( !_quiet ) cerr << "Start cascadic iteration on level " << level+1 << " of " << cascadicLevels << endl;

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
      for( int i = 0; i < pow(2, level); i++ ){
	
	int aux = pow(2, cascadicLevels-level-1);
	if( !_quiet ) cerr << "Solve geodesic ( " << aux * (2*i + 0)<< ", " << aux * (2*i + 1) << ", " << aux * (2*i + 2) << " )... " << endl;
    
        VectorType leftShape, rightShape;
	for( int j = 0; j < ConfiguratorType::Dim; j++ ){
	  leftShape.appendReference(  fullPath[ aux * (2*i + 0) * ConfiguratorType::Dim + j] );
	  rightShape.appendReference( fullPath[ aux * (2*i + 2) * ConfiguratorType::Dim + j] );
	}
	
	Geod3Grad<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, leftShape, rightShape );
        Geod3Hess<ConfiguratorType, MembraneDeformationType, BendingDeformationType  >  dE2( *_topology, _pparser, leftShape, rightShape, _regularizeMatrix );
	if( _fixBoundary ){	  
	  dE.setBoundaryMask( *_bdryMask );
	  dE2.setBoundaryMask( *_bdryMask );
	}

	// initialization
	VectorType InitialGuess( leftShape, aol::STRUCT_COPY );
	if( (level==0) && (i==0) && _pparser.hasVariable("initialShape") ){
	  cerr  << "Load initial shape from " << _pparser.getString("initialShape") << endl;
	  _topolMesh.loadFromPLY( _pparser.getString("initialShape") );
          _topolMesh.toVector( InitialGuess );    
	}
	else{
	  InitialGuess = leftShape;	
	  // linear initialization
	  if( !_pparser.checkAndGetBool("initializeTrivially") ){
	    InitialGuess += rightShape;
	    InitialGuess /= 2.;
	  }
	}
	
	// solving 
	VectorType singleSolution;
	for( int j = 0; j < ConfiguratorType::Dim; j++ )
	  singleSolution.appendReference( fullPath[ aux * (2*i + 1) * ConfiguratorType::Dim + j] );
        NewtonMethodSolver<ConfiguratorType>( _topology->getGrid(), dE, dE2, _solverType, _numOfDescentSteps, _stopCriterion, _quiet ).solve( InitialGuess, singleSolution, _lagrangeSetup );         
      }
      
    } // end cascadic
    
    
    
    //! final relaxation   
    aol::Scalar<RealType> energy;
    GeodesicEnergy<ConfiguratorType, MembraneDeformationType, BendingDeformationType > E( *_topology, _pparser, startShape, endShape, _numOfShells  );

    if( cascadicLevels > 1 ){   
      
      if( !_quiet ) cerr << "=============================================================" << endl;
      if( !_quiet ) cerr << "Final relaxation with " << _numOfShells << " shapes..." << endl;
      if( !_quiet && _lagrangeSetup ) cerr << "Using Lagrange setup." << endl;
      
      // intial energy
      E.apply( Path, energy );
      if( !_quiet ) cerr << "Initial energy = " <<  energy[0] * (_numOfFreeShells+1) << endl;       
    
      // save initialization      
      if( _saving && _pparser.checkAndGetBool("saveInitialization") ){
	if( !_quiet ) cerr <<"Save initialization ... " << endl;
        ostringstream savename;
        savename << _destDirectory << _savenameStem << "_initial_relaxation";  
        saveResults( fullPath, savename.str() );
      }     
      
      // setup gradient and Hessian
      GeodesicGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, startShape, endShape, _numOfShells  );
      GeodesicHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE2( *_topology, _pparser, startShape, endShape, _numOfShells  );
    
      // fix boundary?
      if( _fixBoundary ){
        dE.setBoundaryMask( *_bdryMask );
        dE2.setBoundaryMask( *_bdryMask );
      }
    
      // rigid body motion handler
      typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
      RBMHandlerType CHandler( *_topology, startShape, NumOfSingleConstraints, _numOfFreeShells );

      GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( dE, CHandler );
      GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( dE2, CHandler );
   
      // append Lagrange multipliers
      VectorType Destination;           
      VectorType lagrangeMult( NumOfConstraints, 1 );
      Destination.appendReference( Path );
      if( _lagrangeSetup ) Destination.appendReference( lagrangeMult );
      VectorType Argument( Destination );
    
      // initialize solver and solve
      RealType stopCriterion = max( _stopCriterion, 1e-8 * sqrt( _topology->getNumVertices() * _numOfFreeShells) );
      if( !_quiet ) cerr <<"Start solving with stopping criterion " << stopCriterion << " ... " << endl;
      NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( *_topology, dL, dL2, ConfiguratorType::Dim * _numOfFreeShells, NumOfConstraints, _solverType, _numOfDescentSteps, stopCriterion );
      zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
      zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
      zeroFinder.setQuietMode( _quiet );      
      zeroFinder.apply( Argument, Destination );   
      
      if( !_quiet && _lagrangeSetup ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Path ) << endl;
    }
   
    // save solution    
    if( _saving ){
      if( !_quiet ) cerr <<"Save final solution... " << endl;
      ostringstream savename;
      savename << _destDirectory << _savenameStem;       
      saveResults( fullPath, savename.str() );
    }
              
    // final energy
    E.apply( Path, energy );
    if( !_quiet ) cerr << "Final energy = " <<  energy[0] * (_numOfFreeShells+1) << endl;       
    return energy[0] * (_numOfFreeShells+1);
  }
  
protected:
  void initializePathLinearly( const VectorType& startShape, const VectorType& endShape, VectorType& Path ) const {
    if( !_quiet ) cerr <<"Initialize path linearly... " << endl;
    VectorType difference( endShape );
    difference -= startShape;
    for( int i = 0; i < _numOfFreeShells; i++ )
      for( int j = 0; j < ConfiguratorType::Dim; j++ )
	Path[ ConfiguratorType::Dim*i + j ].addMultiple( difference[j], 1.*(i+1)/(_numOfFreeShells+1) );
  }
  
  // initialization order: try (1) to load path, (2) load single shape, (3) initialize linearly (4) initialize trivially with start shape
  void initializePath( const VectorType& startShape, const VectorType& endShape, VectorType& Path ) const {
    
      // load initialization from file
      if( _pparser.hasVariable("initialShapeStem") ){
        loadPath( Path );
	return;
      }      
      
      // load initial shape ?
      VectorType initialShape( startShape, aol::STRUCT_COPY );
      if( _pparser.hasVariable("initialShape") ){
        if( !_quiet ) cerr <<"Load initial shape " << _pparser.getString("initialShape") << endl;
        _topolMesh.loadFromPLY( _pparser.getString("initialShape") );
        _topolMesh.toVector( initialShape );      
      }
      else
        initialShape = startShape;
    
      // initialize with start shape on coarsest level
      for( int i = 0; i < _numOfFreeShells; i++ )
        for( int j = 0; j < ConfiguratorType::Dim; j++ )
          Path[ ConfiguratorType::Dim*i + j ] = initialShape[ j ];  
	
      if( _pparser.hasVariable("initialShape") || _pparser.checkAndGetBool("initializeTrivially") ){
	if( !_quiet ) cerr <<"Initialized with constant shape." << endl;
	return;
      }
      
      // initialize linearly  
      initializePathLinearly( startShape, endShape, Path );
  }
  
  // saving results
  void saveResults( const MeshType& startShape, const MeshType& endShape, const VectorType& Path, string filename ) const {   
    
    if( Path.numComponents() / ConfiguratorType::Dim != _numOfShells - 2 )
      throw aol::Exception ( "GeodesicOp::saveResults(): sizes do not match(1)!", __FILE__, __LINE__ );
    
    for( int k = 0; k < _numOfShells; k++ ){
      MeshType mesh;
      ostringstream savename;
      savename << filename << ( k < 10 ? "0" : "" ) << k << ".ply" << ends;
      
      // start shape
      if( k == 0 )
        startShape.saveAsPLY( savename.str() );;
      
      // end shape
      if( k == (_numOfShells-1) )
	endShape.saveAsPLY( savename.str() );;
      
      // intermediate shape
      if( (k>0) && (k < _numOfShells-1) ){
	MeshType auxMesh( startShape );
	VectorType temp;
	for( int j = 0; j < ConfiguratorType::Dim; j++ )
	  temp.appendReference( Path[ (k-1)*ConfiguratorType::Dim + j ] );
	auxMesh.fromVector( temp );
        auxMesh.saveAsPLY( savename.str() );
      }      
    }
  }
  
  // saving results
  void saveResults(  const VectorType& Path, string filename ) const {   
    
    if( Path.numComponents() / ConfiguratorType::Dim != _numOfShells )
      throw aol::Exception ( "GeodesicOp::saveResults(): sizes do not match(2)!", __FILE__, __LINE__ );
    
    for( int k = 0; k < _numOfShells; k++ ){
      ostringstream savename;
      savename << filename << ( k < 10 ? "0" : "" ) << k << ".ply" << ends;
	MeshType auxMesh( _topolMesh );
	VectorType temp;
	for( int j = 0; j < ConfiguratorType::Dim; j++ )
	  temp.appendReference( Path[ k*ConfiguratorType::Dim + j ] );
	auxMesh.fromVector( temp );
        auxMesh.saveAsPLY( savename.str() ); 
    }
  }
  
  // build up progressive mesh hierarchy
  void computeDecimation( MeshContainer& startShapes, MeshContainer& endShapes, aol::RandomAccessContainer< ProlongInfo<RealType> >& prolongationInfos ) const {
    
    if( !_quiet ) cerr<<"Original mesh: "<<startShapes[0].getNumVertices()<<" vertices, "<<startShapes[0].getNumFaces()<<" faces.\n";

    //! generate progressive meshes
    aol::StopWatch dectime;
    dectime.start();
    if( !_quiet ) cerr<<"Start decimation...\n";
    for( int level = 0; level < _levels - 1; level++ ){

      SimulProgMesh< MeshType, RealType > progMeshGenerator;
      DecimationInfo<RealType> dInfo( startShapes[level].getNumVertices() );

      // calculate decimations
      MeshContainer coarseMeshes;
      coarseMeshes.pushBack( startShapes[level] );
      coarseMeshes.pushBack( endShapes[level] );
      progMeshGenerator.calcDecimation( coarseMeshes, _theta, dInfo );

      // store coarse meshes
      startShapes.pushBack( coarseMeshes[0] );
      endShapes.pushBack( coarseMeshes[1] );

      // calculate and store prolongations for initial and end shapes
      if( !_quiet ) cerr<<"Prolongate...";
      progMeshGenerator.calcProlongation( coarseMeshes[0], startShapes[level], dInfo, prolongationInfos[2*level], 0 );
      progMeshGenerator.calcProlongation( coarseMeshes[1], endShapes[level],   dInfo, prolongationInfos[2*level+1], 1 );

      if( !_quiet ) cerr<<"done.\n"<<level<<"th level: "<<coarseMeshes[0].getNumVertices()<<" vertices, "<<coarseMeshes[0].getNumFaces()<<" faces.\n";
    }
    dectime.stop();
    if( !_quiet ) cerr<<"\nDecimation done in "<<dectime.elapsedWallClockTime ()<<" seconds.\n\n";
  }
  
  // compute prolongation needed for multilevel method
  void computeProlongation( const MeshContainer& startShapes, 
			    const MeshContainer& /*endShapes*/, 
			    const aol::RandomAccessContainer< ProlongInfo<RealType> >& prolongationInfos, 
			    int coarseLevel, 
			    VectorType& SolutionOnCurrentLevel ) const {
      
      int fineLevel = coarseLevel - 1;			      
      if( !_quiet ) cerr<<"Start prolongation from level " << _levels - coarseLevel << " to " << _levels - fineLevel << endl;
      VectorType prolongFromLeft( ConfiguratorType::Dim * _numOfFreeShells,  startShapes[fineLevel].getNumVertices() ),
                 prolongFromRight( ConfiguratorType::Dim * _numOfFreeShells, startShapes[fineLevel].getNumVertices() );

      MeshType coarseMesh( startShapes[coarseLevel] );
      for( int i = 0; i < _numOfFreeShells; i++ ){
        SimulProgMesh< MeshType, RealType > progMeshGenerator;
        // initialize coarse mesh with solution
        VectorType coarseMV;
        for( int j = 0; j < ConfiguratorType::Dim; j++ )
          coarseMV.appendReference( SolutionOnCurrentLevel[ ConfiguratorType::Dim*i+j ] );
        coarseMesh.fromVector( coarseMV );
        // prolongate from left
        VectorType auxL, auxR;
        for( int j = 0; j < ConfiguratorType::Dim; j++ )
          auxL.appendReference( prolongFromLeft[ ConfiguratorType::Dim*i+j ] );
        progMeshGenerator.prolongate( coarseMesh, auxL, prolongationInfos[2 * fineLevel] );
        // prolongate from right
        for( int j = 0; j < ConfiguratorType::Dim; j++ )
          auxR.appendReference( prolongFromRight[ ConfiguratorType::Dim*i+j ] );
        progMeshGenerator.prolongate( coarseMesh, auxR, prolongationInfos[2 * fineLevel + 1] );

      }
      if( !_quiet ) cerr<<"done." << endl;

      // comupte new initial values as linear interpolation
      if( !_quiet ) cerr<<"Start re-initialisation...";
      SolutionOnCurrentLevel.reallocate( ConfiguratorType::Dim * _numOfFreeShells, startShapes[fineLevel].getNumVertices() );

      // intermediate meshes
      RealType deformations = static_cast<RealType>( _numOfShells - 1 );
      for( int i = 0; i < _numOfFreeShells; i++ )	
        for( int j = 0; j<ConfiguratorType::Dim; j++ ){
          SolutionOnCurrentLevel[ConfiguratorType::Dim*i+j].addMultiple( prolongFromLeft[ConfiguratorType::Dim*i+j],  (deformations - i - 1) / deformations );
          SolutionOnCurrentLevel[ConfiguratorType::Dim*i+j].addMultiple( prolongFromRight[ConfiguratorType::Dim*i+j],  (i + 1) / deformations  );
        }
      
      if( !_quiet ) cerr<<"done." << endl;
  }
  
  //
  void checkAndUpdateBdryMask( const MeshType& mesh ) const {
    if( _bdryMask->size() == mesh.getNumVertices() )
      return;
    delete _bdryMask;
    _bdryMask = new aol::BitVector( mesh.getNumVertices() );
    mesh.fillBoundaryMask( *_bdryMask );
  }
  
  // load path from file
  void loadPath( VectorType& Path ) const {
    for( int k = 0; k < _numOfFreeShells; k++ ){
      // loadname
      ostringstream loadname;
      loadname <<  _pparser.getString("initialShapeStem");
      if( k < 9 ) loadname <<"0";
      loadname << k + 1 <<".ply";      
      cerr << "Load " << loadname.str() << endl;

      // load
      MeshType mesh;
      mesh.loadFromPLY( loadname.str() );
      VectorType shell;
      for( int i = 0; i < 3; i++ )
         shell.appendReference( Path[k*3+i] ); 
      mesh.toVector( shell );   
    }
  }
  
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class to shoot discrete geodesics.
//! \author Heeren
//!
//! Given p as "position" and q as "variational shape" compute EXP^K_p( v ), with v := alpha * (q - p),
//! where K is the number of shooting steps and alpha is an emphasizing factor
template<typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType>
class ShootingOp{
  
   typedef typename ConfiguratorType::RealType RealType;
   typedef aol::MultiVector<RealType> VectorType;
   typedef typename ConfiguratorType::InitType MeshType;
    
   typedef SparseMatrix<RealType> SubMatrixType;
   typedef SparseBlockMatrix<SubMatrixType> BlockMatrixType;
   
   const aol::ParameterParser& _pparser;
   const MeshTopologySaver<MeshType>* _topology;
   mutable MeshType _topolMesh;    
   VectorType *_position, *_variationShape;
   RealType _factor;   
   int _shootingSteps;   
   int _solverType;
   int _numOfDescentSteps;
   RealType _stopCriterion;
   char _destDirectory[1024], _savenameStem[1024];
   bool _fixBoundary;
   aol::BitVector *_bdryMask;  
   bool _regularizeMatrix, _lagrangeSetup;
   bool _quiet, _deleteTopology, _saving;
   
   const static int _NumOfConstraints = 6;
   
public:
  ShootingOp( const aol::ParameterParser& pparser ) 
    : _pparser( pparser), 
    _topolMesh( _pparser.getString("variationShape") ), 
    _factor( _pparser.getDoubleOrDefault( "emphFactorShooting", 1.0 ) ),
    _shootingSteps( pparser.getInt("numOfShootingSteps") ),
    _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
    _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
    _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
    _fixBoundary( _pparser.checkAndGetBool("fixBoundary") ),
    _bdryMask( NULL ),
    _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
    _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
    _quiet( false ),
    _deleteTopology( true ),
    _saving( true ){
      
      // dymanic allocating
      _topology = new MeshTopologySaver<MeshType>( _topolMesh );
      _position = new VectorType( 3, _topolMesh.getNumVertices() );
      _variationShape = new VectorType( 3, _topolMesh.getNumVertices() );
      
      _topolMesh.toVector( *_variationShape );
      _topolMesh.loadFromPLY( _pparser.getString("position") );
      _topolMesh.toVector( *_position );
      
      pparser.getString( "destDirectory", _destDirectory );
      pparser.getString( "savenameStem", _savenameStem );
      
      // boundary mask? 
      if( _fixBoundary ){
        if( pparser.hasVariable("bdryMask") ){
	  cerr << "Load boundary mask..." << endl;	
	  _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
        }
        else{
	  cerr << "Set boundary mask..." << endl;	
	   _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	   _topolMesh.fillBoundaryMask( *_bdryMask );
        }
      } 
      
      printSettings();
    }
  
  ShootingOp( const aol::ParameterParser& pparser, const MeshTopologySaver<MeshType>& Topology, const VectorType& variationShape ) 
    : _pparser( pparser),    
    _topology( &Topology ),
    _topolMesh( Topology.getGrid() ),
    _position( NULL ),
    _variationShape( NULL ),
    _factor( 1.0 ),
    _shootingSteps( -1 ),
    _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
    _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
    _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
    _fixBoundary( _pparser.checkAndGetBool("fixBoundary") ),
    _bdryMask( NULL ),
    _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
    _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
    _quiet( true ),
    _deleteTopology( false ),
    _saving( false ){ 
      _variationShape = new VectorType( variationShape );
      
      // boundary mask? 
      if( _fixBoundary ){
        if( pparser.hasVariable("bdryMask") ){
	  cerr << "Load boundary mask..." << endl;	
	  _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
        }
        else{
	  cerr << "Set boundary mask..." << endl;	
	   _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	   _topolMesh.fillBoundaryMask( *_bdryMask );
        }
      } 
    }
  
  ~ShootingOp() {
    if( _position )
      delete _position;
    if( _variationShape )    
      delete _variationShape;
    if( _deleteTopology && _topology )
      delete _topology;
    if( _fixBoundary ) 
      delete _bdryMask;
  }
  
  void setBoundaryMask( const BitVector& Mask ){
    if( !_quiet ) cerr << "Set boundary mask..." << endl;	
    if( _fixBoundary )
      delete _bdryMask;
    _bdryMask = new aol::BitVector( Mask );
    _fixBoundary = true;
  }
  
  void setSaveMode( const aol::ParameterParser& pparser ){
    pparser.getString( "destDirectory", _destDirectory );
    pparser.getString( "savenameStem", _savenameStem );
    _saving = true;
  }
  
  void printSettings() const {
    cerr << endl << "===========================" << endl;
    cerr << "Solving flags:" << endl;	
    cerr << " - max number of steps = " << _numOfDescentSteps << endl;
    cerr << " - stoping criterion = " << _stopCriterion << endl;
    if( _lagrangeSetup )
      cerr << " - using Lagrange setup" << endl;
    if( _regularizeMatrix )
      cerr << " - regularize matrix by adding eps*Id" << endl;
    if( _fixBoundary )
      cerr << " - regularize matrix by adding eps*Id" << endl;
    cerr << "===========================" << endl << endl;
  }

  void setQuietMode( bool Quiet ){
    _quiet = Quiet;
  }
  
  //! Given p as "position" and q as "variational shape" compute p + K*v, 
  //! where v := alpha * (q - p)  and K is the number of shooting steps  
  void executeLinearly( ) const {
    if( !_position )
      throw aol::Exception ( "ShootingOp::executeLinearly(): Position not set!", __FILE__, __LINE__ );
    if( !_variationShape )
      throw aol::Exception ( "ShootingOp::executeLinearly(): Variation not set!", __FILE__, __LINE__ );
    if( !_quiet ) cerr << "Start to apply linear shooting operator " << _shootingSteps << " times..." << endl;
    
    executeLinearly( *_position, *_variationShape, _shootingSteps, _factor );
  }
  
  //! Given a p as "position" and q as "variational shape" compute p + K*v, 
  //! where v := alpha * (q - p)  and K is the number of shooting steps
  void executeLinearly(  const VectorType& position, const VectorType& variationShape, int shootingSteps, RealType alpha = 1. ) const {
           
    if( !_topology )
      throw aol::Exception ( "ShootingOp::executeLinearly(): Topology not set!", __FILE__, __LINE__ );
    
    // S = S_0 + alpha (S_1 - S_0) = (1 - alpha) S_0 + alpha S_1
    VectorType Variation( variationShape ), Position( position );
    Variation -= position;
    
    // linear shooting 
    for( int k = 0; k < shootingSteps; k++ ){
      
      Position.addMultiple( Variation, alpha );      
      
      // saving
      if( _saving && ( (k+1) % _pparser.getIntOrDefault("saveShootingSteps", 1) == 0 ) ){      
	  _topolMesh.fromVector( Position );
	  int number = k / _pparser.getIntOrDefault("saveShootingSteps", 1);
	  ostringstream name;
          name << _destDirectory << _savenameStem   << number << ".ply" << ends;
	  _topolMesh.saveAsPLY( name.str() );	
      }
    }    
  }
  
  //! Given p as "position" and q as "variational shape" compute EXP^K_p( v ), 
  //! where v := alpha * (q - p)  and K is the number of shooting steps  
  void execute( ) const {
    if( !_position )
      throw aol::Exception ( "ShootingOp::execute(): Position not set!", __FILE__, __LINE__ );
    if( !_variationShape )
      throw aol::Exception ( "ShootingOp::execute(): Variation not set!", __FILE__, __LINE__ );
    if( !_quiet ) cerr << "Start to apply shooting operator " << _shootingSteps << " times..." << endl;
    if( !_quiet ) cerr << "Scale variation by " << _factor << "." << endl << endl;
    
    execute( *_position, *_variationShape, _shootingSteps, _factor );
  }
  
  //! Given p as "position" and q as "variational shape" compute EXP^K_p( v ), 
  //! where v := alpha * (q - p)  and K is the number of shooting steps
  void execute( const VectorType& position, const VectorType& variationShape, int shootingSteps, RealType alpha = 1. ) const {    
         
    if( !_topology )
      throw aol::Exception ( "ShootingOp::execute(): Topology not set!", __FILE__, __LINE__ );
    
    // S = S_0 + alpha (S_1 - S_0) = (1 - alpha) S_0 + alpha S_1
    VectorType Position( position ), nextPosition( position );
    nextPosition *= 1. - alpha;
    nextPosition.addMultiple( variationShape, alpha );
        
    VectorType shootedVector( *_variationShape, aol::FLAT_COPY );
    
    // shooting 
    for( int k = 0; k < shootingSteps; k++ ){
      if( !_quiet ) cerr << "Shoot vector: " << k+1 << " of " << shootingSteps << " steps..." << endl;
      
      // trivial initialization
      shootedVector = nextPosition;
        
      Exp2Energy<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dF( *_topology, _pparser, Position, nextPosition );
      Exp2Deriv<ConfiguratorType, MembraneDeformationType, BendingDeformationType > dF2( *_topology, _pparser, nextPosition, _regularizeMatrix );
      
      // fix boundary?
      if( _fixBoundary ){
	  dF.setBoundaryMask( *_bdryMask );
	  dF2.setBoundaryMask( *_bdryMask );
      }
	
      // solving
      typename aol::NewtonInfo<RealType>::TIMESTEP_CONTROLLER timestepController = aol::NewtonInfo<RealType>::ARMIJO;
      //typename aol::NewtonInfo<RealType>::TIMESTEP_CONTROLLER timestepController = aol::NewtonInfo<RealType>::NEWTON_OPTIMAL;
      NewtonMethodSolver<ConfiguratorType>( *_topology, position, dF, dF2, _solverType, _numOfDescentSteps, _stopCriterion, _quiet, timestepController, 0.1, _pparser.getIntOrDefault("numDerivativeHoldingSteps", 1) ).solve( nextPosition, shootedVector, _lagrangeSetup );	
      
      // saving
      if( _saving && ( (k+1) % _pparser.getIntOrDefault("saveShootingSteps", 1) == 0 ) ){    
	  _topolMesh.fromVector( shootedVector );
	  int number = k / _pparser.getIntOrDefault("saveShootingSteps", 1);
	  ostringstream vtkname;
	  if( _pparser.checkAndGetBool("saveAsVTK") ){
            vtkname << _destDirectory << _savenameStem   << number << ".vtk" << ends;
	    //cerr << vtkname.str() << endl;
	    _topolMesh.saveAsVTK( vtkname.str() );	    
	  }
	  else{
	    vtkname << _destDirectory << _savenameStem   << number << ".ply" << ends;
	    //cerr << vtkname.str() << endl;
	    _topolMesh.saveAsPLY( vtkname.str() );	
	  }
      }
      
      // update
      Position = nextPosition;
      nextPosition = shootedVector;   
      
      if( !_quiet ) cerr << "--------------------------------" << endl;
    }    
    
    if( _saving )
      saveFinalResult();

  }
  
  const VectorType& getVariationShape() const{
    return *_variationShape;
  }
  
  // save as ply (if required: save as obj or vtk)
  void saveFinalResult() const {
    if( _topolMesh.getNumVertices() == 0 )
      throw aol::Exception ( "ShootingOp::saveFinalResult(): topolMesh not set!", __FILE__, __LINE__ );
    
    _topolMesh.fromVector( *_variationShape );
    ostringstream savename;
    savename << _destDirectory << _savenameStem << ".ply" << ends;
    _topolMesh.saveAsPLY( savename.str() );
    if( _pparser.checkAndGetBool("saveAsOBJ") ){
      ostringstream objname;
      objname << _destDirectory << _savenameStem <<".obj" << ends;
      _topolMesh.saveAsOBJ( objname.str() );
    }
    if( _pparser.checkAndGetBool("saveAsVTK") ){
      ostringstream vtkname;
      vtkname << _destDirectory << _savenameStem << ".vtk" << ends;
      _topolMesh.saveAsVTK( vtkname.str() );      
    }
  }
  
};

//!===============================================================================================================================
//!===============================================================================================================================
//! \brief Class to compute discrete parallel transport via Schilds ladder.
//! \author Heeren
//!
//! Class computes parallel transport of a variational shape along a given path; possibly scaled by some epsilon.
//! Schilds ladder is based of a sequence of geodesic parallelograms, where diagonal geodesics in each parallelogram are given by N-point geodesics.
template<typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType>
class ParallelTransportOp{
  
   typedef typename ConfiguratorType::RealType RealType;
   typedef aol::MultiVector<RealType> VectorType;
   typedef typename ConfiguratorType::InitType MeshType;
    
   typedef SparseMatrix<RealType> SubMatrixType;
   typedef SparseBlockMatrix<SubMatrixType> BlockOp;
   
   const aol::ParameterParser& _pparser;
   const MeshTopologySaver<MeshType>* _topology;
   mutable MeshType _topolMesh;   
   int _solverTypeInterpolation, _solverTypeExtrapolation;
   int _numOfDescentSteps;
   RealType _stopCriterion;
   bool _fixBoundary;
   aol::BitVector *_bdryMask;  
   bool _regularizeMatrix, _lagrangeSetup;   
   VectorType* _path;    
   VectorType* _transportedVariationShape;
   RealType _scalingParTransport;   
   bool _saveOnlyFinalResult, _saveIntermediateMeshes, _saveAsOBJ, _quiet;
   char _directory[1024], _savenameStem[1024]; 
   bool _deleteTopology, _saving;
    
public:
  ParallelTransportOp( const aol::ParameterParser& pparser ) 
      : _pparser( pparser ),
      _topolMesh( pparser.getString("variationShape") ),
      _solverTypeInterpolation( pparser.getIntOrDefault("solverTypeInterpolation", LU) ),
      _solverTypeExtrapolation( pparser.getIntOrDefault("solverTypeExtrapolation", LU) ),
      _numOfDescentSteps( pparser.getIntOrDefault("numOfDescentSteps",1000) ),
      _stopCriterion( pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
      _fixBoundary( _pparser.checkAndGetBool("fixBoundary") || _pparser.hasVariable("bdryMask") ),
      _bdryMask(NULL),
      _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
      _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
      _scalingParTransport( pparser.getDoubleOrDefault("scaledParallelTransport", 0.) ),
      _saveOnlyFinalResult( pparser.getIntOrDefault("saveOnlyFinalResult", 0) ),
      _saveIntermediateMeshes( pparser.getIntOrDefault("saveIntermediateMesh", 0) ),
      _saveAsOBJ( pparser.getIntOrDefault("saveAsOBJ", 0) ),
      _quiet( false ),
      _deleteTopology( true ),
      _saving( true ){  
	
	// read in the directory, where results are to be saved
        pparser.getString( "destDirectory", _directory );
	pparser.getString( "savenameStem", _savenameStem );
	
    // load path
    _topology = new MeshTopologySaver<MeshType>( _topolMesh );
    _path = new VectorType( 3 * pparser.getInt("numOfShells"), _topolMesh.getNumVertices() );   
    _transportedVariationShape = new VectorType( 3, _topolMesh.getNumVertices() );
    
    // load or initialize path
    loadPath( pparser.hasVariable("path") );
    
    // boundary mask? 
    if( _fixBoundary ){
      if( pparser.hasVariable("bdryMask") ){
	cerr << "Load boundary mask..." << endl;	
	 _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
      }
      else{
	cerr << "Set boundary mask..." << endl;	
	 _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	 _topolMesh.fillBoundaryMask( *_bdryMask );
      }
    }
    
    // print solving flags
    printSettings();
  }  
  
  ParallelTransportOp( const aol::ParameterParser& pparser, const MeshTopologySaver<MeshType>& Topology, bool Quiet = true ) 
      : _pparser( pparser ),
      _topology( &Topology ),
      _solverTypeInterpolation( pparser.getIntOrDefault("solverTypeInterpolation", LU) ),
      _solverTypeExtrapolation( pparser.getIntOrDefault("solverTypeExtrapolation", LU) ),
      _numOfDescentSteps( pparser.getIntOrDefault("numOfDescentSteps",1000) ),
      _stopCriterion( pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
      _fixBoundary( false ),
      _bdryMask(NULL),
      _regularizeMatrix( pparser.checkAndGetBool("regularizeMatrix") ),
      _lagrangeSetup( pparser.checkAndGetBool("LagrangeSetup") ),
      _path( NULL ),
      _quiet( Quiet ),
      _deleteTopology( false ),
      _saving( false ){   
	_transportedVariationShape = new VectorType( 3, _topology->getNumVertices() );
      }
  
  ~ParallelTransportOp(){ 
    if( _fixBoundary ) 
      delete _bdryMask;     
    if( _path )
      delete _path;    
    if( _transportedVariationShape ) 
      delete _transportedVariationShape;
    if( _deleteTopology && _topology )
      delete _topology;
  }
  
  void setBoundaryMask( const BitVector& Mask ){
    if( !_quiet ) cerr << "Set boundary mask..." << endl;	
    _fixBoundary = true;
    _bdryMask = new aol::BitVector( Mask );
  }
  
  void printSettings() const {
    cerr << endl << "===========================" << endl;
    cerr << "Solving flags:" << endl;	
    cerr << " - max number of steps = " << _numOfDescentSteps << endl;
    cerr << " - stoping criterion = " << _stopCriterion << endl;
    if( _lagrangeSetup )
      cerr << " - using Lagrange setup" << endl;
    if( _regularizeMatrix )
      cerr << " - regularize matrix by adding eps*Id" << endl;
    if( _bdryMask )
      cerr<< " - fixing boundary" << endl;
    cerr << "===========================" << endl << endl;
  }
  
  
  //! computes parallel transport of a variational shape along a given path; possibly scaled by some epsilon.
  //! Parallel transport is realized via Schilds ladder.
  //! Here diagonal geodesics in each geodesic parallelogram are given by 3-point geodesics ( for higher order use execute() ).
  void execute3Point () const {
    _topolMesh.toVector( *_transportedVariationShape );    
    execute3Point( *_path, *_transportedVariationShape, _scalingParTransport );    
  }
    
  //! computes parallel transport of a variational shape along a given path; possibly scaled by some epsilon.
  //! Parallel transport is realized via Schilds ladder.
  //! Here diagonal geodesics in each geodesic parallelogram are given by 3-point geodesics ( for higher order use execute() ).
  void execute3Point ( const VectorType& Path, const VectorType& VariationShape, RealType ScalingParTransport ) const {
    
    if( _topolMesh.getNumVertices() == 0 )
      createMeshFromGeometryAndTopology<MeshType>( *_topology, VariationShape, _topolMesh );
    
    int numShells = Path.numComponents() / ConfiguratorType::Dim ;

    // bring into more convenient form
    aol::RandomAccessContainer<VectorType> shells( numShells  );
    for ( int k = 0; k < numShells; k++ )
      for ( int i = 0; i < 3; i++ )
        shells[k].appendReference( Path[k*3+i] );
      
    // save path?
    if( _saving && _pparser.checkAndGetBool("savePath") ){
      if( !_quiet ) cerr<<"Save path...\n";
      for( int k = 0; k < numShells; k++ ){
        ostringstream savename;
        savename << _directory << "path" << ( k < 10 ? "0" : "" ) << k << ".ply" << ends;
        _topolMesh.fromVector( shells[k] );
        _topolMesh.saveAsPLY( savename.str() );
      }    
    }
      
    if( !_quiet ) cerr << endl << "Perform parallel transport: " << endl;  
    
    // scale variation?
    VectorType scaledVariationShape( *_transportedVariationShape, aol::FLAT_COPY );
    scaledVariationShape = VariationShape;
    if( ScalingParTransport > 0 ){
        // var_eps = shell_0 + eps * (var - shell_0)
        if( !_quiet ) cerr<<"Scale variation with eps = " << ScalingParTransport << "...\n";
        scaledVariationShape *= ScalingParTransport;
	scaledVariationShape.addMultiple( shells[0], 1. - ScalingParTransport );
    }    
    
    // saving
    if( _saving && !_saveOnlyFinalResult )
      saveTransported( shells[0], scaledVariationShape, ScalingParTransport, 0 );
     
    // start parallel transport via Schild`s ladder
    for( int k = 1; k < numShells ; k++ ){

        if( !_quiet ) cerr<<"Start "<<k<<"th step of " << numShells-1 << " of parallel transport:\n";
	
	//! geodesic
        if( !_quiet ) cerr<<"1.) Calculate intermediate shape...\n";
        VectorType intermShell( scaledVariationShape, aol::STRUCT_COPY );
	Geod3Grad<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE( *_topology, _pparser, scaledVariationShape, shells[k] );
        Geod3Hess<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dE2( *_topology, _pparser, scaledVariationShape, shells[k], _regularizeMatrix );
	// fix boundary?
	if( _fixBoundary ){
	  dE.setBoundaryMask( *_bdryMask );
	  dE2.setBoundaryMask( *_bdryMask );
	}
	// solving
        VectorType InitialGuess;
	averageVector<VectorType>( scaledVariationShape, shells[k], InitialGuess );
        NewtonMethodSolver<ConfiguratorType>( *_topology, shells[k], dE, dE2, _solverTypeInterpolation, _numOfDescentSteps, _stopCriterion, _quiet ).solve( InitialGuess, intermShell, _lagrangeSetup );   
	
	//! shooting
        if( !_quiet ) cerr<<"2.) Calculate Exp2...\n";
	Exp2Energy<ConfiguratorType, MembraneDeformationType, BendingDeformationType >  dF( *_topology, _pparser, shells[k-1], intermShell );
        Exp2Deriv<ConfiguratorType, MembraneDeformationType, BendingDeformationType > dF2( *_topology, _pparser, intermShell, _regularizeMatrix );
	// fix boundary?
	if( _fixBoundary ){
	  dF.setBoundaryMask( *_bdryMask );
	  dF2.setBoundaryMask( *_bdryMask );
	}
	// solving
	extrapolateVector<VectorType>( shells[k-1], intermShell, InitialGuess );
	//InitialGuess = intermShell;
	NewtonMethodSolver<ConfiguratorType>( *_topology, shells[k], dF, dF2, _solverTypeExtrapolation, _numOfDescentSteps, _stopCriterion, _quiet ).solve( InitialGuess, scaledVariationShape, _lagrangeSetup ); 

	
        // saving
        if( _saving && !_saveOnlyFinalResult )
	  saveTransported( shells[k], scaledVariationShape, ScalingParTransport, k );
        
        if( _saving && _saveIntermediateMeshes ){
          ostringstream savename;
          savename << _directory << _savenameStem <<"_inter" << ( k < 10 ? "_0" : "_" ) << k << ".ply" << ends;
          _topolMesh.fromVector( intermShell );
          _topolMesh.saveAsPLY( savename.str() );
        }
        
        if( !_quiet ) cerr << "------------------------------------" << endl;
    }    
          
    // scaling?
    if( ScalingParTransport > 0 ){      
      RealType eps_inv = 1./ ScalingParTransport;
      if( !_quiet ) cerr<<"Rescale variation with 1 / eps = " << eps_inv << "...\n";
      scaledVariationShape *= eps_inv ;
      scaledVariationShape.addMultiple( shells[numShells-1], 1. - eps_inv );
    }
    
    if( _saving )
      saveFinalResult();
    
  }
    
  //! computes parallel transport of a variational shape along a given path; possibly scaled by some epsilon.
  //! Parallel transport is realized via Schilds ladder.
  //! Here diagonal geodesics in each geodesic parallelogram are given by N-point geodesics, where N is the level of geodesics used in parallelogram.
  void execute () const {
    _topolMesh.toVector( *_transportedVariationShape );    
    execute( *_path, *_transportedVariationShape, _scalingParTransport, _pparser.getIntOrDefault("cascadicLevelsInParTranp", 1) );    
  }
  
  //! computes parallel transport of a variational shape along a given path; possibly scaled by some epsilon.
  //! Parallel transport is realized via Schilds ladder.
  //! Here diagonal geodesics in each geodesic parallelogram are given by N-point geodesics, where N is the level of geodesics used in parallelogram.
  void execute ( const VectorType& Path, const VectorType& VariationShape, RealType ScalingParTransport, int levelInParallelogram ) const {
    
    if( levelInParallelogram < 2 )
      return execute3Point( Path, VariationShape, ScalingParTransport );

    int lengthOfLongGeodesicInParalelogram = pow( 2, levelInParallelogram ) + 1;
    int lengthOfShortGeodesicInParalelogram = pow( 2, levelInParallelogram - 1 ) + 1;
    int numOfShootingStepsInParallelogram = lengthOfShortGeodesicInParalelogram - 1;
    
    if( _topolMesh.getNumVertices() == 0 )
      createMeshFromGeometryAndTopology<MeshType>( *_topology, VariationShape, _topolMesh );
    
    int numShells = Path.numComponents() / ConfiguratorType::Dim ;

    // bring into more convenient form
    aol::RandomAccessContainer<VectorType> shells( numShells  );
    for ( int k = 0; k < numShells; k++ )
      for ( int i = 0; i < 3; i++ )
        shells[k].appendReference( Path[k*3+i] );
      
    if( !_quiet ) cerr << endl << "Perform parallel transport (fine): " << endl;  
    
    // save path?
    if( _saving && _pparser.checkAndGetBool("savePath") ){
      if( !_quiet ) cerr<<"Save path...\n";
      for( int k = 0; k < numShells; k++ ){
        ostringstream savename;
        savename << _directory << "path" << ( k < 10 ? "0" : "" ) << k << ".ply" << ends;
        _topolMesh.fromVector( shells[k] );
        _topolMesh.saveAsPLY( savename.str() );
      }    
    }
    
    // scale variation?
    VectorType scaledVariationShape( *_transportedVariationShape, aol::FLAT_COPY );
    scaledVariationShape = VariationShape;
    if( ScalingParTransport > 0 ){
        // var_eps = shell_0 + eps * (var - shell_0)
        if( !_quiet ) cerr<<"Scale variation with eps = " << ScalingParTransport << "...\n";
        scaledVariationShape *= ScalingParTransport;
	scaledVariationShape.addMultiple( shells[0], 1. - ScalingParTransport );
    }    
    
    // saving
    if( _saving && !_saveOnlyFinalResult )
      saveTransported( shells[0], scaledVariationShape, ScalingParTransport, 0 );
    
    // initialize geodesic operator
    GeodesicOp<ConfiguratorType, MembraneDeformationType, BendingDeformationType > geodOp( _pparser, *_topology );
    geodOp.setQuietMode( false );
    // fix boundary?
    if( _fixBoundary )
      geodOp.setBoundaryMask( *_bdryMask );	
    
    // initialize shooting op
    ShootingOp<ConfiguratorType, MembraneDeformationType, BendingDeformationType > shootOp( _pparser, *_topology, scaledVariationShape );
    shootOp.setQuietMode( false );
    // fix boundary?
    if( _fixBoundary )
      shootOp.setBoundaryMask( *_bdryMask );	
    
     
    // start parallel transport via Schild`s ladder
    for( int k = 1; k < numShells ; k++ ){
      
        if( !_quiet ) cerr<<"\n\n**************************************************" << endl;
        if( !_quiet ) cerr<<"Start "<<k<<"th step of " << numShells << " of fine parallel transport:\n";
	if( !_quiet ) cerr<<"**************************************************" << endl;
	
	//! long geodesic	
        if( !_quiet ) cerr<<"1.) Calculate intermediate shape of " << lengthOfLongGeodesicInParalelogram << "-point geodesic...\n";
        VectorType intermShell, Path;
	geodOp.computeGeodesicCascadic( shootOp.getVariationShape(), shells[k], levelInParallelogram, Path );	
	// intermediate shell is midpoint of geodesic (caution: only intermediate shells stored in Path!)
	for( int i = 0; i < ConfiguratorType::Dim; i++ )
	  intermShell.appendReference( Path[ConfiguratorType::Dim * ((lengthOfLongGeodesicInParalelogram-2) / 2) + i] );
  
	
	//! short geodesic
	if( !_quiet ) cerr<<"\n**************************************************" << endl;
        if( !_quiet ) cerr<<"2.) Calculate " << lengthOfShortGeodesicInParalelogram << "-point geodesic...\n";
	VectorType shortPath, shootingStartPoint;
	geodOp.computeGeodesicCascadic( shells[k-1], intermShell, levelInParallelogram-1, shortPath );	
	// shootingStartPoint is now last argument (caution: only intermediate shells stored in Path!)
	for( int i = 0; i < ConfiguratorType::Dim; i++ )
	  shootingStartPoint.appendReference( shortPath[ConfiguratorType::Dim * (lengthOfShortGeodesicInParalelogram - 3) + i] );
	
	
	//! shooting
	if( !_quiet ) cerr<<"\n**************************************************" << endl;
	if( !_quiet ) cerr<<"3.) Shoot " << numOfShootingStepsInParallelogram << "...\n";
	shootOp.execute( shootingStartPoint, intermShell, numOfShootingStepsInParallelogram );
	
        // saving
        if( _saving && !_saveOnlyFinalResult )
	  saveTransported( shells[k], shootOp.getVariationShape(), ScalingParTransport, k );
        
        if( _saving && _saveIntermediateMeshes ){
          ostringstream savename;
          savename << _directory << _savenameStem <<"_inter" << ( k < 10 ? "_0" : "_" ) << k << ".ply" << ends;
          _topolMesh.fromVector( intermShell );
          _topolMesh.saveAsPLY( savename.str() );
        }

    }    
    
    //
    if( !_quiet ) cerr << "==============================================================" << endl;
    if( !_quiet ) cerr << "==============================================================" << endl << endl;
          
    // scaling?
    if( ScalingParTransport > 0 ){      
      scaledVariationShape = shootOp.getVariationShape();
      RealType eps_inv = 1./ ScalingParTransport;
      if( !_quiet ) cerr<<"Rescale variation with 1 / eps = " << eps_inv << "...\n";
      scaledVariationShape *= eps_inv ;
      scaledVariationShape.addMultiple( shells[numShells-1], 1. - eps_inv );
    }
    
    if( _saving )
      saveFinalResult();
    
  }
  
  const VectorType& getTransportedShape() const{
    return *_transportedVariationShape;
  }
 
protected:     
  void saveTransported( const VectorType& ShapeOnPath, const VectorType& TransportedShape, RealType ScalingParTransport, int number) const {
    ostringstream savename;
    savename << _directory << _savenameStem << ( number < 10 ? "_0" : "_" ) << number << ".ply" << ends;
    VectorType temp( TransportedShape );
    RealType eps_inv = 1./ ScalingParTransport;
    temp *= eps_inv ;
    temp.addMultiple( ShapeOnPath, 1. - eps_inv );
    _topolMesh.fromVector( temp );
    _topolMesh.saveAsPLY( savename.str() );
  }
  
  void saveFinalResult() const {
    if( _topolMesh.getNumVertices() == 0 )
      throw aol::Exception ( "ParallelTransportOp::saveFinalResult(): topolMesh not set!", __FILE__, __LINE__ );
    _topolMesh.fromVector( *_transportedVariationShape );
    
    ostringstream savename;
    savename << _directory << _savenameStem << ".ply" << ends;
    _topolMesh.saveAsPLY( savename.str() );
    if( _saveAsOBJ ){
      ostringstream objname;
      objname << _directory << _savenameStem <<".obj" << ends;
      _topolMesh.saveAsOBJ( objname.str() );
    }
  }
  
  void saveFinalResult( string savename ) const {
    if( _topolMesh.getNumVertices() == 0 )
      throw aol::Exception ( "ParallelTransportOp::saveFinalResult(): topolMesh not set!", __FILE__, __LINE__ );
    _topolMesh.fromVector( *_transportedVariationShape );
    _topolMesh.saveAsPLY( savename );
  }
  
  void loadPath ( bool loadPath ) {        
    if( !loadPath )
      return initializePathLinearly();
    
    for( int k = 0; k < _pparser.getInt("numOfShells"); k++ ){
      // loadname
      ostringstream loadname;
      loadname << _directory << _pparser.getString("path");
      if(k<10) loadname <<"0";
      loadname << k <<".ply";      
      cerr << "Load " << loadname.str() << endl;

      // load
      MeshType mesh;
      mesh.loadFromPLY( loadname.str() );
      VectorType shell;
      mesh.toVector( shell );   
      for( int i = 0; i < 3; i++ )
        (*_path)[k*3+i] = shell[i]; 
    }
  }
  
  void initializePathLinearly(){
    if(!_quiet) cerr << "Initialize path linearly..." << endl;
    int numOfShells = _pparser.getInt("numOfShells");
    MeshType mesh( _pparser.getString("endShape") );
    VectorType startShape, difference;
    mesh.toVector( difference );
    mesh.loadFromPLY( _pparser.getString("startShape") );
    mesh.toVector( startShape );
    difference -= startShape;
    for( int i = 0; i < numOfShells; i++ )
      for( int j = 0; j < ConfiguratorType::Dim; j++ ){
	(*_path)[ ConfiguratorType::Dim*i + j ] = startShape[j];
	(*_path)[ ConfiguratorType::Dim*i + j ].addMultiple( difference[j], 1.*i / (numOfShells-1) );
      }
  }

 
};  
 

#endif