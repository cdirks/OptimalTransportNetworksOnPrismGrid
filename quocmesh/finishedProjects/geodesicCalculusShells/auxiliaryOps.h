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

#ifndef __AUXILIARYOPS_H
#define __AUXILIARYOPS_H

#include <Newton.h>
#include <trustRegionMethod.h>

#include <aol.h>
#include "geometry.h"

using namespace aol;
using namespace qc;
//using namespace om;

enum SolverType{ Cholesky, BiCholesky, LU, CG }; 
enum SolverMethod{ GRADIENTDESCENT, QUASINEWTON, NEWTON }; 

//!===========================================================================================================================
//!===========================================================================================================================


//!===========================================================================================================================
//!   UNIT CONFIGURATORS
//!   In most FE operators in the Quocmesh library, integration is done in some base class whereas the integrand is evaluated in derived classes.
//!   On the other hand, typical elastic deformation energies of shells (more precisely: elastic deformation densities) are 
//!     (1) constant on faces
//!     (2) require the evaluation of the determinant of the metric
//!   For the numerical quadrature performed in the base class the computation of the integration weights requires the computation of the area of
//!   a triangle, which is proportional to the squareroot of the determinant of the metric.
//!   To avoid computing this quantity twice (in the base class and in the integrand) we remove all geometric information from the base class
//!   and introduce a UnitTriangMeshConfigurator, which handles a (topological) triangular mesh which only consists of unit triangles (embedded in R^2).
//!   The integrand class now takes care of all geometric quantities.
//!===========================================================================================================================

//! Base function set for unit triangle with center quadrature rule.
//! Unit triangle embedded in R^2 is given by the three positions (0,0), (1,0) and (0,1).
template <typename RealType, typename TriangleType>
class UnitTriangMeshBaseFunctionSetCenterQuad 
  : public aol::BaseFunctionSetInterface<RealType, aol::Vec2<RealType>, aol::Vec2<RealType>, 3, aol::CenterQuadrature<RealType>, UnitTriangMeshBaseFunctionSetCenterQuad<RealType, TriangleType> >  {

  static RealType _b1   ( const aol::Vec2<RealType> &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const aol::Vec2<RealType> &c ) { return c[0]; }
  static RealType _b3   ( const aol::Vec2<RealType> &c ) { return c[1]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _basis[3];
  
  const aol::Vec3<RealType> _d1b, _d2b;

public:
  UnitTriangMeshBaseFunctionSetCenterQuad(  ) : _d1b(-1.,1.,0.), _d2b(-1.,0.,1.) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
  }
  
  enum { numBaseFuncs = 3 };

  void setTriangle ( const TriangleType &/*T*/ ) { }
  
  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &/*RefCoord*/, aol::Vec2<RealType> &Gradient ) const {
    Gradient[0] = _d1b[BaseFuncNum];
    Gradient[1] = _d2b[BaseFuncNum];
  }

  inline aol::Vec2<RealType> evaluateGradient ( int BaseFuncNum, int /*QuadPoint*/ ) const {
    return aol::Vec2<RealType>( _d1b[BaseFuncNum], _d2b[BaseFuncNum] );
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
};


//! Mesh configurator for topological triangular mesh, i.e. each triangle of the mesh is the unit triangle embedded in R^2,
//! given by the three positions (0,0), (1,0) and (0,1) with area 0.5
template <typename MeshType>
class UnitTriangMeshConfigurator {
protected:
  const MeshType &_mesh;
  
public:
  typedef typename MeshType::RealType     RealType;
  typedef aol::CenterQuadrature<RealType> QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                                 InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>                    DomVecType;
  typedef aol::Vec3<RealType>                      VecType;
  typedef aol::Mat<2, 2, RealType>                 MatType;
  typedef aol::Vector<RealType>                    VectorType;
  typedef aol::Vector<RealType>                    ArrayType;
  typedef aol::SparseMatrix<RealType>              MatrixType;
  typedef aol::BitVector                           MaskType;
  typedef typename MeshType::ElementType           ElementType;
  
  typedef UnitTriangMeshBaseFunctionSetCenterQuad<RealType,ElementType> BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;


  UnitTriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }

  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 3;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 3;
  }

  int getNumGlobalDofs( ) const {
    return this->_mesh.getNumVertices();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }
  
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return T.globNodeIdx( localIndex );
  }

  inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( T, localIndex0 );
    glob[1] = localToGlobal ( T, localIndex1 );
  }

  RealType vol ( const ElementType &/*T*/ ) const {
    return .5;
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }

};


//!===========================================================================================================================
//! FE OPERATOR INTERFACES
//!===========================================================================================================================
//!
//! Assumes center quadrature rule; only topology information is stored here; all geometric information is provided in derived classes
//! cf. aol::FENonlinVectorDiffOpInterface<> but we here do not need to evaluate volume of triangle (as it is constant!)
template <typename MeshType, typename Imp >
class UnitTriangleCenterQuadFEVectorInterface : public aol::Op< aol::MultiVector<typename MeshType::RealType> > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::ElementType ElementType;
  
  const MeshTopologySaver<MeshType> &_topology;
  const aol::Vec3<RealType> _d1b, _d2b;
  
public:
  explicit UnitTriangleCenterQuadFEVectorInterface( const MeshTopologySaver<MeshType> & Topology ) : _topology( Topology ), _d1b(-1.,1.,0.), _d2b(-1.,0.,1.) { }

  virtual ~UnitTriangleCenterQuadFEVectorInterface( ) {}

  void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    aol::Mat<3,2, RealType> ArgMatrix;
    for ( int faceIdx = 0; faceIdx < _topology.getNumFaces(); faceIdx++ ) {
      
      this->asImp().getNonlinearity ( Arg, faceIdx, ArgMatrix );
      
      for ( int dof = 0; dof < 3; dof++ ) {
        aol::Vec2<RealType> grad( _d1b[dof], _d2b[dof] );
        aol::Vec3<RealType> tmp; 
        ArgMatrix.mult ( grad, tmp );
	// note: no integration weights, area of unit triangle is 0.5!
        for ( int d = 0; d < 3; d++ )
          Dest[d][ _topology.getNodeOfTriangle ( faceIdx, dof ) ] += 0.5 * tmp[d];
      }
    }
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const MultiVector<RealType> &Arg,
                         const int &ElementIdx,
                         aol::Mat<3,2, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( Arg, ElementIdx, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

//! cf. aol::FELinAsymMatrixWeightedStiffInterface<>
//! provides an easy interface to Finite Element operators of the form \f$ \mbox{div}(A(x)\nabla u)\f$,
//! where \f$A\f$ is a ASYMMETRIC coefficient matrix.
template <typename MeshType, typename Imp, GridGlobalIndexMode IndexMode = QUOC_GRID_INDEX_MODE>
class UnitTriangleCenterQuadFEMatrixInterface :
      public aol::FELinOpInterface< typename MeshType::RealType, UnitTriangMeshConfigurator<MeshType>, UnitTriangleCenterQuadFEMatrixInterface<MeshType, Imp, IndexMode>, IndexMode > {

protected:
  typedef typename MeshType::RealType RealType;
  typedef UnitTriangMeshConfigurator<MeshType> ConfiguratorType;

  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::ElementType ElementType;
  
  aol::RandomAccessContainer< aol::Vec2<RealType> > _gradBasis;
  
public:     
  explicit UnitTriangleCenterQuadFEMatrixInterface ( const MeshTopologySaver<MeshType> & Topology, OperatorType OpType = ONTHEFLY )
      : aol::FELinOpInterface< RealType, ConfiguratorType, UnitTriangleCenterQuadFEMatrixInterface<MeshType, Imp, IndexMode>, IndexMode > ( Topology.getGrid(), OpType ){
	_gradBasis.pushBack( aol::Vec2<RealType>(-1., -1.) );
	_gradBasis.pushBack( aol::Vec2<RealType>( 1.,  0.) );
	_gradBasis.pushBack( aol::Vec2<RealType>( 0.,  1.) );
      }



  //! this function has to be provided in the implementation (derived class) of the interface
  inline void getCoeffMatrix ( const int& ElementIdx, MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( ElementIdx, Matrix );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    
    const int numDofs = this->getConfigurator().getNumLocalDofs ( El );
    LocalMatrix.setZero();

    MatType mat;
    aol::Vec2<RealType> matgrad1;
    getCoeffMatrix ( El.globIdx(), mat );
    
    // note: area of unit triangle is 0.5
    for ( int i = 0; i < numDofs; ++i ) {
      mat.mult ( _gradBasis[i], matgrad1 );
      for ( int j = 0; j < numDofs; ++j ) 
        LocalMatrix[j][i] += 0.5 * ( matgrad1 * _gradBasis[j] );
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};


//!===========================================================================================================================
//! CONSTRAINT HANDLER
//!===========================================================================================================================

//! \brief Base class for handling constraints
//! \author Heeren
//! Suppose we have a generic Lagrange function L[x,\lambda], where x = (x_j)_j and \lambda_i, i = 1, ..., N are the Lagrange multipliers.
//! We assume L[x,\lambda] = E[x] + \sum_i \lambda_i f_i[x], where f_i[x] = 0 represents the ith constraint.
//! This class provides for a fixed number N of constraints f_i
//!   (a) f_i[x] via evalConstraint()
//!   (b) \partial_{x_j} \sum_i \lambda_i f_i[x] via getGradOfConstraint()
//!   (c) \partial^2_{\lambda_k, x_j} \sum_i \lambda_i f_i[x] via getGradOfConstraintAndPositionInHessian( )
template <typename ConfiguratorType>
class ConstraintHandlerBase {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;
  
  const int _numOfConstraints;
  
public:
  ConstraintHandlerBase( const int NumOfConstraints ) : _numOfConstraints( NumOfConstraints ) {}
  
  virtual ~ConstraintHandlerBase () {}

  int getNumOfConstraints( ) const { return _numOfConstraints; }
  
  // returns f_i[x]
  virtual RealType evalConstraint( const VectorType& Arg, int i ) const = 0;
  
  // returns \sum_i f_i[x]
  RealType checkConstraints ( const VectorType& Arg ) const {
    RealType res = 0.;
    for( int i = 0; i < _numOfConstraints; i++ )
      res += evalConstraint( Arg, i );
    return res;
  }
  
  //
  virtual void getGradOfConstraint( const VectorType& Arg, int i, VectorType& Constraint, std::vector<int>& relevantBlocks, std::vector<RealType>& scaling ) const = 0;

  // derivative of the i-th constraint and its position in Hessian
  virtual void getGradOfConstraintAndPositionInHessian( const VectorType& Arg, int i, VectorType& Constraint, std::vector< aol::Vec2<int> >& relevantBlocks, std::vector<RealType>& scaling ) const = 0;

};

//! \brief Handles rigid body motion constraints for a single shell.
//! \author Heeren
//! There are usually 2 x Dim = 6 constraints, as we have to fix translations (one in each dimension) and rotations (one in each dimension).
template <typename ConfiguratorType>
class RigidBodyMotionsConstraintHandler : public ConstraintHandlerBase<ConfiguratorType> {
  
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;
  
  aol::Vector<RealType> _allOnes, _massOpAllOnes;
  VectorType _massOpIdentity;
  const VectorType& _identity;
  aol::MassOp<ConfiguratorType> *_massOp; 

public:
  RigidBodyMotionsConstraintHandler( const MeshTopologySaver<GridType>& topology, const VectorType& geometry, int NumOfConstraints ) 
    : ConstraintHandlerBase<ConfiguratorType>( NumOfConstraints ), 
    _allOnes( geometry[0].size(), 1. ),
    _massOpAllOnes( geometry[0].size() ),
    _massOpIdentity( 3, geometry[0].size() ),
    _identity( geometry )
    {
      _massOp = new aol::MassOp<ConfiguratorType>( topology.getGrid(), aol::ONTHEFLY );
      _massOp->apply( _allOnes, _massOpAllOnes );    
      for( int i = 0; i < 3; i++ )
	_massOp->apply( _identity[i], _massOpIdentity[i] );
    }
    
  ~RigidBodyMotionsConstraintHandler(){
    // delete mass op
    if( _massOp ) 
      delete _massOp; 
  }
 
  //
  RealType evalConstraint( const VectorType& Arg, int i ) const {
    VectorType displacement( Arg );
    displacement -= _identity;
    
    if( i < ConfiguratorType::Dim )
      return displacement[i] * _massOpAllOnes;
    else{
      int idx = i - ConfiguratorType::Dim;
      return _massOpIdentity[ (idx + 1)%3 ] * displacement[idx] - _massOpIdentity[idx] * displacement[ (idx + 1)%3 ];
    }    
  } 
  
  // derivative of the i-th constraint
  void getGradOfConstraint( const VectorType& /*Arg*/, int i, VectorType& Constraint, std::vector<int>& relevantBlocks, std::vector<RealType>& scaling  ) const {
    // (M*1[i], 0, 0)
    if( i < ConfiguratorType::Dim ){
      Constraint.clear();
      Constraint.appendReference( _massOpAllOnes );
      relevantBlocks.push_back( i );
      scaling.push_back( 1. );
    }
    // (M*\id[i+1], - M*\id[i], 0)
    else{
      int idx = i - ConfiguratorType::Dim;
      Constraint.clear();
      
      Constraint.appendReference( _massOpIdentity[ ( idx + 1 ) % 3 ] );
      relevantBlocks.push_back( idx );
      scaling.push_back( 1. );
      
      Constraint.appendReference( _massOpIdentity[ idx ] );
      relevantBlocks.push_back( ( idx + 1 ) % 3 );
      scaling.push_back( -1. );
    }    
  }

  // derivative of the i-th constraint and its position in Hessian
  void getGradOfConstraintAndPositionInHessian( const VectorType& Arg, int i, VectorType& Constraint, std::vector< aol::Vec2<int> >& relevantBlocks, std::vector<RealType>& scaling ) const {
    const int dim = ConfiguratorType::Dim;
    std::vector<int> relBlocks;
    getGradOfConstraint( Arg, i, Constraint, relBlocks, scaling );
    if( i < ConfiguratorType::Dim ){
      relevantBlocks.push_back( aol::Vec2<int>(i,i+dim) );
    }
    else{
      const int idx = i - dim;
      relevantBlocks.push_back( aol::Vec2<int>( idx, dim + i ) );
      relevantBlocks.push_back( aol::Vec2<int>( (idx+1)%3, dim + i ) );
    }
  }

};

//! \brief Handles rigid body motion constraints for a sequence of shells.
//! \author Heeren
//! There are usually 2 x Dim = 6 constraints per shell, as we have to fix translations (one in each dimension) and rotations (one in each dimension).
//! See documentation of class RigidBodyMotionsConstraintHandler<> above.
template <typename ConfiguratorType>
class MultipleRigidBodyMotionsConstraintHandler : public ConstraintHandlerBase<ConfiguratorType> {
  
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;

  const int _numOfGeometries, _numOfSingleConstraints;
  const RigidBodyMotionsConstraintHandler<ConfiguratorType> _rbmHandler;

public:
  MultipleRigidBodyMotionsConstraintHandler( const MeshTopologySaver<GridType>& topology, const VectorType& referenceObject, int NumOfSingleConstraints, int NumOfGeometries ) 
    : ConstraintHandlerBase<ConfiguratorType>( NumOfSingleConstraints * NumOfGeometries ), 
    _numOfGeometries( NumOfGeometries ), _numOfSingleConstraints( NumOfSingleConstraints ), _rbmHandler( topology, referenceObject, _numOfSingleConstraints ){  }

  //
  RealType evalConstraint( const VectorType& Arg, int i ) const {
    VectorType SingleArg;
    fillSingleArg( Arg, i, SingleArg );
    return _rbmHandler.evalConstraint( SingleArg, i%_numOfSingleConstraints );
  } 
  
  // derivative of the i-th constraint
  void getGradOfConstraint( const VectorType& Arg, int i, VectorType& Constraint, std::vector<int>& relevantBlocks, std::vector<RealType>& scaling  ) const {
    VectorType SingleArg;
    fillSingleArg( Arg, i, SingleArg );
    _rbmHandler.getGradOfConstraint( SingleArg, i%_numOfSingleConstraints, Constraint, relevantBlocks, scaling );
    // translate relevant blocks
    for( uint j = 0; j < relevantBlocks.size(); j++ )
      relevantBlocks[j] += (i / _numOfSingleConstraints) * ConfiguratorType::Dim;
  }

  // derivative of the i-th constraint and its position in Hessian
  void getGradOfConstraintAndPositionInHessian( const VectorType& Arg, int i, VectorType& Constraint, std::vector< aol::Vec2<int> >& relevantBlocks, std::vector<RealType>& scaling ) const {
    VectorType SingleArg;
    fillSingleArg( Arg, i, SingleArg );
    _rbmHandler.getGradOfConstraintAndPositionInHessian( SingleArg, i%_numOfSingleConstraints, Constraint, relevantBlocks, scaling );
    
    // translate relevant blocks
    for( uint j = 0; j < relevantBlocks.size(); j++ ){
      relevantBlocks[j][0] +=  ConfiguratorType::Dim * (i / _numOfSingleConstraints);
      relevantBlocks[j][1] += ConfiguratorType::Dim * (_numOfGeometries - 1) + (i / _numOfSingleConstraints) * _numOfSingleConstraints;
    }
  }

protected:
  void fillSingleArg( const VectorType& Arg, int i, VectorType& SingleArg ) const {
     int counter = i / _numOfSingleConstraints * ConfiguratorType::Dim;
     for( int j = 0; j < ConfiguratorType::Dim; j++ )
       SingleArg.appendReference( Arg[counter++] );
  }
};


//!===========================================================================================================================
//! LAGRANGE OPERATORS
//!===========================================================================================================================

//! \brief Allocate block matrix ( A B; B^T 0)
//! \author Heeren
template< typename TriMeshConfType, typename BlockMatrixType >
void allocateSystemMatrixForLagrange( BlockMatrixType &SystemMat, int numOfDOFs, int numOfGeometry, int numOfConstraints ) {
    
    if( (SystemMat.getNumRows() != numOfGeometry+numOfConstraints) || (SystemMat.getNumCols() != numOfGeometry+numOfConstraints) )
      throw aol::Exception ( "allocateSystemMatrixForLagrange: wrong number of blocks!!!", __FILE__, __LINE__ );

    // allocate matrix blocks
    for ( int i = 0; i < numOfGeometry; i++ )
      for ( int j = 0; j < numOfGeometry; j++ )
        if( !SystemMat.getPointer(i, j))
          SystemMat.allocateMatrix ( i, j, numOfDOFs, numOfDOFs );

    // allocate matrix blocks for contraints
    for ( int row = 0; row < numOfGeometry; row++ )
      for ( int col = numOfGeometry; col < numOfGeometry+numOfConstraints; col++ ){
        if( !SystemMat.getPointer(row, col) )
          SystemMat.allocateMatrix ( row, col, numOfDOFs, 1 );
        if( !SystemMat.getPointer(col, row) )
          SystemMat.allocateMatrix ( col, row, 1, numOfDOFs );
      }

    // allocate zero matrix blocks
    for ( int i = numOfGeometry; i < numOfGeometry+numOfConstraints; i++ )
      for ( int j = numOfGeometry; j < numOfGeometry+numOfConstraints; j++ )
        SystemMat.allocateMatrix ( i, j, 1, 1 );
}

//! Lagrange function L[x,\lambda] = E[x] + \sum_i \lambda_i F_i[x] for some 'Energy' E
//! where \lambda_i are Lagrange multipliers and F_i constraint functions provided by the 'ConstraintHandler'.
//! \author Heeren
template <typename ConfiguratorType, typename ConstraintHandlerType >
class GenericLagrangeFunction : public aol::Op< aol::MultiVector< typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;

  const Op< VectorType, aol::Scalar<RealType> >& _energy;
  const ConstraintHandlerType& _constraintHandler;
  
public:
  GenericLagrangeFunction ( const Op< VectorType, aol::Scalar<RealType> >& Energy, const ConstraintHandlerType &ConstraintHandler ) : _energy( Energy ), _constraintHandler( ConstraintHandler ){}

  // Arg is (x,y), where x contains the geometric information and y_1,...,y_m are the Lagrange multipliers
  void applyAdd( const VectorType& Arg, aol::Scalar<RealType>& Dest ) const {
        
    const int numOfConstraints = _constraintHandler.getNumOfConstraints();
    const int geometryComponents = Arg.numComponents() - numOfConstraints;
    
    // arg contains only geometric information x
    VectorType arg;
    for( int i = 0; i < geometryComponents; i++ )
      arg.appendReference( Arg[i] );

    _energy.applyAdd( arg, Dest );

    // now account for the y_i
    for( int i = 0; i < numOfConstraints; i++ )
      Dest[0] += Arg[geometryComponents+i][0] * _constraintHandler.evalConstraint( arg, i );

  }
  
};

//! Derivative of the Lagrange function L[x,\lambda] = E[x] + \sum_i \lambda_i F_i[x],
//! i.e. vector of all partial derivatives w.r.t. components of x = {x_j}_j and all \lambda_i
//! \author Heeren
template <typename ConfiguratorType, typename ConstraintHandlerType >
class GenericLagrangeGradient : public aol::Op< aol::MultiVector< typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;

  const Op< VectorType, VectorType >& _gradient;
  const ConstraintHandlerType& _constraintHandler;
  
public:
  GenericLagrangeGradient ( const Op< VectorType, VectorType >& Gradient, const ConstraintHandlerType &ConstraintHandler ) : _gradient( Gradient ), _constraintHandler( ConstraintHandler ){}

  // Arg is (x,y), where x contains the geometric information and y_1,...,y_m are the Lagrange multipliers
  void applyAdd( const VectorType& Arg, VectorType& Dest ) const {
    
    const int numOfConstraints = _constraintHandler.getNumOfConstraints();
    const int geometryComponents = Arg.numComponents() - numOfConstraints;
    
    if( Arg.numComponents() != Dest.numComponents() )
      throw aol::Exception ( "GenericLagrangeGradient::applyAdd: wrong number of blocks!!!", __FILE__, __LINE__ );
    
    // 
    VectorType arg, dest;
    for( int i = 0; i < geometryComponents; i++ ){
      arg.appendReference( Arg[i] );
      dest.appendReference( Dest[i] );
    }
    _gradient.applyAdd( arg, dest );
    
#ifdef _OPENMP
#pragma omp parallel for
#endif  
    for( int i = 0; i < numOfConstraints; i++ ){
      VectorType GradOfConstraint;
      std::vector<int> relevantBlocks;
      std::vector<RealType> scaling;
      _constraintHandler.getGradOfConstraint( arg, i, GradOfConstraint, relevantBlocks, scaling );

      for( unsigned int k = 0; k < relevantBlocks.size(); k++ )
        dest[ relevantBlocks[k] ].addMultiple( GradOfConstraint[k], scaling[k] * Arg[ geometryComponents + i ][0] );

      Dest[ geometryComponents + i ][0] += _constraintHandler.evalConstraint( arg, i );
    }

  }
  
};

//! Hessian of the Lagrange function L[x,\lambda] = E[x] + \sum_i \lambda_i F_i[x],
//! i.e. block matrix ( A B; B^T 0) where A is the Hessian of E[x]
//! \author Heeren
template <typename ConfiguratorType, typename ConstraintHandlerType, typename BlockMatrixType = aol::SparseBlockMatrix< aol::SparseMatrix<typename ConfiguratorType::RealType> > >
class GenericLagrangeHessian : public aol::Op< aol::MultiVector< typename ConfiguratorType::RealType>, BlockMatrixType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType GridType;
  typedef aol::MultiVector<RealType>          VectorType;

  const Op< VectorType, BlockMatrixType >& _hessian;
  const ConstraintHandlerType& _constraintHandler;
  
public:
  GenericLagrangeHessian ( const Op< VectorType, BlockMatrixType >& Hessian, const ConstraintHandlerType &ConstraintHandler ) : _hessian( Hessian ), _constraintHandler( ConstraintHandler ){}

  //
  void applyAdd( const VectorType& Arg, BlockMatrixType& Dest ) const {
    
    const int numOfConstraints = _constraintHandler.getNumOfConstraints();
    const int geometryComponents = Arg.numComponents() - numOfConstraints;
    
    if( Arg.numComponents() != Dest.getNumCols() )
      throw aol::Exception ( "GenericLagrangeHessian::applyAdd: wrong number of cols!!!", __FILE__, __LINE__ );
    if( Arg.numComponents() != Dest.getNumRows() )
      throw aol::Exception ( "GenericLagrangeHessian::applyAdd: wrong number of rows!!!", __FILE__, __LINE__ );

    int numOfDOFs = Arg[0].size();   
    allocateSystemMatrixForLagrange<ConfiguratorType, BlockMatrixType>( Dest, numOfDOFs, geometryComponents, numOfConstraints );
   
    VectorType arg;
    for( int i = 0; i < geometryComponents; i++ )
      arg.appendReference( Arg[i] );

    // get reference on subblock of Dest
    std::vector< std::vector< std::vector< int > > > blockPositions( geometryComponents );    
    for ( int i = 0; i < geometryComponents; ++i ){
      blockPositions[i].resize ( geometryComponents );
      for ( int j = 0; j < geometryComponents; ++j ) {
	blockPositions[i][j].resize ( 2 );
        blockPositions[i][j][0] = i;
        blockPositions[i][j][1] = j;
      }
    }
  
    BlockMatrixType dest(geometryComponents, geometryComponents);
    Dest.getSubBlockMatrix(geometryComponents, geometryComponents, blockPositions, dest );    
    _hessian.applyAdd( arg, dest );
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int c = 0; c < numOfConstraints; c++ ){
      VectorType GradOfConstraint;
      std::vector< aol::Vec2<int> > relevantBlocks;
      std::vector<RealType> scaling;
      _constraintHandler.getGradOfConstraintAndPositionInHessian( arg, c, GradOfConstraint, relevantBlocks, scaling );
      
#ifdef _OPENMP
#pragma omp parallel for
#endif    
      for( unsigned int k = 0; k < relevantBlocks.size(); k++ )
        for( int j = 0; j < numOfDOFs; j++ ){
          Dest.getReference( relevantBlocks[k][0], relevantBlocks[k][1] ).add( j, 0, scaling[k] * GradOfConstraint[k][j] );
          Dest.getReference( relevantBlocks[k][1], relevantBlocks[k][0] ).add( 0, j, scaling[k] * GradOfConstraint[k][j] );
        }  
    }

  }
  
};

//!===========================================================================================================================
//! SOLVER
//!===========================================================================================================================

//! Newton Iteration written by Benedikt Wirth
//!TODO use aol::NewtonIterationSparseBlockMatrix<aol::SparseMatrix<>, aol::CholeskyBiBlockInverseOp<ConfiguratorType, aol::SparseMatrix<>, ConfiguratorType::Dim>, ConfiguratorType::Dim > instead!?!
template <typename ConfiguratorType, typename VectorType, typename SecondDerivativeType>
class NewtonMethod :
  public aol::NewtonIterationBase<typename ConfiguratorType::RealType,VectorType,VectorType,SecondDerivativeType> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename SecondDerivativeType::MatrixType SubMatrixType; 
  typedef typename ConfiguratorType::InitType GridType;
  
  typedef aol::IdentityOp<VectorType> PreconditionType; 
  mutable PreconditionType *_pPrecond;
  
  //const typename ConfiguratorType::InitType &_grid;
  const int _solverType;

  virtual void prepareSolver () const {
    if( this->_pSolver )
      delete this->_pSolver;
 
    if( _solverType == BiCholesky ){
      throw aol::Exception ( "NewtonMethod: BiCholesky not available!!!", __FILE__, __LINE__ );
/*           
      //cerr << "maxNumNonZeroesPerRow = " << this->_pMatDF->maxNumNonZeroesPerRow() << endl;
      //! Input A, b, where A is supposed to be non-singular; solves (A^TA)x = A^Tb for x by Cholesky factorization of A^TA
      aol::CholeskyBiBlockInverseOp<ConfiguratorType,SubMatrixType,_numOfBLocks>* pSolver = new aol::CholeskyBiBlockInverseOp<ConfiguratorType,SubMatrixType,_numOfBLocks>( _grid, *(this->_pMatDF) );
      this->_pSolver = pSolver;
*/        
    }
 
    if( _solverType == Cholesky ){
      //cerr << "maxNumNonZeroesPerRow = " << this->_pMatDF->maxNumNonZeroesPerRow() << endl;
      if ( this->_pSolver == NULL )
        this->_pSolver = new aol::CholeskyBlockInverseOp<RealType,SubMatrixType>( this->_pMatDF->maxNumNonZeroesPerRow() );
      static_cast<aol::CholeskyBlockInverseOp<RealType,SubMatrixType>*>( this->_pSolver )->setMatrix( *(this->_pMatDF) );
    }
    
    if( _solverType == LU ){
      //! Input as above, computs LU-factorization of A and solves Ax=b for x
      aol::UMFPACKBlockInverseOp<RealType,SubMatrixType>* pSolver = new aol::UMFPACKBlockInverseOp<RealType,SubMatrixType>( *(this->_pMatDF) );    
      this->_pSolver = pSolver;
    }
    
    if( _solverType == CG ){
      _pPrecond = new PreconditionType( *(this->_pMatDF) );
      this->_pSolver = new aol::PCGInverse<VectorType>( *(this->_pMatDF), *_pPrecond, this->_pInfo->getSolverInfo() );
      // the solver-accuracy is automatically set to fit to the Newton-minimizer accuracy
      static_cast<aol::PCGInverse<VectorType>*>( this->_pSolver )->setStopping( aol::STOPPING_ABSOLUTE );
    }
  }
  
  RealType computeErrorNorm ( const VectorType& /*x_k*/, const VectorType& delta_x_k, const VectorType& /*F_x_k*/, bool initial ) const {
    if( initial )
      return 1.;
    else
      return delta_x_k.norm();
  }
  
public:
  NewtonMethod( const typename ConfiguratorType::InitType &Grid,
                const aol::Op<VectorType> &F,
                const aol::Op<VectorType, SecondDerivativeType> &DF,
		const int NumOfGeometryBlocks,
		const int NumOfConstraintBlocks,
		const int SolverType,
                const int MaxIterations = 1000,
                const RealType StopEpsilon = 1e-8 ) :
    aol::NewtonIterationBase<RealType,VectorType,VectorType,SecondDerivativeType>( new SecondDerivativeType( NumOfGeometryBlocks+NumOfConstraintBlocks, NumOfGeometryBlocks+NumOfConstraintBlocks ), F, DF, MaxIterations, StopEpsilon, false, NULL ), _solverType( SolverType ){
      
    // check  
    if( (_solverType < 0) || (_solverType > 3) )
      throw aol::Exception ( "NewtonMethod: unknown solver!!!", __FILE__, __LINE__ );
    
    allocateSystemMatrixForLagrange<ConfiguratorType, SecondDerivativeType>( *(this->_pMatDF), Grid.getNumVertices(), NumOfGeometryBlocks, NumOfConstraintBlocks );
  }
  
    NewtonMethod( const MeshTopologySaver<GridType> &Topology,
                const aol::Op<VectorType> &F,
                const aol::Op<VectorType, SecondDerivativeType> &DF,
		const int NumOfGeometryBlocks,
		const int NumOfConstraintBlocks,
		const int SolverType,
                const int MaxIterations = 1000,
                const RealType StopEpsilon = 1e-8 ) :
    aol::NewtonIterationBase<RealType,VectorType,VectorType,SecondDerivativeType>( new SecondDerivativeType( NumOfGeometryBlocks+NumOfConstraintBlocks, NumOfGeometryBlocks+NumOfConstraintBlocks ), F, DF, MaxIterations, StopEpsilon, false, NULL ), _solverType( SolverType ){
      
    // check  
    if( (_solverType < 0) || (_solverType > 3) )
      throw aol::Exception ( "NewtonMethod: unknown solver!!!", __FILE__, __LINE__ );
    
    allocateSystemMatrixForLagrange<ConfiguratorType, SecondDerivativeType>( *(this->_pMatDF), Topology.getNumVertices(), NumOfGeometryBlocks, NumOfConstraintBlocks );
  }
  
  void setNumDerivativeHoldingSteps( int num ) const {
    this->_pInfo->setNumDerivativeHoldingSteps( num );
  }
  
  virtual ~NewtonMethod(){  }
};

//! NewtonMethodSolver
template <typename ConfiguratorType, typename SecondDerivativeType = aol::SparseBlockMatrix< aol::SparseMatrix<typename ConfiguratorType::RealType> > >
class NewtonMethodSolver {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType MeshType;
  typedef typename SecondDerivativeType::MatrixType SubMatrixType; 
  typedef aol::Scalar<RealType> ScalarType;
  typedef aol::MultiVector<RealType> VectorType;
  typedef typename aol::NewtonInfo<RealType>::TIMESTEP_CONTROLLER TimestepControllerType;
  
  //const typename ConfiguratorType::InitType &_grid;
  const MeshTopologySaver<MeshType>* _topology;
  const VectorType* _position;
  bool _deleteTopSaver, _deletePosition;
  
  //const aol::Op<VectorType, ScalarType > &_F;
  const aol::Op<VectorType> &_DF;
  const aol::Op<VectorType, SecondDerivativeType> &_D2F;

  const int _solverType, _maxIterations;
  TimestepControllerType _timestepController;
  int _numDerivativeHoldingSteps;
  mutable RealType _epsilon, _sigma;
  bool _quiet;

  
public:
  NewtonMethodSolver( const MeshType &Grid,
                      const aol::Op<VectorType> &DF,
                      const aol::Op<VectorType, SecondDerivativeType> &D2F,
		      int SolverType,
                      int MaxIterations = 1000,
		      RealType StopEpsilon = 1.e-8,
		      bool Quiet = true,
		      TimestepControllerType TimestepController = aol::NewtonInfo<RealType>::NEWTON_OPTIMAL,                      
		      RealType Sigma = 0.1,
		      int numDerivativeHoldingSteps = 1 ) : _DF( DF ), _D2F( D2F ), _solverType(SolverType), _maxIterations( MaxIterations ), _timestepController( TimestepController ), _numDerivativeHoldingSteps( numDerivativeHoldingSteps ), _epsilon( StopEpsilon ), _sigma( Sigma ), _quiet( Quiet ) {
			_topology = new MeshTopologySaver<MeshType>( Grid );
			_deleteTopSaver = true;
			// set position to be fixed in Lagrange setup
			VectorType tempPosition( 3, Grid.getNumVertices() );
			Grid.toVector( tempPosition );
			_position = new VectorType( tempPosition );
			_deletePosition = true;
		      }
		      
  NewtonMethodSolver( const MeshTopologySaver<MeshType>& Topology,
                      const aol::Op<VectorType> &DF,
                      const aol::Op<VectorType, SecondDerivativeType> &D2F,
		      int SolverType,
                      int MaxIterations = 1000,
		      RealType StopEpsilon = 1.e-8,
		      bool Quiet = true,
		      TimestepControllerType TimestepController = aol::NewtonInfo<RealType>::NEWTON_OPTIMAL,                      
		      RealType Sigma = 0.1,
		      int numDerivativeHoldingSteps = 1 ) : _topology( &Topology ), _deleteTopSaver( false ), _deletePosition( false ),  _DF( DF ), _D2F( D2F ), _solverType(SolverType), _maxIterations( MaxIterations ), _timestepController( TimestepController ), _numDerivativeHoldingSteps( numDerivativeHoldingSteps ), _epsilon( StopEpsilon ), _sigma( Sigma ), _quiet( Quiet ) {}
		      
  NewtonMethodSolver( const MeshTopologySaver<MeshType>& Topology,
		      const VectorType& Position, 
                      const aol::Op<VectorType> &DF,
                      const aol::Op<VectorType, SecondDerivativeType> &D2F,
		      int SolverType,
                      int MaxIterations = 1000,
		      RealType StopEpsilon = 1.e-8,
		      bool Quiet = true,
		      TimestepControllerType TimestepController = aol::NewtonInfo<RealType>::NEWTON_OPTIMAL,                      
		      RealType Sigma = 0.1,
		      int numDerivativeHoldingSteps = 1  ) : _topology( &Topology ), _position( &Position ), _deleteTopSaver( false ), _deletePosition( false ),  _DF( DF ), _D2F( D2F ), _solverType(SolverType), _maxIterations( MaxIterations ), _timestepController( TimestepController ), _numDerivativeHoldingSteps( numDerivativeHoldingSteps ), _epsilon( StopEpsilon ), _sigma( Sigma ), _quiet( Quiet ) {}
                      
                      
  ~NewtonMethodSolver(){
    if( _deleteTopSaver )
      delete _topology;
    if( _deletePosition )
      delete _position;
  }
  
  void solve( const VectorType& InitialGuess, VectorType& Solution, bool Lagrange = false ) const { 
    
    if( _numDerivativeHoldingSteps > 1 ) cerr << "Caution!  NumDerivativeHoldingSteps is set to " <<  _numDerivativeHoldingSteps << endl;
    
    // update stopping criterion
    RealType newEps = 1e-8 * sqrt( InitialGuess.numComponents() * InitialGuess[0].size() );
    if( newEps > _epsilon ){
      _epsilon = newEps;
      if( !_quiet ) cerr << "Stopping criterion set to " <<  _epsilon << endl;
    }
    
    if( Lagrange ){
      
      if( !_position )
	throw aol::Exception ( "NewtonMethodSolver::solve(): position not set!", __FILE__, __LINE__ );
      
      const int NumOfConstraints = 6;
      //RigidBodyMotionsConstraintHandler<ConfiguratorType> CHandler( *_topology, *_position, NumOfConstraints );
      
      typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
      RBMHandlerType CHandler( *_topology, *_position, NumOfConstraints, 1 );
    
      //GenericLagrangeFunction< ConfiguratorType, RBMHandlerType >  L( _F, CHandler );
      GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( _DF, CHandler );
      GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( _D2F, CHandler );
  
      VectorType Arg( InitialGuess, aol::FLAT_COPY ), Sol( Solution, aol::FLAT_COPY );
      VectorType lagrangeArg( NumOfConstraints, 1 ), lagrangeSol( NumOfConstraints, 1 );    
      Arg.appendReference( lagrangeArg );
      Sol.appendReference( lagrangeSol );

      NewtonMethod<ConfiguratorType, VectorType, SecondDerivativeType> zeroFinder( *_topology, dL, dL2, ConfiguratorType::Dim, NumOfConstraints, _solverType, _maxIterations, _epsilon );
      zeroFinder.setTimestepController( _timestepController );
      if( _numDerivativeHoldingSteps > 1 ) zeroFinder.setNumDerivativeHoldingSteps( _numDerivativeHoldingSteps );
      zeroFinder.setSigma( _sigma ); 
      zeroFinder.setQuietMode( _quiet );
      zeroFinder.apply( Arg, Sol );
            
      if( !_quiet ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Solution ) << endl;
    }
    else
    {
      NewtonMethod<ConfiguratorType, VectorType, SecondDerivativeType> zeroFinder( *_topology, _DF, _D2F, ConfiguratorType::Dim, 0, _solverType, _maxIterations, _epsilon );
      zeroFinder.setTimestepController( _timestepController );
      if( _numDerivativeHoldingSteps > 1 ) zeroFinder.setNumDerivativeHoldingSteps( _numDerivativeHoldingSteps );
      zeroFinder.setSigma( _sigma ); 
      zeroFinder.setQuietMode( _quiet );
      zeroFinder.apply( InitialGuess, Solution );
    }    
  }
  
  void setQuietMode( bool Quiet ) {
    _quiet = Quiet;
  }
};

#endif