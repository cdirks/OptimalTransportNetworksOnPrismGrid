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

#ifndef __FELEVELSETOPINTERFACE_H
#define __FELEVELSETOPINTERFACE_H

#include <tpCFEGrid.h>
#include <tpCFEUtils.h>
#include <scalarArray.h>
#include <discreteFunction.h>
#include <gridBase.h>
#include <isolineIterator2d.h>
#include <configurators.h>


namespace qc {

/**
 * \brief This class provides an interface for general linear operators on functions that live on levelsets (given as aol::Vector on a grid).
 * It provides assembly routines for the corresponding matrices and apply routines for direct application to a aol::Vector.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, Dimension Dim, typename Imp>
class FELevelsetOpInterface {};

/**
 * \brief This class provides an interface for general linear operators on functions that live on levelsets (given as aol::Vector on a grid).
 * It provides assembly routines for the corresponding matrices and apply routines for direct application to a aol::Vector.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename Imp>
class FELevelsetOpInterface<ConfiguratorType, QC_2D, Imp> :
      public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const ConfiguratorType _config;
  const typename ConfiguratorType::ArrayType &_levelsetFunction;
  mutable typename ConfiguratorType::MatrixType *_mat;
  const aol::OperatorType _opType;

public:
  /**
   * The grid and the levelset function together define the zero levelset, where the operands of the operator live.
   */
  explicit FELevelsetOpInterface ( const typename ConfiguratorType::InitType &Grid, const typename ConfiguratorType::ArrayType &LevelsetFunction, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      _config ( Grid ),
      _levelsetFunction ( LevelsetFunction ),
      _mat ( NULL ),
      _opType ( OpType ) {}

  virtual ~FELevelsetOpInterface( ) {
    delete _mat;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    switch ( _opType ) {
      case aol::ONTHEFLY:
        multiplyOnTheFly ( Arg, Dest );
        break;
      case aol::ASSEMBLED:
        if ( !_mat )
          assembleMatrix( );
        _mat->applyAdd ( Arg, Dest );
        break;
      default:
        throw aol::UnimplementedCodeException ( "FELevelsetOpInterface::applyAdd: unsupported opType", __FILE__, __LINE__ );
    };
  }

  // clears the assembled matrix
  void reset( ) {
    if ( _mat )
      delete _mat;
    _mat = NULL;
  }

  typename ConfiguratorType::MatrixType& getMatrix( ) {
    if ( !_mat ) {
      _mat = _config.createNewMatrix( );
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

protected:
  void multiplyOnTheFly ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // will contain the global indices of the DOFs
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    // will contain the local system matrix for one single element
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    // traverse all lines on the interface
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _config.getInitializer(), _levelsetFunction );
    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( isoIt, localMatrix );

      const int numLocalDofs = _config.getNumLocalDofs ( isoIt->el );

      // get the global indices to the current Dofs
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[i] = _config.localToGlobal ( isoIt->el, i );

      // finally add the locally computed values to the result
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Dest[globalDofs[i]] += localMatrix [ i ][ j ] * Arg[globalDofs[j]] ;
    }
  }

  void assembleMatrix( ) const {
    if ( _mat )
      delete _mat;
    _mat = _config.createNewMatrix( );
    assembleAddMatrix ( *_mat );
  }

public:
  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {
    // will contain the global indices of the DOFs
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    // will contain the local system matrix for one single element
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    // traverse all lines on the interface
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _config.getInitializer(), _levelsetFunction );
    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {
      // assemble the local matrix for the current element
      this->asImp().prepareLocalMatrix ( isoIt, localMatrix );

      const int numLocalDofs = _config.getNumLocalDofs ( isoIt->el );

      // get the global indices of the local Dofs of the current element
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[i] = _config.localToGlobal ( isoIt->el, i );

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Mat.add ( globalDofs[i], globalDofs[j], Factor * localMatrix [ i ][ j ] );
    }
  }

protected:
  // barton-nackman
inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief This class provides an interface for general linear operators on functions that live on levelsets (given as aol::Vector on a grid).
 * It provides assembly routines for the corresponding matrices and apply routines for direct application to a aol::Vector.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename Imp>
class FELevelsetOpInterface<ConfiguratorType, QC_3D, Imp> :
      public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const ConfiguratorType _config;
  const typename ConfiguratorType::ArrayType &_levelsetFunction;
  tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> _cfeGrid;
  mutable typename ConfiguratorType::MatrixType *_mat;
  const aol::OperatorType _opType;

public:
  /**
   * The grid and the levelset function together define the zero levelset, where the operands of the operator live.
   */
  explicit FELevelsetOpInterface ( const typename ConfiguratorType::InitType &Grid, const typename ConfiguratorType::ArrayType &LevelsetFunction, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      _config ( Grid ),
      _levelsetFunction ( LevelsetFunction ),
      _cfeGrid ( Grid.getGridDepth() ),
      _mat ( NULL ),
      _opType ( OpType ) {
    _cfeGrid.setAdjustLevelset ( 1.0e-10 );
    _cfeGrid.setDomainFrom ( LevelsetFunction );
    _cfeGrid.detectVirtualNodes ();
  }

  virtual ~FELevelsetOpInterface( ) {
    delete _mat;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    switch ( _opType ) {
      case aol::ONTHEFLY:
        multiplyOnTheFly ( Arg, Dest );
        break;
      case aol::ASSEMBLED:
        if ( !_mat )
          assembleMatrix( );
        _mat->applyAdd ( Arg, Dest );
        break;
      default:
        throw aol::UnimplementedCodeException ( "FELevelsetOpInterface::applyAdd: unsupported opType", __FILE__, __LINE__ );
    };
  }

  // clears the assembled matrix
  void reset( ) {
    if ( _mat )
      delete _mat;
    _mat = NULL;
  }

  typename ConfiguratorType::MatrixType& getMatrix() {
    if ( !_mat ) {
      _mat = _config.createNewMatrix( );
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

protected:
  void multiplyOnTheFly ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // will contain the global indices of the DOFs
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    // will contain the local system matrix for one single element
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    // traverse all triangles on the interface
    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( itit, localMatrix );

      typename ConfiguratorType::ElementType el ( itit.getTriangRef().getElement() );
      const int numLocalDofs = _config.getNumLocalDofs ( el );

      // get the global indices to the current Dofs
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[i] = _config.localToGlobal ( el, i );

      // finally add the locally computed values to the result
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Dest[globalDofs[i]] += localMatrix [ i ][ j ] * Arg[globalDofs[j]] ;
    }
  }

  void assembleMatrix() const {
    if ( _mat )
      delete _mat;
    _mat = _config.createNewMatrix();
    assembleAddMatrix ( *_mat );
  }

public:
  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {
    // will contain the global indices of the DOFs
    int globalDofs[ConfiguratorType::maxNumLocalDofs];
    // will contain the local system matrix for one single element
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;

    // traverse all triangles on the interface
    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( itit, localMatrix );

      const typename ConfiguratorType::ElementType el ( itit.getTriangRef().getElement() );
      const int numLocalDofs = _config.getNumLocalDofs ( el );

      // get the global indices to the current Dofs
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[i] = _config.localToGlobal ( el, i );

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Mat.add ( globalDofs[i], globalDofs[j], Factor * localMatrix [ i ][ j ] );
    }
  }

protected:
  // barton-nackman
inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief This class implements the operator (expressed in matrix form)
 * \f$\left(\int_\Gamma w(x)\varphi_i(x)\varphi_j(x) dx\right)_{ij}\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\varphi_k\f$ is the \f$k\f$th
 * finite element basis function, and the weight \f$w\f$ has to be defined in the derived classes as "getCoeff".
 *
 * Currently, only midpoint quadrature on the line segments of \f$\Gamma\f$ is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, Dimension Dim, typename Imp>
class FELinLevelsetScalarWeightedMassInterface {};

/**
 * \brief This class implements the operator (expressed in matrix form)
 * \f$\left(\int_\Gamma w(x)\varphi_i(x)\varphi_j(x) dx\right)_{ij}\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\varphi_k\f$ is the \f$k\f$th
 * finite element basis function, and the weight \f$w\f$ has to be defined in the derived classes as "getCoeff".
 *
 * Currently, only midpoint quadrature on the line segments of \f$\Gamma\f$ is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename Imp>
class FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_2D, Imp> :
      public FELevelsetOpInterface<ConfiguratorType, QC_2D, FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_2D, Imp> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit FELinLevelsetScalarWeightedMassInterface ( const typename ConfiguratorType::InitType &Grid, const typename ConfiguratorType::ArrayType &LevelsetFunction, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      FELevelsetOpInterface<ConfiguratorType, QC_2D, FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_2D, Imp> > ( Grid, LevelsetFunction, OpType ) {}

  // this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, const DomVecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getCoeff ( El, RefCoord );
  }

  // this function computes the numerical quadrature of the bilinear form and saves the values locally
  void prepareLocalMatrix ( qc::IsoLineIterator2d<ConfiguratorType> &IsoIt, aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> &LocalMatrix ) const {
    const int numDofs = this->_config.getNumLocalDofs ( IsoIt->el );
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config.getBaseFunctionSet ( IsoIt->el );

    // compute center of mass for line (in local coordinates, not scaled by h)
    typename ConfiguratorType::VecType lineNodes[2];
    for ( int i = 0; i < 2; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        lineNodes[i][j] = IsoIt->points[i][j] - IsoIt->el[j];
    typename ConfiguratorType::VecType centerOfMass ( lineNodes[0] + lineNodes[1] );
    centerOfMass /= static_cast<RealType> ( 2 );

    // compute values needed for midpoint quadrature
    const RealType length = sqrt ( aol::Sqr ( IsoIt->points[0][0] - IsoIt->points[1][0] ) + aol::Sqr ( IsoIt->points[0][1] - IsoIt->points[1][1] ) ) * this->_config.getInitializer().H();
    RealType coeff = this->asImp().getCoeff ( IsoIt->el, centerOfMass );

    // compute upper half of local matrix for current element
    RealType basisVal[this->_config.maxNumLocalDofs];
    for ( int i = 0; i < numDofs; i++ )
      basisVal[i] = bfs.evaluate ( i, centerOfMass );
    for ( int i = 0; i < numDofs; i++ )
      for ( int j = i; j < numDofs; j++ )
        LocalMatrix[i][j] =  basisVal[i] * basisVal[j] * coeff;

    // compute lower half and scale the local matrix
    for ( int i = 0; i < numDofs; i++ )
      for ( int j = i; j < numDofs; j++ )
        LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * length;
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief This class implements the operator (expressed in matrix form)
 * \f$\left(\int_\Gamma w(x)\varphi_i(x)\varphi_j(x) dx\right)_{ij}\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\varphi_k\f$ is the \f$k\f$th
 * finite element basis function, and the weight \f$w\f$ has to be defined in the derived classes as "getCoeff".
 *
 * Currently, only midpoint quadrature on the line segments of \f$\Gamma\f$ is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename Imp>
class FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_3D, Imp> :
      public FELevelsetOpInterface<ConfiguratorType, QC_3D, FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_3D, Imp> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit FELinLevelsetScalarWeightedMassInterface ( const typename ConfiguratorType::InitType &Grid, const typename ConfiguratorType::ArrayType &LevelsetFunction, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      FELevelsetOpInterface<ConfiguratorType, QC_3D, FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, QC_3D, Imp> > ( Grid, LevelsetFunction, OpType ) {}

  // this function has to be provided in the implementation (derived class) of the interface
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, const DomVecType& RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().getCoeff ( El, RefCoord );
  }

  // this function computes the numerical quadrature of the bilinear form and saves the values locally
  void prepareLocalMatrix ( const tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > &ItIt, aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> &LocalMatrix ) const {
    const typename ConfiguratorType::ElementType el ( ItIt.getTriangRef().getElement() );
    const int numDofs = this->_config.getNumLocalDofs ( el );
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config.getBaseFunctionSet ( el );

    // compute center of mass for triangle (in local coordinates, not scaled by h)
    std::vector< aol::Vec3< RealType > > triangleNodes = ItIt.getTriangRef().getLocalCoordinates();
    typename ConfiguratorType::VecType centerOfMass ( triangleNodes[0] + triangleNodes[1] + triangleNodes[2] );
    centerOfMass /= static_cast<RealType> ( 3 );

    // compute values needed for midpoint quadrature
    const RealType area = ItIt.getTriangRef().getArea ( this->_config.getInitializer().H() );
    RealType coeff = this->asImp().getCoeff ( el, centerOfMass );

    // compute upper half of local matrix for current element
    RealType basisVal[this->_config.maxNumLocalDofs];
    for ( int i = 0; i < numDofs; i++ )
      basisVal[i] = bfs.evaluate ( i, centerOfMass );
    for ( int i = 0; i < numDofs; i++ )
      for ( int j = i; j < numDofs; j++ )
        LocalMatrix[i][j] =  basisVal[i] * basisVal[j] * coeff;

    // compute lower half and scale the local matrix
    for ( int i = 0; i < numDofs; i++ )
      for ( int j = i; j < numDofs; j++ )
        LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * area;
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma \vec{f}\left(\vec\phi(x),\nabla\vec\phi(x),x\right)\cdot\vec\varphi_i(x) dx\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec f\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename qc::Dimension Dim, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorOpInterface {};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma \vec{f}\left(\vec\phi(x),\nabla\vec\phi(x),x\right)\cdot\vec\varphi_i(x) dx\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec f\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorOpInterface<ConfiguratorType, qc::QC_2D, NumCompArg, NumCompDest, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const ConfiguratorType _config;
  const InitType &_grid;
  const ArrayType &_levelSetFunction;

public:
  FENonlinLevelsetVectorOpInterface ( const InitType &Grid,
                                      const ArrayType &LevelSetFunction ) :
      _config ( Grid ),
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // check, whether Arg and Dest have correct size
    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << Dest.numComponents() << " " << NumCompDest << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "nl" (it will contain the nonlinear function value at each quadrature point) as well as "val" (quadrature result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Vec<NumCompDest, RealType> nl, val;
    const RealType h  = _grid.H();

    // traverse all lines on the interface
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _grid, _levelSetFunction );
    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {
      // compute center of mass for line (in local coordinates, not scaled by h)
      typename ConfiguratorType::VecType lineNodes[2];
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          lineNodes[i][j] = isoIt->points[i][j] - isoIt->el[j];
      typename ConfiguratorType::VecType centerOfMass ( lineNodes[0] + lineNodes[1] );
      centerOfMass /= static_cast<RealType> ( 2 );

      // compute values needed for midpoint quadrature
      const qc::Element el ( isoIt->el );
      const RealType length = sqrt ( aol::Sqr ( isoIt->points[0][0] - isoIt->points[1][0] ) + aol::Sqr ( isoIt->points[0][1] - isoIt->points[1][1] ) ) * h;
      asImp().getNonlinearity ( discrFuncsArg, el, centerOfMass, nl );

      // for each hat function compute its product with the nonlinear function and the midpoint quadrature of that
      const int numLocalDofs = _config.getNumLocalDofs ( el );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( el );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        val = nl * ( length * bfs.evaluate ( dof, centerOfMass ) );

        for ( int d = 0; d < NumCompDest; d++ )
          Dest[d][ _config.localToGlobal ( el, dof ) ] += val[d];
      }
    }
  }

  // Interface function; has to be provided in derived classes.
  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Vec<NumCompDest, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, NL );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma \vec{f}\left(\vec\phi(x),\nabla\vec\phi(x),x\right)\cdot\vec\varphi_i(x) dx\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec f\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorOpInterface<ConfiguratorType, qc::QC_3D, NumCompArg, NumCompDest, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const ConfiguratorType _config;
  const InitType &_grid;
  const ArrayType &_levelSetFunction;
  tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> _cfeGrid;

public:
  FENonlinLevelsetVectorOpInterface ( const InitType &Grid,
                                      const ArrayType &LevelSetFunction ) :
      _config ( Grid ),
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ),
      _cfeGrid ( _grid.getGridDepth() ) {
    _cfeGrid.setAdjustLevelset ( 1.0e-10 );
    _cfeGrid.setDomainFrom ( LevelSetFunction );
    _cfeGrid.detectVirtualNodes ();
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // check, whether Arg and Dest have correct size
    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << Dest.numComponents() << " " << NumCompDest << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "nl" (it will contain the nonlinear function value at each quadrature point) as well as "val" (quadrature result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Vec<NumCompDest, RealType> nl, val;
    const RealType h  = _cfeGrid.H();

    // traverse all triangles on the interface
    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      // compute center of mass for triangle (in local coordinates, not scaled by h)
      std::vector< aol::Vec3< RealType > > triangleNodes = itit.getTriangRef().getLocalCoordinates();
      typename ConfiguratorType::VecType centerOfMass ( triangleNodes[0] + triangleNodes[1] + triangleNodes[2] );
      centerOfMass /= static_cast<RealType> ( 3 );

      // compute values needed for midpoint quadrature
      const qc::Element el ( itit.getTriangRef().getElement() );
      const RealType area = itit.getTriangRef().getArea ( h );
      asImp().getNonlinearity ( discrFuncsArg, el, centerOfMass, nl );

      // for each hat function compute its product with the nonlinear function and the midpoint quadrature of that
      const int numLocalDofs = _config.getNumLocalDofs ( el );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( el );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        val = nl * ( area * bfs.evaluate ( dof, centerOfMass ) );

        for ( int d = 0; d < NumCompDest; d++ )
          Dest[d][ _config.localToGlobal ( el, dof ) ] += val[d];
      }
    }
  }

  // Interface function; has to be provided in derived classes.
  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Vec<NumCompDest, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, NL );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma\vec{\vec{f}}\left(\vec\phi(x),\nabla\vec\phi(x),x\right):\nabla \vec\varphi_i(x) dx\!\!\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec{\vec{f}}\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename qc::Dimension Dim, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorDiffOpInterface {};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma\vec{\vec{f}}\left(\vec\phi(x),\nabla\vec\phi(x),x\right):\nabla \vec\varphi_i(x) dx\!\!\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec{\vec{f}}\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorDiffOpInterface<ConfiguratorType, qc::QC_2D, NumCompArg, NumCompDest, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const ConfiguratorType _config;
  const InitType &_grid;
  const ArrayType &_levelSetFunction;

public:
  FENonlinLevelsetVectorDiffOpInterface ( const InitType &Grid,
                                          const ArrayType &LevelSetFunction ) :
      _config ( Grid ),
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // check, whether Arg and Dest have correct size
    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << Dest.numComponents() << " " << NumCompDest << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "nl" (it will contain the nonlinear function value at each quadrature point) as well as "val" (quadrature result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> nl;
    aol::Vec<NumCompDest, RealType> val;
    typename ConfiguratorType::VecType grad;
    const RealType h  = _grid.H();

    // traverse all lines on the interface
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _grid, _levelSetFunction );
    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {
      // compute center of mass for triangle (in local coordinates, not scaled by h)
      typename ConfiguratorType::VecType lineNodes[2];
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          lineNodes[i][j] = isoIt->points[i][j] - isoIt->el[j];
      typename ConfiguratorType::VecType centerOfMass ( lineNodes[0] + lineNodes[1] );
      centerOfMass /= static_cast<RealType> ( 2 );

      // compute values needed for midpoint quadrature
      const qc::Element el ( isoIt->el );
      const RealType length = sqrt ( aol::Sqr ( isoIt->points[0][0] - isoIt->points[1][0] ) + aol::Sqr ( isoIt->points[0][1] - isoIt->points[1][1] ) ) * h;
      asImp().getNonlinearity ( discrFuncsArg, el, centerOfMass, nl );

      // for each hat function compute its diff product with the nonlinear function and the midpoint quadrature of that
      const int numLocalDofs = _config.getNumLocalDofs ( el );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( el );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        bfs.evaluateGradient ( dof, centerOfMass, grad );
        val = nl * grad * length;

        for ( int d = 0; d < NumCompDest; d++ )
          Dest[d][ _config.localToGlobal ( el, dof ) ] += val[d];
      }
    }
  }

  // Interface function; has to be provided in derived classes.
  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, NL );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma\vec{\vec{f}}\left(\vec\phi(x),\nabla\vec\phi(x),x\right):\nabla \vec\varphi_i(x) dx\!\!\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor, \f$\vec\varphi_i\f$ is the vector-valued \f$i\f$th
 * finite element basis function, and \f$\vec{\vec{f}}\f$ has to be defined in the derived classes as "getNonlinearity".
 * \f$\vec\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp>
class FENonlinLevelsetVectorDiffOpInterface<ConfiguratorType, qc::QC_3D, NumCompArg, NumCompDest, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const ConfiguratorType _config;
  const InitType &_grid;
  const ArrayType &_levelSetFunction;
  tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> _cfeGrid;

public:
  FENonlinLevelsetVectorDiffOpInterface ( const InitType &Grid,
                                          const ArrayType &LevelSetFunction ) :
      _config ( Grid ),
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ),
      _cfeGrid ( _grid.getGridDepth() ) {
    _cfeGrid.setAdjustLevelset ( 1.0e-10 );
    _cfeGrid.setDomainFrom ( LevelSetFunction );
    _cfeGrid.detectVirtualNodes ();
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // check, whether Arg and Dest have correct size
    if ( Arg.numComponents() != NumCompArg || Dest.numComponents() != NumCompDest ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << Dest.numComponents() << " " << NumCompDest << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "nl" (it will contain the nonlinear function value at each quadrature point) as well as "val" (quadrature result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> nl;
    aol::Vec<NumCompDest, RealType> val;
    typename ConfiguratorType::VecType grad;
    const RealType h  = _cfeGrid.H();

    // traverse all triangles on the interface
    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      // compute center of mass for triangle
      std::vector< aol::Vec3< RealType > > triangleNodes = itit.getTriangRef().getLocalCoordinates();
      aol::Vec3<RealType> centerOfMass ( triangleNodes[0] + triangleNodes[1] + triangleNodes[2] );
      centerOfMass /= static_cast<RealType> ( 3 );

      // compute values needed for midpoint quadrature
      const qc::Element el ( itit.getTriangRef().getElement() );
      const RealType area = itit.getTriangRef().getArea ( h );
      asImp().getNonlinearity ( discrFuncsArg, el, centerOfMass, nl );

      // for each hat function compute its diff product with the nonlinear function and the midpoint quadrature of that
      const int numLocalDofs = _config.getNumLocalDofs ( el );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( el );
      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        bfs.evaluateGradient ( dof, centerOfMass, grad );
        val = nl * grad * area;

        for ( int d = 0; d < NumCompDest; d++ )
          Dest[d][ _config.localToGlobal ( el, dof ) ] += val[d];
      }
    }
  }

  // Interface function; has to be provided in derived classes.
  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( DiscFuncs, El, NL );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the integral
 * \f$\int_\Gamma f\left(\vec\phi(x),\nabla\vec\phi(x),x\right) dx\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$f\f$ has to be
 * defined in the derived classes as "getNonlinearity".
 * \f$\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, typename qc::Dimension Dim, int NumCompArg, typename Imp>
class FENonlinLevelsetIntegrationVectorInterface {};

/**
 * \brief Via "apply" or "applyAdd" this class computes the integral
 * \f$\int_\Gamma f\left(\vec\phi(x),\nabla\vec\phi(x),x\right) dx\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$f\f$ has to be
 * defined in the derived classes as "getNonlinearity".
 * \f$\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, typename Imp>
class FENonlinLevelsetIntegrationVectorInterface<ConfiguratorType, qc::QC_2D, NumCompArg, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const InitType &_grid;
  const ArrayType &_levelSetFunction;

public:
  FENonlinLevelsetIntegrationVectorInterface ( const InitType &Grid,
                                               const ArrayType &LevelSetFunction ) :
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // check, whether Arg has correct size
    if ( Arg.numComponents() != NumCompArg ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "val" (it will contain the result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Scalar<RealType> val;
    const RealType h  = _grid.H();

    // traverse all lines on the interface
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _grid, _levelSetFunction );
    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {
      // compute center of mass for triangle (in local coordinates, not scaled by h)
      typename ConfiguratorType::VecType lineNodes[2];
      for ( int i = 0; i < 2; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          lineNodes[i][j] = isoIt->points[i][j] - isoIt->el[j];
      typename ConfiguratorType::VecType centerOfMass ( lineNodes[0] + lineNodes[1] );
      centerOfMass /= static_cast<RealType> ( 2 );

      // midpoint quadrature
      const qc::Element el ( isoIt->el );
      const RealType length = sqrt ( aol::Sqr ( isoIt->points[0][0] - isoIt->points[1][0] ) + aol::Sqr ( isoIt->points[0][1] - isoIt->points[1][1] ) ) * h;
      Dest += length * asImp().evaluateIntegrand ( discrFuncsArg, el, centerOfMass );
    }
  }

  // Interface function; has to be provided in derived classes.
  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, RefCoord );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the integral
 * \f$\int_\Gamma f\left(\vec\phi(x),\nabla\vec\phi(x),x\right) dx\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$f\f$ has to be
 * defined in the derived classes as "getNonlinearity".
 * \f$\phi\f$ is the function passed to the apply-routines.
 *
 * Currently, only midpoint quadrature is implemented.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, typename Imp>
class FENonlinLevelsetIntegrationVectorInterface<ConfiguratorType, qc::QC_3D, NumCompArg, Imp> :
      public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;

  const InitType &_grid;
  const ArrayType &_levelSetFunction;
  tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> _cfeGrid;

public:
  FENonlinLevelsetIntegrationVectorInterface ( const InitType &Grid,
                                               const ArrayType &LevelSetFunction ) :
      _grid ( Grid ),
      _levelSetFunction ( LevelSetFunction ),
      _cfeGrid ( _grid.getGridDepth() ) {
    _cfeGrid.setAdjustLevelset ( 1.0e-10 );
    _cfeGrid.setDomainFrom ( LevelSetFunction );
    _cfeGrid.detectVirtualNodes ();
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // check, whether Arg has correct size
    if ( Arg.numComponents() != NumCompArg ) {
      cerr << Arg.numComponents() << " " << NumCompArg << " " << endl;
      throw aol::Exception ( "Mismatching number of vectors.", __FILE__, __LINE__ );
    }

    // make a function of Arg and define "val" (it will contain the result)
    aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> discrFuncsArg ( _grid, Arg );
    aol::Scalar<RealType> val;
    const RealType h  = _cfeGrid.H();

    // traverse all triangles on the interface
    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      // compute center of mass for triangle
      std::vector< aol::Vec3< RealType > > triangleNodes = itit.getTriangRef().getLocalCoordinates();
      aol::Vec3<RealType> centerOfMass ( triangleNodes[0] + triangleNodes[1] + triangleNodes[2] );
      centerOfMass /= static_cast<RealType> ( 3 );

      // midpoint quadrature
      const qc::Element el ( itit.getTriangRef().getElement() );
      const RealType area = itit.getTriangRef().getArea ( h );
      Dest += area * asImp().evaluateIntegrand ( discrFuncsArg, el, centerOfMass );
    }
  }

  // Interface function; has to be provided in derived classes.
  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               const typename ConfiguratorType::DomVecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, RefCoord );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma\vec\phi(x)\cdot\vec\varphi_idx\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$\vec\varphi_i\f$
 * is the \f$i\f$th finite element basis function.
 * \f$\phi\f$ is the function passed to the apply-routines.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumComp>
class LevelsetVectorMassOp :
      public FENonlinLevelsetVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, NumComp, NumComp, LevelsetVectorMassOp<ConfiguratorType, NumComp> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

public:
  LevelsetVectorMassOp ( const qc::GridDefinition &Grid,
                         const ArrayType &LevelSetFunction ) :
      FENonlinLevelsetVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, NumComp, NumComp, LevelsetVectorMassOp<ConfiguratorType, NumComp> > ( Grid, LevelSetFunction ) {}

  // Interface function; has to be provided in derived classes.
  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumComp> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &RefCoord,
                         aol::Vec<NumComp, typename ConfiguratorType::RealType> &NL ) const {
    DiscFuncs.evaluate ( El, RefCoord, NL );
  }
};

/**
 * \brief This class implements the operator (expressed in matrix form)
 * \f$\left(\int_\Gamma \varphi_i(x)\varphi_j(x) dx\right)_{ij}\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$\varphi_k\f$ is the \f$k\f$th
 * finite element basis function.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class LevelsetMassOp :
      public FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, ConfiguratorType::Dim, LevelsetMassOp<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

public:
  LevelsetMassOp ( const qc::GridDefinition &Grid,
                   const ArrayType &LevelSetFunction ) :
      FELinLevelsetScalarWeightedMassInterface<ConfiguratorType, ConfiguratorType::Dim, LevelsetMassOp<ConfiguratorType> > ( Grid, LevelSetFunction ) {}

  // Interface function; has to be provided in derived classes.
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::NumberTrait<RealType>::one;
  }
};

/**
 * \brief Via "apply" or "applyAdd" this class computes the vector
 * \f$\left(\int_\Gamma\nabla\vec\phi(x):\nabla\vec\varphi_idx\right)_i\f$,
 * where \f$\Gamma\f$ is the zero levelset passed to the constructor and \f$\vec\varphi_i\f$
 * is the \f$i\f$th finite element basis function.
 * \f$\phi\f$ is the function passed to the apply-routines.
 *
 * \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest>
class LevelsetVectorStiffOp :
      public FENonlinLevelsetVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, NumCompArg, NumCompDest, LevelsetVectorStiffOp<ConfiguratorType, NumCompArg, NumCompDest> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

public:
  LevelsetVectorStiffOp ( const typename ConfiguratorType::InitType &Grid, const ArrayType &LevelSetFunction ) :
      FENonlinLevelsetVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, NumCompArg, NumCompDest, LevelsetVectorStiffOp<ConfiguratorType, NumCompArg, NumCompDest> > ( Grid, LevelSetFunction ) {}

  void getNonlinearity ( aol::DiscreteVectorFunctionDefault<ConfiguratorType, NumCompArg> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         const typename ConfiguratorType::VecType &RefCoord,
                         aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> &NL ) const {
    DiscFuncs.evaluateGradient ( El, RefCoord, NL );
  }
};


} // end namespace

#endif
