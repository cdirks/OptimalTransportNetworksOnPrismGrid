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

/** \brief Example for using adaptive Grids
 *
 *  \author Paetz
 */

#include <FEOpInterface.h>
#include <configurators.h>
#include <bitVector.h>
#include <solver.h>
#include <quocMatrices.h>
#include <maskedVector.h>
#include <scalarArray.h>
#include <preconditioner.h>

#include "adaptiveGridTP.h"
#include "estimatorTP.h"


template <typename _RealType, typename _MatrixType>
class AdaptiveGridConfiguratorBaseUGB {
public:
  typedef typename AdaptiveGridTP<EstimatorTP2d<double>, double>::ElementIterator ElementIteratorType;
  typedef qc::Element                                                             ElementType;
  typedef AdaptiveGridTP<EstimatorTP2d<double>, double>                           InitType;
  typedef _RealType                                                               RealType;
  typedef _MatrixType                                                             MatrixType;

  static const qc::Dimension DimOfWorld = qc::QC_2D;

  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_ADAPTIVEGRID_INDEX_MODE;

protected:
  const InitType &_grid;   //!< memorize reference to the grid

public:

  const InitType& getGrid() const { return _grid;}

  const aol::Vec3<int> &_size;   //!< memorize grid size

  AdaptiveGridConfiguratorBaseUGB ( const AdaptiveGridTP<EstimatorTP2d<double>, double> &Grid ) :
    _grid ( Grid ), _size ( Grid.getSize() ) {
#ifdef VERBOSE
    cerr << "constructor = " << this->_size;
#endif
  }
  //! returns the begin iterator of the grid
  const typename AdaptiveGridTP<EstimatorTP2d<double>, double>::ElementIterator &begin( ) const {
    return this->_grid.makeElementIterator();
  }

  //! returns the end iterator of the grid
  const inline qc::EndElement &end( ) const {
    return this->_grid._endEl;
  }

  const InitType& getInitializer( ) const { return this->_grid; }

  int getNumGlobalDofs( ) const {
    return _size[0] * _size[1] * _size[2]; //this->_grid.getNumberOfNodes();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<DimOfWorld> ( _grid ) );
    return mat;
  }

  AdaptiveGridConfiguratorBaseUGB& operator= ( const AdaptiveGridConfiguratorBaseUGB& ) {
    throw aol::UnimplementedCodeException ( "qc::AdaptiveGridConfiguratorBaseUGB::operator= not implemented", __FILE__, __LINE__ );
  }

};






template <typename _RealType, typename _QuadType, typename _MatrixType>
class AdaptiveGridConfiguratorUGB : public AdaptiveGridConfiguratorBaseUGB<_RealType, _MatrixType> {
public:
  AdaptiveGridConfiguratorUGB ( const AdaptiveGridTP<EstimatorTP2d<double>, double> &Grid )
    : AdaptiveGridConfiguratorBaseUGB<_RealType, _MatrixType> ( Grid ), _baseFuncSetPointer ( NULL ) {

    _baseFuncSetPointer = new BaseFuncSetType*[Grid.getGridDepth() + 1];
    for ( int i = 0; i <= Grid.getGridDepth(); i++ ) {
      _baseFuncSetPointer[i] = new BaseFuncSetType ( 1.0 / _RealType ( 1 << i ) );
    }

  }

  typedef aol::Vector<_RealType>       VectorType;
  typedef _RealType                    RealType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> ArrayType;
  typedef qc::BitArray<qc::QC_2D>               MaskType;
  typedef aol::FullMatrix<_RealType>   FullMatrixType;

private:

  BaseFuncSetType **_baseFuncSetPointer;

public:
  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;


  static const int maxNumLocalDofs = 4;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 4;
  }

  inline RealType vol ( const qc::Element &El ) const {
    const RealType numEl = 1 << ( El.level() );
    return aol::Sqr ( 1.0 / numEl );
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &El ) const {
    return *_baseFuncSetPointer[El.level()];
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    const int elSize = 1 << ( this->getInitializer().getGridDepth() - El.level() );
    int r =
      ( El.x() + ( ( localIndex & 1 ) != 0 ) * elSize ) * 1       +
      ( El.y() + ( ( localIndex & 2 ) != 0 ) * elSize ) * this->_size[0];
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( El, localIndex0 );
    glob[1] = localToGlobal ( El, localIndex1 );
  }

  AdaptiveGridConfiguratorUGB<_RealType, _QuadType, _MatrixType>& operator= ( const AdaptiveGridConfiguratorUGB<_RealType, _QuadType, _MatrixType>& ) {
    throw aol::UnimplementedCodeException ( "qc::RectangularGridConfiguratorUGB<_RealType, qc::QC_2D, _QuadType,_MatrixType>::operator= not implemented", __FILE__, __LINE__ );
  }

};

#define KB 1024
#define CACHE_SIZE 8192 * KB

int main ( int /*argc*/, char** /*argv*/ ) {
  try {

    const int depth = 7;
    const int N = ( 1 << depth ) + 1;

    aol::Vec3<int>       size ( N, N, 1 );

    typedef AdaptiveGridConfiguratorUGB<double, aol::GaussQuadrature<double, qc::QC_2D, 3>, aol::SparseMatrix<double> > ConfiguratorType;
    typedef ConfiguratorType::InitType GridType;
    typedef ConfiguratorType::MatrixType MatrixType;

    double threshold = 0.01;
    EstimatorTP2d<double> est2 ( size[0] );
    est2.setThreshold ( threshold );

    GridType newGrid ( &est2 );

    int scaleMesh = 5;
    qc::ScalarArray<int, qc::QC_2D> mesh ( scaleMesh * ( size.x() - 1 ) + 1, scaleMesh * ( size.y() - 1 ) + 1 );
    qc::ScalarArray<double, qc::QC_2D> tmpSolution ( size ), boundary_conditions ( size ), rhs ( size ), solution ( size ), rhsTmp ( size );

    for ( int i = 0; i < N; i++ ) {
      boundary_conditions.set ( i, 0, 1.0 );
      for ( int j = 0; j < N; j++ ) {
        if ( i < 3.0 * N / 4.0 && j < 3.0 * N / 4.0 ) tmpSolution.set ( i, j, 0.0 );
        else tmpSolution.set ( i, j, 1.0 );
      }
    }

    newGrid.getEstimator()->makeSaturatedErrorArray ( tmpSolution );

    est2.saveAsVector ( "errorIndicator.dat" );
    est2 *= 255;
    est2.savePNG ( "errorIndicator.png" );

    aol::StiffOp< ConfiguratorType, aol::QUOC_ADAPTIVEGRID_INDEX_MODE > stiffOp ( newGrid, aol::ASSEMBLED );
    newGrid.initVertexArrays ( stiffOp );

    newGrid.getActivityVector ( rhsTmp );
    rhsTmp *= 255;
    rhsTmp.saveAsVector ( "hanging.dat" );
    rhsTmp.savePNG ( "hanging.png" );
    newGrid.drawMesh ( stiffOp, mesh, scaleMesh );

    MatrixType stiffMatrix ( newGrid );
    stiffMatrix.setZero();
    stiffOp.assembleAddMatrix ( stiffMatrix );

    stiffMatrix.apply ( boundary_conditions, rhs );
    rhs *= -1.0;

    for ( int i = 0; i < N; i++ ) {
      stiffMatrix.setRowColToZero ( i );
      stiffMatrix.setRowColToZero ( N * ( N - 1 ) + i );
      stiffMatrix.set ( N * ( N - 1 ) + i, N * ( N - 1 ) + i, 1.0 );
      rhs.set ( i, 0, 0.0 );
      rhs.set ( i, N - 1, 0.0 );
    }

    newGrid.setInactiveNodesFromMatrixToUnity ( stiffMatrix );

    aol::DiagonalPreconditioner<aol::Vector<double> > diagPre ( stiffMatrix );
    aol::PCGInverse< aol::Vector<double>, MatrixType > solver ( stiffMatrix, diagPre, 1e-16, 5000, aol::STOPPING_ABSOLUTE );

    solver.apply ( rhs, solution );
    solution += boundary_conditions;
    newGrid.interpolateResult ( solution );

    mesh.saveAsVector ( "mesh.dat" );
    solution.saveAsVector ( "solution.dat" );
    mesh.savePNG ( "mesh.png" );
    mesh *= 255;
    solution *= 255;
    solution.savePNG ( "solution.png" );
    solution.saveASCII ( "solution.txt" );

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
