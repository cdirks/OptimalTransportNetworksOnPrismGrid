#ifndef __ADAPTIVEPRISMGRIDCONFIGURATOR_H
#define __ADAPTIVEPRISMGRIDCONFIGURATOR_H

#include <sparseMatrices.h>
#include <smallMat.h>
#include <smallVec.h>

#include <adaptiveFEPrismMesh.h>
#include <baseFunctionSetPrismGrid.h>

/**
 * Configurator class for adaptive three-dimensional prism grid used for problem arising from functional lifting 
 * \author Dirks 
 */


template < typename _RealType, typename _GridType, typename _QuadType, typename _PrismType = PrismBaseElement < _RealType, int > >
class AdaptiveFEPrismMeshConfigurator {
    
public: 
    
  typedef _RealType RealType;
  typedef _GridType InitType;
  typedef _QuadType QuadType;  
  typedef aol::Vec3 < RealType > DomVecType;
  typedef aol::Mat < 3, 3, RealType > MatType;
  typedef aol::SparseMatrix < RealType > MatrixType;
  typedef aol::Vec3 < RealType > VecType;
  typedef BaseFunctionSetPrismGrid < _RealType, _QuadType, _PrismType > BaseFuncSetType; 
  typedef typename InitType::ElementIteratorType ElementIteratorType;
  typedef _PrismType ElementType;
  typedef typename InitType::MapType MapType;
  typedef aol::Vector < _RealType > VectorType; 
  
  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DimOfWorld = Dim;
  static const int maxNumLocalDofs = 6;
    
  AdaptiveFEPrismMeshConfigurator ( InitType & grid ) : _grid ( grid ) {}
    
    
  inline int getNumLocalDofs ( const ElementType & ) const { 
    return 6;
  }
  
  inline int getNumGlobalDofs () const { 
    return _grid.getNumberOfNodes() - _grid.getHangingNodeMap().size();  
  }
  
  const ElementIteratorType & begin () const {
    return _grid.begin_it;
  }
  
  inline const ElementIteratorType & end() const {
    return _grid.end_it;
  }
  
  //! index mapping for local dof to global dof
  inline int localToGlobal ( const ElementType & El, const int localIndex ) const { 
    return El.getGlobalNodeIndex ( localIndex );
  }
  
  //! create a new, clean matrix
  MatrixType * createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( this->_grid );
    return mat;
  }
  
  //! return corresponding base function set for element El
  const BaseFuncSetType & getBaseFunctionSet ( const ElementType & El ) const { 
    _baseFuncSet.setPrism ( El );
    return _baseFuncSet;  
  }
 
  RealType vol ( const ElementType & El ) const {
    return El.vol();
  }
  
  //! return grid reference
  const InitType & getInitializer( ) const { 
    return _grid; 
  }
  
  
  mutable BaseFuncSetType _baseFuncSet;
  
  
protected:
    
  InitType &_grid;      //! reference to the grid
    
};

#endif