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

#ifndef __SUBGRIDSPARSEMATRIX_H
#define __SUBGRIDSPARSEMATRIX_H

#include <sparseMatrices.h>
#include <gridSize.h>
#include <eigenvectors.h>
#include <indexTranslationMatrix.h>
#include <qmException.h>

/*    CONTENTS OF THIS FILE:
 *
 *  class SubGridSparseMatrix:  1.1 declaration
 *                              2.1 implementation
 *
 *  class ActionSubGridSparseMatrixEigenvalues: 1.2 declaration
 *                                                  implementation
 *
 */
namespace nb {

// ================================= 1.1 ====================================

/** a general, unstructured sparse matrix
 */
template <typename DataType, typename ConfigType>
class SubGridSparseMatrix : public aol::GenSparseMatrix<DataType> {

public:
  SubGridSparseMatrix ( const ConfigType & config );
  SubGridSparseMatrix ( const typename ConfigType::InitType & grid );
  SubGridSparseMatrix ( int numGlobalDofs );
  SubGridSparseMatrix ( const SubGridSparseMatrix& mat, aol::CopyFlag copyFlag = aol::DEEP_COPY );
  virtual ~SubGridSparseMatrix( );
  virtual aol::SparseRow<DataType>* newDefaultRow( ) const;
  void rebuild ( const ConfigType & config );

  template <class T>
  void addMatrixProduct ( const aol::Matrix<T> &M1, const aol::Matrix<T> &M2 );

  void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const;
  bool isRowSet ( int r ) const;

  int writeIndexTranslationSmallToLargeTo ( aol::Vector<int> & translation ) const;

  SubGridSparseMatrix<DataType, ConfigType> & operator= ( const SubGridSparseMatrix<DataType, ConfigType> & rhs );

private:
  void destroy();
};

//  ================================= 1.2 ===================================
//
//!  This class only works for block operators with equally structured block
//!  entries!
//!
//!  \author von Deylen

template <typename DataType, typename ConfigType>
class ActionBlockSubGridSparseMatrixEigenvalues {

public:
  template <typename BlockEntryType>
  ActionBlockSubGridSparseMatrixEigenvalues ( const aol::BlockOp<DataType, BlockEntryType> & blockOp ) :
      _blockOp ( blockOp.getNumRows(), blockOp.getNumCols() ) {
    for ( int i = 0; i < blockOp.getNumRows(); ++i )
      for ( int j = 0; j < blockOp.getNumCols(); ++j ) {
        SubGridSparseMatrix<DataType, ConfigType> * cast_matrix
        = dynamic_cast<const SubGridSparseMatrix<DataType, ConfigType> * > ( blockOp.getPointer ( i, j ) );
        _blockOp.setPointer ( i, j, cast_matrix );
      }
  }

  template <typename BlockOpType>
  aol::Vec2<int> getTotalSize ( const BlockOpType & blockOp ) {
    int n = blockOp.getNumRows();
    // int m = blockOp.getNumCols();
    aol::Vec2<int> ret;
    for ( int i = 0; i < n; ++i )
      if ( blockOp.getPointer ( i, i ) ) {
        const aol::Matrix<DataType> * M_ij = blockOp.getPointer ( i, i );
        if ( !M_ij )
          throw aol::Exception ( "Op::getTotalSize() : Not all contained operators are of "
                                 "type Matrix<RealType>.", __FILE__, __LINE__ );

        ret[0] += M_ij->getNumRows();
        ret[1] += M_ij->getNumCols();
      }
    return ret;
  }

  DataType getMaxEigenvalue() const {

    aol::Vec2<int> size = getTotalSize ( _blockOp );
    QUOC_ASSERT ( size[0] == size[1] );
    int dimLarge = size[0];

    QUOC_ASSERT ( _blockOp.getNumRows() == _blockOp.getNumCols() );
    int numBlocks = _blockOp.getNumRows();

    // create small matrix
    typedef aol::FullMatrix<DataType> FullType;
    aol::Vector<int> indexTranslation;
    int dimSmallSingleBlock = _blockOp ( 0, 0 ).writeIndexTranslationSmallToLargeTo ( indexTranslation );
    int dimSmall = dimSmallSingleBlock * numBlocks;
    aol::TransparentIndexTranslationMatrix<FullType> smallMatUnblocked ( dimSmall, dimLarge, indexTranslation );

    for ( int i = 0; i < numBlocks; ++i )
      for ( int j = 0; j < numBlocks; ++j ) {
        FullType smallMat ( dimSmallSingleBlock, dimSmallSingleBlock );
        smallMat.FillFromLarge ( _blockOp ( i, j ) );
        smallMatUnblocked.addMultipleAtPosition ( i * dimSmallSingleBlock, j * dimSmallSingleBlock, smallMat );
      }

    // calculate eigenvectors
    aol::MultiVector<DataType> ev_small;
    aol::QREigenvectorOp<FullType> evOp;
    evOp.apply ( smallMatUnblocked, ev_small );
    aol::Vector<DataType> eigenvalues ( ev_small[0] );

    // make eigenvectors on full grid from small eigenvectors
    // by translating all parts separately
    aol::MultiVector<DataType> ev_large ( dimSmall, dimLarge );
    for ( int i = 0; i < dimSmallSingleBlock; ++i )
      for ( int j = 0; j < numBlocks; ++j ) {
        aol::Vector<DataType> smallTmp ( dimSmallSingleBlock );
        aol::Vector<DataType> largeTmp ( _blockOp ( i, i ).getNumCols() );

        ev_small[i + 1].getBlock ( dimSmallSingleBlock * j, smallTmp );
        smallMatUnblocked.TranslateSmallToLarge ( smallTmp, largeTmp );
        ev_large[i].setBlock ( dimSmallSingleBlock * j, largeTmp );
      }

    cout << "Eigenvalues of DF[0][0]: " << endl << endl << eigenvalues << endl << endl;

    return eigenvalues[dimSmall - 1];
  }

protected:
  const aol::BlockOp<DataType, SubGridSparseMatrix<DataType, ConfigType> > _blockOp;
};


// ================================= 2.1 ====================================

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType>::
SubGridSparseMatrix ( const ConfigType &config )
  : aol::GenSparseMatrix<DataType> ( config.getNumGlobalDofs(),
                                     config.getNumGlobalDofs() ) {
  rebuild ( config );
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType>::
SubGridSparseMatrix ( const typename ConfigType::InitType & grid )
  : aol::GenSparseMatrix<DataType> ( qc::GridSize<ConfigType::DomDim> ( grid ) ) {
  rebuild ( ConfigType ( grid ) );
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType>::
SubGridSparseMatrix ( int numGlobalDofs )
  : aol::GenSparseMatrix<DataType> ( numGlobalDofs, numGlobalDofs ) {}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType>::
SubGridSparseMatrix ( const SubGridSparseMatrix& mat, aol::CopyFlag copyFlag )
    : aol::GenSparseMatrix<DataType> ( mat.getNumRows(), mat.getNumCols() ) {
  // copy Matrix members other than M, N separately
  this->format = mat.format;
  this->prettyFormat = mat.prettyFormat;

  switch ( copyFlag ) {
    case aol::STRUCT_COPY: {
      for ( int i = 0; i < this->getNumRows(); ++i )
      if ( mat.rows[i] )
        this->rows[i] = new aol::SparseRow<DataType>;
    }
    break;
    case aol::DEEP_COPY: {
      for ( int i = 0; i < this->getNumRows(); ++i )
      if ( mat.rows[i] )
        this->rows[i] = new aol::SparseRow<DataType> ( static_cast<const aol::SparseRow<DataType> & > ( * ( mat.rows[i] ) ) );
    }
    break;
    default: {
    throw aol::UnimplementedCodeException ( "This CopyFlag is not implemented yet.", __FILE__, __LINE__ );
    }
    break;
  }
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType>::
~SubGridSparseMatrix( ) {
  destroy();
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
aol::SparseRow<DataType>* SubGridSparseMatrix<DataType, ConfigType>::
newDefaultRow( ) const {
  return new aol::SparseRow<DataType>( );
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
void SubGridSparseMatrix<DataType, ConfigType>::
rebuild ( const ConfigType &config ) {
  for ( int i = 0; i < this->getNumRows(); i++ ) {
    if ( this->rows[i] ) {
      delete this->rows[i];
      this->rows[i] = NULL;
    }
  }
  typename ConfigType::MaskType nodeExistMask ( qc::GridSize<ConfigType::Dim> ( config.getInitializer() ) );
  config.writeNodeExistMaskTo ( nodeExistMask );
  for (int i = 0; i < nodeExistMask.size(); ++i)
    if ( nodeExistMask[i] && !this->rows[i] )
      this->rows[i] = newDefaultRow();
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
template <class T>
void SubGridSparseMatrix<DataType, ConfigType>::
addMatrixProduct ( const aol::Matrix<T> &M1, const aol::Matrix<T> &M2 ) {
  // cerr << "calling SubGridSparseMatrix::addMatrixProduct()\n";
  if (    this->getNumRows( ) != M1.getNumRows()
          || this->getNumCols( ) != M2.getNumCols( )
          || M1.getNumCols() != M2.getNumRows() )
    throw aol::Exception ( "aol::Matrix::addMatrixProduct(): incompatible "
                           "matrix sizes.\n", __FILE__, __LINE__ );

  vector<typename aol::Row<T>::RowEntry > vec1, vec2;

  for ( int i = 0; i < this->getNumRows(); i++ )
    if ( this->rows[i] ) {
      M1.makeRowEntries ( vec1, i );
      typename vector<typename aol::Row<T>::RowEntry >::iterator it1;
      for ( it1 = vec1.begin(); it1 != vec1.end(); ++it1 ) {
        M2.makeRowEntries ( vec2, it1->col );
        typename vector<typename aol::Row<T>::RowEntry >::iterator it2;
        for ( it2 = vec2.begin(); it2 != vec2.end(); ++it2 ) {
          if ( !this->rows[i] )
            cerr << it1->value * it2->value << " ";
          else
            add ( i, it2->col, this->getUnsetRowsDiagEntry() * it1->value * it2->value );
        }
      }
    }
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
void SubGridSparseMatrix<DataType, ConfigType>::
makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
  if ( this->rows[RowNum] ) {
    static_cast<aol::SparseRow<DataType>* > ( this->rows[RowNum] )->makeRowEntries ( vec, RowNum );
  } else {
    // nothing saved: return row from identity matrix.
    vec.resize ( 1 );
    vec[0].col = RowNum;
    vec[0].value = this->getUnsetRowsDiagEntry();
  }
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
bool SubGridSparseMatrix<DataType, ConfigType>::
isRowSet ( int r ) const {
  return ( this->rows[ r ] != NULL );
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
int SubGridSparseMatrix<DataType, ConfigType>::
writeIndexTranslationSmallToLargeTo ( aol::Vector<int> & translation ) const {
  translation.resize ( this->getNumRows() );
  int index = 0;
  for ( int i = 0; i < this->getNumRows(); ++i )
    if ( isRowSet ( i ) )
      translation[index++] = i;
  translation.resize ( index );
  return index;
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
SubGridSparseMatrix<DataType, ConfigType> & SubGridSparseMatrix<DataType, ConfigType>::
operator= ( const SubGridSparseMatrix & rhs ) {
  setUnsetRowsDiagEntry ( rhs.getUnsetRowsDiagEntry() );

  if ( this->getNumRows() != rhs.getNumRows() )
    throw aol::Exception ( "SubGridSparseMatrix::operator==(...): "
                           "row count does not match.", __FILE__, __LINE__ );
  if ( this->getNumCols() != rhs.getNumCols() )
    throw aol::Exception ( "SubGridSparseMatrix::operator==(...): "
                           "col count does not match.", __FILE__, __LINE__ );

  size_t n = this->rows.size();
  for ( size_t i = 0; i < n; ++i ) {
    // delete old content
    if ( this->rows[i] )
      delete this->rows[i];

    // copy from rhs
    if ( rhs.rows[i] ) {
      // although this is a downcast, we use static_cast instead
      // of dynamic_cast for reasons of speed. By getting
      // a SubGridSparseMatrix as rhs, we can be totally sure that
      // rhs.rows only contains SparseRows.
      const aol::SparseRow<DataType> * rhs_row
      = static_cast<const aol::SparseRow<DataType> * > ( rhs.rows[i] );
      this->rows[i] = new aol::SparseRow<DataType> ( *rhs_row );
    } else
      this->rows[i] = NULL;
  }
  return *this;
}

//---------------------------------------------------------------------------

template <typename DataType, typename ConfigType>
void SubGridSparseMatrix<DataType, ConfigType>::
destroy() {
  for ( unsigned int i = 0; i < this->rows.size(); i++ ) {
    delete this->rows[ i ];
  }
}

//---------------------------------------------------------------------------

} // end of namespace nb.

#endif
