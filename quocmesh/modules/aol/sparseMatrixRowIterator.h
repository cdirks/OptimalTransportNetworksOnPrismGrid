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

#ifndef __SPARSEMATRIXROWITERATOR_H
#define __SPARSEMATRIXROWITERATOR_H

/*
 * SparseMatrixRowIterator.h
 *
 *  Created on: Jul 19, 2012
 *      Author: Toelkes
 */

#include <rows.h>
#include <sparseMatrices.h>

namespace aol {

//! \brief Iterator for a sparse row of an aol::SparseMatrix < DataType >
//! \author Toelkes
template < typename DataType >
class SparseMatrixRowIterator : public std::vector < typename aol::Row < DataType >::RowEntry >::iterator {
protected:
  typedef typename std::vector < typename aol::Row < DataType >::RowEntry >           RowVectorType;
  typedef typename std::vector < typename aol::Row < DataType >::RowEntry >::iterator IteratorType;

  //! The sparse row
  aol::SparseRow < DataType > &_sparseRow;

  //! The vector containing the row entries (column and value).
  RowVectorType &_row;

public:
  //! Constructs the SparseMatrixRowIterator from a sparse matrix and a row number
  //! \param sparseMatrix Reference to the sparse matrix containing the row to be iterated over
  //! \param row          The row to be iterated over
  SparseMatrixRowIterator ( aol::SparseMatrix < DataType > &sparseMatrix, unsigned int row )
    : _sparseRow ( * ( static_cast < aol::SparseRow < DataType >* > ( sparseMatrix.rows[row] ) ) ), _row ( _sparseRow.row ) {
    IteratorType::operator= ( _row.begin () );
  }

  //! \brief Resets the iterator to the first entry of the row
  void setToFirstRowEntry () {
    *this = _row.begin ();
  }

  //! \brief Determine if the iterator reached the end of the row
  bool atEnd () const {
    return *this == _row.end ();
  }

  //! \brief Determine if the iterator hasn't reached the end of the row yet
  bool notAtEnd () const {
    return *this != _row.end ();
  }
};

} /* namespace aol */
#endif /* __SPARSEMATRIXROWITERATOR_H */
