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

#ifndef __VERYSPARSEMATRIX_H
#define __VERYSPARSEMATRIX_H

#include <sparseMatrices.h>
#include "lookupMap.h"

namespace aol {

/** Class for matrices that are very sparse in the sense that only few rows contain entries at all.
 *  \author Schwen (MEVIS)
 */
template <typename DataType>
class VerySparseMatrix : public Matrix<DataType> {
protected:
  aol::LookupMap< int, aol::SparseRow<DataType> > _rows;

public:
  //! Standard constructor creates matrix of size 0 by 0
  VerySparseMatrix ( );

  //! destructor
  ~VerySparseMatrix ( );

  //! (explicit) copy constructor
  explicit VerySparseMatrix ( const VerySparseMatrix<DataType> &other );

  //! (mathematical) assignment operator
  VerySparseMatrix<DataType>& operator= ( const VerySparseMatrix<DataType> &other );

  //! constructor for given size Rows x Columns
  VerySparseMatrix ( const int Rows, const int Columns );

  //! change size and delete old contents
  void reallocate ( const int Rows, const int Columns );

  DataType get ( int I, int J ) const;

  void set ( int I, int J, DataType Val );

  void add ( int I, int J, DataType Val );

  void setZero ( );

  void apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const;

  void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const;

  void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const;

  void makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const;

  ostream& printVerySparse ( ostream& out = cout, const aol::Format DataFormatter = aol::mixedFormat ) const;

  VerySparseMatrix<DataType>& operator*= ( const DataType Factor );
};

}

#endif
