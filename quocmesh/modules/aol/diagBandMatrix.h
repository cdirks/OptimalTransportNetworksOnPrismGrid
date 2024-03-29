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

#ifndef __DIAGBANDMATRIX_H
#define __DIAGBANDMATRIX_H

#include <bandMatrix.h>

namespace aol {

/** Class for squared banded matrices with fixed number of bands below and above the main diagonal.
 *  \author Schwen (MEVIS)
 */
template< typename DataType, int BLower, int BUpper = BLower >
class DiagBandMatrix : public aol::Matrix<DataType> {
protected:
  DataType*  _pData;
  const int _nBands;

public:
  //! Standard constructor creates matrix of size 0 by 0
  DiagBandMatrix ( );

  //! destructor
  ~DiagBandMatrix ( );

  //! (explicit) copy constructor
  explicit DiagBandMatrix ( const DiagBandMatrix<DataType,BLower,BUpper> &other );

  //! (mathematical) assignment operator
  DiagBandMatrix<DataType,BLower,BUpper>& operator= ( const DiagBandMatrix<DataType,BLower,BUpper> &other );

  //! constructor for given size
  explicit DiagBandMatrix ( const int NRowsCols );

  explicit DiagBandMatrix ( const qc::GridStructure &Size );

  DiagBandMatrix ( const int NRows, const int NCols );

  //! change size and delete old contents
  void reallocate ( const int nRowsCols );

  void reallocate ( const int nRows, const int nCols );

  DataType get ( int I, int J ) const;

  void set ( int I, int J, DataType Val );

  void add ( int I, int J, DataType Val );

  void setZero ( );

  void apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const;

  void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const;

  void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const;

  void makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const;

protected:
  inline int offset ( const int I, const int J ) const {
    return ( J - I );
  }

  inline int dataIndex ( const int I, const int J ) const {
    return ( ( _nBands ) * I + ( offset ( I, J ) + BLower ) );
  }

  inline int dataIndexRO ( const int Row, const int Offset ) const {
    return ( ( _nBands ) * Row + ( Offset + BLower ) );
  }

  inline bool validEntry ( const int I, const int J ) const {
    return ( J >= ( I - BLower ) && ( J <= ( I + BUpper ) ) );
  }
};

}

#endif
