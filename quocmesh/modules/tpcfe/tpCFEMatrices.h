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

#ifndef __TPCFEMATRICES_H
#define __TPCFEMATRICES_H

#include <matrix.h>
#include <bandMatrix.h>
#include <sparseMatrices.h>
#include <quocMatrices.h>

#include <arrayExtensions.h>
#include <iterators.h>

#include <tpCFEGrid.h>
#include <tpCFELookup.h>

namespace tpcfe {

// implement efficient matrices for tpcfe computations here

// A. based on genSparseMatrix (row-wise organized matrices - probably slower but more memory-efficient)

template <typename GridType>
class AffineFE3D15Row : public aol::Row< typename GridType::RealType > {
public:
  typedef typename GridType::RealType DataType;

protected:
  const GridType &_grid;                    // sizeof ( * )
  DataType *_thisRow;                       // 15 * sizeof( DataType )

  mutable int _lastReturnedRowEntry;        // sizeof ( int )

public:
  explicit AffineFE3D15Row ( const GridType &Grid ) : _grid ( Grid ), _thisRow ( NULL ) {
    initGetRowEntries();
    if ( _grid.getDimOfWorld() != qc::QC_3D ) {
      throw aol::UnimplementedCodeException ( "AffineFE15Row not implemented yet for dimension different from  3", __FILE__, __LINE__ );
    }
    _thisRow = new DataType[ NUM_NEIGHBORS ];
    memset ( _thisRow, 0, NUM_NEIGHBORS * sizeof ( DataType ) );
  }

private:
  AffineFE3D15Row ( ); // class has const reference member, so standard constructor should not be called.


public:
  explicit AffineFE3D15Row ( const AffineFE3D15Row<GridType> &other ) : _grid ( other._grid ), _thisRow ( NULL ), _lastReturnedRowEntry ( other._lastReturnedRowEntry ) {
    if ( &_grid != & ( other._grid ) ) {
      throw aol::Exception ( "AffineFE3D15Row<GridType>: copy constructor does not work for different grid", __FILE__, __LINE__ );
    }
    _thisRow = new DataType[ NUM_NEIGHBORS ];
    for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
      _thisRow[i] = other._thisRow[i];
    }
  }

  AffineFE3D15Row<GridType>& operator= ( const AffineFE3D15Row<GridType> &other ) {
    if ( &_grid != & ( other._grid ) ) {
      throw aol::Exception ( "AffineFE3D15Row<GridType>::operator= cannot assign row for different grid", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
      _thisRow[i] = other._thisRow[i];
    }
    _lastReturnedRowEntry = other._lastReturnedRowEntry;
    return ( this );
  }

  ~AffineFE3D15Row() {
    if ( _thisRow ) {
      delete[] _thisRow;
      _thisRow = NULL;
    }
  }

  int numStoredEntries() const {
    return NUM_NEIGHBORS;
  }

  int numNonZeroes ( ) const {
    short nnz = 0;
    for ( int i = 0; i < NUM_NEIGHBORS; ++i )
      if ( _thisRow[i] != aol::NumberTrait<DataType>::zero )
        ++nnz;
    return ( nnz );
  }

  DataType get ( int I, int J ) const {
    for ( short k = 0; k < NUM_NEIGHBORS; ++k )
      if ( J == ( I + _grid.getCfeGridIndexOffset ( k ) ) )
        return _thisRow[k];

    // if not found:
    return 0.0;
  }

  void set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    bool okay = ( Value == static_cast<DataType> ( 0.0 ) );
#endif
    for ( short k = 0; k < NUM_NEIGHBORS; ++k )
      if ( J == I + _grid.getCfeGridIndexOffset ( k ) ) {
        _thisRow[k] = Value;
#ifdef BOUNDS_CHECK
        okay = true;
#endif
      }
#ifdef BOUNDS_CHECK
    if ( !okay ) {
      cerr << I << ", " << J << " " << Value << endl;
      throw aol::Exception ( "tpcfe::AffineFE15Row::set: This position may not be used.", __FILE__, __LINE__ );
    }
#endif
  }

  /** row-vector scalar mulitplication
   */
  DataType mult ( const aol::Vector<DataType> &src, const int rowIndex ) const {
    DataType dst = aol::NumberTrait<DataType>::zero;
    for ( short k = 0; k < NUM_NEIGHBORS; ++k )
      dst += _thisRow[k] * src[ rowIndex + _grid.getCfeGridIndexOffset ( k ) ];
    return ( dst );
  }

  virtual DataType multMaskedFunctorTrue ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorTrue> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorFalse ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorFalse> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorIdentity ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorIdentity> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorNegate ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorNegate> ( Src, Row, Mask );
  }

  /** masked row-vector scalar mulitplication (if desired only include nodes that are (not) masked)
   */
  template <typename BitMaskFunctorType>
  DataType multMasked ( const aol::Vector<DataType> &src, int rowIndex,
                        const aol::BitVector & Mask ) {
    BitMaskFunctorType maskFunctor;
    DataType dst = static_cast<DataType> ( 0.0 );
    for ( short k = 0; k < NUM_NEIGHBORS; ++k )
      if ( maskFunctor ( Mask[rowIndex + _grid.getCfeGridIndexOffset ( k ) ] ) )
        dst += _thisRow[k] * src[ rowIndex + _grid.getCfeGridIndexOffset ( k ) ];
    return ( dst );
  }


  /** row is mapped to scalar * row
   */
  void scale ( DataType factor ) {
    for ( short k = 0; k < NUM_NEIGHBORS; ++k )
      _thisRow[k] *= factor;
  }

  /** compute row sum
   */
  DataType sum ( int /*I*/ ) {
    // I is ignored
    DataType result = 0;
    for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
      result += _thisRow[k];
    }
    return ( result );
  }

  void scale ( int /*I*/, DataType factor ) {  // 1st parameter is ignored
    scale ( factor );
  }

  void add ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    bool okay = ( Value == static_cast<DataType> ( 0.0 ) );
#endif
    for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
      if ( J == ( I + _grid.getCfeGridIndexOffset ( k ) ) ) {
        _thisRow[k] += Value;
#ifdef BOUNDS_CHECK
        okay = true;
#endif
      }
    }
#ifdef BOUNDS_CHECK
    if ( !okay ) {
      cerr << I << ", " << J << " " << Value << endl;
      throw aol::Exception ( "tpcfe::AffineFE15Row::add: This position may not be used.", __FILE__, __LINE__ );
    }
#endif
  }

  void addMultiple ( const int /* RowNum */, const tpcfe::AffineFE3D15Row<GridType> &other, const DataType factor = aol::NumberTrait<DataType>::one ) {
    for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
      _thisRow[i] += factor * other._thisRow[i];
    }
  }

  using aol::Row<DataType>::addMultiple;

  void setZero() {
    memset ( _thisRow, 0, sizeof ( DataType ) * NUM_NEIGHBORS );
  }


  bool checkForNANsAndINFs() const {
    for ( short i = 0; i < NUM_NEIGHBORS; ++i ) {
      if ( !aol::isFinite ( _thisRow[ i ] ) ) {
        return true;
      }
    }
    return false;
  }

  //! use this method only if you know what you are doing!
  DataType* getDataPointer ( ) const {
    return ( _thisRow );
  }

protected:

  virtual void makeSortedRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( NUM_NEIGHBORS );

    for ( short i = 0; i < NUM_NEIGHBORS; ++i ) {
      vec[i].col   = RowNum + _grid.getCfeGridIndexOffset ( i ) ;
#ifdef BOUNDS_CHECK
      if ( vec[i].col < 0 || vec[i].col > _grid.getNumberOfNodes() ) {
        cerr << vec[i].col << endl;
        throw aol::Exception ( "tpcfe::AffineFE15Row::makeRowEntries: Illegal entry in makeRowEntries", __FILE__, __LINE__ );
      }
#endif
      vec[i].value = _thisRow[i];
    }
  }

  virtual void makeRowEntries (  vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    makeSortedRowEntries ( vec, RowNum );
  }

  void initGetRowEntries() const {
    _lastReturnedRowEntry = -1;
  }

  void getNextRowEntry ( typename aol::Row<DataType>::RowEntry &Entry, int I ) const {
    if ( _lastReturnedRowEntry < ( NUM_NEIGHBORS - 1 ) ) {
      ++_lastReturnedRowEntry;
      Entry.col = I + _grid.getCfeGridIndexOffset ( _lastReturnedRowEntry ) ;
      Entry.value = _thisRow[_lastReturnedRowEntry];
#ifdef BOUNDS_CHECK
      if ( Entry.col < 0 || Entry.col > _grid.getNumberOfNodes() ) {
        cerr << Entry.col << endl;
        throw aol::Exception ( "tpcfe::AffineFE15Row::getNextRowEntry: Illegal entry in getNextRowEntry", __FILE__, __LINE__ );
      }
#endif
    } else {
      Entry.col = -1;
      Entry.value = 0.0;
    }
  }

};

/** This matrix type is somewhat similar to FL's old SparseBandMatrix.
 *  Only one default row is stored (and once set it is used for all other default rows that cannot be modified -- all write access is completely ignored), so writing and reading must be done in an appropriate way.
 *  By default, this matrix type only works for constant coefficients. If weighted ops are used with varying weights, the matrix can be forced to store all rows by calling CFEGrid::setCoarsenedFlag ( true ).
 */
template < typename GridType >
class CFEHybridMatrixBase : public aol::Matrix< typename GridType::RealType > {
public:
  typedef typename GridType::RealType           DataType;
  typedef unsigned char                         RowTypeType;

  static const RowTypeType SavedRowBitMask        = static_cast<RowTypeType> ( 128 );
  static const RowTypeType DefaultRowBitMask      = static_cast<RowTypeType> (  64 );
  static const RowTypeType ExplicitRowBitMask     = static_cast<RowTypeType> ( SavedRowBitMask );
  static const RowTypeType ADefaultRowBitMask     = static_cast<RowTypeType> ( SavedRowBitMask | DefaultRowBitMask );
  static const RowTypeType DefaultRowIndexBitMask = static_cast<RowTypeType> ( ~ADefaultRowBitMask );
  static const RowTypeType UnsetRowBitMask        = static_cast<RowTypeType> (   0 );

protected:
  const GridType&                               _grid;

  qc::AArray< aol::Row<DataType>*, qc::QC_3D >  _rows;                  // Pointers to rows

  RowTypeType                                   _numDefaultRowTypes;    // number of different default rows (CFE_CD: 2 (interior/exterior), CFE_LIEHR: 2 (-/+)

  int                                           _numExplicit;           // number of explicit rows
  aol::Vector<int>                              _numDefault;            // number of default rows of the different types

  qc::ScalarArray< RowTypeType, qc::QC_3D >     _rowTypes;              // types of the rows at different positions

protected:
  bool isSavedRow   ( const int i ) const {  return ( ( _rowTypes[i] & SavedRowBitMask    ) != 0 ); } // and this is exactly the criterion for write access.
  bool isDefaultRow ( const int i ) const {  return ( ( _rowTypes[i] & DefaultRowBitMask  ) != 0 ); }

  bool isExplicitRow    ( const int i ) const { return ( isSavedRow ( i ) && !isDefaultRow ( i ) ); }
  bool isADefaultRow    ( const int i ) const { return ( ( _rowTypes[i] & ADefaultRowBitMask ) != 0 ); }
  bool isUnsetRow       ( const int i ) const { return ( _rowTypes[i] == UnsetRowBitMask      ); }

  RowTypeType whichDefaultRow ( const int i ) const { return ( _rowTypes[i] & DefaultRowIndexBitMask ); }

public:
  CFEHybridMatrixBase ( const GridType &Grid, const int numDefaultRowTypes /*should be determined by grid; possibly more for multiply jumping coefficient*/ )
      : aol::Matrix<DataType> ( Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ),
      _grid ( Grid ),
      _rows ( Grid ),
      _numDefaultRowTypes ( numDefaultRowTypes ),
      _numExplicit ( 0 ), _numDefault ( numDefaultRowTypes ),
      _rowTypes ( Grid ) {

    _rows.setAll ( NULL );
    _rowTypes.setAll ( UnsetRowBitMask );

  }

  virtual ~CFEHybridMatrixBase() {}

private:
  CFEHybridMatrixBase ( ); // class contains const reference, thus standard and copy constructor and assignment operator should not be called.
  explicit CFEHybridMatrixBase ( const CFEHybridMatrixBase<GridType> &other );
  CFEHybridMatrixBase<GridType>& operator= ( const CFEHybridMatrixBase<GridType> &other );

private:
  virtual void init ( ) = 0;
  // This function is pure virtual so that it must be implemented in derived classes.
  // The constructor of derived classes must call init because the constructor of the parent class cannot call virtual methods.


public:
  virtual DataType get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
      cerr << I << " " << J << " is out of bounds.\n";
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::get: Index out of bounds", __FILE__, __LINE__ );
    }
    if ( _rows[I] == NULL )
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::get: row unset", __FILE__, __LINE__ );
#endif
    return ( _rows[I]->get ( I, J ) );
  }

  virtual void set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
      cerr << I << " " << J << " is out of bounds.\n";
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::set: Index out of bounds", __FILE__, __LINE__ );
    }
    if ( _rows[I] == NULL )
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::set: row unset", __FILE__, __LINE__ );
#endif

    if ( isSavedRow ( I ) ) {
      _rows[I]->set ( I, J, Value );
    }      // else do nothing
  }

  virtual void add ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
      cerr << I << " " << J << " is out of bounds.\n";
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::add: Index out of bounds", __FILE__, __LINE__ );
    }
    if ( _rows[I] == NULL )
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::add: row unset", __FILE__, __LINE__ );
#endif

    if ( isSavedRow ( I ) ) {
      _rows[I]->add ( I, J, Value );
    } // else do nothing
  }

  void addMultiple ( const tpcfe::CFEHybridMatrixBase<GridType> &other, const DataType factor ) {
#ifdef BOUNDS_CHECK
    if ( _rowTypes != other._rowTypes )
      throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::addMultiple: incompatible matrices", __FILE__, __LINE__ );
    // maybe do more checking
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( isSavedRow ( i ) ) {
        const type_info &currentTypeId = typeid ( * ( _rows[i] ) );

        if ( currentTypeId == typeid ( * ( other._rows[i] ) ) ) {
          static const type_info &AffineFE3D15RowTypeId = typeid ( tpcfe::AffineFE3D15Row<GridType> );
          static const type_info &SparseRowTypeId       = typeid ( aol::SparseRow<DataType> );
          static const type_info &DiagonalRowTypeId     = typeid ( aol::DiagonalRow<DataType> );

          if ( currentTypeId == AffineFE3D15RowTypeId ) {
            dynamic_cast< tpcfe::AffineFE3D15Row<GridType>* > ( _rows[i] )->addMultiple ( i, * ( dynamic_cast< tpcfe::AffineFE3D15Row<GridType>* > ( other._rows[i] ) ), factor );
          } else if ( currentTypeId == SparseRowTypeId ) {
            dynamic_cast< aol::SparseRow<DataType>* >         ( _rows[i] )->addMultiple ( i, * ( dynamic_cast< aol::SparseRow<DataType>* >         ( other._rows[i] ) ), factor );
          } else if ( currentTypeId ==  DiagonalRowTypeId ) {
            dynamic_cast< aol::DiagonalRow<DataType>* >       ( _rows[i] )->addMultiple ( i, * ( dynamic_cast< aol::DiagonalRow<DataType>* >       ( other._rows[i] ) ), factor );
          } else {
            throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::addMultiple: Illegal typeid", __FILE__, __LINE__ );
          }
        } else {
          throw aol::Exception ( "tpcfe::CFEHybridMatrixBase::addMultiple: Different typeids", __FILE__, __LINE__ );
        }
      }
    }
  }

  using aol::Matrix<DataType>::operator+=;

  tpcfe::CFEHybridMatrixBase<GridType>& operator+= ( const tpcfe::CFEHybridMatrixBase<GridType> &other ) {
    addMultiple ( other, aol::NumberTrait<DataType>::one );
    return ( *this );
  }

  void applyAdd ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {

    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      throw ( aol::Exception ( "Wrong dimension.", __FILE__, __LINE__ ) );
    }

    // applyAdd *this to each component of Arg:
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < this->_numRows; ++i ) {
      Dest[i] += _rows[i]->mult ( Arg, i );
    }

  }

  void apply ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      throw ( aol::Exception ( "Wrong dimension.", __FILE__, __LINE__ ) );
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < this->_numRows; ++i ) {
      Dest[ i ] = _rows[i]->mult ( Arg, i );
    }

  }

  void makeRowEntries ( std::vector< typename aol::Row< DataType >::RowEntry > &vec, const int RowNum ) const {
    _rows[ RowNum ]->makeRowEntries ( vec, RowNum );
    // a distinction of different rowtype cases and precomputing row entry vectors did not turn out to be more efficient. double check this!
  }

  void makeSortedRowEntries ( std::vector< typename aol::Row< DataType >::RowEntry > &vec, const int RowNum ) const {
    _rows[ RowNum ]->makeSortedRowEntries ( vec, RowNum );
    // a distinction of different rowtype cases and precomputing row entry vectors did not turn out to be more efficient. double check this!
  }

  virtual void printStatistics ( ) const {
  }

private:
  // methods we wish to hide

  virtual aol::SparseRow<DataType>* newDefaultRow() const {
    // pure virtual function in base class, must be implemented.
    throw aol::Exception ( "CFEHybridMatrixBase::newDefaultRow() does not make sense", __FILE__, __LINE__ );
    return new aol::SparseRow<DataType>();
  }

  void setIdentity() {
    throw aol::UnimplementedCodeException ( "CFEHybridMatrixBase::setIdentity(): Standard setIdentity() incompatible with this type of matrix.", __FILE__, __LINE__ );
  }

  void reallocate ( const int /*M*/, const int /*N*/ ) {
    throw aol::UnimplementedCodeException ( "CFEHybridMatrixBase::reallocate( int, int ) does not make sense", __FILE__, __LINE__ );
  }

  void setZero() {
    throw aol::UnimplementedCodeException ( "CFEHybridMatrixBase::setZero(): setZero() incompatible with this type of matrix.", __FILE__, __LINE__ );
    // might work, but needs to be checked
  }

public:
  CFEHybridMatrixBase ( const int Rows, const int Cols ) : aol::Matrix<DataType> ( Rows, Cols ), _grid ( NULL ) {
    throw aol::UnimplementedCodeException ( "CFEHybridMatrixBase::CFEHybridMatrixBase( int, int ) does not make sense", __FILE__, __LINE__ );
  }

  // end of class CFEHybridMatrixBase
};


/** "Basis" class for template specialization
 */
template < typename GridType >
class CFEHybridMatrix : public CFEHybridMatrixBase< GridType > { };


/** Memory-efficient matrix class for CFE (CD complicated domain) computations storing only one instance of default rows.
 *  An important property is that domain and Dirichlet nodes must be set in the grid before instantiating the matrix and must not be changed afterwards.
 *  ATTENTION: use this matrix only if you really know what you are doing and compile with --IKnowWhatIAmDoing :-)
 */
template < typename DataType, typename NodalCoeffType >
class CFEHybridMatrix< CFEGrid<DataType, CFE_CD, NodalCoeffType > > : public CFEHybridMatrixBase< CFEGrid<DataType, CFE_CD, NodalCoeffType > > {
protected:
  typedef CFEGrid<DataType, CFE_CD, NodalCoeffType > GridType;
  aol::DiagonalRow<DataType>       theDefaultExteriorRow;
  tpcfe::AffineFE3D15Row<GridType> theDefaultInteriorRow;

public:
  explicit CFEHybridMatrix ( const GridType &Grid ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), theDefaultExteriorRow ( ), theDefaultInteriorRow ( Grid ) {
    init();
  }

  CFEHybridMatrix ( const GridType &Grid, const qc::BitArray<qc::QC_3D> &InteriorRows, const qc::BitArray<qc::QC_3D> &ExteriorRows, const qc::BitArray<qc::QC_3D> &ExplicitRows ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), theDefaultExteriorRow ( ), theDefaultInteriorRow ( Grid ) {
    init ( InteriorRows, ExteriorRows, ExplicitRows );
  }

  ~CFEHybridMatrix ( ) {
    for ( int i = 0; i < this->_rows.size(); ++i ) {
      if ( this->isExplicitRow ( i ) ) {
        delete ( this->_rows[i] );
      }
      this->_rows[i] = NULL;
    }
  }

  virtual void printStatistics ( ) const {
    const int64_t we_require = this->_numExplicit             * ( NUM_NEIGHBORS * sizeof ( DataType ) + sizeof ( void* ) + sizeof ( int ) ) +  this->_grid.getNumberOfNodes() * ( sizeof ( typename CFEHybridMatrixBase<GridType>::RowTypeType )  + sizeof ( void* ) ); // + O(1)
    const int64_t fs_require = this->_grid.getNumberOfNodes() * ( NUM_NEIGHBORS * sizeof ( DataType ) + sizeof ( void* ) + sizeof ( int ) ); // + O(1)
    const double saving     = 100 * ( 1.0 - static_cast<double> ( we_require ) / static_cast<double> ( fs_require ) );

    cerr << this->_numExplicit << " rows explicitely stored, " << this->_numDefault.sum() << " default rows, " << endl;
    cerr << "Memory needed: approximately " << we_require << " bytes, " << we_require / 1048576 << " MiB." << endl;
    cerr << "Full 15 band storage would require approximately " << fs_require << " bytes (" << fs_require / 1048576 << " MiB), so we save " << saving << " percent." << endl;
  }

private:
  virtual void init ( ) {
    // first store sparseRows for all boundary nodes (of bounding box, i. e. cube)
    // \todo think about whether to distinguish between non-domain/non-DOF and other
    for ( qc::RectangularBoundaryIterator<qc::QC_3D> fbnit ( this->_grid ); fbnit.notAtEnd(); ++fbnit ) {
      this->_rows.getRef ( *fbnit ) = new aol::SparseRow< DataType >;  // boundary rows must be sparse rows because AffineFE rows have entries outside
      this->_rowTypes.set ( *fbnit, this->ExplicitRowBitMask );
      ++ ( this->_numExplicit );
#ifdef VERBOSE
      cerr << "Setting " << *fbnit << " = " << this->_grid.getIndexMapperRef() ( *fbnit ) << " to explicit sparse Row" << endl;
#endif
    }

    // then find rowTypes for all interior nodes:
    // all Dirichlet nodes and their neighbors as well as all near-interface domain nodes must be stored explicitly as AffineFE3D rows.
    // note that cube-boundary nodes were already set and are ignored below.
    qc::BitArray<qc::QC_3D> storeExplicitly ( qc::GridSize<qc::QC_3D> ( this->_grid ) );
    storeExplicitly.setAll ( false );

    if ( this->_grid.isCoarsened() ) {

      storeExplicitly.setAll ( true );

    } else {

      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );

        if ( ( this->_grid.isDOFNode ( rowIndex ) &&  ! ( this->_grid.isInsideDomain ( *fnit ) ) )   ||   this->_grid.isDirichletNode ( rowIndex ) ) {
          // neighbors of DOF nodes outside the domain and neighbors of Dirichlet nodes must be stored explicitely
          for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
            qc::CoordType nbpos = ( *fnit ) + tpcfe::CFETopoLookup::_tetNbOffsets[i];

            if ( this->_grid.isAdmissibleNode ( nbpos ) ) {
              storeExplicitly.set ( nbpos, true );
            }

          }
        }
      }

    }

    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
      if ( this->_rows[ rowIndex ] == NULL ) {
        // only treat those not set yet (i. e. only treat cube-interior)

        if ( storeExplicitly[ rowIndex ] ) {

          // all nodes marked for explicit storage: store explicit AffineFE3D15Row

          this->_rows[ rowIndex ] =  new tpcfe::AffineFE3D15Row<GridType> ( this->_grid );
          this->_rowTypes[ rowIndex ] = this->ExplicitRowBitMask;
          ++ ( this->_numExplicit );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to explicit AffineFE3D15 Row" << endl;
#endif

        } else if ( this->_grid.isInsideDomain ( *fnit ) ) {

          // cube-interior, domain-interior, far from interface and Dirichlet: DefaultInterior

          this->_rows[ rowIndex ] = &theDefaultInteriorRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 0; // 0 for interior
          ++ ( this->_numDefault[0] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultInterior Row" << endl;
#endif

        } else {

          // cube-interior domain-exterior and far from interface (i. e. non-DOF), far from Dirichlet: DefaultExterior

          this->_rows[ rowIndex ] = &theDefaultExteriorRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 1; // 1 for exterior
          ++ ( this->_numDefault[1] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultExterior Row" << endl;
#endif

        }

      }
    }

    // detect TheDefaultRows:
    for ( unsigned char defaultType = 0; defaultType < 2; ++defaultType ) {
      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
        if ( this->_rowTypes[ rowIndex ] == ( this->DefaultRowBitMask | defaultType ) ) {
          this->_rowTypes[ rowIndex ] |= this->SavedRowBitMask;
          break;
        }
      }
    }

#ifdef VERBOSE
    for ( int i = 0; i < this->getNumRows(); ++i )
      cerr << "Row " << aol::intFormat ( i ) << " has rowType " << static_cast<int> ( this->_rowTypes[i] ) << "; " << this->_rows[i] << endl;
#endif

  }

  void init ( const qc::BitArray<qc::QC_3D> &/*InteriorRows*/, const qc::BitArray<qc::QC_3D> &/*ExteriorRows*/, const qc::BitArray<qc::QC_3D> &/*ExplicitRows*/ ) {
    throw aol::UnimplementedCodeException ( "init ( 3 bitArrays ) not implemented yet", __FILE__, __LINE__ );
  }

public:
  void setRowToZero ( const int Row ) {

#ifdef BOUNDS_CHECK
    if ( Row < 0 && Row >= this->_numCols ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEHybridMatrix (CFE_CD): setRowToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    if ( this->isSavedRow ( Row ) ) {
      this->_rows[ Row ]->setZero();
    } else if ( this->isDefaultRow ( Row ) && this->whichDefaultRow ( Row ) == 1 ) {
      // do nothing
    } else {
      cerr << Row << endl;
      throw aol::Exception ( "CFEHybridMatrix (CFE_CD): Row may not be modified", __FILE__, __LINE__ );
    }

  }

  void setColToZero ( const int Col ) {
#ifdef BOUNDS_CHECK
    if ( Col < 0 && Col >= this->_numRows ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEHybridMatrix (CFE_CD)::setColToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
      int row = Col + this->_grid.getCfeGridIndexOffset ( k );
      if ( row >= 0 && row < this->_numCols ) { // else neighbor does not exist: do nothing.
        this->set ( row, Col, aol::NumberTrait<DataType>::zero ); // we rely on row doing the correct thing when setting to zero.
      }
    }
  }

  // end of class CFEHybridMatrix<DataType, CFE_CD, NodalCoeffType>
};


/** Memory-efficient matrix class for CFE (LIEHR jumping coefficient) computations storing only one instance of default rows.
 *  An important property is that Dirichlet nodes must be set in the grid before instantiating the matrix and must not be changed afterwards.
 *  ATTENTION: use this matrix only if you really know what you are doing and compile with --IKnowWhatIAmDoing :-)
 */
template < typename DataType >
class CFEHybridMatrix< CFEGrid<DataType, CFE_LIEHR > > : public CFEHybridMatrixBase< CFEGrid<DataType, CFE_LIEHR > > {
protected:
  typedef CFEGrid<DataType, CFE_LIEHR > GridType;
  tpcfe::AffineFE3D15Row<GridType> theDefaultPlusRow, theDefaultMinusRow;

public:
  explicit CFEHybridMatrix ( const GridType &Grid ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), theDefaultPlusRow ( Grid ), theDefaultMinusRow ( Grid ) {
    init();
  }

  ~CFEHybridMatrix ( ) {
    for ( int i = 0; i < this->_rows.size(); ++i ) {
      if ( this->isExplicitRow ( i ) ) {
        delete ( this->_rows[i] );
      }
      this->_rows[i] = NULL;
    }
  }

private:
  virtual void init ( ) {

    // first store sparseRows for all boundary nodes (of bounding box, i. e. cube)

    for ( qc::RectangularBoundaryIterator<qc::QC_3D> fbnit ( this->_grid ); fbnit.notAtEnd(); ++fbnit ) {
      this->_rows.getRef ( *fbnit ) = new aol::SparseRow< DataType >;  // boundary rows must be sparse rows because AffineFE rows have entries outside
      this->_rowTypes.set ( *fbnit, this->ExplicitRowBitMask );
      ++ ( this->_numExplicit );
#ifdef VERBOSE
      cerr << "Setting " << *fbnit << " = " << this->_grid.getIndexMapperRef() ( *fbnit ) << " to explicit sparse Row" << endl;
#endif
    }

    // then find rowTypes for all interior nodes:
    // all Dirichlet nodes and their neighbors as well as all near-interface domain nodes must be stored explicitly as AffineFE3D rows.
    // note that cube-boundary nodes were already set and are ignored below.
    qc::BitArray<qc::QC_3D> storeExplicitly ( qc::GridSize<qc::QC_3D> ( this->_grid ) );
    storeExplicitly.setAll ( false );

    if ( this->_grid.isCoarsened() ) {

        storeExplicitly.setAll ( true );

    } else {

      switch ( this->_grid.getNumStructures() ) {

      case 1: // standard case
        for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
          // explicit storage for all nodes that have a Dirichlet neighbor or a neighbor on the opposite side of the interface

          const DataType strucVal = this->_grid.getStructureRef ( 0 ).getValue ( *fnit );

          for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
            qc::CoordType nbpos = ( *fnit ) + tpcfe::CFETopoLookup::_tetNbOffsets[i];

            if ( this->_grid.isAdmissibleNode ( nbpos ) ) {
              if ( this->_grid.isDirichletNode ( nbpos )  ||  ( strucVal * this->_grid.getStructureRef ( 0 ).getValue ( nbpos ) < 0 ) )
                storeExplicitly.set ( *fnit, true );
            }
          }
        }
        break;

      default:
        throw aol::Exception ( "CFEHybridMatrix (CFE_LIEHR): illegal number of structures", __FILE__, __LINE__ );
      }

    }


    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
      if ( this->_rows[ rowIndex ] == NULL ) {
        // only treat those not set yet (i. e. only treat cube-interior)

        if ( storeExplicitly[ rowIndex ] ) {

          // all nodes marked for explicit storage: store explicit AffineFE3D15Row

          this->_rows[ rowIndex ] =  new tpcfe::AffineFE3D15Row<GridType> ( this->_grid );
          this->_rowTypes[ rowIndex ] = this->ExplicitRowBitMask;
          ++ ( this->_numExplicit );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to explicit AffineFE3D15 Row" << endl;
#endif

        } else if ( this->_grid.getStructureRef ( 0 ).getValue ( *fnit ) < 0 ) {

          this->_rows[ rowIndex ] = &theDefaultMinusRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 0; // 0 for minus
          ++ ( this->_numDefault[0] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultMinus Row" << endl;
#endif

        } else if ( this->_grid.getStructureRef ( 0 ).getValue ( *fnit ) > 0 ) {

          this->_rows[ rowIndex ] = &theDefaultPlusRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 1; // 1 for plus
          ++ ( this->_numDefault[1] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultPlus Row" << endl;
#endif

        } else {
          throw aol::Exception ( "Strucutre value == 0. This should not happen", __FILE__, __LINE__ );
        }

      }
    }

    // detect TheDefaultRows:
    for ( unsigned char defaultType = 0; defaultType < 2; ++defaultType ) {
      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
        if ( this->_rowTypes[ rowIndex ] == ( this->DefaultRowBitMask | defaultType ) ) {
          this->_rowTypes[ rowIndex ] |= this->SavedRowBitMask;
          break;
        }
      }
    }

#ifdef VERBOSE
    for ( int i = 0; i < this->getNumRows(); ++i )
      cerr << "Row " << aol::intFormat ( i ) << " has rowType " << static_cast<int> ( this->_rowTypes[i] ) << "; " << this->_rows[i] << endl;
#endif
  }

  void setRowToZero ( const int Row ) {
#ifdef BOUNDS_CHECK
    if ( Row < 0 && Row >= this->_numCols ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEHybridMatrix (CFE_LIEHR) setRowToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    if ( this->isSavedRow ( Row ) ) {
      this->_rows[ Row ]->setZero();
    } else {
      cerr << Row << endl;
      throw aol::Exception ( "CFEHybridMatrix (CFE_LIEHR): Row may not be modified", __FILE__, __LINE__ );
    }
  }

  void setColToZero ( const int Col ) {
#ifdef BOUNDS_CHECK
    if ( Col < 0 && Col >= this->_numRows ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEHybridMatrix (CFE_LIEHR)::setColToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
      int row = Col + this->_grid.getCfeGridIndexOffset ( k );
      if ( row >= 0 && row < this->_numCols ) { // else neighbor does not exist: do nothing.
        this->set ( row, Col, aol::NumberTrait<DataType>::zero ); // we rely on row doing the correct thing when setting to zero.
      }
    }
  }

  // end of class CFEHybridMatrix<DataType, CFE_LIEHR>
};


/** Common basis class for CFE_TPOS and CFE_TPOSELAST. Constructor is protected so it cannot be used on its own
 */
template< typename GridType >
class CFETPOSHybridMatrix : public CFEHybridMatrixBase< GridType > {
public:
  typedef typename GridType::RealType DataType;

  virtual void printStatistics ( ) const {
    int64_t weRequire = this->_grid.getNumberOfNodes() * ( sizeof ( typename CFEHybridMatrixBase<GridType>::RowTypeType )  + sizeof ( void* ) );
    for ( int i = 0; i < this->_rows.size(); ++i )
      if ( this->isSavedRow ( i ) )
        weRequire += this->_rows[i]->numStoredEntries() * ( sizeof ( int ) + sizeof ( DataType ) );
    // + O(1)
    cerr << this->_numExplicit << " rows explicitely stored, " << this->_numDefault.sum() << " default rows, " << endl;
    cerr << "Memory needed: approximately " << weRequire << " bytes, " << ( weRequire >> 20 ) << " MiB." << endl;
  }

protected:
  tpcfe::AffineFE3D15Row<GridType> _theDefaultMinusRow, _theDefaultPlusRow;

  explicit CFETPOSHybridMatrix ( const GridType &Grid ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), _theDefaultMinusRow ( Grid ), _theDefaultPlusRow ( Grid ) {
    init();
  }

  ~CFETPOSHybridMatrix ( ) {
    for ( int i = 0; i < this->_rows.size(); ++i ) {
      if ( this->isExplicitRow ( i ) ) {
        delete ( this->_rows[i] );
      }
      this->_rows[i] = NULL;
    }
  }

private:
  virtual void init ( ) {
    qc::BitArray<qc::QC_3D> typeExplicitSparse ( qc::GridSize<qc::QC_3D> ( this->_grid ) );
    typeExplicitSparse.setAll ( false );

    if ( this->_grid.isCoarsened() ) {

      typeExplicitSparse.setAll ( true );

    } else {

      switch ( this->_grid.getNumStructures() ) {

      case 1: {
        // first store sparseRows for all boundary nodes (of bounding box, i. e. cube)
        for ( qc::RectangularBoundaryIterator<qc::QC_3D> fbnit ( this->_grid ); fbnit.notAtEnd(); ++fbnit ) {
          typeExplicitSparse.set ( *fbnit, true );
        }

        for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
          if ( this->_grid.isNearInterface ( *fnit ) ) {
            for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( *fnit, 2, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) {
              // 2 nbhd is sufficient even though maximal distance to neighboring point is 3.
              // If N is nb with distande 3, there is another near-interface node (on other side of interface) with distance 2 to N.
              typeExplicitSparse.set ( *bit, true );
            }
          }

          // standard neighbors of Dirichlet nodes must also be stored explcitly (if there are jumping-coefficient-additional neighbors, they will be treated by the near-interface loop)
          if ( this->_grid.isDirichletNode ( *fnit ) ) {
            for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
              qc::CoordType nbpos = ( *fnit ) + tpcfe::CFETopoLookup::_tetNbOffsets[i];
              if ( this->_grid.isAdmissibleNode ( nbpos ) ) {
                typeExplicitSparse.set ( nbpos, true );
              }
            }
          }
        }
      } break;

      default:
        throw aol::Exception ( "CFETPOSHybridMatrix: illegal number of structures", __FILE__, __LINE__ );

      }
    }

    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      const int rowIndex = this->_grid.getIndexMapperRef().getGlobalIndex ( *fnit );
      if ( typeExplicitSparse.get ( *fnit ) == true ) {
        this->_rows.getRef ( *fnit ) = new aol::SparseRow< DataType >;
        this->_rowTypes.set ( *fnit, this->ExplicitRowBitMask );
        ++ ( this->_numExplicit );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to explicit Sparse Row" << endl;
#endif
      } else if ( this->_grid.getStructureRef ( 0 ).getValue ( *fnit ) < 0 ) {
        this->_rows[ rowIndex ] = &_theDefaultMinusRow;
        this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 0; // 0 for minus
        ++ ( this->_numDefault[0] );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultMinus Row" << endl;
#endif
      } else {
        this->_rows[ rowIndex ] = &_theDefaultPlusRow;
        this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 1; // 1 for plus
        ++ ( this->_numDefault[1] );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultPlus Row" << endl;
#endif
      }
    }

    for ( unsigned char defaultType = 0; defaultType < 2; ++defaultType ) {
      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
        if ( this->_rowTypes[ rowIndex ] == ( this->DefaultRowBitMask | defaultType ) ) {
          this->_rowTypes[ rowIndex ] |= this->SavedRowBitMask;
          break;
        }
      }
    }


  }

public:
  void setRowToZero ( const int Row ) {
#ifdef BOUNDS_CHECK
    if ( Row < 0 && Row >= this->_numCols ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFETPOSHybridMatrix: setRowToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif
    this->_rows[ Row ]->setZero();
  }

  void setColToZero ( const int Col ) {
#ifdef BOUNDS_CHECK
    if ( Col < 0 && Col >= this->_numRows ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFETPOSHybridMatrix::setColToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    qc::CoordType pos = this->_grid.getIndexMapperRef().splitGlobalIndex ( Col );
    for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( pos, 3, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) {
      const int row = this->_grid.getIndexMapperRef().getGlobalIndex ( *bit );
      this->set ( row, Col, aol::NumberTrait<DataType>::zero );
    }
  }

}; // end class CFETPOSHybridMatrix


/** Warning: first version, not efficient yet!
 */
template < typename DataType >
class CFEHybridMatrix< CFEGrid<DataType, CFE_TPOS > > : public CFETPOSHybridMatrix< CFEGrid<DataType, CFE_TPOS > > {
public:
  explicit CFEHybridMatrix ( const CFEGrid<DataType, CFE_TPOS > &Grid ) : CFETPOSHybridMatrix< CFEGrid<DataType, CFE_TPOS > > ( Grid ) {
  }
};

template < typename DataType, typename NodalCoeffType >
class CFEHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > : public CFETPOSHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > {
public:
  explicit CFEHybridMatrix ( const CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > &Grid ) : CFETPOSHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > ( Grid ) {
  }
};

template < typename DataType, typename NodalCoeffType >
class CFEHybridMatrix< CFEGrid< DataType, CFE_CDWI_ELAST, NodalCoeffType > > : public CFETPOSHybridMatrix< CFEGrid< DataType, CFE_CDWI_ELAST, NodalCoeffType > > {
public:
  explicit CFEHybridMatrix ( const CFEGrid< DataType, CFE_CDWI_ELAST, NodalCoeffType > &Grid ) : CFETPOSHybridMatrix< CFEGrid< DataType, CFE_CDWI_ELAST, NodalCoeffType > > ( Grid ) {
    if ( Grid.isCoarsened() == false ) {
      throw ( aol::Exception ( "CFEHybridMatrix< DataType, CFE_CDWI_ELAST > so far only works for grids with coarsenedFlag set", __FILE__, __LINE__ ) );
    }
  }
};


// B. based on genBandMatrix (band structure - might faster but wastes memory due to explicitely stored zeroes)

/** "Basis" class for template specialization
 */
template < typename GridType >
class CFEBandMatrix { };

/** CFEBandMatrix for complicated domains has 15-band structure of affine FE in 3D. */
// this probably works exactly in the same way for CFE_CDWI_LIEHR
template < typename _DataType, typename NodalCoeffType >
class CFEBandMatrix< CFEGrid<_DataType, CFE_CD, NodalCoeffType> > : public qc::AffineFEBandMatrix< _DataType, qc::QC_3D > {
public:
  typedef _DataType DataType;

  explicit CFEBandMatrix ( const CFEGrid<DataType, CFE_CD, NodalCoeffType> &grid ) : qc::AffineFEBandMatrix< DataType, qc::QC_3D > ( grid ) {}

  CFEBandMatrix ( const int rows, const int cols ) : qc::AffineFEBandMatrix< DataType, qc::QC_3D > ( rows, cols ) {
    throw aol::Exception ( "CFEBandMatrix( const int rows, const int cols ) constructor is not useful.", __FILE__, __LINE__ );
  }

  void printStatistics ( ) {
  }
};


/** CFEBandMatrix for jumping coefficients (LIEHR extrapolation) has 15-band structure of affine FE in 3D */
template < typename _DataType >
class CFEBandMatrix< CFEGrid<_DataType, CFE_LIEHR> > : public qc::AffineFEBandMatrix< _DataType, qc::QC_3D > {
public:
  typedef _DataType DataType;

  explicit CFEBandMatrix ( const CFEGrid<DataType, CFE_LIEHR> &grid ) : qc::AffineFEBandMatrix< DataType, qc::QC_3D > ( grid ) {}

  CFEBandMatrix ( const int rows, const int cols ) : qc::AffineFEBandMatrix< DataType, qc::QC_3D > ( rows, cols ) {
    throw aol::Exception ( "CFEBandMatrix( const int rows, const int cols ) constructor is not useful.", __FILE__, __LINE__ );
  }

  void printStatistics ( ) {}
};


/** Common basis class for CFE_TPOS and CFE_TPOSELAST. Constructor is protected so it cannot be used on its own
 */
template< typename GridType >
class CFETPOSBandMatrix : public aol::GenBandMatrix< typename GridType::RealType > {
public:
  typedef typename GridType::RealType DataType;

protected:
  explicit CFETPOSBandMatrix ( const GridType &grid ) : aol::GenBandMatrix<DataType> ( grid.getNumberOfNodes(), grid.getNumberOfNodes() ) {
    if ( grid.getGridDepth() < 3 ) {
      throw aol::Exception ( "CFEBandMatrix<DataType, CFE_TPOS>: grid level too small. Will not work and would not be useful.", __FILE__, __LINE__ );
    }

    aol::Vector<int> offsetVector ( aol::Cub ( 7 ) );
    short i = 0;
    for ( qc::RectangularIterator<qc::QC_3D> lit ( qc::CoordType ( -3, -3, -3 ), qc::CoordType ( 4, 4, 4 ) ); lit.notAtEnd(); ++lit, ++i ) {
      offsetVector[i] = qc::ILexCombine3 ( ( *lit ) [0], ( *lit ) [1], ( *lit ) [2], grid.getNumX(), grid.getNumY() );
    }

    this->reallocate ( grid.getNumberOfNodes(), grid.getNumberOfNodes(), offsetVector );
  }


  CFETPOSBandMatrix ( const int rows, const int cols ) : aol::GenBandMatrix< DataType > ( rows, cols ) {
    throw aol::Exception ( "CFEBandMatrix( const int rows, const int cols ) constructor is not useful.", __FILE__, __LINE__ );
  }

public:
  void printStatistics ( ) { }

  //! makeRowEntries, here only copying nonzero entries (this makes sense as we have lots of zero entries stored)
  void makeRowEntries ( std::vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.clear();
    vec.reserve ( this->_endApply[ RowNum ] - this->_startApply[ RowNum ] );
    for ( int j = this->_startApply[RowNum]; j < this->_endApply[RowNum]; ++j ) {
      const DataType value = this->_pData[ this->map_index ( j, RowNum ) ];
      if ( value != aol::NumberTrait<DataType>::zero ) {
        vec.push_back ( typename aol::Row<DataType>::RowEntry ( RowNum + this->_localToGlobal[j], value ) );
      }
    }
  }

}; // end class CFETPOSBandMatrix


/** CFEBandMatrix for jumping coefficients (TPOS extrapolation) has at most 89-band structure on the level on which it is assembled;
 *  coarsening can produce bigger band structure. Thus we store the generous upper bound of 7^3 bands. This is probably a waste of memory
 *  and also time during apply, but column-wise operations are faster than with an aol::SparseMatrix.
 *  \author Schwen
 */
template< typename DataType >
class CFEBandMatrix< CFEGrid<DataType, CFE_TPOS> > : public CFETPOSBandMatrix< CFEGrid<DataType, CFE_TPOS> > {
public:
  explicit CFEBandMatrix ( const CFEGrid<DataType, CFE_TPOS> &grid ) : CFETPOSBandMatrix< CFEGrid<DataType, CFE_TPOS> > ( grid ) {
  }

  CFEBandMatrix ( const int rows, const int cols ) : CFETPOSBandMatrix< CFEGrid<DataType, CFE_TPOS> > ( rows, cols ) {
  }
};

template< typename DataType >
class CFEBandMatrix< CFEGrid<DataType, CFE_CDWI_TPOS> > : public CFETPOSBandMatrix< CFEGrid<DataType, CFE_CDWI_TPOS> > {
public:
  explicit CFEBandMatrix ( const CFEGrid<DataType, CFE_CDWI_TPOS> &grid ) : CFETPOSBandMatrix< CFEGrid<DataType, CFE_CDWI_TPOS> > ( grid ) {
  }

  CFEBandMatrix ( const int rows, const int cols ) : CFETPOSBandMatrix< CFEGrid<DataType, CFE_CDWI_TPOS> > ( rows, cols ) {
  }
};

template< typename DataType, typename NodalCoeffType  >
class CFEBandMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > : public CFETPOSBandMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > {
public:
  explicit CFEBandMatrix ( const CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > &grid ) : CFETPOSBandMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > ( grid ) {
  }

  CFEBandMatrix ( const int rows, const int cols ) : CFETPOSBandMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > ( rows, cols ) {
  }
};



// C. matrix classes for periodic boundary values with additional (remapped) entries, makes sense only based on hybrid matrix

/** "Basis" class for template specialization
 */
template < typename GridType >
class CFEPeriodicHybridMatrix : public CFEHybridMatrixBase< GridType > { };

/** Memory-efficient matrix class for CFE (CD complicated domain) computations storing only one instance of default rows for periodic boundary conditions on all cube faces
 *  An important property is that domain and Dirichlet nodes must be set in the grid before instantiating the matrix and must not be changed afterwards.
 *  ATTENTION: use this matrix only if you really know what you are doing and compile with --IKnowWhatIAmDoing :-)
 */
template < typename DataType >
class CFEPeriodicHybridMatrix< CFEGrid<DataType, CFE_CD > > : public CFEHybridMatrixBase< CFEGrid<DataType, CFE_CD > > {
protected:
  typedef CFEGrid<DataType, CFE_CD > GridType;
  aol::DiagonalRow<DataType>       theDefaultExteriorRow;
  tpcfe::AffineFE3D15Row<GridType> theDefaultInteriorRow;

public:
  explicit CFEPeriodicHybridMatrix ( const GridType &Grid ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), theDefaultExteriorRow ( ), theDefaultInteriorRow ( Grid ) {
    init();
  }

  CFEPeriodicHybridMatrix ( const GridType &Grid, const qc::BitArray<qc::QC_3D> &InteriorRows, const qc::BitArray<qc::QC_3D> &ExteriorRows, const qc::BitArray<qc::QC_3D> &ExplicitRows ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), theDefaultExteriorRow ( ), theDefaultInteriorRow ( Grid ) {
    init ( InteriorRows, ExteriorRows, ExplicitRows );
  }

  ~CFEPeriodicHybridMatrix ( ) {
    for ( int i = 0; i < this->_rows.size(); ++i ) {
      if ( this->isExplicitRow ( i ) ) {
        delete ( this->_rows[i] );
      }
      this->_rows[i] = NULL;
    }
  }

  virtual void printStatistics ( ) const {
    const int64_t we_require = this->_numExplicit             * ( NUM_NEIGHBORS * sizeof ( DataType ) + sizeof ( void* ) + sizeof ( int ) ) +  this->_grid.getNumberOfNodes() * ( sizeof ( typename CFEHybridMatrixBase<GridType>::RowTypeType )  + sizeof ( void* ) ); // + O(1)
    const int64_t fs_require = this->_grid.getNumberOfNodes() * ( NUM_NEIGHBORS * sizeof ( DataType ) + sizeof ( void* ) + sizeof ( int ) ); // + O(1)
    const float saving     = 100 * ( 1.0 - static_cast<double> ( we_require ) / static_cast<double> ( fs_require ) );

    cerr << this->_numExplicit << " rows explicitely stored, " << this->_numDefault.sum() << " default rows, " << endl;
    cerr << "Memory needed: approximately " << we_require << " bytes, " << we_require / 1048576 << " MiB." << endl;
    cerr << "Full storage (15 bands plus periodic entries) would require approximately " << fs_require << " bytes (" << fs_require / 1024 << " KiB, " << fs_require / 1048576 << " MiB), so we save " << saving << " percent." << endl;
  }

private:
  virtual void init ( ) {

    // first store sparseRows for all boundary nodes and those in second-to-last layer (due to periodicity, they will get additional neighbors)
    // \todo think about whether to distinguish between non-domain/non-DOF and other
    const int lastX = this->_grid.getNumX() - 1, lastY = this->_grid.getNumY() - 1, lastZ = this->_grid.getNumZ() - 1;
    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      if ( ( *fnit ) [0] == 0 || ( *fnit ) [0] == lastX - 1 || ( *fnit ) [0] == lastX ||
           ( *fnit ) [1] == 0 || ( *fnit ) [1] == lastY - 1 || ( *fnit ) [1] == lastY ||
           ( *fnit ) [2] == 0 || ( *fnit ) [2] == lastZ - 1 || ( *fnit ) [2] == lastZ    ) {
        this->_rows.getRef ( *fnit ) = new aol::SparseRow< DataType >;  // boundary rows must be sparse rows because AffineFE rows have entries outside
        this->_rowTypes.set ( *fnit, this->ExplicitRowBitMask );
        ++ ( this->_numExplicit );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << this->_grid.getIndexMapperRef() ( *fnit ) << " to explicit sparse Row" << endl;
#endif
      }
    }

    // then find rowTypes for all interior nodes:
    // all Dirichlet nodes and their neighbors as well as all near-interface domain nodes must be stored explicitly as AffineFE3D rows.
    // note that cube-boundary nodes were already set and are ignored below.
    qc::BitArray<qc::QC_3D> storeExplicitly ( qc::GridSize<qc::QC_3D> ( this->_grid ) );
    storeExplicitly.setAll ( false );

    if ( this->_grid.isCoarsened() ) {

      storeExplicitly.setAll ( true );

    } else {

      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );

        if ( ( this->_grid.isDOFNode ( rowIndex ) &&  ! ( this->_grid.isInsideDomain ( *fnit ) ) )   ||   this->_grid.isDirichletNode ( rowIndex ) ) {
          // neighbors of DOF nodes outside the domain and neighbors of Dirichlet nodes must be stored explicitely
          for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
            qc::CoordType nbpos = ( *fnit ) + tpcfe::CFETopoLookup::_tetNbOffsets[i];

            if ( this->_grid.isAdmissibleNode ( nbpos ) ) {
              storeExplicitly.set ( nbpos, true );
            }

          }
        }
      }

    }

    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
      if ( this->_rows[ rowIndex ] == NULL ) {
        // only treat those not set yet (i. e. only treat cube-interior)

        if ( storeExplicitly[ rowIndex ] ) {

          // all nodes marked for explicit storage: store explicit AffineFE3D15Row

          this->_rows[ rowIndex ] =  new tpcfe::AffineFE3D15Row<GridType> ( this->_grid );
          this->_rowTypes[ rowIndex ] = this->ExplicitRowBitMask;
          ++ ( this->_numExplicit );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to explicit AffineFE3D15 Row" << endl;
#endif

        } else if ( this->_grid.isInsideDomain ( *fnit ) ) {

          // cube-interior, domain-interior, far from interface and Dirichlet: DefaultInterior

          this->_rows[ rowIndex ] = &theDefaultInteriorRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 0; // 0 for interior
          ++ ( this->_numDefault[0] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultInterior Row" << endl;
#endif

        } else {

          // cube-interior domain-exterior and far from interface (i. e. non-DOF), far from Dirichlet: DefaultExterior

          this->_rows[ rowIndex ] = &theDefaultExteriorRow;
          this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 1; // 1 for exterior
          ++ ( this->_numDefault[1] );
#ifdef VERBOSE
          cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultExterior Row" << endl;
#endif

        }

      }
    }

    // detect TheDefaultRows:
    for ( unsigned char defaultType = 0; defaultType < 2; ++defaultType ) {
      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
        if ( this->_rowTypes[ rowIndex ] == ( this->DefaultRowBitMask | defaultType ) ) {
          this->_rowTypes[ rowIndex ] |= this->SavedRowBitMask;
          break;
        }
      }
    }

#ifdef VERBOSE
    for ( int i = 0; i < this->getNumRows(); ++i )
      cerr << "Row " << aol::intFormat ( i ) << " has rowType " << static_cast<int> ( this->_rowTypes[i] ) << "; " << this->_rows[i] << endl;
#endif

  }

  void init ( const qc::BitArray<qc::QC_3D> &/*InteriorRows*/, const qc::BitArray<qc::QC_3D> &/*ExteriorRows*/, const qc::BitArray<qc::QC_3D> &/*ExplicitRows*/ ) {
    throw aol::UnimplementedCodeException ( "init ( 3 bitArrays ) not implemented yet", __FILE__, __LINE__ );
  }

public:
  void setRowToZero ( const int Row ) {
#ifdef BOUNDS_CHECK
    if ( Row < 0 && Row >= this->_numCols ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEPeriodicHybridMatrix (CFE_CD): setRowToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    if ( this->isSavedRow ( Row ) ) {
      this->_rows[ Row ]->setZero();
    } else if ( this->isDefaultRow ( Row ) && this->whichDefaultRow ( Row ) == 1 ) {
      // do nothing
    } else {
      cerr << Row << endl;
      throw aol::Exception ( "CFEPeriodicHybridMatrix (CFE_CD): Row may not be modified", __FILE__, __LINE__ );
    }

  }

  void setColToZero ( const int Col ) {
#ifdef BOUNDS_CHECK
    if ( Col < 0 && Col >= this->_numRows ) {
      throw aol::OutOfBoundsException ( "tpcfe::CFEPeriodicHybridMatrix (CFE_CD)::setColToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
      int row = Col + this->_grid.getCfeGridIndexOffset ( k );
      if ( row >= 0 && row < this->_numCols ) { // else neighbor does not exist: do nothing.
        this->set ( row, Col, aol::NumberTrait<DataType>::zero ); // we rely on row doing the correct thing when setting to zero.
      }
    }

    // if at periodic face, also need to consider neighbors of facing present nodes
    qc::CoordType pos, ppos;
    this->_grid.getIndexMapperRef().splitGlobalIndex ( Col, pos );
    ppos = pos;
    const aol::Vec3<int> si = this->_grid.getSize();

    for ( int c = 0; c < 3; ++c ) {
      if ( pos[c] == 0 ) {
        ppos[c] = si[c]-1;
      } else if ( pos[c] == si[c]-1 ) {
        ppos[c] = 0;
      } else {
        continue;
      }

      if ( pos != ppos ) { // just to make sure ...
        for ( short k = 0; k < NUM_NEIGHBORS; ++k ) {
          int row = this->_grid.getIndexMapperRef().getGlobalIndex ( ppos ) + this->_grid.getCfeGridIndexOffset ( k );
          if ( row >= 0 && row < this->_numCols ) { // else neighbor does not exist: do nothing.
            this->set ( row, Col, aol::NumberTrait<DataType>::zero ); // we rely on row doing the correct thing when setting to zero.
            if ( this->get ( row, Col ) != aol::NumberTrait<DataType>::zero && ! ( this->isSavedRow ( row ) || ( this->isDefaultRow ( row ) && this->whichDefaultRow ( row ) == 1 ) ) )
              cerr << "setColToZero: cannot set entry " << Col << "in row " << row << endl;
          }
        }
      }

    }

  }

  // end of class CFEPeriodicHybridMatrix<DataType, CFE_CD>
};


template < typename GridType >
class CFETPOSPeriodicHybridMatrix : public CFEHybridMatrixBase< GridType > {
public:
  typedef typename GridType::RealType DataType;

protected:
  tpcfe::AffineFE3D15Row<GridType> _theDefaultMinusRow, _theDefaultPlusRow;

  explicit CFETPOSPeriodicHybridMatrix ( const GridType &Grid ) : CFEHybridMatrixBase< GridType > ( Grid, 2 ), _theDefaultMinusRow ( Grid ), _theDefaultPlusRow ( Grid ) {
    init();
  }

  ~CFETPOSPeriodicHybridMatrix ( ) {
    for ( int i = 0; i < this->_rows.size(); ++i ) {
      if ( this->isExplicitRow ( i ) ) {
        delete ( this->_rows[i] );
      }
      this->_rows[i] = NULL;
    }
  }

private:
  virtual void init ( ) {
    qc::BitArray<qc::QC_3D> typeExplicitSparse ( qc::GridSize<qc::QC_3D> ( this->_grid ) );
    typeExplicitSparse.setAll ( false );

    if ( this->_grid.isCoarsened() ) {

      typeExplicitSparse.setAll ( true );

    } else {

      switch ( this->_grid.getNumStructures() ) {

      case 1: {
        const aol::Vec3<int> upTo ( 3, 3, 3 ), stFr ( this->_grid.getNumX() - 4, this->_grid.getNumY() - 4, this->_grid.getNumZ() - 4 );
        for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {

          // store sparse rows for all nodes possibly "affected" by the boundary, i.e. distance <= 3 from boundary:
          if ( ( *fnit ) [0] <= upTo[0] || ( *fnit ) [0] >= stFr[0] ||
               ( *fnit ) [1] <= upTo[1] || ( *fnit ) [1] >= stFr[1] ||
               ( *fnit ) [2] <= upTo[2] || ( *fnit ) [2] >= stFr[2] ) {
            typeExplicitSparse.set ( *fnit, true );

            continue; // remaining possibilities inside for loop not relevant
          }

          // near the interface, have extended supports
          if ( this->_grid.isNearInterface ( *fnit ) ) {
            for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( *fnit, 2, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) {
              typeExplicitSparse.set ( *bit, true );
            }
          }

          // standard neighbors of Dirichlet nodes must also be stored explcitly (if there are jumping-coefficient-additional neighbors, they will be treated by the near-interface loop)
          if ( this->_grid.isDirichletNode ( *fnit ) ) {
            for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
              qc::CoordType nbpos = ( *fnit ) + tpcfe::CFETopoLookup::_tetNbOffsets[i];
              if ( this->_grid.isAdmissibleNode ( nbpos ) ) {
                typeExplicitSparse.set ( nbpos, true );
              }
            }
          }

        }

      } break;

      default:
        throw aol::Exception ( "CFETPOSPeriodicHybridMatrix: illegal number of structures", __FILE__, __LINE__ );

      }
    }

    for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
      const int rowIndex = this->_grid.getIndexMapperRef().getGlobalIndex ( *fnit );
      if ( typeExplicitSparse.get ( *fnit ) == true ) {
        this->_rows.getRef ( *fnit ) = new aol::SparseRow< DataType >;
        this->_rowTypes.set ( *fnit, this->ExplicitRowBitMask );
        ++ ( this->_numExplicit );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to explicit Sparse Row" << endl;
#endif
      } else if ( this->_grid.getStructureRef ( 0 ).getValue ( *fnit ) < 0 ) {
        this->_rows[ rowIndex ] = &_theDefaultMinusRow;
        this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 0; // 0 for minus
        ++ ( this->_numDefault[0] );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultMinus Row" << endl;
#endif
      } else {
        this->_rows[ rowIndex ] = &_theDefaultPlusRow;
        this->_rowTypes[ rowIndex ] = this->DefaultRowBitMask | 1; // 1 for plus
        ++ ( this->_numDefault[1] );
#ifdef VERBOSE
        cerr << "Setting " << *fnit << " = " << rowIndex << " to DefaultPlus Row" << endl;
#endif
      }
    }

    for ( unsigned char defaultType = 0; defaultType < 2; ++defaultType ) {
      for ( qc::RectangularIterator<qc::QC_3D> fnit ( this->_grid ); fnit.notAtEnd(); ++fnit ) {
        const int rowIndex = this->_grid.getIndexMapperRef() ( *fnit );
        if ( this->_rowTypes[ rowIndex ] == ( this->DefaultRowBitMask | defaultType ) ) {
          this->_rowTypes[ rowIndex ] |= this->SavedRowBitMask;
          break;
        }
      }
    }


  }

public:
  void setRowToZero ( const int Row ) {
#ifdef BOUNDS_CHECK
    if ( Row < 0 && Row >= this->_numCols ) {
      throw aol::OutOfBoundsException ( "CFETPOSPeriodicHybridMatrix::setRowToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif
    this->_rows[ Row ]->setZero();
  }

  void setColToZero ( const int Col ) {
#ifdef BOUNDS_CHECK
    if ( Col < 0 && Col >= this->_numRows ) {
      throw aol::OutOfBoundsException ( "CFETPOSPeriodicHybridMatrix::setColToZero: index out of bounds", __FILE__, __LINE__ );
    }
#endif

    qc::CoordType pos = this->_grid.getIndexMapperRef().splitGlobalIndex ( Col );
    for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( pos, 3, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) {  // here, need 3-neighborhood in all directions
      const int row = this->_grid.getIndexMapperRef().getGlobalIndex ( *bit );
      this->set ( row, Col, aol::NumberTrait<DataType>::zero );
    }

    // if at periodic face, also need to consider neighbors of facing present nodes
    qc::CoordType ppos;

    ppos = pos;
    const aol::Vec3<int> si = this->_grid.getSize();

    for ( short c = 0; c < 3; ++c ) {
      if ( pos[c] == 0 ) {
        ppos[c] = si[c] - 1;
      } else if ( pos[c] == si[c] - 1 ) {
        ppos[c] = 0;
      } else {
        continue;
      }

      if ( pos != ppos ) { // just to make sure ...
        for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( ppos, 3, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) {  // again, need 3-neighborhood in all directions
          const int row = this->_grid.getIndexMapperRef().getGlobalIndex ( *bit );
          this->set ( row, Col, aol::NumberTrait<DataType>::zero );
        }

      }
    }

  }

  // end of class CFETPOSPeriodicHybridMatrix
};

template < typename DataType >
class CFEPeriodicHybridMatrix< CFEGrid<DataType, CFE_TPOS > > : public CFETPOSPeriodicHybridMatrix< CFEGrid<DataType, CFE_TPOS > > {
public:
  explicit CFEPeriodicHybridMatrix ( const CFEGrid<DataType, CFE_TPOS > &Grid ) : CFETPOSPeriodicHybridMatrix< CFEGrid<DataType, CFE_TPOS > > ( Grid ) {
  }
};

template < typename DataType, typename NodalCoeffType >
class CFEPeriodicHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > : public CFETPOSPeriodicHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > {
public:
  explicit CFEPeriodicHybridMatrix ( const CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > &Grid ) : CFETPOSPeriodicHybridMatrix< CFEGrid< DataType, CFE_TPOSELAST, NodalCoeffType > > ( Grid ) {
  }
};

} // end namespace

#endif
