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

#ifndef __TPCFESTANDARDOP_H
#define __TPCFESTANDARDOP_H

#include <op.h>
#include <progressBar.h>
#include <quoc.h>

#include <tpCFEElement.h>
#include <tpCFEGrid.h>
#include <tpCFEMatrices.h>


namespace { // nameless
// the following workaround is necessary because depending on constraint type, the constructor must or must not be given the grid.

template< typename MatrixType, tpcfe::ConstraintType CT, typename NodalCoeffType >
class CFENewMatrixForConfiguratorCreator {
public:
  static MatrixType* create ( const tpcfe::CFEGrid< typename MatrixType::DataType, CT, NodalCoeffType > &Grid ) {
    return ( new MatrixType ( Grid ) );
  }
};

template< typename MatrixType, typename NodalCoeffType>
class CFENewMatrixForConfiguratorCreator< MatrixType, tpcfe::CFE_DOMAIN, NodalCoeffType > {
public:
  static MatrixType* create ( const tpcfe::CFEGrid< typename MatrixType::DataType, tpcfe::CFE_DOMAIN, NodalCoeffType > &Grid ) {
    if ( Grid.hasDomain() ) {
      const int numDofs = Grid.largestDofIndex();
      return ( new MatrixType ( numDofs, numDofs ) );
    } else {
      return ( new MatrixType ( Grid ) );
    }
  }
};

} // end of nameless namespace


namespace tpcfe {

/** Utility class that is used for on-the-fly multiplication of CFE-operators.
 */
template < typename RealType >
class CFEOnTheFlyMultiplicator {
protected:
  const aol::Vector < RealType > &_arg;
  aol::Vector < RealType > &_dest;
public:
  CFEOnTheFlyMultiplicator ( const aol::Vector < RealType > &Arg,
                             aol::Vector < RealType > &Dest ) : _arg ( Arg ), _dest ( Dest ) {}

  void processContribution ( const int row, const int col, const RealType value ) {
    _dest[row] += value * _arg[col];
  }

  const char* getProgressBarText ( ) const {
    return ( "Applying" );
  }
};

/** Utility class that is used for the assembly of CFE-matrices
 */
template < typename RealType, typename MatrixType, ConstraintType CT >
class CFEMatrixAssembler {
protected:
  MatrixType &_mat;

public:
  template< typename GridType >
  CFEMatrixAssembler ( MatrixType &Mat, const GridType & ) : _mat ( Mat ) {}

  void processContribution ( const int row, const int col, const RealType value ) {
    if ( value != aol::NumberTrait<RealType>::zero )
      _mat.add ( row, col, value );
  }

  const char* getProgressBarText ( ) const {
    return ( "Assembling" );
  }

};

/** Utility class that is used for the assembly of CFE-matrices in
 *  the case of the CFE_DOMAIN mode, where a redistribution of the degrees
 *  of freedom has been done.
 */
template < typename RealType, typename MatrixType >
class CFEMatrixAssembler< RealType, MatrixType, CFE_DOMAIN > {
protected:
  typedef tpcfe::CFEGrid<RealType, CFE_DOMAIN> GridType;

  MatrixType  &_mat;
  const GridType    &_grid;

public:
  CFEMatrixAssembler ( MatrixType &Mat, const GridType &Grid ) : _mat ( Mat ), _grid ( Grid ) {}

  void processContribution ( const int row, const int col, const RealType value ) {
    int rowIndex, colIndex;

    // if the given a row or col index is negative, the dof is a virtual node
    // the global index is obtained by removing the "-" sign
    // for all other dofs the global index is stored in _dofMap
    if ( row < 0 ) rowIndex = -row;
    else rowIndex = _grid.remappedDof ( row );
    if ( col < 0 ) colIndex = -col;
    else colIndex = _grid.remappedDof ( col );

    if ( value != static_cast<RealType> ( 0.0 ) ) _mat.add ( rowIndex, colIndex, value );
  }

  const char* getProgressBarText ( ) const {
    return ( "Assembling" );
  }
};


/** CFE configurator. This class provides access to the grid and the matrices
 *  which shall be used by the CFE. Here we have the realization for 3D only
 */
template < typename _GridType, typename _MatrixType >
class CFEConfigurator {
public:
  typedef typename _GridType::RealType                       RealType;
  typedef _GridType                                          GridType;
  typedef _MatrixType                                        MatrixType;

  typedef tpcfe::CFEElement < RealType >                     ElementType;

  typedef aol::Mat<4, 4, RealType>                           TetraMatrixType;

  static const qc::Dimension                                 DimOfWorld = qc::QC_3D;
  static const short                                         maxNumLocalDofs = 8;
  static const tpcfe::ConstraintType                         CT = GridType::CT;

  /** A standard constructor
   */
  explicit CFEConfigurator ( const GridType & Grid ) : _grid ( Grid ) {}

  /** Give out a new matrix
   */
  MatrixType *createNewMatrix () const {
    return ( CFENewMatrixForConfiguratorCreator< MatrixType, CT, typename GridType::NodalCoeffType >::create ( _grid ) );
  }

  /** Return the grid
   */
  const GridType &grid() const { return _grid; }

protected:
  const GridType & _grid;
};

/** Standard Interface for CFE operators. It is an extension of aol::Op and
 *  thus can be plugged into iterative solvers.
 *  New operators should not directly derive from CFELinOpInterface, because for
 *  most purposes the class CFEStandardOp already provides more functionality.
 */
template < typename T_ConfiguratorType, typename Imp >
class CFELinOpInterface : public aol::Op < aol::Vector < typename T_ConfiguratorType::RealType > > {
public:
  typedef T_ConfiguratorType                                   ConfiguratorType;
  typedef aol::Vector < typename ConfiguratorType::RealType >  VectorType;

private:
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MatrixType                MatrixType;
  typedef typename ConfiguratorType::GridType                  GridType;

protected:
  aol::OperatorType      _opType;
  const ConfiguratorType _config;

public:
  bool  _quietMode;

protected:
  mutable MatrixType * _mat;

public:
  /** A constructor which takes a configurator
   */
  explicit CFELinOpInterface ( ConfiguratorType Config, aol::OperatorType OpType = aol::ONTHEFLY ) :
      _opType ( OpType ), _config ( Config ), _quietMode ( false ), _mat ( NULL ) { }

  /** This is a constructor which builds its own configurator with the given grid
   */
  explicit CFELinOpInterface ( const GridType &Grid, aol::OperatorType OpType = aol::ONTHEFLY ) :
      _opType ( OpType ), _config ( Grid ), _quietMode ( false ), _mat ( NULL ) { }

  /** Return the matrix which is stored for this configurator
   */
  MatrixType &getMatrixRef() const {
    if ( _opType == aol::ASSEMBLED && !_mat )
      assembleMatrix();

    if ( _mat ) {
      return *_mat;
    } else {
      throw aol::Exception ( "tpcfe::CFELinOpInterface::getMatrixRef: _mat not set", __FILE__, __LINE__ );
    }
  }

  /** Write an image which shows the non-zero pattern of the associated matrix
   */
  void writeNonZeroPattern ( const char *fileName ) const {
    _mat->writeNonZeroPattern ( fileName );
  }

  /** The standard applyAdd does an "on-the-fly" traversal or a matrix multiplication.
   *  If the matrix has not been assembled yet, it will be assembled now
   */
  void applyAdd ( const aol::Vector <RealType > &Arg,
                  aol::Vector < RealType > &Dest ) const {
    switch ( _opType ) {
      case aol::ONTHEFLY: {
          tpcfe::CFEOnTheFlyMultiplicator<RealType> otfm ( Arg, Dest );
          loopOverElements<tpcfe::CFEOnTheFlyMultiplicator<RealType> > ( otfm );
          break;
        }
      case aol::ASSEMBLED : {
          if ( !_mat ) assembleMatrix();
          _mat->applyAdd ( Arg, Dest );
          break;
        }
      default:
        throw aol::UnimplementedCodeException ( "CFELinOpInterface::applyAdd: unknown opType", __FILE__, __LINE__ );
    }
  }

  /** Compute the contribution of the given element.
   *  This is a standard method which calls @see tpcfe::CFEStandardOp::addSubTetraContrib()
   *  or @see tpcfe::CFEStandardOp::computeNonInterfacedHexaContrib()
   *  on the derived subclasses depending on whether this element is
   *  interfaced or not.
   *  What is done with the contribution is determined by the given instance of
   *  ContribProcessorType. In the standard cases this
   *  will be either CFEOnTheFlyMultiplicator or CFEMatrixAssembler
   */
  template <typename ContribProcessorType>
  void computeElementContribution ( tpcfe::CFEElement<RealType> &e, ContribProcessorType &cpt ) const {

    asImp().preProcessLocalContribution ( e, cpt );

    if ( ! ( e.cfeType().representsInterfaced() ) ) {
      // Element is not interfaced
      asImp().computeNonInterfacedHexaContrib ( e, cpt );
    } else {
      // Prepare all data needed for this element
      e.computeAssembleData ( _config.grid() );

      // Iterate over tetrahedra of this element
      for ( tpcfe::CFETetraInElementIterator<RealType> it ( e ); it.notAtEnd(); ++it ) {
        switch ( it->getVolType() ) {
          case CFETopoTetra::P:
          case CFETopoTetra::N:
            asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::PPPN:
          case CFETopoTetra::NNNP:
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::NNNPPP:
          case CFETopoTetra::PPPNNN:
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          default:
            break;
        }
      }
    }
    asImp().postProcessLocalContribution ( e, cpt );
  }

  /** Compute the contribution of the given element. Works exactly as
   *  computeElementContribution, but considers elements with the correct
   *  sign only
   */
  template <typename ContribProcessorType>
  void computeSignElementContribution ( tpcfe::CFEElement<RealType> &e,
                                        ContribProcessorType &cpt,
                                        const signed char sign = -1             ) const {

    const unsigned char pureCFEType = e.cfeType()._pureType;

    if ( ( sign == -1 && pureCFEType == 0x00 ) || ( sign == + 1 && pureCFEType == 0xFF ) ) return;

    asImp().preProcessLocalContribution ( e, cpt );

    if ( ( sign == + 1 && pureCFEType == 0x00 ) ||
         ( sign == -1 && pureCFEType == 0xFF ) ) {
      // Element is not interfaced
      asImp().computeNonInterfacedHexaContrib ( e, cpt );
    } else {
      // Prepare all data needed for this element
      e.computeAssembleData ( _config.grid() );
      // Iterate over tetrahedra of this element
      for ( tpcfe::CFETetraInElementIterator<RealType> it ( e ); it.notAtEnd(); ++it ) {
        switch ( it->getVolType() ) {
          case CFETopoTetra::P:
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::N:
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::PPPN:
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::NNNP:
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::NNNPPP:
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          case CFETopoTetra::PPPNNN:
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == + 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt ); ++it;
            if ( sign == - 1 ) asImp().addSubTetraContrib ( e, *it, cpt );
            break;
          default:
            break;
        }
      }
    }
    asImp().postProcessLocalContribution ( e, cpt );
  }


  /** Assemble the matrix into mat
   */
  template < typename MatrixType >
  void assembleAddMatrix ( MatrixType & mat ) const {
    tpcfe::CFEMatrixAssembler<RealType, MatrixType, ConfiguratorType::CT>  assembler ( mat, _config.grid() );
    loopOverElements<tpcfe::CFEMatrixAssembler<RealType, MatrixType, ConfiguratorType::CT> > ( assembler );
  }


private:
  template <typename ASSEMBLER_TYPE>
  void loopOverElements ( ASSEMBLER_TYPE &ma ) const {
    aol::ProgressBar<> pb ( ma.getProgressBarText() );
    pb.start ( _config.grid().getNumberOfElements() );

    if ( _config.grid().hasDomain() ) { // the grid has a domain set

      const qc::GridSize<qc::QC_3D> gridSize ( _config.grid() );
      for ( typename GridType::FullElementIterator it ( _config.grid() ); it.notAtEnd(); ++it ) {
        if ( !_quietMode ) pb++;
        tpcfe::CFEElement<RealType> el ( *it, gridSize, _config.grid().getElType ( *it ) );
        if ( _config.grid().getCT() == CFE_CDWI_TPOS || _config.grid().getCT() == CFE_CDWI_LIEHR ) {
          const int structureNo = el.cfeType()._structureNo;

          // integrate not interfaced elements in domain
          if ( ! ( el.cfeType().representsInterfaced() ) ) {
            asImp().computeSignElementContribution ( el, ma, -1 );
          } else {
            // integrate correct part of domain-cutted elements
            if ( structureNo == MAX_STRUCT_ID ) {
              asImp().computeSignElementContribution ( el, ma, -1 );
            }
            // integrate interface-cutted elements
            else {
              asImp().computeElementContribution ( el, ma );
            }
          }
        } else {
          asImp().computeSignElementContribution ( el, ma, -1 );
        }
      }

    } else {

      const qc::GridSize<qc::QC_3D> gridSize ( _config.grid() );
      for ( typename GridType::FullElementIterator it ( _config.grid() ); it.notAtEnd(); ++it ) {
        if ( !_quietMode ) pb++;
        tpcfe::CFEElement<RealType> el ( *it, gridSize, _config.grid().getElType ( *it ) );
        asImp().computeElementContribution ( el, ma );
      }

    }
    if ( !_quietMode ) pb.finish();
  }


public:
  /** Make the (possibly already assembled) matrix invalid by deleting it
   */
  void invalidateMatrix() {
    cerr << "Invalidating matrix " << endl;
    if ( _mat != NULL ) delete _mat;
    _mat = NULL;
  }

  /** Print the matrix object
   */
  void dump ( ostream & ) {
    if ( !_mat ) throw aol::Exception ( "CFELinOpInterface::dump Matrix object does not exist", __FILE__, __LINE__ );
    //          _mat->dump(out);
    _mat->dump();
  }

  /** Print one line of the matrix object
   */
  void dumpLine ( const int which, ostream &out ) {
    if ( !_mat ) throw aol::Exception ( "CFELinOpInterface::dump Matrix object does not exist", __FILE__, __LINE__ );
    _mat->dumpLine ( which, out );
  }

  /** (Re)allocate matrix for this operator
   */
  void allocateMatrix () const {
    if ( _mat ) delete _mat;
    _mat = _config.createNewMatrix ();
  }

  /** Ask the configurator to build a new matrix and then assemble it
   */
  void assembleMatrix () const {
    allocateMatrix();
    assembleAddMatrix ( *_mat );
  }

  /** Checks whether the matrices of the operators are strictly equal
   */
  bool operator== ( CFELinOpInterface &other ) const {
    return ( *_mat ) == other.getMatrixRef();
  }

  /** Checks whether the matrices of the operators are approximately equal, i.e. fails if
   *  the entries differ by more than epsilon
   */
  bool isApproxEqual ( const CFELinOpInterface &other, RealType epsilon ) const {
    return _mat->isApproxEqual ( other.getMatrixRef(), epsilon );
  }

  /** Return the configurator
   */
  const ConfiguratorType &getConfig() const { return _config; }

  /** call makeRowEntries on Matrix, for compatibility with GaussSeidelInverse
   */
  void makeRowEntries ( std::vector<typename aol::Row<RealType>::RowEntry > &vec, const int RowNum ) const {
    if ( _mat ) {
      _mat->makeRowEntries ( vec, RowNum );
    } else {
      throw aol::Exception ( "tpcfe::CFELinOpInterface: Cannot makeRowEntries: _mat not set", __FILE__, __LINE__ );
    }
  }

protected:
  /** Barton-Nackman-Trick avoids virtual functions because of
   *  the static cast
   */
  inline Imp & asImp () {
    return static_cast < Imp & > ( *this );
  }

  inline const Imp & asImp () const {
    return static_cast < const Imp & > ( *this );
  }

  ~CFELinOpInterface() {
    if ( _mat )
      delete ( _mat );
  }

  CFELinOpInterface ( const CFELinOpInterface<T_ConfiguratorType, Imp> &other );

};


/** A base class for standard operators which use precomputed matrices stored in
 *  lookup tables
 *  \ingroup tpCFE
 */
template < typename ConfiguratorType, typename ImpType >
class CFEStandardOp : public CFELinOpInterface < ConfiguratorType, ImpType > {
protected:
  typedef typename ConfiguratorType::RealType                                                          RealType;
  typedef typename ConfiguratorType::GridType                                                          GridType;
  typedef typename ConfiguratorType::MatrixType                                                        MatrixType;

  typedef typename GridType::VNType                                                                    VNType;
  typedef aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType>     LocalMatrixType;
  typedef typename ConfiguratorType::TetraMatrixType                                                   TetraMatrixType;

protected:
  LocalMatrixType  _defaultLocalHexaMatrix;       //!< contains the hexahedron element matrix for a non-interfaced element. This matrix is obtained by adding up all _defaultLocalTetraMatrices for the 6 standard tetrahedra
  TetraMatrixType  _defaultLocalTetraMatrix[6];   //!< stores the matrix for standard tetrahedra which have volume 1./6.

  mutable TetraMatrixType _localTetraMatrix;      //!< stores the local sub-tetra matrix. is filled in prepareTetraMatrix

  const GridType &_grid;

public:
  explicit CFEStandardOp ( const typename ConfiguratorType::GridType & Grid, const aol::OperatorType OpType = aol::ONTHEFLY ) : CFELinOpInterface< ConfiguratorType, ImpType > ( Grid, OpType ), _grid ( Grid ) {}

protected:
  /** Create the default matrices
   */
  void init() {
    this->asImp().createDefaultTetraMatrices();
    this->asImp().buildHexaMatrix();
  }

  /** Assemble the local matrix for a whole hexahedron
   *  from the matrices of the standard tetrahedra
   */
  void buildHexaMatrix() {

    _defaultLocalHexaMatrix.setZero();

    for ( CFETopoTetraIterator it ( CFEType ( -1, 0 ) ); it.notAtEnd(); ++it ) {
      const int parent = it->getParent();

      for ( short i = 0; i < 4; ++i ) {
        const int i1 = ( *it ) ( i, 0 );
        for ( short j = 0; j < 4; ++j ) {
          const int j1 = ( *it ) ( j, 0 );
          _defaultLocalHexaMatrix[i1][j1] += _defaultLocalTetraMatrix[parent][i][j];
        }
      }
    }
  }

public:

  /** Print out the local hexa matrix
   */
  void dumpLocalHexaMatrix ( ostream &out ) const {
    for ( short i = 0; i < 8; ++i ) {
      for ( short j = 0; j < 8; ++j ) {
        out << _defaultLocalHexaMatrix[i][j] << " ";
      }
      out << endl;
    }
  }

  /** Preprocess any computations. The term "local" means Hexahedron not Tetrahedron or Subtetrahedron.
   */
  template <typename ContribProcessorType>
  inline void preProcessLocalContribution ( const CFEElement< RealType>& /*e*/,
                                            ContribProcessorType& /*cpt*/ ) const {}

  /** Postprocess any computations. The term "local" means Hexahedron not Tetrahedron or Subtetrahedron.
   */
  template <typename ContribProcessorType>
  inline void postProcessLocalContribution ( const CFEElement< RealType>& /*e*/,
                                             ContribProcessorType& /*cpt*/ ) const {}

  /** copy the local matrix for an element
   */
  template <typename ContribProcessorType>
  inline void computeNonInterfacedHexaContrib ( const CFEElement< RealType> &e,
                                                ContribProcessorType &cpt ) const {
    for ( short i = 0; i < 8; ++i ) {
      const int i1 = e.globalIndex ( i );
      for ( short j = 0; j < 8; ++j ) {
        const int j1 = e.globalIndex ( j );
        cpt.processContribution ( i1, j1, _defaultLocalHexaMatrix[i][j] );
      }
    }
  }

  /** Return the matrix contribution the coupling i, j would give for standard
   *  basis functions i and j.
   *  In this default implementation it is assumed that the sub-tetra-matrix is
   *  just a scaled version of the parental std-tetra-matrix. Note that this is not
   *  satisfied in all cases. It holds for the standard mass matrices but not for the
   *  standard stiffness matrices!
   */
  inline const RealType getTetraContrib ( const CFEElement< RealType >& /*e*/, const CFETetra<RealType>& /*t*/,
                                          const int /*parent*/, const int i, const int j ) const {
    return _localTetraMatrix[i][j];
  }

  /** Add the constribution of a sub tetrahedron according to the construction of cfe weights
   */
  template <typename ContribProcessorType>
  inline void addSubTetraContrib ( const CFEElement< RealType > &e,
                                   const CFETetra<RealType> & t,
                                   ContribProcessorType &cpt        ) const {

    const short parent = t.getParent();
    this->asImp().preprocessLocalTetraMatrix ( e, t );

    for ( short i = 0; i < 4; ++i ) {
      unsigned char sI = t.isVirtualNode ( i );
      for ( short j = 0; j < 4; ++j ) {
        unsigned char sJ = t.isVirtualNode ( j );
        sJ <<= 1;

        const RealType contrib = this->asImp().getTetraContrib ( e, t, parent, i, j );
        switch ( sI | sJ ) {
          case 0: { // Both nodes are no virtual nodes
            const int i1 = e.globalIndex ( t ( i, 0 ) );
            const int j1 = e.globalIndex ( t ( j, 0 ) );  // Global index for node j

            cpt.processContribution ( i1, j1, contrib );
            break;
          }
          case 1: { // Node i is a virtual node
            const int j1 = e.globalIndex ( t ( j, 0 ) );  // Global index for node j

            // Global indices of nodes defining the edge on which the virtual node lies
            const int en1 = e.globalIndex ( t ( i, 0 ) );
            const int en2 = e.globalIndex ( t ( i, 1 ) );

            // Get the virtual node object from the grid
            const VNType& n = this->_grid.getVirtualNodeRef ( en1, en2 );

            //Ignore basis function if virtual node is a dirichlet node
            if ( n._isDirichlet ) { continue; }

            const int numIConstraints = n.numConstraints();
            for ( int l = 0; l < numIConstraints; ++l ) {
              const RealType coeff = n.weight ( l );
              const int i1 = n.constraint ( l );

              cpt.processContribution ( i1, j1, contrib*coeff );
            }
            //cerr << "1";
            break;
          }
          case 2: { // Node j is a virtual node
            const int i1 = e.globalIndex ( t ( i, 0 ) );

            // Global indices of nodes defining the edge on which the virtual node lies
            const int en1 = e.globalIndex ( t ( j, 0 ) );
            const int en2 = e.globalIndex ( t ( j, 1 ) );

            // Get the virtual node object from the grid
            const VNType& n = this->_grid.getVirtualNodeRef ( en1, en2 );

            //Ignore basis function if virtual node is a dirichlet node
            if ( n._isDirichlet ) { continue; }

            const int numJConstraints = n.numConstraints();
            for ( int l = 0; l < numJConstraints; ++l ) {
              const RealType coeff = n.weight ( l );
              const int j1 = n.constraint ( l );

              cpt.processContribution ( i1, j1, contrib*coeff );
            }
            //cerr << "2";
            break;
          }
          case 3: { // Both nodes are virtual nodes
            // Global indices of nodes defining the edge on which the virtual node lies
            const int eIn1 = e.globalIndex ( t ( i, 0 ) );
            const int eIn2 = e.globalIndex ( t ( i, 1 ) );

            const int eJn1 = e.globalIndex ( t ( j, 0 ) );
            const int eJn2 = e.globalIndex ( t ( j, 1 ) );

            // Get the virtual node objects from the grid
            const VNType& nI = this->_grid.getVirtualNodeRef ( eIn1, eIn2 );
            const VNType& nJ = this->_grid.getVirtualNodeRef ( eJn1, eJn2 );

            //Ignore basis function if virtual node is a dirichlet node
            if ( nI._isDirichlet ) continue;
            if ( nJ._isDirichlet ) continue;

            const int numIConstraints = nI.numConstraints();
            const int numJConstraints = nJ.numConstraints();
            for ( int k = 0; k < numIConstraints; ++k ) {
              const RealType coeffI = nI.weight ( k );
              const int i1 = nI.constraint ( k );
              for ( int l = 0; l < numJConstraints; ++l ) {
                const RealType coeffJ = nJ.weight ( l );
                const int j1 = nJ.constraint ( l );

                cpt.processContribution ( i1, j1, contrib*coeffI*coeffJ );
              }
            }
            //cerr << "3";
            break;
          }
          default:
            cerr << static_cast< int > ( sI | sJ ) << endl;
            throw aol::Exception ( "CFELinOpInterface::addSubTetraContrib: Cannot be!", __FILE__, __LINE__ );
        }
      }
    }
  }

  /** Do any computations which have to be performed before the assembly loop.
   *  E.g. it is more efficient to assemble the local tetra matrix here and
   *  distribute it only in the assembly loop.
   *  Not virtual because we are using Barton-Nackman
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType > &, const CFETetra<RealType> & ) const {}

  /** Do any computations -- whatever that may be -- which have to be performed after the assembly loop.
   *  Not virtual because we are using Barton-Nackman
   */
  void postprocessLocalTetraMatrix ( const CFEElement< RealType > &, const CFETetra<RealType> & ) const {}

  // end class CFEStandardOp
};



/** Composite Finite Element Mass Matrix
 */
template < typename _ConfiguratorType >
class CFEMassOp : public CFEStandardOp < _ConfiguratorType, CFEMassOp < _ConfiguratorType > > {
public:
  typedef CFEStandardOp < _ConfiguratorType, CFEMassOp < _ConfiguratorType > > BaseType;
  typedef typename BaseType::ConfiguratorType ConfiguratorType;
  typedef typename BaseType::VectorType       VectorType;

protected:
  typedef typename ConfiguratorType::RealType RealType;


public:
  explicit CFEMassOp ( const typename ConfiguratorType::GridType & Grid, aol::OperatorType OpType = aol::ONTHEFLY ) :
      CFEStandardOp < ConfiguratorType,  CFEMassOp < ConfiguratorType > > ( Grid, OpType ) {
    this->init();
  }

  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    const RealType stdTetraVolFactor = aol::Cub ( this->_grid.H() );

    for ( short l = 0; l < 6; ++l )
      for ( short i = 0; i < 4; ++i )
        for ( short j = 0; j < 4; ++j )
          this->_defaultLocalTetraMatrix[l][i][j] = tpcfe::CFELookup<RealType>::_localMassMatrixRefTetra[i][j] * stdTetraVolFactor;
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType > &/*e*/, const CFETetra<RealType> & t ) const {
    const RealType volumeFactor = t.getVolume() / tpcfe::CFELookup<RealType>::_stdTetraVolume * aol::Cub ( this->_grid.H() );

    for ( short i = 0; i < 4; ++i )
      for ( short j = 0; j < 4; ++j )
        this->_localTetraMatrix[i][j] = tpcfe::CFELookup<RealType>::_localMassMatrixRefTetra[i][j] * volumeFactor;
  }

  // end class CFEMassOp
};



/** Composite Finite Element Stiffness Matrix
 */
template < typename _ConfiguratorType >
class CFEStiffOp : public CFEStandardOp < _ConfiguratorType, CFEStiffOp < _ConfiguratorType > > {
public:
  typedef CFEStandardOp < _ConfiguratorType, CFEStiffOp < _ConfiguratorType > > BaseType;
  typedef typename BaseType::ConfiguratorType ConfiguratorType;
  typedef typename BaseType::VectorType       VectorType;

protected:
  typedef typename ConfiguratorType::RealType RealType;

public:
  explicit CFEStiffOp ( const typename ConfiguratorType::GridType & Grid, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      CFEStandardOp < ConfiguratorType,  CFEStiffOp< ConfiguratorType > > ( Grid, OpType ) {
    this->init();
  }

  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    const RealType scale = this->_grid.H();

    for ( short l = 0; l < 6; ++l ) {
      for ( short i = 0; i < 4; ++i ) {
        for ( short j = 0; j < 4; ++j ) {
          this->_defaultLocalTetraMatrix[l][i][j] =
            tpcfe::CFELookup<RealType>::_localStiffnessMatrixStdTetra[l][i][j] * scale;
        }
      }
    }
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType > &/*e*/, const CFETetra<RealType> & t ) const {
    // Prepare stiffness matrix
    t.computeInverseTransformation();
    t.computeStiffnessMatrix();

    // Compute scaling factor
    const RealType det = t.determinant();
    const RealType scale = this->_grid.H() / ( 6 * det ); // note: t._lsm is the local stiffness matrix scaled by det^2, hence this factor needs to be divided out

    for ( short i = 0; i < 4; ++i ) {
      for ( short j = 0; j < 4; ++j ) {
        this->_localTetraMatrix[i][j] = t._lsm.get ( i, j ) * scale;
      }
    }

  }

  // end class CFEStiffOp
};


/** This class performs the computation of average (e.g. diffusivity or elasticity) coefficients for one (cubic) CFEElement.
 *  \author Schwen (based on code by Preusser)
 */
template < typename RealType, typename NodalCoeffType >
class CFEWeightProvider {
protected:
  const qc::AArray< NodalCoeffType, qc::QC_3D > &_weights;
  aol::RandomAccessContainer<NodalCoeffType>     _localWeights;
  NodalCoeffType                                 _plusWeight, _minusWeight;
  aol::Vector < signed char >                     _sign;

public:
  CFEWeightProvider ( const qc::AArray< NodalCoeffType, qc::QC_3D > &Weights, const CFEElement<RealType> &Element )
    : _weights ( Weights ), _localWeights ( 8 ), _plusWeight( aol::ZTrait<NodalCoeffType>::zero ), _minusWeight( aol::ZTrait<NodalCoeffType>::zero ), _sign ( 8 ) {

    updateForElement ( Element );
  }

  //! reuse this weightProvider for a different element (of the same grid)
  void updateForElement ( const CFEElement<RealType> &Element ) {
    const CFEType cfeType = Element.cfeType();
    int numMinusNodes = 0, numPlusNodes = 0;
    NodalCoeffType pWeight = aol::ZTrait<NodalCoeffType>::zero, mWeight = aol::ZTrait<NodalCoeffType>::zero;

    for ( int i = 0, bitMask = 1; i < 8; ++i, bitMask <<= 1 ) {
      const int index = Element.globalIndex ( i );
      const NodalCoeffType v = _weights [ index ];
      _sign[i] = ( cfeType._pureType & bitMask ) ? -1 : + 1;
      if ( _sign[i] == -1 ) {
        mWeight += v;
        ++numMinusNodes;
      } else {
        pWeight += v;
        ++numPlusNodes;
      }
      _localWeights[i] = v;
    }

    _minusWeight = mWeight;
    _minusWeight /= static_cast<RealType>(numMinusNodes); // setting to NaN if no MinusNode found is intended here

    _plusWeight = pWeight;
    _plusWeight /= static_cast<RealType>(numPlusNodes);
  }

  // no implementation of destructor, copy constructor necessary
  // assignment operator should fail due to const reference member

  //! Return the weight value at the given global index
  inline NodalCoeffType nodalWeight ( const int globalIndex ) const {
    return ( _weights[ globalIndex ] );
  }

  //! Return the weight value at the given local index for the current element
  inline NodalCoeffType getLocalWeight ( const int localIndex ) const {
    return ( _localWeights[localIndex] );
  }

  //! Return the average weight at the given side of the interface
  inline NodalCoeffType meanWeight ( const signed char sign ) const {
    return ( sign == -1 ) ? _minusWeight : _plusWeight;
  }

  //! Return the average weight on the inside
  inline NodalCoeffType meanMinusWeight() const {
    return ( _minusWeight );
  }

  //! Return the average weight on the outside
  inline NodalCoeffType meanPlusWeight() const {
    return ( _plusWeight );
  }

  //! Return the sign of the level set function at a local node
  inline signed char getSign ( const int index ) const {
    return ( _sign[index] );
  }

protected:
  //! standard constructor cannot work due to const reference member
  CFEWeightProvider();
};


/** Implements a standard operator for bilinear forms with non-constant coefficients
 *  \author Preusser
 */
template < typename ConfiguratorType, typename ImpType >
class CFEWeightedStandardOp : public CFEStandardOp < ConfiguratorType, ImpType > {
protected:
  typedef typename ConfiguratorType::RealType   RealType;

  const qc::AArray<RealType, qc::QC_3D> &_nodalCoeffs;
  mutable CFEWeightProvider< RealType, RealType > _weightProvider; //!< store weight data for current element

public:
  CFEWeightedStandardOp ( const qc::AArray<RealType, qc::QC_3D> &NodalCoeffs,
                          const typename ConfiguratorType::GridType & Grid,
                          const aol::OperatorType OpType = aol::ONTHEFLY ) : CFEStandardOp< ConfiguratorType, ImpType > ( Grid, OpType ), _nodalCoeffs ( NodalCoeffs ), _weightProvider ( NodalCoeffs, CFEElement<RealType>() ) {
    if ( qc::GridSize<qc::QC_3D> ( Grid ) != qc::GridSize<qc::QC_3D>( _nodalCoeffs ) )
      throw aol::Exception ( "CFEWeightedStandardOp::CFEWeightedStandardOp dimensions of grid and weights do not match", __FILE__, __LINE__ );
  }

  /** Clear all matrix entries in the local matrix and gather the local weights
   */
  template <typename ContribProcessorType>
  inline void preProcessLocalContribution ( const CFEElement< RealType> &e,
                                            ContribProcessorType &/*cpt*/ ) const {
    _weightProvider.updateForElement ( e );
  }

  /** Assemble the local matrix for a non-interfaced element
   */
  template <typename ContribProcessorType>
  inline void computeNonInterfacedHexaContrib ( const CFEElement< RealType> &e,
                                                ContribProcessorType &cpt ) const {
    const RealType elWeight = _weightProvider.meanWeight ( e.getNonInterfacedSign() );
    for ( short i = 0; i < 8; ++i ) {
      const int i1 = e.globalIndex ( i );
      for ( short j = 0; j < 8; ++j ) {
        const int j1 = e.globalIndex ( j );
        cpt.processContribution ( i1, j1, elWeight * this->_defaultLocalHexaMatrix[i][j] );
      }
    }
  }

  /** Multiply the local tetra matrix by the weight of the parent node of the tetra
   */

  inline const RealType getTetraContrib ( const CFEElement< RealType >& /*e*/,
                                          const CFETetra<RealType> & t,
                                          const int /*parent*/,
                                          const int i, const int j ) const {

    const RealType weight = _weightProvider.meanWeight ( t.getSign() );
    return weight*this->_localTetraMatrix[i][j];
  }

  // end class CFEWeightedStandardOp
};



/** Provides a mass matrix which works with real cfe-basis functions (with jumping gradients)
 */
template < typename _ConfiguratorType >
class CFEMassOpWI : public CFEWeightedStandardOp < _ConfiguratorType, CFEMassOpWI < _ConfiguratorType > > {
public:
  typedef CFEStandardOp < _ConfiguratorType, CFEMassOpWI < _ConfiguratorType > > BaseType;
  typedef typename BaseType::ConfiguratorType ConfiguratorType;
  typedef typename BaseType::VectorType       VectorType;

protected:
  typedef typename ConfiguratorType::RealType RealType;

public:
  CFEMassOpWI ( const qc::AArray<RealType, qc::QC_3D> &NodalCoeffs,
                const typename ConfiguratorType::GridType &Grid,
                aol::OperatorType OpType = aol::ONTHEFLY ) :
      CFEWeightedStandardOp < ConfiguratorType,  CFEMassOpWI < ConfiguratorType > > ( NodalCoeffs, Grid, OpType ) {
    this->init();
  }


  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    const RealType stdTetraVolFactor = aol::Cub ( this->_grid.H() );

    for ( int l = 0; l < 6; ++l )
      for ( int i = 0; i < 4; ++i )
        for ( int j = 0; j < 4; ++j )
          this->_defaultLocalTetraMatrix[l][i][j] = tpcfe::CFELookup<RealType>::_localMassMatrixRefTetra[i][j] * stdTetraVolFactor;
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType >& /*e*/, const CFETetra<RealType> & t ) const {
    const RealType volumeFactor = t.getVolume() / tpcfe::CFELookup<RealType>::_stdTetraVolume * aol::Cub ( this->_grid.H() );

    for ( int i = 0; i < 4; ++i )
      for ( int j = 0; j < 4; ++j )
        this->_localTetraMatrix[i][j] = tpcfe::CFELookup<RealType>::_localMassMatrixRefTetra[i][j] * volumeFactor;
  }

  // end class CFEMassOpWI
};



/** Class for a stiffness-matrix with a spatially dependent coefficient
 *  Provides assymbly of matrices for bilinear forms like
 *  \f[  (a(x) \nabla\phi_i, \nabla\phi_j)_{0,2,E} \f]
 *  \author Preusser
 */
template < typename _ConfiguratorType >
class CFEStiffOpWI : public CFEWeightedStandardOp < _ConfiguratorType, CFEStiffOpWI < _ConfiguratorType > > {
public:
  typedef CFEStandardOp < _ConfiguratorType, CFEStiffOpWI < _ConfiguratorType > > BaseType;
  typedef typename BaseType::ConfiguratorType ConfiguratorType;
  typedef typename BaseType::VectorType       VectorType;

protected:
  typedef typename ConfiguratorType::RealType RealType;

public:
  CFEStiffOpWI ( const qc::AArray<RealType, qc::QC_3D> &NodalCoeffs,
                 const typename ConfiguratorType::GridType & Grid, const aol::OperatorType OpType = aol::ONTHEFLY ) :
      CFEWeightedStandardOp < ConfiguratorType,  CFEStiffOpWI< ConfiguratorType > > ( NodalCoeffs, Grid, OpType ) {
    this->init();
  }

  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    const RealType scale = this->_grid.H();

    for ( short l = 0; l < 6; ++l )
      for ( short i = 0; i < 4; ++i )
        for ( short j = 0; j < 4; ++j )
          this->_defaultLocalTetraMatrix[l][i][j] = tpcfe::CFELookup<RealType>::_localStiffnessMatrixStdTetra[l][i][j] * scale;
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType >& /*e*/, const CFETetra<RealType> & t ) const {
    // Prepare stiffness matrix
    t.computeInverseTransformation();
    t.computeStiffnessMatrix();

    // Compute scaling factor
    const RealType det  = t.determinant();
    const RealType scale = this->_grid.H() / ( 6 * det ) ;

    for ( short i = 0; i < 4; ++i )
      for ( short j = 0; j < 4; ++j )
        this->_localTetraMatrix[i][j] = t._lsm.get ( i, j ) * scale;
  }

  // end class CFEStiffOpWI
};

}

#endif
