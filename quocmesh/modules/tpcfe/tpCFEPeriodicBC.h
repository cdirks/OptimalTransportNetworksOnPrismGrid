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

#ifndef __TPCFEPERIODICBC_H
#define __TPCFEPERIODICBC_H

#include <periodicBC.h>

#include <tpCFEGrid.h>
#include <tpCFEMultigrid.h>

namespace tpcfe {

/** Handling of periodicity in all three space directions for CFE for complicated domains.
 *  Basis class from which we will derive by template specialization.
 *  \author Schwen
 */
template< typename GridType >
class CFEPeriodicityHandlerBase : public qc::PeriodicityHandlerBase< GridType, typename GridType::RealType, qc::QC_3D > {

protected:
  typedef typename GridType::RealType DataType;

  const aol::Vec3<int> _size;
  qc::BitArray<qc::QC_3D> _presentDOFNodes;    //!< Locations at which data is stored, i. e. either node is DOF or node is facingPresentNode of DOF

public:
  explicit CFEPeriodicityHandlerBase ( const GridType &Grid ) : qc::PeriodicityHandlerBase< GridType, typename GridType::RealType, qc::QC_3D > ( Grid ), _size ( Grid.getSize() ), _presentDOFNodes ( qc::GridSize<qc::QC_3D> ( Grid ) ) {
  }

  // default copy constructor, assignment operator and destructor are correct.

  virtual bool isPeriodicNode ( const qc::CoordType &pos ) const {
    return ( ( pos[0] == ( _size[0] - 1 ) ) || ( pos[1] == ( _size[1] - 1 ) ) || ( pos[2] == ( _size[2] - 1 ) ) );
  }

  virtual qc::CoordType facingPresentNode ( const qc::CoordType &pos ) const {
    qc::CoordType ret ( pos );

    for ( int i = 0; i < 3; ++i ) {
      if ( ret[i] == ( _size[i] - 1 ) ) {
        ret[i] = 0;
      }
    }

    return ( ret );
  }

  bool isPresentDOFNode ( const qc::CoordType &pos ) const {
    return ( _presentDOFNodes.get ( pos ) );
  }

  bool isPresentDOFNode ( const int index ) const {
    return ( _presentDOFNodes.get ( index ) );
  }

  const qc::BitArray<qc::QC_3D>& getPresentDOFNodeMask ( ) const {
    return ( _presentDOFNodes );
  }

  void restrictToPresentDOFs ( aol::Vector<DataType> &Arg ) const {
    for ( int i = 0; i < Arg.size(); ++i )
      if ( ! _presentDOFNodes.get ( i ) )
        Arg[i] = aol::NumberTrait<DataType>::zero;
  }

  void restrictToPresentDOFs ( aol::MultiVector<DataType> &Arg ) const {
    for ( int i = 0; i < Arg.numComponents(); ++i )
      restrictToPresentDOFs ( Arg[i] );
  }

  void restrictNonPresentDOFEntries ( aol::Matrix<DataType> &Arg, const DataType DiagValue = aol::NumberTrait<DataType>::one ) const {
    for ( int i = 0; i < this->_grid.getNumberOfNodes(); ++i ) {
      if ( ! _presentDOFNodes.get ( i ) ) {
        Arg.setRowColToZero ( i );
        Arg.set ( i, i, DiagValue );
      }
    }
  }

  template< typename BlockMatrixType >
  void restrictNonPresentDOFEntries ( BlockMatrixType &Arg ) const {
    for ( int i = 0; i < Arg.getNumRows(); ++i )
      for ( int j = 0; j < Arg.getNumCols(); ++j )
        restrictNonPresentDOFEntries ( Arg.getReference ( i, j ), ( i == j ? aol::NumberTrait<DataType>::one : aol::NumberTrait<DataType>::zero ) );
  }

};


/** Handling of periodicity in all three space directions for CFE for complicated domains.
 *  Basis class for template specialization.
 *  \author Schwen
 */
template< typename GridType >
class CFEPeriodicityHandler;


template< typename DataType >
class CFEPeriodicityHandler< CFEGrid< DataType, CFE_CD > > : public CFEPeriodicityHandlerBase< CFEGrid< DataType, CFE_CD > > {
  typedef CFEGrid< DataType, CFE_CD > GridType;

public:
  explicit CFEPeriodicityHandler ( const GridType &Grid ) : CFEPeriodicityHandlerBase< GridType > ( Grid ) {

    this->_presentDOFNodes.setAll ( false );

    for ( qc::RectangularIterator<qc::QC_3D> bit ( this->_grid ); bit.notAtEnd(); ++bit )
      if ( this->_grid.isDOFNode ( *bit ) )
        this->_presentDOFNodes.set ( this->facingPresentNode ( *bit ), true );

  }

  // default copy constructor, assignment operator and destructor are correct.


  //! faster CFE variant; make sure setDOFMaskFromDirichletAndDomainNodeMask() has been called.
  virtual void periodicallyCollapseMatrix ( aol::Matrix<DataType> &matrix, const DataType diagValue = aol::NumberTrait<DataType>::one ) const {
#ifdef VERBOSE
    aol::ProgressBar<> pb ( "Collapsing matrix" );
    pb.start ( this->_grid.getNumberOfBoundaryNodes() );
#endif
    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bit ( this->_grid ); bit.notAtEnd(); ++bit ) {
#ifdef VERBOSE
      pb++;
#endif
      const qc::CoordType pos ( *bit ), fpos ( this->facingPresentNode ( pos ) );

      if ( pos != fpos ) { // note: fpos need not be DOF!
        const int gP = this->_map ( pos ), gF = this->_map ( fpos );

        // fast eliminate col
        for ( int nb = 0; nb < tpcfe::NUM_NEIGHBORS; ++nb ) {
          const qc::CoordType Npos ( pos + tpcfe::CFELookup<DataType>::_tetNbOffsets[nb] );  // neighbor

          if ( tpcfe::CFEMG_is_inside ( Npos, this->_size ) ) {
            const qc::CoordType fNpos ( this->facingPresentNode ( Npos ) );                        // facing present node of neighbor (this is symmetric; may be useless if Npos is not inside, but then it is not used).

            const int gN = this->_map ( Npos ), gFN = this->_map ( fNpos );

            matrix.add ( gFN, gF, matrix.get ( gN, gP ) );
            matrix.set ( gN, gP, aol::NumberTrait<DataType>::zero );

            matrix.add ( gF, gFN, matrix.get ( gP, gN ) );
            matrix.set ( gP, gN, aol::NumberTrait<DataType>::zero );
          }
        }
        matrix.set ( gP, gP, diagValue );
      }
    }
#ifdef VERBOSE
    pb.finish();
#endif
  }

  //! copies from facingPresentNodes to periodicNodes, then restricts to grid's DOFs
  virtual void extendPeriodicBC ( aol::Vector<DataType> &Arg ) const {
    qc::PeriodicityHandlerBase< GridType, DataType, qc::QC_3D >::extendPeriodicBC ( Arg );
    this->_grid.restrictToDofs ( Arg );
  }

  using qc::PeriodicityHandlerBase< GridType, DataType, qc::QC_3D >::extendPeriodicBC;

};


template< typename GridType >
class CFETPOSPeriodicityHandler : public CFEPeriodicityHandlerBase< GridType > {

protected:
  typedef typename GridType::RealType RealType;

  explicit CFETPOSPeriodicityHandler ( const GridType &Grid ) : CFEPeriodicityHandlerBase< GridType > ( Grid ) {

    this->_presentDOFNodes.setAll ( false );

    for ( qc::RectangularIterator<qc::QC_3D> bit ( this->_grid ); bit.notAtEnd(); ++bit )
      this->_presentDOFNodes.set ( this->facingPresentNode ( *bit ), true );

  }

  // default copy constructor, assignment operator and destructor are correct.


public:
  //! faster CFE variant; make sure setDOFMaskFromDirichletAndDomainNodeMask() has been called.
  virtual void periodicallyCollapseMatrix ( aol::Matrix<RealType> &matrix, const RealType diagValue = aol::NumberTrait<RealType>::one ) const {
#ifdef VERBOSE
    aol::ProgressBar<> pb ( "Collapsing matrix" );
    pb.start ( this->_grid.getNumberOfBoundaryNodes() );
#endif
    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bit ( this->_grid ); bit.notAtEnd(); ++bit ) {
#ifdef VERBOSE
      pb++;
#endif
      const qc::CoordType pos ( *bit ), fpos ( this->facingPresentNode ( pos ) );

      if ( pos != fpos ) { // note: fpos need not be DOF!
        const int gP = this->_map ( pos ), gF = this->_map ( fpos );

        for ( qc::LocalLInfBoxIterator< qc::QC_3D > bit ( pos, 3, aol::Vec3<int> ( 0, 0, 0 ), this->_grid.getSize() ); bit.notAtEnd(); ++bit ) { // 2 or 3 neighborhood?
          const qc::CoordType Npos ( *bit );  // neighbor
          const qc::CoordType fNpos ( this->facingPresentNode ( Npos ) );                        // facing present node of neighbor (this is symmetric; may be useless if Npos is not inside, but then it is not used).

          const int gN = this->_map ( Npos ), gFN = this->_map ( fNpos );

          matrix.add ( gFN, gF, matrix.get ( gN, gP ) );
          matrix.set ( gN, gP, aol::NumberTrait<RealType>::zero );

          matrix.add ( gF, gFN, matrix.get ( gP, gN ) );
          matrix.set ( gP, gN, aol::NumberTrait<RealType>::zero );
        }
        matrix.set ( gP, gP, diagValue );
      }
    }
#ifdef VERBOSE
    pb.finish();
#endif
  }

  //! copies from facingPresentNodes to periodicNodes
  virtual void extendPeriodicBC ( aol::Vector<RealType> &Arg ) const {
    qc::PeriodicityHandlerBase< GridType, RealType, qc::QC_3D >::extendPeriodicBC ( Arg );
  }

  using qc::PeriodicityHandlerBase< GridType, RealType, qc::QC_3D >::extendPeriodicBC;

};


template< typename RealType >
class CFEPeriodicityHandler< CFEGrid< RealType, CFE_TPOS > > : public CFETPOSPeriodicityHandler< CFEGrid< RealType, CFE_TPOS > > {
  typedef CFEGrid< RealType, CFE_TPOS > GridType;
public:
  CFEPeriodicityHandler ( const GridType &Grid ) : CFETPOSPeriodicityHandler < GridType > ( Grid ) {
  }
};


template < typename RealType, typename NodalCoeffType >
class CFEPeriodicityHandler< CFEGrid< RealType, CFE_TPOSELAST, NodalCoeffType > > : public CFETPOSPeriodicityHandler< CFEGrid< RealType, CFE_TPOSELAST, NodalCoeffType > > {
  typedef CFEGrid< RealType, CFE_TPOSELAST, NodalCoeffType > GridType;
public:
  CFEPeriodicityHandler ( const GridType &Grid ) : CFETPOSPeriodicityHandler < GridType > ( Grid ) {
  }
};



/** CFE Prolongation operator respecting periodic boundary conditions
 *  \author Schwen
 */
template< typename GridType >
class CFEPeriodicProlongOp : public aol::BiOp< aol::Vector< typename GridType::RealType > > {

protected:
  typedef typename GridType::RealType       RealType;
  typedef CFEPeriodicityHandler< GridType > PHType;

  const GridType                       &_coarseGrid, &_fineGrid;
  const PHType                         &_coarsePH, &_finePH;
  const int                            _wi_co, _wi_fi;
  aol::OperatorType                    _opType;
  mutable aol::SparseMatrix<RealType>  *_pmat;


public:
  CFEPeriodicProlongOp ( const GridType &coarseGrid, const GridType &fineGrid, const PHType &coarsePH, const PHType &finePH, aol::OperatorType OpType = aol::ONTHEFLY )
    : _coarseGrid ( coarseGrid ), _fineGrid ( fineGrid ), _coarsePH ( coarsePH ), _finePH ( finePH ), _wi_co ( coarseGrid.getWidth() ), _wi_fi ( fineGrid.getWidth() ), _opType ( OpType ), _pmat ( NULL ) {
    if ( _coarseGrid.getGridDepth() != _fineGrid.getGridDepth() - 1 ) {
      throw aol::Exception ( "Cannot construct CFE Prolongation operator: incorrect grid dephts", __FILE__, __LINE__ );
    }

        // CFE_TPOS different

    if ( _opType == aol::ASSEMBLED ) {
      assembleMatrix();
    }

  }

public:
  ~CFEPeriodicProlongOp ( ) {
    if ( _pmat ) {
      delete ( _pmat );
    }
  }

public:
  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    if ( _opType == aol::ASSEMBLED ) {

      _pmat->apply ( Arg, Dest );

    } else {

      if ( ( Arg.size() != _coarseGrid.getNumberOfNodes() ) || Dest.size() != _fineGrid.getNumberOfNodes() ) {
        throw aol::Exception ( "Cannot apply CFE Prolongation operator: incorrect vector size.", __FILE__, __LINE__ );
      }

      Dest.setZero();

      switch ( _coarseGrid.CT ) {
      case CFE_CD :

        for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( _coarseGrid ); cPosIt.notAtEnd(); ++cPosIt ) {
          const qc::CoordType cPos = *cPosIt, cPosP = _coarsePH.facingPresentNode ( cPos );
          if ( _coarsePH.isPresentDOFNode ( cPos ) || _coarsePH.isPresentDOFNode ( cPosP ) ) {
            for ( int i = 0; i < tpcfe::NUM_NEIGHBORS; ++i ) {
              const qc::CoordType fNPos = static_cast<short> ( 2 ) * cPos + CFELookup<RealType>::_tetNbOffsets[i];
              if ( CFEMG_is_inside ( fNPos, _wi_fi ) && _finePH.isPresentDOFNode ( fNPos ) ) {
                Dest.add ( _fineGrid.getIndexMapperRef() ( fNPos ), ( i == NEIGHBOR_SELF ? 1.0 : 0.5 ) * Arg.get ( _coarseGrid.getIndexMapperRef() ( cPosP ) ) );
              }
            }
          }
        }

        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEPeriodicProlongOp::applyAdd is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

      }
    }

  }

  using aol::BiOp< aol::Vector<RealType> >::apply;

  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    throw aol::Exception ( "Are you sure you want to applyAdd the CFEPeriodicProlongOp?", __FILE__, __LINE__ );
    aol::Vector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp< aol::Vector<RealType> >::applyAdd;

protected:
  virtual void assembleMatrix ( ) const {
    if ( _pmat ) delete _pmat ;
    _pmat = new aol::SparseMatrix<RealType> ( _fineGrid.getNumberOfNodes() , _coarseGrid.getNumberOfNodes() );
    assembleMatrix ( *_pmat );
  }

public:
  virtual void assembleMatrix ( aol::SparseMatrix<RealType> &mat ) const {
    mat.setZero();

    switch ( _coarseGrid.CT ) {
    case CFE_CD :
      for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( _coarseGrid ); cPosIt.notAtEnd(); ++cPosIt ) {
        const qc::CoordType cPos = *cPosIt, cPosP = _coarsePH.facingPresentNode ( cPos );
        if ( _coarsePH.isPresentDOFNode ( cPos ) || _coarsePH.isPresentDOFNode ( cPosP ) ) {
          for ( int i = 0; i < tpcfe::NUM_NEIGHBORS; ++i ) {
            const qc::CoordType fNPos = static_cast<short> ( 2 ) * cPos + CFELookup<RealType>::_tetNbOffsets[i];
            if ( CFEMG_is_inside ( fNPos, _wi_fi ) && _finePH.isPresentDOFNode ( fNPos ) ) {
              mat.set ( _fineGrid.getIndexMapperRef() ( fNPos ),
                        _coarseGrid.getIndexMapperRef() ( cPosP ),
                        ( i == NEIGHBOR_SELF ? 1.0 : 0.5 ) );
            }
          }
        }
      }

      break;

    default:
      throw aol::UnimplementedCodeException ( "CFEPeriodicProlongOp::assembleMatrix is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

    }
  }

private:
  CFEPeriodicProlongOp ( );
  CFEPeriodicProlongOp ( const CFEPeriodicProlongOp<GridType>& );
  CFEPeriodicProlongOp<GridType>& operator= ( const CFEPeriodicProlongOp<GridType>& );

  // end class CFEPeriodicProlongOp
};



/** CFE Prolongation operator respecting periodic boundary conditions
 *  \author Schwen
 */
template< typename GridType >
class CFEPeriodicRestrictOp : public aol::BiOp< aol::Vector< typename GridType::RealType > > {

protected:
  typedef typename GridType::RealType       RealType;
  typedef CFEPeriodicityHandler< GridType > PHType;

  const GridType                       &_coarseGrid, &_fineGrid;
  const PHType                         &_coarsePH, &_finePH;
  const int                            _wi_co, _wi_fi;
  aol::OperatorType                    _opType;
  mutable aol::SparseMatrix<RealType>  *_pmat;


public:
  CFEPeriodicRestrictOp ( const GridType &coarseGrid, const GridType &fineGrid, const PHType &coarsePH, const PHType &finePH, aol::OperatorType opType = aol::ONTHEFLY )
    : _coarseGrid ( coarseGrid ), _fineGrid ( fineGrid ), _coarsePH ( coarsePH ), _finePH ( finePH ), _wi_co ( coarseGrid.getWidth() ), _wi_fi ( fineGrid.getWidth() ), _opType ( opType ), _pmat ( NULL ) {

    if ( _coarseGrid.getGridDepth() != _fineGrid.getGridDepth() - 1 ) {
      throw aol::Exception ( "Cannot construct CFE Restriction operator: incorrect grid dephts", __FILE__, __LINE__ );
    }

    // CFE_TPOS different

    if ( _opType == aol::ASSEMBLED ) {
      this->assembleMatrix();
    }

  }

public:
  ~CFEPeriodicRestrictOp ( ) {
    if ( _pmat ) {
      delete ( _pmat );
    }
  }

public:
  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    if ( _opType == aol::ASSEMBLED ) {

      _pmat->apply ( Arg, Dest );

    } else {

      if ( ( Arg.size() != _fineGrid.getNumberOfNodes() ) || Dest.size() != _coarseGrid.getNumberOfNodes() ) {
        throw aol::Exception ( "Cannot apply CFE Restriction operator: incorrect vector size.", __FILE__, __LINE__ );
      }

      Dest.setZero();

      switch ( _coarseGrid.CT ) {
      case CFE_CD :
        for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( _coarseGrid ); cPosIt.notAtEnd(); ++cPosIt ) {
          const qc::CoordType cPos = *cPosIt, cPosP = _coarsePH.facingPresentNode ( cPos );
          if ( _coarsePH.isPresentDOFNode ( cPos ) || _coarsePH.isPresentDOFNode ( cPosP ) ) {
            RealType res = aol::NumberTrait<RealType>::zero;
            for ( int i = 0; i < tpcfe::NUM_NEIGHBORS; ++i ) {
              const qc::CoordType fNPos = static_cast<short> ( 2 ) * cPos + CFELookup<RealType>::_tetNbOffsets[i];
              if ( CFEMG_is_inside ( fNPos, _wi_fi ) && _finePH.isPresentDOFNode ( fNPos ) ) {

                res += ( ( i == NEIGHBOR_SELF ? 1.0 : 0.5 ) *
                         Arg.get ( _fineGrid.getIndexMapperRef() ( fNPos ) ) );
              }
            }
            Dest.add ( _coarseGrid.getIndexMapperRef() ( cPosP ), res );
          }
        }

        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEPeriodicRestrictOp::applyAdd is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

      }

    }

  }

  using aol::BiOp< aol::Vector< RealType > >::apply;


  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    throw aol::Exception ( "Are you sure you want to applyAdd the CFEPeriodicRestrictOp?", __FILE__, __LINE__ );
    aol::Vector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp< aol::Vector< RealType > >::applyAdd;

protected:
  virtual void assembleMatrix ( ) const { // obviously, constness is nonsense. but necessary ...
    if ( _pmat ) delete _pmat ;
    _pmat = new aol::SparseMatrix<RealType> ( _coarseGrid.getNumberOfNodes(), _fineGrid.getNumberOfNodes() );
    assembleMatrix ( *_pmat );
  }

public:
  virtual void assembleMatrix ( aol::SparseMatrix<RealType> &mat ) const {
    mat.setZero();

    switch ( _coarseGrid.CT ) {
    case CFE_CD :
      for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( _coarseGrid ); cPosIt.notAtEnd(); ++cPosIt ) {
        const qc::CoordType cPos = *cPosIt, cPosP = _coarsePH.facingPresentNode ( cPos );
        if ( _coarsePH.isPresentDOFNode ( cPos ) || _coarsePH.isPresentDOFNode ( cPosP ) ) {
          for ( int i = 0; i < tpcfe::NUM_NEIGHBORS; ++i ) {
            const qc::CoordType fNPos = static_cast<short> ( 2 ) * cPos + CFELookup<RealType>::_tetNbOffsets[i];
            if ( CFEMG_is_inside ( fNPos, _wi_fi ) && _finePH.isPresentDOFNode ( fNPos ) ) {
              mat.set ( _coarseGrid.getIndexMapperRef() ( cPosP ),
                        _fineGrid.getIndexMapperRef() ( fNPos ),
                        ( i == NEIGHBOR_SELF ? 1.0 : 0.5 ) );
            }
          }
        }
      }

      break;

    default:
      throw aol::UnimplementedCodeException ( "CFEPeriodicRestrictOp::applyAdd is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

    }
  }

private:
  CFEPeriodicRestrictOp ( );
  CFEPeriodicRestrictOp ( const CFEPeriodicRestrictOp<GridType>& );
  CFEPeriodicRestrictOp<GridType>& operator= ( const CFEPeriodicRestrictOp<GridType>& );

   // end class CFEPeriodicRestrictOp
};



/** CFE Multigrid solver with projection to average-0 constraints for periodic boundary conditions.
 *  DOF/DomainNode masks are managed by the grid (as usual), periodicity is handled via periodicityHandlers.
 *  \attention only works for periodic geometries
 *  \author Schwen
 */
template< class CFEOpType, class CFEMassOpType, mg::CoarseningMode CMode = mg::MATRIXMULT_COARSENING >
class CFEMultigridProjectAvgConstr : public mg::ExplicitOperatorHierarchyMultigrid < typename CFEOpType::VectorType,
                                                                                     typename CFEOpType::ConfiguratorType::MatrixType,
                                                                                     mg::GaussSeidelSelectiveSmoother< typename CFEOpType::VectorType, typename CFEOpType::ConfiguratorType::MatrixType >,
                                                                                     tpcfe::CFEPeriodicRestrictOp< typename CFEOpType::ConfiguratorType::GridType >,
                                                                                     tpcfe::CFEPeriodicProlongOp< typename CFEOpType::ConfiguratorType::GridType >, typename CFEOpType::ConfiguratorType::GridType > {

protected:
  typedef typename CFEOpType::ConfiguratorType::RealType   RealType;
  typedef typename CFEOpType::VectorType                   VectorType;
  typedef typename CFEOpType::ConfiguratorType::MatrixType MatrixType;
  typedef typename CFEOpType::ConfiguratorType::GridType   GridType;
  typedef CFEPeriodicityHandler< GridType >                PHType;

  RealType                                            _relax;
  aol::RandomAccessContainer<VectorType>*             _pConstrVecLev;
  aol::RandomAccessContainer<VectorType>*             _pProjDirsVecLev;
  RealType                                            _projectThreshold;
  std::vector< const PHType* >                        _pPeriodicityHandlers;

public:
  // more options: coarsening stopping criteria!
  CFEMultigridProjectAvgConstr ( const GridType &Grid,
                                 const MatrixType &fineMatrix,
                                 const CFEPeriodicityHandler < GridType > &PeriodicityHandler,
                                 const int ExplicitLevel = 2,
                                 const int PreSmooth = 3, const int PostSmooth = 3,
                                 const RealType Relax = 1.0,
                                 const mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                                 const RealType CycleEps = 1.0e-16,
                                 const int MaxCycle = 100,
                                 const aol::StoppingMode stop = aol::STOPPING_UNSET )
      : mg::ExplicitOperatorHierarchyMultigrid < VectorType,
                                                 MatrixType,
                                                 mg::GaussSeidelSelectiveSmoother< VectorType, MatrixType >,
                                                 tpcfe::CFEPeriodicRestrictOp< GridType >,
                                                 tpcfe::CFEPeriodicProlongOp< GridType >, GridType > ( fineMatrix,
                                                                                                       Grid,
                                                                                                       ExplicitLevel,
                                                                                                       PreSmooth, PostSmooth,
                                                                                                       Relax,
                                                                                                       CycleMode, CycleEps, MaxCycle, stop ),
      _relax ( Relax ),
      _pConstrVecLev ( NULL ),
      _pProjDirsVecLev ( NULL ),
      _projectThreshold ( 1e-7 ),
      _pPeriodicityHandlers ( this->_depth + 1 ) {

    if ( this->_explicitLevel >= this->_depth )
      throw aol::Exception ( "CFEMultigridEqConstr: explicit level must be < finest level.", __FILE__, __LINE__ );

    if ( this->_explicitLevel <= 0 )
      throw aol::Exception ( "CFEMultigridEqConstr: explicit level must be > 0.", __FILE__, __LINE__ );

    if ( Grid.CT != tpcfe::CFE_CD )
      throw aol::UnimplementedCodeException ( "CFEMultigridEqConstr does not work yet for AssemblyTypes other than CFE_CD", __FILE__, __LINE__ );

    _pConstrVecLev = new aol::RandomAccessContainer<VectorType>[ this->_depth + 1 ];
    _pProjDirsVecLev = new aol::RandomAccessContainer<VectorType>[ this->_depth + 1 ];

    for ( int i = 0; i < this->_depth; ++i )
      _pPeriodicityHandlers[i] = NULL;

    _pPeriodicityHandlers[ this->_depth ] = &PeriodicityHandler;

    init();

    this->setCoarseSolverSteps ( 5000 );
    this->setCoarseSolverEps ( 1.0e-16 );

    initConstraints();
  }

  ~CFEMultigridProjectAvgConstr ( ) {
    delete[] ( _pConstrVecLev );    _pConstrVecLev = NULL;
    delete[] ( _pProjDirsVecLev );  _pProjDirsVecLev = NULL;

    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      if ( _pPeriodicityHandlers[i] ) {
        delete ( _pPeriodicityHandlers[i] );
        _pPeriodicityHandlers[i] = NULL;
      }
    }
  }

public:
  virtual void coarsenGrid ( const int coarseLevel ) {

    GridType *pCoarseGrid = new GridType ( coarseLevel );

    coarsenDomainNodeMask ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // only for CFE_CD, exists here
    /*  coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid );  */ // does NOT exist here
    coarsenDOFMask        ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // exists here.

    pCoarseGrid->setCoarsenedFlag ( true );

    this->setpGrid ( coarseLevel, pCoarseGrid );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( Dest, _pConstrVecLev[ this->_depth ], _pProjDirsVecLev[ this->_depth ], _projectThreshold * Arg.norm() );
    mg::AbstractMultigrid< VectorType >::apply ( Arg, Dest );
  }


  void setProjectThreshold ( const RealType newThres ) {
    _projectThreshold = newThres;
  }


protected:

  virtual void init ( ) {
    cerr << "Setting up list of grids ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenGrid ( i );
      cerr << i << " ";
    }

    cerr << " creating periodicity handlers ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      _pPeriodicityHandlers[i] = new CFEPeriodicityHandler< GridType > ( this->getGridRef ( i ) );
      cerr << i << " ";
    }

    cerr << " setting up restriction and prolongation operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      setProlongationOperator ( i );
      setRestrictionOperator ( i );
      cerr << i << " ";
    }

    cerr << " coarsening operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenOperator ( i );
      cerr << i << " ";
    }

    cerr << " setting up smoothOps ... ";
    for ( int i = ( this->_depth - 0 ); i > this->_explicitLevel; --i ) {
      // need smoothOp on finest level but not on _explicitLevel
      setSmoother ( i );
      cerr << i << " ";
    }

    setCoarsePreconditioner();

    cerr << " done." << endl;

  }

  void initConstraints ( ) {
    // constraints to be set for all levels greater than this->_explicitLevel
    const int DIMENSION_ONE = 1;

    std::vector< MatrixType* > pMassMats ( this->_depth + 1 );
    for ( int i = 0; i < this->_depth + 1; ++i ) {
      pMassMats[i] = NULL;
    }

    // resize, alloc, ...
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      pMassMats[ i ] = new MatrixType ( this->getGridRef ( i ) );

      _pProjDirsVecLev[ i ].reallocate ( DIMENSION_ONE );
      _pProjDirsVecLev[ i ][0].reallocate ( this->getGridRef ( i ) );

      _pConstrVecLev[ i ].reallocate ( DIMENSION_ONE );
      _pConstrVecLev[ i ][0].reallocate ( this->getGridRef ( i ) );
    }

    // projection directions = all-one vectors
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      _pProjDirsVecLev[ i ][0].setAll ( aol::NumberTrait<RealType>::one );
      _pPeriodicityHandlers[ i ]->restrictToPresentDOFs ( _pProjDirsVecLev[ i ][0] );
      _pPeriodicityHandlers[ i ]->restrictPeriodicBC ( _pProjDirsVecLev[ i ][0] );
    }

    // finest mass matrix
    {
      CFEMassOpType finestMassOp ( this->getGridRef ( this->_depth ), aol::ONTHEFLY );
      finestMassOp.assembleAddMatrix ( * ( pMassMats[ this->_depth] ) );
    }
    _pPeriodicityHandlers[ this->_depth ]->periodicallyCollapseMatrix ( * ( pMassMats[ this->_depth] ) );

    // coarsen mass matrix
    for ( int coarseLevel = ( this->_depth - 1 ); coarseLevel >= this->_explicitLevel; --coarseLevel ) {
      switch ( CMode ) {
      case mg::MATRIXMULT_COARSENING:
        {
          const int
            fineLevel   = coarseLevel + 1,
            coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes(),
            fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

          cerr << " left-multiplying by restriction matrix ... ";

          typedef aol::SparseMatrix<RealType> SparseMatrixType;
          SparseMatrixType *re_mat = new SparseMatrixType ( coarseSize, fineSize );
          this->getRestrictionOpRef ( coarseLevel ).assembleMatrix ( ( *re_mat ) );

          SparseMatrixType *re_X_fineop = new SparseMatrixType ( coarseSize, fineSize );

          ( *re_X_fineop ).addMatrixProduct ( ( *re_mat ), * ( pMassMats[ fineLevel ] ) );

          delete ( re_mat );

          cerr << " right-multiplying by prolongation matrix ... ";

          SparseMatrixType *pr_mat  = new SparseMatrixType ( fineSize, coarseSize );
          this->getProlongationOpRef ( coarseLevel ).assembleMatrix ( ( *pr_mat ) );

          ( * ( pMassMats[ coarseLevel ] ) ).addMatrixProduct ( ( *re_X_fineop ), ( *pr_mat ) );

          delete ( re_X_fineop );
          delete ( pr_mat );
        }
        break;

      case mg::ONTHEFLY_COARSENING:
        throw aol::UnimplementedCodeException ( "CFEMultigridProjectAvgConstr: ONTHEFLY_COARSENING not implemented", __FILE__, __LINE__ );
        break;

      default:
        throw aol::Exception ( "tpcfe::CFEMultigridProjectAvgConstr: illegal RePrCoOpType", __FILE__, __LINE__ );
      }
    }

    // constraints = gimel vectors
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      pMassMats[ i ]->apply ( _pProjDirsVecLev[ i ][0], _pConstrVecLev[ i ][0] );
      const RealType volumeFactor = _pProjDirsVecLev[ i ][0] * _pConstrVecLev[ i ][0];
      _pConstrVecLev[ i ][0] /= volumeFactor;


#ifdef DEBUG
      // maybe the threshold need to be chosen bigger?
      smallOrDie ( 1.0 - ( _pProjDirsVecLev[ i ][0] * _pConstrVecLev[ i ][0] ), 1e-15 * _pProjDirsVecLev[ i ][0].size(), "constraint times projection direction == 1?", __FILE__, __LINE__ );
      smallOrDie ( aol::ProjectEqConstrSolver<VectorType>::checkCorrectionResiduumNeutrality ( this->getOperatorRef ( i ), _pProjDirsVecLev[ i ][0] ), 1e-15 * _pProjDirsVecLev[ i ][0].size(), "Correction direction neutral for residuum?", __FILE__, __LINE__ );
#endif
    }

    // delete mass matrices created here
    for ( int i = 0; i <= this->_depth; ++i ) {
      if ( pMassMats[i] ) {
        delete ( pMassMats[i] );
        pMassMats[i] = NULL;
      }
    }

  }


  virtual void setRestrictionOperator ( const int level ) {
    this->setpRestrictionOperator ( level, new tpcfe::CFEPeriodicRestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), * ( _pPeriodicityHandlers[level] ), * ( _pPeriodicityHandlers[level+1] ) ) );
  }

  virtual void setProlongationOperator ( const int level ) {
    this->setpProlongationOperator ( level, new tpcfe::CFEPeriodicProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), * ( _pPeriodicityHandlers[level] ), * ( _pPeriodicityHandlers[level+1] ) ) );
  }


  virtual void coarsenOperator ( const int coarseLevel ) {
    // this will become the operator on coarseLevel
    MatrixType *re_X_op_X_pr = new MatrixType ( this->getGridRef ( coarseLevel ) );

    switch ( CMode ) {
    case mg::MATRIXMULT_COARSENING:
      {
        const int
          fineLevel   = coarseLevel + 1,
          coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes(),
          fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

        cerr << " left-multiplying by restriction matrix ... ";

        typedef aol::SparseMatrix<RealType> SparseMatrixType;
        SparseMatrixType *re_mat = new SparseMatrixType ( coarseSize, fineSize );
        this->getRestrictionOpRef ( coarseLevel ).assembleMatrix ( ( *re_mat ) );

        SparseMatrixType *re_X_fineop = new SparseMatrixType ( coarseSize, fineSize );

        ( *re_X_fineop ).addMatrixProduct ( ( *re_mat ), this->getOperatorRef ( fineLevel ) ); // i. e. get matrix

        delete ( re_mat );

        cerr << " right-multiplying by prolongation matrix ... ";

        SparseMatrixType *pr_mat  = new SparseMatrixType ( fineSize, coarseSize );
        this->getProlongationOpRef ( coarseLevel ).assembleMatrix ( ( *pr_mat ) );

        ( *re_X_op_X_pr ).addMatrixProduct ( ( *re_X_fineop ), ( *pr_mat ) );

        delete ( re_X_fineop );
        delete ( pr_mat );

        for ( int i = 0; i < coarseSize; ++i ) {
          if ( ! ( _pPeriodicityHandlers[ coarseLevel ]->isPresentDOFNode ( i ) ) ) {
            re_X_op_X_pr->set ( i, i, aol::NumberTrait< RealType >::one );
          }
        }

      }
      break;

    case mg::ONTHEFLY_COARSENING:
      throw aol::UnimplementedCodeException ( "CFEMultigridProjectAvgConstr: ONTHEFLY_COARSENING not implemented", __FILE__, __LINE__ );
      break;

    default:
      throw aol::Exception ( "tpcfe::CFEMultigrid: illegal RePrCoOpType", __FILE__, __LINE__ );
    }

    cerr << " repairing zero eigenspaces of operators ... ";
    // modify coarsened operators so that constant vectors are exactly in zero eigenspace
    // note that this does not work if there is only one DOF left on coarseLevel 0: would create singular matrix.
    for ( int i = 0; i < (*re_X_op_X_pr).getNumRows(); ++i ) {
      std::vector< typename aol::Row< RealType >::RowEntry > RowEntries;
      (*re_X_op_X_pr).makeRowEntries ( RowEntries, i );
      RealType rowSum = aol::NumberTrait<RealType>::zero;
      for ( typename std::vector< typename aol::Row<RealType>::RowEntry >::const_iterator it = RowEntries.begin(); it != RowEntries.end(); ++it ) {
        rowSum += it->value;
      }

      if ( _pPeriodicityHandlers[ coarseLevel ]->isPresentDOFNode ( i ) ) {
        (*re_X_op_X_pr).add ( i, i, - rowSum );
      }
    }

    this->setpOperator ( coarseLevel, re_X_op_X_pr );
  }

  virtual void setSmoother ( const int level ) {
    this->setpSmoother ( level, new mg::GaussSeidelSelectiveSmoother< VectorType, MatrixType > ( this->getOperatorRef ( level ), _pPeriodicityHandlers[level]->getPresentDOFNodeMask(), 0, _relax, aol::GAUSS_SEIDEL_FORWARD ) );
  }

  virtual void setCoarsePreconditioner ( ) {
    aol::DiagonalMatrix< RealType >* diagMat = new aol::DiagonalMatrix< RealType > ( this->getGridRef ( this->_explicitLevel ) );
    for ( int i = 0; i < this->getGridRef( this->_explicitLevel ).getNumberOfNodes(); ++i ) {
      const RealType omatii = this->getOperatorRef( this->_explicitLevel ).getDiag ( i );
      if ( omatii != aol::ZOTrait<RealType>::zero ) {
        diagMat->set( i, i, aol::ZOTrait<RealType>::one / aol::Sqr ( omatii ) );
      } else {
        diagMat->set( i, i, aol::ZOTrait<RealType>::one );
      }
    }
    this->_coarsePreconditioner = diagMat;
  }


  // pre-smoothing
  virtual void preSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    this->smooth ( Rhs, X, Level, this->_preSmoothSteps );
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );
  }


  // post-smoothing
  virtual void postSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    // this should not be necessary
    // aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );

    this->smooth ( Rhs, X, Level, this->_postSmoothSteps );

    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );
  }

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    VectorType modifiedRhs ( Rhs ), atrhs ( Rhs, aol::STRUCT_COPY );
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( modifiedRhs, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );

    aol::CompositeOp< VectorType > ata;
    ata.appendReference ( this->getOperatorRef ( Level ) );
    ata.appendReference ( this->getOperatorRef ( Level ) );
    this->getOperatorRef( Level ).apply ( modifiedRhs, atrhs );

    aol::PCGInverse< VectorType > CoarseSolver ( ata, *(this->_coarsePreconditioner), this->coarseSolverEps, this->coarseSolverSteps, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM, cerr );
    CoarseSolver.setQuietMode( this->_verbose < 4 );

    CoarseSolver.apply ( atrhs, X );

    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );
  }

  // end class CFEMultigridProjectAvgConstr
};


/** CFE Multigrid solver with projection to average-0 constraints for periodic boundary conditions in the vector-valued (elasticity) case.
 *  DOF/DomainNode masks are managed by the grid (as usual), periodicity is handled via periodicityHandlers.
 *  \attention only works for periodic geometries
 *  \author Schwen
 */
template< class CFEBlockOpType, class CFEMassOpType, mg::CoarseningMode CMode = mg::MATRIXMULT_COARSENING >
class CFEBlockMultigridProjectAvgConstr : public mg::ExplicitOperatorHierarchyMultigrid < typename CFEBlockOpType::VectorType,
                                                                                          aol::BlockMatrix< typename CFEBlockOpType::ConfiguratorType::MatrixType >,
                                                                                          mg::BlockGaussSeidelSelectiveSmoother < typename CFEBlockOpType::VectorType,
                                                                                          aol::BlockMatrix< typename CFEBlockOpType::ConfiguratorType::MatrixType > > ,
                                                                                          tpcfe::CFEPeriodicRestrictOp< typename CFEBlockOpType::ConfiguratorType::GridType >,
                                                                                          tpcfe::CFEPeriodicProlongOp< typename CFEBlockOpType::ConfiguratorType::GridType >, typename CFEBlockOpType::ConfiguratorType::GridType > {

protected:
  typedef typename CFEBlockOpType::ConfiguratorType::RealType          RealType;
  typedef typename CFEBlockOpType::VectorType                          VectorType;
  typedef typename CFEBlockOpType::ConfiguratorType::MatrixType        MatrixType;
  typedef typename CFEBlockOpType::ConfiguratorType::GridType          GridType;
  typedef aol::BlockMatrix< MatrixType >             BlockMatrixType;
  typedef CFEPeriodicityHandler< GridType >                            PHType;

  RealType                                                             _relax;
  aol::RandomAccessContainer<VectorType>*                              _pConstrVecLev;
  aol::RandomAccessContainer<VectorType>*                              _pProjDirsVecLev;
  RealType                                                             _projectThreshold;
  std::vector< const PHType* >                                         _pPeriodicityHandlers;

public:
  // more options: coarsening stopping criteria!
  CFEBlockMultigridProjectAvgConstr ( const GridType &Grid,
                                      const BlockMatrixType &fineBlockMatrix,
                                      const CFEPeriodicityHandler < GridType > &PeriodicityHandler,
                                      const int ExplicitLevel = 2,
                                      const int PreSmooth = 3, const int PostSmooth = 3,
                                      const RealType Relax = 1.0,
                                      const mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                                      const RealType CycleEps = 1.0e-16,
                                      const int MaxCycle = 100,
                                      const aol::StoppingMode stop = aol::STOPPING_UNSET )
      : mg::ExplicitOperatorHierarchyMultigrid < VectorType,
                                                 BlockMatrixType,
                                                 mg::BlockGaussSeidelSelectiveSmoother< VectorType, BlockMatrixType >,
                                                 tpcfe::CFEPeriodicRestrictOp< GridType >,
                                                 tpcfe::CFEPeriodicProlongOp< GridType >, GridType > ( fineBlockMatrix,
                                                                                                       Grid,
                                                                                                       ExplicitLevel,
                                                                                                       PreSmooth, PostSmooth,
                                                                                                       Relax,
                                                                                                       CycleMode, CycleEps, MaxCycle, stop ),
      _relax ( Relax ),
      _pConstrVecLev ( NULL ),
      _pProjDirsVecLev ( NULL ),
      _projectThreshold ( 1e-7 ),
      _pPeriodicityHandlers ( this->_depth + 1 ) {

    if ( this->_explicitLevel >= this->_depth )
      throw aol::Exception ( "CFEBlockMultigridEqConstr: explicit level must be < finest level.", __FILE__, __LINE__ );

    if ( this->_explicitLevel <= 0 )
      throw aol::Exception ( "CFEMultigridEqConstr: explicit level must be > 0.", __FILE__, __LINE__ );

    if ( Grid.CT != tpcfe::CFE_CD )
      throw aol::UnimplementedCodeException ( "CFEBlockMultigridEqConstr does not work yet for AssemblyTypes other than CFE_CD", __FILE__, __LINE__ );

    _pConstrVecLev = new aol::RandomAccessContainer<VectorType>[ this->_depth + 1 ];
    _pProjDirsVecLev = new aol::RandomAccessContainer<VectorType>[ this->_depth + 1 ];

    for ( int i = 0; i < this->_depth; ++i )
      _pPeriodicityHandlers[i] = NULL;

    _pPeriodicityHandlers[ this->_depth ] = &PeriodicityHandler;

    init();

    this->setCoarseSolverSteps ( 5000 );
    this->setCoarseSolverEps ( 1.0e-16 );

    initConstraints();
  }

  ~CFEBlockMultigridProjectAvgConstr ( ) {
    delete[] ( _pConstrVecLev );    _pConstrVecLev = NULL;
    delete[] ( _pProjDirsVecLev );  _pProjDirsVecLev = NULL;

    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      if ( _pPeriodicityHandlers[i] ) {
        delete ( _pPeriodicityHandlers[i] );
        _pPeriodicityHandlers[i] = NULL;
      }
    }

  }

public:
  virtual void coarsenGrid ( const int coarseLevel ) {

    GridType *pCoarseGrid = new GridType ( coarseLevel );

    coarsenDomainNodeMask ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // only for CFE_CD, exists here
    /*  coarsenDirichletMask  ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid );  */ // does NOT exist here
    coarsenDOFMask        ( this->getGridRef ( coarseLevel + 1 ), *pCoarseGrid ); // exists here.

    pCoarseGrid->setCoarsenedFlag ( true );

    this->setpGrid ( coarseLevel, pCoarseGrid );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( Dest, _pConstrVecLev[ this->_depth ], _pProjDirsVecLev[ this->_depth ], _projectThreshold * Arg.norm() );
    mg::AbstractMultigrid< VectorType >::apply ( Arg, Dest );
  }


  void setProjectThreshold ( const RealType newThres ) {
    _projectThreshold = newThres;
  }


protected:

  virtual void init ( ) {
    cerr << "Setting up list of grids ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenGrid ( i );
      cerr << i << " ";
    }

    cerr << " creating periodicity handlers ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      _pPeriodicityHandlers[i] = new CFEPeriodicityHandler< GridType > ( this->getGridRef ( i ) );
      cerr << i << " ";
    }

    cerr << " setting up restriction and prolongation operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      setProlongationOperator ( i );
      setRestrictionOperator ( i );
      cerr << i << " ";
    }

    cerr << " coarsening operators ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenOperator ( i );
      cerr << i << " ";
    }

    cerr << " setting up smoothOps ... ";
    for ( int i = ( this->_depth - 0 ); i > this->_explicitLevel; --i ) {
      // need smoothOp on finest level but not on _explicitLevel
      setSmoother ( i );
      cerr << i << " ";
    }

    setCoarsePreconditioner();

    cerr << " done." << endl;
  }

  void initConstraints ( ) {
    // constraints to be set for all levels greater than this->_explicitLevel
    const int DIMENSION_THREE = 3;

    std::vector< MatrixType* > pMassMats ( this->_depth + 1 );
    for ( int i = 0; i < this->_depth + 1; ++i ) {
      pMassMats[i] = NULL;
    }

    // resize, alloc, ...
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      pMassMats[ i ] = new MatrixType ( this->getGridRef ( i ) );

      _pProjDirsVecLev[ i ].reallocate ( DIMENSION_THREE );
      for ( int comp = 0; comp < 3; ++comp ) {
        _pProjDirsVecLev[ i ][comp].reallocate ( this->getGridRef ( i ) );
      }

      _pConstrVecLev[ i ].reallocate ( DIMENSION_THREE );
      for ( int comp = 0; comp < 3; ++comp ) {
        _pConstrVecLev[ i ][comp].reallocate ( this->getGridRef ( i ) );
      }
    }

    // projection directions = all-one vectors
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      for ( int comp = 0; comp < 3; ++comp ) {
        _pProjDirsVecLev[ i ][comp][comp].setAll ( aol::NumberTrait<RealType>::one );
        _pPeriodicityHandlers[ i ]->restrictToPresentDOFs ( _pProjDirsVecLev[ i ][comp][comp] ); // the other entries of the multivector are zero anyway.
        _pPeriodicityHandlers[ i ]->restrictPeriodicBC ( _pProjDirsVecLev[ i ][comp][comp] );
      }
    }

    // finest mass matrix
    {
      CFEMassOpType finestMassOp ( this->getGridRef ( this->_depth ), aol::ONTHEFLY );
      finestMassOp.assembleAddMatrix ( * ( pMassMats[ this->_depth] ) );
    }
    _pPeriodicityHandlers[ this->_depth ]->periodicallyCollapseMatrix ( * ( pMassMats[ this->_depth] ) );

    // coarsen mass matrix
    for ( int coarseLevel = ( this->_depth - 1 ); coarseLevel >= this->_explicitLevel; --coarseLevel ) {
      switch ( CMode ) {
      case mg::MATRIXMULT_COARSENING:
        {
          const int
            fineLevel   = coarseLevel + 1,
            coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes(),
            fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

          cerr << " left-multiplying by restriction matrix ... ";

          typedef aol::SparseMatrix<RealType> SparseMatrixType;
          SparseMatrixType *re_mat = new SparseMatrixType ( coarseSize, fineSize );
          this->getRestrictionOpRef ( coarseLevel ).assembleMatrix ( ( *re_mat ) );

          SparseMatrixType *re_X_fineop = new SparseMatrixType ( coarseSize, fineSize );

          ( *re_X_fineop ).addMatrixProduct ( ( *re_mat ), * ( pMassMats[ fineLevel ] ) );

          delete ( re_mat );

          cerr << " right-multiplying by prolongation matrix ... ";

          SparseMatrixType *pr_mat  = new SparseMatrixType ( fineSize, coarseSize );
          this->getProlongationOpRef ( coarseLevel ).assembleMatrix ( ( *pr_mat ) );

          ( * ( pMassMats[ coarseLevel ] ) ).addMatrixProduct ( ( *re_X_fineop ), ( *pr_mat ) );

          delete ( re_X_fineop );
          delete ( pr_mat );
        }
        break;

      case mg::ONTHEFLY_COARSENING:
        throw aol::UnimplementedCodeException ( "CFEBlockMultigridProjectAvgConstr: ONTHEFLY_COARSENING not implemented", __FILE__, __LINE__ );
        break;

      default:
        throw aol::Exception ( "tpcfe::CFEBlockMultigridProjectAvgConstr: illegal RePrCoOpType", __FILE__, __LINE__ );
      }
    }

    // constraints = gimel vectors
    for ( int i = ( this->_depth ); i >= this->_explicitLevel; --i ) {
      for ( int comp = 0; comp < 3; ++comp ) {
        pMassMats[ i ]->apply ( _pProjDirsVecLev[ i ][comp][comp], _pConstrVecLev[ i ][comp][comp] ); // component comp of RandomAccessContainer is MultiVector that only has entries in component Comp
        const RealType volumeFactor = _pProjDirsVecLev[ i ][comp][comp] * _pConstrVecLev[ i ][comp][comp];
        _pConstrVecLev[ i ][comp][comp] /= volumeFactor;

#ifdef DEBUG
        // maybe the threshold need to be chosen bigger?
        smallOrDie ( 1.0 - ( _pProjDirsVecLev[ i ][comp][comp] * _pConstrVecLev[ i ][comp][comp] ), 1e-15 * _pConstrVecLev[ i ][comp][comp].size(), "constraint times projection direction == 1?", __FILE__, __LINE__ );
        smallOrDie ( aol::ProjectEqConstrSolver< VectorType >::checkCorrectionResiduumNeutrality ( this->getOperatorRef ( i ), _pProjDirsVecLev[ i ][comp] ), 1e-15 * _pConstrVecLev[ i ][comp][comp].size() ,
                     "Correction direction neutral for residuum?", __FILE__, __LINE__ );
#endif
      }
    }

    // delete mass matrices created here
    for ( int i = 0; i <= this->_depth; ++i ) {
      if ( pMassMats[i] ) {
        delete ( pMassMats[i] );
        pMassMats[i] = NULL;
      }
    }

  }


  virtual void setRestrictionOperator ( const int level ) {
    this->setpRestrictionOperator ( level, new tpcfe::CFEPeriodicRestrictOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), * ( _pPeriodicityHandlers[level] ), * ( _pPeriodicityHandlers[level+1] ) ) );
  }

  virtual void setProlongationOperator ( const int level ) {
    this->setpProlongationOperator ( level, new tpcfe::CFEPeriodicProlongOp< GridType > ( this->getGridRef ( level ), this->getGridRef ( level + 1 ), * ( _pPeriodicityHandlers[level] ), * ( _pPeriodicityHandlers[level+1] ) ) );
  }


  virtual void coarsenOperator ( const int coarseLevel ) {
    // this will become the operator on coarseLevel
    BlockMatrixType *newpBlockMat = new BlockMatrixType ( this->getGridRef ( coarseLevel ) );

    switch ( CMode ) {
    case mg::MATRIXMULT_COARSENING:
      {
        const int
          fineLevel   = coarseLevel + 1,
          coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes(),
          fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

        cerr << " left-multiplying by restriction matrix ... ";

        typedef aol::BlockMatrix< aol::SparseMatrix<RealType> > BlockSparseMatrixType;
        BlockSparseMatrixType *re_mat = new BlockSparseMatrixType ( 3, 3, coarseSize, fineSize );

        for ( int i = 0; i < 3; ++i )
          this->getRestrictionOpRef ( coarseLevel ).assembleMatrix ( re_mat->getReference ( i, i ) );

        BlockSparseMatrixType *re_X_fineop = new BlockSparseMatrixType ( 3, 3, coarseSize, fineSize );

        ( *re_X_fineop ).addMatrixProduct ( ( *re_mat ), this->getOperatorRef ( fineLevel ) ); // i. e. get matrix

        delete ( re_mat );

        cerr << " right-multiplying by prolongation matrix ... ";

        BlockSparseMatrixType *pr_mat  = new BlockSparseMatrixType ( 3, 3, fineSize, coarseSize );
        for ( int i = 0; i < 3; ++i )
          this->getProlongationOpRef ( coarseLevel ).assembleMatrix ( pr_mat->getReference ( i, i ) );

        ( *newpBlockMat ).addMatrixProduct ( ( *re_X_fineop ), ( *pr_mat ) );

        delete ( re_X_fineop );
        delete ( pr_mat );

        for ( int j = 0; j < 3; ++j ) {
          for ( int i = 0; i < coarseSize; ++i ) {
            if ( ! ( _pPeriodicityHandlers[ coarseLevel ]->isPresentDOFNode ( i ) ) ) {
              newpBlockMat->getReference ( j, j ).set ( i, i, aol::NumberTrait< RealType >::one );
              newpBlockMat->getReference ( j, ( j + 1 ) % 3 ).set ( i, i, aol::NumberTrait< RealType >::zero );
              newpBlockMat->getReference ( j, ( j + 2 ) % 3 ).set ( i, i, aol::NumberTrait< RealType >::zero );
            }
          }
        }
      }
      break;

    case mg::ONTHEFLY_COARSENING:
      throw aol::UnimplementedCodeException ( "CFEBlockMultigridProjectAvgConstr: ONTHEFLY_COARSENING not implemented", __FILE__, __LINE__ );
      break;

    default:
      throw aol::Exception ( "tpcfe::CFEBlockMultigridProjectAvgConstr: illegal RePrCoOpType", __FILE__, __LINE__ );
    }

    cerr << " repairing zero eigenspaces of operators ... ";
    for ( int i = 0; i < (*newpBlockMat).getNumRows(); ++i ) {
      for ( int j = 0; j < (*newpBlockMat).getReference ( i, 0 ).getNumRows(); ++j ) {
        std::vector< typename aol::Row< RealType >::RowEntry > RowEntries;
        (*newpBlockMat).makeUnblockedRowEntries ( RowEntries, i, j );
        RealType rowSum = aol::NumberTrait<RealType>::zero;
        for ( typename std::vector< typename aol::Row<RealType>::RowEntry >::const_iterator it = RowEntries.begin(); it != RowEntries.end(); ++it ) {
          rowSum += it->value;
        }

        if ( _pPeriodicityHandlers[ coarseLevel ]->isPresentDOFNode ( j ) )
          (*newpBlockMat).getReference ( i, i ).add ( j, j, - rowSum );

      }
    }

    this->setpOperator ( coarseLevel, newpBlockMat );
  }

  virtual VectorType* createNewVector ( const int Level ) const {
    // creation of MultiVectors works differently ...
    VectorType* pNewVector = new VectorType ( this->getGridRef ( Level ).getDimOfWorld(), this->getGridRef ( Level ).getNumberOfNodes() );
    return ( pNewVector );
  }

  virtual int getVectorSize ( int const Level ) const {
    // total vector size (for L2-norms) is dimension times vector size
    return ( this->getGridRef ( Level ).getDimOfWorld() *  this->getGridRef ( Level ).getNumberOfNodes() );
  }

  virtual void setSmoother ( const int level ) {
    this->setpSmoother ( level, new mg::BlockGaussSeidelSelectiveSmoother< VectorType, BlockMatrixType > ( this->getOperatorRef ( level ), _pPeriodicityHandlers[level]->getPresentDOFNodeMask(), 0, _relax, aol::GAUSS_SEIDEL_FORWARD ) );
  }

  virtual void setCoarsePreconditioner ( ) {
    aol::BlockMatrix< aol::DiagonalMatrix<RealType> > * diagP = new aol::BlockMatrix< aol::DiagonalMatrix<RealType> > ( this->getGridRef ( this->_explicitLevel ) );

    for ( int i = 0; i < this->getGridRef( this->_explicitLevel ).getNumberOfNodes(); ++i ) {
      aol::Matrix33<RealType> diagBlock;
      for ( short j = 0; j < 3; ++j )
        for ( short k = 0; k < 3; ++k )
          diagBlock.set ( j, k, this->getOperatorRef ( this->_explicitLevel ).getReference ( j, k ).get ( i, i ) );

      aol::Matrix33<RealType> dBInv;
      if ( diagBlock.det() != 0.0 )
        dBInv = diagBlock.inverse();
      else
        dBInv.setIdentity();

      for ( short j = 0; j < 3; ++j )
        for ( short k = 0; k < 3; ++k )
          diagP->getReference ( j, k ).set ( i, i, dBInv.get ( j, k ) );
    }

    this->_coarsePreconditioner = diagP;
  }


  // pre-smoothing
  virtual void preSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    this->smooth ( Rhs, X, Level, this->_preSmoothSteps );
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold );
  }


  // post-smoothing
  virtual void postSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    // this should not be necessary
    // aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold );

    this->smooth ( Rhs, X, Level, this->_postSmoothSteps );

    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold );
  }

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    VectorType modifiedRhs ( Rhs ), atrhs ( Rhs, aol::STRUCT_COPY );
    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( modifiedRhs, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );

    aol::CompositeOp< VectorType > ata;
    ata.appendReference ( this->getOperatorRef ( Level ) );
    ata.appendReference ( this->getOperatorRef ( Level ) );
    this->getOperatorRef( Level ).apply ( modifiedRhs, atrhs );

    aol::PCGInverse< VectorType > CoarseSolver ( ata, *(this->_coarsePreconditioner), this->coarseSolverEps, this->coarseSolverSteps, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM, cerr );
    CoarseSolver.setQuietMode( this->_verbose < 4 );

    CoarseSolver.apply ( atrhs, X );

    aol::ProjectEqConstrSolver<VectorType>::projectToEqConstr ( X, _pConstrVecLev[ Level ], _pProjDirsVecLev[ Level ], _projectThreshold * Rhs.norm() );
  }

  // end class CFEBlockMultigridProjectAvgConstr
};


}


#endif
