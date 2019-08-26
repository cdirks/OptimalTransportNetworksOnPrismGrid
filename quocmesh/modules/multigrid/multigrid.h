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

#ifndef __MULTIGRID_H
#define __MULTIGRID_H

#include <solver.h>

namespace mg {

//! Cycle modes for the Multigrid solver
enum MGCycleMode {
  MULTIGRID_V_CYCLE,
  MULTIGRID_W_CYCLE
};

//! How to compute coarsened operators.
enum CoarseningMode {
  ONTHEFLY_COARSENING = aol::ONTHEFLY, // it makes sense to have the same values here.
  MATRIXMULT_COARSENING = aol::ASSEMBLED
};

/** A general multigrid basis class.
 *  \todo maybe derive from IterativeInverseOp (requires restructuring)
 *  \author Schwen (based on code by others)
 */
template <typename VectorType>
class AbstractMultigrid : public aol::InverseOp< VectorType > {
protected:
  typedef typename VectorType::DataType DataType;

  const int
    _depth,               // finest grid level
    _explicitLevel,       // level up to which we coarsen, on this level, we use an "exact" solver
    _preSmoothSteps,      // number of presmoothing steps
    _postSmoothSteps;     // number of postsmoothing steps

  const MGCycleMode
    _cycleMode;           // cycle mode

  mutable int
    _cycle,
    _verbose;

  mutable DataType
    _cycleEps,           // solver accuracy
    _resNorm,             // norm of residuum
    _resQuot;             // = _resNormNew / _resNormOld

  const int
    _maxNumCycles;      // maximal number of multigrid cycles

  aol::StoppingMode _stopping;

public:

  AbstractMultigrid ( const int Depth,
                      const int ExplicitLevel = 2,
                      const int PreSmooth = 3, const int PostSmooth = 3,
                      const MGCycleMode CycleMode = MULTIGRID_V_CYCLE,
                      const DataType CycleEps = 1.0e-16,
                      const int MaxCycle = 100,
                      const aol::StoppingMode stop = aol::STOPPING_UNSET )
      : _depth ( Depth ), _explicitLevel ( ExplicitLevel ),
      _preSmoothSteps ( PreSmooth ), _postSmoothSteps ( PostSmooth ),
      _cycleMode ( CycleMode ), _cycle ( 0 ),
      _verbose ( 1 ),
      _cycleEps ( CycleEps ),
      _resNorm ( 100.0 ),  _resQuot ( 1.0 ),
      _maxNumCycles ( MaxCycle ), _stopping ( stop ) {

    // constructor does nothing

  }

public:
  virtual ~AbstractMultigrid ( ) {
    // destructor doees nothing
  }

public:
  /** @see setVerboseMode ( int )
   */
  void setQuietMode ( const bool Quiet = true ) {
    _verbose = Quiet ? 0 : 2 ;
  }

  //! greater Value for Verbose causes more output ( 0 - totally quiet )
  void setVerboseMode ( const int Verbose = 2 ) {
    _verbose = Verbose;
  }

  void apply ( const VectorType &Arg, VectorType &Dest ) const;

  void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType tmp ( Dest, aol::DEEP_COPY ); // deep or struct copy? original value as initial guess is no worse than zero, but does this make sense in an applyAdd?
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void setAccuracy ( const DataType cycleEps ) {
    _cycleEps = cycleEps;
  }

  int getCount ( ) const {
    return ( _cycle );
  }

  //! Select the type of stopping criterion
  void setStopping ( const aol::StoppingMode stop ) {
    _stopping = stop;
  }

protected:
  // the following functions must be implemented in derived classes.

  virtual void mgProlongate ( const VectorType &Coarse, VectorType &Fine, const int CoarseLevel ) const = 0 ;

  virtual void mgRestrict ( const VectorType &Fine, VectorType &Coarse, const int CoarseLevel ) const = 0 ;

  virtual void smooth ( const VectorType &RHS, VectorType &Smoothened, const int Level, const int NumIterations ) const = 0 ;

  virtual void applyOperator ( const VectorType &Arg, VectorType &Dest, const int Level ) const = 0 ;

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const = 0 ;

  //! get l2 norm squared of residuum
  virtual DataType getResiduum ( const VectorType &Arg, const VectorType &Dest, const int Level ) const {
    VectorType res ( Arg, aol::STRUCT_COPY );
    applyOperator ( Dest, res, Level );
    res -= Arg;
    return ( res * res );
  }

  //! pre-smoothing
  virtual void preSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    smooth ( Rhs, X, Level, _preSmoothSteps );
  }

  //! coarse-grid correction
  virtual void coarseGridCorrect ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    VectorType
      dummy        ( Rhs, aol::STRUCT_COPY ),
      residuum     ( Rhs, aol::STRUCT_COPY ),
      *prhs_coarse =  createNewVector ( Level - 1 ) , // delete at end of method
      *px_coarse   =  createNewVector ( Level - 1 )  ;

    // Step 2: Compute residuum
    applyOperator ( X, dummy, Level );
    residuum = Rhs;
    residuum -= dummy;

    if ( _verbose >= 3 ) {
      DataType resnormsqr = residuum * residuum;
      cerr << "       Level " << aol::mixedFormat ( Level ) << " down,   Residuum (l2 norm ^2) is " << aol::mixedFormat ( resnormsqr )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resnormsqr / getVectorSize ( Level ) ) ) << " ] " << aol::color::reset << endl;
    }


    // Step 3: Restrict the residuum
    mgRestrict ( residuum, *prhs_coarse, Level - 1 );


    // Step 4: recursive call one level up: compute coarse grid correction
    betweenSmooth ( *prhs_coarse, *px_coarse, Level - 1 );


    // Step 5: Prolongate coarse-grid-correction to fine grid
    mgProlongate ( *px_coarse, dummy, Level - 1 );


    // Step 6: Add correction
    X += dummy;

    delete ( prhs_coarse );
    delete ( px_coarse );

  }

  //! recursive call of the multigrid solver on coarse level
  virtual void betweenSmooth ( const VectorType &coarseRHS, VectorType &coarseX, const int coarseLevel ) const {
    coarseX.setZero( );

    // One recursion for V cycle, two for W cycle
    int numLoops = 0;

    switch ( _cycleMode ) {
    case MULTIGRID_V_CYCLE: numLoops = 1; break;
    case MULTIGRID_W_CYCLE: numLoops = 2; break;
    default: throw aol::Exception ( "unknown MGCycleMode", __FILE__, __LINE__ );
    }

    for ( int i = 0; i < numLoops; ++i ) {
      solve ( coarseRHS, coarseX, coarseLevel );
    }

  }

  //! post-smoothing
  virtual void postSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    smooth ( Rhs, X, Level, _postSmoothSteps );
  }

  //! This is the normalization factor for L2 residuum, i. e. dimension * #DOF in one dimension
  virtual int getVectorSize ( int const Level ) const = 0 ;

  void solve ( const VectorType &Rhs, VectorType &X, const int Level ) const ; // implementation below.

  //! overload this method if creating new vector requires more information, e. g. for MultiVector.
  //! make sure to delete all vectors created by this method!
  virtual VectorType* createNewVector ( const int Level ) const  = 0;

  // end of class AbstractMultigrid
};


template <typename VectorType>
void AbstractMultigrid<VectorType>::solve ( const VectorType &Rhs, VectorType &X, const int Level ) const {
  if ( Level == _explicitLevel ) {

    solveCoarseProblem ( Rhs, X, Level );

  } else {

    // For computing residuum improvement on intermediate levels ( if verboseMode >= 4 )
    DataType
      resNormBeforeCycle  = -1.0,
      resNormAfterCycle   = -1.0,
      resNormBeforeCGC    = -1.0,
      resNormAfterCGC     = -1.0,
      resNormBeforeSmooth = -1.0,
      resNormAfterSmooth  = -1.0;

    if ( _verbose >= 4 ) {
      resNormBeforeCycle = resNormBeforeSmooth = getResiduum ( Rhs, X, Level );
      cerr << "       Level " << aol::mixedFormat ( Level ) << " down, before Smooth   " << aol::mixedFormat ( resNormBeforeSmooth )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resNormBeforeSmooth / getVectorSize ( Level ) ) ) << " ] " << aol::color::reset  << endl;
    }


    // Step 1: Presmoothing
    preSmooth ( Rhs, X, Level );

    if ( _verbose >= 4 ) {
      resNormBeforeCGC = resNormAfterSmooth = getResiduum ( Rhs, X, Level );
      cerr << "       Level " << aol::mixedFormat ( Level ) << " down, after  Smooth   " << aol::mixedFormat ( resNormAfterSmooth )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resNormAfterSmooth / getVectorSize ( Level ) ) ) << " ] " << aol::color::reset
           << "      (~" << aol::mixedFormat ( resNormAfterSmooth / resNormBeforeSmooth ) << " )" << endl;
    }


    // Steps 2 to 6: Coarse grid correction

    coarseGridCorrect ( Rhs, X, Level );


    if ( _verbose >= 4 ) {
      resNormAfterCGC = resNormBeforeSmooth = getResiduum ( Rhs, X, Level );
      cerr << "       Level " << aol::mixedFormat ( Level ) << " up, after  CGC        " << aol::mixedFormat ( resNormAfterCGC )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resNormAfterCGC / getVectorSize ( Level ) ) ) << " ] " << aol::color::reset
           << "      (~" << aol::mixedFormat ( resNormAfterCGC / resNormBeforeCGC ) << " )" << endl;
    }


    // Step 7: Postsmoothing

    postSmooth ( Rhs, X, Level );

    if ( _verbose >= 4 ) {
      resNormAfterCycle = resNormAfterSmooth = getResiduum ( Rhs, X, Level );
      cerr << "       Level " << aol::mixedFormat ( Level ) << " up, after Smooth      "
           << aol::color::cyan
           << aol::mixedFormat ( resNormAfterSmooth )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resNormAfterSmooth / getVectorSize ( Level ) ) ) << " ] " << aol::color::cyan
           << "      (~" << aol::mixedFormat ( resNormAfterSmooth / resNormBeforeSmooth ) << " )"
           << aol::color::reset << endl

           << "       Level " << aol::mixedFormat ( Level ) << " up, Final             "
           << aol::color::red
           << aol::mixedFormat ( resNormAfterCycle )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( resNormAfterCycle / getVectorSize ( Level ) ) ) << " ] " << aol::color::red
           << "      (~" << aol::mixedFormat ( resNormAfterCycle / resNormBeforeCycle ) << " )"
           << aol::color::reset << endl;
    }

  }

  if ( Level == _depth ) {
    _resQuot = _resNorm;
    _resNorm = getResiduum ( Rhs, X, _depth );
    _resQuot = _resNorm / _resQuot ;
  }

}

template <typename VectorType>
void AbstractMultigrid<VectorType>::apply ( const VectorType &Arg, VectorType &Dest ) const {
  // no Dest.clear() to use initial guess.

  _cycle = 0;

  aol::StopWatch timer;

  if ( &Arg == &Dest ) {
    throw aol::Exception ( "Arg and Dest must not be the same object. ", __FILE__, __LINE__ );
  }
  _resNorm = getResiduum ( Arg, Dest, _depth );

  if ( _stopping == aol::STOPPING_UNSET )
    cerr << aol::color::red << "Stopping criterion not set. Will use absolute stopping." << aol::color::reset << endl;

  if ( _verbose >= 1 ) {
    cerr << "---> " << aol::mixedFormat ( 0 ) << "Multigrid Cycle(s), Residuum (l_2 norm ^2) is "
         << aol::color::residuum << aol::mixedFormat ( _resNorm )
         << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( _resNorm / getVectorSize ( _depth ) ) ) << " ] "
         << aol::color::reset << "   \r";
  }

  if ( _verbose >= 2 ) {
    cerr << endl;
  }

  DataType stopEps = _cycleEps;
  if ( _stopping == aol::STOPPING_UNSET )
    cerr << aol::color::red << "Stopping criterion not set. Will use absolute stopping." << aol::color::reset << endl;
  if ( _stopping == aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM )
    stopEps *= _resNorm;
  if ( _stopping == aol::STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE )
    stopEps  *= Arg.normSqr ();
  if ( stopEps == 0 )
    stopEps = _cycleEps;

  do {
    timer.start();
    solve ( Arg, Dest, _depth );
    timer.stop();
    ++_cycle;
    if ( _verbose >= 1 ) {
      cerr << "---> " << aol::mixedFormat ( _cycle ) << "Multigrid Cycle(s), Residuum (l_2 norm ^2) is "
           << aol::color::residuum << aol::mixedFormat ( _resNorm )
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( _resNorm / getVectorSize ( _depth ) ) ) << " ] " << aol::color::reset
           << "    (~ " << aol::mixedFormat ( _resQuot ) << " ) "
           << aol::color::green << "     [ " << aol::mixedFormat ( sqrt ( _resQuot ) ) << " ] " << aol::color::reset
           << " (~ " << aol::shortFormat ( timer.elapsedWallClockTime() ) << " sec WC, " << aol::shortFormat ( timer.elapsedCpuTime() ) << " sec CPU time )         \r";
    }
    if ( _verbose >= 2 ) {
      cerr << endl;
    }

  } while ( ( _resNorm > stopEps ) && ( _cycle < _maxNumCycles ) && ( _resQuot < 1.0 ) );

  if ( _verbose >= 1 ) {

    cerr << endl;

    if ( _verbose >= 1 && _cycle == _maxNumCycles && _resNorm > stopEps ) {
      cerr << aol::color::error << "Too many Multigrid Cycles! " << _maxNumCycles << endl << aol::color::reset ;
    }
    if ( _resQuot >= 1.0 ) {
      cerr << aol::color::error << "Residuum has increased!" << endl << aol::color::reset ;
    }
  }
}



/** A multigrid basis class that stores operators on all grid levels and prolongation/restriction operators.
 *  \author Schwen (based on code by others)
 */

template < typename VectorType, typename OpType, typename SmoothOpType, typename ReOpType, typename PrOpType >
class OperatorHierarchyMultigrid : public AbstractMultigrid<VectorType> {
  // general multigrid solver for a hierarchy of coarse operators

  typedef typename VectorType::DataType DataType;

protected:
  int coarseSolverSteps;
  DataType coarseSolverEps;

public:
  OperatorHierarchyMultigrid ( const int Depth,
                               const int ExplicitLevel = 2,
                               const int PreSmooth = 3, const int PostSmooth = 3,
                               const MGCycleMode CycleMode = MULTIGRID_V_CYCLE,
                               const DataType CycleEps = 1.0e-16,
                               const int MaxCycle = 100,
                               const aol::StoppingMode stop = aol::STOPPING_UNSET  )
      : AbstractMultigrid<VectorType> ( Depth, ExplicitLevel, PreSmooth, PostSmooth, CycleMode, CycleEps, MaxCycle, stop ), coarseSolverSteps ( 2000 ), coarseSolverEps ( 1e-20 ) {

    // nothing to do

  }


public:
  virtual ~OperatorHierarchyMultigrid ( ) {

    // nothing to do

  }

public:

  void setCoarseSolverSteps ( const int cSS ) {
    coarseSolverSteps = cSS;
  }

  void setCoarseSolverEps ( const DataType cSE ) {
    coarseSolverEps = cSE;
  }

protected:

  // the following functions must be implemented in derived classes

  virtual const PrOpType& getProlongationOpRef ( const int CoarseLevel ) const = 0;

  virtual const ReOpType& getRestrictionOpRef ( const int CoarseLevel ) const = 0;

  virtual const OpType& getOperatorRef ( const int Level ) const = 0;

  virtual SmoothOpType& getSmootherRef ( const int Level ) const = 0;

  virtual void mgProlongate ( const VectorType &Coarse, VectorType &Fine, const int CoarseLevel ) const  {
    getProlongationOpRef ( CoarseLevel ).apply ( Coarse, Fine );
  }

  virtual void mgRestrict ( const VectorType &Fine, VectorType &Coarse, const int CoarseLevel ) const  {
    getRestrictionOpRef ( CoarseLevel ).apply ( Fine, Coarse );
  }

  virtual void smooth ( const VectorType &RHS, VectorType &Smoothened, const int Level, const int NumIterations ) const   {
    SmoothOpType& smoother = getSmootherRef ( Level );
    smoother.setNumIterations ( NumIterations );
    smoother.apply ( RHS, Smoothened );
  }

  virtual void applyOperator ( const VectorType &Arg, VectorType &Dest, const int Level ) const {
    getOperatorRef ( Level ).apply ( Arg, Dest );
  }

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    aol::CGInverse<VectorType> inv ( getOperatorRef ( Level ), coarseSolverEps, coarseSolverSteps );
    inv.setQuietMode ( ( this->_verbose < 4 ) );
    inv.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    inv.apply ( Rhs, X );
  }

  // end of class OperatorHierarchyMultigrid
};


/** A multigrid basis class that, given the operator on the finest grid, computes the coarsened operators and explicitely stores them.
 *  \author Schwen (based on code by others)
 */

template < typename VectorType, typename OperatorType, typename SmoothOpType, typename ReOpType, typename PrOpType, typename GridType >
class ExplicitOperatorHierarchyMultigrid : public OperatorHierarchyMultigrid< VectorType, OperatorType, SmoothOpType, ReOpType, PrOpType > {
protected:
  typedef typename VectorType::DataType DataType;

private:
  std::vector< const OperatorType* >                 _operatorList;      // i = operator on level i
  std::vector< SmoothOpType* >                       _smoothOpList;      // i = iterative smoother on level i; cannot be const because we may change number of iteratons
  std::vector< const ReOpType* >                     _restrictOpList;    // i = restriction   i+1 -> i
  std::vector< const PrOpType* >                     _prolongOpList;     // i = prolongation  i -> i+1
  std::vector< const GridType* >                     _gridList;          // i = grid on level i
  DataType                                           _relax;

protected:
  aol::Op<VectorType> *_coarsePreconditioner;

public:
  ExplicitOperatorHierarchyMultigrid ( const OperatorType &fineOperator,
                                       const GridType &fineGrid,
                                       const int ExplicitLevel = 2,
                                       const int PreSmooth = 3, const int PostSmooth = 3,
                                       const DataType Relax = 1.0,
                                       const MGCycleMode CycleMode = MULTIGRID_V_CYCLE,
                                       const DataType CycleEps = 1.0e-16,
                                       const int MaxCycle = 100,
                                       const aol::StoppingMode stop = aol::STOPPING_UNSET  ) : OperatorHierarchyMultigrid < VectorType, OperatorType, SmoothOpType, ReOpType, PrOpType > ( fineGrid.getGridDepth(), ExplicitLevel, PreSmooth, PostSmooth, CycleMode, CycleEps, MaxCycle, stop ),
      _operatorList ( this->_depth + 1 ),
      _smoothOpList ( this->_depth + 1 ),
      _restrictOpList ( this->_depth + 1 ),
      _prolongOpList ( this->_depth + 1 ),
      _gridList ( this->_depth + 1 ),
      _relax ( Relax ),
      _coarsePreconditioner ( NULL ) {

    // initialize vectors with null pointers
    for ( int i = 0; i < this->_depth + 1; ++i ) {
      _operatorList[i]   = NULL;
      _smoothOpList[i]   = NULL;
      _restrictOpList[i] = NULL;
      _prolongOpList[i]  = NULL;
      _gridList[i]       = NULL;
    }

    _operatorList[this->_depth] = &fineOperator;
    _gridList[this->_depth] = &fineGrid;

    // Do not call init() here, call it only in derived classes!
    // Calling init means calling this class's init, regardless of virtual qualifier!

  }

public:
  ~ExplicitOperatorHierarchyMultigrid ( ) {
    if ( this->_verbose >= 1 )
      cerr << "Destroying multigrid solver: freeing memory ... ";

    // delete all operators allocated in this class: all except for finest operator
    for ( unsigned int i = 0; i < _operatorList.size() - 1;  ++i )
      if ( _operatorList[i] ) {
        delete ( _operatorList[i] );
        _operatorList[i] = NULL;
      }

    if ( this->_verbose >= 1 )
      cerr << "o";


    // delete all smooth, prolong and restrict operators allocated in this class.
    for ( unsigned int i = 0; i < _smoothOpList.size();  ++i )
      if ( _smoothOpList[i] ) {
        delete ( _smoothOpList[i] );
        _smoothOpList[i] = NULL;
      }

    if ( this->_verbose >= 1 )
      cerr << "s";


    // delete all restrict ops allocated in this class
    for ( unsigned int i = 0; i < _restrictOpList.size();  ++i )
      if ( _restrictOpList[i] ) {
        delete ( _restrictOpList[i] );
        _restrictOpList[i] = NULL;
      }

    if ( this->_verbose >= 1 )
      cerr << "r";


    // delete all prolong ops allocated in this class
    for ( unsigned int i = 0; i < _prolongOpList.size();  ++i )
      if ( _prolongOpList[i] ) {
        delete ( _prolongOpList[i] );
        _prolongOpList[i] = NULL;
      }

    if ( this->_verbose >= 1 )
      cerr << "p";


    // delete all grids allocated in this class: all except for finest grid
    for ( unsigned int i = 0; i < _gridList.size() - 1;  ++i )
      if ( _gridList[i] ) {
        delete ( _gridList[i] );
        _gridList[i] = NULL;
      }

    if ( this->_verbose >= 1 )
      cerr << "g";


    if ( _coarsePreconditioner ) {
      delete ( _coarsePreconditioner );
      _coarsePreconditioner = NULL;
    }

    if ( this->_verbose >= 1 )
      cerr << " done." << endl;
  }

protected:
  //! Method to set up restriction, prolongation, coarse grid and smoothing operators. Call it in constructors of derived classes.
  virtual void init ( ) {

    cerr << "Setting up list of grids ... ";
    for ( int i = ( this->_depth - 1 ); i >= this->_explicitLevel; --i ) {
      coarsenGrid ( i );
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

protected:
  virtual void coarsenGrid ( const int coarseLevel ) = 0;
  virtual void coarsenOperator ( const int coarseLevel ) = 0;

public:
  virtual const OperatorType& getOperatorRef ( const int Level ) const {
#ifdef BOUNDS_CHECK
    if ( ( Level < this->_explicitLevel ) || ( Level > this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::getOperatorRef: illegal Level", __FILE__, __LINE__ ) );
#endif
    return * ( _operatorList[ Level ] );
  }

  virtual SmoothOpType& getSmootherRef ( const int Level ) const {
#ifdef BOUNDS_CHECK
    if ( ( Level <= this->_explicitLevel ) || ( Level > this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::getSmootherRef: illegal Level", __FILE__, __LINE__ ) );
#endif
    return * ( _smoothOpList[ Level ] );
  }

  virtual const PrOpType& getProlongationOpRef ( const int CoarseLevel ) const {
#ifdef BOUNDS_CHECK
    if ( ( CoarseLevel < this->_explicitLevel ) || ( CoarseLevel >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::getProlongationOpRef: illegal Level", __FILE__, __LINE__ ) );
#endif
    return * ( _prolongOpList[ CoarseLevel ] );
  }

  virtual const ReOpType& getRestrictionOpRef ( const int CoarseLevel ) const {
#ifdef BOUNDS_CHECK
    if ( ( CoarseLevel < this->_explicitLevel ) || ( CoarseLevel >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::getRestrictionOpRef: illegal Level", __FILE__, __LINE__ ) );
#endif
    return * ( _restrictOpList[ CoarseLevel ] );
  }

  const GridType& getGridRef ( const int Level ) const {
#ifdef BOUNDS_CHECK
    if ( ( Level < this->_explicitLevel ) || ( Level > this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::getGridRef: illegal Level", __FILE__, __LINE__ ) );
#endif
    return* ( _gridList[Level] );
  }

  //! Overload for MultiVectors!
  virtual int getVectorSize ( const int Level ) const {
    return ( getGridRef ( Level ).getNumberOfNodes() );
  }

  virtual const OperatorType& getOpReference ( ) {
    return ( getOperatorRef ( this->_depth ) );
  }

protected:
  //! Returns non-constant reference to operator (to avoid accidental use). "Usually", derived classes should not need this method.
  OperatorType& getNonConstOperatorRef ( int Level ) {
    return * ( _operatorList[ Level ] );
  }

  void setpOperator ( const int Level, const OperatorType *pOperator ) {
#ifdef BOUNDS_CHECK
    if ( ( Level < this->_explicitLevel ) || ( Level >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::setpOperator: illegal Level", __FILE__, __LINE__ ) );
#endif
    if ( _operatorList[Level] ) {
      delete ( _operatorList[Level] );
    }
    _operatorList[Level] = pOperator;
  }

  void setpGrid ( const int Level, const GridType *pGrid ) {
#ifdef BOUNDS_CHECK
    if ( ( Level < this->_explicitLevel ) || ( Level >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::setpGrid: illegal Level", __FILE__, __LINE__ ) );
#endif
    if ( _gridList[Level] ){
      delete ( _gridList[Level] );
    }
    _gridList[Level] = pGrid;
  }

  void setpRestrictionOperator ( const int level, const ReOpType *rOp ) {
#ifdef BOUNDS_CHECK
    if ( ( level < this->_explicitLevel ) || ( level >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::setpRestrictionOperator: illegal Level", __FILE__, __LINE__ ) );
#endif

    if ( _restrictOpList[level] ) {
      delete ( _restrictOpList[level] );
    }
    _restrictOpList[level] = rOp;
  }

  void setpProlongationOperator ( const int level, const PrOpType *pOp ) {
#ifdef BOUNDS_CHECK
    if ( ( level < this->_explicitLevel ) || ( level >= this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::setProlongationOperator: illegal Level", __FILE__, __LINE__ ) );
#endif
    if ( _prolongOpList[level] ) {
      delete ( _prolongOpList[level] );
    }
    _prolongOpList[level] = pOp;
  }

  virtual void setRestrictionOperator ( const int level ) = 0;

  virtual void setProlongationOperator ( const int level ) = 0;

  void setpSmoother ( const int level, SmoothOpType *pSmoother ) {
#ifdef BOUNDS_CHECK
    if ( ( level <= this->_explicitLevel ) || ( level > this->_depth ) )
      throw ( aol::Exception ( "mg::ExplicitOperatorHierarchyMultigrid::setpSmoother: illegal Level", __FILE__, __LINE__ ) );
#endif
    if ( _smoothOpList[level] ) {
      delete ( _smoothOpList[level] );
    }
    _smoothOpList[level] = pSmoother;
  }

  virtual void setSmoother ( const int level ) {
    setpSmoother ( level, new SmoothOpType ( *_operatorList[level], 0, _relax ) );
  }

  virtual VectorType* createNewVector ( const int Level ) const {
    return ( new VectorType ( getGridRef ( Level ) ) );
  }

  virtual void setCoarsePreconditioner ( ) {
    // default: do nothing!
  }

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    if ( this->getResiduum ( Rhs, X, Level ) == aol::NumberTrait<DataType>::zero ) // exact float comparison intended
      return;

    if ( _coarsePreconditioner ) {
      aol::PCGInverse<VectorType> inv ( this->getOperatorRef ( Level ), *_coarsePreconditioner, this->coarseSolverEps, this->coarseSolverSteps, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM, cerr );
      inv.setQuietMode ( ( this->_verbose < 4 ) );
      inv.apply ( Rhs, X );
    } else {
      aol::CGInverse<VectorType> inv ( this->getOperatorRef ( Level ), this->coarseSolverEps, this->coarseSolverSteps, aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM, cerr  );
      inv.setQuietMode ( ( this->_verbose < 4 ) );
      inv.apply ( Rhs, X );
    }
  }

  // end of class ExplicitOperatorHierarchyMultigrid
};


}

#endif
