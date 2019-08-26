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

#ifndef __SMOOTHER_H
#define __SMOOTHER_H

#include <op.h>
#include <vec.h>
#include <GaussSeidel.h>

namespace mg {

/** Basis class for Jacobi and Gauss-Seidel smoother (and maybe others)
 *  In contrast to InverseOps (solvers), Smoothers always perform a given number of iterations and do not care about residuum
 *  \author Schwen (Droske)
 */
template < typename VectorType, typename OpType >
class Smoother : public aol::Op<VectorType> {
  typedef typename VectorType::DataType       DataType;

protected:
  const OpType&  _op;
  int            _numIter; // in contrast to InverseOps, this is a fixed number of iterations and not a maximal one
  DataType       _relax;

public:
  Smoother ( const OpType &Op,
             const int NumIter,
             const DataType Relax )
      : _op ( Op ),
      _numIter ( NumIter ),
      _relax ( Relax ) {
  }

public:

  void setNumIterations ( const int NumIter ) {
    _numIter = NumIter;
  }

  virtual void setRelax ( const DataType Relax ) {
    _relax = Relax;
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType tmp ( Dest, aol::STRUCT_COPY ); // should not need deep copy
    apply ( Arg, tmp );
    Dest += tmp;
  }

  //! implement this in derived classes.
  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const = 0;

private:
  Smoother ( ) {
    throw aol::UnimplementedCodeException ( "mg::Smoother standard constructor not implemented", __FILE__, __LINE__ );
  }

  Smoother ( const Smoother& ){
    throw aol::UnimplementedCodeException ( "mg::Smoother copy constructor not implemented", __FILE__, __LINE__ );
  }

  Smoother& operator= ( const Smoother& ) {
    throw aol::UnimplementedCodeException ( "mg::Smoother:operator= not implemented", __FILE__, __LINE__ );
  }


};

/** Jacobi-Type Smoothing Operator which can be used for multigrid residual-smoothing.
 *  \author Schwen (Droske)
 */
template < typename VectorType, typename OpType >
class JacobiSmoother : public Smoother<VectorType, OpType> {
  typedef typename VectorType::DataType       DataType;

public:
  JacobiSmoother ( const OpType &Op,
                   const int NumIter = 50,
                   const DataType Relax = 1. )
      : Smoother<VectorType, OpType> ( Op, NumIter, Relax ) {
  }

public:
  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType dummy ( Dest, aol::STRUCT_COPY );

    std::vector<typename aol::Row<DataType>::RowEntry > vec;
    for ( int iter = 0; iter < this->_numIter; ++iter ) {
      this->_op.apply ( Dest, dummy );

      for ( int i = 0; i < Dest.size(); ++i ) {
        this->_op.makeRowEntries ( vec, i );

        for ( typename vector<typename aol::Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
          if ( it->col == i ) {
            Dest[i] += this->_relax * ( Arg[i] - dummy[i] ) / it->value;
            break;
          }
        }
      }
    }
  }

};


/** Gauss-Seidel-Type Smoothing Operator which can be used for multigrid residual-smoothing, scalar and vector valued case.
 *  The default behavior is "symmetric" Gauss-Seidel meaning that we alternate between forward and backward loop through rows (for historical reasons).
 *  Red-Black Gauss-Seidel might be useful for parallel computations.
 *  \author Schwen (Droske)
 */
template < typename VectorType, typename OpType >
class GaussSeidelSmoother : public Smoother< VectorType, OpType > {

public:
  typedef typename VectorType::DataType       DataType;

protected:
  aol::GaussSeidelSweepingMode _gss;
  aol::GaussSeidelSweeper<DataType, VectorType, OpType> _sweeper;

public:
  GaussSeidelSmoother ( const OpType &Op,
                        const int Iter = 50,
                        const DataType Relax = 1.0,
                        aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
      : Smoother<VectorType, OpType> ( Op, Iter, Relax ), _gss ( Gss ), _sweeper ( Op, Relax ) {}

public:

  virtual void setRelax ( const DataType Relax ) {
    Smoother<VectorType, OpType>::setRelax ( Relax );
    _sweeper.setRelax ( Relax );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    for ( int iter = 0; iter < this->_numIter; ++iter ) {
      switch ( _gss ) {
      case aol::GAUSS_SEIDEL_FORWARD:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 ? aol::GAUSS_SEIDEL_BACKWARD : aol::GAUSS_SEIDEL_FORWARD ) );
        break;
      case aol::GAUSS_SEIDEL_RED_BLACK:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_EVEN_ONLY );
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_ODD_ONLY );
        break;
      default:
        throw aol::Exception("mg::GaussSeidelSmoother::apply: No Gauss-Seidel sweeping mode selected", __FILE__, __LINE__ );
        break;
      }
    }
  }

}; // end class GaussSeidelSmoother



/** Gauss-Seidel smoother that only operates on subset of vector entries specified by a bitField.
 *  \author Schwen
 */
template < typename VectorType, typename OpType >
class GaussSeidelSelectiveSmoother : public Smoother< VectorType, OpType > {

public:
  typedef typename VectorType::DataType       DataType;

protected:
  aol::GaussSeidelSweepingMode                                    _gss;
  aol::GaussSeidelSelectiveSweeper<DataType, VectorType, OpType>  _sweeper;

public:
  GaussSeidelSelectiveSmoother ( const OpType &Op,
                                 const aol::BitVector &SweepMask,
                                 const int Iter = 50,
                                 const DataType Relax = 1.0,
                                 aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
    : Smoother<VectorType, OpType> ( Op, Iter, Relax ), _gss ( Gss ), _sweeper ( Op, SweepMask, Relax ) {}

  GaussSeidelSelectiveSmoother ( const OpType &Op,
                                 const int Iter = 50,
                                 const DataType Relax = 1.0,
                                 aol::GaussSeidelSweepingMode = aol::GAUSS_SEIDEL_SYMMETRIC )
    : Smoother<VectorType, OpType> ( Op, Iter, Relax ) {
    throw aol::Exception( "mg::GaussSeidelSelectiveSmoother: constructor without sweepMask exists for compatibility only and must not be called", __FILE__, __LINE__ );
  }


public:

  virtual void setRelax ( const DataType Relax ) {
    Smoother<VectorType, OpType>::setRelax ( Relax );
    _sweeper.setRelax ( Relax );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    for ( int iter = 0; iter < this->_numIter; ++iter ) {
      switch ( _gss ) {
      case aol::GAUSS_SEIDEL_FORWARD:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 ? aol::GAUSS_SEIDEL_BACKWARD : aol::GAUSS_SEIDEL_FORWARD ) );
        break;
      case aol::GAUSS_SEIDEL_RED_BLACK:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_EVEN_ONLY );
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_ODD_ONLY );
        break;
      default:
        throw aol::Exception("mg::GaussSeidelSelectiveSmoother::apply: No Gauss-Seidel sweeping mode selected", __FILE__, __LINE__ );
        break;
      }
    }
  }

}; // end class GaussSeidelSelectiveSmoother


/** Gauss-Seidel-Type Smoothing Operator that performs the implicit reordering of unknowns as described for BlockGaussSeidelSmoother.
 *  \author Schwen
 *  \todo   think about how do reduce duplicate code here
 */
template < typename VectorType, typename OpType >
class BlockGaussSeidelSmoother : public Smoother< VectorType, OpType > {

public:
  typedef typename VectorType::DataType       DataType;

protected:
  aol::GaussSeidelSweepingMode _gss;
  aol::BlockGaussSeidelSweeper<DataType, OpType> _sweeper;

public:
  BlockGaussSeidelSmoother ( const OpType &Op,
                             const int Iter = 50,
                             const DataType Relax = 1.0,
                             aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
      : Smoother<VectorType, OpType> ( Op, Iter, Relax ), _gss ( Gss ), _sweeper ( Op, Relax ) {}

public:

  virtual void setRelax ( const DataType Relax ) {
    Smoother<VectorType, OpType>::setRelax ( Relax );
    _sweeper.setRelax ( Relax );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    for ( int iter = 0; iter < this->_numIter; ++iter ) {
      switch ( _gss ) {
      case aol::GAUSS_SEIDEL_FORWARD:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 == 1 ? aol::GAUSS_SEIDEL_BACKWARD : aol::GAUSS_SEIDEL_FORWARD ) );
        break;
      case aol::GAUSS_SEIDEL_RED_BLACK:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_EVEN_ONLY );
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_ODD_ONLY );
        break;
      case aol::GAUSS_SEIDEL_ZEBRA2:
        _sweeper.apply ( Arg, Dest, aol::GAUSS_SEIDEL_ZEBRA2_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_ZEBRA2_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 == 1 ? aol::GAUSS_SEIDEL_ZEBRA2_BACKWARD : aol::GAUSS_SEIDEL_ZEBRA2_FORWARD ) );
        break;
      default:
        throw aol::Exception("mg::GaussSeidelSmoother::apply: No Gauss-Seidel sweeping mode selected", __FILE__, __LINE__ );
        break;
      }
    }
  }

}; // end class BlockGaussSeidelSmoother


/** Block Gauss-Seidel smoother that only operates on subset of vector entries specified by a bitField.
 *  \author Schwen
 */
template < typename VectorType, typename BlockOpType >
class BlockGaussSeidelSelectiveSmoother : public Smoother< VectorType, BlockOpType > {

public:
  typedef typename VectorType::DataType       DataType;

protected:
  aol::GaussSeidelSweepingMode _gss; // currently ignored.
  aol::BlockGaussSeidelSelectiveSweeper<DataType, BlockOpType> _sweeper;

public:
  BlockGaussSeidelSelectiveSmoother ( const BlockOpType &Op,
                                      const int Iter = 50,
                                      const DataType Relax = 1.0,
                                      aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
      : Smoother<VectorType, BlockOpType> ( Op, Iter, Relax ), _gss ( Gss ), _sweeper ( Op, Relax ) {
    throw aol::Exception ( "mg::BlockGaussSeidelSelectiveSmoother: this constructor exists only for compatibility but must not be called. Please specify sweepMask! ", __FILE__, __LINE__ );
  }


  BlockGaussSeidelSelectiveSmoother ( const BlockOpType &Op,
                                      const aol::BitVector &SweepMask,
                                      const int Iter = 50,
                                      const DataType Relax = 1.0,
                                      aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
      : Smoother<VectorType, BlockOpType> ( Op, Iter, Relax ), _gss ( Gss ), _sweeper ( Op, Relax, SweepMask ) {}

public:

  virtual void setRelax ( const DataType Relax ) {
    Smoother<VectorType, BlockOpType>::setRelax ( Relax );
    _sweeper.setRelax ( Relax );
  }


  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    for ( int iter = 0; iter < this->_numIter; ++iter ) {
      switch ( _gss ) {
      case aol::GAUSS_SEIDEL_FORWARD:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 ? aol::GAUSS_SEIDEL_BACKWARD : aol::GAUSS_SEIDEL_FORWARD ) );
        break;
      case aol::GAUSS_SEIDEL_RED_BLACK:
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_EVEN_ONLY );
        _sweeper.apply( Arg, Dest, aol::GAUSS_SEIDEL_ODD_ONLY );
        break;
      case aol::GAUSS_SEIDEL_ZEBRA2:
        _sweeper.apply ( Arg, Dest, aol::GAUSS_SEIDEL_ZEBRA2_FORWARD );
        break;
      case aol::GAUSS_SEIDEL_ZEBRA2_SYMMETRIC:
        _sweeper.apply( Arg, Dest, ( iter % 2 == 1 ? aol::GAUSS_SEIDEL_ZEBRA2_BACKWARD : aol::GAUSS_SEIDEL_ZEBRA2_FORWARD ) );
        break;
      default:
        throw aol::Exception("mg::GaussSeidelSmoother::apply: No Gauss-Seidel sweeping mode selected", __FILE__, __LINE__ );
        break;
      }
    }
  }

}; // end class BlockGaussSeidelSelectiveSmoother


}

#endif
