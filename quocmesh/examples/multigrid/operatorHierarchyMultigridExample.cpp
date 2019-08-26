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

/** \file
 *  \brief using OperatorHierarchyMultigrid
 *
 *   Example for using the OperatorHierarchyMultigrid class, a
 *   multigrid solver for which the coarsened operators need to be
 *   constructed.
 *
 *  \author Schwen
 */

#include <multigrid.h>
#include <smoother.h>
#include <FEOpInterface.h>
#include <configurators.h>
#include <preconditioner.h>
#include <solver.h>
#include <restriction.h>
#include <prolongation.h>
#include <quocMatrices.h>

/** This example illustrates the use of mg::ExplicitOperatorHierarchyMultigrid.
 *  We need to construct the coarse operators ourselves which is shown here.
 *  \author Schwen
 */
typedef aol::Vector<double>                                VectorType;
typedef qc::UniformGridSparseMatrix<double>                OpType;
typedef mg::JacobiSmoother < VectorType, OpType >          SmoothOpType;  // Gauss-Seidel might be more efficient in real applications
typedef qc::RestrictOp< double, qc::STD_MG_RESTRICT >      RestrictOpType;
typedef qc::ProlongOp< double >                            ProlongOpType;
typedef qc::GridDefinition                                 GridType;

class HeatSolver3D : public mg::OperatorHierarchyMultigrid< VectorType, OpType, SmoothOpType, RestrictOpType, ProlongOpType > {
public:
  HeatSolver3D( int level, double tau ) : mg::OperatorHierarchyMultigrid< VectorType, OpType, SmoothOpType, RestrictOpType, ProlongOpType > ( level, 1 ), _level ( level ){ // 1 is the level where we stop coarsening
    _gridList.resize( _level+1 );
    _operatorList.resize( _level+1);
    _smoothOpList.resize( _level+1 );
    _restrictOpList.resize( _level+1 );
    _prolongOpList.resize( _level+1 );

    for ( int l = _level; l >= 1; --l ) {
      _gridList[l] = new qc::GridDefinition ( l, qc::QC_3D );
    }
    // to be safe, one should set the unused entries in the lists to NULL, also after deleting. to keep this example simpler, we omitted that and "know" which entries are set.

    for ( int l = _level; l >= 1; --l ) {
      _operatorList[l] = new OpType( *_gridList[l] );
      aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( *_gridList[l], aol::ONTHEFLY );
      aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( *_gridList[l], aol::ONTHEFLY );
      stiffOp.assembleAddMatrix( *_operatorList[l] );
      (*_operatorList[l]) *= tau;
      massOp.assembleAddMatrix( *_operatorList[l] );
    }

    for ( int l = _level; l >= 1; --l ) {
      _smoothOpList[l] = new SmoothOpType( *_operatorList[l], 3 ); // number of iterations
    }


    for ( int l = _level-1; l >= 1; --l ) {
      _restrictOpList[l] = new RestrictOpType ( *_gridList[l], *_gridList[l+1] );
    }

    for ( int l = _level-1; l >= 1; --l ) {
      _prolongOpList[l] = new ProlongOpType ( *_gridList[l], *_gridList[l+1] );
    }

  }

  ~HeatSolver3D ( ) {
    // delete allocated:
    for ( int l = _level; l >= 1; --l ) {
      delete ( _gridList[l] );
    }

    for ( int l = _level; l >= 1; --l ) {
      delete ( _operatorList[l] );
    }

    for ( int l = _level; l >= 1; --l ) {
      delete ( _smoothOpList[l] );
    }


    for ( int l = _level-1; l >= 1; --l ) {
      delete ( _restrictOpList[l] );
    }

    for ( int l = _level-1; l >= 1; --l ) {
      delete ( _prolongOpList[l] );
    }

  }

protected:
  virtual VectorType* createNewVector ( const int Level ) const {
    return( new VectorType( _gridList[Level]->getNumberOfNodes() ) );
  }

  virtual const ProlongOpType& getProlongationOpRef ( int CoarseLevel ) const {
    return *( _prolongOpList[CoarseLevel] );
  }

  virtual const RestrictOpType& getRestrictionOpRef ( int CoarseLevel ) const {
    return *( _restrictOpList[CoarseLevel] );
  }

  virtual const OpType& getOperatorRef ( int Level ) const {
    return *( _operatorList[Level] );
  }

  virtual SmoothOpType& getSmootherRef ( int Level ) const {
    return *( _smoothOpList[Level] );
  }

  virtual int getVectorSize ( int const Level ) const {
    return ( _gridList[Level]->getNumberOfNodes() );
  }

private:

  std::vector< const GridType* >                     _gridList;          // i = grid on level i
  std::vector< OpType* >                             _operatorList;      // i = operator on level i
  std::vector< SmoothOpType* >                       _smoothOpList;      // i = iterative smoother on level i
  std::vector< RestrictOpType* >                     _restrictOpList;    // i = restriction   i+1 -> i
  std::vector< ProlongOpType* >                      _prolongOpList;     // i = prolongation  i -> i+1
  int                                                _level;

};

int main ( int /*argc*/, char** /*argv*/ ) {
  try {

    // change me!
    const int depth = 3;

    const int steps = 500;
    const double eps = 1.0e-16;

    qc::GridDefinition grid ( depth, qc::QC_3D );

    aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
    aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( grid, aol::ONTHEFLY );

    qc::UniformGridSparseMatrix<double> heatMat ( grid );

    const double tau = 0.5 *  grid.H();

    // set up the matrix we need to invert: M + tau L

    stiffOp.assembleAddMatrix ( heatMat );
    heatMat *= tau;
    massOp.assembleAddMatrix ( heatMat );



    qc::ScalarArray<double, qc::QC_3D>
      img ( grid ),
      rhs ( grid ),
      mg_soln ( grid );


    // we fill our image with some random data to make things interesting ...
    for ( int i = 0; i < img.size(); ++i ) {
      img[i] = static_cast<double> ( rand() ) / static_cast<double> ( RAND_MAX );
    }

    // set up right hand side
    massOp.apply ( img, rhs );

    HeatSolver3D heat_solver( depth, tau );
    heat_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    heat_solver.setVerboseMode ( 2 ); // verbosity mode can be increased up to 6, which produces a lot of output

    // apply solver
    heat_solver.apply ( rhs, mg_soln );


    // We now compare the result to a solution obtained by a standard diagonally preconditioned conjugate gradient solver:

    qc::ScalarArray<double, qc::QC_3D> dpcg_soln ( grid );

    aol::DiagonalPreconditioner< aol::Vector<double> > diag_prec ( heatMat );
    aol::PCGInverse< aol::Vector<double> > dpcg_solver ( heatMat, diag_prec, eps, steps );
    dpcg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    dpcg_solver.apply ( rhs, dpcg_soln );


    // finally, compare the results:
    qc::ScalarArray<double, qc::QC_3D> tmp( grid );

    tmp = mg_soln;
    tmp -= dpcg_soln;

    cerr << "DPCG vs MG:    " << dpcg_soln.norm() << " " << mg_soln.norm() << " " << tmp.norm() <<  endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
