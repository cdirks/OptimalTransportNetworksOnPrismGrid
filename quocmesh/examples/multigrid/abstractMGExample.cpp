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
 *  \brief using AbstractMultigrid
 *
 *   Example for using the AbstractMultigrid class, an abstract
 *   multigrid framework which does not store coarse operators but
 *   implements a method applyOperator that must handle each level
 *   (e. g. on the fly).
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

/** This example illustrates the use of mg::AbstractMultigrid
 *  No coarse grid operators and grid transfer operators are stored, so they need to be constructed each time they are used. This is quite time-consuming.
 *  \author Schwen
 */
typedef aol::Vector<double>                                VectorType;
typedef qc::UniformGridSparseMatrix<double>                OpType;
typedef mg::JacobiSmoother < VectorType, OpType >          SmoothOpType;  // Gauss-Seidel might be more efficient in real applications
typedef qc::RestrictOp< double, qc::STD_MG_RESTRICT >      RestrictOpType;
typedef qc::ProlongOp< double >                            ProlongOpType;
typedef qc::GridDefinition                                 GridType;

class HeatSolver3D : public mg::AbstractMultigrid< VectorType > {
public:
  HeatSolver3D( int level, double tau ) : mg::AbstractMultigrid< VectorType > ( level, 1 ), _tau ( tau ) { // 1 is the level where we stop coarsening
  }

protected:

  virtual void mgProlongate ( const VectorType &Coarse, VectorType &Fine, int CoarseLevel ) const {
    qc::GridDefinition coarse_grid ( CoarseLevel, qc::QC_3D ), fine_grid ( CoarseLevel+1, qc::QC_3D );
    qc::ProlongOp<double> prolong_op ( coarse_grid, fine_grid, aol::ONTHEFLY );
    prolong_op.apply ( Coarse, Fine );
  }

  virtual void mgRestrict ( const VectorType &Fine, VectorType &Coarse, int CoarseLevel ) const {
    qc::GridDefinition coarse_grid ( CoarseLevel, qc::QC_3D ), fine_grid ( CoarseLevel+1, qc::QC_3D );
    RestrictOpType restrict_op ( coarse_grid, fine_grid, aol::ONTHEFLY );
    restrict_op.apply ( Fine, Coarse );
  }

  virtual void smooth ( const VectorType &RHS, VectorType &Smoothened, int Level, int NumIterations ) const {
    qc::GridDefinition grid ( Level, qc::QC_3D );
    aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
    aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( grid, aol::ONTHEFLY );

    OpType heatMat( grid );
    stiffOp.assembleAddMatrix( heatMat );
    heatMat *= _tau;
    massOp.assembleAddMatrix( heatMat );

    SmoothOpType smoothOp ( heatMat, 1 );

    for ( int i = 0; i < NumIterations; ++i ) {
      smoothOp.apply( RHS, Smoothened );
    }
  }

  virtual void applyOperator ( const VectorType &Arg, VectorType &Dest, int Level ) const {
    qc::GridDefinition grid ( Level, qc::QC_3D );
    aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
    aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( grid, aol::ONTHEFLY );

    aol::LinCombOp< VectorType > heatOp;
    heatOp.appendReference( massOp );
    heatOp.appendReference( stiffOp, _tau );

    heatOp.apply ( Arg, Dest );
  }

  virtual void solveCoarseProblem ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    qc::GridDefinition grid ( Level, qc::QC_3D );
    aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > > stiffOp ( grid, aol::ONTHEFLY );
    aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > >  massOp ( grid, aol::ONTHEFLY );

    aol::LinCombOp< VectorType > heatOp;
    heatOp.appendReference( massOp );
    heatOp.appendReference( stiffOp, _tau );

    aol::CGInverse<VectorType> coarseSolver ( heatOp );
    coarseSolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    coarseSolver.setQuietMode();
    coarseSolver.apply( Rhs, X );
  }

  virtual VectorType* createNewVector ( const int level ) const {
    return( new VectorType( getVectorSize(level) ) );
  }

  virtual int getVectorSize ( int const level ) const {
    return ( aol::Cub( (1 << level) + 1 ) ); // ( 2^level + 1 )^3
  }

private:
  double _tau;

};


int main ( int /*argc*/, char** /*argv*/ ) {
  try {

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
