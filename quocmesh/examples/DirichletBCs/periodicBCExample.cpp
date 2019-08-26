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
 *  \brief using periodic boundary conditions
 *
 *  This example shows how periodic boundary conditions can be treated in a FE computation:
 *  \f[     - \Delta u = f    \mbox{ in } \Omega                   \f]
 *  \f[       u(1,y,z) = u(0,y,z)    \mbox{ on } \partial \Omega   \f]
 *  \f[       u(x,1,z) = u(x,0,z)    \mbox{ on } \partial \Omega   \f]
 *  \f[       u(x,y,1) = u(x,y,0)    \mbox{ on } \partial \Omega   \f]
 *  \f[ \int_\Omega u  = 0                                         \f]
 *  where the last constraint makes the solution unique.
 *
 *  \author Schwen
 */


#include <configurators.h>
#include <FEOpInterface.h>
#include <quocMatrices.h>

#include <periodicBC.h>

#include <preconditioner.h>
#include <solver.h>

typedef double RealType;

typedef qc::GridDefinition                                                                                GridType;
typedef aol::SparseMatrix< RealType >                                                                     MatrixType;
typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > ConfiguratorType;
typedef aol::MassOp< ConfiguratorType >                                                                   MassOpType;
typedef aol::StiffOp< ConfiguratorType >                                                                  StiffOpType;

int main ( int, char** ) {
  try {

    GridType grid ( 3, qc::QC_3D );

    qc::QuocPeriodicityHandler < GridType, RealType, qc::QC_3D > periodicityHandler( grid );

    // set MassOp and StiffOp
    MatrixType massMat ( grid ), stiffMat ( grid );
    {
      MassOpType massOp ( grid, aol::ONTHEFLY );
      massOp.assembleAddMatrix ( massMat );
      StiffOpType stiffOp ( grid, aol::ONTHEFLY );
      stiffOp.assembleAddMatrix ( stiffMat );
    }


    qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), source ( grid );


    // set source term (periodic)
    qc::CoordType p1 ( grid.getWidth() / 3, grid.getWidth() / 3, grid.getWidth() / 3 ), p2 = static_cast<short>( 2 ) * p1;
    source.set ( p1, +1.0 );
    source.set ( p2, -1.0 );
    massMat.applyAdd ( source, rhs );
    periodicityHandler.collapsePeriodicBC ( rhs ); // for general rhs, this has an effect

    // peroidize massMat and stiffMat.
    periodicityHandler.periodicallyCollapseMatrix( massMat );
    periodicityHandler.periodicallyCollapseMatrix( stiffMat );

    // set neutral functions
    aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
    allOnes[0].reallocate ( grid.getNumberOfNodes() );
    allOnes[0].setAll ( 1.0 );
    periodicityHandler.restrictPeriodicBC ( allOnes[0] );


    // set average constraint
    aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
    constrVec[0].resize( grid.getNumberOfNodes() );
    massMat.apply ( allOnes[0], constrVec[0] );
    const RealType volumeFactor = constrVec[0] * allOnes[0];
    constrVec[0] /= volumeFactor;
    // massMat is periodized and allOnes restricted, so constrVec[0] only lives on present DOFs.


    // PCG solver with projection
    aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
    aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    solver.setProjectThreshold ( 1.0e-6 );
    solver.apply ( rhs, soln );


    // check whether constraint is satisfied
    if ( constrVec[0] * soln / soln.size() > 1e-10 )
      abort();


    // check whether residuum is really small
    {
      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      stiffMat.apply ( soln, dummy ); // misuse of source as dummy vector
      dummy -= rhs;

      if ( dummy.norm() / dummy.size() > 1e-10 )
        abort();
    }


    // periodic extension: copying values to nodes that are periodic copies of the other nodes
    periodicityHandler.extendPeriodicBC ( soln );

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
