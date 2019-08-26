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
 *  \brief usage of periodic boundary conditions for Composite Finite Elements
 *
 *  Introductory example of usage of CFE for complicated domains and
 *  periodic boundary conditions.
 *
 *  \author Schwen
 */

#include <tpCFEPeriodicBC.h>
#include <tpCFEGrid.h>
#include <tpCFELevelsets.h>
#include <tpCFEStandardOp.h>
#include <linearSmoothOp.h>
#include <shapeLevelsetGenerator.h>


typedef double RealType;

typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >        GridType;
typedef tpcfe::CFEPeriodicHybridMatrix< GridType >        MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
typedef tpcfe::CFEMassOp< ConfiguratorType >              MassOpType;
typedef tpcfe::CFEStiffOp< ConfiguratorType >             StiffOpType;


int main ( int, char** ) {
  try {

    // set up grid
    GridType grid ( 3 );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
    qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    // set up periodicity handler
    tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

    // set up mass and stiffness matrices
    MatrixType massMat ( grid ), stiffMat ( grid );
    {
      MassOpType massOp ( grid, aol::ONTHEFLY );
      massOp.assembleAddMatrix ( massMat );
      StiffOpType stiffOp ( grid, aol::ONTHEFLY );
      stiffOp.assembleAddMatrix ( stiffMat );
    }

    qc::ScalarArray<RealType, qc::QC_3D>  uSmooth ( grid ), rhs ( grid ), soln ( grid );

    // set macroscopic part
    for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
      if ( grid.isDOFNode ( *bit ) ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[2] );
      }
    }

    stiffMat.applyAdd ( uSmooth, rhs );   // note: NONperiodized stiffMat is used here
    periodicityHandler.collapsePeriodicBC ( rhs );

    rhs *= - aol::NumberTrait<RealType>::one;


    // peroidize massMat and stiffMat.
    periodicityHandler.periodicallyCollapseMatrix( massMat );
    periodicityHandler.periodicallyCollapseMatrix( stiffMat );

    periodicityHandler.restrictNonPresentDOFEntries ( stiffMat, 1.0 );


    // set neutral functions
    aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
    allOnes[0].reallocate ( grid.getNumberOfNodes() );
    allOnes[0].setAll ( 1.0 );
    periodicityHandler.restrictToPresentDOFs ( allOnes[0] );

    // check whether adding allOnes[0] to a given vector does not change the residuum, i.e. whether it spans the 0-eigenspace of stiffMat
    tpcfe::smallOrDie ( aol::ProjectEqConstrSolver< aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


    // set average-zero constraint
    aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
    constrVec[0].reallocate( grid.getNumberOfNodes() );
    massMat.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat is used here and allOnes[0] is already restricted to presentDOFs
    const RealType volumeFactor = constrVec[0] * allOnes[0];
    constrVec[0] /= volumeFactor;
    // constrVec[0] only lives on presentDOFs, so no restriction necessary.


    // PCG solver with projection
    aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
    aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    solver.setProjectThreshold ( 1.e-9 );
    solver.apply ( rhs, soln );

    // test whether constraint is satisfied and soln solves system.
    tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

    aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
    stiffMat.apply ( soln, dummy );
    dummy -= rhs;
    tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "Solves system?", __FILE__, __LINE__ );


    // periodic extension and addition of macroscopic part
    periodicityHandler.extendPeriodicBC ( soln );
    soln += uSmooth;

    // soln.saveSlices ( "soln_%03d.pgm", qc::QC_Y, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, soln.getMinValue(), soln.getMaxValue() );
    // Note that values in the top and bottom layer are not constant one or zero at the boundary.
    // This is not a bug - it is the interpolation along edges to virtual nodes on the face of the bounding box that matters.

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
