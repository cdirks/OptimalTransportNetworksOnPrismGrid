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
 *  \brief usage of Composite Finite Element code
 *
 *  Introductory example of usage of CFE in the cases of complicated
 *  domains (with Dirichlet boundary values) or discontinuous
 *  coefficients (time steppinng).
 *
 *  \author Schwen
 */

#include <tpCFELevelsets.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>

#include <solver.h>
#include <preconditioner.h>
#include <shapeLevelsetGenerator.h>

typedef double RealType;

int main ( int, char** ) {
  try {

    { // complicated domains, elliptic problem

      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >       GridType;
      typedef tpcfe::CFEBandMatrix < GridType >                MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;
      typedef tpcfe::CFEStiffOp  < ConfiguratorType >          StiffOpType;

      GridType grid ( 3 );

      qc::ScalarArray< RealType, qc::QC_3D > levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();


      // set Dirichlet mask and appropriate boundary values
      qc::BitArray< qc::QC_3D > DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) ); // initialized with false
      for ( qc::RectangularBoundaryIterator< qc::QC_3D > bit ( grid ); bit.notAtEnd(); ++bit ) {
        DirichletMask.set ( *bit, true );
      }

      // grid must know about Dirichlet nodes because it keeps track of DOFs and non-DOFs for complicated domain
      // This is particularly important when using CFEHybridMatrices, these must be assembled after setting the masks correctly.
      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();


      // set boundary conditions
      qc::ScalarArray< RealType, qc::QC_3D > bcond ( grid );
      for ( qc::RectangularBoundaryIterator< qc::QC_3D > bit ( grid ); bit.notAtEnd(); ++bit ) {
        bcond.set ( *bit, static_cast< RealType >( (*bit)[2] ) );
      }

      // set up stiffOp
      MatrixType stiffMat ( grid );
      {
        StiffOpType stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix ( stiffMat );
        // instead, we could use aol::ASSEMBLED operator, call stiffOp.assembleMatrix(); and use stiffOp.getMatrixRef() instead of using separate stiffMat
      }

      grid.restrictToDomain ( bcond ); // vector should not contain entries corresponding to non-domain nodes

      qc::ScalarArray< RealType, qc::QC_3D > Lbcond ( grid );
      stiffMat.apply ( bcond, Lbcond );


      qc::ScalarArray< RealType, qc::QC_3D > rhs ( grid ), soln ( grid );
      rhs -= Lbcond;


      // transform to zero Dirichlet boundary conditions; matrix should have 1.0 on these diagonal entries
      tpcfe::restrictDirichletEntries ( grid, stiffMat, 1.0 );
      grid.restrictDirichletNodes ( rhs );

      // solve
      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat ); // small virtual tetrahedra can lead to big condition number; diagonal preconditioning remedies this
      aol::PCGInverse< aol::Vector<RealType> > solver ( stiffMat, prec, 1.0e-16, 500 );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.apply ( rhs, soln );

      // do not forget to transform back to nonzero Dirichlet BCs
      soln += bcond;

    }


    { // discontinuous coefficients, time-dependent problem

      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_TPOS >       GridType;
      typedef tpcfe::CFEHybridMatrix < GridType >              MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;
      typedef tpcfe::CFEMassOp < ConfiguratorType >            MassOpType;
      typedef tpcfe::CFEStiffOpWI < ConfiguratorType >         WStiffOpType;

      GridType grid ( 3 );

      qc::ScalarArray< RealType, qc::QC_3D > levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
      grid.addStructureFrom ( levelset );

      // attention: if temperature is considered as the unknown
      // quantity, mass-specific heat capacity and heat conductivity
      // need to be considered separately and also a weighted MassOp
      // must be used

      qc::AArray<RealType, qc::QC_3D> coefficients ( grid );
      tpcfe::setCoeffForLevelset ( coefficients, levelset, 2.0 /* interior */, 1.0 /* exterior */ );
      grid.detectAndInitVirtualNodes ( coefficients );

      MassOpType massOp ( grid, aol::ASSEMBLED );
      MatrixType sysMat ( grid );
      {
        WStiffOpType stiffOp ( coefficients, grid, aol::ONTHEFLY );
        RealType tau = 0.5 * grid.H();
        stiffOp.assembleAddMatrix ( sysMat );
        sysMat *= tau;
        massOp.assembleAddMatrix ( sysMat );
      }

      qc::ScalarArray<RealType, qc::QC_3D> uOld ( grid ), uNew ( grid ), rhs ( grid );
      aol::NoiseOperator<RealType> noiseOp;
      noiseOp.applySingle ( uOld );

      massOp.apply ( uOld, rhs );

      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( sysMat );
      aol::PCGInverse< aol::Vector<RealType> > solver ( sysMat, prec, 1.0e-16, 500 );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

      solver.apply ( uOld, uNew );

    }

  } catch ( aol::Exception &exc ) {

    exc.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );

}
