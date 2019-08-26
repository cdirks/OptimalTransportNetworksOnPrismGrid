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

#include <tpCFEPeriodicBC.h>

#include <configurators.h>
#include <FEOpInterface.h>
#include <quocMatrices.h>
#include <anisoStiffOps.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEStandardOp.h>
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>


typedef double RealType;


int main ( int, char** ) {
  try {

    cerr << "This program does not do what the name suggests" << endl;
    abort();

      typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >        GridType;
      typedef tpcfe::CFEPeriodicHybridMatrix< GridType >        MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
      typedef tpcfe::CFEMassOp< ConfiguratorType >              MassOpType;
      typedef tpcfe::CFEStiffOp< ConfiguratorType >             StiffOpType;

      // set up grid
      GridType grid ( 7 );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<RealType>::generateZLevelset ( levelset );


      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();


      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), uSmooth ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[1] );
      }
      grid.restrictToDomain( uSmooth ); // here, use grid's restriction

      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      for ( qc::RectangularBoundaryIterator<qc::QC_3D> bbit ( grid ); bbit.notAtEnd(); ++bbit ) {
        DirichletMask.set ( *bbit, true );
      }

      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      // set MassOp and StiffOp
      MatrixType massMat ( grid ), stiffMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        StiffOpType stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix ( stiffMat );
      }

      soln -= uSmooth;

      tpcfe::restrictNonDomainEntries ( grid, stiffMat );
      grid.restrictToDomain ( soln );

      stiffMat.apply ( soln, rhs );

      tpcfe::restrictDirichletEntries ( grid, stiffMat );
      grid.restrictDirichletNodes ( rhs );

      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
      aol::PCGInverse< aol::Vector<RealType> > solver ( stiffMat, prec, 1e-16, 1000 );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.apply ( rhs, soln );

      soln.save ( "out/ZOscillatory.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      soln += uSmooth;

      soln.save ( "out/ZSum.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      uSmooth.save ( "out/ZSmooth.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
