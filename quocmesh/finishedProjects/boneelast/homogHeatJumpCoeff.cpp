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


// todo: remove duplicate code!

template < typename GridType >
void getMacroscopicHeatFlux ( const GridType &grid, const qc::AArray<typename GridType::RealType, qc::QC_3D> &coeff,
                              const qc::ScalarArray<typename GridType::RealType, qc::QC_3D> &temperature, const typename GridType::RealType aleph,
                              aol::Vec3<typename GridType::RealType> &averageFlux ) {

  typedef typename GridType::RealType RealType;

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    el.computeAssembleData ( grid );

    tpcfe::CFEWeightProvider< RealType, RealType > weightProvider ( coeff, el );

    for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {

      const tpcfe::CFETetra<RealType> &tetra = *tit;
      aol::Vec3<RealType> VertexCoords[4];

      for ( int vtx = 0; vtx < 4; vtx++ )
        tetra.computeGlobalCoordinate ( VertexCoords[ vtx ], el, vtx );

      aol::Vec3<RealType> directions[3];
      for ( int i = 0; i < 3; ++i ) {
        directions[i] = aleph * grid.H() * ( VertexCoords[ i ] - VertexCoords[ 3 ] );
      }

      aol::Matrix33<RealType> directionsMatrix;
      for ( int i = 0; i < 3; ++i )
        directionsMatrix.setRow ( i, directions[i] );

      aol::Matrix33<RealType> TetraMatrix = directionsMatrix.inverse();  // unit is m^{-1}

      aol::Vec<4,RealType> temperatureAtVertex;

      for ( int vtx = 0; vtx < 4; ++vtx ) {

        const int
          lIdx0 = tetra( vtx, 0 ),
          lIdx1 = tetra( vtx, 1 ),
          gIdx0 = el.globalIndex(lIdx0);

        if ( lIdx1 == 11 ) { // at regular nodes, we have data
          temperatureAtVertex[vtx] = temperature[gIdx0];
        } else {             // at virtual nodes, we need to interpolate
          const int gIdx1 = el.globalIndex(lIdx1);
          const tpcfe::CFEVirtualNode< RealType, GridType::CT, RealType > &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
          temperatureAtVertex[vtx] = vn.extrapolate ( temperature );  // extrapolation preserves unit
        }
      }

      aol::Vec3<RealType> differenceTemperature, gradientTemperature;

      for ( int i = 0; i < 3; ++i )
        differenceTemperature[i] = temperatureAtVertex[i] - temperatureAtVertex[3];

      gradientTemperature = TetraMatrix * differenceTemperature;

      const RealType Vol = tetra.getVolume() * aol::Cub ( grid.H() );

      averageFlux += Vol * weightProvider.meanWeight( tetra.getSign() ) * gradientTemperature;

    } // end tetra in element loop
  } // end element loop

  averageFlux /= aol::Cub ( aleph ); // volume of bounding box
}



typedef double RealType;

int main ( int, char** ) {
  try {

    const int
      level = 8,
      nRods = 2;

    // const RealType xTh = 0.38, yTh = 1./3., zTh = 0.24; // A, ganz neu -> b
    // const RealType xPc = 0.0, yPc = 0.0, zPc = 0.0;

    // const RealType xTh = 1./3., yTh = 1./3., zTh = 1./3.; // B, neu -> a
    // const RealType xPc = 0.0, yPc = 0.0, zPc = 0.0;

    // const RealType xTh = 1./3., yTh = 1./3., zTh = 1./3.; // C, neu -> c
    // const RealType xPc = 0.1, yPc = 0.1, zPc = 0.1;

    // const RealType xTh = 1./3., yTh = 1./3., zTh = 1./3.; // D, neu -> d
    // const RealType xPc = 0.3, yPc = 0.0, zPc = 0.0;

    const RealType xTh = 0.38, yTh = 0.33, zTh = 0.24; // E
    const RealType xPc = 0.0, yPc = 0.0, zPc = 0.0;

    cerr << nRods << " rods with  d/l ratios " << xTh << ", " << yTh << ", " << zTh
         << ", removal percentages " << xPc << " " << yPc << " " << zPc << ", on level " << level << endl;

    const tpcfe::ConstraintType CT = tpcfe::CFE_TPOS;
    typedef tpcfe::CFEGrid < RealType, CT >                   GridType;
    typedef tpcfe::CFEPeriodicHybridMatrix< GridType >        MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
    typedef tpcfe::CFEMassOp< ConfiguratorType >              MassOpType;
    typedef tpcfe::CFEStiffOpWI< ConfiguratorType >           WStiffOpType;

    aol::Matrix33<RealType> macroTensor;

    for ( int gradientDirection = 0; gradientDirection < 3; ++gradientDirection ) {

      // set up grid
      GridType grid ( level );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

      qc::ShapeLevelsetGenerator<RealType>::generatePeriodicAnisoRandom3DRodsLevelset ( levelset, nRods, xTh /nRods, yTh / nRods, zTh / nRods, xPc, yPc, zPc );

      grid.addStructureFrom ( levelset );
      qc::AArray<RealType, qc::QC_3D> coeff ( grid );
      tpcfe::setCoeffForLevelset ( coeff, levelset, 237.0, 0.19); // Al-PMMA
      grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

      // set MassOp and StiffOp
      MatrixType massMat ( grid ), wStiffMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        WStiffOpType wStiffOp ( coeff, grid, aol::ONTHEFLY );
        wStiffOp.assembleAddMatrix ( wStiffMat );
      }

      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), uSmooth ( grid ), uOsci ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[ gradientDirection ] );
      }

      wStiffMat.applyAdd ( uSmooth, rhs );                       // note: unperiodized wStiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      cerr << "Periodizing ... ";

      // peroidize massMat and wStiffMat.
      periodicityHandler.periodicallyCollapseMatrix( massMat );
      periodicityHandler.periodicallyCollapseMatrix( wStiffMat );

      periodicityHandler.restrictNonPresentDOFEntries ( wStiffMat, 1.0 );

      cerr << "done." << endl;

      // set neutral functions
      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( allOnes[0] );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver<aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( wStiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


      // set average constraint
      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].reallocate( grid.getNumberOfNodes() );
      massMat.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;


      // PCG solver with projection
      aol::SSORPreconditioner< aol::Vector<RealType>, MatrixType > prec ( wStiffMat );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( wStiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.e-9 );
      solver.apply ( rhs, soln );


      tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-9, "Constraint satisfied?", __FILE__, __LINE__ );

      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      wStiffMat.apply ( soln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-9, "Solves system?", __FILE__, __LINE__ );


      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( soln );

      uOsci = soln;

      //       uOsci.save ( "out/TempOscillatory.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      soln += uSmooth;

      //       soln.save ( "out/TempSum.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );
      //       soln.saveSlices ( "out/TempSumJC_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, soln.getMinValue(), soln.getMaxValue() );

      //       uSmooth.save ( "out/TempSmooth.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      // get macroscopic heat flux for macroscopic tensor!
      aol::Vec3< RealType > averageFlux;
      getMacroscopicHeatFlux<GridType> ( grid, coeff, soln, 1.0, averageFlux );

      for ( int i = 0; i < 3; ++i )
        macroTensor[ i ][ gradientDirection ] = averageFlux[i];

      cerr << "average flux = "
           << aol::longScientificFormat ( averageFlux[0] ) << " "
           << aol::longScientificFormat ( averageFlux[1] ) << " "
           << aol::longScientificFormat ( averageFlux[2] ) << endl;
    }

    cerr << "Macroscopic Tensor: " << endl
         << aol::longScientificFormat ( macroTensor[0][0] ) << " " << aol::longScientificFormat ( macroTensor[0][1] ) << " " << aol::longScientificFormat ( macroTensor[0][2] ) << endl
         << aol::longScientificFormat ( macroTensor[1][0] ) << " " << aol::longScientificFormat ( macroTensor[1][1] ) << " " << aol::longScientificFormat ( macroTensor[1][2] ) << endl
         << aol::longScientificFormat ( macroTensor[2][0] ) << " " << aol::longScientificFormat ( macroTensor[2][1] ) << " " << aol::longScientificFormat ( macroTensor[2][2] ) << endl;

    cerr << "Symmetrized macroscopic tensor: " << endl
         << aol::longScientificFormat ( macroTensor[0][0] ) << " " << aol::longScientificFormat ( 0.5 * ( macroTensor[0][1] + macroTensor[1][0] ) ) << " " << aol::longScientificFormat ( 0.5 * ( macroTensor[0][2] + macroTensor[2][0] ) ) << endl
         << aol::longScientificFormat ( 0.5 * ( macroTensor[0][1] + macroTensor[1][0] ) ) << " " << aol::longScientificFormat ( macroTensor[1][1] ) << " " << aol::longScientificFormat ( 0.5 * ( macroTensor[1][2] + macroTensor[2][1] ) ) << endl
         << aol::longScientificFormat ( 0.5 * ( macroTensor[0][2] + macroTensor[2][0] ) ) << " " << aol::longScientificFormat ( 0.5 * ( macroTensor[1][2] + macroTensor[2][1] ) ) << " " << aol::longScientificFormat ( macroTensor[2][2] ) << endl;

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
