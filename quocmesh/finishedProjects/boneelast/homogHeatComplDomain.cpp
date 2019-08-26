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

template < typename GridType >
void getMacroscopicHeatFlux ( const GridType &grid, const qc::ScalarArray<typename GridType::RealType, qc::QC_3D> &temperature, const typename GridType::RealType aleph, aol::Vec3<typename GridType::RealType> &averageFlux ) {

  typedef typename GridType::RealType RealType;

  RealType vol = aol::NumberTrait<RealType>::zero;

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    el.computeAssembleData ( grid );

    for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); tit.notAtEnd(); ++tit ) {

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
          const tpcfe::CFEVirtualNode<RealType, tpcfe::CFE_CD, RealType> &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
          temperatureAtVertex[vtx] = vn.extrapolate ( temperature );  // extrapolation preserves unit
        }
      }

      aol::Vec3<RealType> differenceTemperature, gradientTemperature;

      for ( int i = 0; i < 3; ++i )
        differenceTemperature[i] = temperatureAtVertex[i] - temperatureAtVertex[3];

      gradientTemperature = TetraMatrix * differenceTemperature;

      const RealType Vol = tetra.getVolume() * aol::Cub ( grid.H() );
      vol += Vol;

      averageFlux += Vol * gradientTemperature;

    } // end tetra in element loop
  } // end element loop

  cerr << "Volume = " << vol << endl;

  averageFlux /= aol::Cub ( aleph ); // volume of bounding box
}



typedef double RealType;

int main ( int, char** ) {
  try {

    typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD >        GridType;
    typedef tpcfe::CFEPeriodicHybridMatrix< GridType >        MatrixType;
    typedef tpcfe::CFEConfigurator < GridType, MatrixType >   ConfiguratorType;
    typedef tpcfe::CFEMassOp< ConfiguratorType >              MassOpType;
    typedef tpcfe::CFEStiffOp< ConfiguratorType >             StiffOpType;        // change to StiffOp with anisotropic coefficient!

    aol::Matrix33<RealType> macroTensor;

    for ( int gradientDirection = 0; gradientDirection < 3; ++gradientDirection ) {

      // set up grid
      GridType grid ( 6 );
      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
      // qc::ShapeLevelsetGenerator<RealType>::generateZ3Levelset ( levelset );
      //  aol::Vec3<RealType> diam ( 0.4/6, 0.3/6, 0.2/6 );
      //  qc::ShapeLevelsetGenerator<RealType>::generateXRotated3DRodsLevelset ( levelset, 3, diam, atan ( 1.0 / 3.0 ) );
      // qc::ShapeLevelsetGenerator<RealType>::generateBallLevelset ( levelset );
      qc::ShapeLevelsetGenerator<RealType>::generateDrunkenZebraLevelset ( levelset );
      // qc::ShapeLevelsetGenerator<RealType>::generateAnisoRandomRodsRemoved3DRodsLevelset ( levelset, 1, 0.4, 0.3, 0.2, 0.0, 0.0, 0.0 );
      levelset.save ( "out/levelset.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();

      grid.setDOFMaskFromDirichletAndDomainNodeMask();
      tpcfe::CFEPeriodicityHandler < GridType > periodicityHandler( grid );

      // set MassOp and StiffOp
      MatrixType massMat ( grid ), stiffMat ( grid );
      {
        MassOpType massOp ( grid, aol::ONTHEFLY );
        massOp.assembleAddMatrix ( massMat );
        StiffOpType stiffOp ( grid, aol::ONTHEFLY );
        stiffOp.assembleAddMatrix ( stiffMat );
      }

      qc::ScalarArray<RealType, qc::QC_3D> rhs ( grid ), soln ( grid ), uSmooth ( grid ), uOsci ( grid );

      // set macroscopic part
      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        uSmooth.set ( *bit, 1.0 * grid.H() * (*bit)[ gradientDirection ] );
      }
      grid.restrictToDomain( uSmooth ); // here, use grid's restriction

      stiffMat.applyAdd ( uSmooth, rhs );                       // note: unperiodized stiffMat.
      periodicityHandler.collapsePeriodicBC ( rhs );

      rhs *= - aol::NumberTrait<RealType>::one;

      cerr << "Periodizing ... ";

      // peroidize massMat and stiffMat.
      periodicityHandler.periodicallyCollapseMatrix( massMat );
      periodicityHandler.periodicallyCollapseMatrix( stiffMat );

      periodicityHandler.restrictNonPresentDOFEntries ( stiffMat, 1.0 );

      cerr << "done." << endl;

      // set neutral functions
      aol::RandomAccessContainer< aol::Vector<RealType> > allOnes ( 1 );
      allOnes[0].reallocate ( grid.getNumberOfNodes() );
      allOnes[0].setAll ( 1.0 );
      periodicityHandler.restrictToPresentDOFs ( allOnes[0] );
      periodicityHandler.restrictPeriodicBC ( allOnes[0] );

      tpcfe::smallOrDie ( aol::ProjectEqConstrSolver<aol::Vector<RealType> >::checkCorrectionResiduumNeutrality ( stiffMat, allOnes[0] ), 1e-10, "Correction direction neutral for residuum?", __FILE__, __LINE__ );


      // set average constraint
      aol::RandomAccessContainer< aol::Vector<RealType> > constrVec ( 1 );
      constrVec[0].reallocate( grid.getNumberOfNodes() );
      massMat.apply ( allOnes[0], constrVec[0] );                 // note: periodized massMat
      const RealType volumeFactor = constrVec[0] * allOnes[0];
      constrVec[0] /= volumeFactor;
      periodicityHandler.restrictToPresentDOFs ( constrVec[0] );  // this should not be necessary
      periodicityHandler.restrictPeriodicBC ( constrVec[0] );     // this should not be necessary


      // PCG solver with projection
      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( stiffMat );
      aol::PCGInverseProjectEqConstr< aol::Vector<RealType> > solver ( stiffMat, prec, constrVec, allOnes, 1e-16, 1000 );

      // solver.setStopping ( aol::STOPPING_ABSOLUTE );
      solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      solver.setProjectThreshold ( 1.e-9 );
      solver.apply ( rhs, soln );


      tpcfe::smallOrDie ( constrVec[0] * soln / soln.size(), 1e-10, "Constraint satisfied?", __FILE__, __LINE__ );

      aol::Vector<RealType> dummy ( rhs, aol::STRUCT_COPY );
      stiffMat.apply ( soln, dummy );
      dummy -= rhs;
      tpcfe::smallOrDie ( dummy.norm() / dummy.size(), 1e-10, "Solves system?", __FILE__, __LINE__ );


      // periodic extension and addition of macroscopic part
      periodicityHandler.extendPeriodicBC ( soln );

      uOsci = soln;

      uOsci.save ( "out/TempOscillatory.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      soln += uSmooth;

      soln.save ( "out/TempSum.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );
      soln.saveSlices ( "out/TempSumCD_%03d.pgm", qc::QC_Z, qc::PGM_UNSIGNED_CHAR_ASCII, NULL, aol::CLIP_THEN_SCALE, soln.getMinValue(), soln.getMaxValue() );

      uSmooth.save ( "out/TempSmooth.dat.bz2", qc::SaveTypeTrait<RealType>::BinarySaveType );

      // get macroscopic heat flux for macroscopic tensor!
      aol::Vec3< RealType > averageFlux;
      getMacroscopicHeatFlux ( grid, soln, 1.0, averageFlux );

      for ( int i = 0; i < 3; ++i )
        macroTensor[ i ][ gradientDirection ] = averageFlux[i];

      cerr << "average flux = " << averageFlux << endl;
    }

    cerr << "Macroscopic Tensor: " << endl << macroTensor << endl;

#if 0
      GridType bigGrid ( 8, qc::QC_3D );
      qc::ScalarArray<RealType, qc::QC_3D> bigLevelset ( bigGrid ), bigSmooth ( bigGrid ), bigOsci ( bigGrid );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        bigLevelset.set ( *bit + qc::CoordType(   0,   0,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType(  64,   0,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType( 128,   0,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType(   0,  64,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType(  64,  64,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType( 128,  64,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType(   0, 128,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType(  64, 128,  0 ), levelset.get(*bit) );
        bigLevelset.set ( *bit + qc::CoordType( 128, 128,  0 ), levelset.get(*bit) );
      }
      bigLevelset.save ( "out/BigZLevelset.dat.bz2", qc::PGM_DOUBLE_BINARY );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        bigOsci.set ( *bit + qc::CoordType(   0,   0,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType(  64,   0,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType( 128,   0,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType(   0,  64,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType(  64,  64,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType( 128,  64,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType(   0, 128,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType(  64, 128,  0 ), uOsci.get(*bit) );
        bigOsci.set ( *bit + qc::CoordType( 128, 128,  0 ), uOsci.get(*bit) );
      }
      bigOsci.save ( "out/BigZOscillatory.dat.bz2", qc::PGM_DOUBLE_BINARY );

      for ( qc::RectangularIterator<qc::QC_3D> bit ( grid ); bit.notAtEnd(); ++bit ) {
        bigSmooth.set ( *bit + qc::CoordType(   0,   0,  0 ), 0 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType(  64,   0,  0 ), 0 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType( 128,   0,  0 ), 0 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType(   0,  64,  0 ), 1 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType(  64,  64,  0 ), 1 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType( 128,  64,  0 ), 1 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType(   0, 128,  0 ), 2 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType(  64, 128,  0 ), 2 + uSmooth.get(*bit) );
        bigSmooth.set ( *bit + qc::CoordType( 128, 128,  0 ), 2 + uSmooth.get(*bit) );
      }
      bigSmooth.save ( "out/BigZSmooth.dat.bz2", qc::PGM_DOUBLE_BINARY );
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
