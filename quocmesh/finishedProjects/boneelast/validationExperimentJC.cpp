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

#include <bitVector.h>
#include <multiArray.h>
#include <parameterParser.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>

#include "computeForce.h"
#include "thresholding.h"

typedef double RealType;
typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_TPOSELAST, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;


int main ( int argc, char **argv ) {
  if ( argc != 2 ) {
    cerr << "usage: validationExperiments PARAMETERFILE" << endl;
    return( EXIT_FAILURE );
  }

  try {
    qc::ScalarArray<RealType, qc::QC_2D>::quietMode = 1;

    aol::ParameterParser parser( argv[1] );

    parser.dump();

    const int
      solverSteps = parser.getInt ( "solverSteps" );

    const RealType
      origRes     = parser.getDouble( "origRes" ),                   // scanning resolution of full-size dataset in m
      shift       = parser.getDouble( "shift" ),                     // in meters or degrees (torsion)
      EAl         = 70.0e9,                                          // in Pa = N / m^2
      nuAl        = 0.35,                                            // in 1
      EPMMA       = 3.0e9,
      nuPMMA      = 0.38;

    char loadFilename[1024], deformationsFNRoot[1024];
    parser.getString ( "loadFilename"  , loadFilename   );
    parser.getString ( "deformationsFNRoot", deformationsFNRoot );

    qc::ScalarArray< RealType, qc::QC_3D > levelsetOrig ( loadFilename );

    const RealType voxelThresholdFinest = RiCaThreshold<RealType> ( levelsetOrig, 0.0, 255.0, aol::Cub ( levelsetOrig.getNumZ() * origRes ) );

    const int
      resFactorFinest = 2,
      resFactorSteps = 1;

    for ( int resF = resFactorFinest + resFactorSteps - 1; resF >= resFactorFinest; --resF ) {
      qc::ScalarArray< RealType, qc::QC_3D > levelsetTmp;
      resampleDatasetByIntegerFactor<RealType> ( levelsetOrig, levelsetTmp, resF  );
      const RealType aleph = levelsetTmp.getNumZ() * origRes * resF;

      qc::ScalarArray< RealType, qc::QC_3D > levelset ( levelsetTmp, aol::STRUCT_COPY );
      thresholdAndDepoempelize< RealType > ( levelsetTmp, voxelThresholdFinest, aol::Cub ( aleph ), levelset, true );

      cerr << "Current resolution: " << levelset.getSize() << ", aleph = " << aleph << endl;

      GridType grid ( qc::GridSize<qc::QC_3D>::createFrom ( levelset ) );

      grid.addStructureFrom ( levelset );

      qc::AArray< GridType::NodalCoeffType, qc::QC_3D > elastCoeff ( grid );
      const GridType::NodalCoeffType ENuMinus ( EAl, nuAl ), ENuPlus ( EPMMA, nuPMMA );
      tpcfe::setCoeffForLevelset ( elastCoeff, levelset, ENuMinus, ENuPlus );

      grid.relaxedDetectAndInitVirtualNodes ( elastCoeff, 1.0e-11, 1.0e-11 );

      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      qc::MultiArray<RealType, 3> dirichletBCs ( grid ), u ( grid ), rhs ( grid );

      setZShift  ( grid, dirichletBCs, DirichletMask, - shift ); // compression

      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << endl;
      cerr << "Total inner volume = " << aol::detailedFormat ( grid.getTotalInnerVolume() * aol::Cub ( aleph ) ) << " m^3." << endl;

      cerr << "Setting up elasticity operator (based on mixed derivative ops)" << endl;
      tpcfe::CFEJCElastOp< ConfiguratorType > elastop ( grid, elastCoeff );

      u -= dirichletBCs;

      elastop.apply ( u, rhs );

      elastop.restrictDirichletEntries();
      grid.restrictDirichletNodes ( rhs );

      elastop.getBlockMatrixRef().getReference(0,0).printStatistics();

      aol::Vector<RealType> solverEpsValues;
      if ( parser.hasVariable ( "solverEpsValues" ) ) {
        parser.getRealVec<RealType> ( "solverEpsValues", solverEpsValues );
      } else {
        solverEpsValues.resize ( 1 );
        solverEpsValues[0] = parser.getDouble( "solverEps" );
      }

      for ( int se = 0; se < solverEpsValues.size(); ++se ) {
        qc::MultiArray<RealType, 3> soln ( grid );
        aol::StopWatch timerSolve;
        timerSolve.start();

        aol::SSORPreconditioner < aol::MultiVector<RealType>, ConfiguratorType::MatrixType > prec ( elastop.getBlockMatrixRef() );
        aol::PCGInverse< aol::MultiVector<RealType> > pcgsolver ( elastop.getBlockMatrixRef(), prec, solverEpsValues[se], solverSteps );
        pcgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE );

        cerr << "memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;
        cerr << "solverEpsilon = " << aol::longScientificFormat ( solverEpsValues[se] ) << endl;

        pcgsolver.apply( rhs, soln );

        timerSolve.stop();
        cerr << "PCG(non-block SSOR) solver took " << timerSolve.elapsedCpuTime()  << endl;

        soln += dirichletBCs;

        if ( solverEpsValues.size() == 1 ) {
          char deformationsFNMask[1024];
          sprintf ( deformationsFNMask, "out/%s_24_resF%02d_def%%d.dat.bz2", deformationsFNRoot, resF );
          cerr << "Starting save to " << deformationsFNMask << endl;
          soln.save ( deformationsFNMask, qc::PGM_DOUBLE_BINARY );
          cerr << "Saved." << endl;
        } // otherwise do not really want it ...

        // force analysis
        const unsigned int surfaceDirection = 2;
        const int surfacePosition = -1;
        const RealType geomTolerance = 1.0e-8f;
        aol::Vec3<RealType> surfaceOuterNormal; surfaceOuterNormal[2] = 1.0f;

        aol::Vector<double> dummy;
        tpcfe::computeForce ( grid, soln, surfaceDirection, surfaceOuterNormal, geomTolerance, surfacePosition, aleph, elastCoeff, cerr, dummy );
      }

    }

  } catch(aol::Exception &ex) {
    ex.dump();
    return( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
