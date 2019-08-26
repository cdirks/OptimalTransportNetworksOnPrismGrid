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
typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_CD, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;

const bool DO_FORCE_ANALYSIS = true;

// #define LOAD_DONT_COMPUTE 1

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
      solverSteps = parser.getInt ( "solverSteps" ),
      shiftType   = parser.getInt ( "shiftType" ); // -2 = compression, 2 = tension, 10 = torsion

    const RealType
      origRes     = parser.getDouble( "origRes" ),                   // scanning resolution of full-size dataset in m
      shift       = parser.getDouble( "shift" ),                     // in meters or degrees (torsion)
      E           = parser.getDouble( "EModulus" ),                  // in Pa = N / m^2
      nu          = parser.getDouble( "nu" ),                        // in 1
      lambda      = tpcfe::computeLambda ( E, nu ),
      mu          = tpcfe::computeMu ( E, nu );

    char loadFilename[1024], deformationsFNRoot[1024];
    parser.getString ( "loadFilename"  , loadFilename   );
    parser.getString ( "deformationsFNRoot", deformationsFNRoot );

    cerr << "lambda = " << lambda << ", mu = " << mu << endl;

    qc::ScalarArray< RealType, qc::QC_3D > levelsetOrig ( loadFilename );

    // const RealType voxelThresholdFinest = RiCaThreshold<RealType> ( levelsetOrig, parser.getDouble ( "HistoAnalyzeMin" ), parser.getDouble ( "HistoAnalyzeMax" ), aol::Cub ( levelsetOrig.getNumZ() * origRes ) );
    const RealType voxelThresholdLoaded = parser.getDouble ( "DeviceThreshold" );

    const int
      resFactorFinest = parser.getInt ( "resFactorFinest" ),
      resFactorSteps = parser.getInt ( "resFactorSteps" );

    for ( int resF = resFactorFinest + resFactorSteps - 1; resF >= resFactorFinest; --resF ) {
      qc::ScalarArray< RealType, qc::QC_3D > levelsetTmp;
      resampleDatasetByIntegerFactor<RealType> ( levelsetOrig, levelsetTmp, resF  );
      const RealType aleph = levelsetTmp.getNumZ() * origRes * resF;

      qc::ScalarArray< RealType, qc::QC_3D > levelset ( levelsetTmp, aol::STRUCT_COPY );
      thresholdAndDepoempelize< RealType > ( levelsetTmp, voxelThresholdLoaded, aol::Cub ( aleph ), levelset, true );

      cerr << "Current resolution: " << levelset.getSize() << ", aleph = " << aleph << endl;

#if 0
      const short slicNo = 60 / resF;
      char fnroot[1024];
      sprintf ( fnroot, "out/visSlice%02d", resF );
      saveSlice<RealType> ( levelset, slicNo, fnroot );

      char fname[1024];
      sprintf ( fname, "out/visVolume%02d", resF );
      levelset.save ( fname, qc::PGM_DOUBLE_BINARY );
#endif

      GridType grid ( qc::GridSize<qc::QC_3D>::createFrom ( levelset ) );

      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();

      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      DirichletMask.setAll ( false );

      qc::MultiArray<RealType, 3> dirichletBCs ( grid ), u ( grid ), rhs ( grid );

      switch ( shiftType ) {
      case -2:
        setZShift  ( grid, dirichletBCs, DirichletMask, - shift ); // compression
        break;
      case 2:
        setZShift  ( grid, dirichletBCs, DirichletMask,   shift ); // tension
        break;
      case 10:
        setTorsion ( grid, dirichletBCs, DirichletMask,   shift ); // torsion
        dirichletBCs *= aleph;
        break;
      default:
        throw aol::Exception ( "Illegal deformation mode", __FILE__, __LINE__ );
      }

      grid.setDirichletMask ( DirichletMask );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

#ifndef LOAD_DONT_COMPUTE
      cerr << "Using 3 x " << grid.getNumDOFs() << " = " << 3 * grid.getNumDOFs() << " degrees of freedom, total number of nodes in cube: " << grid.getNumberOfNodes() << endl;
      cerr << "Total inner volume = " << aol::detailedFormat ( grid.getTotalInnerVolume() * aol::Cub ( aleph ) ) << " m^3." << endl;

      cerr << "Setting up elasticity operator (based on mixed derivative ops)" << endl;
      tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );

      u -= dirichletBCs;

      elastop.restrictNonDomainEntries();

      grid.restrictToDomain ( u );

      elastop.apply ( u, rhs );

      elastop.restrictDirichletEntries();
      grid.restrictDirichletNodes ( rhs );

      elastop.getBlockMatrixRef().getReference(0,0).printStatistics();
#endif
      aol::Vector<RealType> solverEpsValues;
      if ( parser.hasVariable ( "solverEpsValues" ) ) {
        parser.getRealVec<RealType> ( "solverEpsValues", solverEpsValues );
      } else {
        solverEpsValues.resize ( 1 );
        solverEpsValues[0] = parser.getDouble( "solverEps" );
      }

      for ( int se = 0; se < solverEpsValues.size(); ++se ) {
        qc::MultiArray<RealType, 3> soln ( grid );
#ifndef LOAD_DONT_COMPUTE
        aol::StopWatch timerSolve;
        timerSolve.start();

        aol::SSORPreconditioner < aol::MultiVector<RealType>, ConfiguratorType::MatrixType > prec ( elastop.getBlockMatrixRef() );
        aol::PCGInverse< aol::MultiVector<RealType> > pcgsolver ( elastop.getBlockMatrixRef(), prec, solverEpsValues[se], solverSteps );
        pcgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE );

        cerr << "memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;
        cerr << "solverEpsilon = " << aol::longScientificFormat ( solverEpsValues[se] ) << endl;

        pcgsolver.apply( rhs, soln );

        timerSolve.stop();
        cerr << "PCG(non-block SSOR) solver took " << timerSolve.elapsedWallClockTime() << " s (wall clock time)" << endl;

        soln += dirichletBCs;

        if ( solverEpsValues.size() == 1 ) {
          char deformationsFNMask[1024];
          sprintf ( deformationsFNMask, "out/%s_24_resF%02d_def%%d.dat.bz2", deformationsFNRoot, resF );
          cerr << "Starting save to " << deformationsFNMask << endl;
          soln.save ( deformationsFNMask, qc::PGM_DOUBLE_BINARY );
          cerr << "Saved." << endl;
        } // otherwise do not really want it ...
#else
        char deformationsFNMask[1024];
        sprintf ( deformationsFNMask, "out/%s_24_resF%02d_def%%d.dat.bz2", deformationsFNRoot, resF );
        cerr << "Loading from " << deformationsFNMask << endl;
        soln.load ( deformationsFNMask );
#endif
        if ( DO_FORCE_ANALYSIS ) {
          // force analysis
          const unsigned int surfaceDirection = 2;
          const int surfacePosition = -1;
          const RealType geomTolerance = 1.0e-8;
          aol::Vec3<RealType> surfaceOuterNormal; surfaceOuterNormal[2] = 1.0;

          tpcfe::IsotropicElasticityCoefficient<RealType> ENuMinus ( E, nu ), ENuPlus ( aol::NumberTrait<RealType>::NaN, aol::NumberTrait<RealType>::NaN );
          qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > isoCoeff ( grid );
          tpcfe::setCoeffForLevelset ( isoCoeff, levelset, ENuMinus, ENuPlus );

          aol::Vector<double> dummy;
          tpcfe::computeForce ( grid, soln, surfaceDirection, surfaceOuterNormal, geomTolerance, surfacePosition, aleph, isoCoeff, cerr, dummy );
        }
      }

    }

  } catch(aol::Exception &ex) {
    ex.dump();
    return( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
