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

/**
 * \brief Denoising of Time-of-Flight (ToF) data.
 *
 * Options:
 * Data term: L1, L2
 * Regularization: TV, Quadratic Variation
 *
 * \author Bauer
 */

#include <parameterParser.h>
#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <imageTools.h>
#include <triangMesh.h>
#include <signedDistanceOp.h>
#include <signedDistanceSweeping.h>
#include <gradientDescent.h>
#include <gradientflow.h>
#include <hyperelastic.h>
#include <multiArray.h>
#include <generator.h>
#include <RudinOsherFatemi.h>
#include <convolution.h>
#include <meshWithData.h>

#include <iostream>

#include "TOF.h"
#include "SurfaceDenoising.h"


typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_3D, aol::GaussQuadrature<RType, qc::QC_3D, 3> > ConfType3D;

int main ( int argc, char **argv ) {
  try {
    // read parameterfile
    char parameterfilename[1024];
    if ( argc == 2 ) {
      sprintf ( parameterfilename, "%s",  argv[1] );
    } else {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }

    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser ( parameterfilename );

    char templateDataFile[1024];
    char referenceDataFile[1024];
    char saveDirectory[1024];
    int tofSizeX = 0, tofSizeY = 0;
    RType focalLength = 0.;
    RType scaleFactor = 0., shiftOffsetX = 0., shiftOffsetY = 0., shiftOffsetZ = 0.;
    int solverIterations  = 0;
    int smoothingType  = 1;
    int differenceType = 1;
    RType kappa = 0., lambda = 0., mu = 0.;
    RType delta1 = 0., delta2 = 0.;
    int validateDerivative = 0;

    if ( parser.hasVariable ( "templateDataFile" ) ) {
      parser.getString ( "templateDataFile", templateDataFile );
    }
    if ( parser.hasVariable ( "referenceDataFile" ) ) {
      parser.getString ( "referenceDataFile", referenceDataFile );
    }
    if ( parser.hasVariable ( "saveDirectory" ) ) {
      parser.getString ( "saveDirectory", saveDirectory );
    }
    if ( parser.hasVariable ( "tofSizeX" ) ) {
      tofSizeX = parser.getInt ( "tofSizeX" );
    }
    if ( parser.hasVariable ( "tofSizeY" ) ) {
      tofSizeY = parser.getInt ( "tofSizeY" );
    }
    if ( parser.hasVariable ( "focalLength" ) ) {
      focalLength = parser.getReal<RType> ( "focalLength" );
    }
    if ( parser.hasVariable ( "scaleFactor" ) ) {
      scaleFactor = parser.getReal<RType> ( "scaleFactor" );
    }
    if ( parser.hasVariable ( "shiftOffsetX" ) ) {
      shiftOffsetX = parser.getReal<RType> ( "shiftOffsetX" );
    }
    if ( parser.hasVariable ( "shiftOffsetY" ) ) {
      shiftOffsetY = parser.getReal<RType> ( "shiftOffsetY" );
    }
    if ( parser.hasVariable ( "shiftOffsetZ" ) ) {
      shiftOffsetZ = parser.getReal<RType> ( "shiftOffsetZ" );
    }
    if ( parser.hasVariable ( "solverIterations" ) ) {
      solverIterations = parser.getInt ( "solverIterations" );
    }
    if ( parser.hasVariable ( "smoothingType" ) ) {
      smoothingType = parser.getInt ( "smoothingType" );
    }
    if ( parser.hasVariable ( "differenceType" ) ) {
      differenceType = parser.getInt ( "differenceType" );
    }
    if ( parser.hasVariable ( "kappa" ) ) {
      kappa = parser.getReal<RType> ( "kappa" );
    }
    if ( parser.hasVariable ( "lambda" ) ) {
      lambda = parser.getReal<RType> ( "lambda" );
    }
    if ( parser.hasVariable ( "mu" ) ) {
      mu = parser.getReal<RType> ( "mu" );
    }
    if ( parser.hasVariable ( "delta1" ) ) {
      delta1 = parser.getReal<RType> ( "delta1" );
    }
    if ( parser.hasVariable ( "delta2" ) ) {
      delta2 = parser.getReal<RType> ( "delta2" );
    }
    if ( parser.hasVariable ( "validateDerivative" ) ) {
      validateDerivative = parser.getInt ( "validateDerivative" );
    }

    // load template, TOF range data
    std::ifstream InputStreamTOFTemplate;
    InputStreamTOFTemplate.open ( templateDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> templateRangeData ( tofSizeX, tofSizeY );
    templateRangeData.loadRaw ( InputStreamTOFTemplate, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    InputStreamTOFTemplate.close();

    // range to world coord transformation
    RangeToWorldCoordTransformation<RType> transformation ( focalLength, tofSizeX, tofSizeY, scaleFactor, shiftOffsetX, shiftOffsetY, shiftOffsetZ );

    // load reference, CT SDF
    qc::ScalarArray<RType, qc::QC_3D> referenceSDF ( referenceDataFile );

    // generate 2D, 3D grid
    qc::GridSize<qc::QC_2D> grid2DSize ( tofSizeX, tofSizeY );
    const ConfType2D::InitType grid2D ( grid2DSize );
    const ConfType3D::InitType grid3D ( qc::GridSize<qc::QC_3D>::createFrom ( referenceSDF ) );

    // set up energy functionals

    // L1 data term (similarity to noisy ToF data)
    L1DifferenceEnergy<ConfType2D> EL1Difference ( grid2D, templateRangeData, transformation, delta1 );
    VariationOfL1DifferenceEnergy <ConfType2D> DEL1Difference ( grid2D, templateRangeData, transformation, delta1 );

    // L2 data term
    L2DifferenceEnergy<ConfType2D> EL2Difference ( grid2D, templateRangeData, transformation, delta1 );
    VariationOfL2DifferenceEnergy <ConfType2D> DEL2Difference ( grid2D, templateRangeData, transformation, delta1 );

    // TV
    const aol::IsoEnergyOp<ConfType2D> EIso ( grid2D, delta2 );
    const aol::VariationOfIsoEnergyOp<ConfType2D> VarOfEIso ( grid2D, delta2 );

    // Quadratic variation
    // Assemble stiffness matrix L
    qc::FastUniformGridMatrix<RType, ConfType2D::Dim> mat ( qc::GridSize<qc::QC_2D>::createFrom ( grid2D ) );
    aol::StiffOp<ConfType2D> stiff ( grid2D );
    stiff.assembleAddMatrix ( mat );
    DirichletEnergyOp< aol::Vector<RType> > EnergyDenReg ( mat );

    // Matching data term (to reference surface)
    // identity transformation, deformation initialized with 0s
    qc::MultiArray<RType, 2, 3> deformation ( grid2D );
    MatchingEnergyOfRange<ConfType2D, ConfType3D>             ESim (     grid2D, grid3D, referenceSDF, transformation, deformation, 0. );
    FirstPartOfVariationOfMatchingEnergyWRTRange<ConfType2D, ConfType3D>  ESimVariation (  grid2D, grid3D, referenceSDF, transformation, deformation );
    //SecondPartOfVariationOfMatchingEnergyWRTRange<ConfType2D, ConfType3D> ESimVariation1( grid2D, grid3D, deformedReferenceSDF, transformation, deformation, 0. );
    //ThirdPartOfVariationOfMatchingEnergyWRTRange<ConfType2D, ConfType3D>  ESimVariation2( grid2D, grid3D, deformedReferenceSDF, transformation, deformation, 0. );

    // Linear combination of functionals
    aol::LinCombOp<aol::Vector<RType>, aol::Scalar<RType> > E;
    aol::LinCombOp<aol::Vector<RType> > DE;

    std::cout << "Data term: " << differenceType;
    switch ( differenceType ) {
      case 1: // L1
        std::cout << "L1." << std::endl;
        E.appendReference ( EL1Difference, kappa );
        DE.appendReference ( DEL1Difference, kappa );
        break;

      case 2: // L2
        std::cout << "L2." << std::endl;
        E.appendReference ( EL2Difference, kappa );
        DE.appendReference ( DEL2Difference, kappa );
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unsupported differenceType", __FILE__, __LINE__ );
    }

    std::cout << "Regularization: " << smoothingType;
    switch ( smoothingType ) {
      case 1: // Total variation
        std::cout << "TV." << std::endl;
        E.appendReference ( EIso, lambda );
        DE.appendReference ( VarOfEIso, lambda );
        break;

      case 2: // Dirichlet
        std::cout << "Quadratic variation." << std::endl;
        E.appendReference ( EnergyDenReg, lambda );
        DE.appendReference ( stiff, lambda );
        break;
      default:
        throw aol::UnimplementedCodeException ( "Unsupported smoothingType", __FILE__, __LINE__ );
    }

    E.appendReference ( ESim, mu );
    DE.appendReference ( ESimVariation, mu );
    //DE.appendReference( ESimVariation1, mu );
    //DE.appendReference( ESimVariation2, mu );

    // range data to compute
    qc::ScalarArray<RType, qc::QC_2D> newRangeData ( tofSizeX, tofSizeY );
    qc::ScalarArray<RType, qc::QC_2D> tmpRangeData ( tofSizeX, tofSizeY );

    // Validate first derivative
    if ( validateDerivative == 1 ) {
      tmpRangeData = templateRangeData;
      tmpRangeData *= 0.9;

      aol::FirstDerivativeValidator<aol::Vector<RType> > tester ( E, DE, grid2D.H(), aol::FirstDerivativeValidator<aol::Vector<RType> >::LINEAR, 1 );

      char derivativeValidation[1024];
      sprintf ( derivativeValidation, "%sDerivativeValidation.ply", saveDirectory );
      tester.testDirection ( tmpRangeData, derivativeValidation );
      exit ( 0 );
    }

    const bool NoConsoleOutput = false;

    //aol::SimpleGradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, DE, 100);
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, GradientDescentType> GDType;
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, aol::GradientDescent<ConfType, aol::MultiVector<RType> >> GDType;
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::Vector<RType>, aol::H1GradientDescent< ConfType2D, aol::Vector<RType> > > GDType;
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::Vector<RType>, aol::H1GradientDescent< ConfType2D, aol::Vector<RType>, qc::ConvolutionBasedInverseH1Metric<ConfType2D> > > GDType;

    // CG, no smoothing
    typedef aol::GradientDescent<ConfType2D, aol::Vector<RType> > GDType;

    //double tau = 0.01;
    //GDType gradientDescent_solver ( grid, E, DE, 1000 , tau );

    GDType gradientDescent_solver ( grid2D, E, DE, solverIterations );
    //gradientDescent_solver.setFilterWidth( 0.025 );
    gradientDescent_solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG | GDType::DO_NOT_SMOOTH_DESCENT_DIRECTION | ( NoConsoleOutput ? GDType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );

    std::cout << "starting solver..." << std::endl;
    newRangeData = templateRangeData;
    gradientDescent_solver.apply ( newRangeData, tmpRangeData );
    newRangeData = tmpRangeData;
    std::cout << "solver finished." << std::endl;

    //tau = gradientDescent_solver.getStartTau();

    char denDataFileName[1024];
    sprintf ( denDataFileName, "%s/data/range_sm_%i_df_%i_k_%f_l_%f_m_%f_d1_%f_d2_%f_it_%i.raw", saveDirectory, smoothingType, differenceType, kappa, lambda, mu, delta1, delta2, solverIterations );
    newRangeData.save ( denDataFileName, qc::PGM_DOUBLE_BINARY );


    // save surface of raw and denoised data as PLY mesh
    // save distance to reference surface as VTK mesh with scalar field

    aol::TriangMesh<RType> MeshRaw;
    aol::TriangMesh<RType> MeshDen;

    aol::Vector<RType> DenSDF;

    for ( int y = 0; y < tofSizeY; y++ ) {
      for ( int x = 0; x < tofSizeX; x++ ) {
        //RType rNoisy = templateRangeData.get( x, y );
        //aol::Vec3<RType> gNoisy = transformation.evaluateWorldCoord( x, y, rNoisy );

        RType rDen = newRangeData.get ( x, y );
        aol::Vec3<RType> gDen = transformation.evaluateWorldCoord ( x, y, rDen );

        //MeshNoisy.pushBackVertex( gNoisy );
        MeshDen.pushBackVertex ( gDen );

        if ( gDen[0] > 1. || gDen[1] > 1. || gDen[2] > 1. || gDen[0] < 0. || gDen[1] < 0. || gDen[2] < 0. ) {
          DenSDF.pushBack ( 0. );
        } else {
          RType denSDFValue = referenceSDF.interpolate_on01 ( gDen );
          DenSDF.pushBack ( denSDFValue );
        }
      }
    }

    for ( int y = 0; y < tofSizeY - 1; y++ ) {
      for ( int x = 0; x < tofSizeX - 1; x++ ) {
        aol::Vec3<int> UpperFace ( y * tofSizeX + x + 1, y * tofSizeX + x, ( y + 1 ) * ( tofSizeX ) + x );
        aol::Vec3<int> LowerFace ( ( y + 1 ) * ( tofSizeX ) + x, ( y + 1 ) * ( tofSizeX ) + x + 1, y * tofSizeX + x + 1 );

        //MeshNoisy.pushBackTriang( UpperFace );
        //MeshNoisy.pushBackTriang( LowerFace );

        MeshDen.pushBackTriang ( UpperFace );
        MeshDen.pushBackTriang ( LowerFace );
      }
    }

    aol::MeshWithData< aol::TriangMesh<RType> > MeshDenSDF ( MeshDen );
    MeshDenSDF.addData ( DenSDF, "sdf", aol::VERTEX_DATA );

    // save meshes

    //char meshNoisyDataFile[1024];
    //sprintf ( meshNoisyDataFile, "%smeshNoisy.ply", saveDirectory );
    //MeshNoisy.saveAsPLY( meshNoisyDataFile );

    char meshDenDataFile[1024];
    char meshDenSDFDataFile[1024];
    sprintf ( meshDenDataFile, "%s/den/den_sm_%i_df_%i_k_%f_l_%f_m_%f_d1_%f_d2_%f_it_%i.ply", saveDirectory, smoothingType, differenceType, kappa, lambda, mu, delta1, delta2, solverIterations );
    sprintf ( meshDenSDFDataFile, "%s/eval_den/eval_den_sm_%i_df_%i_k_%f_l_%f_m_%f_d1_%f_d2_%f_it_%i.vtk", saveDirectory, smoothingType, differenceType, kappa, lambda, mu, delta1, delta2, solverIterations );
    MeshDen.saveAsPLY ( meshDenDataFile );
    MeshDenSDF.saveAsLegacyVTK ( meshDenSDFDataFile );
  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}

