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
 * \brief Joint Surface Registration and Denoising
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
#include <meshWithData.h>

#include <iostream>

#include "TOF.h"
#include "SurfaceDenoising.h"
#include "SurfaceJointRegistrationDenoising.h"


typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType2D;
typedef qc::RectangularGridConfigurator<RType, qc::QC_3D, aol::GaussQuadrature<RType, qc::QC_3D, 3> > ConfType3D;

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
    char templateNoNoiseDataFile[1024];
    char referenceDataFile[1024];
    char deformedReferenceDataFile[1024];
    char saveDirectory[1024];
    int tofSizeX = 0, tofSizeY = 0;
    RType focalLength = 0.;
    RType scaleFactor = 0., shiftOffsetX = 0., shiftOffsetY = 0., shiftOffsetZ = 0.;
    int solverIterations = 0;
    int smoothingType = 1;
    int differenceType = 1;
    RType alpha = 0.;
    RType kappa = 0., lambda = 0., mu = 0.;
    RType delta1 = 0., delta2 = 0.;
    int validateDerivative = 0;
    int gradientDescent = 0;

    if ( parser.hasVariable ( "templateDataFile" ) ) {
      parser.getString ( "templateDataFile", templateDataFile );
    }
    if ( parser.hasVariable ( "templateNoNoiseDataFile" ) ) {
      parser.getString ( "templateNoNoiseDataFile", templateNoNoiseDataFile );
    }
    if ( parser.hasVariable ( "referenceDataFile" ) ) {
      parser.getString ( "referenceDataFile", referenceDataFile );
    }
    if ( parser.hasVariable ( "deformedReferenceDataFile" ) ) {
      parser.getString ( "deformedReferenceDataFile", deformedReferenceDataFile );
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
    if ( parser.hasVariable ( "alpha" ) ) {
      alpha = parser.getReal<RType> ( "alpha" );
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
    if ( parser.hasVariable ( "gradientDescent" ) ) {
      gradientDescent = parser.getInt ( "gradientDescent" );
    }

    // load template, TOF range data
    std::ifstream InputStreamTOFTemplate;
    InputStreamTOFTemplate.open ( templateDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> templateRangeData ( tofSizeX, tofSizeY );
    templateRangeData.loadRaw ( InputStreamTOFTemplate, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    InputStreamTOFTemplate.close();

    // load template (no noise), TOF range data
    std::ifstream InputStreamTOFTemplateNoNoise;
    InputStreamTOFTemplateNoNoise.open ( templateNoNoiseDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> templateNoNoiseRangeData ( tofSizeX, tofSizeY );
    templateNoNoiseRangeData.loadRaw ( InputStreamTOFTemplateNoNoise, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    InputStreamTOFTemplateNoNoise.close();

    // range to world coord transformation
    RangeToWorldCoordTransformation<RType> transformation ( focalLength, tofSizeX, tofSizeY, scaleFactor, shiftOffsetX, shiftOffsetY, shiftOffsetZ );

    // load reference, CT SDF
    qc::ScalarArray<RType, qc::QC_3D> referenceSDF ( referenceDataFile );

    // generate 2D, 3D grid
    qc::GridSize<qc::QC_2D> grid2DSize ( tofSizeX, tofSizeY );
    const ConfType2D::InitType grid2D ( grid2DSize );
    const ConfType3D::InitType grid3D ( qc::GridSize<qc::QC_3D>::createFrom ( referenceSDF ) );

    // deform reference (SDF)
    qc::DataGenerator<ConfType3D> deformer ( grid3D );
    aol::MultiVector<RType> syntheticDeformation ( ConfType3D::Dim, grid3D.getNumberOfNodes() );
    RType deformationStrength = 0.1;
    deformer.generateNonLinearDeformation ( deformationStrength, syntheticDeformation );
    //qc::ScalarArray<RType, qc::QC_3D> deformedReferenceSDF( referenceDataFile );
    //qc::DeformImage<ConfType3D> ( referenceSDF, grid3D, deformedReferenceSDF, syntheticDeformation );

    qc::ScalarArray<RType, qc::QC_3D> deformedReferenceSDF ( deformedReferenceDataFile );

    // save deformed reference SDF
    char deformedReferenceSDFDataFile[1024];
    sprintf ( deformedReferenceSDFDataFile, "%sdeformedReferenceSDF.dat", saveDirectory );
    deformedReferenceSDF.save ( deformedReferenceSDFDataFile, qc::PGM_DOUBLE_BINARY );

    // set up energy functionals

    // Assemble stiffness matrix L
    //aol::StiffOp<ConfType2D> stiff( grid2D, aol::ASSEMBLED );

    CompleteRegistrationDenoisingEnergy<ConfType2D, ConfType3D> E (
      grid2D,
      grid3D,
      deformedReferenceSDF,
      templateRangeData,
      transformation,
      smoothingType,
      differenceType,
      alpha,
      kappa,
      lambda,
      mu,
      delta1,
      delta2 );

    VariationOfCompleteRegistrationDenoisingEnergy<ConfType2D, ConfType3D> DE (
      grid2D,
      grid3D,
      deformedReferenceSDF,
      templateRangeData,
      transformation,
      smoothingType,
      differenceType,
      alpha,
      kappa,
      lambda,
      mu,
      delta1,
      delta2 );

    //// range data to compute
    qc::ScalarArray<RType, qc::QC_2D> tmpRangeData ( tofSizeX, tofSizeY );
    qc::ScalarArray<RType, qc::QC_2D> newRangeData ( tofSizeX, tofSizeY );

    // Deformation to compute
    // <RType, 2, 3> integrieren auf 2D, 3 Unbekannte (Deformation)
    qc::MultiArray<RType, 2, 3> mtmp ( grid2D );
    qc::MultiArray<RType, 2, 3> phi ( grid2D );

    // initialize range estimate
    newRangeData = templateRangeData;

    aol::MultiVector<RType> deformation_r_tmp ( 0, 0 );
    deformation_r_tmp.appendReference ( mtmp[0] );
    deformation_r_tmp.appendReference ( mtmp[1] );
    deformation_r_tmp.appendReference ( mtmp[2] );
    deformation_r_tmp.appendReference ( tmpRangeData );

    aol::MultiVector<RType> deformation_r ( 0, 0 );
    deformation_r.appendReference ( phi[0] );
    deformation_r.appendReference ( phi[1] );
    deformation_r.appendReference ( phi[2] );
    deformation_r.appendReference ( newRangeData );


    // Validate first derivative
    if ( validateDerivative == 1 ) {
      qc::DataGenerator<ConfType2D> generator ( grid2D );
      generator.generateNonLinearDeformation ( 0.05, deformation_r_tmp );

      deformation_r_tmp[0] = templateRangeData;
      deformation_r_tmp[1] = templateRangeData;
      deformation_r_tmp[0] *= 0.88;
      deformation_r_tmp[1] *= 0.98;

      deformation_r_tmp[2].setAll ( 0.05 );
      deformation_r_tmp[3] = templateRangeData;
      deformation_r_tmp[3] *= 0.9;

      aol::FirstDerivativeValidator<aol::MultiVector<RType> > tester ( E, DE, grid2D.H(), aol::FirstDerivativeValidator<aol::MultiVector<RType> >::LINEAR, 1 );

      char derivativeValidation[1024];
      sprintf ( derivativeValidation, "%sDerivativeValidation.ply", saveDirectory );
      tester.testDirection ( deformation_r_tmp, derivativeValidation );
      exit ( 0 );
    }

    //double tau = 0.01;
    const bool NoConsoleOutput = false;

    //aol::SimpleGradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, DE, 100);
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, GradientDescentType> GDType;
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, aol::GradientDescent<ConfType, aol::MultiVector<RType> >> GDType;
    typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::MultiVector<RType>, aol::H1GradientDescent< ConfType2D, aol::MultiVector<RType> > > GDType;
    //typedef aol::GradientDescent<ConfType2D, aol::MultiVector<RType> > GDType;

    typedef aol::GradientDescentComponentWiseTimestepControlled<ConfType2D, GDType> GDComponentWiseType;
    //typedef aol::GradientDescentComponentWiseTimestepControlled<ConfType2D, aol::H1GradientDescent< ConfType2D, aol::MultiVector<RType> > > GDComponentWiseType;


    //GDType gradientDescent_solver ( grid, E, DE, 1000 , tau );
    //GDType gradientDescent_solver ( grid2D, E, DE, solverIterations );

    aol::Vector<int> tauBlocks ( 2 );
    tauBlocks[0] = 3;
    tauBlocks[1] = 1;

    GDComponentWiseType gradientDescent_solver ( grid2D, E, DE, tauBlocks, solverIterations );
    //gradientDescent_solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG | GDType::DO_NOT_SMOOTH_DESCENT_DIRECTION | ( NoConsoleOutput ? GDType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );
    //gradientDescent_solver.setConfigurationFlags ( GDComponentWiseType::USE_NONLINEAR_CG | GDComponentWiseType::DO_NOT_SMOOTH_DESCENT_DIRECTION | ( NoConsoleOutput ? GDComponentWiseType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );
    gradientDescent_solver.setConfigurationFlags ( GDComponentWiseType::USE_NONLINEAR_CG | ( NoConsoleOutput ? GDComponentWiseType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );

    //gradientDescent_solver.setSaveDirectory(saveDirectory);
    //gradientDescent_solver.activateSaving();
    //gradientDescent_solver.setSaveTimestepOffset(5);

    std::cout << "starting solver..." << std::endl;
    gradientDescent_solver.apply ( deformation_r, deformation_r_tmp );
    std::cout << "solver finished." << std::endl;

    //std::cout << "Energy._cntIn = " << Energy._cntIn << std::endl;
    //std::cout << "Energy._cntOut = " << Energy._cntOut << std::endl;

    deformation_r = deformation_r_tmp;

    //tau = gradientDescent_solver.getStartTau();

    //char filenameResult[1024];
    //sprintf ( filenameResult, "%s/resultRangeData.png", saveDirectory );
    //newRangeData.setOverflowHandling ( aol::CLIP_THEN_SCALE, newRangeData.getMinValue(), newRangeData.getMaxValue() );
    //newRangeData.savePNG ( filenameResult );


    // apply the estimated deformation on the TOF 3D data and save result as PLY mesh
    // in addition, for convenience, save original TOF mesh
    aol::TriangMesh<RType> MeshRaw;
    aol::TriangMesh<RType> MeshDen;
    aol::TriangMesh<RType> MeshDef;
    aol::TriangMesh<RType> MeshDefDen;

    aol::TriangMesh<RType> MeshRawNoNoiseDef;
    aol::TriangMesh<RType> MeshRawNoNoiseGTDef;

    aol::Vector<RType> DefDenSDF;
    aol::Vector<RType> DenSDF;
    aol::Vector<RType> DefSDF;
    aol::Vector<RType> DefL1Diff;

    qc::ScalarArray<RType, qc::QC_3D> synthDefC0 ( syntheticDeformation[0], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );
    qc::ScalarArray<RType, qc::QC_3D> synthDefC1 ( syntheticDeformation[1], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );
    qc::ScalarArray<RType, qc::QC_3D> synthDefC2 ( syntheticDeformation[2], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );

    //synthDefC0.save( "../../../Data/synthDefC0.raw" );
    //synthDefC1.save( "../../../Data/synthDefC1.raw" );
    //synthDefC2.save( "../../../Data/synthDefC2.raw" );

    //qc::ScalarArray<RType, qc::QC_2D> synthDef2DC0( grid2D );
    //qc::ScalarArray<RType, qc::QC_2D> synthDef2DC1( grid2D );
    //qc::ScalarArray<RType, qc::QC_2D> synthDef2DC2( grid2D );

    //qc::ScalarArray<RType, qc::QC_2D> estDef2DC0( grid2D );
    //qc::ScalarArray<RType, qc::QC_2D> estDef2DC1( grid2D );
    //qc::ScalarArray<RType, qc::QC_2D> estDef2DC2( grid2D );

    qc::ScalarArray<RType, qc::QC_2D> phi0 ( phi[0], grid2D.getNumX(), grid2D.getNumY() );
    qc::ScalarArray<RType, qc::QC_2D> phi1 ( phi[1], grid2D.getNumX(), grid2D.getNumY() );
    qc::ScalarArray<RType, qc::QC_2D> phi2 ( phi[2], grid2D.getNumX(), grid2D.getNumY() );



    for ( int y = 0; y < tofSizeY; y++ ) {
      for ( int x = 0; x < tofSizeX; x++ ) {
        RType rRaw = templateRangeData.get ( x, y );
        aol::Vec3<RType> gRaw = transformation.evaluateWorldCoord ( x, y, rRaw );

        RType rRawNoNoise = templateNoNoiseRangeData.get ( x, y );
        aol::Vec3<RType> gRawNoNoise = transformation.evaluateWorldCoord ( x, y, rRawNoNoise );

        RType rDen = newRangeData.get ( x, y );
        aol::Vec3<RType> gDen = transformation.evaluateWorldCoord ( x, y, rDen );

        aol::Vec3<RType> gDef;
        gDef[0] = gRaw[0] + phi[0].get ( x, y );
        gDef[1] = gRaw[1] + phi[1].get ( x, y );
        gDef[2] = gRaw[2] + phi[2].get ( x, y );

        aol::Vec3<RType> gDefDen;
        gDefDen[0] = gDen[0] + phi[0].get ( x, y );
        gDefDen[1] = gDen[1] + phi[1].get ( x, y );
        gDefDen[2] = gDen[2] + phi[2].get ( x, y );

        RType defDenSDFValue = deformedReferenceSDF.interpolate_on01 ( gDefDen );
        DefDenSDF.pushBack ( defDenSDFValue );

        RType denSDFValue = referenceSDF.interpolate_on01 ( gDen );
        DenSDF.pushBack ( denSDFValue );

        aol::Vec3<RType> gRawNoNoiseDef;
        gRawNoNoiseDef[0] = gRawNoNoise[0] + phi[0].get ( x, y );
        gRawNoNoiseDef[1] = gRawNoNoise[1] + phi[1].get ( x, y );
        gRawNoNoiseDef[2] = gRawNoNoise[2] + phi[2].get ( x, y );

        RType defSDFValue = deformedReferenceSDF.interpolate_on01 ( gRawNoNoiseDef );
        DefSDF.pushBack ( defSDFValue );

        aol::Vec3<RType> gRawNoNoiseGTDef;
        gRawNoNoiseGTDef[0] = gRawNoNoise[0] - synthDefC0.interpolate_on01 ( gRawNoNoise );
        gRawNoNoiseGTDef[1] = gRawNoNoise[1] - synthDefC1.interpolate_on01 ( gRawNoNoise );
        gRawNoNoiseGTDef[2] = gRawNoNoise[2] - synthDefC2.interpolate_on01 ( gRawNoNoise );

        //synthDef2DC0.set(x,y, synthDefC0.interpolate_on01( gRawNoNoise ) );
        //synthDef2DC1.set(x,y, synthDefC1.interpolate_on01( gRawNoNoise ) );
        //synthDef2DC2.set(x,y, synthDefC2.interpolate_on01( gRawNoNoise ) );
        //
        //estDef2DC0.set(x,y, synthDefC0.interpolate_on01( gRawNoNoise ) );
        //estDef2DC1.set(x,y, synthDefC1.interpolate_on01( gRawNoNoise ) );
        //estDef2DC2.set(x,y, synthDefC2.interpolate_on01( gRawNoNoise ) );

        RType defL1DiffC0 = - synthDefC0.interpolate_on01 ( gRawNoNoise ) - phi[0].get ( x, y );
        RType defL1DiffC1 = - synthDefC1.interpolate_on01 ( gRawNoNoise ) - phi[1].get ( x, y );
        RType defL1DiffC2 = - synthDefC2.interpolate_on01 ( gRawNoNoise ) - phi[2].get ( x, y );

        RType defL1DiffValue = std::sqrt ( aol::Sqr ( defL1DiffC0 ) + aol::Sqr ( defL1DiffC1 ) + aol::Sqr ( defL1DiffC2 ) );
        DefL1Diff.pushBack ( defL1DiffValue );

        //if( x == 100 && y == 100)
        ////if( defL1DiffValue > 0.001 )
        //{
        //  std::cout << "x,y = " << x << "," << y << std::endl;

        //  std::cout << "synthDefC0 = " << synthDefC0.interpolate_on01( gRawNoNoise ) << std::endl;
        //  std::cout << "synthDefC1 = " << synthDefC1.interpolate_on01( gRawNoNoise ) << std::endl;
        //  std::cout << "synthDefC2 = " << synthDefC2.interpolate_on01( gRawNoNoise ) << std::endl;
        //
        //  std::cout << "synthDefC0/129. = " << synthDefC0.interpolate_on01( gRawNoNoise )/129. << std::endl;
        //  std::cout << "synthDefC1/129. = " << synthDefC1.interpolate_on01( gRawNoNoise )/129. << std::endl;
        //  std::cout << "synthDefC2/129. = " << synthDefC2.interpolate_on01( gRawNoNoise )/129. << std::endl;
        //
        //  std::cout << "gRaw[0] = " << gRaw[0] << std::endl;
        //  std::cout << "gRaw[1] = " << gRaw[1] << std::endl;
        //  std::cout << "gRaw[2] = " << gRaw[2] << std::endl;

        //  std::cout << "gRawNoNoise[0] = " << gRawNoNoise[0] << std::endl;
        //  std::cout << "gRawNoNoise[1] = " << gRawNoNoise[1] << std::endl;
        //  std::cout << "gRawNoNoise[2] = " << gRawNoNoise[2] << std::endl;

        //  std::cout << "gDen[0] = " << gDen[0] << std::endl;
        //  std::cout << "gDen[1] = " << gDen[1] << std::endl;
        //  std::cout << "gDen[2] = " << gDen[2] << std::endl;

        //  std::cout << "phi[0].get( x, y ) = " << phi[0].get( x, y ) << std::endl;
        //  std::cout << "phi[1].get( x, y ) = " << phi[1].get( x, y ) << std::endl;
        //  std::cout << "phi[2].get( x, y ) = " << phi[2].get( x, y ) << std::endl;
        //
        //  std::cout << "defL1DiffC0 = " << defL1DiffC0 << std::endl;
        //  std::cout << "defL1DiffC1 = " << defL1DiffC1 << std::endl;
        //  std::cout << "defL1DiffC2 = " << defL1DiffC2 << std::endl;
        //  std::cout << "defL1DiffValue = " << defL1DiffValue << std::endl;
        //}

        MeshRawNoNoiseDef.pushBackVertex ( gRawNoNoiseDef );
        MeshRawNoNoiseGTDef.pushBackVertex ( gRawNoNoiseGTDef );

        MeshRaw.pushBackVertex ( gRaw );
        MeshDen.pushBackVertex ( gDen );
        MeshDef.pushBackVertex ( gDef );
        MeshDefDen.pushBackVertex ( gDefDen );
      }
    }

    for ( int y = 0; y < tofSizeY - 1; y++ ) {
      for ( int x = 0; x < tofSizeX - 1; x++ ) {
        aol::Vec3<int> UpperFace ( y * tofSizeX + x + 1, y * tofSizeX + x, ( y + 1 ) * ( tofSizeX ) + x );
        aol::Vec3<int> LowerFace ( ( y + 1 ) * ( tofSizeX ) + x, ( y + 1 ) * ( tofSizeX ) + x + 1, y * tofSizeX + x + 1 );

        MeshRaw.pushBackTriang ( UpperFace );
        MeshRaw.pushBackTriang ( LowerFace );

        MeshDen.pushBackTriang ( UpperFace );
        MeshDen.pushBackTriang ( LowerFace );

        MeshDef.pushBackTriang ( UpperFace );
        MeshDef.pushBackTriang ( LowerFace );

        MeshDefDen.pushBackTriang ( UpperFace );
        MeshDefDen.pushBackTriang ( LowerFace );

        MeshRawNoNoiseDef.pushBackTriang ( UpperFace );
        MeshRawNoNoiseDef.pushBackTriang ( LowerFace );

        MeshRawNoNoiseGTDef.pushBackTriang ( UpperFace );
        MeshRawNoNoiseGTDef.pushBackTriang ( LowerFace );
      }
    }


    aol::MeshWithData< aol::TriangMesh<RType> > MeshDefDenSDF ( MeshDefDen );
    MeshDefDenSDF.addData ( DefDenSDF, "sdf", aol::VERTEX_DATA );

    aol::MeshWithData< aol::TriangMesh<RType> > MeshDenSDF ( MeshDen );
    MeshDenSDF.addData ( DenSDF, "sdf", aol::VERTEX_DATA );

    aol::MeshWithData< aol::TriangMesh<RType> > MeshDefL1Diff ( MeshDefDen );
    MeshDefL1Diff.addData ( DefL1Diff, "sdf", aol::VERTEX_DATA );

    aol::MeshWithData< aol::TriangMesh<RType> > MeshRawNoNoiseDefSDF ( MeshRawNoNoiseDef );
    MeshRawNoNoiseDefSDF.addData ( DefSDF, "sdf", aol::VERTEX_DATA );



    char phi0FileName[1024];
    char phi1FileName[1024];
    char phi2FileName[1024];

    char phiFileName[1024];
    char rangeFileName[1024];

    sprintf ( phi0FileName, "%s/data/phi0_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.raw", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( phi1FileName, "%s/data/phi1_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.raw", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( phi2FileName, "%s/data/phi2_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.raw", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );

    sprintf ( phiFileName, "%s/data/phi_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.raw", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( rangeFileName, "%s/data/range_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.raw", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );

    phi0.save ( phi0FileName, qc::PGM_DOUBLE_BINARY );
    phi1.save ( phi1FileName, qc::PGM_DOUBLE_BINARY );
    phi2.save ( phi2FileName, qc::PGM_DOUBLE_BINARY );

    phi.save ( phiFileName, qc::PGM_DOUBLE_BINARY );
    newRangeData.save ( rangeFileName, qc::PGM_DOUBLE_BINARY );


    // save meshes
    char meshRawDataFile[1024];
    sprintf ( meshRawDataFile, "%s/raw/raw_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshRaw.saveAsPLY ( meshRawDataFile );

    char meshRawNoNoiseDefDataFile[1024];
    char MeshRawNoNoiseDefSDFFile[1024];
    sprintf ( meshRawNoNoiseDefDataFile, "%s/nonoisedef/nonoisedef_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( MeshRawNoNoiseDefSDFFile, "%s/eval_nonoisedef/eval_nonoisedef_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.vtk", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshRawNoNoiseDef.saveAsPLY ( meshRawNoNoiseDefDataFile );
    MeshRawNoNoiseDefSDF.saveAsLegacyVTK ( MeshRawNoNoiseDefSDFFile );

    char meshRawNoNoiseGTDefDataFile[1024];
    sprintf ( meshRawNoNoiseGTDefDataFile, "%s/nonoisegtdef/nonoisegtdef_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshRawNoNoiseGTDef.saveAsPLY ( meshRawNoNoiseGTDefDataFile );



    char meshDenDataFile[1024];
    char meshDenSDFDataFile[1024];
    sprintf ( meshDenDataFile, "%s/den/den_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( meshDenSDFDataFile, "%s/eval_den/eval_den_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.vtk", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshDen.saveAsPLY ( meshDenDataFile );
    MeshDenSDF.saveAsLegacyVTK ( meshDenSDFDataFile );

    char meshDefDataFile[1024];
    char meshDefL1DiffDataFile[1024];
    sprintf ( meshDefDataFile, "%s/def/def_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( meshDefL1DiffDataFile, "%s/eval_def/eval_def_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.vtk", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshDef.saveAsPLY ( meshDefDataFile );
    MeshDefL1Diff.saveAsLegacyVTK ( meshDefL1DiffDataFile );

    char meshDefDenDataFile[1024];
    char meshDefDenSDFDataFile[1024];
    sprintf ( meshDefDenDataFile, "%s/defden/defden_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.ply", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    sprintf ( meshDefDenSDFDataFile, "%s/eval_defden/eval_defden_sm_%i_df_%i_a_%f_k_%f_l_%f_m_%f_d1_%f_d2_%f_gd_%i_it_%i.vtk", saveDirectory, smoothingType, differenceType, alpha, kappa, lambda, mu, delta1, delta2, gradientDescent, solverIterations );
    MeshDefDen.saveAsPLY ( meshDefDenDataFile );
    MeshDefDenSDF.saveAsLegacyVTK ( meshDefDenSDFDataFile );
  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}
