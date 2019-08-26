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
 * \brief Denoises ToF data.
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

#include <timestepSaver.h>
#include <firstOrderTVAlgos.h>

#include <iostream>

#include "TOF.h"
#include "SurfaceDenoising.h"


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
    char referenceDataFile[1024];
    char saveDirectory[1024];
    int tofSizeX = 0, tofSizeY = 0;
    RType focalLength = 0.;
    RType scaleFactor = 0., shiftOffset = 0.;
    int solverIterations  = 0;
    int smoothingType  = 1;
    int differenceType = 1;
    RType kappa = 0., lambda = 0., mu = 0.;
    RType delta1 = 0., delta2 = 0.;
    //int validateDerivative = 0;

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
    if ( parser.hasVariable ( "shiftOffset" ) ) {
      shiftOffset = parser.getReal<RType> ( "shiftOffset" );
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
    //if( parser.hasVariable( "validateDerivative" ))
    //{
    //  validateDerivative = parser.getInt( "validateDerivative" );
    //}

    typedef RType RealType;
    typedef ConfType2D ConfiguratorType;

    const RealType gamma = 1. / kappa;

    //ConfiguratorType::InitType grid ( qc::getSizeFromArrayFile( parser.getString( "input-image" ) ) );
    // generate 2D, 3D grid
    qc::GridSize<qc::QC_2D> grid2DSize ( tofSizeX, tofSizeY );
    const ConfType2D::InitType grid2D ( grid2DSize );

    //ConfiguratorType::ArrayType r0 ( grid2D );
    //u0.load( parser.getString( "input-image" ).c_str() );
    //u0 /= u0.getMaxValue();

    // load template, TOF range data
    std::ifstream InputStreamTOFTemplate;
    InputStreamTOFTemplate.open ( templateDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> templateRangeData ( tofSizeX, tofSizeY );
    templateRangeData.loadRaw ( InputStreamTOFTemplate, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    RType scaling = templateRangeData.getMaxValue();
    templateRangeData /= scaling;
    InputStreamTOFTemplate.close();

    qc::FirstOrderPrimalDualROFMinimizer<ConfiguratorType> tvAlgo ( grid2D, gamma, templateRangeData, solverIterations );
    ConfiguratorType::ArrayType newRangeData ( templateRangeData );
    tvAlgo.minimize ( newRangeData );

    newRangeData *= scaling;

    // range to world coord transformation
    RangeToWorldCoordTransformation<RType> transformation ( focalLength, scaleFactor, shiftOffset );


    // apply the estimated deformation on the TOF 3D data and save result as PLY mesh
    // in addition, for convenience, save original TOF mesh
    aol::TriangMesh<RType> MeshNoisy;
    aol::TriangMesh<RType> MeshDenoised;

    for ( int y = 0; y < tofSizeY; y++ ) {
      for ( int x = 0; x < tofSizeX; x++ ) {
        //RType rNoisy = templateRangeData.get( x, y );
        //aol::Vec3<RType> gNoisy = transformation.evaluateWorldCoord( x, y, rNoisy );

        RType rDenoised = newRangeData.get ( x, y );
        aol::Vec3<RType> gDenoised = transformation.evaluateWorldCoord ( x, y, rDenoised );

        //MeshNoisy.pushBackVertex( gNoisy );
        MeshDenoised.pushBackVertex ( gDenoised );
      }
    }

    for ( int y = 0; y < tofSizeY - 1; y++ ) {
      for ( int x = 0; x < tofSizeX - 1; x++ ) {
        aol::Vec3<int> UpperFace ( y * tofSizeX + x + 1, y * tofSizeX + x, ( y + 1 ) * ( tofSizeX ) + x );
        aol::Vec3<int> LowerFace ( ( y + 1 ) * ( tofSizeX ) + x, ( y + 1 ) * ( tofSizeX ) + x + 1, y * tofSizeX + x + 1 );

        //MeshNoisy.pushBackTriang( UpperFace );
        //MeshNoisy.pushBackTriang( LowerFace );

        MeshDenoised.pushBackTriang ( UpperFace );
        MeshDenoised.pushBackTriang ( LowerFace );
      }
    }

    //std::cout << Mesh.getNumVertices() << std::endl;
    //std::cout << Mesh.getNumTriangs() << std::endl;

    // save meshes
    //char meshNoisyDataFile[1024];
    //sprintf ( meshNoisyDataFile, "%smeshNoisy.ply", saveDirectory );
    //MeshNoisy.saveAsPLY( meshNoisyDataFile );

    char meshDenoisedDataFile[1024];
    sprintf ( meshDenoisedDataFile, "%s/den_sm_%i_df_%i_k_%f_l_%f_m_%f_d1_%f_d2_%f_it_%i.ply", saveDirectory, smoothingType, differenceType, kappa, lambda, mu, delta1, delta2, solverIterations );
    MeshDenoised.saveAsPLY ( meshDenoisedDataFile );
  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}

