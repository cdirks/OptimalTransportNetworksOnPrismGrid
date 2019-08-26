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
 * \brief Registers two surfaces (e.g. TOF/CT) using signed distance map as similarity measure
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
    char deformedReferenceDataFile[1024];
    char saveDirectory[1024];
    int tofSizeX = 0, tofSizeY = 0;
    RType focalLength = 0.;
    RType scaleFactor = 0., shiftOffsetX = 0., shiftOffsetY = 0., shiftOffsetZ = 0.;
    int solverIterations = 0;
    RType kappa = 0., lambda = 0.;
    int validateDerivative = 0;
    int gradientDescent = 0;

    if ( parser.hasVariable ( "templateDataFile" ) ) {
      parser.getString ( "templateDataFile", templateDataFile );
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
    if ( parser.hasVariable ( "kappa" ) ) {
      kappa = parser.getReal<RType> ( "kappa" );
    }
    if ( parser.hasVariable ( "lambda" ) ) {
      lambda = parser.getReal<RType> ( "lambda" );
    }
    if ( parser.hasVariable ( "validateDerivative" ) ) {
      validateDerivative = parser.getInt ( "validateDerivative" );
    }
    if ( parser.hasVariable ( "gradientDescent" ) ) {
      gradientDescent = parser.getInt ( "gradientDescent" );
    }

    //// load template, TOF range data
    std::ifstream RangeDataInputStream;
    RangeDataInputStream.open ( templateDataFile, std::ios::in | std::ios::binary );
    qc::ScalarArray<RType, qc::QC_2D> templateRangeData ( tofSizeX, tofSizeY );
    templateRangeData.loadRaw ( RangeDataInputStream, qc::PGM_FLOAT_BINARY, tofSizeX, tofSizeY );
    RangeDataInputStream.close();

    // Benefit joint approach
    //qc::ScalarArray<RType, qc::QC_2D> templateRangeData( templateDataFile );

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

    // similarity measure
    MatchingEnergyOfPhi<ConfType2D, ConfType3D> Energy ( grid2D, grid3D, deformedReferenceSDF, transformation, templateRangeData, 0. );
    VariationOfMatchingEnergyWRTPhi<ConfType2D, ConfType3D> EnergyVariation ( grid2D, grid3D, deformedReferenceSDF, transformation, templateRangeData );

    //MatchingEnergyOfPhi<ConfType2D, ConfType3D> Energy( grid2D, grid3D, referenceSDF, transformation, templateRangeData );
    //VariationOfMatchingEnergyWRTPhi<ConfType2D, ConfType3D> EnergyVariation( grid2D, grid3D, referenceSDF, transformation, templateRangeData );

    // Assemble stiffness matrix L
    //qc::FastUniformGridMatrix<RType, ConfType2D::Dim> mat( qc::GridSize<qc::QC_2D>::createFrom ( grid2D ) );
    aol::StiffOp<ConfType2D> stiff ( grid2D, aol::ASSEMBLED );
    //stiff.assembleAddMatrix( mat );

    // regularization
    const qc::DisplacementLengthEnergy<ConfType2D> Regularizer ( stiff );
    const aol::DiagonalBlockOp<RType> RegularizerVariation ( stiff );

    // symmetric
    //const qc::SymmetricLengthEnergy<ConfType2D> Regularizer( grid2D );
    //const qc::VariationOfSymmetricLengthEnergy<ConfType2D> RegularizerVariation( grid2D );

    // combine energies
    aol::LinCombOp<aol::MultiVector<RType>, aol::Scalar<RType> > E;
    aol::LinCombOp<aol::MultiVector<RType> > DE;

    E.appendReference ( Energy, kappa );
    DE.appendReference ( EnergyVariation, kappa );

    E.appendReference ( Regularizer, lambda );
    DE.appendReference ( RegularizerVariation, lambda );


    // Deformation to compute
    // <RType, 2, 3> integrieren auf 2D, 3 Unbekannte (Deformation)
    qc::MultiArray<RType, 2, 3> mtmp ( grid2D );
    // Initial deformation for solver initialization
    // Phi is initialized with 0's by the aol::MultiVector constructor
    qc::MultiArray<RType, 2, 3> Phi ( grid2D );

    // Validate first derivative
    if ( validateDerivative == 1 ) {
      qc::DataGenerator<ConfType2D> generator ( grid2D );
      generator.generateNonLinearDeformation ( 0.05, mtmp );
      aol::FirstDerivativeValidator<aol::MultiVector<RType> > tester ( E, DE, grid2D.H(), aol::FirstDerivativeValidator<aol::MultiVector<RType> >::LINEAR, 1 );

      char derivativeValidation[1024];
      sprintf ( derivativeValidation, "%sDerivativeValidation.ply", saveDirectory );
      tester.testDirection ( mtmp, derivativeValidation );
      exit ( 0 );
    }

    //solver
    //double tau = 0.01;
    const bool NoConsoleOutput = false;

    //aol::SimpleGradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, DE, 100);
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, GradientDescentType> GDType;
    //typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType, aol::MultiVector<RType>, aol::GradientDescent<ConfType, aol::MultiVector<RType> >> GDType;

    typedef aol::GradientDescent<ConfType2D, aol::MultiVector<RType> > GDType1;
    typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::MultiVector<RType>, aol::H1GradientDescent< ConfType2D, aol::MultiVector<RType> > > GDType2;

    //GDType gradientDescent_solver ( grid, E, DE, 1000 , tau );
    GDType1 gradientDescent_solver1 ( grid2D, E, DE, solverIterations );
    GDType2 gradientDescent_solver2 ( grid2D, E, DE, solverIterations );


    std::cout << "starting solver..." << std::endl;

    if ( gradientDescent == 1 ) {
      std::cout << "GDType1." << std::endl;
      gradientDescent_solver1.setConfigurationFlags ( GDType1::USE_NONLINEAR_CG | GDType1::DO_NOT_SMOOTH_DESCENT_DIRECTION | ( NoConsoleOutput ? GDType1::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );
      gradientDescent_solver1.apply ( Phi, mtmp );
    } else if ( gradientDescent == 2 ) {
      std::cout << "GDType2." << std::endl;
      gradientDescent_solver2.setConfigurationFlags ( GDType2::USE_NONLINEAR_CG | ( NoConsoleOutput ? GDType2::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );
      gradientDescent_solver2.apply ( Phi, mtmp );
    }
    std::cout << "solver finished." << std::endl;

    Phi = mtmp;
    //tau = gradientDescent_solver.getStartTau();


    Phi.save ( "../../../Data/test2.raw", qc::PGM_DOUBLE_BINARY );

    // apply the estimated deformation on the TOF 3D data and save result as PLY mesh
    // in addition, for convenience, save original TOF mesh
    aol::TriangMesh<RType> Mesh;
    aol::TriangMesh<RType> MeshDef;
    aol::TriangMesh<RType> MeshDefGT;

    aol::Vector<RType> DefL2Diff;
    aol::Vector<RType> DefSDF;


    qc::ScalarArray<RType, qc::QC_3D> synthDefC0 ( syntheticDeformation[0], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );
    qc::ScalarArray<RType, qc::QC_3D> synthDefC1 ( syntheticDeformation[1], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );
    qc::ScalarArray<RType, qc::QC_3D> synthDefC2 ( syntheticDeformation[2], grid3D.getNumX(), grid3D.getNumY(), grid3D.getNumZ() );

    for ( int y = 0; y < tofSizeY; y++ ) {
      for ( int x = 0; x < tofSizeX; x++ ) {
        RType r = templateRangeData.get ( x, y );
        aol::Vec3<RType> g = transformation.evaluateWorldCoord ( x, y, r );

        aol::Vec3<RType> gDef;
        gDef[0] = g[0] + Phi[0].get ( x, y );
        gDef[1] = g[1] + Phi[1].get ( x, y );
        gDef[2] = g[2] + Phi[2].get ( x, y );

        aol::Vec3<RType> gDefGT;
        gDefGT[0] = g[0] - synthDefC0.interpolate_on01 ( g );
        gDefGT[1] = g[1] - synthDefC1.interpolate_on01 ( g );
        gDefGT[2] = g[2] - synthDefC2.interpolate_on01 ( g );

        Mesh.pushBackVertex ( g );
        MeshDef.pushBackVertex ( gDef );
        MeshDefGT.pushBackVertex ( gDefGT );

        RType defL2DiffC0 = - synthDefC0.interpolate_on01 ( g ) - Phi[0].get ( x, y );
        RType defL2DiffC1 = - synthDefC1.interpolate_on01 ( g ) - Phi[1].get ( x, y );
        RType defL2DiffC2 = - synthDefC2.interpolate_on01 ( g ) - Phi[2].get ( x, y );

        if ( x == 25 && y == 150 ) {
          std::cout << "synthDefC0.interpolate_on01( g ) = " << synthDefC0.interpolate_on01 ( g ) << std::endl;
          std::cout << "synthDefC1.interpolate_on01( g ) = " << synthDefC1.interpolate_on01 ( g ) << std::endl;
          std::cout << "synthDefC2.interpolate_on01( g ) = " << synthDefC2.interpolate_on01 ( g ) << std::endl;

          std::cout << "synthDefC0.interpolate( g ) = " << synthDefC0.interpolate_on01 ( g ) << std::endl;
          std::cout << "synthDefC1.interpolate( g ) = " << synthDefC1.interpolate_on01 ( g ) << std::endl;
          std::cout << "synthDefC2.interpolate( g ) = " << synthDefC2.interpolate_on01 ( g ) << std::endl;

          std::cout << "Phi[0].get( x, y ) = " << Phi[0].get ( x, y ) << std::endl;
          std::cout << "Phi[1].get( x, y ) = " << Phi[1].get ( x, y ) << std::endl;
          std::cout << "Phi[2].get( x, y ) = " << Phi[2].get ( x, y ) << std::endl;

          std::cout << "defL2DiffC0 = " << defL2DiffC0 << std::endl;
          std::cout << "defL2DiffC1 = " << defL2DiffC1 << std::endl;
          std::cout << "defL2DiffC2 = " << defL2DiffC2 << std::endl;
        }

        RType defL2DiffValue = std::sqrt ( aol::Sqr ( defL2DiffC0 ) + aol::Sqr ( defL2DiffC1 ) + aol::Sqr ( defL2DiffC2 ) );
        DefL2Diff.pushBack ( defL2DiffValue );

        RType defSDFValue = deformedReferenceSDF.interpolate_on01 ( gDef );
        DefSDF.pushBack ( defSDFValue );
      }
    }

    for ( int y = 0; y < tofSizeY - 1; y++ ) {
      for ( int x = 0; x < tofSizeX - 1; x++ ) {
        aol::Vec3<int> UpperFace ( y * tofSizeX + x + 1, y * tofSizeX + x, ( y + 1 ) * ( tofSizeX ) + x );
        aol::Vec3<int> LowerFace ( ( y + 1 ) * ( tofSizeX ) + x, ( y + 1 ) * ( tofSizeX ) + x + 1, y * tofSizeX + x + 1 );

        Mesh.pushBackTriang ( UpperFace );
        Mesh.pushBackTriang ( LowerFace );

        MeshDef.pushBackTriang ( UpperFace );
        MeshDef.pushBackTriang ( LowerFace );

        MeshDefGT.pushBackTriang ( UpperFace );
        MeshDefGT.pushBackTriang ( LowerFace );
      }
    }

    aol::MeshWithData< aol::TriangMesh<RType> > MeshDefL2Diff ( MeshDef );
    MeshDefL2Diff.addData ( DefL2Diff, "L2Diff", aol::VERTEX_DATA );


    aol::MeshWithData< aol::TriangMesh<RType> > MeshDefSDF ( MeshDef );
    MeshDefSDF.addData ( DefSDF, "sdf", aol::VERTEX_DATA );

    //std::cout << Mesh.getNumVertices() << std::endl;
    //std::cout << Mesh.getNumTriangs() << std::endl;

    // save meshes
    char deformedMeshDataFile[1024];
    sprintf ( deformedMeshDataFile, "%s/def/def_k_%f_l_%f_gd_%i_it_%i.ply", saveDirectory, kappa, lambda, gradientDescent, solverIterations );
    MeshDef.saveAsPLY ( deformedMeshDataFile );

    char deformedMeshGTDataFile[1024];
    sprintf ( deformedMeshGTDataFile, "%s/defGT_k_%f_l_%f_gd_%i_it_%i.ply", saveDirectory, kappa, lambda, gradientDescent, solverIterations );
    MeshDefGT.saveAsPLY ( deformedMeshGTDataFile );

    //char originalMeshDataFile[1024];
    //sprintf ( originalMeshDataFile, "%soriginalMesh.ply", saveDirectory );
    //Mesh.saveAsPLY( originalMeshDataFile );

    char meshDefL2DiffDataFile[1024];
    sprintf ( meshDefL2DiffDataFile, "%s/eval_def_L2Diff/eval_def_L2Diff_k_%f_l_%f_gd_%i_it_%i.vtk", saveDirectory, kappa, lambda, gradientDescent, solverIterations );
    MeshDefL2Diff.saveAsLegacyVTK ( meshDefL2DiffDataFile );

    char meshDefSDFDataFile[1024];
    sprintf ( meshDefSDFDataFile, "%s/eval_def_SDF/eval_def_SDF_k_%f_l_%f_gd_%i_it_%i.vtk", saveDirectory, kappa, lambda, gradientDescent, solverIterations );
    MeshDefSDF.saveAsLegacyVTK ( meshDefSDFDataFile );
  } catch ( aol::Exception &el ) {
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  return 0;
}

