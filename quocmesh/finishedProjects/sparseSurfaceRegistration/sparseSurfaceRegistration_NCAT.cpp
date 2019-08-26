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
 * \brief  Sparse surface registration
 * \author Bauer
 *
 * \details
 * In this project, a sparse set of 3-D measurements, denoted \Y (template), is 
 * registered with a dense surface \G (reference).
 * 
 */

#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

// AOL
#include <aol.h>
#include <gradientDescent.h>
#include <gradientflow.h>
#include <FEOpInterface.h>
#include <meshWithData.h>
#include <parameterParser.h> 
#include <triangMesh.h>
#include <Newton.h>

// QUOC
#include <configurators.h>
#include <fastUniformGridMatrix.h>
#include <hyperelastic.h>
#include <generator.h>
#include <imageTools.h>
#include <linearSmoothOp.h>
#include <multiArray.h>
#include <scalarArray.h>
#include <solver.h>

// EIKONAL
#include <signedDistanceOp.h>

// TPCFE
#include <signedDistanceSweeping.h>


#include "sparseSurfaceRegistration.h"

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
void testFENonlinSumInterface ( ) {
  typedef typename ConfiguratorType::RealType RealType;
  const int numPositions = 2;
  aol::MultiVector<RealType> positions ( 2, numPositions );
  aol::Vector<RealType> values ( numPositions );
  positions[0][0] = 0.5;
  positions[1][0] = 0.25;
  positions[0][1] = 0.5;
  positions[1][1] = 0.75;
  values[0] = 1;
  values[1] = 0;

  typename ConfiguratorType::InitType grid ( aol::Vec3<int> ( 17, 17, 1 ) );
  TestFittingEnergy <ConfiguratorType> EFit ( grid, positions, values );
  FENonlinSumInterfaceTest<ConfiguratorType> DEFit ( grid, positions, values );

  aol::SparseMatrix<RealType> regMat;
  qc::assembleLaplaceSquareRegMatrix<ConfiguratorType> ( grid, regMat );
  aol::DiagonalBlockOp<RealType> DEReg ( regMat );
  aol::QuadraticFormOp<aol::MultiVector<RealType> > EReg ( DEReg );

  aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
  const RealType lambda = 0.01;
  E.appendReference ( EFit );
  E.appendReference ( EReg, lambda );

  aol::LinCombOp<aol::MultiVector<RealType> > DE;
  DE.appendReference ( DEFit );
  DE.appendReference ( DEReg, lambda );

  // Test if DE fits as derivative of E
  qc::ScalarArray<RealType, qc::QC_2D> u ( grid );
  aol::MultiVector<RealType> uMVec;
  uMVec.appendReference ( u );
  aol::makeDirectory ( "test" );
  aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, grid.H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.0001 );
  tester.testDirection ( uMVec, "test/test" );
  tester.setSkipPlottingWithSmallError( true );
  tester.setSkippingThreshold( 1e-8 );
  tester.testAllDirections( uMVec, "test/test" );

  // Find a minimizer of E with a gradient descent.
  qc::ScalarArray<RealType, qc::QC_2D> tmp ( grid );
  aol::MultiVector<RealType> tmpMVec;
  tmpMVec.appendReference ( tmp );
  typedef aol::GradientDescent<ConfiguratorType, aol::MultiVector<RealType> > GDType;
  GDType gradientDescent_solver ( grid, E, DE, 1000 );
  gradientDescent_solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG );
  gradientDescent_solver.apply ( uMVec, tmpMVec );
  tmp.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
  tmp.savePNG ( "test/out.png" );
}


// typedefs
typedef double RType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType2D;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_3D, aol::GaussQuadrature<RType,qc::QC_3D,3> > ConfType3D;
typedef ConfType3D::DomVecType Vec3D;

struct Point3D
{
	float X, Y, Z;
};


int main( int argc, char **argv )
{	
try 
{
	// 0: Validation
	// 1: NCAT Evaluation
	int MODE = 1;

	// read parameterfile
    char ParameterFilename[1024];
	RType Kappa = 0., Lambda = 0., Mu = 0.;
	RType Eps = 0.;
    if ( argc == 6 ) 
	{
		sprintf( ParameterFilename, "%s",  argv[1] );
		Kappa = atof(argv[2]);
		Lambda = atof(argv[3]);
		Mu = atof(argv[4]);
		Eps = atof(argv[5]);
    }
	else
	{
		cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
		return EXIT_FAILURE;
	}

	cerr << "Reading parameters from " << ParameterFilename << endl;
	aol::ParameterParser Parser( ParameterFilename );
	
	char TemplateDataFile[1024]; 
	char ReferenceDataFile[1024];
	char ReferenceDataSamplingFile[1024];
	//char DeformedReferenceDataFile[1024];
	char GTDeformedReferenceDataSDFFile[1024];
	char GTDeformedReferenceDataMeshFile[1024];	
	char SaveDirectory[1024];
	char ReferenceDataMeshFile[1024];
	int GridSize2DX = 0, GridSize2DY = 0;
	int ValidateDerivative = 0;
	//RType Kappa = 0., Lambda = 0., Mu = 0.;
	int SolverIterations = 0;

	if( Parser.hasVariable( "TemplateDataFile" ))
	{
		Parser.getString( "TemplateDataFile", TemplateDataFile );
	}	
	if( Parser.hasVariable( "ReferenceDataFile" ))
	{
		Parser.getString( "ReferenceDataFile", ReferenceDataFile );
	}
	if( Parser.hasVariable( "ReferenceDataSamplingFile" ))
	{
		Parser.getString( "ReferenceDataSamplingFile", ReferenceDataSamplingFile );
	}
	//if( Parser.hasVariable( "DeformedReferenceDataFile" ))
	//{
	//	Parser.getString( "DeformedReferenceDataFile", DeformedReferenceDataFile );
	//}
	if( Parser.hasVariable( "GTDeformedReferenceDataSDFFile" ))
	{
		Parser.getString( "GTDeformedReferenceDataSDFFile", GTDeformedReferenceDataSDFFile );
	}
	if( Parser.hasVariable( "GTDeformedReferenceDataMeshFile" ))
	{
		Parser.getString( "GTDeformedReferenceDataMeshFile", GTDeformedReferenceDataMeshFile );
	}
	if( Parser.hasVariable( "SaveDirectory" ))
	{
		Parser.getString( "SaveDirectory", SaveDirectory );
	}
	if( Parser.hasVariable( "ReferenceDataMeshFile" ))
	{
		Parser.getString( "ReferenceDataMeshFile", ReferenceDataMeshFile );
	}
	if( Parser.hasVariable( "GridSize2DX" ))
	{
		GridSize2DX = Parser.getInt( "GridSize2DX" );		
	}
	if( Parser.hasVariable( "GridSize2DY" ))
	{
		GridSize2DY = Parser.getInt( "GridSize2DY" );
	}
	if( Parser.hasVariable( "ValidateDerivative" ))
	{
		ValidateDerivative = Parser.getInt( "ValidateDerivative" );
	}
	//if( Parser.hasVariable( "Kappa" ))
	//{
	//	Kappa = Parser.getReal<RType>( "Kappa" );
	//}
	//if( Parser.hasVariable( "Lambda" ))
	//{
	//	Lambda = Parser.getReal<RType>( "Lambda" );
	//}
	//if( Parser.hasVariable( "Mu" ))
	//{
	//	Lambda = Parser.getReal<RType>( "Mu" );
	//}
	if( Parser.hasVariable( "SolverIterations" ))
	{
		SolverIterations = Parser.getInt( "SolverIterations" );
	}

	// create directory for results
	//char ResultsDir[1024];
	//sprintf ( ResultsDir, "%s/k_%f_l_%f_m_%f_it_%i", SaveDirectory, Kappa, Lambda, Mu, SolverIterations );
	//aol::makeDirectory( ResultsDir );

	// load sparse template data S
	std::ifstream InputStreamS;
	InputStreamS.open( ReferenceDataSamplingFile, std::ios::in | std::ios::binary );
	if ( !InputStreamS.is_open() || InputStreamS.fail() )
	{
		std::cout << "Failed to read S." << std::endl;
		exit(0);
	}
	std::vector<Point3D> ReferenceDataSamplingPoints;
	Point3D Point;
	while(!(InputStreamS.eof()))
	{
		InputStreamS >> Point.X;
		InputStreamS >> Point.Y;
		InputStreamS >> Point.Z;
		ReferenceDataSamplingPoints.push_back(Point);
	}
	InputStreamS.close();

	// setup MultiVector for sparse template data Y
	aol::MultiVector<RType> S( ConfType3D::Dim, ReferenceDataSamplingPoints.size() );	
	for( unsigned int i = 0; i < ReferenceDataSamplingPoints.size(); i++ )
	{
		S[0][i] = ReferenceDataSamplingPoints[i].X;
		S[1][i] = ReferenceDataSamplingPoints[i].Y;
		S[2][i] = ReferenceDataSamplingPoints[i].Z;
	}

	// load sparse template data Y
	std::ifstream InputStream;
	InputStream.open( TemplateDataFile, std::ios::in | std::ios::binary );
	if ( !InputStream.is_open() || InputStream.fail() )
	{
		std::cout << "Failed to read Y." << std::endl;
		exit(0);
	}
	std::vector<Point3D> TemplateDataPoints;
	while(!(InputStream.eof()))
	{
		InputStream >> Point.X;
		InputStream >> Point.Y;
		InputStream >> Point.Z;
		TemplateDataPoints.push_back(Point);
	}
	InputStream.close();

	// setup MultiVector for sparse template data Y
	aol::MultiVector<RType> Y( ConfType3D::Dim, TemplateDataPoints.size() );	
	for( unsigned int i = 0; i < TemplateDataPoints.size(); i++ )
	{
		Y[0][i] = TemplateDataPoints[i].X;
		Y[1][i] = TemplateDataPoints[i].Y;
		Y[2][i] = TemplateDataPoints[i].Z;
	}
	
	// load reference data SDF
	qc::ScalarArray<RType, qc::QC_3D> ReferenceDataSDF( ReferenceDataFile );

	// generate 2D, 3D grid
	qc::GridSize<qc::QC_2D> Grid2DSize( GridSize2DX, GridSize2DY );
	const ConfType2D::InitType Grid2D ( Grid2DSize );
	const ConfType3D::InitType Grid3D ( qc::GridSize<qc::QC_3D>::createFrom ( ReferenceDataSDF ) );

	// generate nonlinear deformation
	qc::DataGenerator<ConfType3D> Deformer( Grid3D );
	aol::MultiVector<RType> SyntheticDeformation( ConfType3D::Dim, Grid3D.getNumberOfNodes() );	
	RType DeformationStrength = 0.1;
	Deformer.generateNonLinearDeformation ( DeformationStrength, SyntheticDeformation );
	qc::ScalarArray<RType, qc::QC_3D> SyntheticDeformationC0( SyntheticDeformation[0], Grid3D.getNumX(), Grid3D.getNumY(), Grid3D.getNumZ());
	qc::ScalarArray<RType, qc::QC_3D> SyntheticDeformationC1( SyntheticDeformation[1], Grid3D.getNumX(), Grid3D.getNumY(), Grid3D.getNumZ());
	qc::ScalarArray<RType, qc::QC_3D> SyntheticDeformationC2( SyntheticDeformation[2], Grid3D.getNumX(), Grid3D.getNumY(), Grid3D.getNumZ());
	
	// deform reference data SDF
	qc::ScalarArray<RType, qc::QC_3D> DeformedReferenceDataSDF (ReferenceDataFile);
	if ( MODE == 0 )
	{
		qc::DeformImage<ConfType3D> ( ReferenceDataSDF, Grid3D, DeformedReferenceDataSDF, SyntheticDeformation );
		char DeformedReferenceDataFile[1024];
		sprintf ( DeformedReferenceDataFile, "%s/Deformed_SDF.raw", SaveDirectory );
		DeformedReferenceDataSDF.save( DeformedReferenceDataFile, qc::PGM_DOUBLE_BINARY );
	}
	else if ( MODE == 1 )
	{
		qc::ScalarArray<RType, qc::QC_3D> GTDeformedReferenceDataSDF( GTDeformedReferenceDataSDFFile );
		DeformedReferenceDataSDF = GTDeformedReferenceDataSDF;
	}
	

	if ( MODE == 0 )
	{
		// deform template data
		for( unsigned int i = 0; i < TemplateDataPoints.size(); i++ )
		{
			Vec3D yi;
			Y.getTo( i, yi );
			yi[0] -= SyntheticDeformationC0.interpolate_on01( yi );
			yi[1] -= SyntheticDeformationC1.interpolate_on01( yi );
			yi[2] -= SyntheticDeformationC2.interpolate_on01( yi );
			Y.set( i, yi );
		}
	}

	// u
	qc::MultiArray<RType, 2, 3> DispU( Grid2D );
	qc::MultiArray<RType, 2, 3> DispU_tmp( Grid2D );
	DispU.setAll(0.);
	DispU_tmp.setAll(0.);

	// w
	aol::MultiVector<RType> DispW( 3, TemplateDataPoints.size() );
	aol::MultiVector<RType> DispW_tmp( 3, TemplateDataPoints.size() );
	DispW.setAll(0.);
	DispW_tmp.setAll(0.);
	
	// concatenate u, w
	aol::MultiVector<RType> DispUDispW( 0, 0);
	DispUDispW.appendReference( DispU[0] );
	DispUDispW.appendReference( DispU[1] );
	DispUDispW.appendReference( DispU[2] );
	DispUDispW.appendReference( DispW[0] );
	DispUDispW.appendReference( DispW[1] );
	DispUDispW.appendReference( DispW[2] );

	aol::MultiVector<RType> DispUDispW_tmp( 0, 0);
	DispUDispW_tmp.appendReference( DispU_tmp[0] );
	DispUDispW_tmp.appendReference( DispU_tmp[1] );
	DispUDispW_tmp.appendReference( DispU_tmp[2] );
	DispUDispW_tmp.appendReference( DispW_tmp[0] );
	DispUDispW_tmp.appendReference( DispW_tmp[1] );
	DispUDispW_tmp.appendReference( DispW_tmp[2] );

	/*
	 * set up energy functionals
	 */

	//! build the regularization matrix
	aol::SparseMatrix<RType> tempRegMat ( Grid2D );
	qc::assembleLaplaceSquareRegMatrix<ConfType2D> ( Grid2D, tempRegMat );
	// Since the matrix will stay the same (and in particular keep its sparsity pattern), we can
	// use a aol::CSR_Matrix here that can apply itself noticeably faster than a aol::SparseMatrix.
	aol::CSR_Matrix<> regMat ( tempRegMat );
		
	SparseSurfaceRegistrationEnergy<ConfType2D, ConfType3D, aol::CSR_Matrix<> > E( Grid2D, Grid3D, ReferenceDataSDF, Y, regMat, Kappa, Lambda, Mu, DispUDispW_tmp );
	VariationOfSparseSurfaceRegistrationEnergy<ConfType2D, ConfType3D, aol::CSR_Matrix<> > DE( Grid2D, Grid3D, ReferenceDataSDF, Y, regMat, Kappa, Lambda, Mu, DispUDispW_tmp );

	// Validate first derivative
	if ( ValidateDerivative == 1 ) 
	{
		qc::DataGenerator<ConfType2D> Generator ( Grid2D );
		Generator.generateNonLinearDeformation ( 0.05, DispU );

		aol::MultiVector<RType> TestDispW;
		TestDispW.appendReference(Y);
		//TestDispW *= 0.05;
		
		aol::FirstDerivativeValidator<aol::MultiVector<RType> > tester ( E, DE, Grid2D.H(), aol::FirstDerivativeValidator<aol::MultiVector<RType> >::LINEAR, 1 );
		char derivativeValidation[1024];
		sprintf ( derivativeValidation, "%sDerivativeValidation", SaveDirectory );
		std::cout << "validate..." << std::endl;
		tester.testDirection ( DispUDispW, derivativeValidation );	
		//tester.testDirection ( DispW, derivativeValidation );	
		std::cout << "done." << std::endl;
		exit(0);
	}

	// \E_{match}
	//MatchingEnergy<ConfType3D> EMatch( Grid3D, ReferenceDataSDF, TemplateData );
	//DerivativeOfMatchingEnergy<ConfType3D> DEMatch( Grid3D, ReferenceDataSDF, TemplateData );
	
	// \E_{corr}
	//CorrelationEnergy<ConfType2D, ConfType3D> ECorr( Grid2D, Grid3D, DeformedReferenceDataSDF, TemplateData );
	//VariationOfCorrelationEnergyWRTDispU<ConfType2D, ConfType3D> VECorrWRTDispU( Grid2D, Grid3D, DeformedReferenceDataSDF, TemplateData, DispW );
	//DerivativeOfCorrelationEnergyWRTDispW<ConfType2D, ConfType3D> DECorrWRTDispW( Grid2D, Grid3D, DeformedReferenceDataSDF, TemplateData, DispU );

	//aol::SparseMatrix<RType> regMat ( Grid2D );
	//assembleLaplaceSquareRegMatrix<ConfType2D> ( Grid2D, regMat );
	//aol::DiagonalBlockOp<RType> DEReg ( regMat );
	//aol::QuadraticFormOp<aol::MultiVector<RType> > EReg ( DEReg );


	// Validate first derivative E_{match}
#if 0
	if ( ValidateDerivative == 1 ) 
	{
		qc::DataGenerator<ConfType2D> Generator ( Grid2D );
		Generator.generateNonLinearDeformation ( 0.05, DispU );
		aol::FirstDerivativeValidator<aol::MultiVector<RType> > tester ( EMatch, DEMatch, Grid2D.H(), aol::FirstDerivativeValidator<aol::MultiVector<RType> >::LINEAR, 1 );
		
		char derivativeValidation[1024];
		sprintf ( derivativeValidation, "%sDerivativeValidation", SaveDirectory );
		tester.testDirection ( DispU, derivativeValidation );	
		exit(0);
	}
#endif 

	std::cout << "Kappa  = " << Kappa  << std::endl;
	std::cout << "Lambda = " << Lambda << std::endl;
	std::cout << "Mu     = " << Mu     << std::endl;

	const bool NoConsoleOutput = false;
	aol::Vector<int> tauBlocks(2);
	tauBlocks[0] = 3;
	tauBlocks[1] = 3;

	
	unsigned int OptimizationMethod = 0;
	switch(OptimizationMethod)
	{
	case 0: //GradientDescent	
		{
			//typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::MultiVector<RType>, aol::GradientDescent< ConfType2D, aol::MultiVector<RType> > > GDType;	
			typedef aol::GradientDescent<ConfType2D, aol::MultiVector<RType> > GDType;
			typedef aol::GradientDescentComponentWiseTimestepControlled<ConfType2D, GDType> GDComponentWiseType;		
			
			GDComponentWiseType gradientDescent_solver ( Grid2D, E, DE, tauBlocks, SolverIterations, aol::ZOTrait< RType >::one, Eps );
			//gradientDescent_solver.setConfigurationFlags ( GDComponentWiseType::DO_NOT_SMOOTH_DESCENT_DIRECTION | GDComponentWiseType::USE_NONLINEAR_CG | GDComponentWiseType::CALC_STOPPING_ENERGY_BEFORE_UPDATING_CURRENT_POSITION | ( NoConsoleOutput ? GDComponentWiseType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );	
			gradientDescent_solver.setConfigurationFlags ( GDComponentWiseType::USE_NONLINEAR_CG | GDComponentWiseType::CALC_STOPPING_ENERGY_BEFORE_UPDATING_CURRENT_POSITION | ( NoConsoleOutput ? GDComponentWiseType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );	
			// only smooth the u part of the unknowns.
			gradientDescent_solver.setSmoothFirstNComponents ( 3 );
			
			aol::StopWatch watch;
			watch.start();
			
			gradientDescent_solver.apply ( DispUDispW, DispUDispW_tmp );
			
			watch.stop();
			printf("Solver Runtime: %ld milliseconds\n", static_cast<long int> ( watch.elapsedWallClockTime()*1000 ) );
			
			break;
		}
	case 1: // QuasiNewton
		{
			typedef aol::QuasiNewtonIterationComponentWiseTimestepControlled<RType> QNComponentWiseType;
			aol::NewtonInfo<RType> newtonInfo ( Eps, SolverIterations, Eps, SolverIterations, aol::STOPPING_ABSOLUTE );
			newtonInfo.setTimestepController ( aol::NewtonInfo<RType>::ARMIJO );  
			//QNComponentWiseType quasiNewton_solver( E, DE, tauBlocks, newtonInfo );
			QNComponentWiseType quasiNewton_solver( E, DE, tauBlocks, SolverIterations, Eps );
			
			std::cout << "starting solver..." << std::endl;
			double solverStart = aol::getRuntimeSoFar();
			quasiNewton_solver.apply ( DispUDispW, DispUDispW_tmp );
			double solverStop = aol::getRuntimeSoFar();
			std::cout << "solver runtime: " << ((solverStop - solverStart)*1000.) << std::endl;
		
			break;
		}
	case 2: //QuasiNewton BFGS
		{
			typedef aol::QuasiNewtonBFGS<RType,aol::MultiVector<RType>,aol::MultiVector<RType> > QNBFGSType;
			aol::NewtonInfo<RType> newtonInfo ( Eps, SolverIterations, Eps, SolverIterations, aol::STOPPING_ABSOLUTE );
			newtonInfo.setTimestepController ( aol::NewtonInfo<RType>::ARMIJO );  
			QNBFGSType quasiNewton_solver( E, DE, newtonInfo );
			
			std::cout << "starting solver..." << std::endl;
			double solverStart = clock();
			quasiNewton_solver.apply ( DispUDispW, DispUDispW_tmp );
			double solverStop = clock();
			std::cout << "solver runtime: " << ((solverStop - solverStart)*1000.) << std::endl;
			
			break;
		}
	default:
		throw aol::UnimplementedCodeException ( "Unknown optimization method", __FILE__, __LINE__ );
	}
	
	DispUDispW = DispUDispW_tmp;

	// save u
	char DispU_Dir[1024];
	sprintf ( DispU_Dir, "%s/u/", SaveDirectory );
	aol::makeDirectory( DispU_Dir );
	char DispUFile[1024];
	sprintf ( DispUFile, "%s/u_k_%f_l_%f_m_%f_it_%i_%f.bin", DispU_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	DispU.saveASCII( DispUFile );

	// save w
	char DispW_Dir[1024];
	sprintf ( DispW_Dir, "%s/w/", SaveDirectory );
	aol::makeDirectory( DispW_Dir );
	char DispWFile[1024];
	sprintf ( DispWFile, "%s/w_k_%f_l_%f_m_%f_it_%i_%f.bin", DispW_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	DispW.saveASCII( DispWFile );


	/*
	 * Evaluation
	 */

	aol::TriangMesh<RType> Mesh_G;
	Mesh_G.loadFromPLY( ReferenceDataMeshFile );

	aol::TriangMesh<RType> Mesh_phiGTG;
	if ( MODE == 1 )
	{
		Mesh_phiGTG.loadFromPLY( GTDeformedReferenceDataMeshFile );
	}

	aol::TriangMesh<RType> Mesh_phiG;
	aol::Vector<RType> d_phiG_phiGTG;
	aol::Vector<RType> d_G_phiGTG;

	aol::MultiVector<RType> phiGlyph( 3, GridSize2DX*GridSize2DY+1 );
	phiGlyph.setAll( 0. );

	int i = 0;
	for( int x = 0; x < GridSize2DX; x++ )
	{
		for( int y = 0; y < GridSize2DY; y++ )
		{
			Vec3D Gi = Mesh_G.getVertex( i );
			Vec3D phiGi = Gi + DispU.get( i );
			Mesh_phiG.pushBackVertex( phiGi );
			d_phiG_phiGTG.pushBack( DeformedReferenceDataSDF.interpolate_on01( phiGi ) );

			d_G_phiGTG.pushBack( DeformedReferenceDataSDF.interpolate_on01( Gi ) );

			if ( MODE == 0 )
			{
				Gi[0] -= SyntheticDeformationC0.interpolate_on01( Gi );
				Gi[1] -= SyntheticDeformationC1.interpolate_on01( Gi );
				Gi[2] -= SyntheticDeformationC2.interpolate_on01( Gi );
				Mesh_phiGTG.pushBackVertex( Gi );
			}

			Vec3D u_i = DispU.get( i );
			phiGlyph.set(i, u_i);

			i++;
		}
	}

	// find out what's going wrong before
	d_G_phiGTG.pushBack( 0. );

	// triangulation
	for( int x = 0; x < GridSize2DX-1; x++ )
	{
		for( int y = 0; y < GridSize2DY-1; y++ )
		{
			aol::Vec3<int> UpperFace(y*GridSize2DX + x+1, y*GridSize2DX + x, (y+1)*(GridSize2DX) + x);
			aol::Vec3<int> LowerFace((y+1)*(GridSize2DX) + x, (y+1)*(GridSize2DX) + x+1, y*GridSize2DX + x+1);

			Mesh_phiG.pushBackTriang( UpperFace );
			Mesh_phiG.pushBackTriang( LowerFace );
			
			if ( MODE == 0 )
			{
				Mesh_phiGTG.pushBackTriang( UpperFace );
				Mesh_phiGTG.pushBackTriang( LowerFace );
			}
		}
	}

	aol::TriangMesh<RType> Mesh_S;
	aol::TriangMesh<RType> Mesh_Y;
	aol::TriangMesh<RType> Mesh_psiY;
	aol::Vector<RType> d_psiY_G;

	for( unsigned int i = 0; i < TemplateDataPoints.size(); i++ )
	{
		aol::Vec3<RType> w_i;
		DispW.getTo( i, w_i );

		aol::Vec3<RType> s_i;
		S.getTo( i, s_i );		
		Mesh_S.pushBackVertex( s_i );
		
		aol::Vec3<RType> y_i;
		Y.getTo( i, y_i );		
		Mesh_Y.pushBackVertex( y_i );

		Mesh_psiY.pushBackVertex( y_i + w_i ); 
		d_psiY_G.pushBack( ReferenceDataSDF.interpolate_on01( y_i + w_i ) );
	}

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_phiG_TEX_d_phiG_phiGTG( Mesh_phiG );
	Mesh_phiG_TEX_d_phiG_phiGTG.addData( d_phiG_phiGTG, "d_phiG_phiGTG", aol::VERTEX_DATA );

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_phiGTG_TEX_d_G_phiGTG( Mesh_phiGTG );
	Mesh_phiGTG_TEX_d_G_phiGTG.addData( d_G_phiGTG, "d_G_phiGTG", aol::VERTEX_DATA );

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_psiY_TEX_d_psiY_G( Mesh_psiY );
	Mesh_psiY_TEX_d_psiY_G.addData( d_psiY_G, "d_psiY_G", aol::VERTEX_DATA );

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_G_phiGlyph( Mesh_G );
	Mesh_G_phiGlyph.addData( phiGlyph, "phiGlyph", aol::VERTEX_DATA );
	


	// save meshes

	char Mesh_phiGTG_FileName[1024];
	sprintf ( Mesh_phiGTG_FileName, "%s/phiGT(G).vtk", SaveDirectory );
	Mesh_phiGTG.saveAsLegacyVTK( Mesh_phiGTG_FileName );

	char Mesh_phiG_Dir[1024];
	sprintf ( Mesh_phiG_Dir, "%s/phi(G)/", SaveDirectory );
	aol::makeDirectory( Mesh_phiG_Dir );
	char Mesh_phiG_FileName[1024];
	sprintf ( Mesh_phiG_FileName, "%s/phi(G)_k_%f_l_%f_m_%f_it_%i_%f.vtk", Mesh_phiG_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	Mesh_phiG.saveAsLegacyVTK( Mesh_phiG_FileName );

	char Mesh_phiG_TEX_d_phiG_phiGTG_Dir[1024];
	sprintf ( Mesh_phiG_TEX_d_phiG_phiGTG_Dir, "%s/phi(G)_TEX_d(phi(G),phiGT(G))/", SaveDirectory );
	aol::makeDirectory( Mesh_phiG_TEX_d_phiG_phiGTG_Dir );
	char Mesh_phiG_TEX_d_phiG_phiGTG_FileName[1024];
	sprintf ( Mesh_phiG_TEX_d_phiG_phiGTG_FileName, "%s/phi(G)_TEX_d(phi(G),phiGT(G))_k_%f_l_%f_m_%f_it_%i_%f.vtk", Mesh_phiG_TEX_d_phiG_phiGTG_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	Mesh_phiG_TEX_d_phiG_phiGTG.saveAsLegacyVTK( Mesh_phiG_TEX_d_phiG_phiGTG_FileName );

	char Mesh_phiGTG_TEX_d_G_phiGTG_FileName[1024];
	sprintf ( Mesh_phiGTG_TEX_d_G_phiGTG_FileName, "%s/phiGT(G)_TEX_d(G,phiGT(G)).vtk", SaveDirectory );
	Mesh_phiGTG_TEX_d_G_phiGTG.saveAsLegacyVTK( Mesh_phiGTG_TEX_d_G_phiGTG_FileName );

	char Mesh_G_phiGlyph_Dir[1024];
	sprintf ( Mesh_G_phiGlyph_Dir, "%s/G_phiGlyph/", SaveDirectory );
	aol::makeDirectory( Mesh_G_phiGlyph_Dir );
	char Mesh_G_phiGlyph_FileName[1024];
	sprintf ( Mesh_G_phiGlyph_FileName, "%s/G_phiGlyph_k_%f_l_%f_m_%f_it_%i_%f.vtk", Mesh_G_phiGlyph_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	Mesh_G_phiGlyph.saveAsLegacyVTK( Mesh_G_phiGlyph_FileName );

	char Mesh_S_FileName[1024];
	sprintf ( Mesh_S_FileName, "%s/S.vtk", SaveDirectory );
	Mesh_S.saveAsLegacyVTK( Mesh_S_FileName );

	char Mesh_Y_FileName[1024];
	sprintf ( Mesh_Y_FileName, "%s/Y.vtk", SaveDirectory );
	Mesh_Y.saveAsLegacyVTK( Mesh_Y_FileName );

	char Mesh_psiY_Dir[1024];
	sprintf ( Mesh_psiY_Dir, "%s/psi(Y)/", SaveDirectory );
	aol::makeDirectory( Mesh_psiY_Dir );
	char Mesh_psiY_FileName[1024];
	sprintf ( Mesh_psiY_FileName, "%s/psi(Y)_k_%f_l_%f_m_%f_it_%i_%f.ply", Mesh_psiY_Dir, Kappa, Lambda, Mu, SolverIterations, Eps );
	Mesh_psiY.saveAsPLY( Mesh_psiY_FileName );

}
catch ( aol::Exception &el ) 
{
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
}
	return 0;
}

