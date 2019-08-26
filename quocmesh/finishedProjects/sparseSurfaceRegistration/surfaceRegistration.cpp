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
 * \brief  Surface registration
 * \author Bauer
 *
 * \details
 * In this project, two surfaces are registered
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


#include "surfaceRegistration.h"


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
	RType Kappa = 0., Lambda = 0., Eps = 0.;
    if ( argc == 5 ) 
	{
		sprintf( ParameterFilename, "%s",  argv[1] );
		Kappa = atof(argv[2]);
		Lambda = atof(argv[3]);
		Eps = atof(argv[4]);
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
	char GTDeformedReferenceDataSDFFile[1024];
	char GTDeformedReferenceDataMeshFile[1024];	
	char SaveDirectory[1024];
	char ReferenceDataMeshFile[1024];
	int GridSize2DX = 0, GridSize2DY = 0;
	//int ValidateDerivative = 0;
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
	//if( Parser.hasVariable( "ValidateDerivative" ))
	//{
	//	ValidateDerivative = Parser.getInt( "ValidateDerivative" );
	//}
	if( Parser.hasVariable( "SolverIterations" ))
	{
		SolverIterations = Parser.getInt( "SolverIterations" );
	}

	// load template data Y (here a dense set of points)
	aol::TriangMesh<RType> TemplateMesh;
	TemplateMesh.loadFromPLY( GTDeformedReferenceDataMeshFile );

	aol::MultiVector<RType> Y( ConfType3D::Dim, GridSize2DX*GridSize2DY );	
	for( int x = 0; x < GridSize2DX; x++ )
	{
		for( int y = 0; y < GridSize2DY; y++ )
		{
			Vec3D Gi = TemplateMesh.getVertex( y*GridSize2DX+x );
			Y[0][y*GridSize2DX+x] = Gi[0];
			Y[1][y*GridSize2DX+x] = Gi[1];
			Y[2][y*GridSize2DX+x] = Gi[2];
		}
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
	

	// u
	qc::MultiArray<RType, 2, 3> DispU( Grid2D );
	qc::MultiArray<RType, 2, 3> DispU_tmp( Grid2D );
	DispU.setAll(0.);
	DispU_tmp.setAll(0.);


	/*
	 * set up energy functionals
	 */

	//! build the regularization matrix
	aol::SparseMatrix<RType> tempRegMat ( Grid2D );
	qc::assembleLaplaceSquareRegMatrix<ConfType2D> ( Grid2D, tempRegMat );
	// Since the matrix will stay the same (and in particular keep its sparsity pattern), we can
	// use a aol::CSR_Matrix here that can apply itself noticeably faster than a aol::SparseMatrix.
	aol::CSR_Matrix<> regMat ( tempRegMat );
		
	SurfaceRegistrationEnergy<ConfType2D, ConfType3D, aol::CSR_Matrix<> > E( Grid2D, Grid3D, ReferenceDataSDF, Y, regMat, Kappa, Lambda );
	VariationOfSurfaceRegistrationEnergy<ConfType2D, ConfType3D, aol::CSR_Matrix<> > DE( Grid2D, Grid3D, ReferenceDataSDF, Y, regMat, Kappa, Lambda );

	std::cout << "Kappa  = " << Kappa  << std::endl;
	std::cout << "Lambda = " << Lambda << std::endl;

	const bool NoConsoleOutput = false;
	typedef aol::GradientDescentWithAutomaticFilterWidth<ConfType2D, aol::MultiVector<RType>, aol::GradientDescent< ConfType2D, aol::MultiVector<RType> > > GDType;	
	GDType gradientDescent_solver ( Grid2D, E, DE, SolverIterations, aol::ZOTrait< RType >::one, Eps );
	gradientDescent_solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG | ( NoConsoleOutput ? GDType::DO_NOT_WRITE_CONSOLE_OUTPUT : 0 ) );	
	std::cout << "starting solver..." << std::endl;
	gradientDescent_solver.apply ( DispU, DispU_tmp );
	std::cout << "done." << std::endl;
	DispU = DispU_tmp;

	// save u
	char DispU_Dir[1024];
	sprintf ( DispU_Dir, "%s/u/", SaveDirectory );
	aol::makeDirectory( DispU_Dir );
	char DispUFile[1024];
	sprintf ( DispUFile, "%s/u_k_%f_l_%f_it_%i_%f.bin", DispU_Dir, Kappa, Lambda, SolverIterations, Eps );
	DispU.saveASCII( DispUFile );

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

	aol::Vector<RType> d_psiY_G;

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_phiG_TEX_d_phiG_phiGTG( Mesh_phiG );
	Mesh_phiG_TEX_d_phiG_phiGTG.addData( d_phiG_phiGTG, "d_phiG_phiGTG", aol::VERTEX_DATA );

	aol::MeshWithData< aol::TriangMesh<RType> > Mesh_phiGTG_TEX_d_G_phiGTG( Mesh_phiGTG );
	Mesh_phiGTG_TEX_d_G_phiGTG.addData( d_G_phiGTG, "d_G_phiGTG", aol::VERTEX_DATA );

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
	sprintf ( Mesh_phiG_FileName, "%s/phi(G)_k_%f_l_%f_it_%i_%f.vtk", Mesh_phiG_Dir, Kappa, Lambda, SolverIterations, Eps );
	Mesh_phiG.saveAsLegacyVTK( Mesh_phiG_FileName );

	char Mesh_phiG_TEX_d_phiG_phiGTG_Dir[1024];
	sprintf ( Mesh_phiG_TEX_d_phiG_phiGTG_Dir, "%s/phi(G)_TEX_d(phi(G),phiGT(G))/", SaveDirectory );
	aol::makeDirectory( Mesh_phiG_TEX_d_phiG_phiGTG_Dir );
	char Mesh_phiG_TEX_d_phiG_phiGTG_FileName[1024];
	sprintf ( Mesh_phiG_TEX_d_phiG_phiGTG_FileName, "%s/phi(G)_TEX_d(phi(G),phiGT(G))_k_%f_l_%f_it_%i_%f.vtk", Mesh_phiG_TEX_d_phiG_phiGTG_Dir, Kappa, Lambda, SolverIterations, Eps );
	Mesh_phiG_TEX_d_phiG_phiGTG.saveAsLegacyVTK( Mesh_phiG_TEX_d_phiG_phiGTG_FileName );

	char Mesh_phiGTG_TEX_d_G_phiGTG_FileName[1024];
	sprintf ( Mesh_phiGTG_TEX_d_G_phiGTG_FileName, "%s/phiGT(G)_TEX_d(G,phiGT(G)).vtk", SaveDirectory );
	Mesh_phiGTG_TEX_d_G_phiGTG.saveAsLegacyVTK( Mesh_phiGTG_TEX_d_G_phiGTG_FileName );

	char Mesh_G_phiGlyph_Dir[1024];
	sprintf ( Mesh_G_phiGlyph_Dir, "%s/G_phiGlyph/", SaveDirectory );
	aol::makeDirectory( Mesh_G_phiGlyph_Dir );
	char Mesh_G_phiGlyph_FileName[1024];
	sprintf ( Mesh_G_phiGlyph_FileName, "%s/G_phiGlyph_k_%f_l_%f_it_%i_%f.vtk", Mesh_G_phiGlyph_Dir, Kappa, Lambda, SolverIterations, Eps );
	Mesh_G_phiGlyph.saveAsLegacyVTK( Mesh_G_phiGlyph_FileName );

}
catch ( aol::Exception &el ) 
{
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
}
	return 0;
}

