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


#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>


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


#include <projector.h>

// typedefs
typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType2D;
typedef qc::RectangularGridConfigurator<RType, qc::QC_3D, aol::GaussQuadrature<RType,qc::QC_3D,3> > ConfType3D;
typedef ConfType3D::DomVecType Vec3D;

struct Point3D
{
	float X, Y, Z;
};


int main( int argc, char **argv )
{	
try 
{
	// read parameterfile
    char ParameterFilename[1024];
    if ( argc == 2 ) 
	{
		sprintf( ParameterFilename, "%s",  argv[1] );		
    }
	else
	{
		cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
		return EXIT_FAILURE;
	}

	cerr << "Reading parameters from " << ParameterFilename << endl;
	aol::ParameterParser Parser( ParameterFilename );
	
	char SDFFile[1024];
	char SDFLevelSetMeshFile[1024];
	char SaveDirectory[1024];
	
	if( Parser.hasVariable( "SDFFile" ))
	{
		Parser.getString( "SDFFile", SDFFile );
	}
	if( Parser.hasVariable( "SDFLevelSetMeshFile" ))
	{
		Parser.getString( "SDFLevelSetMeshFile", SDFLevelSetMeshFile );
	}
	if( Parser.hasVariable( "SaveDirectory" ))
	{
		Parser.getString( "SaveDirectory", SaveDirectory );
	}

	// load reference data SDF
	qc::ScalarArray<RType, qc::QC_3D> SDF( SDFFile );
	const ConfType3D::InitType Grid3D ( qc::GridSize<qc::QC_3D>::createFrom ( SDF ) );
	const aol::DiscreteFunctionDefault<ConfType3D> SDFDiscFuncs( Grid3D, SDF );

	aol::TriangMesh<RType> SDFLevelSetMesh;
	SDFLevelSetMesh.loadFromPLY( SDFLevelSetMeshFile );

	Projector<ConfType3D> Proj( Grid3D, SDFDiscFuncs );

	ConfType3D Configurator3D( Grid3D );
	typedef ConfType3D::ElementType ElementType3D;
	typedef ConfType3D::DomVecType DomVecType3D;
	ElementType3D PEl, POnGEl;
	DomVecType3D PCoord, POnGCoord;

	char DistancesFileName[1024];
	sprintf(DistancesFileName, "%s_dist.txt", SDFLevelSetMeshFile);
	std::ofstream OutputStream(DistancesFileName);
	

	int c = 0;
	RType acc = 0;
	for ( int i = 0; i < SDFLevelSetMesh.getMaxNumVertices(); i++ )
	{
		Vec3D P = SDFLevelSetMesh.getVertex( i );
		//if (	P[0] > 0.2 && P[0] < 0.8 &&
		//		P[1] > 0.2 && P[1] < 0.8 &&
		//		P[2] > 0.2 && P[2] < 0.8 )
		{
			if ( Configurator3D.getLocalCoords( P, PEl, PCoord ) )
			{
				// project
				DomVecType3D POnG = Proj.projectOnG( P, PEl, PCoord );
				// evaluate
				if ( Configurator3D.getLocalCoords( POnG, POnGEl, POnGCoord ) )
				{
					//RType dist1 = SDFDiscFuncs.evaluate( PEl, PCoord );
					RType dist2 = SDFDiscFuncs.evaluate( POnGEl, POnGCoord );
					//std::cout << dist1 << "\t->\t" << dist2 << std::endl;
					c++;
					acc += fabs(dist2);

					OutputStream << dist2 << "\n";
				}
			}
		}
	}

	//// twofold projection
	//ElementType3D POnGEl2;
	//DomVecType3D POnGCoord2;
	//for ( int i = 0; i < SDFLevelSetMesh.getMaxNumVertices(); i++ )
	//{
	//	Vec3D P = SDFLevelSetMesh.getVertex( i );
	//	//if (	P[0] > 0.2 && P[0] < 0.8 &&
	//	//		P[1] > 0.2 && P[1] < 0.8 &&
	//	//		P[2] > 0.2 && P[2] < 0.8 )
	//	{
	//		if ( Configurator3D.getLocalCoords( P, PEl, PCoord ) )
	//		{
	//			// project
	//			DomVecType3D POnG = Proj.projectOnG( P, PEl, PCoord );
	//			// evaluate
	//			if ( Configurator3D.getLocalCoords( POnG, POnGEl, POnGCoord ) )
	//			{
	//				// project a second time
	//				DomVecType3D POnG2 = Proj.projectOnG( POnG, POnGEl, POnGCoord );

	//				if ( Configurator3D.getLocalCoords( POnG2, POnGEl2, POnGCoord2 ) )
	//				{
	//					//RType dist1 = SDFDiscFuncs.evaluate( PEl, PCoord );
	//					RType dist2 = SDFDiscFuncs.evaluate( POnGEl2, POnGCoord2 );
	//					//std::cout << dist1 << "\t->\t" << dist2 << std::endl;
	//					c++;
	//					acc += fabs(dist2);

	//					OutputStream << dist2 << "\n";
	//				}
	//			}
	//		}
	//	}
	//}

	OutputStream.close();

	std::cout << "mean = " << acc/c << std::endl;		

}
catch ( aol::Exception &el ) 
{
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
}
	return 0;
}

