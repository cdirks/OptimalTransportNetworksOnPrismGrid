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
 * This tool reads a set of 3-D points from a *.bin file (assuming quadratic 
 * grid sampling structure) and generates a 3-D triangle mesh.
 *
 * author: Bauer
 */

#include <parameterParser.h> 
#include <aol.h>
#include <scalarArray.h>
#include <configurators.h>
#include <imageTools.h>
#include <triangMesh.h>

#include <iostream>
#include <vector>


typedef double RType;
 
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

	char PointCloudDataFile[1024]; 
	char PLYMeshFile[1024];
	char VTKMeshFile[1024];
	int CloseSurface = 0;

	if( Parser.hasVariable( "PointCloudDataFile" ))
	{
		Parser.getString( "PointCloudDataFile", PointCloudDataFile );
	}	
	if( Parser.hasVariable( "PLYMeshFile" ))
	{
		Parser.getString( "PLYMeshFile", PLYMeshFile );
	}
	if( Parser.hasVariable( "VTKMeshFile" ))
	{
		Parser.getString( "VTKMeshFile", VTKMeshFile );
	}
	if( Parser.hasVariable( "CloseSurface" ))
	{
		CloseSurface = Parser.getInt( "CloseSurface" );
	}


	// read points (assuming quadratic grid sampling structure)
	std::ifstream InputStream;
	InputStream.open( PointCloudDataFile, std::ios::in | std::ios::binary );
	if ( !InputStream.is_open() || InputStream.fail() ) exit(0);
	aol::Vec3<RType> Point;
	aol::TriangMesh<RType> Mesh;
	while(!(InputStream.eof()))
	{
		InputStream >> Point[0];
		InputStream >> Point[1];
		InputStream >> Point[2];

		Mesh.pushBackVertex( Point );
	}
	InputStream.close();

	int DimX = std::sqrt(static_cast<float> ( Mesh.getNumVertices() ));
	int DimY = DimX;
	//std::cout << "DimX = " << DimX << "; DimY = " << DimY << std::endl;

	for( int y = 0; y < DimY-1; y++ )
	{
		for( int x = 0; x < DimX-1; x++ )
		{	
			aol::Vec3<int> UpperFace( y*DimX + x+1, y*DimX + x, (y+1)*(DimX) + x );
			aol::Vec3<int> LowerFace( (y+1)*(DimX) + x, (y+1)*(DimX) + x+1, y*DimX + x+1 );

			Mesh.pushBackTriang( UpperFace );
			Mesh.pushBackTriang( LowerFace );
		}
	}
	
	for( int y = 0; y < DimY; y++ )
	{
		for( int x = 0; x < DimX; x++ )
		{	
			if( CloseSurface == 1 )
			{
				if( x == 0 || x == DimX-1 || y == 0 || y == DimY-1 )
				{
					aol::Vec3<RType> Vertex = Mesh.getVertex(y*DimX+x);
					Vertex[2] = 1.;
					Mesh.setVertex(y*DimX+x, Vertex);
				}
			}
		}
	}
	
	//std::cout << Mesh.getNumVertices() << std::endl;
	//std::cout << Mesh.getNumTriangs() << std::endl;
	
	Mesh.saveAsPLY( PLYMeshFile );
	Mesh.saveAsLegacyVTK( VTKMeshFile );
}	
catch ( aol::Exception &el ) 
{
    el.dump();
    aol::callSystemPauseIfNecessaryOnPlatform();
}
}





