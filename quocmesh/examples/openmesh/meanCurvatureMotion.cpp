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
 * \file
 * \brief Executes several steps of mean curvature motion.
 *
 * Executes several steps of mean curvature motion. Parameter files have the format:
 * \begin{verbatim}
 * inputFilename      [...].ply
 * outputFilename     [...].ply
 * fileFormat         { UDply | ply }       (optional)
 * numSteps           [int]
 * tau                [float]
 * \end{verbatim}
 *
 * \author Heeren, von Deylen
 */

#include <triMesh.h>
#include <triangMeshConfigurators.h>
#include <preconditioner.h>
#include <solver.h>

typedef double RealType;
typedef om::TriMesh<double> MeshType;
typedef aol::TriangMeshConfigurator<RealType, MeshType, aol::CenterQuadrature<RealType> > ConfType;

void mcmStep ( MeshType & mesh, RealType tau ) {

  const int size = mesh.getNumVertices();

  aol::LumpedMassOp<ConfType> lumpedMass ( mesh, aol::DO_NOT_INVERT  );
  aol::StiffOp<ConfType> stiffOp ( mesh, aol::ONTHEFLY  );

  aol::MultiVector<RealType> surf ( 3, size );
  aol::MultiVector<RealType> surf_old ( 3, size );
  aol::MultiVector<RealType> rhs ( 3, size );

  mesh.toVector ( surf );
  surf_old = surf;

  aol::SparseMatrix<RealType> mat ( size, size );
  stiffOp.assembleAddMatrix ( mat );
  mat *= ( tau );
  lumpedMass.assembleAddMatrix ( mat );

  aol::DiagonalBlockOp<RealType> ( lumpedMass ).applyAdd ( surf, rhs );
  aol::SSORPreconditioner<aol::Vector<RealType>, aol::SparseMatrix<RealType> > precond ( mat );

  aol::PCGInverse<aol::Vector<RealType> > cg ( mat, precond, 1e-24, 1000 );
  cg.setStopping ( aol::STOPPING_ABSOLUTE );
  for ( int i = 0; i < 3; i++ ) {
    cg.apply ( rhs[i], surf_old[i] );
  }
  surf_old -= surf;
  //om::NormalProjection<ConfType> nProj ( mesh );
  //nProj.applyAdd( surf_old, surf );
  surf += surf_old;

  mesh.fromVector ( surf );
}

int main ( int argc, char * argv[] ) {

  try {
    // input from user
    string parameterFileName;
    if ( argc < 2 ) {
        cout << "Enter 'parmfile.txt' or create your own parameterfile (by using 'parmfile.txt' as template): ";
        cin >> parameterFileName;
    } else parameterFileName = argv[1];

    aol::ParameterParser parser( parameterFileName.c_str() );

    string inputFilename  = parser.getString ( "inputFilename" );
    string outputFilename = parser.getString ( "outputFilename" );
    int numSteps          = parser.getInt    ( "numSteps" );
    RealType tau          = parser.getDouble ( "tau" );

    string fileFormat = "UDply";
    if ( parser.hasVariable ( "fileFormat" ) )
      fileFormat          = parser.getString ( "fileFormat" );

    MeshType mesh;
    clog << "Loading " << inputFilename << " ... ";
    if ( fileFormat == "UDply" )
      mesh.loadFromUDPLY ( inputFilename );
    else if ( fileFormat == "ply" )
      mesh.loadFromPLY ( inputFilename );
    else
      throw aol::Exception ( ( string ( "TriMesh load function only implemented "
              "for file types UDply and ply, not for " ) + fileFormat ).c_str() );
    clog << "done." << endl;

    clog << "Processing MCM steps ";

    for ( int i = 0; i < numSteps; ++i ) {
      mcmStep ( mesh, tau );
      clog << ".";
    }
    clog << endl;

    clog << "Saving as " << outputFilename << " ... ";
    if ( fileFormat == "UDply" )
      mesh.saveAsUDPLY ( outputFilename );
    else
      mesh.saveAsPLY ( outputFilename );
    clog << "done." << endl;
  }
  catch ( aol::Exception & exc ) {

    cerr << aol::color::red << "Caught aol::Exception(): ";
    exc.dump();
    cerr << endl;
  }
}
