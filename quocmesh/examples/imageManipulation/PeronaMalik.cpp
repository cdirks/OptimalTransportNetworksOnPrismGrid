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
 * \brief Applies Perona–Malik diffusion to an image.
 *
 * \author Berkels
 */

#include <solver.h>
#include <configurators.h>
#include <imageTools.h>
#include <anisoStiffOps.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
//typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3>, qc::CellCenteredCubicGrid<DimensionChoice> > ConfType;

int main ( int argc, char **argv)
{
  try {
    if ( argc != 4 ) {
      cerr << "USAGE: " << argv[0] << " <input_file> <lambda> <tau>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputFileNameBase = aol::getBaseFileName ( inputFileName );
    const RType lambda = atof(argv[2]);
    const RType tau = atof(argv[3]);

    ConfType::ArrayType uOld( inputFileName );

    const ConfType::InitType grid ( uOld.getSize() );

    ConfType::ArrayType uNew( grid );
    ConfType::ArrayType rhs( grid );

    uOld.scaleValuesTo01();

    aol::MassOp<ConfType> mass( grid );

    ConfType::MatrixType mat( qc::GridSize<ConfType::Dim>::createFrom ( grid ) );
    if ( lambda > 0 ) {
      qc::PeronaMalikStiffOp<ConfType> stiff( grid, aol::ONTHEFLY, lambda );
      stiff.setImageReference ( uOld );
      stiff.assembleAddMatrix( mat );
    }
    else if ( lambda < 0 ) {
      qc::AnisotropicDiffusion2DStiffOp<ConfType> stiff( grid, aol::ONTHEFLY, -lambda );
      stiff.computeStructureTensor ( uOld );
      stiff.assembleAddMatrix( mat );
    }
    else {
      aol::StiffOp<ConfType> stiff( grid, aol::ONTHEFLY );
      stiff.assembleAddMatrix( mat );
    }
    mat *= ( tau );
    mass.assembleAddMatrix( mat );

    mass.apply( uOld, rhs );

    aol::CGInverse<aol::Vector<RType> > inv( mat, 1e-16, 10000 );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    inv.apply( rhs, uNew );

    qc::writeImage<RType> ( grid, uNew, aol::strprintf( "%s_%lf_%lf", outputFileNameBase.c_str(), lambda, tau ).c_str(), 1 );
  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
