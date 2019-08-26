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

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <imageTools.h>
#include <generator.h>
#include <convolution.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

int main(int /*argc*/, char ** /*argv*/)
{
  try {
    RType timeLoad,
          timeAssemble,
          timeRHS,
          timeInvert,
          timeWrite,
          timeInvertOTF;

    aol::StopWatch watch;

    watch.start();
    // Create a grid with 257x257 nodes.
    const ConfType::InitType grid( aol::Vec3<int> ( 257, 257, 1 ) );

    // Create arrays corresponding to the grid.
    ConfType::ArrayType u_old( grid );
    ConfType::ArrayType u_new( grid );
    ConfType::ArrayType rhs( grid );

    qc::DataGenerator<ConfType> generator ( grid );
    generator.generateRectangular ( u_old );

    watch.stop();
    u_old /= u_old.getMaxValue();
    timeLoad = watch.elapsedCpuTime();
    cerr << "Loading time: " << watch.timeToString(timeLoad) << endl;
    watch.start();

    RType tau = 0.01;

    // define the operators
    aol::StiffOp<ConfType> stiff( grid );
    aol::MassOp<ConfType> mass( grid );

    // assembling the system matrix
    qc::FastUniformGridMatrix<RType, ConfType::Dim> mat( qc::GridSize<ConfType::Dim>::createFrom ( grid ) );
    cerr << "Assembling matrices... ";
    stiff.assembleAddMatrix( mat );
    mat *= ( tau );
    mass.assembleAddMatrix( mat );
    cerr << "done\n";
    watch.stop();
    timeAssemble = watch.elapsedCpuTime();
    cerr << "Assembling time: " << watch.timeToString(timeAssemble) << endl;
    watch.start();

    // right hand vector
    mass.apply( u_old, rhs );
    watch.stop();
    timeRHS = watch.elapsedCpuTime();
    cerr << "RHS building time: " << watch.timeToString(timeRHS) << endl;
    watch.start();

    // inverse system matrix
    aol::CGInverse<aol::Vector<RType> > inv( mat );
    inv.setStopping( aol::STOPPING_ABSOLUTE );
    inv.apply( rhs, u_new );
    watch.stop();
    timeInvert = watch.elapsedCpuTime();
    cerr << "Inverting time: " << watch.timeToString(timeInvert) << endl;
    watch.start();

    // output
    qc::writeImage<RType> ( grid, u_new, "1", 1 );
    watch.stop();
    timeWrite = watch.elapsedCpuTime();
    cerr << "Write to disk time: " << watch.timeToString(timeWrite) << endl;
    watch.start();

    // constructing on-the-fly op
    aol:: LinCombOp<aol::Vector<RType> > linCombOp;
    linCombOp.appendReference( stiff, tau );
    linCombOp.appendReference( mass );

    // inverse system matrix
    aol::CGInverse<aol::Vector<RType> > invOTF( linCombOp );
    invOTF.setStopping( aol::STOPPING_ABSOLUTE );
    u_new.setZero();
    invOTF.apply( rhs, u_new );
    watch.stop();
    timeInvertOTF = watch.elapsedCpuTime();
    cerr << "OTF inverting time: " << watch.timeToString(timeInvertOTF) << endl;

    // output
    qc::writeImage<RType> ( grid, u_new, "1OTF", 1 );

    // The continuous solution of the heat equation is the convolution with the Gauss kernel.
    // Note: This produces different results than the two methods above because this gives
    // a solution of the time continous problem whereas the two methods above solve the
    // time discrete problem. You can see this if you approximate the time continious solution
    // by solving the time discrete problem multiple times withs a smaller tau.
    // There are still some differences at the boundary though since LinearConvolutionSmoothOp
    // uses zero extension.
    u_new.setZero();
    qc::LinearConvolutionSmoothOp<RType, ConfType::Dim> smooth ( grid.getNumX(), grid.getNumY() );
    smooth.setTau ( tau );
    smooth.apply ( u_old, u_new );
    qc::writeImage<RType> ( grid, u_new, "1conv", 1 );

  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
