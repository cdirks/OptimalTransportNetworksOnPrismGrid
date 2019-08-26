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

// include all eikonal header files
#include<chanVeseDescent.h>
#include<eikonal3d.h>
#include<eikonal.h>
#include<eikonalNA3d.h>
#include<eikonalNA.h>
#include<functionOnLevelSet.h>
#include<segDefs.h>
#include<segment3d.h>
#include<segTrialNode.h>
#include<signedDistanceOp.h>

#include <iostream>

#include <configurators.h>
#include <deformations.h>
#include <signedDistanceSweeping.h>


void testSuccess( const bool success, const char *errorText ) {
  if ( success ) {
    cerr << "OK." << endl;
  } else {
    cerr << aol::color::red << errorText << endl << aol::color::reset;
  }
}

template <typename ConfiguratorType>
bool testSignedDistanceOp ( const typename ConfiguratorType::InitType &Grid, const string ConfigName ) {
  cerr << "--- Testing qc::SignedDistanceOp (" << ConfigName << ") ... " ;
  qc::DataGenerator<ConfiguratorType> generator ( Grid );

  // Generate a signed distance function of a circle.
  typename ConfiguratorType::ArrayType cirlceSignedDistFunc ( Grid );
  generator.generateCircleLevelset ( cirlceSignedDistFunc, 0.15 );

  typename ConfiguratorType::ArrayType temp ( cirlceSignedDistFunc );
  temp *= 20.;

  qc::SignedDistanceOp<ConfiguratorType> signedDistOp( Grid );
  signedDistOp.apply( temp, temp );
  // Calculate the L2-norm of the difference between the original signed distance function and the one obtained by redistancing.
  temp -= cirlceSignedDistFunc;
  temp *= 1./temp.size();

  if( temp.norm() < 1e-5 ){
    cerr << "OK." << endl;
    return true;
  }
  else {
    cerr << "FAILED!" << endl;
    return false;
  }
}

int main( int, char** ) {

  try {
    bool success = true;

    {
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D, 3> > ConfType;
      qc::GridDefinition grid ( 8, ConfType::Dim );
      success &= testSignedDistanceOp<ConfType>( grid, "QuocConfiguratorTraitMultiLin" );
    }

    {
      typedef qc::RectangularGridConfigurator<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfType;

      // To ensure that SignedDistanceOp works with all kind of RectangularGrids, test it
      // with two different non-quadratic grids. One with numY > numX and one with numX > numY.
      aol::Vec3<int> sizeVec1 ( 250, 350, 1 );
      ConfType::InitType  grid1 ( sizeVec1 );
      success &= testSignedDistanceOp<ConfType>( grid1, "RectangularGridConfigurator-Y>X" );

      aol::Vec3<int> sizeVec2 ( 350, 250, 1 );
      ConfType::InitType  grid2 ( sizeVec2 );
      success &= testSignedDistanceOp<ConfType>( grid2, "RectangularGridConfigurator-X>Y" );
    }

    {
      cerr << "--- Testing qc::SignedDistanceOp3D with a plane ... " ;
      typedef qc::RectangularGridConfigurator<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D, 3> > ConfType3D;
      ConfType3D::InitType grid( aol::Vec3<int> ( 35, 20, 25 ) );
      qc::ScalarArray<double, qc::QC_3D> phi0( grid );                 // the computed sphere (exact)
      qc::ScalarArray<double, qc::QC_3D> phiSD( grid );                // the computed signed distance function

      // generate a plane
      for ( int i=0; i<phi0.getNumX(); i++ ) {
        for ( int j=0; j<phi0.getNumY(); j++ ) {
          for ( int k=0; k<phi0.getNumZ(); k++ ) {
            double x = (i * grid.H()) - 0.33;
            phi0.set( i,j,k, x );                               // 0-isosurface is a plane
          }
        }
      }
      // now compute the signed distance function and the error between the original and the computed image
      qc::SignedDistanceOpTrait<ConfType3D>::OpType( grid ).apply( phi0, phiSD );
      // cerr << "phi0.norm():  " << phi0.norm() << ", phiSD.norm: " << phiSD.norm() << ", difference: ";
      phi0 -= phiSD;
      phi0 *= 1./phi0.size();
      success &= phi0.norm() < 1e-10;
      testSuccess( success, "Computing signed distance function of a plane failed!" );

      cerr << "--- Testing qc::SignedDistanceOp3D with a sphere ... " ;
      // generate a sphere with radius 0.247
      const double sphereRadius = 0.247;
      for ( int i=0; i<phi0.getNumX(); i++ ) {
        for ( int j=0; j<phi0.getNumY(); j++ ) {
          for ( int k=0; k<phi0.getNumZ(); k++ ) {
            double x = (i * grid.H()) - 0.5;
            double y = (j * grid.H()) - 0.5;
            double z = (k * grid.H()) - 0.5;
            double value = sqrt( x*x + y*y + z*z );
            phi0.set( i,j,k, value-sphereRadius );               // 0-isosurface is sphere with radius sphereRadius
          }
        }
      }

      // now recalculate the signed distance function
      qc::SignedDistanceOpTrait<ConfType3D>::OpType( grid ).apply( phi0, phiSD );

      phi0 -= phiSD;
      phi0 *= 1./phi0.size();
      // The quality of the signed distance function function obviously has to depend on h.
      success &= phi0.norm() < aol::Sqr ( grid.H() );
      if(!success) cerr << "norm: " << phi0.norm() << endl;
      testSuccess( success, "Computing signed distance function of a sphere failed!" );

      cerr << "--- Testing qc::SignedDistanceOp3D and SignedDistanceSweepingOp3D with an ellipsoid  ... " ;
      // generate an ellipsoid
      // phi0 is no longer a signed distance function to the interface, because isosurfaces different from 0-isosurface are not ellipsoid shaped
      for ( int i=0; i<phi0.getNumX(); i++ ) {
        for ( int j=0; j<phi0.getNumY(); j++ ) {
          for ( int k=0; k<phi0.getNumZ(); k++ ) {
            double x = (i * grid.H()) - 0.5;
            double y = (j * grid.H()) - 0.5;
            double z = (k * grid.H()) - 0.5;
            double value = sqrt( 0.75 * x*x + 0.6* y*y + z*z );
            phi0.set( i,j,k, value - 0.247 );
          }
        }
      }
      // calculate the signed distance function using Fast Marching
      qc::SignedDistanceOpTrait<ConfType3D>::OpType( grid ).apply( phi0, phiSD );

      //claculate signed distance function using sweeping
      qc::SignedDistanceSweepingOp3D<double>( ).applySingle( phi0 );

      phi0 -= phiSD;
      phi0 *= 1./phi0.size();
      // The quality of the signed distance function function obviously has to depend on h.
      success &= phi0.norm() < aol::Sqr ( grid.H() );
      if(!success) cerr << "norm: " << phi0.norm() << endl;
      testSuccess( success, "Computing signed distance function of an ellipsoid failed!" );

    }

    {
      cerr << "--- Testing qc::SignedDistanceSweepingOp3D with a degenerate case  ... " ;

      const int nx = 25, ny = 25, nz = 25;
      const double binExtent = 2.0 / nx;

      qc::ScalarArray<double, qc::QC_3D> data ( nx, ny, nz );
      for ( int i = 0; i < nz; ++i ) {
        for ( int j = 0; j < ny; ++j ) {
          for ( int k = 0; k < nx; ++k ) {
            const double distanceFromCenter = sqrt ( ( i * binExtent - 2 ) * ( i * binExtent - 2 ) + ( j * binExtent ) * ( j * binExtent ) + ( k * binExtent ) * ( k * binExtent ) );
            const double distanceFromCenter2 = sqrt ( ( i * binExtent ) * ( i * binExtent ) + ( j * binExtent - 1 ) * ( j * binExtent - 1 ) + ( k * binExtent ) * ( k * binExtent ) );
            if ( distanceFromCenter2 <= 1.0 || distanceFromCenter <= 1.0 )
              data.set ( k, j, i, 1.0 );
            else
              data.set ( k, j, i, 0 );
          }
        }
      }

      data.set ( 20, 20, 20, 1.0 );
      data.set ( 20, 21, 20 , 1.0 );
      data.set ( 21, 20, 20, 1.0 );
      data.set ( 21, 21, 20, 1.0 );
      data.addToAll ( -0.5 );

      qc::SignedDistanceSweepingOp3D<double> dist;
      dist.applySingle ( data );
      // This test just checks whether using qc::SignedDistanceSweepingOp3D on this data
      // leads to an exception.
      testSuccess ( true, "" );
    }

    if(success) {
      aol::printSelfTestSuccessMessage ( "--                    EIKONAL Self Test Successful                            --" );
      aol::callSystemPauseIfNecessaryOnPlatform();
      return 0;
    }

  } catch(std::exception &ex){
    cerr << "\n\nstd::exception caught:\n";
    cerr << ex.what () << endl;
  } catch(aol::Exception &ex){
    cerr << "\n\naol::Exception caught:\n";
    ex.dump ();
  } catch (...){
    cerr << "\n\nUnknown exception caught.\n";
  }
  aol::printSelfTestFailureMessage ( "!!                    EIKONAL Self Test FAILED                                !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
