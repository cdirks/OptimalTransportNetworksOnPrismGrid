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

#include <configurators.h>
// -----------------------------------------------------------------------------
// smoothVectorField.cpp
// function smoothes a loaded vector-field, saves it and shows it in GRAPE
// -----------------------------------------------------------------------------

#include <quoc.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <preconditioner.h>

#include <mcm.h>
#include <fastUniformGridMatrix.h>

#include <anisotropies.h>
#include <Willmore.h>


#include <qmException.h>
#include <aol.h>
#include <smallMat.h>
#include <parameterParser.h>
#include "grapeInterface3d.h"


typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;


// ------------------- adjustIndex -----------------------------------------
// tests, whether an index i is contained in the regular area ([0,N])
// otherwise it is mirrored into the regular area
inline int adjustInd(int i, int N)
{
  if (i >= 0 && i<N) return i;
  else {
    if (i<0) return -i;
    else return ( 2 * (N-1) - i );
  }
}

// ------------ one method for all three indices -----------------------------
inline void adjustIndices(int &iAdj, int &jAdj, int &kAdj, int i, int j, int k,
                      int Nx, int Ny, int Nz)
{
  iAdj = adjustInd(i, Nx);
  jAdj = adjustInd(j, Ny);
  kAdj = adjustInd(k, Nz);
}





// ----------------------------------------------------------------------------
// method smoothes the picture by solving the heat equation
// ----------------------------------------------------------------------------

void smoothHeatEquation(qc::ScalarArray<double, qc::QC_3D> &img, int N, int timeSteps, double tau)
{

    int d = qc::logBaseTwo (N);
    qc::GridDefinition grid( d, qc::QC_3D );

    tau *= 0.5 * grid.H();

    // Operators
    aol::StiffOp< ConfigType > L( grid );
    aol::LumpedMassOp< ConfigType > M( grid, aol::DO_NOT_INVERT );  // false: do not inverse matrix

    qc::FastUniformGridMatrix<double, qc::QC_3D> mat( grid );
    aol::SSORPreconditioner<aol::Vector<double>, qc::FastUniformGridMatrix<double,qc::QC_3D> > precond( mat );
    aol::PCGInverse<aol::Vector<double> > solver( mat, precond, 1e-16, 1000 );

    qc::ScalarArray<double, qc::QC_3D> rhs( N,N,N );

    cerr << "assembling and applying (M + Tau*lambda*L).";
    mat.setZero();
    L.assembleAddMatrix( mat );  //
    cerr << ".";
    mat *= tau;
    M.assembleAddMatrix( mat );      // this adds to the already assembled matrix
    cerr << ".done, finishing rhs ...";

    // ---------------------------------------------------------------------

    for ( int iter=0; iter<timeSteps; ++iter ) {

      M.apply( img, rhs );
      cerr << "done.\nSolving ...";

      solver.apply( rhs, img );
      cerr << "done!\n";

    }
}


// ----------------------------------------------------------------------------
// method smoothes the picture, with a gaussian mask
// ----------------------------------------------------------------------------

void smoothGaussMask(qc::ScalarArray<double, qc::QC_3D> &feldX,
                 qc::ScalarArray<double, qc::QC_3D> &feldY,
                 qc::ScalarArray<double, qc::QC_3D> &feldZ,
                 const qc::ScalarArray<double, qc::QC_3D> &feldX_old,
                 const qc::ScalarArray<double, qc::QC_3D> &feldY_old,
                 const qc::ScalarArray<double, qc::QC_3D> &feldZ_old, int N)
{
  int i_pre, i_post, j_pre, j_post, k_pre, k_post;
  double val;
  // the weight for the actual element
  double weight = 6.;
  double divisor = 6. + weight;

  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        // get the indices of the circumjacent points
        i_pre  = adjustInd(i-1, N);
        i_post = adjustInd(i+1, N);
        j_pre  = adjustInd(j-1, N);
        j_post = adjustInd(j+1, N);
        k_pre  = adjustInd(k-1, N);
        k_post = adjustInd(k+1, N);

        // now evaluate the mean value of each component
        val = feldX_old.get(i_pre, j,k) + feldX_old.get(i_post, j,k)
              + feldX_old.get(i, j_pre, k) + feldX_old.get(i, j_post, k)
              + feldX_old.get(i,j, k_pre)  + feldX_old.get(i,j, k_post)
              + weight * feldX_old.get( i,j,k );
        feldX.set( i,j,k, val / divisor );

        val = feldY_old.get(i_pre, j,k) + feldY_old.get(i_post, j,k)
              + feldY_old.get(i, j_pre, k) + feldY_old.get(i, j_post, k)
              + feldY_old.get(i,j, k_pre)  + feldY_old.get(i,j, k_post)
              + weight * feldY_old.get( i,j,k );
        feldY.set( i,j,k, val / divisor );

        val = feldZ_old.get(i_pre, j,k) + feldZ_old.get(i_post, j,k)
              + feldZ_old.get(i, j_pre, k) + feldZ_old.get(i, j_post, k)
              + feldZ_old.get(i,j, k_pre)  + feldZ_old.get(i,j, k_post)
              + weight * feldZ_old.get( i,j,k );
        feldZ.set( i,j,k, val / divisor );

      }
}



// ----------------------------------------------------------------------------
// method smoothes the picture with a median filter
// ----------------------------------------------------------------------------


void arrangeVals(double* &values, int length)
{
  // use a simple bubble-sort (the array is short, so I think it's ok ;-)
  for (int i=0; i<length-1; i++) {
    for (int j=0; j<length-1-i; j++)
      if (values[j+1] < values[j]) {    /* compare the two neighbors           */
        double tmp = values[j];         /* swap values[j] and values[j+1]      */
        values[j] = values[j+1];
        values[j+1] = tmp;
    }
  }
}

void smoothMedianComponent(qc::ScalarArray<double, qc::QC_3D> &feld,
                           const qc::ScalarArray<double, qc::QC_3D> &feld_old, int N)
{
  // width of the neighbourhood
  int width = 1;
  int length = (2*width + 1) * (2*width + 1) * (2*width + 1);
  double *values = new double[length];
  int ind_i, ind_j, ind_k;    // indices

  // walk through the whole field
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      for (int k=0; k<N; ++k)
      {
        // now get the values from the neighbours and of the actual point (i,j,k) (9 elements)
        int count = 0;
        for (int is = i-width; is <= i+width; ++is)
          for (int js = j-width; js <= j+width; ++js)
            for (int ks = k-width; ks <= k+width; ++ks)
            {
              // get validated indices
              adjustIndices(ind_i, ind_j, ind_k, is, js, ks, N, N, N);
              values[count] = feld_old.get( ind_i, ind_j, ind_k );
              count ++;
            }

        if (count != length) cerr<<"ACHTUNG: count != length\n";
        // now arrange the values by size
        arrangeVals(values, length);
        // and take the average element
        int ind = static_cast<int>( length / 2.) + 1;
        feld.set( i,j,k, values[ind] );
      }
}


void smoothMedian(qc::ScalarArray<double, qc::QC_3D> &feldX,
                 qc::ScalarArray<double, qc::QC_3D> &feldY,
                 qc::ScalarArray<double, qc::QC_3D> &feldZ,
                 const qc::ScalarArray<double, qc::QC_3D> &feldX_old,
                 const qc::ScalarArray<double, qc::QC_3D> &feldY_old,
                 const qc::ScalarArray<double, qc::QC_3D> &feldZ_old, int N)
{
  // smooth each component
  smoothMedianComponent(feldX, feldX_old, N);
  smoothMedianComponent(feldY, feldY_old, N);
  smoothMedianComponent(feldZ, feldZ_old, N);
}





// this represents Tolga's method: First the Tensorproduct N*N^T is calculated,
// smoothed and then get the biggest eigenvector.
void smoothEVofTensorProduct( qc::ScalarArray<double, qc::QC_3D> *feld, int N,
                              int timesteps, double tau )
{
  // this will be the tensor-product
  qc::ScalarArray<double, qc::QC_3D> Tensor[3][3];
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      Tensor[i][j].reallocate ( N, N, N );


  cerr<<"\nCalculating tensor-product...";

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      // now traverse the whole scalararray in this place of the tensor-matrix
      for (int k=0; k<N; k++)
        for (int l=0; l<N; l++)
          for (int m=0; m<N; m++)
            (Tensor[i][j]).set(k,l,m, feld[i].get(k,l,m) * feld[j].get(k,l,m) );
    }


    cerr<<"smoothing the tensor-matrix...\n";
    // now smooth this rank 1 - tensor-matrix component-wise (smooth upper right matrix)

//     for (int i=0; i<3; i++)
//       for (int j=0; j<3; j++)
//       {
//         cerr<<"Smoothing component ["<<i<<"]["<<j<<"]...\n";
//         smoothHeatEquation( Tensor[i][j], N, timesteps, tau );
//       }
//

    for (int i=0; i<3; i++)
      for (int j=0; j<i+1; j++)
      {
        cerr<<"Smoothing component ["<<i<<"]["<<j<<"]...\n";
        smoothHeatEquation( Tensor[i][j], N, timesteps, tau );
      }
    // finally copy lower left matrix
    for (int k=0; k<N; k++)
      for (int l=0; l<N; l++)
        for (int m=0; m<N; m++) {
          (Tensor[0][1]).set(k,l,m, (Tensor[1][0]).get(k,l,m) );
          (Tensor[0][2]).set(k,l,m, (Tensor[2][0]).get(k,l,m) );
          (Tensor[1][2]).set(k,l,m, (Tensor[2][1]).get(k,l,m) );
        }

    cerr<<"getting the eigenvalues...";
    // ok, normally now we've got a rank 3 - matrix => get the biggest
    // eigenvector as the new direction
    for (int k=0; k<N; k++)
      for (int l=0; l<N; l++)
        for (int m=0; m<N; m++)
        {
          // now build a matrix
          aol::Matrix33Symm<double> smoothedMat(
            Tensor[0][0].get(k,l,m), Tensor[0][1].get(k,l,m), Tensor[0][2].get(k,l,m),
            Tensor[1][0].get(k,l,m), Tensor[1][1].get(k,l,m), Tensor[1][2].get(k,l,m),
            Tensor[2][0].get(k,l,m), Tensor[2][1].get(k,l,m), Tensor[2][2].get(k,l,m) );


          if (smoothedMat.infinityNorm() > 0.1)
          {
            // and get the eigenvectors
            aol::Vec3<double> eigenVals;
            aol::Matrix33<double> eigenVecs;
            smoothedMat.eigenVectors(eigenVals, eigenVecs);

            // sort the eigenvalues (mini-bubble) => ev[index[0]] is the biggest one
            int index[3] = { 0,1,2 };

            for (int i=0; i<2; i++) {
              if ( eigenVals[index[i]] < eigenVals[index[i+1]] ) {
                int temp = index[i]; index[i] = index[i+1]; index[i+1] = temp;
              }
            }
            if ( eigenVals[index[0]] < eigenVals[index[1]] ) {
              int temp = index[0]; index[0] = index[1]; index[1] = temp;
            }

            // save the eigenvector belonging to the biggest eigenvalue
            feld[0].set( k,l,m, eigenVecs[0][ index[0] ] );
            feld[1].set( k,l,m, eigenVecs[1][ index[0] ] );
            feld[2].set( k,l,m, eigenVecs[2][ index[0] ] );
          }
          else
          {
            feld[0].set( k,l,m, 0. );
            feld[1].set( k,l,m, 0. );
            feld[2].set( k,l,m, 0. );
          }

        }
}



// now the main program
int main( int argc, char **argv ) {

  // read the parameters with an parameter-file
  if ( argc != 2 ) {
    string s = "USAGE: ";
    s += argv[0];
    s += " <parameterfile>";
    throw aol::Exception( s.c_str(), __FILE__, __LINE__  );
  }

  try {
    aol::ParameterParser parser( argv[1] );

    // load image into scalar array
    char loadName[ 1024 ];
    char filename[ 1024 ];
    parser.getString( "loadName", loadName );

    cerr<<"Restoring Vektorfeld...";

    sprintf( filename, "%sX.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> feldX_old( filename );
    sprintf( filename, "%sY.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> feldY_old( filename );
    sprintf( filename, "%sZ.bz2", loadName);
    qc::ScalarArray<double, qc::QC_3D> feldZ_old( filename );

    int N = feldX_old.getNumX();
    qc::ScalarArray<double, qc::QC_3D> feldX( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> feldY( N,N,N );
    qc::ScalarArray<double, qc::QC_3D> feldZ( N,N,N );


    if (feldY_old.getNumX() != N || feldZ_old.getNumX() != N)
    {
      cerr<<"\n\nError: The vector-field doesn't have the same size as the image !!!!\n\n)";
      exit(43);
    }

    cerr<<"Restoring finished! \n\n";



    // -------------------------------------------------------------------------------
    // now the smoothing itself ......

    aol::StopWatch watch;
    watch.start();

    int steps = parser.getInt( "numberSteps" );

    for (int i=0; i < steps; i++) {

      cerr.precision(4);
      cerr<<"\rSmoothing vector-field, step "<<i+1<<" of "<<steps<<"      ";

      int method = parser.getInt( "method" );

      switch( method ) {
        // 0 = Median
        case 0 : smoothMedian(feldX, feldY, feldZ, feldX_old, feldY_old, feldZ_old, N);
                 break;

        // 1 = Gaussian Mask
        case 1 : smoothGaussMask(feldX, feldY, feldZ, feldX_old, feldY_old, feldZ_old, N);
                 break;

        // 2 = solving the heat equation
        case 2 : {
                  feldX = feldX_old;
                  feldY = feldY_old;
                  feldZ = feldZ_old;
                  double tau = aol::Sqr( parser.getDouble("tau") );
                  int timeSteps = parser.getInt( "timesteps");
                  cerr << "\nsmoothing x-component...\n";
                  smoothHeatEquation( feldX, N, timeSteps, tau );
                  cerr << "\nsmoothing y-component...\n";
                  smoothHeatEquation( feldY, N, timeSteps, tau );
                  cerr << "\nsmoothing z-component...\n";
                  smoothHeatEquation( feldZ, N, timeSteps, tau );
                 }
                 break;


        // 3 = Tolga's method using the biggest eigenvector of the smoothed tensor-product
        case 3 : {
                  feldX = feldX_old;
                  feldY = feldY_old;
                  feldZ = feldZ_old;
                  cerr<<"vor tau...";
                  double tau = aol::Sqr( parser.getDouble("tau") );
                  cerr<<"vor timesteps...";
                  int timeSteps = parser.getInt( "timesteps");
                  qc::ScalarArray<double, qc::QC_3D> feld[3];
                  feld[0].reallocate ( feldX ); feld[0] = feldX;
                  feld[1].reallocate ( feldY ); feld[1] = feldY;
                  feld[2].reallocate ( feldZ ); feld[2] = feldZ;

                  smoothEVofTensorProduct( feld, N, timeSteps, tau );
                  feldX = feld[0]; feldY = feld[1]; feldZ = feld[2];
                 }
                 break;

      };

    }
    cerr<<"\rSmoothing vector-field...100 %          \n";

    watch.stop();
    cerr<<"Finished! Took "<<watch.elapsedCpuTime()<<" seconds!\nATTENTION: NOT REGISTERED YET!\n";

    // -------------------------------------------------------------------------------



    // now save the smoothed vector field
    cerr<<"Saving the smoothed vector field ... \n";
    char saveName[ 1024 ];
    parser.getString( "saveName", saveName );

    sprintf( filename, "%sX.bz2", saveName);
    feldX.save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sY.bz2", saveName);
    feldY.save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sZ.bz2", saveName);
    feldZ.save( filename, qc::PGM_DOUBLE_BINARY );

    cerr<<"... sucessfully finished! \n\n";


    // finally view the result in GRAPE

    aol::MultiVector<double> smooth (0, feldX.getNumX() * feldY.getNumY() * feldZ.getNumZ());
    smooth.appendReference (feldX);
    smooth.appendReference (feldY);
    smooth.appendReference (feldZ);
    // view also the original data to compare them
    aol::MultiVector<double> original (0, feldX_old.getNumX() * feldY_old.getNumY()
                                          * feldZ_old.getNumZ() );
    original.appendReference (feldX_old);
    original.appendReference (feldY_old);
    original.appendReference (feldZ_old);
#ifdef USE_EXTERNAL_GRAPE
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&smooth, "smoothed VF");
    addVectorData(mesh, &original, "Original VF");
    // and then start GRAPE, thats it!
    initStartGrape(mesh, "smoothed VF");
#else
    cerr << "cannot display in GRAPE without grape external" << endl;
#endif
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

  return EXIT_SUCCESS;

}

