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
#include <qmException.h>
#include <aol.h>
#include <quoc.h>
#include <parameterParser.h>
#include <smallMat.h>

//#include "grapeInterface3d.h"
using namespace std;


/**************************************************************
 * class for calculating the mass, center of gravity and
 * the moments of inertia (Massenträgheitsmoment) of given
 * data in 3D
 * @author Nemitz
 */


template <typename DataType>
class moments {
 public:
  qc::ScalarArray<DataType, qc::QC_3D> _image, _mass;
  qc::ScalarArray<DataType, qc::QC_3D> *_centerOfMass[3];        // Schwerpunkte
  qc::ScalarArray<DataType, qc::QC_3D> *_massCentroidAxis[3];    // Hauptträgheitsachsen
  qc::ScalarArray<DataType, qc::QC_3D> *_moments[3][3];
  DataType _delta;                  // avoiding singularities
  DataType _epsilon;                // filter-width
  DataType _radiusFactor;           // for multiplying with the radii of the structures
  aol::StopWatch _watch;

 public:
  // constructor
  moments( const char *imgname, DataType delta, DataType epsilon )
    : _image(imgname), _mass(_image.getNumX(), _image.getNumY(), _image.getNumZ() )
  {
    _delta = delta;
    _epsilon = epsilon;

    for (int i=0; i<3; i++)
      _centerOfMass[i] = new qc::ScalarArray<DataType, qc::QC_3D> (_image.getNumX(), _image.getNumY(),
                                                          _image.getNumZ());

    for (int i=0; i<3; i++)
      _massCentroidAxis[i] = new qc::ScalarArray<DataType, qc::QC_3D> (_image.getNumX(), _image.getNumY(),
                                                          _image.getNumZ());

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        _moments[i][j] = new qc::ScalarArray<DataType, qc::QC_3D> (_image.getNumX(), _image.getNumY(),
                                                          _image.getNumZ());
  }

  ~moments( ) {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        delete _moments[i][j];

    for (int i=0; i<3; i++)
      delete _massCentroidAxis[i];

    for (int i=0; i<3; i++)
      delete _centerOfMass[i];
  }

  // ---------------------------------- methods -------------------------------------------------

  void setRadiusFactor( DataType f ) { _radiusFactor = f; }


  // auxiliary method
  DataType sqr(DataType x) { return(x*x); }

  // load already calculated moments
  void loadMoments(const char *loadName)
  {
    cerr<<"Restoring moments ... \n";
    char filename[ 1024 ];

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++) {
        sprintf( filename, "%s_%03d%03d.bz2", loadName, i,j );
        _moments[i][j] = new qc::ScalarArray<DataType, qc::QC_3D> ( filename );
      }
    cerr<<"Restoring finished! \n\n";
  }

  // load already calculated moments
  void saveMassCentroidAxis(const char *saveName)
  {
    cerr<<"Saving the Mass Centroid Axis ... \n";
    char filename[ 1024 ];

    sprintf( filename, "%sX.bz2", saveName);
    _massCentroidAxis[0]->save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sY.bz2", saveName);
    _massCentroidAxis[1]->save( filename, qc::PGM_DOUBLE_BINARY );
    sprintf( filename, "%sZ.bz2", saveName);
    _massCentroidAxis[2]->save( filename, qc::PGM_DOUBLE_BINARY );

    cerr<<"... sucessfully finished! \n\n";
  }

  // ------------------- adjustIndex -----------------------------------------
  // tests, whether an index i is contained in the regular area ([0,N])
  // otherwise it is mirrored into the regular area
  inline int adjustInd(int i, int N)
  {
    if (i >= 0 && i<N) return i;
    else {
      if (i<0) return -i;
      else return ( 2 * (N-1) - i );
//       if (i<0) return 0;
//       else return ( N-1 );
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


  // -------------------------------- MASS -----------------------------------------

  // calc the mass-matrix (for the whole picture)
  void calcMass(void)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();

    double prozent = 0;
    _watch.start();
    cerr<<endl;

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Masse..."<<prozent<<" %            ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          _mass.set(x,y,z, 0);

          // now loops over the epsilon-environment
          for (int i = x-radius; i < x+radius; ++i)
            for (int j = y-radius; j < y+radius; ++j)
              for (int k = z-radius; k < z+radius; ++k)
                if ( (sqr(x-i) + sqr(y-j) + sqr(z-k)) < sqr(_epsilon) )
                  _mass.add( x,y,z, _image.get( adjustInd(i,Nx), adjustInd(j,Ny), adjustInd(k,Nz) ) );
        }
      }
    }
    _watch.stop();
    cerr<<"\rBerechnung der Masse abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n";
  }


  // calc the mass-matrix (for the whole picture) with radii given in a scalar-array
  void calcMassWithAdjustedRadii(qc::ScalarArray<DataType, qc::QC_3D> &Radii)
  {
    double radius;
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();

    double prozent = 0;
    _watch.start();
    cerr<<endl;

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Masse mit angepassten Radien..."<<prozent<<" %            ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          _mass.set(x,y,z, 0);

          // only if there is a radius != 0
          radius = _radiusFactor * Radii.get(x,y,z);
          int intRadius = static_cast<int>(radius);
          if (intRadius != 0) {
            // now loops over the epsilon-environment
            for (int i = x-intRadius; i < x+intRadius; ++i)
              for (int j = y-intRadius; j < y+intRadius; ++j)
                for (int k = z-intRadius; k < z+intRadius; ++k)
                  if ( (sqr(x-i) + sqr(y-j) + sqr(z-k)) < sqr(radius) )
                    _mass.add( x,y,z, _image.get( adjustInd(i,Nx), adjustInd(j,Ny), adjustInd(k,Nz) ) );
          }
        }
      }
    }
    _watch.stop();
    cerr<<"\rAngepasste Berechnung der Masse abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n";
  }


  // calc the mass only in one point
  DataType calcMass(int xk, int yk, int zk)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();

    DataType res = 0;

    // now loops over the epsilon-environment
    for (int i = xk-radius; i < xk+radius; ++i)
      for (int j = yk-radius; j < yk+radius; ++j)
        for (int k = zk-radius; k < zk+radius; ++k)
          if ( (sqr(xk-i) + sqr(yk-j) + sqr(zk-k)) < sqr(_epsilon) )
            res += _image.get( adjustInd(i,Nx), adjustInd(j,Ny), adjustInd(k,Nz) );

    return res;
  }


  // ---------------------- CENTER OF MASS -----------------------------

  // calc the centers of gravity for the whole picture
  void calcCenterOfMass(void)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType nenner;
    DataType imDat1, imDat2, imDat3;

    double prozent = 0;
    _watch.start();

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Schwerpunkte..."<<prozent<<" %                ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          _centerOfMass[0]->set(x,y,z, 0);
          _centerOfMass[1]->set(x,y,z, 0);
          _centerOfMass[2]->set(x,y,z, 0);

          nenner = _mass.get(x,y,z) + _delta;    // mass is non-negative

          // now loops over the epsilon-environment
          for (int i = x-radius; i < x+radius; ++i)
            for (int j = y-radius; j < y+radius; ++j)
              for (int k = z-radius; k < z+radius; ++k)
                if ( (sqr(x-i) + sqr(y-j) + sqr(z-k) ) < sqr(_epsilon) )
                {
                  adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
                  imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-x) / nenner;
                  imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-y) / nenner;
                  imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-z) / nenner;

                  _centerOfMass[0]->add( x,y,z, imDat1);
                  _centerOfMass[1]->add( x,y,z, imDat2);
                  _centerOfMass[2]->add( x,y,z, imDat3);
                }
        }
      }
    }
    _watch.stop();
    cerr<<"\rBerechnung der Schwerpunkte abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n";
  }

  // calc the centers of gravity for the whole picture with radii given in a scalar-array
  void calcCenterOfMassWithAdjustedRadii(qc::ScalarArray<DataType, qc::QC_3D> &Radii)
  {
    double radius;
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType nenner;
    DataType imDat1, imDat2, imDat3;

    double prozent = 0;
    _watch.start();

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Schwerpunkte mit angepassten Radien..."<<prozent<<" %                ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          _centerOfMass[0]->set(x,y,z, 0);
          _centerOfMass[1]->set(x,y,z, 0);
          _centerOfMass[2]->set(x,y,z, 0);

          // only if there is a radius != 0
          radius = _radiusFactor * Radii.get(x,y,z);
          int intRadius = static_cast<int>(radius);
          if (intRadius != 0) {

            nenner = _mass.get(x,y,z) + _delta;    // mass is non-negative

            // now loops over the epsilon-environment
            for (int i = x-intRadius; i < x+intRadius; ++i)
              for (int j = y-intRadius; j < y+intRadius; ++j)
                for (int k = z-intRadius; k < z+intRadius; ++k)
                  if ( (sqr(x-i) + sqr(y-j) + sqr(z-k) ) < sqr(radius) )
                  {
                    adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
                    imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-x) / nenner;
                    imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-y) / nenner;
                    imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-z) / nenner;

                    _centerOfMass[0]->add( x,y,z, imDat1);
                    _centerOfMass[1]->add( x,y,z, imDat2);
                    _centerOfMass[2]->add( x,y,z, imDat3);
                  }
          }
        }
      }
    }
    _watch.stop();
    cerr<<"\rAngepasste Berechnung der Schwerpunkte abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n";
  }


  // calc the center of gravity only in one point
  aol::Vec2<DataType> calcCenterOfMass(int xk, int yk, int zk)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType nenner;
    DataType imDat1, imDat2, imDat3;

    aol::Vec3<DataType> res(0,0,0);

    nenner = calcMass(xk,yk,zk) + _delta;

    // now loops over the epsilon-environment
    for (int i = xk-radius; i < xk+radius; ++i)
      for (int j = yk-radius; j < yk+radius; ++j)
        for (int k = zk-radius; k < zk+radius; ++k)
          if ( (sqr(xk-i) + sqr(yk-j) + sqr(zk-k)) < sqr(_epsilon) )
          {
            adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
            imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-xk) / nenner;
            imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-yk) / nenner;
            imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-zk) / nenner;

            res[0] += imDat1;
            res[1] += imDat2;
            res[2] += imDat3;
          }

    return res;
  }


  // ------------------------- MOMENTS -----------------------------------

  // calc the moments for the whole picture and save them
  void calcMoments(const char *saveName)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType imDat1, imDat2, imDat3, nenner;

    double prozent = 0;
    _watch.start();

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Momente..."<<prozent<<" %      ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          // clear the 3x3-tensor-matrix
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                _moments[i][j]->set(x,y,z, 0);

          nenner = _mass.get(x,y,z) + _delta;

          // now loops over the epsilon-environment
          for (int i = x-radius; i < x+radius; ++i)
            for (int j = y-radius; j < y+radius; ++j)
              for (int k = z-radius; k < z+radius; ++k)
                if ( (sqr(x-i) + sqr(y-j) + sqr(z-k) ) < sqr(_epsilon) )
                {
                  adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
                  imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-x-(_centerOfMass[0]->get(x,y,z)));
                  imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-y-(_centerOfMass[1]->get(x,y,z)));
                  imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-z-(_centerOfMass[2]->get(x,y,z)));

                  // tensor-product
                  _moments[0][0]->add( x,y,z, imDat1*imDat1/nenner);
                  _moments[0][1]->add( x,y,z, imDat1*imDat2/nenner);
                  _moments[0][2]->add( x,y,z, imDat1*imDat3/nenner);
                  _moments[1][0]->add( x,y,z, imDat2*imDat1/nenner);
                  _moments[1][1]->add( x,y,z, imDat2*imDat2/nenner);
                  _moments[1][2]->add( x,y,z, imDat2*imDat3/nenner);
                  _moments[2][0]->add( x,y,z, imDat3*imDat1/nenner);
                  _moments[2][1]->add( x,y,z, imDat3*imDat2/nenner);
                  _moments[2][2]->add( x,y,z, imDat3*imDat3/nenner);
                }
        }
      }
    }
    _watch.stop();
    cerr<<"\rBerechnung der Momente abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n\n";
    if ( *saveName )  {
      cerr<<"Saving moments ... \n";
      char filename[ 1024 ];
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)  {
          sprintf( filename, "%s_%03d%03d.bz2", saveName, i,j );
          _moments[i][j]->save( filename, qc::PGM_DOUBLE_BINARY );
        }
    }
  }


  // calc the moments for the whole picture with radii given in a scalar-array and save them
  void calcMomentsWithAdjustedRadii(const char *saveName, qc::ScalarArray<DataType, qc::QC_3D> &Radii)
  {
    double radius;
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType imDat1, imDat2, imDat3, nenner;

    double prozent = 0;
    _watch.start();

    for (int x=0; x<Nx; ++x)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Momente mit angepassten Radien..."<<prozent<<" %      ";

      for (int y=0; y<Ny; ++y)
      {
        for (int z=0; z<Nz; ++z)
        {
          // clear the 3x3-tensor-matrix
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                _moments[i][j]->set(x,y,z, 0);

          // only if there is a radius != 0
          radius = _radiusFactor * Radii.get(x,y,z);
          int intRadius = static_cast<int>(radius);
          if (intRadius != 0) {

            nenner = _mass.get(x,y,z) + _delta;

            // now loops over the epsilon-environment
            for (int i = x-intRadius; i < x+intRadius; ++i)
              for (int j = y-intRadius; j < y+intRadius; ++j)
                for (int k = z-intRadius; k < z+intRadius; ++k)
                  if ( (sqr(x-i) + sqr(y-j) + sqr(z-k) ) < sqr(radius) )
                  {
                    adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
                    imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-x-(_centerOfMass[0]->get(x,y,z)));
                    imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-y-(_centerOfMass[1]->get(x,y,z)));
                    imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-z-(_centerOfMass[2]->get(x,y,z)));

                    // tensor-product
                    _moments[0][0]->add( x,y,z, imDat1*imDat1/nenner);
                    _moments[0][1]->add( x,y,z, imDat1*imDat2/nenner);
                    _moments[0][2]->add( x,y,z, imDat1*imDat3/nenner);
                    _moments[1][0]->add( x,y,z, imDat2*imDat1/nenner);
                    _moments[1][1]->add( x,y,z, imDat2*imDat2/nenner);
                    _moments[1][2]->add( x,y,z, imDat2*imDat3/nenner);
                    _moments[2][0]->add( x,y,z, imDat3*imDat1/nenner);
                    _moments[2][1]->add( x,y,z, imDat3*imDat2/nenner);
                    _moments[2][2]->add( x,y,z, imDat3*imDat3/nenner);
                  }
          }
        }
      }
    }
    _watch.stop();
    cerr<<"\rAngepasste Berechnung der Momente abgeschlossen! ("<<_watch.elapsedCpuTime()<<" sec)\n\n";
    if ( *saveName )  {
      cerr<<"Saving moments ... \n";
      char filename[ 1024 ];
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)  {
          sprintf( filename, "%s_%03d%03d.bz2", saveName, i,j );
          _moments[i][j]->save( filename, qc::PGM_DOUBLE_BINARY );
        }
    }
  }


  // calc the moment-matrix only in one point
  aol::Matrix33<DataType> calcMoments(int xk, int yk, int zk, aol::Vec3<DataType> &centerOfMass)
  {
    int radius = static_cast<int>(_epsilon);
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();
    int iAdj, jAdj, kAdj;
    DataType imDat1, imDat2, imDat3, nenner;

    aol::Matrix33<DataType> res(0,0,0, 0,0,0, 0,0,0);

    nenner = calcMass(xk,yk,zk) + _delta;

    // now loops over the epsilon-environment
    for (int i = xk-radius; i < xk+radius; ++i)
      for (int j = yk-radius; j < yk+radius; ++j)
        for (int k = zk-radius; k < zk+radius; ++k)
          if ( ( sqr(xk-i) + sqr(yk-j) +sqr(zk-k) ) < sqr(_epsilon) )
          {
            adjustIndices(iAdj, jAdj, kAdj, i,j,k, Nx, Ny, Nz);
            imDat1 = _image.get( iAdj, jAdj, kAdj ) * (i-xk-(centerOfMass[0]));
            imDat2 = _image.get( iAdj, jAdj, kAdj ) * (j-yk-(centerOfMass[1]));
            imDat3 = _image.get( iAdj, jAdj, kAdj ) * (k-zk-(centerOfMass[2]));

            // tensor-product
            res[0][0] += imDat1*imDat1/nenner;
            res[0][1] += imDat1*imDat2/nenner;
            res[0][2] += imDat1*imDat3/nenner;
            res[1][0] += imDat2*imDat1/nenner;
            res[1][1] += imDat2*imDat2/nenner;
            res[1][2] += imDat2*imDat3/nenner;
            res[2][0] += imDat3*imDat1/nenner;
            res[2][1] += imDat3*imDat2/nenner;
            res[2][2] += imDat3*imDat3/nenner;
          }

    return res;
  }




  // -------------------------------------------------------------------------------------
  // --------------- Calculations of the EIGENVECTOR of the biggest Eigenvalue -----------
  // ----------------- entspricht Hauptträgheitsachsen (mass centroid axis) --------------
  // -------------------------------------------------------------------------------------

  void calcMassCentroidAxis( DataType barrier )
  {
    int Nx = _image.getNumX();
    int Ny = _image.getNumY();
    int Nz = _image.getNumZ();

    double prozent = 0;
    _watch.start();

    for (int x=0; x<Nx; x++)
    {
      prozent = static_cast<double> ( x * 100. / Nx );
      cerr.precision(4);
      cerr<<"\rBerechne Hauptträgheitsachsen..."<<prozent<<" %              ";

      for (int y=0; y<Ny; y++)
        for (int z=0; z<Nz; z++)
        {
          // the tensor-matrix
          aol::Matrix33Symm<DataType> moment(
            _moments[0][0]->get(x,y,z), _moments[0][1]->get(x,y,z), _moments[0][2]->get(x,y,z),
            _moments[1][0]->get(x,y,z), _moments[1][1]->get(x,y,z), _moments[1][2]->get(x,y,z),
            _moments[2][0]->get(x,y,z), _moments[2][1]->get(x,y,z), _moments[2][2]->get(x,y,z) );

          //cerr<<"\nMatrix: \n"<<moment<<endl;

          aol::Vec3<DataType> eigenVals;
          aol::Matrix33<DataType> eigenVecs;
          moment.eigenVectors(eigenVals, eigenVecs);

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

          // now calc the metric for linearly anisotropy
          DataType sum = eigenVals[index[0]] + eigenVals[index[1]] + eigenVals[index[2]];
          if (sum != 0)
          {
            DataType cl = ( eigenVals[index[1]] + eigenVals[index[2]] ) / sum;
            // now test, whether the mass-structure is linear (cl < eps) or not (cl >= eps)

//             cerr<<"Eigenwerte: "<<eigenVals[0]<<", "<<eigenVals[1]<<", "<<eigenVals[2]<<" - grösste: ";
//             cerr<<eigenVals[index[0]]<<", kleinste: "<<eigenVals[index[2]]<<", cl: "<<cl<<endl;

            if (cl < barrier)    // that means ev2 = ev1 = 0 << ev0 (ca.) => mass like a stab
            {                    // now the direction of the ellipsoid is the belonging eigenvector
              _massCentroidAxis[0]->set( x,y,z, eigenVecs[0][index[0]]);
              _massCentroidAxis[1]->set( x,y,z, eigenVecs[1][index[0]]);
              _massCentroidAxis[2]->set( x,y,z, eigenVecs[2][index[0]]);
            } else {
              _massCentroidAxis[0]->set( x,y,z, 0.);
              _massCentroidAxis[1]->set( x,y,z, 0.);
              _massCentroidAxis[2]->set( x,y,z, 0.);
            }

          } else    // 0-matrix => set axis to 0
          {
            _massCentroidAxis[0]->set( x,y,z, 0.);
            _massCentroidAxis[1]->set( x,y,z, 0.);
            _massCentroidAxis[2]->set( x,y,z, 0.);
          }
        }    // end for z
    }  // end for x
    _watch.stop();
     cerr<<"\rBerechnung der Hauptträgheitsachsen abgeschlossen ("<<_watch.elapsedCpuTime()<<" sec)\n\n";
  }



  // ------------------------------ VISUALIZATION ----------------------------------------


  void showInGrape(void)
  {
#ifdef USE_EXTERNAL_GRAPE
    // convert this to a genmesh (it will automatically be tested
    // whether the data is quadratic or not)

    // the image itself and the masses
    GENMESH3D* mesh = quocmesh_convert_to_gmesh3d (&_image, "Bilddaten");
//     _mass /= _mass.getMaxValue();
//     addScalarData(mesh, &_image, "Daten");
//     addScalarData(mesh, &_mass, "Masse");
//
//     // center of mass
//     aol::MultiVector<double> SP (0, _image.getNumX() * _image.getNumY() * _image.getNumZ());
//     SP.append ((*_centerOfMass[0]));
//     SP.append ((*_centerOfMass[1]));
//     SP.append ((*_centerOfMass[2]));
//     addVectorData(mesh, &SP, "Schwerpunkte");

    // center of mass
    aol::MultiVector<double> TA (0, _image.getNumX() * _image.getNumY() * _image.getNumZ());
    TA.appendReference ((*_massCentroidAxis[0]));
    TA.appendReference ((*_massCentroidAxis[1]));
    TA.appendReference ((*_massCentroidAxis[2]));
    addVectorData(mesh, &TA, "Traegheitsachsen");

    // and then start GRAPE, thats it!
    initStartGrape(mesh, "Pfad");
#else
    cerr << "cannot display in GRAPE without grape external" << endl;
#endif
  }




  protected:

};
