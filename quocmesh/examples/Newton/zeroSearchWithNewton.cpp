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
 * \brief Searches root of the system $F(x,y) = (x + \cos y, y - \sin x)$.
 *
 * Searches root of the system $F(x,y) = (x + \cos y, y - \sin x)$.
 * The configuration is read from a parameter file.
 * Note that NewtonIterationBase cannot be endowed with a simple Op as solver,
 * but requires an IterativeInverseOp.
 *
 * \author von Deylen
 */


#include <vec.h>
#include <matrix.h>
#include <Newton.h>

#include <cmath>

using namespace std;

typedef long double                 RealType;
typedef aol::Vector<RealType>       VecType;
typedef aol::FullMatrix<RealType>   MatType;

class F : public aol::Op<VecType, VecType> {
public:
  void applyAdd ( const VecType & x, VecType & Fx ) const {
    Fx[0] += x[0] + cos ( x[1] );
    Fx[1] += x[1] - sin ( x[0] );
  }
};

class DF : public aol::Op<VecType, MatType> {
public:
  void applyAdd ( const VecType & x, MatType & DFx ) const {
    DFx.set ( 0, 0,     1.            );
    DFx.set ( 0, 1,     -sin ( x[1] ) );
    DFx.set ( 1, 0,     -cos ( x[0] ) );
    DFx.set ( 1, 1,     1.            );
  }
};

class Newton2d : public aol::NewtonIterationBase<RealType, VecType, VecType, MatType> {
public:
  Newton2d ( MatType * pMatDF,
             const aol::Op<VecType, VecType> & f,
             const aol::Op<VecType, MatType> & df,
             aol::NewtonInfo<RealType> & newtonInfo )
    : aol::NewtonIterationBase<RealType, VecType, VecType, MatType> ( pMatDF, f, df, newtonInfo )
    , _pMatDFT ( NULL )
  {}

  ~Newton2d () {
    delete _pMatDFT;
  }

  virtual void prepareSolver () const {
    if (!_pMatDFT)
      _pMatDFT = new MatType ( *_pMatDF );
    _pMatDF->transposeTo ( *_pMatDFT );

    delete _pSolver;
    _pSolver = new aol::BiCGInverse<VecType> ( *_pMatDF, *_pMatDFT, this->_pInfo->getSolverInfo() );
  }
  mutable MatType * _pMatDFT;
};

int main( int, char** ) {

  try {
    F f;
    DF df;
    MatType * dfMatrix = new MatType ( 2, 2 );

    aol::NewtonInfo<RealType> * newtonInfo
      = aol::NewtonInfo<RealType>::createFromParameterFile ( "NewtonInfo.par" );
    Newton2d newton ( dfMatrix, f, df, *newtonInfo );

    VecType x0 ( 2 ), x ( 2 );
    x0[0] = 1000.;
    x0[1] = 1000.;

    newton.apply ( x0, x );
    cout << x << endl << endl;
  }
  catch(std::exception &ex){
    cerr << aol::color::error << endl;
    cerr << "\n\nstd::exception caught:\n";
    cerr << ex.what () << endl;
  }
  catch(aol::Exception &ex){
    cerr << aol::color::error << endl;
    cerr << "\n\naol::Exception caught:\n";
    ex.dump ();
  }
  catch (...){
    cerr << aol::color::error << endl;
    cerr << "\n\nUnknown exception caught.\n";
  }
  cerr << aol::color::reset;
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
