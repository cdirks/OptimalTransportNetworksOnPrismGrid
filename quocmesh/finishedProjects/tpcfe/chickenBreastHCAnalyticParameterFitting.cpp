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

#include <vec.h>
#include <randomGenerator.h>
#include <scalarArray.h>


template< typename RealType >
class My1DOptimizer {
protected:
  const aol::Vector<RealType> &_measuredValues;
  const RealType _startBeta, _stopBeta;
  const int _numBeta, _numSteps;

public:
  My1DOptimizer ( const aol::Vector<RealType> &MeasuredValues, const RealType StartBeta, const RealType StopBeta, const int NumBeta, const int NumSteps )
    : _measuredValues ( MeasuredValues ),
      _startBeta ( StartBeta ),
      _stopBeta ( StopBeta ),
      _numBeta ( NumBeta ),
      _numSteps ( NumSteps ) {
  }

  virtual ~My1DOptimizer ( ) { }

  std::pair<RealType, RealType> optimize ( ) {
    RealType currentRange = _stopBeta - _startBeta, currentOpt = 0.5 * ( _startBeta + _stopBeta ), currentOptValue = evaluateFunction ( currentOpt );

    for ( int i = 0; i < _numSteps; ++i ) {
      aol::Vector<RealType> betaValues ( _numBeta ), fctValues ( _numBeta );
      for ( int j = 0; j < _numBeta; ++j ) {
        betaValues[j] = currentOpt - currentRange / 2 + ( 1.0 * j ) / ( _numBeta - 1) * currentRange;
        fctValues[j] = evaluateFunction ( betaValues[j] );
      }

      cerr << betaValues << endl << fctValues << endl;

      const std::pair < int, RealType > indMin = fctValues.getMinIndexAndValue();
      currentOpt = betaValues[ indMin.first ];
      currentOptValue = indMin.second;
      currentRange /= ( _numBeta - 1 );
    }

    return ( std::pair<RealType, RealType> ( currentOpt, currentOptValue ) );
  }

protected:
  virtual RealType evaluateFunction ( const RealType arg ) const = 0;
};


#if 0
class MyTestOptimizer : public My1DOptimizer<double> {
public:
  MyTestOptimizer ( const aol::Vector<double> &MeasuredValues, const double StartBeta, const double StopBeta, const int NumBeta, const int NumSteps )
    : My1DOptimizer<double> ( MeasuredValues, StartBeta, StopBeta, NumBeta, NumSteps ) {
  }

protected:
  double evaluateFunction ( const double arg ) const {
    double ret = 0.0;
    for ( int i = 0; i < _measuredValues.size(); ++i ) {
      ret += aol::Sqr ( arg * i * i - _measuredValues[i] );
    }
    return ( ret );
  }
};


int main ( int, char** ) {
  aol::Vector<double> myVec ( 100 );
  for ( int i = 0; i < 100; ++i )
    myVec[i] = 1.2453336897365247 * i * i;

  MyTestOptimizer mto ( myVec, 1.0, 2.0, 20, 10 );
  std::pair< double, double > minimum = mto.optimize();
  cerr << minimum.first << " " << minimum.second << endl;
}
#endif

class MyAHOptimizer : public My1DOptimizer<double> {
public:
  MyAHOptimizer ( const aol::Vector<double> &MeasuredValues, const double StartBeta, const double StopBeta, const int NumBeta, const int NumSteps )
    : My1DOptimizer<double> ( MeasuredValues, StartBeta, StopBeta, NumBeta, NumSteps ) {
  }

  double SmirnovValue ( const double beta, const double time ) const {
    static const double uZero = -34;
    return ( uZero * exp ( beta * beta * time ) * ( 1 - aol::Erf ( beta * sqrt( time ) ) ) );
  }

protected:
  // midpoint integration assuming measured temperature at center of time interval
  double evaluateFunction ( const double arg ) const {
    double ret = 0.0;
    const double
      timeStep = 120;
    for ( int i = 0; i < _measuredValues.size(); ++i ) {
      const double
        midTime = timeStep * ( i + 0.5 ),
        SmirnovVal = SmirnovValue ( arg, midTime );
      ret += timeStep * aol::Sqr ( SmirnovVal - _measuredValues[i] );
    }
    return ( ret );
  }
};


class MyHOptimizer : public My1DOptimizer<double> {
protected:
  double _betaKnown;

public:
  MyHOptimizer ( const aol::Vector<double> &MeasuredValues, const double StartBeta, const double StopBeta, const int NumBeta, const int NumSteps, const double BetaKnown )
    : My1DOptimizer<double> ( MeasuredValues, StartBeta, StopBeta, NumBeta, NumSteps ), _betaKnown ( BetaKnown ) {
  }

  double SmirnovValue ( const double h, const double time ) const {
    static const double uZero = -34, xZero = 0.05;
    return (  uZero * aol::Erf ( xZero / ( 2 * _betaKnown / h * sqrt ( time ) ) ) + uZero * exp ( _betaKnown * _betaKnown * time + h * xZero ) * ( 1 - aol::Erf ( xZero / ( 2 * _betaKnown / h * sqrt ( time ) ) + _betaKnown * sqrt ( time ) ) )  );
  }

protected:
  // midpoint integration assuming measured temperature at center of time interval
  double evaluateFunction ( const double arg ) const {
    double ret = 0.0;
    const double
      timeStep = 120;
    for ( int i = 0; i < _measuredValues.size(); ++i ) {
      const double
        midTime = timeStep * ( i + 0.5 ),
        SmirnovVal = SmirnovValue ( arg, midTime );
      ret += timeStep * aol::Sqr ( SmirnovVal - _measuredValues[i] );
    }
    return ( ret );
  }

};


int main ( int, char** ) {

  double optimizedAHValue = -1.0, optimizedHValue = -1.0;

  {
    aol::Vector<double> measuredu0 ( 117 ), simulatedu0 ( 117 );
    ifstream measuredu0loader ( "u0measured.dat" );
    measuredu0loader >> measuredu0;
    MyAHOptimizer ahOptimizer ( measuredu0, 0.0, 0.1, 21, 10 );
    std::pair< double, double > minimumu0 = ahOptimizer.optimize();
    cerr << minimumu0.first << " " << minimumu0.second << endl;

    for ( int i = 0; i < simulatedu0.size(); ++i ) {
      simulatedu0[i] = ahOptimizer.SmirnovValue ( minimumu0.first, 120 * i );
    }

    ofstream comparisonu0 ( "comparisonu0.dat" );
    for ( int i = 0; i < simulatedu0.size(); ++i ) {
      comparisonu0 << 120 * i << " " << measuredu0[i] << " " << simulatedu0[i] << endl;
    }

    optimizedAHValue = minimumu0.first;
  }

  {
    aol::Vector<double> measuredu5 ( 117 ), simulatedu5 ( 117 );
    ifstream measuredu5loader ( "u5measured.dat" );
    measuredu5loader >> measuredu5;
    MyHOptimizer hOptimizer ( measuredu5, 0.0, 100.0, 21, 10, optimizedAHValue );
    std::pair< double, double > minimumu5 = hOptimizer.optimize();
    cerr << minimumu5.first << " " << minimumu5.second << endl;

    for ( int i = 0; i < simulatedu5.size(); ++i ) {
      simulatedu5[i] = hOptimizer.SmirnovValue ( minimumu5.first, 120 * i );
    }

    ofstream comparisonu5 ( "comparisonu5.dat" );
    for ( int i = 0; i < simulatedu5.size(); ++i ) {
      comparisonu5 << 120 * i << " " << measuredu5[i] << " " << simulatedu5[i] << endl;
    }

    optimizedHValue = minimumu5.first;
  }

  cerr << "Parameter fitting yields:" << endl
       << "a * h = " << optimizedAHValue << endl
       << "h = " << optimizedHValue << endl;

  return ( EXIT_SUCCESS );
}

