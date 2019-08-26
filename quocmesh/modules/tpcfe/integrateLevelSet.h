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

#ifndef __INTEGRATELEVELSET_H
#define __INTEGRATELEVELSET_H

// ***************************************************************
// class for integrating a given function over a given level-set.
// works in 2d (based on Droskes Mumford2d) and in 3d (based on
// Liehrs CFE-Ops).
// The Level Set function has to be implemented in the derived
// class (put into practice with Barton-Nackman)
// ***************************************************************

#include <tpCFEGrid.h>
#include <tpCFEUtils.h>
#include <scalarArray.h>
#include <discreteFunction.h>
#include <gridBase.h>
#include <isolineIterator2d.h>
#include <configurators.h>


namespace qc {

template <typename RealType, Dimension Dim/*, typename Imp*/>
class IntegrateFunctionOverLevelSet {
};

template <typename RealType, Dimension Dim/*, typename Imp*/>
class IntegrateFunctionSqrOverLevelSet {
};

template <typename RealType, Dimension Dim/*, typename Imp*/>
class IntegrateTangentGradientOfFunctionSqrOverLevelSet {
};

// ------------------------ class for integrating in 2D ----------------------


// !!!!!!!!
// ******** ATTENTION: class not really complete yet
// IntegrateFunctionOverLevelSet has to be converted to Barton-Nackman
// Thereby someone is able to implement arbitrary function-evaluate methods
// in an easy way
// !!!!!!!!

template <typename RealType/*, typename Imp*/>
class IntegrateFunctionOverLevelSet<RealType, qc::QC_2D> {
protected:
  const qc::GridDefinition &_grid;
  const qc::ScalarArray<RealType, qc::QC_2D> &_ls;

  static RealType _b1   ( RealType X, RealType Y ) { return ( 1. - X ) * ( 1 - Y ); }
  static RealType _b2   ( RealType X, RealType Y ) { return X* ( 1. - Y ); }
  static RealType _b3   ( RealType X, RealType Y ) { return ( 1. - X ) *Y; }
  static RealType _b4   ( RealType X, RealType Y ) { return X*Y; }
  typedef RealType ( *BASIS_FUNC_TYPE ) ( RealType, RealType );
  BASIS_FUNC_TYPE basis[ 4 ];

  qc::FastILexMapper<qc::QC_2D> _mapper;

  typedef qc::QuocConfiguratorTraitMultiLin< RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;

public:
  IntegrateFunctionOverLevelSet ( const qc::GridDefinition &Grid,
                                  const qc::ScalarArray<RealType, qc::QC_2D> &LevelSetFunction )
      : _grid ( Grid ), _ls ( LevelSetFunction ), _mapper ( Grid ) {
    basis[ 0 ] = _b1;
    basis[ 1 ] = _b2;
    basis[ 2 ] = _b3;
    basis[ 3 ] = _b4;
  }

  RealType integrate ( const qc::Array<RealType> &Function ) {
    qc::IsoLineManager2d<ConfiguratorType> isoManager ( _grid, _ls );

    const RealType alpha[2] = { 0.21132487, 0.78867513 };

    RealType integral = 0.;

    for ( qc::IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ ) {

      RealType v[2], len;
      RealType w[2] = { isoIt->weights[0], isoIt->weights[1] };
      aol::Vec3<RealType> pt ( isoIt->el.x(), isoIt->el.y(), 0. );
      aol::Vec3<RealType> qpts[2];

      aol::Vec3<RealType> tmp;
      tmp =   isoIt->points[0];
      tmp -=  isoIt->points[1];
      len = tmp.norm();

      for ( int j = 0; j < 2; j++ ) {
        // compute speed at end-points
        v[j] = w[j] * Function.get ( isoIt->intersect[j][0] )
               + ( 1. - w[j] ) * Function.get ( isoIt->intersect[j][1] );

        // compute quadrature points on line
        qpts[j]  = isoIt->points[0];
        qpts[j] *= alpha[j] / ( 1. - alpha[j] );
        qpts[j] += isoIt->points[1];
        qpts[j] *= 1. - alpha[j];

        qpts[j] -= pt; // make reference coords
      }

      // compute line integrals here
      // quadrature
      RealType a = 0.;
      for ( int j = 0; j < 2; j++ ) {
        a += 0.5 * ( v[0] * alpha[j] + v[1] * ( 1. - alpha[j] ) );
      }
      a *= len * _grid.H();

      integral += a;
    }
    return integral;
  }

protected:

//   barton-nackman
//   inline Imp& asImp() { return static_cast<Imp&>(*this); }
//   inline const Imp& asImp() const { return static_cast<const Imp&>(*this); }
};


// ------------------------ class for integrating in 3D ----------------------


template <typename RealType, typename Imp>
class IntegrateFunctionOverLevelSetBase3d {
public:
  const qc::GridDefinition                    &_grid;
  const qc::ScalarArray<RealType, qc::QC_3D>  &_levelSetFunction;
  tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE>   _cfeGrid;


  IntegrateFunctionOverLevelSetBase3d ( const qc::GridDefinition &Grid,
                                        const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction )
      : _grid ( Grid ), _levelSetFunction ( LevelSetFunction ), _cfeGrid ( _grid.getGridDepth() ) {
    _cfeGrid.setAdjustLevelset ( 1.0e-10 );
    _cfeGrid.setDomainFrom ( LevelSetFunction );
    _cfeGrid.detectVirtualNodes ();
  }

  RealType integrate ( bool Verbose = 0 )  {
    RealType val = aol::ZTrait<RealType>::zero;
    const RealType h  = _cfeGrid.H();

    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      std::vector< aol::Vec3< RealType > > dest = itit.getTriangRef().getLocalCoordinates();
      aol::Vec3<RealType> centerOfMass;
      centerOfMass = ( dest[0] + dest[1] + dest[2] ) / ( static_cast<RealType> ( 3 ) );

      // midpoint quadrature
      const double area    = itit.getTriangRef().getArea ( h );
      const double valTemp = asImp().evaluateFunction ( itit.getTriangRef().getElement(), centerOfMass );

      if ( aol::isNaN ( valTemp ) )
        cout << "found NaN. Will continue without counting this triangle." << endl;
      else
        val += area * valTemp;
      if ( Verbose )
        cerr << "(" <<  valTemp << ", " << area << "), ";

    }

    return ( val );
  }


  RealType maxValueOnMidpoints () {
    RealType val = aol::ZTrait<RealType>::zero;

    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      std::vector< aol::Vec3< RealType > > dest = itit.getTriangRef().getLocalCoordinates();
      aol::Vec3<RealType> centerOfMass;
      centerOfMass = ( dest[0] + dest[1] + dest[2] ) / ( static_cast<RealType> ( 3 ) );
      RealType valTemp = asImp().evaluateFunction ( itit.getTriangRef().getElement(), centerOfMass );
      if ( valTemp > val )
        val = valTemp;
    }
    return ( val );
  }


  void integrateSaveEachValue ( qc::ScalarArray<RealType, qc::QC_3D> & values, bool Verbose = 0 )  {
    const RealType h  = _cfeGrid.H();

    for ( tpcfe::CFEInterfaceTriangleIterator< tpcfe::CFEGrid<RealType, tpcfe::CFE_NONE> > itit ( _cfeGrid ); itit.notAtEnd(); ++itit ) {
      std::vector< aol::Vec3< RealType > > dest = itit.getTriangRef().getLocalCoordinates();
      aol::Vec3<RealType> centerOfMass;
      centerOfMass = ( dest[0] + dest[1] + dest[2] ) / ( static_cast<RealType> ( 3 ) );

      // midpoint quadrature
      const double area    = itit.getTriangRef().getArea ( h );
      const double valTemp = asImp().evaluateFunction ( itit.getTriangRef().getElement(), centerOfMass );

      values.add ( itit.getElementRef().x(), itit.getElementRef().y(), itit.getElementRef().z(),  area * valTemp );
      if ( Verbose )
        cerr << "(" <<  valTemp << ", " << area << "), ";

    }
  }


  //! ---------- method to integrate only over an euclidian epsilon-ball ---------
  //! arg Function is the function to be integrated,
  //! params are the epicenter of the ball and its radius
  RealType integrateBall ( double epsilon, aol::Vec3<double> centre, bool Ausgabe = 0 ) /*const*/ {

    // the value of the integral
    RealType val = 0;
    const qc::GridSize<qc::QC_3D> gridSize = _cfeGrid.getSize();
    const double h  = _cfeGrid.H();

    qc::Ball epsBall ( _grid, epsilon );
    qc::Ball::iteratorBall ballIt;
    qc::Ball::iteratorBall ballItEnd ( epsBall.end ( centre ) );

    for ( ballIt = epsBall.begin ( centre ); ballIt != ballItEnd; ballIt++ ) {
      tpcfe::CFEElement<RealType> cfeEl ( gridSize, ballIt->x(), ballIt->y(), ballIt->z(), _cfeGrid.getElType ( ballIt->x(), ballIt->y(), ballIt->z() ) );
      if ( cfeEl.cfeType().representsInterfaced() )
        _cfeGrid.getCutRelations ( cfeEl );

      if ( _cfeGrid.isInterfaced ( cfeEl ) ) {
        for ( tpcfe::CFETopoTetraIterator tit ( cfeEl.cfeType(), -1 ); tit.notAtEnd(); ++tit ) {
          const tpcfe::CFETopoTetra &tet = *tit;
          if ( tet.getVirtualNodeNum() == 3 ) {                       // tetra with three virtual nodes has one interface-face

            std::vector< aol::Vec3<RealType> > dest;
            dest.reserve ( 3 );
            for ( int i = 0; i < 4; ++i ) {
              if ( tet.isVirtualNode ( i ) ) {
                aol::Vec3<RealType> tmp;
                tet.computeLocalCoordinate ( tmp, cfeEl, i );
                dest.push_back ( tmp );
              }
            }                                                         // now dest contains the global coordinates of the interface triangle

            // dest[i] are the coordinates of the triangle => calc its center of mass
            aol::Vec3<RealType> centerOfMass;
            centerOfMass = ( dest[0] + dest[1] + dest[2] ) / ( static_cast<RealType> ( 3 ) );

            // evaluate the function in the center of mass and add the value to the integral
            double valTemp = asImp().evaluateFunction ( cfeEl, centerOfMass );

            // area of the triangle
            // it's not necessary to compute world coords in the element because we only consider differences => *h
            for ( int i = 0; i < 3; ++i ) dest[i] *= h;
            aol::Vec3<RealType> Temp ( dest[2] - dest[0] );
            double area = 0.5 * ( Temp.crossProduct ( ( dest[2] - dest[1] ) ) ).norm() ;

            val += area * valTemp;
            if ( Ausgabe ) cerr << "(" <<  valTemp << ", " << area << "), ";
          }
        }
      }
    }

    return val;

  }



protected:
  //   barton-nackman
inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};



// ------------------------------------------------------------------------------------------



// ---------- derived class for functions given by a scalarArray3d ------------

template <typename RealType>
class IntegrateFunctionOverLevelSet<RealType, qc::QC_3D> : public qc::IntegrateFunctionOverLevelSetBase3d
      <RealType, IntegrateFunctionOverLevelSet<RealType, qc::QC_3D> > {

public:
  // TODO: only temporarely
  typedef qc::QuocConfiguratorTraitMultiLin < RealType, qc::QC_3D,
  aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfigType;

protected:
  aol::DiscreteFunctionDefault<ConfigType> _discFunc;

public:
  // ---------------------- constructor --------------------------------------------
  IntegrateFunctionOverLevelSet ( const qc::GridDefinition &Grid,
                                  const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction,
                                  const qc::ScalarArray<RealType, qc::QC_3D> &Function  ) :
      qc::IntegrateFunctionOverLevelSetBase3d < RealType,
      IntegrateFunctionOverLevelSet<RealType, qc::QC_3D> > ( Grid, LevelSetFunction ),
      _discFunc ( Grid, Function ) {}

  // ----------------------------methods --------------------------------------------
  // evaluate the function by using the level set data

  RealType evaluateFunction ( const qc::Element &El,
                              const aol::Vec3<RealType> &refCoords )  const {
    return _discFunc.evaluate ( El, refCoords );
  }

};

// ---------- derived class for calculating the square -----------------------
// ---------- of a function over the given level set. ------------------------

//! ATTENTION: Be aware, that this method delivers other results than just integrating
//! over the function with square nodal values, because here you first interpolate and
//! then compute the square, whereas otherwise you interpolate the squares.

template <typename RealType>
class IntegrateFunctionSqrOverLevelSet<RealType, qc::QC_3D> :
      public qc::IntegrateFunctionOverLevelSetBase3d <RealType, IntegrateFunctionSqrOverLevelSet<RealType, qc::QC_3D> > {

public:
  // TODO: only temporarely
  typedef qc::QuocConfiguratorTraitMultiLin < RealType, qc::QC_3D,
  aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfigType;

protected:
  aol::DiscreteFunctionDefault<ConfigType> _discFunc;

public:
  // ---------------------- constructor --------------------------------------------
  IntegrateFunctionSqrOverLevelSet ( const qc::GridDefinition &Grid,
                                     const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction,
                                     const qc::ScalarArray<RealType, qc::QC_3D> &Function  ) :
      qc::IntegrateFunctionOverLevelSetBase3d < RealType,
      IntegrateFunctionSqrOverLevelSet<RealType, qc::QC_3D> > ( Grid, LevelSetFunction ),
      _discFunc ( Grid, Function ) {}

  // ----------------------------methods --------------------------------------------
  // evaluate the function by using the level set data

  RealType evaluateFunction ( const qc::Element &El,
                              const aol::Vec3<RealType> &refCoords )  const {
    return aol::Sqr ( _discFunc.evaluate ( El, refCoords ) );
  }

};


// ---------- derived class for calculating the square ---------------------------------------
// ---------- of the gradient of a function over the given level set (for H^1-norms). --------

template <typename RealType>
class IntegrateTangentGradientOfFunctionSqrOverLevelSet<RealType, qc::QC_3D> :
      public qc::IntegrateFunctionOverLevelSetBase3d <RealType, IntegrateTangentGradientOfFunctionSqrOverLevelSet<RealType, qc::QC_3D> > {

public:
  // TODO: only temporarely
  typedef qc::QuocConfiguratorTraitMultiLin < RealType, qc::QC_3D,
  aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfigType;

protected:
  aol::DiscreteFunctionDefault<ConfigType> _discFunc;

public:
  // ---------------------- constructor --------------------------------------------
  IntegrateTangentGradientOfFunctionSqrOverLevelSet ( const qc::GridDefinition &Grid,
                                                      const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction,
                                                      const qc::ScalarArray<RealType, qc::QC_3D> &Function  ) :
      qc::IntegrateFunctionOverLevelSetBase3d < RealType,
      IntegrateTangentGradientOfFunctionSqrOverLevelSet<RealType, qc::QC_3D> > ( Grid, LevelSetFunction ),
      _discFunc ( Grid, Function ) {}

  // ----------------------------methods --------------------------------------------
  // evaluate the function by using the level set data

  RealType evaluateFunction ( const qc::Element &El,
                              const aol::Vec3<RealType> &refCoords )  const {
    aol::Vec3<RealType> gradPhi,                // gradient of the level-set-function
    gradient,               // euclidian gradient of the function
    tangentialGradient;     // tangential gradient of the function
    _discFunc.evaluateGradient ( El, refCoords, gradient );

    aol::DiscreteFunctionDefault<ConfigType> _levelSetDiscrFunc ( this->_grid, this->_levelSetFunction );
    _levelSetDiscrFunc.evaluateGradient ( El, refCoords, gradPhi );

    // compute the projection matrix to project the gradient to the tangent space
    aol::Matrix33<RealType> PMat;
    for ( int i = 0; i < ConfigType::VecType::dim; i++ )
      for ( int j = 0; j < ConfigType::VecType::dim; j++ )
        PMat[i][j] = -gradPhi[i] * gradPhi[j];

    for ( int i = 0; i < ConfigType::VecType::dim; i++ )
      PMat[i][i] += 1.;

    PMat.mult ( gradient, tangentialGradient );

    return tangentialGradient.normSqr();
  }

};


// ---------- derived class for calculating the volume -----------------------
// ---------- enclosed by the level set. -------------------------------------

template <typename RealType>
class LevelSetVolume : public qc::IntegrateFunctionOverLevelSetBase3d
      <RealType, LevelSetVolume<RealType> > {

public:
//    // TODO: only temporarely
  typedef qc::QuocConfiguratorTraitMultiLin < RealType, qc::QC_3D,
  aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfigType;

protected:
  aol::DiscreteFunctionDefault<ConfigType> _discFuncLevelSet;
  RealType h;

public:
  // ---------------------- constructor --------------------------------------------
  LevelSetVolume ( const qc::GridDefinition &Grid,
                   const qc::ScalarArray<RealType, qc::QC_3D> &LevelSetFunction ) :
      qc::IntegrateFunctionOverLevelSetBase3d < RealType,
      LevelSetVolume<RealType> > ( Grid, LevelSetFunction ),
      _discFuncLevelSet ( Grid, LevelSetFunction ) {
    h = Grid.H();
  }

  // ----------------------------methods --------------------------------------------
  // evaluate the function by using the level set data

  RealType evaluateFunction ( const qc::Element &El,
                              const aol::Vec3<RealType> &refCoords )  const {
    // use Gauss and get the volume by integrating over x*n
    // first calc n = grad u / || grad u ||

    aol::Vec3<RealType> grad;
    _discFuncLevelSet.evaluateGradient ( El, refCoords, grad );
    grad /= grad.norm();

    aol::Vec3<RealType> x ( refCoords );

    RealType res = grad * x;
    res /= 3.;                  // because div x = 3 (at least in 3D...)
    return res;
  }


  RealType getVolume ( bool Ausgabe = 0 ) {
    return this->integrate ( Ausgabe );
  }

};



} // end namespace



#endif
