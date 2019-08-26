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
#include <scalarArray.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <quocTimestepSaver.h>
#include <smallVec.h>

#include <FEOpInterface.h>
#include <gradientDescent.h>
#include <deformations.h>
#include <signedDistanceOp.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>


template< typename DataType >
void generate_ball_levelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, DataType radius = 0.3, DataType ctr_x = 0.5, DataType ctr_y = 0.5, DataType ctr_z = 0.5 ) {
  const int width = levelset.getNumX();
  aol::ProgressBar<> pb ( "Creating ball levelset." );
  pb.start ( aol::Cub ( width ) );

  for ( int k = 0; k < width; ++k ) {
    for ( int j = 0; j < width; ++j ) {
      for ( int i = 0; i < width; ++i ) {
        const DataType
          x = ( 1.0 * i ) / ( width - 1.0 ),
          y = ( 1.0 * j ) / ( width - 1.0 ),
          z = ( 1.0 * k ) / ( width - 1.0 ),

          value = aol::Sqr ( x - ctr_x ) + aol::Sqr ( y - ctr_y ) + aol::Sqr ( z - ctr_z ) - aol::Sqr ( radius );

        levelset.set ( i, j, k, value );

      };
    };
  };
}

/**
 * \brief This class is used to compute the vector of sums \f$(\sum_ic_i\zeta(x_i)\phi_j(x_i))_j\f$,
 * where the \f$\phi_j\f$ are the finite element basis functions and the discrete function \f$\zeta\f$ is passed to the
 * method "apply(...)" as a Vector. The points \f$x_i\f$ and weights \f$c_i\f$ are passed to the constructor as a std::vector
 * and aol::Vector, respectively.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class DiscreteScaledMassOp :
  public aol::FELinOpInterface<typename ConfiguratorType::RealType, ConfiguratorType, DiscreteScaledMassOp<ConfiguratorType> > {

private:
  // abbreviations for types
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType  MatType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::ElementType ElementType;

  // define the computation-settings and the grid
  ConfiguratorType _config;
  // contains all points $x_i$ (in local coordinates), ordered by the element in which they appear
  // (allowed, since ElementType has an operator<)
  mutable std::map<ElementType, std::vector<VecType> > _points;
  // contains all weights $c_i$ in the same order as the points
  mutable std::map<ElementType, std::vector<RealType> > _c;

public:
  DiscreteScaledMassOp( const typename ConfiguratorType::InitType &Grid,
                        const std::vector<VecType> &Points,
                        const aol::Vector<RealType> &C,
                        typename aol::OperatorType OpType = aol::ONTHEFLY ) :
      aol::FELinOpInterface<typename ConfiguratorType::RealType,ConfiguratorType,DiscreteScaledMassOp<ConfiguratorType> > ( _config, OpType ),
      _config ( Grid ) {
      // definition/declaration of the number of $x_i$ and of $x_i$ in local/global coordinates
      int numPoints = Points.size();
      ElementType el;
      VecType localCoord;
      // for each $x_i$ do...
      for ( int i = 0; i < numPoints; i++ ) {
        // obtain local coordinates of $x_i$
        _config.getLocalCoords ( Points[i], el, localCoord );
        // save local coordinates
        if ( _points.find( el ) == _points.end() ) {
          _points[el] = std::vector<VecType>();
          _c[el] = std::vector<RealType>();
        }
        _points[el].push_back( localCoord );
        _c[el].push_back( C[i] );
      }
      /*// check, whether points have correctly been saved
      cerr<<_points.size()<<endl;
      typename std::map<ElementType, std::vector<VecType> >::iterator iter;
      VecType glob;
      for( iter = _points.begin(); iter != _points.end(); ++iter ) {
        _config.getGlobalCoords(iter->first,iter->second[0],glob);
        cerr << "key: " << iter->first << ", value: " << iter->second << ", " << glob << endl;
      }
      */
  }

  void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    // local number of degrees of freedom on the element
    const int numDofs = _config.getNumLocalDofs ( El );
    // initialise the local operator matrix with zero
    for ( int i = 0; i < numDofs; i++ )
      for ( int j = 0; j < numDofs; j++ )
        LocalMatrix[i][j] = 0.;
    // if there are points $x_i$ within the element, do...
    if ( _points.find( El ) != _points.end() ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
      const int numPoints = _points[El].size();
      // for each $x_p$ on the element do...
      for ( int p = 0; p < numPoints; p++ )
        // fill the upper half of the discrete local mass matrix
        for ( int i = 0; i < numDofs; i++ ) {
          RealType basisi = bfs.evaluate( i, _points[El][p] );
          for ( int j = i; j < numDofs; j++ )
            LocalMatrix[i][j] +=  basisi * bfs.evaluate( j, _points[El][p] ) * _c[El][p];
        }
      // fill the lower half of the matrix symmetrically
      for ( int i = 0; i < numDofs; i++ )
        for ( int j = i; j < numDofs; j++ )
          LocalMatrix[j][i] = LocalMatrix[i][j] = LocalMatrix[i][j] * _config.vol( El );
    }
  }
};


/**
 * \brief This class computes via "apply(...)": \f$\int_\Omega(||\nabla\zeta||^2-1)^2dx+\sum_ic_i\zeta(x_i)^2+\epsilon\int_\Omega\nabla H(\zeta)dx\f$,
 * where the domain \f$\Omega\f$ is passed to the constructor as a grid "Grid", the points \f$x_i\f$ are passed as
 * a std::vector of coordinate vectors, and \f$\zeta\f$ is passed as an aol::Vector to "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class Energy : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType,Energy<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;

  // for the terms $\sum_ic_i\zeta(x_i)^2$
  const DiscreteScaledMassOp<ConfiguratorType> _discreteMassOp;
  const RealType _discreteMassWeight;
  const RealType _epsilon;
  const aol::ArcTanHeavisideFunction<RealType> _heaviside;
  const aol::HeavisideLevelsetLengthEnergy< ConfiguratorType, aol::ArcTanHeavisideFunction<RealType> > _lengthPenalty;

public:
  Energy( const typename ConfiguratorType::InitType &Grid,
          const std::vector<VecType> &Points,
          const aol::Vector<RealType> &C,
          const RealType DiscreteMassWeight,
          const RealType Epsilon ) :
    // load the grid
    aol::FENonlinIntegrationScalarInterface<ConfiguratorType,Energy<ConfiguratorType> >( Grid ),
    // initialise the mass operator
    _discreteMassOp( Grid, Points, C ),
    _discreteMassWeight( DiscreteMassWeight ),
    _epsilon( Epsilon ),
    _heaviside( Grid.H() ),
    _lengthPenalty( Grid, _heaviside ) {
  }

  /**
   * \brief Returns \f$(||\nabla\zeta||^2-1)^2)\f$ at the point specified by "El", "QuadPoint", "RefCoord". \f$\zeta\f$ is contained in "DiscFuncs".
   */
  inline RealType evaluateIntegrand( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::DomVecType& /*RefCoord*/) const {
    VecType grad;
    DiscFuncs.evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    return aol::Sqr( grad.normSqr()-1 );
  }

  /**
   * \brief Returns \f$\int_\Omega(||\nabla\zeta||^2-1)^2dx+\sum_ic_i\zeta(x_i)^2+\epsilon\int_\Omega\nabla H(\zeta)dx\f$ in "Dest",
   * where the domain \f$\Omega\f$ and the points \f$x_i\f$ were passed to the constructor and \f$\zeta\f$ is given by "Arg".
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // auxiliary vector
    aol::Vector<RealType> aux( Arg );
    //! compute discrete weighted mass operation
    _discreteMassOp.apply( Arg, aux );
    Dest += aux * Arg * _discreteMassWeight;
    //! compute the eikonal energy
    // careful: apply of the parent class would not work; it would call the applyAdd of this class
    aol::FENonlinIntegrationScalarInterface<ConfiguratorType,Energy<ConfiguratorType> >::applyAdd( Arg, Dest );
    //! compute the length penalty
    aol::MultiVector<RealType> tmpArg;
    tmpArg.appendReference( Arg );
    Dest /= _epsilon;
    _lengthPenalty.applyAdd( tmpArg, Dest );
    Dest *= _epsilon;
  }
};

/**
 * \brief This class computes via "apply(...)": \f$(\int_\Omega4(||\nabla\zeta||^2-1)\nabla\zeta\nabla\phi_jdx+\sum_ic_i2\zeta(x_i)\phi_j(x_i)+\epsilon\int_\Omega\nabla^2H(\zeta)\phi_jdx)_j\f$,
 * where the domain \f$\Omega\f$ is passed to the constructor as a grid "Grid", the points \f$x_i\f$ are passed as
 * a std::vector of coordinate vectors, \f$\zeta\f$ is passed as an aol::Vector to "apply(...)", and the \f$\phi_j\f$
 * denote the finite element basis functions belonging to the mesh.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class EnergyVariation : public aol::FENonlinDiffOpInterface<ConfiguratorType,EnergyVariation<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;

  // for the terms $\sum_ic_i\zeta(x_i)^2$
  const DiscreteScaledMassOp<ConfiguratorType> _discreteMassOp;
  const RealType _discreteMassWeight;
  const RealType _epsilon;
  const aol::ArcTanHeavisideFunction<RealType> _heaviside;
  const aol::FullVariationOfHeavisideLevelsetLengthEnergy< ConfiguratorType, aol::ArcTanHeavisideFunction<RealType> > _lengthPenaltyVariation;

public:
  EnergyVariation( const typename ConfiguratorType::InitType &Grid,
                   const std::vector<VecType> &Points,
                   const aol::Vector<RealType> &C,
                   const RealType DiscreteMassWeight,
                   const RealType Epsilon ) :
    // load the grid
    aol::FENonlinDiffOpInterface<ConfiguratorType,EnergyVariation<ConfiguratorType> >( Grid ),
    // initialise the mass and stiffness operator
    _discreteMassOp( Grid, Points, C ),
    _discreteMassWeight( DiscreteMassWeight ),
    _epsilon( Epsilon ),
    _heaviside( Grid.H() ),
    _lengthPenaltyVariation( Grid, _heaviside ) {
  }

  /**
   * \brief Returns \f$4(||\nabla\zeta||^2-1)\nabla\zeta\f$ at the point specified by "El", "QuadPoint", "RefCoord", in "NL". \f$\zeta\f$ is contained in "DiscFuncs".
   */
  inline void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                                    const typename ConfiguratorType::ElementType &El,
                                    int QuadPoint,
                                    const typename ConfiguratorType::DomVecType& /*RefCoord*/,
                                    VecType &NL ) const {
    DiscFuncs.evaluateGradientAtQuadPoint( El, QuadPoint, NL );
    NL *= 4 * ( NL.normSqr()-1 );
  }

  /**
   * \brief Returns \f$(\int_\Omega4(||\nabla\zeta||^2-1)\nabla\zeta\nabla\phi_jdx+\sum_ic_i2\zeta(x_i)\phi_j(x_i)+\epsilon\int_\Omega\nabla^2H(\zeta)\phi_jdx)_j\f$ in "Dest",
   * where the domain \f$\Omega\f$ and the points \f$x_i\f$ were passed to the constructor and \f$\zeta\f$ is given by "Arg".
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    Dest /= 2*_discreteMassWeight;
    //! compute discrete weighted mass operator
    _discreteMassOp.applyAdd( Arg, Dest );
    Dest *= 2*_discreteMassWeight;
    //! compute the eikonal operator
    // careful: apply of the parent class would not work; it would call the applyAdd of this class
    aol::FENonlinDiffOpInterface<ConfiguratorType,EnergyVariation<ConfiguratorType> >::applyAdd( Arg, Dest );    //! compute the length penalty
    //! compute the length penalty term
    aol::MultiVector<RealType> tmpArg, tmpDest;
    tmpArg.appendReference( Arg );
    tmpDest.appendReference( Dest );
    Dest /= _epsilon;
    _lengthPenaltyVariation.applyAdd( tmpArg, tmpDest );
    Dest *= _epsilon;
  }
};

template<typename ConfiguratorType>
void findLevelSet( const char PointsFileName[], const int id = 1 ) {
  // read in point cloud $x_i$ of dxf-file
  std::vector<typename ConfiguratorType::VecType> points;
  char line[256];
  ifstream dxfFile( PointsFileName );
  while ( !dxfFile.eof() ) {
    dxfFile.getline( line, 256 );
    if ( !strcmp( line, "POINT\r" ) ) {
      typename ConfiguratorType::VecType point;
      for ( int i = 0; i < 4; i++ )
        dxfFile.getline ( line, 256 );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
        dxfFile.getline ( line, 256 );
        dxfFile.getline ( line, 256 );
        point[i] = strtod( line, NULL );
      }
      // scale point into the interval $[0,1]^n$ (best by a constant factor to obtain reproducible results)
      point /= 400.;
      point += typename ConfiguratorType::VecType( .5, .5, .1 );
      points.push_back( point );
    }
  }
  dxfFile.close();
  /*// read in point cloud $x_i$ of vtk-polydata file
  std::vector<typename ConfiguratorType::VecType> points;
  char line[256], word[256] = "";
  ifstream vtkFile( PointsFileName );
  // neglect header lines until points start
  while ( strcmp( word, "POINTS" ) && !vtkFile.eof() ) {
    vtkFile.getline( line, 256 );
    strncpy( word, line, 6 );
  }
  while ( !vtkFile.fail() && !vtkFile.eof() ) {
    typename ConfiguratorType::VecType point;
    vtkFile >> point[0];
    if ( vtkFile.fail() )
      break;
    for ( int i = 1; i < ConfiguratorType::Dim; i++ )
      vtkFile >> point[i];
    points.push_back( point );
  }
  vtkFile.close();*/
  // rescale the point cloud to fit into [0,1]^d
  typename ConfiguratorType::RealType minX = points[0][0], minY = points[0][1], minZ = points[0][2], maxX = minX, maxY = minY, maxZ = minZ;
  int numPoints = points.size();
  for ( int i = 0; i < numPoints; i++ ) {
    // compute range of values
    minX = aol::Min( minX, points[i][0] );
    minY = aol::Min( minY, points[i][1] );
    minZ = aol::Min( minZ, points[i][2] );
    maxX = aol::Max( maxX, points[i][0] );
    maxY = aol::Max( maxY, points[i][1] );
    maxZ = aol::Max( maxZ, points[i][2] );
  }
  typename ConfiguratorType::VecType offset( ( minX + maxX ) / 2, ( minY + maxY ) / 2, ( minZ + maxZ ) / 2 );
  typename ConfiguratorType::VecType shift( .5, .5, .5 );
  typename ConfiguratorType::RealType scale = 1.5 * aol::Max( aol::Max( maxX - minX, maxY - minY ), maxZ - minZ );
  scale = 1;
  for ( int i = 0; i < numPoints; i++ ) {
    points[i] -= offset;
    points[i] /= scale;
    points[i] += shift;
  }

  // set weights $c_i$
  aol::Vector<typename ConfiguratorType::RealType> c( points.size() );
  c.setAll( 1. );

  int lowestLevel = 5, highestLevel = 7;

  // initialise the level set function $\zeta$
  typename ConfiguratorType::InitType coarsestGrid( lowestLevel, ConfiguratorType::Dim );
  typename ConfiguratorType::ArrayType initialZeta( coarsestGrid );
  /*if ( ConfiguratorType::Dim == 2 ) {
    qc::DataGenerator<ConfiguratorType> levelSetFunctionInitialiser( coarsestGrid );
    levelSetFunctionInitialiser.generateCircleLevelset( initialZeta, .25, .5, .5 );
  } else
    generate_ball_levelset<typename ConfiguratorType::RealType> ( initialZeta, .1, .5, .5, .85 );
  */
  initialZeta.setAll( 1. );

  /*// save points as ply-file
  sm::SurfMesh<typename ConfiguratorType::RealType> surf;
  surf.reserve(numPoints,0);
  for (int i=0;i<numPoints;i++) surf.push_back_vertex(points[i][0],points[i][1],points[i][2]);
  surf.writePlyFile("surf.ply");*/

  // perform multilevel energy minimisation
  qc::MultilevelArray<typename ConfiguratorType::RealType> zeta( highestLevel, ConfiguratorType::Dim );

  /*// load initial data
  typename ConfiguratorType::InitType finestGrid( highestLevel, ConfiguratorType::Dim );
  typename ConfiguratorType::ArrayType(zeta.current(),aol::FLAT_COPY).load("./feet/feet_example/China/levelset20.bz2");
  zeta.current().addToAll(-zeta.current().getMinValue());
  zeta.current() /= zeta.current().getMaxValue();
  zeta.current().addToAll(-.5);
  qc::SignedDistanceOp3D<typename ConfiguratorType::RealType> signedDistOp( finestGrid );
  signedDistOp.apply( zeta.current(), zeta.current() );*/

  zeta.levRestrict( lowestLevel, highestLevel );
  zeta.setCurLevel( lowestLevel );
  zeta.current() = initialZeta;
  for ( int level = lowestLevel; level <= highestLevel; level++ ) {
    // produce a corresponding grid and choose correct level
    typename ConfiguratorType::InitType grid( level, ConfiguratorType::Dim );
    zeta.setCurLevel( level );
    typename ConfiguratorType::ArrayType curZeta( zeta.current(), aol::DEEP_COPY );

    // perform the energy minimisation, finding the level set
    typename ConfiguratorType::RealType discreteMassWeight = 10000., surfacePenaltyWeight = 1e-3;
    Energy<ConfiguratorType> E( grid, points, c, discreteMassWeight, surfacePenaltyWeight );
    EnergyVariation<ConfiguratorType> DE( grid, points, c, discreteMassWeight, surfacePenaltyWeight );
    aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType,aol::Vector<typename ConfiguratorType::RealType> > gradientDescent( grid, E, DE, 200 );
    gradientDescent.apply( curZeta, zeta.current() );

    // reinitialize as a signed distance function
    if ( ConfiguratorType::Dim == 2 ) {
      qc::SignedDistanceOp<ConfiguratorType> signedDistOp( grid );
      signedDistOp.apply( zeta.current(), zeta.current() );
    } else {
      qc::SignedDistanceOp3D<typename ConfiguratorType::RealType> signedDistOp( grid );
      signedDistOp.apply( zeta.current(), zeta.current() );
    }

    // prolongate
    if ( level < highestLevel )
      zeta.levProlongate();

    // store the results
    char fileName[1024];
    // due to compression, file extension should be bz2
    sprintf( fileName, "%slevelset_%d.%s", "./", id, "bz2" );
    (typename ConfiguratorType::ArrayType(zeta.current(),aol::FLAT_COPY)).save( fileName, qc::PGM_DOUBLE_BINARY );
    //qc::QuocTimestepSaver<ConfiguratorType> timeStepSaver( 1, 0, "levelS", true );
    //timeStepSaver.setSaveDirectory( "./" );
    //timeStepSaver.saveTimestepPNGDrawRedAndBlueIsoline( 0, 0, 0., 0.1, zeta, zeta, grid );
  }
}

//////////////////////////////////////////////////////////////////////////////////////
//! ---------------------------- main program --------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype,
// finite element type, used quadrature rules
typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;

int main( int argc, char *argv[] ) {

  try {
   // read in file names of data files
    if ( argc > 50 ) {     // wrong number of input files
      cerr << "Please specify between 1 and 50 input files." << endl;
      return EXIT_FAILURE;
    } else if ( argc < 2 )// no file specified
      findLevelSet<ConfiguratorType>( "JEREMYLLascii.dxf" );
    else                  // read in filenames of data files and perform the algorithm
      for ( int i=1; i<argc; i++ )
        findLevelSet<ConfiguratorType>( argv[i], i );
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
