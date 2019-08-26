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
#include <smallVec.h>

#include <FEOpInterface.h>
#include <Newton.h>
#include <deformations.h>
#include <signedDistanceOp.h>
#include <gradientDescent.h>

#include <prolongation.h>

/**
 * \brief This class computes via "apply(...)": \f$ \int_\Omega (u-v\circ\phi)^2 dx \f$ in "Dest",
 * where the domain \f$ \Omega \f$ and the images \f$ u,v \f$ are passed to the constructor, and \f$ \phi(x) = Ax+b \f$
 * with \f$ b \f$ and \f$ A \f$ are passed to "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class AffineImageAlignmentEnergy :
  public aol::Op< aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType VectorType;

  // image grid
  const typename ConfiguratorType::InitType _grid;
  // auxiliary operator
  const aol::MassOp<ConfiguratorType> _massOp;
  // reference image
  const VectorType _u;
  // template image
  const VectorType _v;

public:
  AffineImageAlignmentEnergy( const typename ConfiguratorType::InitType &Grid,
                              const VectorType &U,
                              const VectorType &V ) :
    _grid( Grid ),
    _massOp( Grid ),
    _u( U ),
    _v( V ) {
  }

  /**
   * \brief Returns \f$ \int_\Omega (u-v\circ\phi)^2 dx \f$ in "Dest",
   * where the domain \f$ \Omega \f$ and the images \f$ u,v \f$ were passed to the constructor and \f$ \phi(x) = Ax+b \f$
   * with \f$ b \f$ being the last column of "Arg" and \f$ A \f$ the rest.
   */
  void applyAdd( const aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // initialize the vector $v\circ\phi$ $u-v\circ\phi$ and an auxiliary vector
    VectorType deformedV( _v, aol::STRUCT_COPY ), difference( _u, aol::DEEP_COPY ), tmp( _u, aol::STRUCT_COPY );
    // obtain matrix $A$ and vector $b$
    typename ConfiguratorType::VecType shift;
    typename ConfiguratorType::MatType transformMatrix;
    Arg.getCol( ConfiguratorType::Dim, shift );
    Arg.getSubMatrix( 0, 0, transformMatrix );
    // compute $u-v\circ\phi$
    aol::MultiVector<RealType> displacement( ConfiguratorType::Dim, _v.size() );
    qc::DataGenerator<ConfiguratorType> generator( _grid );
    generator.generateAffineDisplacement( transformMatrix, shift, displacement );
    qc::DeformImage<ConfiguratorType>( _v, _grid, deformedV, displacement, difference );
    difference -= deformedV;
    // return $\int_\Omega (u-v\circ\phi$)^2 dx$
    _massOp.apply( difference, tmp );
    Dest += difference * tmp;
  }
};

/**
 * \brief This class computes via "apply(...)" the gradient of \f$ \int_\Omega (u-v\circ\phi)^2 dx \f$ wrt \f$ A \f$ and \f$ b \f$ in "Dest",
 * where the domain \f$ \Omega \f$ and the images \f$ u,v \f$ are passed to the constructor, and \f$ \phi(x) = Ax+b \f$
 * with \f$ b \f$ and \f$ A \f$ are passed to "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class AffineImageAlignmentGradient :
  public aol::Op< aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType VectorType;

  // image grid
  const typename ConfiguratorType::InitType _grid;
  // reference image
  const VectorType _u;
  // template image
  const VectorType _v;

public:
  AffineImageAlignmentGradient( const typename ConfiguratorType::InitType &Grid,
                                const VectorType &U,
                                const VectorType &V ) :
    _grid( Grid ),
    _u( U ),
    _v( V ) {
  }

  /**
   * \brief Returns gradient of \f$ \int_\Omega (u-v\circ\phi)^2 dx \f$ in "Dest",
   * where the domain \f$ \Omega \f$ and the images \f$ u,v \f$ were passed to the constructor and \f$ \phi(x) = Ax+b \f$
   * with \f$ b \f$ being the last column of "Arg" and \f$ A \f$ the rest.
   */
  void applyAdd( const aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Arg, aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Dest ) const {
    // initialize the vector $v\circ\phi$ $u-v\circ\phi$ and an auxiliary vector
    VectorType deformedV( _v, aol::STRUCT_COPY ), difference( _u, aol::DEEP_COPY ), tmp( _u, aol::STRUCT_COPY );
    // obtain matrix $A$ and vector $b$
    typename ConfiguratorType::VecType shift;
    typename ConfiguratorType::MatType transformMatrix;
    Arg.getCol( ConfiguratorType::Dim, shift );
    Arg.getSubMatrix( 0, 0, transformMatrix );
    // compute $u-v\circ\phi$
    aol::MultiVector<RealType> displacement( ConfiguratorType::Dim, _v.size() );
    qc::DataGenerator<ConfiguratorType> generator( _grid );
    generator.generateAffineDisplacement( transformMatrix, shift, displacement );
    qc::DeformImage<ConfiguratorType>( _v, _grid, deformedV, displacement, difference );
    difference -= deformedV;
    // return $\int_\Omega -2(u-v\circ\phi)\nabla v dx$ and $\int_\Omega -2(u-v\circ\phi)\nabla v x^T dx$
    VectorType ones( _u, aol::STRUCT_COPY );
    aol::MultiVector<RealType> id( ConfiguratorType::Dim, _u.size() );
    ones.addToAll( 1 );
    generator.generateIdentity( id );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      aol::Scalar<RealType> deriv;
      aol::WeightedSemiDiffOp<ConfiguratorType> derivOp( _grid, ones, i );
      derivOp.apply( deformedV, tmp );
      Dest( i, ConfiguratorType::Dim ) += tmp * _v;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        aol::Scalar<RealType> deriv;
        aol::WeightedSemiDiffOp<ConfiguratorType> derivOp( _grid, id[j], i );
        derivOp.apply( deformedV, tmp );
        Dest( i, j ) += tmp * _v;
      }
    }
    Dest *= -2;
  }
};

template <typename ConfiguratorType>
class WeightedDiffVectorIntegrator :
  public aol::FENonlinIntegrationVectorGeneralInterface< ConfiguratorType, WeightedDiffVectorIntegrator<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  // the weight
  aol::DiscreteFunctionDefault<ConfiguratorType> _w;
  // the _i-th row and _j-th column of the gradient shall be integrated
  int _i, _j;

public:
  WeightedDiffVectorIntegrator( typename ConfiguratorType::InitType &Grid, typename ConfiguratorType::ArrayType &W, int I, int J ) :
    aol::FENonlinIntegrationVectorGeneralInterface< ConfiguratorType, WeightedDiffVectorIntegrator<ConfiguratorType> >( Grid ),
    _w( Grid, W ),
    _i( I ),
    _j( J ) {}

  RealType evaluateIntegrand ( const aol::auto_container< ConfiguratorType::Dim, aol::DiscreteFunctionDefault< ConfiguratorType > > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/) const {
    typename ConfiguratorType::VecType df;
    DiscFuncs[_i].evaluateGradientAtQuadPoint ( El, QuadPoint, df );
    return _w.evaluateAtQuadPoint( El, QuadPoint ) * df[_j];
  }
};

template <typename ConfiguratorType>
class WeightedVectorIntegrator :
  public aol::FENonlinIntegrationVectorGeneralInterface< ConfiguratorType, WeightedVectorIntegrator<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  // the weight
  aol::DiscreteFunctionDefault<ConfiguratorType> _w;
  // the _i-th row shall be integrated
  int _i;

public:
  WeightedVectorIntegrator( typename ConfiguratorType::InitType &Grid, typename ConfiguratorType::ArrayType &W, int I ) :
    aol::FENonlinIntegrationVectorGeneralInterface< ConfiguratorType, WeightedVectorIntegrator<ConfiguratorType> >( Grid ),
    _w( Grid, W ),
    _i( I ) {}

  RealType evaluateIntegrand ( const aol::auto_container< ConfiguratorType::Dim, aol::DiscreteFunctionDefault< ConfiguratorType > > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/) const {
    return _w.evaluateAtQuadPoint( El, QuadPoint ) * DiscFuncs[_i].evaluateAtQuadPoint( El, QuadPoint );
  }
};


//////////////////////////////////////////////////////////////////////////////////////
//! ---------------------------- main program --------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////

// define the settings/configuration of the program, i.e. computation accuracy, dimension, gridtype,
// finite element type, used quadrature rules
typedef double RealType;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;

int main( int argc, char *argv[] ) {

  /*try*/ {
    if ( argc < 5 ) {     // wrong number of input files
      cerr << "Usage: " << argv[0] << " <image number> <image wordstem> <displacement wordstem> <result directory>" << endl;
      return EXIT_FAILURE;
    } else {

      /*// the reference image
      ConfiguratorType::ArrayType u( argv[1] );
      u.addToAll( -.5 );
      signedDistOp.apply( u, u );

      for ( int i = 3; i <= argc; i++ ) {
        // the template image
        ConfiguratorType::ArrayType v( argv[i] ), deformedV( v, aol::STRUCT_COPY );
        v.addToAll( -v.getMinValue() );
        v /= v.getMaxValue();
        v.addToAll( -.5 );
        signedDistOp.apply( v, v );
        // compute the affine alignment transformation
        aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> initialX, x;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          initialX( j, j ) = 1;
        AffineImageAlignmentEnergy<ConfiguratorType> e( grid, u, v );
        AffineImageAlignmentGradient<ConfiguratorType> de( grid, u, v );
        aol::GridlessGradientDescent< RealType, aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> > gradientDescent( e, de, 100 );
        gradientDescent.apply( initialX, x );

        // obtain matrix $A$ and vector $b$
        ConfiguratorType::VecType column;
        x.getCol( ConfiguratorType::Dim, column );
        const ConfiguratorType::VecType shift( column );
        ConfiguratorType::MatType transformMatrix;
        for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
          x.getCol( i, column );
          transformMatrix.setCol( i, column );
        }
        // transform the image
        aol::MultiVector<RealType> displacement( ConfiguratorType::Dim, v.size() );
        qc::DataGenerator<ConfiguratorType> generator( grid );
        generator.generateAffineDisplacement( transformMatrix, shift, displacement );
        qc::DeformImage<ConfiguratorType>( v, grid, deformedV, displacement, true );

        // save the transformed image
        deformedV.save( "ImgTem.bz2", qc::PGM_DOUBLE_BINARY );
        u.save( "ImgRef.bz2", qc::PGM_DOUBLE_BINARY );
      }*/

      /*aol::Vector<RealType> defShift( ConfiguratorType::Dim * ( ConfiguratorType::Dim + 1 ) );
      aol::Scalar<RealType> aux;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          (WeightedDiffVectorIntegrator<ConfiguratorType>( grid, vFine, i, j )).apply( d, aux );
          defShift[i*ConfiguratorType::Dim+j] = aux;
        }
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        defShift[i*(ConfiguratorType::Dim+1)] += 1.0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
        (WeightedVectorIntegrator<ConfiguratorType>( grid, vFine, i )).apply( d, aux );
        defShift[ConfiguratorType::Dim*ConfiguratorType::Dim+i-1] = aux;
      }
      cerr << defShift << endl;
      sprintf( filename, "%s/DefShift%d.bz2", argv[4], image + 1 + initialImage );
      defShift.saveAsVector( filename );*/


      // find the rigid body motion for all input shapes
      int initialImage = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule( static, 4 )
#endif
      for ( int image = 0; image < atoi( argv[1] ); image++ ) {

        // load the original volumetric shape
        char filename[1024];
        sprintf( filename, "%s%d.bz2", argv[2], image + 1 + initialImage );
        cerr << "load " << filename << endl;
        ConfiguratorType::ArrayType sh( filename );
        /*ConfiguratorType::ArrayType shAux( filename ), sh( 257, 257, 257 );
        ConfiguratorType::InitType coarseGrid( 7, ConfiguratorType::Dim );
        ConfiguratorType::InitType fineGrid( 8, ConfiguratorType::Dim );
        qc::ProlongOp< RealType >( coarseGrid, fineGrid ).apply( shAux, sh );*/

        // define an appropriate grid structure
        int depth = 0, width = sh.getNumX();
        while ( width > 1 ) {
          depth++;
          width = width >> 1;
        }
        ConfiguratorType::InitType grid( depth, ConfiguratorType::Dim );

        // load the deformation
        aol::MultiVector<RealType> d( ConfiguratorType::Dim, sh.size() );
        for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
          sprintf( filename, "%s%d_%d.bz2", argv[3], image + initialImage, i );
          cerr << "load " << filename << endl;
          ( ConfiguratorType::ArrayType( d[i], grid, aol::FLAT_COPY ) ).load( filename );
        }

        // put all data into a multilevel structure
        qc::MultilevelArray<RealType> shape( grid );
        qc::MultiDimMultilevelArray<RealType> deformation( grid, ConfiguratorType::Dim );
        shape.current() = sh;
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          deformation.getArray( i ) = d[i];

        // restrict to the initial level of the multilevel approach
        shape.levRestrict( depth / 2, depth );
        deformation.levRestrict( depth / 2, depth );

        // compute the rigid alignment transformation
        aol::Vector<RealType> initialX( 6 ), x( 6 );
        ConfiguratorType::VecType offset;
        offset.setAll( .5 );
        for ( int level = depth / 2; level <= depth; level++ ) {
          // obtain the current shape and deformation
          ConfiguratorType::InitType currentGrid( level, ConfiguratorType::Dim );
          shape.setCurLevel( level );
          deformation.setCurLevel( level );
          aol::MultiVector<RealType> d;
          for ( int dir = 0; dir < ConfiguratorType::Dim; dir++ )
            d.appendReference( deformation.getArray( dir ) );

          // find the best fit rigid body motion
          qc::RigidDisplacementProjectionEnergy<ConfiguratorType> e( currentGrid, d, shape.current(), offset );
          qc::RigidDisplacementProjectionGradient<ConfiguratorType> de( currentGrid, d, shape.current(), offset );
          aol::NewtonInfo<RealType> newtonInfo ( 1E-10, 300, 1E-20, 1000, aol::STOPPING_ABSOLUTE );
          aol::QuasiNewtonIteration< RealType, aol::Vector<RealType>, aol::Vector<RealType> > descent( e, de, newtonInfo );
          descent.apply( initialX, x );
          initialX = x;
        }

        // obtain matrix $A$ and vector $b$
        int index1[3] = { 0, 0, 1 };
        int index2[3] = { 1, 2, 2 };
        int numberAngles = x.size() - ConfiguratorType::Dim;
        ConfiguratorType::MatType transformMatrix, rotationMatrix;
        transformMatrix.setIdentity();
        for ( int i = 0; i < numberAngles; i++ ) {
          RealType alpha = x[ConfiguratorType::Dim + i];
          rotationMatrix.setIdentity();
          rotationMatrix[index1[i]][index1[i]] = cos( alpha );
          rotationMatrix[index2[i]][index2[i]] = cos( alpha );
          rotationMatrix[index1[i]][index2[i]] = -sin( alpha );
          rotationMatrix[index2[i]][index1[i]] = sin( alpha );
          transformMatrix *= rotationMatrix;
        }

        ConfiguratorType::VecType shift, invShift;
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          shift[i] = x[i];
        shift += offset;
        shift *= -1;
        transformMatrix.multAdd( offset, shift );

        // transform the image affinely
        rotationMatrix = transformMatrix.inverse();
        rotationMatrix.mult( shift, invShift );
        aol::MultiVector<RealType> displacement( ConfiguratorType::Dim, shape.current().size() );
        ConfiguratorType::ArrayType deformedShape( shape.current(), aol::STRUCT_COPY );
        ( qc::DataGenerator<ConfiguratorType>( grid ) ).generateAffineDisplacement( rotationMatrix, invShift, displacement );
        qc::DeformImage<ConfiguratorType>( shape.current(), grid, deformedShape, displacement, true );

        // compute norm of difference between displacement and projection into affine displacement space
        ConfiguratorType::ArrayType dispNorm( shape.current(), aol::STRUCT_COPY );
        ( qc::DataGenerator<ConfiguratorType>( grid ) ).generateAffineDisplacement( transformMatrix, shift, displacement );
        displacement -= d;
        for ( int i = 0; i < d[0].size(); i++ ) {
          for ( int j = 0; j < ConfiguratorType::Dim; j++ )
            dispNorm[i] += aol::Sqr( displacement[j][i] );
          dispNorm[i] = sqrt( dispNorm[i] );
        }

        // save the transformed image
        sprintf( filename, "%s/deformedImage%d.bz2", argv[4], image + 1 + initialImage );
        deformedShape.save( filename, qc::PGM_DOUBLE_BINARY );
        sprintf( filename, "%s/displacement%d.bz2", argv[4], image + 1 + initialImage );
        dispNorm.save( filename, qc::PGM_DOUBLE_BINARY );
        sprintf( filename, "%s/RotShift%d.bz2", argv[4], image + 1 + initialImage );
        x.saveAsVector( filename );
      }
    }
  }
  /*catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }*/

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
