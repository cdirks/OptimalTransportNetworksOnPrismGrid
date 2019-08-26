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

#ifndef __UTILITIES_H
#define __UTILITIES_H


#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <auxiliary.h>
#include <deformations.h>
#include <imageTools.h>

//(h(x)-g^{-1}(x))|\det(D[g^{-1}](x))|
template <typename ConfiguratorType>
class SymTransSubDetOp:
      public aol::FENonlinVectorOpInterface<ConfiguratorType,
      ConfiguratorType::Dim,
      ConfiguratorType::Dim,
  SymTransSubDetOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the input g transform
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,  ConfiguratorType::Dim> _InputTrans;
  const qc::GridDefinition &_grid;
public:
  // initialization constructor
  SymTransSubDetOp( const qc::GridDefinition &Grid,
                    const aol::MultiVector<RealType> &InputTransDofs)
      : aol::FENonlinVectorOpInterface<ConfiguratorType,
      ConfiguratorType::Dim,
      ConfiguratorType::Dim,
      SymTransSubDetOp<ConfiguratorType> > ( Grid ),
      _InputTrans( Grid,InputTransDofs),
  _grid( Grid ) {
  }

  // compute
  void getNonlinearity( aol::auto_container< ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::VecType &/*RefCoord*/,
                        aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {

    typename ConfiguratorType::MatType mat;
    typename ConfiguratorType::VecType tran_1;
    typename ConfiguratorType::VecType tran_2;
    typename ConfiguratorType::VecType grad;

    // compute h(x)
    for ( int i=0; i<ConfiguratorType::Dim; i++ ) {
      tran_1[i] = DiscFuncs[i].evaluateAtQuadPoint( El, QuadPoint );
    }

    //compute g(x)
    _InputTrans.evaluateAtQuadPoint(El, QuadPoint,tran_2);

    // compute D[g]
    _InputTrans.evaluateGradientAtQuadPoint( El, QuadPoint, mat );
    for ( int i = 0; i < ConfiguratorType::Dim; i ++ ) {
      mat[i][i] += 1.;
    }

    //(h(x)-g^{-1}(x))
    // The identity mapping, that would have to be added is canceled out by the subtraction
    // VERY bad approximation of the inverse! Should be changed!
    NL = tran_1 + tran_2;

    //(h(x)-g^{-1}(x))|\det(D[g^{-1}](x))|
    NL /= fabs(mat.det());
  }
};

namespace qc {

/** @brief The function that exports 1D histogram (1d vector) to data file, that can be visualized via gnuplot
 */
template <typename RealType>
void writeHistogram(const char *fileName, const aol::Vector<RealType> &vec) {

  // number of bin and step size
  const int numOfBin = vec.getSize();
  const RealType step = 1.0/static_cast<RealType>(numOfBin);
  RealType curStep = 0.0;

  // open file
  ofstream out( fileName );
  if (!out) {
    throw aol::FileException( "cannot open file for output.", __FILE__, __LINE__ );
  }

  int i;
  for(i=0;i<numOfBin;i++,curStep+=step) {
    out << curStep << " " << vec[i]<<endl;
  }
}

/** @brief The function that exports 1D histogram (1d vector) to data file, that can be visualized via gnuplot
 */
template <typename RealType>
void writeHistogram(const char *fileName, qc::ScalarArray<RealType, qc::QC_2D> &Table ,int fr) {

  // number of bin and step size
  const int numOfBin = Table.getNumX();
  const RealType step = 1.0/static_cast<RealType>(numOfBin);
  RealType curStep = 0.0;

  if(fr<0 || fr>1) {
    throw aol::Exception( "fr needs to be equal to zero or one", __FILE__, __LINE__);
  }
  // open file
  ofstream out( fileName );
  if (!out) {
    throw aol::FileException( "cannot open file for output.", __FILE__, __LINE__ );
  }

  int i;
  if(fr==0) {
    for(i=0;i<numOfBin;i++,curStep+=step) {
      out << curStep << " " <<Table.getClip(i,0)<<endl;
    }
  } else {
    for(i=0;i<numOfBin;i++,curStep+=step) {
      out << curStep << " " <<Table.getClip(0,i)<<endl;
    }
  }
}

/** @brief Compute absoluted difference between two images
*/
template <typename ConfiguratorType>
double ComputeDifference( const qc::GridDefinition &Grid,
                          const qc::Array<typename ConfiguratorType::RealType> &Image1,
                          const qc::Array<typename ConfiguratorType::RealType> &Image2,
                          qc::Array<typename ConfiguratorType::RealType> &ImageDiff,
                          int style=0) {
  //const int num=Grid.getWidth();
  //ConfiguratorType::RealType h = Grid.H();
  GridDefinition::OldFullNodeIterator fnit;
  double Sum_diff=0.;
  typename ConfiguratorType::RealType diff=0.;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper(Grid);
  if(style==0) {
    for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {
      diff =aol::Abs(Image1.get(*fnit)-Image2.get(*fnit));
      Sum_diff+=diff;
      ImageDiff.set(*fnit,diff);
      //ImageDiff.set(*fnit, aol::Abs(Image1.get(*fnit)-Image2.get(*fnit)));
    }
  } else {
    for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {
      diff =aol::Sqr(Image1.get(*fnit)-Image2.get(*fnit));
      Sum_diff +=diff;
      ImageDiff.set(*fnit,diff);
      //ImageDiff.set(*fnit,aol::Abs(Image1.get(*fnit)-Image2.get(*fnit)));
    }
  }
  return Sum_diff;
}

/** @brief Compute difference between ref images and deformed image
*/
template <typename ConfiguratorType>
double ComputeDifference( const qc::GridDefinition &Grid,
                          const qc::GridDefinition &Cur_grid,
                          const qc::Array<typename ConfiguratorType::RealType> &ImageRef,
                          const qc::Array<typename ConfiguratorType::RealType> &ImageTem,
                          const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                          qc::Array<typename ConfiguratorType::RealType> &ImageDiff,
                          int style=0) {
  typename qc::Array<typename ConfiguratorType::RealType> DeformedImage(Grid);
  deformImageWithCoarseDeformation<ConfiguratorType>( ImageTem, Grid, Cur_grid, DeformedImage,Phi);
  //DeformImage<ConfiguratorType>(ImageTem,Grid,DeformedImage, Phi);
  return(ComputeDifference<ConfiguratorType>(Grid,ImageRef,DeformedImage,ImageDiff,style));

}
/** @brief Output the diff image
*/
template<typename ConfiguratorType>
double writeDiffImage(const qc::GridDefinition &Grid,
                      const qc::GridDefinition &Cur_grid,
                      const char *Filename,
                      const qc::Array<typename ConfiguratorType::RealType> &ImageRef,
                      const qc::Array<typename ConfiguratorType::RealType> &ImageTem,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                      int style=0) {
  typename qc::Array<typename ConfiguratorType::RealType> DiffImage(Grid);
  double diff= ComputeDifference<ConfiguratorType>( Grid,Cur_grid,ImageRef,ImageTem,Phi,DiffImage,style) ;
  // double diff= ComputeDifference<ConfiguratorType>( Grid,ImageRef,ImageTem,DiffImage,style) ;
  qc::writeImage<typename ConfiguratorType::RealType> ( Grid, DiffImage, Filename );
  return diff;
}

/**@brief Compute the max/min of determine of Jacobin matrix of transform*/
template<typename ConfiguratorType>
void JaocibiEvaluation(const qc::GridDefinition &Grid,
                       const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                       double& JacMin, double& JacMax) {
  qc::GridDefinition::OldFullElementIterator  it;
  qc::Element el;
  int QuadPoint=0;
  double JacDet;
  JacMin=1;
  JacMax=1;
  typename ConfiguratorType::VecType nulVec;
  typename ConfiguratorType::MatType Jacobi;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> funPhi( Grid, Phi );

  for (it = Grid.begin(); it != Grid.end(); ++it ) {


    funPhi.evaluateGradientAtQuadPoint(*it, QuadPoint ,Jacobi);

    for ( int i = 0; i < ConfiguratorType::Dim; i ++ ) {
      Jacobi[i][i] += 1.;
    }

    JacDet = Jacobi.det();
    JacMin = (JacMin>JacDet)?JacDet:JacMin;
    JacMax = (JacMax>JacDet)?JacMax:JacDet;
  }
}

}  // end namespace qc

template <typename ConfiguratorType>
class ConsistencyEnergy
: public aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                                    ConsistencyEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const qc::GridDefinition &_grid;
  aol::auto_container< ConfiguratorType::Dim, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrPhiInverse;
protected:
public:
  ConsistencyEnergy( const qc::GridDefinition &Grid,
                     const aol::MultiVector<RealType> &PhiInverseDofs)
      : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                                   ConsistencyEnergy<ConfiguratorType> > ( Grid ),
  _grid( Grid ) {
    for ( int c=0; c<PhiInverseDofs.numComponents(); c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( this->getConfigurator(), PhiInverseDofs[c] );
      _discrPhiInverse.set_copy( c, temp );
    }
  };

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    typename ConfiguratorType::VecType offset, transformed_local_coord;
    qc::Element transformed_el;
    RealType energy = 0.;

    // calculate \Abs{g(h(x))-x}^2
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      offset[i] = discrFuncs[i].evaluateAtQuadPoint( El, QuadPoint );
    }
    if( !qc::transformCoord<ConfiguratorType> ( _grid, El, RefCoord, offset, transformed_el, transformed_local_coord ) ){
      return 0;
    }

    // Since _discrPhiInverse only contains the displacement, we have to omit "-x",
    // but add the offset (comes from evaluating x at the deformed position), i.e.
    // g( h(x) ) - x = ( x + u(x) ) + v( x + u(x) ) - x
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      energy += aol::Sqr(_discrPhiInverse[i].evaluate( transformed_el, transformed_local_coord ) + offset[i]);
    }

    // calculate \Abs{h(g(x))-x}^2
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      offset[i] = _discrPhiInverse[i].evaluateAtQuadPoint( El, QuadPoint );
    }
    if( !qc::transformCoord<ConfiguratorType> ( _grid, El, RefCoord, offset, transformed_el, transformed_local_coord ) ){
      return 0;
    }
    // Since discrFuncs only contains the displacement, we have to omit "-x",
    // but add the offset (comes from evaluating x at the deformed position), i.e.
    // g( h(x) ) - x = ( x + u(x) ) + v( x + u(x) ) - x
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      energy += aol::Sqr(discrFuncs[i].evaluate( transformed_el, transformed_local_coord ) + offset[i]);
    }
    return 0.5 * energy;
  }
};

template <typename RealType>
void fillWithBox( qc::ScalarArray<RealType, qc::QC_2D> &Array, const int EndColumn, const RealType Angle ) {
  Array.setZero();
  if( EndColumn >= Array.getNumX() )
    throw aol::Exception( "EndColumn >= Array.getNumX() !", __FILE__, __LINE__);
  if( EndColumn < 0 )
    throw aol::Exception( "EndColumn < 0 !", __FILE__, __LINE__);
  for( int j = 0; j < Array.getNumY(); j++ ) {
    for( int i = 0; i < Array.getNumX(); i++ ) {
      if( ( cos( Angle )*(i-EndColumn)) > (sin( Angle ) *(j-Array.getNumY()/2.)) ) {
        Array.set(i, j, 1.);
      }
    }
  }
}

template <typename RealType>
void addFilledBox( qc::ScalarArray<RealType, qc::QC_2D> &Array, const int StartX, const int EndX, const int StartY, const int EndY, const RealType Value = aol::ZOTrait<RealType>::one ) {
  if( EndX >= Array.getNumX() || StartX < 0 || EndY >= Array.getNumY() || StartX < 0 )
    throw aol::Exception( "Arguments out of range!", __FILE__, __LINE__);
  for( int j = StartY; j <= EndY ; j++ ) {
    for( int i = StartX; i <= EndX; i++ ) {
      Array.set(i, j, Value);
    }
  }
}

template <typename RealType>
void clipAndLogaritmicScaleTo01( qc::ScalarArray<RealType, qc::QC_2D> &Array, const RealType Min, const RealType Max, const RealType LogarithmBase = 10.) {
  // Higher values of LogarithmBase will pronounce values near Min
  RealType range = Max - Min;
  for( int i = 0; i < Array.getNumX(); i++ ) {
    for( int j = 0; j < Array.getNumY(); j++ ) {
      RealType tmp;
      const RealType value = Array.get(i,j);
      if(value < Min)
        tmp = Min;
      else {
        if(value > Max)
          tmp = Max;
        else
          tmp = value;
      }
      // tmp is scaled to a value in [1,LogarithmBase]
      tmp = (tmp - Min)/(1.0*range)*(LogarithmBase-1)+1;
      Array.set(i,j, ( log10(tmp)/log10(LogarithmBase)) );
    }
  }
}

template <typename RealType>
void joinTwoArrays2D( const qc::ScalarArray<RealType, qc::QC_2D> &Arg1, const qc::ScalarArray<RealType, qc::QC_2D> &Arg2, qc::ScalarArray<RealType, qc::QC_2D> &Dest, const int Spacing = 0 ) {
  if( Arg1.getNumY() != Arg2.getNumY() )
    throw aol::Exception( "Arg1.getNumY() != Arg2.getNumY() !", __FILE__, __LINE__);
  Dest.reallocate(Arg1.getNumX() + Spacing + Arg2.getNumX(), Arg1.getNumY());
  for( int j = 0; j < Arg1.getNumY(); j++ ) {
    for( int i = 0; i < Arg1.getNumX(); i++ ) {
      Dest.set( i, j, Arg1.get(i,j) );
    }
    for( int i = 0; i < Spacing; i++ ) {
      Dest.set( i + Arg1.getNumX(), j, 0.5 );
    }
    for( int i = 0; i < Arg2.getNumX(); i++ ) {
      Dest.set( i + Arg1.getNumX() + Spacing, j, Arg2.get(i,j) );
    }
  }
}

template <typename RealType>
void joinThreeArrays2D( const qc::ScalarArray<RealType, qc::QC_2D> &Arg1, const qc::ScalarArray<RealType, qc::QC_2D> &Arg2, const qc::ScalarArray<RealType, qc::QC_2D> &Arg3, qc::ScalarArray<RealType, qc::QC_2D> &Dest, const int Spacing = 0  ) {
  qc::ScalarArray<RealType, qc::QC_2D> temp;
  joinTwoArrays2D( Arg1, Arg2, temp, Spacing );
  joinTwoArrays2D( temp, Arg3, Dest, Spacing );
}

template <typename ConfiguratorType>
void resampleScalarArray( const qc::GridDefinition &Grid, const typename ConfiguratorType::ArrayType &Arg, typename ConfiguratorType::ArrayType &Dest ) {
  const typename ConfiguratorType::RealType h = Grid.H();
  typename ConfiguratorType::RealType scaling[3];
  scaling[0] = (Arg.getNumX() - 1) * h;
  scaling[1] = (Arg.getNumY() - 1) * h;
  scaling[2] = (Arg.getNumZ() - 1) * h;

  typename ConfiguratorType::VecType coords;
  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper(Grid);
  for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {
    for( int i = 0; i < ConfiguratorType::Dim; i++) {
      coords[i] = (*fnit)[i] * scaling[i];
    }
    Dest[ mapper.getGlobalIndex(*fnit) ] = Arg.interpolate( coords );
  }
}

#endif//__UTILITIES_H
