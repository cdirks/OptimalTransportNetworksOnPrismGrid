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

#include <FEOpInterface.h>
#include <deformations.h>
#include <hyperelastic.h>

#ifdef _OPENMP
#include <omp.h>
#endif

template<typename ConfiguratorType>
class AmbrosioTortorelliEdgeDetector :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType,AmbrosioTortorelliEdgeDetector<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _image;

public:
  AmbrosioTortorelliEdgeDetector ( const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &Image ) :
    aol::FELinScalarWeightedMassInterface<ConfiguratorType,AmbrosioTortorelliEdgeDetector<ConfiguratorType> >( Grid ),
    _image( Grid, Image ) {}

  RealType getCoeff ( const typename ConfiguratorType::ElementType &El,
                      int QuadPoint,
                      const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType imageGradient;
    _image.evaluateGradientAtQuadPoint( El, QuadPoint, imageGradient );
    return imageGradient.normSqr();
  }
};

template<typename ConfiguratorType>
class AmbrosioTortorelliEdgeDetectorMixedHessian :
  public aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType,AmbrosioTortorelliEdgeDetectorMixedHessian<ConfiguratorType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _image, _phaseField;

public:
  AmbrosioTortorelliEdgeDetectorMixedHessian ( const typename ConfiguratorType::InitType &Grid,
                                               const aol::Vector<RealType> &Image,
                                               const aol::Vector<RealType> &PhaseField,
                                               const bool Transpose) :
    aol::FELinVectorWeightedSemiDiffInterface<ConfiguratorType,AmbrosioTortorelliEdgeDetectorMixedHessian<ConfiguratorType> >( Grid, Transpose ),
    _image( Grid, Image ),
    _phaseField( Grid, PhaseField ) {}

  void getCoefficientVector ( const typename ConfiguratorType::ElementType &El,
                              int QuadPoint,
                              const typename ConfiguratorType::VecType &/*LocalCoord*/,
                              typename ConfiguratorType::VecType &Vector ) const {
    _image.evaluateGradientAtQuadPoint( El, QuadPoint, Vector );
    Vector *= 4. * _phaseField.evaluateAtQuadPoint( El, QuadPoint );
  }
};

template<typename ConfiguratorType>
class AmbrosioTortorelliSegmentationEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  AmbrosioTortorelliSegmentationEnergy( const typename ConfiguratorType::InitType &Grid,
                                        const aol::Vector<RealType> &Image,
                                        const RealType Alpha,
                                        const RealType Nu,
                                        const RealType Epsilon ) :
    _grid( Grid ),
    _massOp( Grid, aol::ASSEMBLED ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: smoothed image; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY );
    // fitting term
    aux1 -= Arg[0];
    _massOp.apply( aux1, aux2 );
    Dest += ( aux1 * aux2 ) * _alpha;
    // smoothing term
    AmbrosioTortorelliEdgeDetector<ConfiguratorType>( _grid, Arg[0] ).apply( Arg[1], aux2 );
    Dest += ( Arg[1] * aux2 );
    // perimeter
    _stiffOp.apply( Arg[1], aux2 );
    Dest += ( Arg[1] * aux2 ) * _nu * _epsilon / 2.;
    aux1 = Arg[1];
    aux1.addToAll( -1. );
    _massOp.apply( aux1, aux2 );
    Dest += ( aux1 * aux2 ) * _nu / _epsilon / 2.;
  }
};

template<typename ConfiguratorType>
class AmbrosioTortorelliSegmentationDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  AmbrosioTortorelliSegmentationDerivative( const typename ConfiguratorType::InitType &Grid,
                                            const aol::Vector<RealType> &Image,
                                            const RealType Alpha,
                                            const RealType Nu,
                                            const RealType Epsilon ) :
    _grid( Grid ),
    _massOp( Grid, aol::ASSEMBLED ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: smoothed image; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY );

    // derivative wrt smoothed image
    // fitting term
    aux1 -= Arg[0];
    _massOp.apply( aux1, aux2 );
    Dest[0].addMultiple( aux2, -2. * _alpha );
    // smoothing term
    aol::SquaredWeightStiffOp<ConfiguratorType>( _grid, Arg[1] ).apply( Arg[0], aux2 );
    Dest[0].addMultiple( aux2, 2. );

    // derivative wrt phase field
    // smoothing term
    AmbrosioTortorelliEdgeDetector<ConfiguratorType>( _grid, Arg[0] ).apply( Arg[1], aux2 );
    Dest[1].addMultiple( aux2, 2. );
    // perimeter
    _stiffOp.apply( Arg[1], aux2 );
    Dest[1].addMultiple( aux2, _nu * _epsilon );
    aux1 = Arg[1];
    aux1.addToAll( -1. );
    _massOp.apply( aux1, aux2 );
    Dest[1].addMultiple( aux2, _nu / _epsilon );
  }
};

template<typename ConfiguratorType, typename SubMatrixType>
class AmbrosioTortorelliSegmentationHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::BlockOp<SubMatrixType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  AmbrosioTortorelliSegmentationHessian ( const typename ConfiguratorType::InitType &Grid,
                                          const aol::Vector<RealType> &Image,
                                          const RealType Alpha,
                                          const RealType Nu,
                                          const RealType Epsilon ) :
    _grid( Grid ),
    _massOp( Grid, aol::ASSEMBLED ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: smoothed image; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::BlockOp<SubMatrixType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY );

    // Hessian wrt smoothed image
    // fitting term
    _massOp.assembleAddMatrix( Dest.getReference( 0, 0 ), 2. * _alpha );
    // smoothing term
    aol::SquaredWeightStiffOp<ConfiguratorType>( _grid, Arg[1] ).assembleAddMatrix( Dest.getReference( 0, 0 ), 2. );

    // Hessian wrt phase field
    // smoothing term
    AmbrosioTortorelliEdgeDetector<ConfiguratorType>( _grid, Arg[0] ).assembleAddMatrix( Dest.getReference( 1, 1 ), 2. );
    // perimeter
    _stiffOp.assembleAddMatrix( Dest.getReference( 0, 0 ), _nu * _epsilon );
    _massOp.assembleAddMatrix( Dest.getReference( 0, 0 ), _nu / _epsilon );

    // mixed derivative (smoothing term)
    AmbrosioTortorelliEdgeDetectorMixedHessian<ConfiguratorType>( _grid, Arg[0], Arg[1], false ).assembleAddMatrix( Dest.getReference( 0, 1 ) );
    AmbrosioTortorelliEdgeDetectorMixedHessian<ConfiguratorType>( _grid, Arg[0], Arg[1], true ).assembleAddMatrix( Dest.getReference( 1, 0 ) );
  }
};

template<typename ConfiguratorType>
class ModicaMortolaSegmentationEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  ModicaMortolaSegmentationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                    const aol::Vector<RealType> &Image,
                                    const RealType Alpha,
                                    const RealType Nu,
                                    const RealType Epsilon ) :
    _grid( Grid ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: both gray values; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY ), aux3( Arg[1], aol::DEEP_COPY );
    // fitting term
    aux1.addToAll( -Arg[0][0] );
    aux3.addToAll( 1. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux3 ).apply( aux1, aux2 );
    Dest += ( aux1 * aux2 ) * _alpha / 4.;
    aux1.addToAll( Arg[0][0] - Arg[0][1] );
    aux3.addToAll( -2. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux3 ).apply( aux1, aux2 );
    Dest += ( aux1 * aux2 ) * _alpha / 4.;
    // perimeter
    _stiffOp.apply( Arg[1], aux2 );
    Dest += ( Arg[1] * aux2 ) * _nu * _epsilon / 2.;
    aux1 = Arg[1];
    aux1.addToAll( 1. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux3 ).apply( aux1, aux2 );
    Dest += ( aux1 * aux2 ) * aol::Sqr( .75 ) * _nu / _epsilon / 2.;
  }
};

template<typename ConfiguratorType>
class ModicaMortolaSegmentationDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  ModicaMortolaSegmentationDerivative ( const typename ConfiguratorType::InitType &Grid,
                                        const aol::Vector<RealType> &Image,
                                        const RealType Alpha,
                                        const RealType Nu,
                                        const RealType Epsilon ) :
    _grid( Grid ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _massOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: both gray values; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY ), aux3( Arg[1], aol::DEEP_COPY ),
                          ones( Arg[1], aol::STRUCT_COPY );
    ones.setAll( 1. );

    // derivative wrt gray values
    // fitting term
    aux1.addToAll( -Arg[0][0] );
    aux3.addToAll( 1. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux3 ).apply( aux1, aux2 );
    Dest[0][0] -= ( aux2 * ones ) * _alpha / 2.;
    aux1.addToAll( Arg[0][0] - Arg[0][1] );
    aux3.addToAll( -2. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux3 ).apply( aux1, aux2 );
    Dest[0][1] -= ( aux2 * ones ) * _alpha / 2.;

    // derivative wrt phase field
    // fitting term
    aux1 = _image;
    aux3 = Arg[1];
    aux1.addToAll( -Arg[0][0] );
    aux3.addToAll( 1. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux1 ).apply( aux3, aux2 );
    Dest[1].addMultiple( aux2, _alpha / 2. );
    aux1.addToAll( Arg[0][0] - Arg[0][1] );
    aux3.addToAll( -2. );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux1 ).apply( aux3, aux2 );
    Dest[1].addMultiple( aux2, _alpha / 2. );
    // perimeter
    _stiffOp.apply( Arg[1], aux2 );
    Dest[1].addMultiple( aux2, _nu * _epsilon );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, Arg[1] ).apply( Arg[1], aux2 );
    _massOp.apply( Arg[1], aux1 );
    aux2 -= aux1;
    Dest[1].addMultiple( aux2, 2.25 * _nu / _epsilon / 2. );
  }
};

template<typename ConfiguratorType,typename SubMatrixType>
class ModicaMortolaSegmentationHessian :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::BlockOp<SubMatrixType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::Vector<RealType> &_image;
  const RealType _alpha, _nu, _epsilon;

public:
  ModicaMortolaSegmentationHessian( const typename ConfiguratorType::InitType &Grid,
                                    const aol::Vector<RealType> &Image,
                                    const RealType Alpha,
                                    const RealType Nu,
                                    const RealType Epsilon ) :
    _grid( Grid ),
    _stiffOp( Grid, aol::ASSEMBLED ),
    _massOp( Grid, aol::ASSEMBLED ),
    _image( Image ),
    _alpha( Alpha ),
    _nu( Nu ),
    _epsilon( Epsilon ) {}

  //! 1st component of Arg: both gray values; 2nd component of Arg: phase field
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::BlockOp<SubMatrixType> &Dest ) const {
    aol::Vector<RealType> aux1( _image, aol::DEEP_COPY ), aux2( _image, aol::STRUCT_COPY ), aux3( Arg[1], aol::DEEP_COPY );

    // Hessian wrt gray values
    // fitting term
    aux3.addToAll( 1. );
    _massOp.apply( aux3, aux2 );
    Dest.getReference( 0, 0 ).add( 0, 0, ( aux3 * aux2 ) / 2. );
    aux3.addToAll( -2. );
    _massOp.apply( aux3, aux2 );
    Dest.getReference( 0, 0 ).add( 1, 1, ( aux3 * aux2 ) / 2. );

    // Hessian wrt phase field
    // fitting term
    aux1.addToAll( -Arg[0][0] );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux1 ).assembleAddMatrix( Dest.getReference( 1, 1 ), .5 );
    aux1.addToAll( Arg[0][0] - Arg[0][1] );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, aux1 ).assembleAddMatrix( Dest.getReference( 1, 1 ), .5 );
    // perimeter
    _stiffOp.assembleAddMatrix( Dest.getReference( 1, 1 ), _nu * _epsilon );
    _massOp.assembleAddMatrix( Dest.getReference( 1, 1 ), -1.125 * _nu / _epsilon );
    aol::SquaredWeightMassOp<ConfiguratorType>( _grid, Arg[1] ).assembleAddMatrix( Dest.getReference( 1, 1 ), 3.375 * _nu / _epsilon );

    // mixed derivative
    aux1 = _image;
    aux3 = Arg[1];
    aux1.addToAll( -Arg[0][0] );
    aux3.addToAll( 1. );
    aol::WeightedMassOp<ConfiguratorType>( _grid, aux1 ).apply( aux3, aux2 );
    for ( int i = 0; i < aux2.size(); i++ ) {
      Dest.getReference( 0, 1 ).add( 0, i, aux2[i] );
      Dest.getReference( 1, 0 ).add( i, 0, aux2[i] );
    }
    aux1.addToAll( Arg[0][0] - Arg[0][1] );
    aux3.addToAll( -2. );
    aol::WeightedMassOp<ConfiguratorType>( _grid, aux1 ).apply( aux3, aux2 );
    for ( int i = 0; i < aux2.size(); i++ ) {
      Dest.getReference( 0, 1 ).add( 1, i, aux2[i] );
      Dest.getReference( 1, 0 ).add( i, 1, aux2[i] );
    }
  }
};

template<typename ConfiguratorType>
class MismatchEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;

public:
  explicit MismatchEnergy ( const typename ConfiguratorType::InitType &Grid ) :
    _grid( Grid ) {}

  //! 1st, 2nd component of Arg: input images; remaining components of Arg: displacement
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> disp;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      disp.appendReference( Arg[2+i] );
    qc::MismatchEnergyWRTDisp<ConfiguratorType>( _grid, Arg[1], Arg[0] ).applyAdd( disp, Dest );
  }
};

template<typename ConfiguratorType>
class MismatchDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;

public:
  explicit MismatchDerivative ( const typename ConfiguratorType::InitType &Grid ) :
    _grid( Grid ),
    _massOp( Grid, aol::ASSEMBLED ) {}

  //! 1st, 2nd component of Arg: input images; remaining components of Arg: displacement
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> disp, dispVar;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      disp.appendReference( Arg[2+i] );
      dispVar.appendReference( Dest[2+i] );
    }
    // derivative wrt displacement
    qc::MismatchEnergyVariationWRTDisp<ConfiguratorType>( _grid, Arg[1], Arg[0] ).applyAdd( disp, dispVar );
    // derivative wrt non-deformed image
    aol::Vector<RealType> aux1( Arg[1], aol::STRUCT_COPY );
    _massOp.apply( Arg[1], aux1 );
    Dest[1].addMultiple( aux1, 2. );
    qc::DeformMassOp<ConfiguratorType>( _grid, disp, aol::ONTHEFLY, qc::RIGHT ).apply( Arg[0], aux1 );
    Dest[1].addMultiple( aux1, -2. );
    // derivative wrt deformed image
    qc::DeformMassOp<ConfiguratorType>( _grid, disp, aol::ONTHEFLY, qc::BOTH ).apply( Arg[0], aux1 );
    Dest[0].addMultiple( aux1, 2. );
    qc::DeformMassOp<ConfiguratorType>( _grid, disp, aol::ONTHEFLY, qc::LEFT ).apply( Arg[1], aux1 );
    Dest[0].addMultiple( aux1, -2. );
  }
};

template<typename ConfiguratorType,typename SegmentationEnergyType>
class JointSegmentationRegistrationEnergy :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const SegmentationEnergyType _segmEnergy0, _segmEnergy1;
  const MismatchEnergy<ConfiguratorType> _mismatchEnergy;
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _elastDensity;
  const qc::HyperelasticEnergy<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > _elastEnergy;
  const RealType _a;

public:
  JointSegmentationRegistrationEnergy ( const typename ConfiguratorType::InitType &Grid,
                                        const aol::Vector<RealType> &Image0,
                                        const aol::Vector<RealType> &Image1,
                                        const RealType Nu,
                                        const RealType Epsilon,
                                        const RealType Alpha,
                                        const RealType A,
                                        const RealType LengthChangeWeight,
                                        const RealType AreaChangeWeight,
                                        const RealType VolumeChangeWeight ) :
    _grid( Grid ),
    _segmEnergy0( _grid, Image0, Alpha, Nu, Epsilon ),
    _segmEnergy1( _grid, Image1, Alpha, Nu, Epsilon ),
    _mismatchEnergy( Grid ),
    _elastDensity( LengthChangeWeight, AreaChangeWeight, VolumeChangeWeight ),
    _elastEnergy( _grid, _elastDensity ),
    _a( A ) {}

  //! 1st, 2nd component of Arg: smoothed images; 3rd, 4th component: phase fields; remaining components: displacement
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> arg1, arg2, arg3, disp;
    // segmentation energy
    arg1.appendReference( Arg[0] );
    arg1.appendReference( Arg[2] );
    _segmEnergy0.applyAdd( arg1, Dest );
    arg2.appendReference( Arg[1] );
    arg2.appendReference( Arg[3] );
    _segmEnergy1.applyAdd( arg2, Dest );
    // deformation regularization
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      disp.appendReference( Arg[4+i] );
    _elastEnergy.applyAdd( disp, Dest );
    // matching energy
    arg3.appendReference( Arg[2] );
    arg3.appendReference( Arg[3] );
    arg3.appendReference( disp );
    aol::Scalar<RealType> dest;
    _mismatchEnergy.applyAdd( arg3, dest );
    Dest.addMultiple( dest, _a );
  }
};

template<typename ConfiguratorType, typename SegmentationDerivativeType>
class JointSegmentationRegistrationDerivative :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const SegmentationDerivativeType _segmDeriv0, _segmDeriv1;
  const MismatchDerivative<ConfiguratorType> _mismatchDeriv;
  const qc::HyperelasticEnergyDensityDefault<ConfiguratorType> _elastDensity;
  const qc::HyperelasticGradient<ConfiguratorType,qc::HyperelasticEnergyDensityDefault<ConfiguratorType> > _elastDeriv;
  const RealType _a;
  const int _mode; // compute gradient only wrt: 2 - displacement, 1 - segmentation variables, 0 - both

public:
  JointSegmentationRegistrationDerivative ( const typename ConfiguratorType::InitType &Grid,
                                            const aol::Vector<RealType> &Image0,
                                            const aol::Vector<RealType> &Image1,
                                            const RealType Nu,
                                            const RealType Epsilon,
                                            const RealType Alpha,
                                            const RealType A,
                                            const RealType LengthChangeWeight,
                                            const RealType AreaChangeWeight,
                                            const RealType VolumeChangeWeight,
                                            const int Mode = 0 ) :
    _grid( Grid ),
    _segmDeriv0( _grid, Image0, Alpha, Nu, Epsilon ),
    _segmDeriv1( _grid, Image1, Alpha, Nu, Epsilon ),
    _mismatchDeriv( Grid ),
    _elastDensity( LengthChangeWeight, AreaChangeWeight, VolumeChangeWeight ),
    _elastDeriv( _grid, _elastDensity ),
    _a( A ),
    _mode( Mode ) {}

  //! 1st, 2nd component of Arg: smoothed images; 3rd, 4th component: phase fields; remaining components: displacement
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> arg1, dest1, arg2, dest2, arg3, dest3, disp, destDisp;
    // segmentation energy
    arg1.appendReference( Arg[0] );
    arg1.appendReference( Arg[2] );
    dest1.appendReference( Dest[0] );
    dest1.appendReference( Dest[2] );
    _segmDeriv0.applyAdd( arg1, dest1 );
    arg2.appendReference( Arg[1] );
    arg2.appendReference( Arg[3] );
    dest2.appendReference( Dest[1] );
    dest2.appendReference( Dest[3] );
    _segmDeriv1.applyAdd( arg2, dest2 );
    // deformation regularization
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      disp.appendReference( Arg[4+i] );
      destDisp.appendReference( Dest[4+i] );
    }
    _elastDeriv.applyAdd( disp, destDisp );
    // matching energy
    arg3.appendReference( Arg[2] );
    arg3.appendReference( Arg[3] );
    arg3.appendReference( disp );
    aol::MultiVector<RealType> dest( arg3, aol::STRUCT_COPY );
    _mismatchDeriv.applyAdd( arg3, dest );
    for ( int i = 0; i < arg3.numComponents(); i++ )
      Dest[2+i].addMultiple( dest[i], _a );
    // delete part of the derivative for alternating descent
    if ( _mode == 2 )
      for ( int i = 0; i < 4; i++ )
        Dest[i].setZero();
    if ( _mode == 1 )
      for ( int i = 4; i < 4 + ConfiguratorType::Dim; i++ )
        Dest[i].setZero();
  }
};




typedef double RealType;
const qc::Dimension Dim = qc::QC_2D;
typedef qc::QuocConfiguratorTraitMultiLin<RealType, Dim, aol::Quadrature2D<RealType,aol::SimpsonQuadrature<RealType> > > ConfiguratorType;

int main( int /*argc*/, char *argv[] ) {
#ifdef _OPENMP
  omp_set_num_threads( 4 );
#endif

  try {
    typedef ConfiguratorType::ArrayType ArrayType;

    // read in parameters
    const aol::ParameterParser parameterParser( argv[1] );
    char saveDir[1024];
    parameterParser.getString( "destDirectory", saveDir );
    const bool flagAT = ( parameterParser.getInt( "segmMethod" ) );
    aol::Vector<int> level, numIter, alternationMode;
    parameterParser.getIntVec( "level", level );
    parameterParser.getIntVec( "numIter", numIter );
    if ( parameterParser.hasVariable( "alternationMode" ) )
      parameterParser.getIntVec( "alternationMode", alternationMode );
    else
      alternationMode.resize( level.size() );
    aol::Vector<RealType> a, nu, alpha, epsilon, elastWeight1, elastWeight2, elastWeight3;
    parameterParser.getRealVec( "a", a );
    parameterParser.getRealVec( "nu", nu );
    parameterParser.getRealVec( "alpha", alpha );
    parameterParser.getRealVec( "epsilon", epsilon );
    parameterParser.getRealVec( "elastWeight1", elastWeight1 );
    parameterParser.getRealVec( "elastWeight2", elastWeight2 );
    parameterParser.getRealVec( "elastWeight3", elastWeight3 );
    const bool onlyProlongateDisplacement = parameterParser.hasVariable( "onlyProlongateDisplacement" );

    // read in images
    const ArrayType img0( parameterParser.getString( "inputImage0" ).c_str() ), img1( parameterParser.getString( "inputImage1" ).c_str() );

    // initialize variables
    const qc::GridSize<Dim> gridSize( img0 );
    const ConfiguratorType::InitType grid( gridSize );
    const int fineLevel = grid.getGridDepth();
    qc::MultiDimMultilevelArray<RealType,ArrayType> inputImages( fineLevel, Dim, 2 ),
                                                    phaseFields( fineLevel, Dim, 2 ),
                                                    smoothImages( fineLevel, Dim, 2 ),
                                                    displacement( fineLevel, Dim, Dim );
    inputImages.getArray( 0 ) = img0;
    inputImages.getArray( 1 ) = img1;
    smoothImages.getArray( 0 ) = img0;
    smoothImages.getArray( 1 ) = img1;
    if ( flagAT )
      for ( int i = 0; i < 2; i++ )
        phaseFields.getArray( i ).setAll( 1. );
    aol::MultiVector<RealType> grayVals( 2, 2 );
    grayVals[0][1] = 256;
    grayVals[1][1] = 256;
    inputImages.levRestrict( 2, fineLevel );
    phaseFields.levRestrict( 2, fineLevel );
    smoothImages.levRestrict( 2, fineLevel );
    displacement.levRestrict( 2, fineLevel );

    // perform multiscale method
    for ( int iter = 0; iter < level.size(); iter++ ) {
      // set current level
      inputImages.setCurLevel( level[iter] );
      phaseFields.setCurLevel( level[iter] );
      smoothImages.setCurLevel( level[iter] );
      displacement.setCurLevel( level[iter] );
      const ConfiguratorType::InitType grid( level[iter], Dim );
      // perform minimization
      aol::MultiVector<RealType> arg;
      if ( flagAT )
        smoothImages.appendReferencesTo( arg );
      else
        arg.appendReference( grayVals );
      phaseFields.appendReferencesTo( arg );
      displacement.appendReferencesTo( arg );
      aol::MultiVector<RealType> initArg( arg, aol::DEEP_COPY );
      if ( flagAT ) {               // Ambrosio-Tortorelli method
        JointSegmentationRegistrationEnergy<ConfiguratorType,AmbrosioTortorelliSegmentationEnergy<ConfiguratorType> >
          e( grid, inputImages.getArray( 0 ), inputImages.getArray( 1 ), nu[iter], epsilon[iter], alpha[iter], a[iter], elastWeight1[iter], elastWeight2[iter], elastWeight3[iter] );
        JointSegmentationRegistrationDerivative<ConfiguratorType,AmbrosioTortorelliSegmentationDerivative<ConfiguratorType> >
          de( grid, inputImages.getArray( 0 ), inputImages.getArray( 1 ), nu[iter], epsilon[iter], alpha[iter], a[iter], elastWeight1[iter], elastWeight2[iter], elastWeight3[iter], alternationMode[iter] );
        aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimizer( e, de, numIter[iter], 1.e-4, 50, false );
        minimizer.setTimestepController( aol::NewtonInfo<RealType>::ARMIJO );
        minimizer.apply( initArg, arg );
      } else {                      // Modica-Mortola method
        JointSegmentationRegistrationEnergy<ConfiguratorType,ModicaMortolaSegmentationEnergy<ConfiguratorType> >
          e( grid, inputImages.getArray( 0 ), inputImages.getArray( 1 ), nu[iter], epsilon[iter], alpha[iter], a[iter], elastWeight1[iter], elastWeight2[iter], elastWeight3[iter] );
        JointSegmentationRegistrationDerivative<ConfiguratorType,ModicaMortolaSegmentationDerivative<ConfiguratorType> >
          de( grid, inputImages.getArray( 0 ), inputImages.getArray( 1 ), nu[iter], epsilon[iter], alpha[iter], a[iter], elastWeight1[iter], elastWeight2[iter], elastWeight3[iter], alternationMode[iter] );
        aol::QuasiNewtonIteration<RealType,aol::MultiVector<RealType>,aol::MultiVector<RealType> > minimizer( e, de, numIter[iter], 1.e-4, 50, false );
        minimizer.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );
        minimizer.apply( initArg, arg );
      }
      // save result
      char filename[1024];
      aol::MultiVector<RealType> disp;
      displacement.appendReferencesTo( disp );
      for ( int i = 0; i < Dim; i++ ) {
        sprintf( filename, "%sdisp%d.bz2", saveDir, i );
        displacement.getArray( i ).save( filename, qc::PGM_DOUBLE_BINARY );
      }
      ArrayType defImg( inputImages.getArray( 0 ), aol::STRUCT_COPY );
      qc::DeformImage<ConfiguratorType>( inputImages.getArray( 0 ), grid, defImg, disp );
      sprintf( filename, "%sdefImg.pgm", saveDir );
      defImg.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
      if ( flagAT ) {               // Ambrosio-Tortorelli method
        for ( int i = 0; i < 2; i++ ) {
          sprintf( filename, "%simg%d.pgm", saveDir, i );
          smoothImages.getArray( i ).save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
          ArrayType phf( phaseFields.getArray( i ), aol::DEEP_COPY );
          phf *= 256;
          sprintf( filename, "%sphf%d.pgm", saveDir, i );
          phf.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
        }
      } else {                      // Modica-Mortola method
        for ( int i = 0; i < 2; i++ ) {
          ArrayType colPhf( phaseFields.getArray( i ), aol::DEEP_COPY );
          colPhf.addToAll( -1. );
          colPhf *= -.5 * ( grayVals[i][1] - grayVals[i][0] );
          colPhf.addToAll( grayVals[i][0] );
          sprintf( filename, "%sphf%d.pgm", saveDir, i );
          colPhf.save( filename, qc::PGM_UNSIGNED_CHAR_BINARY );
        }
      }
      // prolongate result
      if ( !onlyProlongateDisplacement ) {
        phaseFields.levProlongate();
        smoothImages.levProlongate();
      }
      displacement.levProlongate();
    }
  }
  catch ( aol::Exception &el ) {
    // print out error message
    el.dump();
  }

  // pause before closing the program (to let the user first read the output)
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_SUCCESS;
}
