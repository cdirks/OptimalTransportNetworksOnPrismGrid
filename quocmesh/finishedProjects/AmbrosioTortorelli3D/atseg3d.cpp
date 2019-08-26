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
 * \brief  The class implements Ambrosio-Tortorelli based segmentation
 * \author Berkels, Han
 */

// quocmesh defined
#include <parameterParser.h>
#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <preconditioner.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <convolution.h>
#include <AmbrosioTortorelli.h>
#include <gradientDescent.h>
#include <Newton.h>

// user defined
#include "registration.h"

// type of image pixel
typedef double RType;

// configuration
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;

template <typename ConfiguratorType>
class DerivativeOfL2ConvolutionFidelityTermWithRegardToU : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  const typename ConfiguratorType::InitType &_grid;
  qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> _massMatrix;
  const aol::MassOp<ConfiguratorType> _mass;
  const qc::Convolution<qc::QC_2D> &_conv;
  qc::ScalarArray<RealType, qc::QC_2D> _convolvedKernel;
public:
  DerivativeOfL2ConvolutionFidelityTermWithRegardToU( const typename ConfiguratorType::InitType &Initializer,
                                                      const ArrayType &Kernel,
                                                      const ArrayType &KernelMinus,
                                                      qc::Convolution<qc::QC_2D> &Conv )
  : _grid( Initializer ),
    _massMatrix( _grid ),
    _mass( _grid ),
    _conv( Conv ),
    _convolvedKernel( _grid )
  {
    _conv.convolve( Kernel, KernelMinus, _convolvedKernel );
    _mass.assembleAddMatrix( _massMatrix );
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    ArrayType tmp( _grid );
    ArrayType argArray ( Arg, _grid );
    _conv.convolve( argArray, _convolvedKernel, tmp );
    _massMatrix.apply( tmp, Dest );
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};


/**@brief The operator computing the the energy respect ot phase field function v
*/
template <typename ConfiguratorType, typename AnisotropyType>
class AnisoATEnergyWithRegardToV: public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> >
{
  /// type define
  typedef typename ConfiguratorType::RealType RealType;

  /// internal grid
  const typename ConfiguratorType::InitType &_grid;

  /// The type of anisotroy
  const AnisotropyType &_anisotropy;

  /// the matrix
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &_matLinearPart;

  /// mass operator
  const aol::MassOp<ConfiguratorType> _mass;

  /// The parameter alpha and epsilon;
  const RealType _alpha, _epsilon;

  /// The constant of mass
  RealType _massconstant;

  /// temp vector
  mutable aol::Vector<RealType> _tmp,_tmp2;
public:

  /// Initialization constructor
  AnisoATEnergyWithRegardToV(const typename ConfiguratorType::InitType &Initializer,
                             const AnisotropyType &Anisotropy,
                             const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &matLinearPart,
                             const RealType Alpha,
                             const RealType Epsilon):
  _grid(Initializer),
  _anisotropy(Anisotropy),
  _matLinearPart(matLinearPart),
  _mass(Initializer),
  _alpha(Alpha),
  _epsilon(Epsilon),
  _tmp(_grid.getNumberOfNodes()),
  _tmp2(_grid.getNumberOfNodes())
  {
    _tmp.setAll(1.); // massconstant = alpha/(8.*epsilon)*M*1vec*1vec
    _mass.apply(_tmp, _tmp2);
    _massconstant = _alpha/(8.*_epsilon)*(_tmp* _tmp2);
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest) const{
    RealType Energy = _massconstant; // Energy = (alpha/(8*epsilon))*M*1vec*1vec
    _matLinearPart.apply(Arg, _tmp); //tmp = (alpha/(4*epsilon))*M*V + M[GradU2]*V
    Energy += 0.5*(Arg * _tmp); // Energy = (alpha/(8*epsilon))*M*V*V + 1/2*M[GradU2]*V*V

    _tmp.setAll(1.); // Energy += -alpha/(4*epsilon)*M*V*1vec
    _mass.apply(Arg, _tmp2);
    Energy += _alpha/(-4.*_epsilon)*(_tmp2 * _tmp);

    qc::ATAnisoGradSqrOp<ConfiguratorType, AnisotropyType> atAnisoGradSqrOp( _grid, Arg, _anisotropy);
    atAnisoGradSqrOp.apply(_tmp, _tmp2); //Argument tmp is ignored
    Energy += (_tmp * _tmp2);
    Dest[0] = Energy;
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/) const{
  }
};
/**@brief The operator computing the first variation of energy respect to the phase field function v
*/
template <typename ConfiguratorType, typename AnisotropyType>
class VariationOfAnisoATEnergyWithRegardToV : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const AnisotropyType &_anisotropy;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &_matLinearPart;
  const aol::MassOp<ConfiguratorType> _mass;
  const RealType _alpha, _epsilon;
public:
  VariationOfAnisoATEnergyWithRegardToV(const typename ConfiguratorType::InitType &Initializer, const AnisotropyType &Anisotropy, const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &matLinearPart, const RealType Alpha, const RealType Epsilon):
  _grid(Initializer), _anisotropy(Anisotropy), _matLinearPart(matLinearPart), _mass(Initializer), _alpha(Alpha), _epsilon(Epsilon){
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    aol::Vector<RealType> tmp(Dest.getSize());

    tmp.setAll( -1.*_alpha / ( _epsilon * 4. ) ); // Dest = -(alpha/(4*epsilon))*M(1,...,1)
    _mass.apply( tmp, Dest );

    _matLinearPart.applyAdd( Arg, Dest); // Dest += M[GradU2]*v + (alpha/(4*epsilon))*M*v

    Dest *= 1./(_alpha*_epsilon);  // Dest += alpha*epsilon \int b(v)\cdot \grad phi_i
    qc::ATAnisotropyDiffOp<ConfiguratorType, AnisotropyType> atAnisotropyDiffOp(_grid, Arg, _anisotropy);
    atAnisotropyDiffOp.applyAdd(tmp, Dest); //Argument tmp is not used by atAnisotropyDiffOpForU
    Dest *= (_alpha*_epsilon);
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};
/**@brief The operator computing the first variation of energy respect to the function U
*/
template <typename ConfiguratorType, typename AnisotropyType>
class VariationOfAnisoATEnergyWithRegardToU : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const AnisotropyType &_anisotropy;
  const aol::Vector<RealType> &_v;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &_matLinearPart;
  const aol::MassOp<ConfiguratorType> _mass;
  aol::Vector<RealType> _u0constantVec; // stores -beta*M*u_0
public:
  VariationOfAnisoATEnergyWithRegardToU(const typename ConfiguratorType::InitType &Initializer,
                                        const AnisotropyType &Anisotropy,
                                        const aol::Vector<RealType> &U0,
                                        const aol::Vector<RealType> &V,
                                        const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &matLinearPart,
                                        const RealType Beta):
  _grid(Initializer), _anisotropy(Anisotropy), _v(V), _matLinearPart(matLinearPart), _mass(Initializer), _u0constantVec(Initializer.getNumberOfNodes()){
    _mass.apply(U0, _u0constantVec);
    _u0constantVec *= (-1.*Beta);
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    Dest = _u0constantVec; // Dest=-beta*M*u_0
    _matLinearPart.applyAdd(Arg, Dest); // Dest += beta*M*u

    qc::ATAnisotropyDiffOpForU<ConfiguratorType, AnisotropyType> atAnisotropyDiffOpForU(_grid, _v, _anisotropy);
    atAnisotropyDiffOpForU.applyAdd(Arg, Dest);
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};
/**@brief The operator computing the second variation of energy respect to the phase field function v
*/
template <typename ConfiguratorType, typename AnisotropyType>
class SecondVariationOfAnisoATEnergyWithRegardToV : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, qc::FastUniformGridMatrix<typename ConfiguratorType::RealType,ConfiguratorType::Dim> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const AnisotropyType &_anisotropy;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &_matLinearPart;
  const RealType _alpha, _epsilon;
public:
  SecondVariationOfAnisoATEnergyWithRegardToV(const typename ConfiguratorType::InitType &Initializer, const AnisotropyType &Anisotropy, const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &matLinearPart, const RealType Alpha, const RealType Epsilon):
  _grid(Initializer), _anisotropy(Anisotropy), _matLinearPart(matLinearPart), _alpha(Alpha), _epsilon(Epsilon){
  }
  virtual void apply( const aol::Vector<RealType> &Arg, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &Dest ) const{
    Dest.setZero(); //Dest = alpha*epsilon L[A(v)]
    qc::ATAnisotropyStiffOp<ConfiguratorType, AnisotropyType> atAnisotropyStiffOp(_grid, Arg, _anisotropy);
    atAnisotropyStiffOp.assembleAddMatrix(Dest);
    Dest.scale(_alpha*_epsilon);
    Dest += _matLinearPart; //Dest +=  M[GradU2] + (alpha/(4*epsilon))*M

  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};
/**@brief The operator computing the second variation of energy with anisotropy respect to the function U
*/
template <typename ConfiguratorType, typename AnisotropyType>
class SecondVariationOfAnisoATEnergyWithRegardToU : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, qc::FastUniformGridMatrix<typename ConfiguratorType::RealType,ConfiguratorType::Dim> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const AnisotropyType &_anisotropy;
  const aol::Vector<RealType> &_v;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &_matLinearPart;
public:
  SecondVariationOfAnisoATEnergyWithRegardToU(const typename ConfiguratorType::InitType &Initializer,
                                              const AnisotropyType &Anisotropy,
                                              const aol::Vector<RealType> &V,
                                              const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &matLinearPart):
  _grid(Initializer), _anisotropy(Anisotropy), _v(V), _matLinearPart(matLinearPart){
  }
  virtual void apply( const aol::Vector<RealType> &Arg, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &Dest ) const{
    Dest.setZero(); //Dest = L[v^2*A(u)]
    qc::ATAnisotropyStiffOpForU<ConfiguratorType, AnisotropyType> atAnisotropyStiffOpForU(_grid, Arg, _v, _anisotropy);
    atAnisotropyStiffOpForU.assembleAddMatrix(Dest);
    Dest += _matLinearPart; //Dest += beta*M
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};
/**@brief The operator computing the first variation of energy with Double well v function
*/
template <typename ConfiguratorType>
class VariationOfDoubleWellATEnergy : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> _matLinearPart;
  const RealType _alpha, _epsilon;
public:
  VariationOfDoubleWellATEnergy(const typename ConfiguratorType::InitType &Initializer, const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> matLinearPart, const RealType Alpha, const RealType Epsilon):
  _grid(Initializer), _matLinearPart(matLinearPart), _alpha(Alpha), _epsilon(Epsilon){
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    qc::ATMassV2Op<ConfiguratorType> atMassV2Op( _grid, Arg); // Dest = alpha/(2*epsilon)M[v^2]*v
    atMassV2Op.apply( Arg, Dest);
    Dest *= _alpha/(2.*_epsilon);

    _matLinearPart.applyAdd( Arg, Dest); // Dest += M[GradU2]*v - (alpha/(2*epsilon))*M*v + L*v
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};
/**@brief The operator computing the second variation of double well energy respect ot phase field function v (theor)
*/
template <typename ConfiguratorType>
class SecondVariationOfDoubleWellATEnergy : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, qc::FastUniformGridMatrix<typename ConfiguratorType::RealType,ConfiguratorType::Dim> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> _matLinearPart;
  const RealType _alpha, _epsilon;
public:
  SecondVariationOfDoubleWellATEnergy(const typename ConfiguratorType::InitType &Initializer, const qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> matLinearPart, const RealType Alpha, const RealType Epsilon):
  _grid(Initializer), _matLinearPart(matLinearPart), _alpha(Alpha), _epsilon(Epsilon){
  }
  virtual void apply( const aol::Vector<RealType> &Arg, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &Dest ) const{
    Dest.setZero(); //Dest = (3*alpha)/(2*epsilon)M[v^2]
    qc::ATMassV2Op<ConfiguratorType> atMassV2Op( _grid, Arg);
    atMassV2Op.assembleAddMatrix(Dest);
    Dest.scale( (3.*_alpha)/(2.*_epsilon) );
    Dest += _matLinearPart; //Dest +=  M[GradU2] - (alpha/(2*epsilon))*M + L
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

/**@brief Segmentation class for single level*/
template <typename ConfiguratorType>
class ExtendedATSegmentation : public qc::ATSegmentation<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

public:
  ExtendedATSegmentation ( const RealType Epsilon,
                           const RealType Alpha,
                           const RealType Beta,
                           const qc::Array<RealType> &U0 )
    : qc::ATSegmentation<ConfiguratorType> ( Epsilon, Alpha, Beta, U0 ) { }

  void updateUwithAnisotropy ( aol::Vector<RealType> &U,
                              const aol::Vector<RealType> &V,
                              const int InnerIterations = 0,
                              const int OuterIterations = 0) {
    // use Newton Iteration to solve F(v)=0

    cerr << "Solving for u" << endl;
    this->_mat.setZero(); //mat = beta*M
    this->_mass.assembleAddMatrix( this->_mat );
    this->_mat.scale( this->_beta );

    const RealType delta = 0.1; //sqrt(0.1);
    //typedef qc::ATAnisotropy2Norm<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta , 1., 1.);
    typedef qc::ATAnisotropy1Norm<ConfiguratorType> AnisoType;
    AnisoType anisotropy( delta );
    //typedef qc::ATAnisotropyOktaeder<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta );
    //AnisotropyShapeWriter<RealType, AnisoType> anisotropyShapeWriter(anisotropy);
    //anisotropyShapeWriter.writeWulffShape( "wulffshape.pgm" );
    //anisotropyShapeWriter.writeFrankDiagram( "frankdiagram.pgm" );

    VariationOfAnisoATEnergyWithRegardToU<ConfiguratorType, AnisoType> F(this->_grid, anisotropy, this->_u0, V, this->_mat, this->_beta);
    SecondVariationOfAnisoATEnergyWithRegardToU<ConfiguratorType, AnisoType> DF(this->_grid, anisotropy, V, this->_mat);

    aol::NewtonIteration<ConfiguratorType> newton_solver(this->_grid, F, DF, 2, InnerIterations, OuterIterations);
    newton_solver.apply(U, this->_rhs);
    U = this->_rhs;
  }

  void updateVwithAnisotropyOnlyInV(const aol::Vector<RealType> &U,
                             aol::Vector<RealType> &V,
                             const int InnerIterations = 0,
                             const int OuterIterations = 0){
    // use Newton Iteration to solve F(v)=0

    cerr << "Solving for v" << endl;
    this->_mat.setZero(); //mat = (alpha/(4*epsilon))*M
    this->_mass.assembleAddMatrix( this->_mat );
    this->_mat.scale( this->_alpha / ( 4. * this->_epsilon) );

    qc::ATMassGradU2Op<ConfiguratorType> atMassGradU2( this->_grid, U );
    atMassGradU2.assembleAddMatrix( this->_mat ); //mat += M[GradU2]

    const RealType delta = 0.1; //sqrt(0.1);
    //typedef qc::ATAnisotropy2Norm<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta , 1., 1.);
    typedef qc::ATAnisotropy1Norm<ConfiguratorType> AnisoType;
    AnisoType anisotropy( delta );
    //typedef qc::ATAnisotropyOktaeder<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta );
    //AnisotropyShapeWriter<RealType, AnisoType> anisotropyShapeWriter(anisotropy);
    //anisotropyShapeWriter.writeWulffShape( "wulffshape.pgm" );
    //anisotropyShapeWriter.writeFrankDiagram( "frankdiagram.pgm" );

    VariationOfAnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> F(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);
    SecondVariationOfAnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> DF(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);

    aol::NewtonIteration<ConfiguratorType> newton_solver(this->_grid, F, DF, 2, InnerIterations, OuterIterations);
    newton_solver.apply(V, this->_rhs);
    V = this->_rhs;

    //AnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> E(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);

    //GradientDescent<ConfiguratorType> gradientDescent_solver(this->_grid, E, F, InnerIterations, OuterIterations);
    //gradientDescent_solver.apply(V, this->_rhs);
    //V = this->_rhs;
  }
  void updateVwithAnisotropyInUAndV(const aol::Vector<RealType> &U,
                             aol::Vector<RealType> &V,
                             const int InnerIterations = 0,
                             const int OuterIterations = 0){
    // use Newton Iteration to solve F(v)=0

    cerr << "Solving for v" << endl;
    const RealType delta = 0.1; //sqrt(0.1);
    //typedef qc::ATAnisotropy2Norm<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta , 1., 1.);
    typedef qc::ATAnisotropy1Norm<ConfiguratorType> AnisoType;
    AnisoType anisotropy( delta );
    //typedef qc::ATAnisotropyOktaeder<ConfiguratorType> AnisoType;
    //AnisoType anisotropy( delta );
    //AnisotropyShapeWriter<RealType, AnisoType> anisotropyShapeWriter(anisotropy);
    //anisotropyShapeWriter.writeWulffShape( "wulffshape.pgm" );
    //anisotropyShapeWriter.writeFrankDiagram( "frankdiagram.pgm" );

    this->_mat.setZero(); //mat = (alpha/(4*epsilon))*M
    this->_mass.assembleAddMatrix( this->_mat );
    this->_mat.scale( this->_alpha / ( 4. * this->_epsilon) );

    qc::ATMassAnisoU2Op<ConfiguratorType, AnisoType> atMassAnisoU2Op( this->_grid, U, anisotropy );
    atMassAnisoU2Op.assembleAddMatrix( this->_mat ); //mat += M[GradU2]

    VariationOfAnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> F(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);
    SecondVariationOfAnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> DF(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);

    aol::NewtonIteration<ConfiguratorType> newton_solver(this->_grid, F, DF, 2, InnerIterations, OuterIterations);
    newton_solver.apply(V, this->_rhs);
    V = this->_rhs;

    //AnisoATEnergyWithRegardToV<ConfiguratorType, AnisoType> E(this->_grid, anisotropy, this->_mat, this->_alpha, this->_epsilon);

    //GradientDescent<ConfiguratorType> gradientDescent_solver(this->_grid, E, F, InnerIterations, OuterIterations);
    //gradientDescent_solver.apply(V, this->_rhs);
    //V = this->_rhs;
  }
/*
  void updateVwithDoubleWell(const aol::Vector<RealType> &U,
                             aol::Vector<RealType> &V,
                             const int InnerIterations = 0,
                             const int OuterIterations = 0){
    // use Newton Iteration to solve F(v)=0

    cerr << "Solving for v" << endl;
    qc::ATMassGradU2Op<ConfiguratorType> atMassGradU2( this->_grid, U );

    this->_mat.setZero();

    _stiff.assembleAddMatrix( this->_mat ); //mat = alpha*epsilon*L
    this->_mat.scale( this->_alpha * this->_epsilon );

    atMassGradU2.assembleAddMatrix( this->_mat ); //mat += M[GradU2]

    this->_mat.scale( (-2. * this->_epsilon) / this->_alpha); // mat += -(alpha/(2*epsilon))*M
    this->_mass.assembleAddMatrix( this->_mat );
    this->_mat.scale( this->_alpha / ( -2. * this->_epsilon) );
    VariationOfDoubleWellATEnergy<ConfiguratorType> F(this->_grid, this->_mat, this->_alpha, this->_epsilon);
    SecondVariationOfDoubleWellATEnergy<ConfiguratorType> DF(this->_grid, this->_mat, this->_alpha, this->_epsilon);
    NewtonIteration<ConfiguratorType> newton_solver(this->_grid, F, DF, 2, InnerIterations, OuterIterations);
    newton_solver.apply(V, this->_rhs);
    V = this->_rhs;
  }
*/
  void segmentWithAnisotropy(){

    ArrayType u( this->_grid );
    ArrayType v( this->_grid );
    ArrayType tmp( this->_grid );

    //char fn[1024];

    v.setAll( 0. );

    for ( int outer=0; outer < 5; outer++ ) {
      cerr << ">>>>>>>> epsilon = " << this->_epsilon << endl;

      for (int iter = 0; iter < 10; iter++)
      {
// First, given v solve for u:

        //updateU( u, v, iter, outer);
        updateUwithAnisotropy( u, v, iter, outer);

// write output
/*
        u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
        sprintf( fn, "results/u_%02d_%03d.pgm", outer, iter );
        u.save( fn );
*/
        u.saveSlices( "results/u_old%2d.pgm", qc::QC_X, 5, NULL, aol::CLIP_THEN_SCALE, 0, 1);

// now solve for v:

        updateVwithAnisotropyInUAndV(u, v, iter, outer);
        //updateV(u, v, iter, outer);
        //updateVwithDoubleWell(u, v, iter, outer);

// write output...
/*
    v.setOverflowHandling( aol::CLIP_THEN_SCALE, 0., 1. );
    sprintf( fn, "results/v_%02d_%03d.pgm", outer, iter );
    v.save( fn );
*/
        v.saveSlices( "results/u_old%2d.pgm", qc::QC_X, 5, NULL, aol::CLIP_THEN_SCALE, 0, 1);
      }
      this->_epsilon /= 2.;
    }
  }

  // basic segmentation without anisotropi
  void BasicSegment( const int MaxIterations = 1 ) {
    ArrayType u( this->_grid );
    ArrayType v( this->_grid );

    segment ( u, v, MaxIterations );

    qc::writeImage<RealType> ( this->_grid, u, "results/u", true );
    qc::writeImage<RealType> ( this->_grid, v, "results/v", true );
  }

  void updateUwithKernel( aol::Vector<RealType> &U,
                          const aol::Vector<RealType> &V,
                          const int InnerIterations = 0,
                          const int /*OuterIterations*/ = 0 ){

    qc::Convolution<qc::QC_2D> conv( this->_grid.getWidth(), this->_grid.getHeight() );

    ArrayType kernel( this->_grid );
    ArrayType kernelMinus( this->_grid );
    typename ConfiguratorType::VecType v (10., 10.);
    qc::generateMotionBlurKernel<RealType>( v, kernel );
    kernelMinus.pointMirrorAtOrigin( kernel );

    DerivativeOfL2ConvolutionFidelityTermWithRegardToU<ConfiguratorType> DE( this->_grid, kernel, kernelMinus, conv );

    aol::LinCombOp<aol::Vector<RealType> > op;

    cerr << "Solving for u at iteration " << InnerIterations << endl;
    cerr << "Assembling matrices ";
    // \beta \int (u \phi) +  \int (v^2 \nabla u \cdot \nabla \phi) = \beta \int (u_0 \phi)
    qc::ATStiffV2Op<ConfiguratorType> atStiffV2( this->_grid, V );

    this->_mat.setZero();

    op.appendReference( DE, this->_beta );
    this->_mass.assembleAddMatrix( this->_mat ); // op = beta * Mconv

    atStiffV2.assembleAddMatrix( this->_mat ); // op += L[v^2]
    op.appendReference( this->_mat );
    cerr << "done\n";

    // rhs could be cached to get a minor speed increase
    ArrayType tmp( this->_grid );
    conv.convolve( this->_u0, kernelMinus, tmp );
    this->_mass.apply( tmp, this->_rhs ); //rhs = beta * M (u_0 conv h^-)
    this->_rhs *= this->_beta;

    aol::CGInverse<aol::Vector<RealType> > solver( op, 1e-10, 1000 );
    solver.setStopping( aol::STOPPING_ABSOLUTE );
    solver.apply( this->_rhs, U );
  }

};

template <typename ConfiguratorType>
class ATSegmentationMultilevelDescent {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  const RealType _alpha, _beta, _presmoothFactor, _presmooth;
  const int _max_depth;
  const qc::GridDefinition &_grid;
  const qc::GridDefinition *_cur_grid;
  int _cur_level;
  qc::MultilevelArray<RealType> _reference, _org_reference, _reference0, _v;
  vector<const qc::GridDefinition*> _grids;
  RealType _epsilon;

public:
  ATSegmentationMultilevelDescent( aol::ParameterParser &Parser, qc::GridDefinition &Grid )
    :_alpha(Parser.getDouble( "alpha" )),
     _beta(Parser.getDouble( "beta" )), _presmoothFactor(Parser.getDouble("presmoothFactor")),
     _presmooth(Parser.getDouble( "presmooth" )),
     _max_depth( Parser.getInt( "highestLevel" ) ),
     _grid( Grid ),
     _cur_grid( NULL ),
     _cur_level( _max_depth ),
     _reference( _grid ),
     _org_reference( _grid ),
     _reference0( _grid ),
     _v( _grid )
  {
    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
      _grids.push_back( new qc::GridDefinition( level, ConfiguratorType::Dim ) );
    }

    _cur_grid = _grids[ _cur_level ];


    char fn_reference[1024];
    Parser.getString( "reference", fn_reference );
    /************************************ BEGIN: loading image *************************************************/
    //cerr << "reading image from: " << fn_reference << endl;

    ArrayType ref( _reference.current( ), aol::FLAT_COPY );

    ref.load( fn_reference );
    ref /= ref.getMaxValue();
    _org_reference.current( ) = _reference.current( );
/*
    // presmooth the images
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( *_cur_grid );
    linSmooth.setSigma( _presmooth * _cur_grid->H() );
    linSmooth.apply( ref, ref );
*/
    _reference.levRestrict( 0, _cur_level );
    _org_reference.levRestrict( 0, _cur_level );
 }

  virtual ~ATSegmentationMultilevelDescent( ) {
    for ( int level = 0; level <= _grid.getGridDepth(); level++ ) {
     delete _grids[ level ];
   }
  }

  RealType H() {
    return _cur_grid->H();
  }

  void setLevel( int Level ) {
    _cur_level = Level;
    _cur_grid = _grids[ _cur_level ];
    _reference.setCurLevel( Level );
    _org_reference.setCurLevel( Level );
    _reference0.setCurLevel( Level );
    _v.setCurLevel( Level );
  }

  int getLevel( ) const {
    return _cur_level;
  }

  void prolongate( ) {
    if ( _cur_level < _grid.getGridDepth() ) {
      _v.levProlongate( );
      _reference.levProlongate( );
      setLevel( _cur_level + 1 );
    }
  }
  virtual void descentOnCurrentGrid( RealType epsilon );

};

template <typename ConfiguratorType>
void ATSegmentationMultilevelDescent<ConfiguratorType>::descentOnCurrentGrid( RealType Epsilon ) {

  _epsilon = Epsilon;
/*
  qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> mat( *_cur_grid );

  aol::SSORPreconditioner<aol::Vector<RealType>, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> > ssor( mat );
  aol::PCGInverse<aol::Vector<RealType> > solver( mat, ssor, 1e-10, 1000 );
  solver.setQuietMode( false );

  qc::ScalarArray<RealType, qc::QC_2D> rhs( *_cur_grid );
  qc::ScalarArray<RealType, qc::QC_2D> tmp( *_cur_grid );
*/
  qc::LinearSmoothOp<RealType> linSmooth0;
  linSmooth0.setCurrentGrid( *_cur_grid );
  linSmooth0.setSigma( _presmoothFactor *(*_cur_grid).H() );
  linSmooth0.apply( _org_reference[ _cur_level ], _reference0[ _cur_level ] );

  qc::writeImage<RealType> ( *_cur_grid, _reference0[ _cur_level ], "results/smoothed_reference0" );

  ExtendedATSegmentation<ConfiguratorType> atSegmentation( Epsilon, _alpha, _beta, _reference0[ _cur_level ]);
  ArrayType differenceU( *_cur_grid );
  ArrayType differenceV( *_cur_grid );

  char fn[1024];
  ///////////////////////////////////// BEGIN ITERATION ////////////////////////////////////

  cerr << "\n\n----------------------------------------------------\n";
  cerr << "Segmentation on level " << _cur_level << " started\n";
  cerr << "----------------------------------------------------\n\n";
  int iter = 0;
  do{
    differenceU = _reference[ _cur_level ];
    differenceV = _v[ _cur_level ];
    atSegmentation.updateU( _reference[ _cur_level ], _v[ _cur_level ], iter);
    //atSegmentation.updateUwithKernel( _reference[ _cur_level ], _v[ _cur_level ], iter);
    //atSegmentation.updateUwithAnisotropy( _reference[ _cur_level ], _v[ _cur_level ], iter);

    //sprintf( fn, "results/u_%f_%03d", Epsilon, iter );
    sprintf( fn, "results/u_%02d_%03d", _cur_level, iter );
    qc::writeImage<RealType> ( *_cur_grid, _reference[ _cur_level ], fn );

    atSegmentation.updateV( _reference[ _cur_level ], _v[ _cur_level ], iter);
    //atSegmentation.updateVwithAnisotropyInUAndV( _reference[ _cur_level ], _v[ _cur_level ], iter);

    //sprintf( fn, "results/v_%f_%03d", Epsilon, iter );
    sprintf( fn, "results/v_%02d_%03d", _cur_level, iter );
    qc::writeImage<RealType> ( *_cur_grid, _v[ _cur_level ], fn );

/*
    qc::ATEnergyOp<ConfiguratorType> atEnergyOp( *_cur_grid,
                                             _reference0[_cur_level ],
                                             _beta,
                                             1.,
                                             _alpha,
                                             _epsilon );
    aol::MultiVector<RealType> uvReference(0,0);
    uvReference.appendReference( _reference[_cur_level] );
    uvReference.appendReference( _v[_cur_level] );

    aol::Scalar<RealType> scalarTemp;
    atEnergyOp.apply( uvReference, scalarTemp );
    cerr << "ATCost = " << scalarTemp[0] << endl;
*/

    differenceU -= _reference[ _cur_level ];
    differenceV -= _v[ _cur_level ];
    iter++;
  } while ( ((differenceU.norm() > 0.1) || (differenceV.norm() > 0.1)) && ( iter < 500));
}

int main( int argc, char **argv ) {
  try {
    char parameterfilename[1024];

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf( parameterfilename, "%s",  argv[1] );
    }
    if ( argc == 1 ) {
      if( ConfType::Dim == 3)
        sprintf( parameterfilename, "ParameterFiles/Seg_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/Seg_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    aol::StopWatch watch;
    watch.start();

/*
    ExtendedATSegmentation<ConfType> atSegmentation( parser );
    atSegmentation.BasicSegment();
*/
    const RType epsilon = parser.getDouble( "epsilon" );
    //const RType epsilonFactor = parser.getDouble( "epsilonFactor" );
    const int highestLevel = parser.getInt( "highestLevel" );
    const int startLevel = parser.getInt( "startLevel" );
    const int stopLevel = parser.getInt( "stopLevel" );
    qc::GridDefinition Grid( highestLevel, ConfType::Dim );
    ATSegmentationMultilevelDescent<ConfType> mld( parser, Grid );
    mld.setLevel( startLevel );
    for ( int level = startLevel + 1; level <= stopLevel; level++ ) {
      mld.prolongate( );
      //RType epsilon = mld.H() * epsilonFactor;
      cerr << endl << "epsilon = " << epsilon << endl << endl;
      mld.descentOnCurrentGrid( epsilon );
    }

    watch.stop();
    watch.printReport( cerr );
  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}

