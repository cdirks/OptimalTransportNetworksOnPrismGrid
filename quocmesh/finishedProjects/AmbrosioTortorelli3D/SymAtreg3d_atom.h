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

#ifndef __SYMATREG3D_ATOM_H
#define __SYMATREG3D_ATOM_H

/**
 * \file
 * \brief  The classes of symmetric registration on single level
 * \author Berkels, Han
 */

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <gradientDescent.h>
#include <AmbrosioTortorelli.h>
#include <registration.h>
#include "utilities.h"

// E
template <typename ConfiguratorType>
void APlus1MinusB( const qc::GridDefinition &Grid,
                   const aol::Vector<typename ConfiguratorType::RealType> &A,
                   const aol::Vector<typename ConfiguratorType::RealType> &B,
                   aol::Vector<typename ConfiguratorType::RealType> &Dest){
  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper( Grid );
  for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {
    const int index = mapper.getGlobalIndex(*fnit);
    Dest[ index ] = A[ index ] + 1. - B[ index ];
  }
}

template <typename ConfiguratorType>
class SymATRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,
                                                          aol::Scalar<typename ConfiguratorType::RealType> >
{
public:
  // type define
  typedef typename ConfiguratorType::RealType RealType;

  // internal grid
  const typename ConfiguratorType::InitType &_grid;

  //member variables
#ifdef INCLUDE_WEIGHTS
  const aol::Vector<RealType>               &_weight;
#endif
#ifdef INCLUDE_WEIGHTS2
  const aol::Vector<RealType>               &_weight_other;
#endif
  const aol::Vector<RealType>               &_U;
  const aol::Vector<RealType>               &_V;
  const aol::MultiVector<RealType>          &_invT;
  const RealType                            _mu;
  const RealType                            _lambda;
  const RealType                            _kappa;

  const aol::StiffOp<ConfiguratorType>      &_stiff;
  // store the energy component
  mutable RealType                          _lastCoupledEnergy;
  mutable RealType                          _lastInternalEnergy;
  mutable RealType                          _lastConsistencyEnergy;
  mutable RealType                          _lastGlobalEnergy;
public:

//  // initialization constructor
  SymATRegistrationEnergyWithRegardToPhi(const qc::GridDefinition &Grid,
                                         const aol::StiffOp<ConfiguratorType> &Stiff,
#ifdef INCLUDE_WEIGHTS
                                         const aol::Vector<RealType> &Weight,
#endif
#ifdef INCLUDE_WEIGHTS2
                                         const aol::Vector<RealType> &Weight_other,
#endif
                                         const aol::Vector<RealType> &u,
                                         const aol::Vector<RealType> &v,
                                         const aol::MultiVector<RealType> &invT,
                                         const RealType mu,
                                         const RealType lambda,
                                         const RealType kappa):
    _grid(Grid),
#ifdef INCLUDE_WEIGHTS
    _weight(Weight),
#endif
#ifdef INCLUDE_WEIGHTS2
    _weight_other(Weight_other),
#endif
    _U(u),
    _V(v),
    _invT(invT),
    _mu(mu),
    _lambda(lambda),
    _kappa(kappa),
    _stiff(Stiff),
    _lastCoupledEnergy(0.),
    _lastInternalEnergy(0.),
    _lastConsistencyEnergy(0.),
    _lastGlobalEnergy(0.)
  {
  }

  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest) const{
    aol::Scalar<RealType> temp;

    RealType E_ext = CoupledEnergy(_grid,_U,_V,Arg);

    qc::DisplacementLengthEnergy<ConfiguratorType> displacementLenghtEnergy(_stiff);
    displacementLenghtEnergy(Arg, temp);
    RealType E_int = temp[0];

    ConsistencyEnergy<ConfiguratorType> consisEnergy(_grid, _invT);
    consisEnergy.apply(Arg,temp);
    RealType E_con = temp[0];

    _lastCoupledEnergy = E_ext;
    _lastInternalEnergy = E_int;
    _lastConsistencyEnergy = E_con;
    Dest[0] = _lastGlobalEnergy = _mu*_lastCoupledEnergy + _lambda*_lastInternalEnergy + _kappa*_lastConsistencyEnergy;
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
  void getLastEnergy( aol::Vector<RealType> &Energy ) const{
    Energy[0] = _lastGlobalEnergy;
    Energy[1] = _lastCoupledEnergy;
    Energy[2] = _lastInternalEnergy;
    Energy[3] = _lastConsistencyEnergy;
  }
private:
  RealType CoupledEnergy( const qc::GridDefinition &Grid,
                          const aol::Vector<RealType>  &U,
                          const aol::Vector<RealType>  &V,
                          const aol::MultiVector<RealType> &phi) const
  {
    // predefine the variables
    RealType E = 0.;
    aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );

    // E *= L[v^2\circ\phi]R*R
#ifdef INCLUDE_WEIGHTS
#ifdef INCLUDE_WEIGHTS2
    aol::Vector<RealType> tmp2( Grid.getNumberOfNodes() );
    APlus1MinusB<ConfiguratorType>(Grid, V, _weight_other, tmp2);
    ATWeightedStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( Grid, _weight, tmp2, phi);
#else
    ATWeightedStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( Grid, _weight, V, phi);
#endif
#else
    qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( Grid, V, phi);
#endif
    atStiffV2Deform.apply(U, tmp);
    E += U * tmp;

    E /= 2.;

    return E;
  }
};

template <typename ConfiguratorType>
class VariationOfSymATRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> >
{
  // type define
  typedef typename ConfiguratorType::RealType RealType;

  // member variables
  const typename ConfiguratorType::InitType &_grid;
#ifdef INCLUDE_WEIGHTS
  const aol::Vector<RealType>               &_weight;
#endif
#ifdef INCLUDE_WEIGHTS2
  const aol::Vector<RealType>               &_weight_other;
#endif
  const aol::Vector<RealType>               &_U;
  const aol::Vector<RealType>               &_V;
  const aol::MultiVector<RealType>          &_extT;
  const RealType _mu, _lambda, _kappa;

  const aol::DiagonalBlockOp<RealType> _blockStiff;

public:
  // initialization constructor
  VariationOfSymATRegistrationEnergyWithRegardToPhi(const typename ConfiguratorType::InitType &grid,
                                                    const aol::StiffOp<ConfiguratorType> &Stiff,
#ifdef INCLUDE_WEIGHTS
                                                    const aol::Vector<RealType> &Weight,
#endif
#ifdef INCLUDE_WEIGHTS2
                                                    const aol::Vector<RealType> &Weight_other,
#endif
                                                    const aol::Vector<RealType> &u,
                                                    const aol::Vector<RealType> &v,
                                                    const aol::MultiVector<RealType>  &extT,
                                                    const RealType mu,
                                                    const RealType lambda,
                                                    const RealType kappa):
    _grid(grid),
#ifdef INCLUDE_WEIGHTS
    _weight(Weight),
#endif
#ifdef INCLUDE_WEIGHTS2
    _weight_other(Weight_other),
#endif
    _U(u),
    _V(v),
    _extT(extT),
    _mu(mu),
    _lambda(lambda),
    _kappa(kappa),
    _blockStiff(Stiff)
  {
  }

  // compute the energy
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest) const
  {
    // predefine the operators
#ifdef INCLUDE_WEIGHTS
#ifdef INCLUDE_WEIGHTS2
    aol::Vector<RealType> tmp2( _grid.getNumberOfNodes() );
    APlus1MinusB<ConfiguratorType>(_grid, _V, _weight_other, tmp2);
    ATWeightedDeformationGradient<ConfiguratorType,ConfiguratorType::Dim> defGradOp( _grid, _weight, tmp2, _U);
#else
    ATWeightedDeformationGradient<ConfiguratorType,ConfiguratorType::Dim> defGradOp( _grid, _weight, _V, _U);
#endif
#else
    qc::ATDeformationGradient<ConfiguratorType,ConfiguratorType::Dim> defGradOp( _grid, _V, _U);
#endif
    SymTransSubDetOp<ConfiguratorType> TransSubDetOp(_grid,_extT);
    qc::VariationOfConsistencyEnergyOp<ConfiguratorType> TransSubJobJobOp(_grid,_extT);

    //link the operators
    aol::LinCombOp<aol::MultiVector<RealType> > LinkOp;
    LinkOp.appendReference(defGradOp,_mu);
    LinkOp.appendReference(_blockStiff,_lambda);
    LinkOp.appendReference(TransSubJobJobOp,_kappa);
    LinkOp.appendReference(TransSubDetOp,_kappa);

    //solve
    LinkOp.apply(Arg,Dest);
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/) const
  {
      throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
private:

};

template <typename ConfiguratorType>
class SymATRegistrationEnergy: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,
                                              aol::Scalar<typename ConfiguratorType::RealType> >
{
public:
  // type define
  typedef typename ConfiguratorType::RealType RealType;

  // internal grid
  const typename ConfiguratorType::InitType &_grid;

  //member variables
  const aol::Vector<RealType>               &_uR;
  const aol::Vector<RealType>               &_uT;
  const aol::Vector<RealType>               &_vR;
  const aol::Vector<RealType>               &_vT;
  const RealType                            _mu;
  const RealType                            _lambda;
  const RealType                            _kappa;

  const aol::StiffOp<ConfiguratorType>      &_stiff;
  // store the energy component
  mutable RealType                          _lastCoupledEnergy_R2T;
  mutable RealType                          _lastInternalEnergy_R2T;
  mutable RealType                          _lastConsistencyEnergy_R2T;
  mutable RealType                          _lastCoupledEnergy_T2R;
  mutable RealType                          _lastInternalEnergy_T2R;
  mutable RealType                          _lastConsistencyEnergy_T2R;
  mutable RealType                          _lastGlobalEnergy;
public:

//  // initialization constructor
  SymATRegistrationEnergy( const qc::GridDefinition &Grid,
                           const aol::StiffOp<ConfiguratorType> &Stiff,
                           const aol::Vector<RealType> &UR,
                           const aol::Vector<RealType> &UT,
                           const aol::Vector<RealType> &VR,
                           const aol::Vector<RealType> &VT,
                           const RealType Mu,
                           const RealType Lambda,
                           const RealType Kappa):
    _grid(Grid),
    _uR(UR),
    _uT(UT),
    _vR(VR),
    _vT(VT),
    _mu(Mu),
    _lambda(Lambda),
    _kappa(Kappa),
    _stiff(Stiff),
    _lastCoupledEnergy_R2T(0.),
    _lastInternalEnergy_R2T(0.),
    _lastConsistencyEnergy_R2T(0.),
    _lastCoupledEnergy_T2R(0.),
    _lastInternalEnergy_T2R(0.),
    _lastConsistencyEnergy_T2R(0.),
    _lastGlobalEnergy(0.)
  {
  }

  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest) const{
    aol::MultiVector<RealType> phi_T2R(0,0);
    aol::MultiVector<RealType> phi_R2T(0,0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      phi_T2R.appendReference( Arg[i] );
      phi_R2T.appendReference( Arg[i+ConfiguratorType::Dim] );
    }

    aol::Scalar<RealType> temp;

    _lastCoupledEnergy_T2R = CoupledEnergy(_grid,_uR,_vT,phi_T2R);
    _lastCoupledEnergy_R2T = CoupledEnergy(_grid,_uT,_vR,phi_R2T);

    qc::DisplacementLengthEnergy<ConfiguratorType> displacementLenghtEnergy(_stiff);
    displacementLenghtEnergy.apply(phi_T2R, temp);
    _lastInternalEnergy_T2R = temp[0];
    displacementLenghtEnergy.apply(phi_R2T, temp);
    _lastInternalEnergy_R2T = temp[0];

    //ConsistencyEnergy<ConfiguratorType> consisEnergy(_grid, phi_R2T);
    //consisEnergy.apply(phi_T2R,temp);
    //RealType temp0 = temp[0];
    qc::ConsistencyEnergyOp<ConfiguratorType> consisEnergyPhi_T2R(_grid, phi_R2T);
    consisEnergyPhi_T2R.apply(phi_T2R,temp);
    _lastConsistencyEnergy_T2R = temp[0];
    qc::ConsistencyEnergyOp<ConfiguratorType> consisEnergyPhi_R2T(_grid, phi_T2R);
    consisEnergyPhi_R2T.apply(phi_R2T,temp);
    _lastConsistencyEnergy_R2T = temp[0];

    //cerr << temp0 << endl;
    //cerr << _lastConsistencyEnergy_R2T + _lastConsistencyEnergy_T2R << endl;

    Dest[0] = _lastGlobalEnergy = _mu*(_lastCoupledEnergy_R2T+_lastCoupledEnergy_T2R)
                               + _lambda*(_lastInternalEnergy_R2T + _lastInternalEnergy_T2R)
                               + _kappa*(_lastConsistencyEnergy_R2T + _lastConsistencyEnergy_T2R);
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
  void getLastEnergy( aol::Vector<RealType> &Energy ) const{
    Energy.resize(7);
    Energy[0] = _lastGlobalEnergy;
    Energy[1] = _lastCoupledEnergy_T2R;
    Energy[2] = _lastCoupledEnergy_R2T;
    Energy[3] = _lastInternalEnergy_T2R;
    Energy[4] = _lastInternalEnergy_R2T;
    Energy[5] = _lastConsistencyEnergy_T2R;
    Energy[6] = _lastConsistencyEnergy_R2T;


  }
private:
  RealType CoupledEnergy( const qc::GridDefinition &Grid,
                          const aol::Vector<RealType>  &U,
                          const aol::Vector<RealType>  &V,
                          const aol::MultiVector<RealType> &phi) const
  {
    // predefine the variables
    RealType E = 0.;
    aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );

    // E *= L[v^2\circ\phi]R*R
    qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( Grid, V, phi);
    atStiffV2Deform.apply(U, tmp);
    E += U * tmp;

    E /= 2.;

    return E;
  }
};

template <typename ConfiguratorType>
class VariationOfSymATRegistrationEnergy: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> >
{
  // type define
  typedef typename ConfiguratorType::RealType RealType;

  // member variables
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType>               &_uR;
  const aol::Vector<RealType>               &_uT;
  const aol::Vector<RealType>               &_vR;
  const aol::Vector<RealType>               &_vT;
  const RealType _mu, _lambda, _kappa;

  const aol::StiffOp<ConfiguratorType>      &_stiff;

public:
  // initialization constructor
  VariationOfSymATRegistrationEnergy( const typename ConfiguratorType::InitType &Grid,
                                      const aol::StiffOp<ConfiguratorType> &Stiff,
                                      const aol::Vector<RealType> &UR,
                                      const aol::Vector<RealType> &UT,
                                      const aol::Vector<RealType> &VR,
                                      const aol::Vector<RealType> &VT,
                                      const RealType Mu,
                                      const RealType Lambda,
                                      const RealType Kappa):
    _grid(Grid),
    _uR(UR),
    _uT(UT),
    _vR(VR),
    _vT(VT),
    _mu(Mu),
    _lambda(Lambda),
    _kappa(Kappa),
    _stiff(Stiff)
  {
  }

  // compute the energy
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest) const
  {
    aol::MultiVector<RealType> phi_T2R(0,0);
    aol::MultiVector<RealType> phi_R2T(0,0);
    aol::MultiVector<RealType> phi_T2R_Dest(0,0);
    aol::MultiVector<RealType> phi_R2T_Dest(0,0);
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      phi_T2R.appendReference( Arg[i] );
      phi_R2T.appendReference( Arg[i+ConfiguratorType::Dim] );
      phi_T2R_Dest.appendReference( Dest[i] );
      phi_R2T_Dest.appendReference( Dest[i+ConfiguratorType::Dim] );
    }

    VariationOfSymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> DEPhi_T2R(_grid, _stiff, _uR, _vT, phi_R2T, _mu, _lambda, _kappa);
    DEPhi_T2R.apply( phi_T2R, phi_T2R_Dest );
    VariationOfSymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> DEPhi_R2T(_grid, _stiff, _uT, _vR, phi_T2R, _mu, _lambda, _kappa);
    DEPhi_R2T.apply( phi_R2T, phi_R2T_Dest );


  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/) const
  {
      throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
private:

};


template <typename ConfiguratorType>
class SymATRegistration : public qc::RegistrationInterface<ConfiguratorType> {

  /**Type define*/
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  typedef qc::RegistrationInterface<ConfiguratorType>  SuperClass;


  const RealType _alpha;
  const RealType _beta;
  const RealType _nu;
  const RealType _mu;
  const RealType _kappa;
  const RealType _epsilon;
  const RealType _sigma;
  const int _numGradientSteps;
  RealType _filterWidth;

    /** Predefined mass operator*/
  const aol::MassOp<ConfiguratorType> _mass;

  /** Predefined stiff operator*/
  const aol::StiffOp<ConfiguratorType> _stiff;

  /** Predefined Lumped mass inversion operator*/
  const aol::LumpedMassOp<ConfiguratorType> _lumpedMassInv;

  /** The system matrix*/
  qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> _mat;

  /** The sor precondition operator*/
  aol::SSORPreconditioner<aol::Vector<RealType>, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> > _ssor;

  /** The CG solver*/
  aol::PCGInverse<aol::Vector<RealType> > _solver;

public:
  SymATRegistration(RealType alpha,
                    RealType beta,
                    RealType nu,
                    RealType mu,
                    RealType kappa,
                    RealType lambda,
                    RealType epsilon,
                    RealType sigma,
                    const int NumGradientSteps,
                    const ArrayType &Reference0,
                    const ArrayType &Template0)
  :SuperClass(Reference0,Template0,lambda),
  _alpha(alpha),
  _beta(beta),
  _nu(nu),
  _mu(mu),
  _kappa(kappa),
  _epsilon(epsilon),
  _sigma(sigma),
  _numGradientSteps(NumGradientSteps),
  _filterWidth( sqrt(10.*this->_grid.H()) ),
  _mass (this->_grid ),
  _stiff (this->_grid, aol::ASSEMBLED),
  _lumpedMassInv(this->_grid, aol::INVERT ),
  _mat( this->_grid ),
  _ssor( _mat ),
  _solver( _mat, _ssor, 1e-10, 1000 )
  {
    _solver.setQuietMode( false );
    _solver.setStopping( aol::STOPPING_ABSOLUTE );
  }


  void update_U(const aol::Vector<RealType> &Image,
#ifdef INCLUDE_WEIGHTS
                const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS2
                const aol::Vector<RealType> &Weight_other,
#endif
                const aol::Vector<RealType> &V_own,
                const aol::Vector<RealType> &V_other,
                const aol::MultiVector<RealType> &Tran,
                aol::Vector<RealType> &U);

  void update_V(const aol::Vector<RealType> &U_own,
#ifdef INCLUDE_WEIGHTS2
                const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS
                const aol::Vector<RealType> &Weight_other,
#endif
                const aol::Vector<RealType> &U_other,
                //const aol::Vector<RealType> &V_other,
                const aol::MultiVector<RealType> &Tran,
                aol::Vector<RealType> &V_own);

  void update_T(const aol::Vector<RealType> &U,
                const aol::Vector<RealType> &V,
                const aol::MultiVector<RealType> &Tran_in,
                int ID,
                std::ofstream& os,
                aol::MultiVector<RealType> &Tran_out);

  void update_T_AM(const aol::Vector<RealType>&U,
                   const aol::Vector<RealType>&V,
                   const aol::MultiVector<RealType> &Tran_in,
#ifdef INCLUDE_WEIGHTS
                   const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS2
                   const aol::Vector<RealType> &Weight_other,
#endif
                   int ID,
                   std::ofstream& osEnergy,
                   std::ofstream& osTau,
                   aol::MultiVector<RealType> &Tran_out);


  void updateBothTransformations( const aol::Vector<RealType> &UR,
                                  const aol::Vector<RealType> &UT,
                                  const aol::Vector<RealType> &VR,
                                  const aol::Vector<RealType> &VT,
                                  aol::MultiVector<RealType> &Phi_T2R,
                                  aol::MultiVector<RealType> &Phi_R2T,
                                  int ID,
                                  std::ofstream &osEnergy );

};

template <typename ConfiguratorType>
void SymATRegistration<ConfiguratorType>
::update_U( const aol::Vector<RealType>&Image,
#ifdef INCLUDE_WEIGHTS
            const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS2
            const aol::Vector<RealType> &Weight_other,
#endif
            const aol::Vector<RealType>&V_own,
            const aol::Vector<RealType>&V_other,
            const aol::MultiVector<RealType> &Tran,
            aol::Vector<RealType> &U)
{
/*
  \alpha\int_\Omega u_R\ \vartheta\;dx +\beta\int_\Omega v_R^2\nabla
  u_R\cdot\nabla\vartheta\;dx +\mu\int_{\Omega}(v_T\circ h)^2\nabla
  u_R \cdot \nabla\vartheta\;dx = \alpha\int_\Omega u_R^0 \
  \vartheta\;dx.
*/

  // user-defined arraies
  ArrayType  rhs( this->_grid );

  // user defined operator
  qc::ATStiffV2Op<ConfiguratorType> atStiffV2( this->_grid, V_own );
#ifdef INCLUDE_WEIGHTS
#ifdef INCLUDE_WEIGHTS2
  aol::Vector<RealType> tmp2( this->_grid.getNumberOfNodes() );
  APlus1MinusB<ConfiguratorType>(this->_grid, V_other, Weight_other, tmp2);
  ATWeightedStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( this->_grid, Weight_own, tmp2, Tran);
#else
  ATWeightedStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( this->_grid, Weight_own, V_other, Tran);
#endif
#else
  qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( this->_grid, V_other, Tran);
#endif

  // clear the system matrix
  _mat.setZero();

  // mat =  M
  _mass.assembleAddMatrix( _mat);

  // mat = \frac{\alpha}{\beta}M
  _mat *= ( _alpha / _beta );

  // mat = \frac{\alpha}{\beta}M + S[v_R^2]
  atStiffV2.assembleAddMatrix(_mat);

  // mat = \frac{\alpha}{\mu}M + \frac{\beta}{\mu}S[v_R^2]
  _mat *= ( _beta / _mu );

  // mat = \frac{\alpha}{\mu}M +\frac{\beta}{\mu}S[v_R^2] + S[(v_T\circ h)^2]
  atStiffV2Deform.assembleAddMatrix(_mat);

  // right side vector
  // rhs = \frac{\alpha}{\mu}M^0
  _mass.apply( Image, rhs );
  // Intead of multiplying the matrix with _mu, multiply the rhs with 1./_mu
  rhs *= _alpha/_mu;

  // solve the linear equation and save the result in U
  _solver.apply( rhs, U );
}

template <typename ConfiguratorType>
void SymATRegistration<ConfiguratorType>
::update_V( const aol::Vector<RealType> &U_own,
#ifdef INCLUDE_WEIGHTS2
            const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS
            const aol::Vector<RealType> &Weight_other,
#endif
            const aol::Vector<RealType> &U_other,
            //const aol::Vector<RealType> &V_other,
            const aol::MultiVector<RealType> &Tran,
            aol::Vector<RealType> &V_own)
{
/*
  \beta\int_\Omega\Abs{\nabla u_R}^2v_R\vartheta\;dx
  +\frac{\nu}{4\epsilon}\int_\Omega v_R\vartheta \;dx
  +\nu\epsilon\int_\Omega\nabla v_R\cdot\nabla \vartheta\;dx \\+
  \mu\int_{\Omega}\|\nabla u_T\|^2\circ g^{-1}v_R\cdot\vartheta
  |det(D[g^{-1}](x))| \;dx =\frac{\nu}{4\epsilon}\int_\Omega
  \vartheta\;dx
*/

  ArrayType tmp( this->_grid );
  ArrayType rhs( this->_grid );
  //user defined operator
  qc::ATMassGradU2Op<ConfiguratorType> atMassGradU2(this->_grid, U_own );

  qc::ATGradSqrOp<ConfiguratorType> atGradSqrOp( this->_grid );
  atGradSqrOp.apply( U_other, rhs ); //rhs is just used as temp vector here to store gradRSqrDeform
#ifdef INCLUDE_WEIGHTS
  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper( this->_grid );
  for(fnit = this->_grid.begin(); fnit != this->_grid.end(); fnit++){
    const int index = mapper.getGlobalIndex(*fnit);
    rhs[ index ] = rhs[ index ] * Weight_other[ index ];
  }
#endif
  _lumpedMassInv.apply( rhs, tmp );

  qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( this->_grid );
  transform.setDeformation( Tran );
  transform.apply( tmp, rhs );

  qc::ATMassGradU2DeformOp<ConfiguratorType> atMassGradR2Deform( this->_grid, rhs, Tran );

  // clear system matrix
  _mat.setZero();

  // mat = M[\|\nabla u_R\|^2]
  atMassGradU2.assembleAddMatrix(_mat);

#ifdef INCLUDE_WEIGHTS2
  aol::Vector<RealType> tmp2( this->_grid.getNumberOfNodes() );
  aol::Vector<RealType> tmp3( this->_grid.getNumberOfNodes() );
  tmp2.setAll( -1. );
  tmp2 += Weight_own;
  _mat.apply( tmp2, tmp3 );
#endif

  // mat = \frac{4\beta\epsilon}{\nu}M[\|\nabla u_R\|^2]
  _mat *= ( 4.* _epsilon * _beta / _nu);

  // mat = \frac{4\beta\epsilon}{\nu}M[\|\nabla u_R\|^2]+ M
  _mass.assembleAddMatrix(_mat);

  // mat = \frac{\beta}{\nu\epsilon} M[\|\nabla u_R\|^2] + \frac{1}{4\epsilon^2} M
  _mat *= ( 0.25/(_epsilon*_epsilon));

  // mat = \frac{\beta}{\nu\epsilon} M[\|\nabla u_R\|^2] + \frac{1}{4\epsilon^2} M + S
  _stiff.assembleAddMatrix(_mat);

  // mat = \frac{\beta}{\mu}M[\|\nabla u_R\|^2] + \frac{\nu}{4\epsilon\mu}M + \frac{\nu\epsilon}{\mu}S
  _mat *= (_nu * _epsilon / _mu);

  // mat = \frac{\beta}{\mu}M[\|\nabla u_R\|^2] + \frac{\nu}{4\epsilon\mu}M + \frac{\nu\epsilon}{\mu}S+ M^{'}
  //simMassGradR2Deform.assembleAddMatrix(_mat);
  atMassGradR2Deform.assembleAddMatrix(_mat);

  // Intead of multiplying the matrix with _mu, multiply the rhs with 1./_mu
  tmp.setAll( 0.25*_nu/(_epsilon*_mu) );
  _mass.apply( tmp, rhs );
#ifdef INCLUDE_WEIGHTS2
  rhs += tmp3;
#endif
  // the content of rhs is not used anymore (atMassGradR2Deform is unusable from now on)

  _solver.apply( rhs, V_own);
}

template <typename ConfiguratorType>
void SymATRegistration<ConfiguratorType>
::update_T( const aol::Vector<RealType>&U,
            const aol::Vector<RealType>&V,
            const aol::MultiVector<RealType> &Tran_in,
            int ID,
            std::ofstream& os,
            aol::MultiVector<RealType> &Tran_out)
{
  cerr << "Deperacted. If you don't want to use time step with control, replace GradientDescent with SimpleGradientDescent in update_T_AM\n";
  /*
   \var{\partial_{h}E_{\text{G}},\psi} =
   \mu\int_\Omega\Abs{\nabla u_{R}}^2(v_{T}\circ h)\nabla(v_{T}\circ h )\cdot \psi \;dx
   + \lambda\int_\Omega D h:D\varphi \;dx \qquad \qquad\qquad\qquad+\kappa\int_{\Omega}
   (h(x)-g^{-1}(x))|\det(D[g^{-1}](x))|\cdot \psi(x) \;dx +\kappa \int_\Omega
   ([g\circ h](x) -x) D[g\circ h](x) \cdot \psi(x)\;dx
  */

  // predefine the operators
  qc::ATDeformationGradient<ConfiguratorType,ConfiguratorType::Dim> defGradOp( this->_grid, V, U);
  SymTransSubDetOp<ConfiguratorType> TransSubDetOp(this->_grid,Tran_in);
  qc::VariationOfConsistencyEnergyOp<ConfiguratorType> TransSubJobJobOp(this->_grid,Tran_in);

  /** Predefined Diagonal block stiff operator*/
  aol::DiagonalBlockOp<RealType> _blockStiff(_stiff);

  //link the operators
  aol::LinCombOp<aol::MultiVector<RealType> > LinkOp;
  LinkOp.appendReference(defGradOp,_mu);
  LinkOp.appendReference(_blockStiff,this->_lambda);
  LinkOp.appendReference(TransSubJobJobOp,_kappa);
  LinkOp.appendReference(TransSubDetOp,_kappa);

  // tmp
  aol::MultiVector<RealType> TempUnknown( ConfiguratorType::Dim,this->_grid.getNumberOfNodes() );
  aol::MultiVector<RealType> TempUnknown_2( ConfiguratorType::Dim,this->_grid.getNumberOfNodes() );
  //solve
  LinkOp.apply(Tran_out,TempUnknown);

  // inverse mass matrix
  for(int i= 0 ; i<ConfiguratorType::Dim;i++)
  {
    _lumpedMassInv.apply( TempUnknown[i], TempUnknown_2[i] );
  }
  TempUnknown_2*= -1.;

  // smoothing transform
  qc::LinearSmoothOp<RealType> linSmooth;
  linSmooth.setCurrentGrid(this->_grid);
  linSmooth.setTau( this->_tau * this->_grid.H() );
  for(int i= 0 ; i<ConfiguratorType::Dim;i++)
  {
    linSmooth.apply( TempUnknown_2[i], TempUnknown[i] );
  }


  // normalization
  RealType max = 0.0;
  for(int i= 0 ; i<ConfiguratorType::Dim;i++)
  {
     max = aol::Max(RealType(max),aol::Abs( TempUnknown[i].getMaxValue()));
     max = aol::Max(RealType(max),aol::Abs( TempUnknown[i].getMinValue()));
  }
  TempUnknown *= this->_grid.H() / max;


   Tran_out    += TempUnknown;

   ////////////////////////////////////////////
   // print out energy
   ////////////////////////////////////////////
   SymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> E(this->_grid,_stiff,U,V,Tran_in,_mu,this->_lambda,_kappa);
   aol::Vector<RealType> energy(4);
   aol::Scalar<RealType> temp;
   E.apply(Tran_out,temp);
   E.getLastEnergy(energy);

   char temString[1000];
   sprintf(temString,"%3d\t%e\t%e\t%e\t%e",ID,energy[0],energy[1],energy[2],energy[3]);
   os<<temString<<endl;
   cerr<<"T: "<<temString<<endl;

   return;


}


template <typename ConfiguratorType>
void SymATRegistration<ConfiguratorType>
::update_T_AM( const aol::Vector<RealType> &U,
               const aol::Vector<RealType> &V,
               const aol::MultiVector<RealType> &Tran_in,
#ifdef INCLUDE_WEIGHTS
               const aol::Vector<RealType> &Weight_own,
#endif
#ifdef INCLUDE_WEIGHTS2
               const aol::Vector<RealType> &Weight_other,
#endif
               int ID,
               std::ofstream &osEnergy,
               std::ofstream &osTau,
               aol::MultiVector<RealType> &Tran_out)
{
#ifdef INCLUDE_WEIGHTS
#ifdef INCLUDE_WEIGHTS2
  SymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> E(this->_grid,_stiff,Weight_own,Weight_other,U,V,Tran_in,_mu,this->_lambda,_kappa);
  VariationOfSymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> DE(this->_grid,_stiff,Weight_own,Weight_other,U,V, Tran_in,_mu,this->_lambda,_kappa);
#else
  SymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> E(this->_grid,_stiff,Weight_own,U,V,Tran_in,_mu,this->_lambda,_kappa);
  VariationOfSymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> DE(this->_grid,_stiff,Weight_own,U,V, Tran_in,_mu,this->_lambda,_kappa);
#endif
#else
  SymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> E(this->_grid,_stiff,U,V,Tran_in,_mu,this->_lambda,_kappa);
  VariationOfSymATRegistrationEnergyWithRegardToPhi<ConfiguratorType> DE(this->_grid,_stiff, U,V, Tran_in,_mu,this->_lambda,_kappa);
#endif
  aol::MultiVector<RealType> mtmp( ConfiguratorType::Dim, this->_grid.getNumberOfNodes() );

  aol::RecordedGradientDescent<ConfiguratorType> gradientDescent_solver( osTau,ID,this->_grid, E, DE, 1, this->_tau );
  //GradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( _grid, E, DE, 1 );
  gradientDescent_solver.apply(Tran_out, mtmp);
  this->_tau = gradientDescent_solver.getStartTau();
  Tran_out = mtmp;

  ////////////////////////////////////////////
  // print out energy
  ////////////////////////////////////////////

  aol::Vector<RealType> energy(4);
  E.getLastEnergy(energy);

  char temString[1000];
  sprintf(temString,"%3d\t%e\t%e\t%e\t%e",ID,energy[0],energy[1],energy[2],energy[3]);
  osEnergy<<temString<<endl;
  cerr<<"T: "<<temString<<endl;
}

template <typename ConfiguratorType>
void SymATRegistration<ConfiguratorType>
::updateBothTransformations( const aol::Vector<RealType> &UR,
                             const aol::Vector<RealType> &UT,
                             const aol::Vector<RealType> &VR,
                             const aol::Vector<RealType> &VT,
                             aol::MultiVector<RealType> &Phi_T2R,
                             aol::MultiVector<RealType> &Phi_R2T,
                             int ID,
                             std::ofstream &osEnergy )
{
  SymATRegistrationEnergy<ConfiguratorType> E( this->_grid, _stiff, UR, UT, VR, VT, _mu, this->_lambda, _kappa );
  VariationOfSymATRegistrationEnergy<ConfiguratorType> DE( this->_grid, _stiff, UR, UT, VR, VT, _mu, this->_lambda, _kappa );

  aol::MultiVector<RealType> mtmp( 2*ConfiguratorType::Dim, this->_grid.getNumberOfNodes() );

  //RecordedGradientDescent<ConfiguratorType> gradientDescent_solver( osTau,ID,this->_grid, E, DE, 1, this->_tau );
  aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, E, DE, _numGradientSteps, this->_tau, aol::ZOTrait<RealType>::zero);
  gradientDescent_solver.setFilterWidth( _filterWidth );

  aol::MultiVector<RealType> phiAndPsi( 0, 0 );
  for( int i = 0; i < ConfiguratorType::Dim; i++ )
    phiAndPsi.appendReference( Phi_T2R[i] );
  for( int i = 0; i < ConfiguratorType::Dim; i++ )
    phiAndPsi.appendReference( Phi_R2T[i] );

  gradientDescent_solver.apply(phiAndPsi, mtmp);
  this->_tau = gradientDescent_solver.getStartTau();
  _filterWidth = gradientDescent_solver.getFilterWidth();
  phiAndPsi = mtmp;

  ////////////////////////////////////////////
  // print out energy
  ////////////////////////////////////////////

  aol::Vector<RealType> energy(7);
  E.getLastEnergy(energy);

  char temString[1000];
  sprintf(temString,"%3d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",ID,energy[0],energy[1],energy[2],energy[3],
                                                                    energy[4],energy[5],energy[6],this->_tau);
  osEnergy<<temString<<endl;
  cerr<<"T: "<<temString<<endl;
}
#endif // __SYMATREG3D_ATOM_H
