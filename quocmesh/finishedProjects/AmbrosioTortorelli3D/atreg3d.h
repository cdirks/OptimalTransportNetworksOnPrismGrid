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

#ifndef __ATREG3D_H
#define __ATREG3D_H

#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <preconditioner.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <gradientDescent.h>
#include <AmbrosioTortorelli.h>
#include <registration.h>
#include <mutualInformation.h>
#include "ambrotelli.h"

// E
template <typename ConfiguratorType>
class ATRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> >
{
  /// type define
  typedef typename ConfiguratorType::RealType RealType;

  /// internal grid
  const typename ConfiguratorType::InitType &_grid;

  //aol::Vector<RealType> _phaseField;
  const aol::Vector<RealType> &_phaseField;
  const aol::Vector<RealType> &_reference;

  /// temp vector
  mutable aol::Vector<RealType> _tmp;

public:

  /// Initialization constructor
  ATRegistrationEnergyWithRegardToPhi(const typename ConfiguratorType::InitType &Initializer,
                                      const aol::Vector<RealType> &PhaseField,
                                      const aol::Vector<RealType> &Reference):
  _grid(Initializer),
  _phaseField(PhaseField),
  _reference(Reference),
  _tmp(_grid.getNumberOfNodes())
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest) const{
    RealType Energy = 0.;
    qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( _grid, _phaseField, Arg ); // E = 1/2*L[v^2\circ\phi]R*R
    atStiffV2Deform.apply(_reference, _tmp);
    Energy = (_reference * _tmp)/2.;
    Dest[0] = Energy;
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

template <typename ConfiguratorType>
class VariationOfATRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> >
{
  /// type define
  typedef typename ConfiguratorType::RealType RealType;

  /// internal grid
  const typename ConfiguratorType::InitType &_grid;

  const aol::Vector<RealType> &_phaseField;
  const aol::Vector<RealType> &_reference;

  /// The parameter lambda;
  const RealType _lambda;

  const aol::DiagonalBlockOp<RealType> _blockStiff;

public:

  /// Initialization constructor
  VariationOfATRegistrationEnergyWithRegardToPhi(const typename ConfiguratorType::InitType &Initializer,
                                                 const aol::StiffOp<ConfiguratorType> &Stiff,
                                                 const aol::Vector<RealType> &PhaseField,
                                                 const aol::Vector<RealType> &Reference,
                                                 const RealType Lambda):
  _grid(Initializer),
  _phaseField(PhaseField),
  _reference(Reference),
  _lambda(Lambda),
  _blockStiff(Stiff)
  {
  }
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest) const{

    qc::ATDeformationGradient<ConfiguratorType, ConfiguratorType::Dim> defGradOp( _grid, _phaseField , _reference );

    aol::LinCombOp<aol::MultiVector<RealType> > gradOp;

    gradOp.appendReference( defGradOp );
    //gradOp.appendReference( _blockStiff, _lambda );

    gradOp.apply( Arg, Dest );

    // Generate Identity
/*
    aol::MultiVector<RealType> identity( _dim, _grid.getNumberOfNodes() );

    generateIdentity<RealType, _dim>(_grid, identity);

    Dest *= 1./lambda;
    blockStiff.applyAdd( identity, Dest);
    Dest *= lambda;
*/
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

template <typename ConfiguratorType>
class ATRegistrationBase : public qc::RegistrationInterface<ConfiguratorType> {
protected:
  /**Type define*/
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  /**Computation parameters*/
  const RealType _epsilon;
  const RealType _alpha;
  const RealType _beta;

  /**The right hand side vector*/
  ArrayType _rhs;

  /** Predefined mass operator*/
  const aol::MassOp<ConfiguratorType> _mass;

  /** Predefined stiff operator*/
  const aol::StiffOp<ConfiguratorType> _stiff;

  const aol::LumpedMassOp<ConfiguratorType> _lumpedMassInv;
  ATRegistrationBase(const RealType Epsilon,
                     const RealType Alpha,
                     const RealType Beta,
                     const RealType Lambda,
                     const ArrayType &Reference0,
                     const ArrayType &Template0,
                     aol::OperatorType StiffOpType = aol::ONTHEFLY )
  : qc::RegistrationInterface<ConfiguratorType>( Reference0, Template0, Lambda),
    _epsilon(Epsilon),
    _alpha(Alpha),
    _beta(Beta),
    _rhs(this->_grid),
    _mass (this->_grid ),
    _stiff (this->_grid, StiffOpType),
    _lumpedMassInv( this->_grid, aol::INVERT )
  {
  }
  ///compute \phi
public:
  void updateTransformation(const aol::Vector<RealType> &Reference,
                            const aol::Vector<RealType> &PhaseField,
                            aol::MultiVector<RealType> &Phi){
    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
    ATRegistrationEnergyWithRegardToPhi<ConfiguratorType> registrationEnergy(this->_grid, PhaseField, Reference);
    qc::DisplacementLengthEnergy<ConfiguratorType> regularizationEnergy(_stiff);

    E.appendReference( registrationEnergy );
    E.appendReference( regularizationEnergy, this->_lambda );

    aol::LinCombOp<aol::MultiVector<RealType> > DE;
    qc::ATDeformationGradient<ConfiguratorType, ConfiguratorType::Dim> defGradOp( this->_grid, PhaseField , Reference );
    const aol::DiagonalBlockOp<RealType> _blockStiff(_stiff);

    DE.appendReference( defGradOp );
    DE.appendReference( _blockStiff, this->_lambda );
/*
    const RealType _beta = 1.0;
    const int _numberOfIntensityLevels = 6;
    typedef aol::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfTypeFeatureDomain;
    MIRegistrationEnergyWithRegardToPhi<ConfiguratorType, ConfTypeFeatureDomain> miEnergy(this->_grid, this->_reference0, this->_template0, _numberOfIntensityLevels, _beta);
    E.appendReference( miEnergy );
    VariationOfMIRegistrationEnergyWithRegardToPhi<ConfiguratorType> variationOfMIEnergy(this->_grid, this->_reference0, this->_template0, _numberOfIntensityLevels, _beta);
    DE.appendReference( variationOfMIEnergy );
*/
    aol::MultiVector<RealType> mtmp( this->_dim, this->_grid.getNumberOfNodes() );

    aol::GradientDescent<ConfiguratorType, aol::MultiVector<RealType> > gradientDescent_solver( this->_grid, E, DE, 1 , this->_tau);
    gradientDescent_solver.apply(Phi, mtmp);
    this->_tau = gradientDescent_solver.getStartTau();
    Phi = mtmp;
  }
};

template <typename ConfiguratorType>
class ATRegistration : public ATRegistrationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  /** The system matrix*/
  qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> _mat;

  /** The sor precondition operator*/
  aol::SSORPreconditioner<aol::Vector<RealType>, qc::FastUniformGridMatrix<RealType,ConfiguratorType::Dim> > _ssor;

  /** The CG solver*/
  aol::PCGInverse<aol::Vector<RealType> > _solver;

public:
  ATRegistration(const RealType Epsilon,
                 const RealType Alpha,
                 const RealType Beta,
                 const RealType Lambda,
                 const ArrayType &Reference0,
                 const ArrayType &Template0)
  : ATRegistrationBase<ConfiguratorType>(Epsilon, Alpha, Beta, Lambda, Reference0, Template0, aol::ASSEMBLED),
    _mat( this->_grid ),
    _ssor( _mat ),
    _solver( _mat, _ssor, 1e-10, 1000 )
  {
    _solver.setQuietMode( false );
    _solver.setStopping( aol::STOPPING_ABSOLUTE );
  }
  /// compute U_T
  void updateTemplate(aol::Vector<RealType> &Template,
                      const aol::Vector<RealType> &PhaseField){
    // \beta \int (T \phi) +  \int (v^2 \nabla T \cdot \nabla \phi) = \beta \int (T_0 \phi)
    cerr << "Assembling matrices ";
    qc::ATStiffV2Op<ConfiguratorType> atStiffV2( this->_grid, PhaseField );

    _mat.setZero();

    this->_mass.assembleAddMatrix( _mat ); // mat = beta * M
    _mat *= ( this->_beta );

    atStiffV2.assembleAddMatrix( _mat ); // mat += S[v^2]
    cerr << "done\n";

    cerr << "Assembling rhs ";
    this->_mass.apply( this->getTemplImageReference(), this->_rhs ); //rhs = beta * M T_0
    this->_rhs *= this->_beta;
    cerr << "done\n";

    _solver.apply( this->_rhs, Template );
  }
  /// compute U_R
  void updateReference(aol::Vector<RealType> &Reference,
                       const aol::Vector<RealType> &PhaseField,
                       const aol::MultiVector<RealType> &Phi){
    // \beta \int (R \phi) +  \int (v^2 \circ \Phi \nabla R \cdot \nabla \phi) = \beta \int (R_0 \phi)
    cerr << "Assembling matrices ";
    qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( this->_grid, PhaseField, Phi);

    _mat.setZero();

    this->_mass.assembleAddMatrix( _mat ); // mat = beta * M
    _mat *= ( this->_beta );

    atStiffV2Deform.assembleAddMatrix( _mat ); // mat += S[v^2 \circ \Phi]
    cerr << "done\n";

    cerr << "Assembling rhs ";
    this->_mass.apply( this->getRefImageReference(), this->_rhs ); //rhs = beta * M R_0
    this->_rhs *= this->_beta;
    cerr << "done\n";

    _solver.apply( this->_rhs, Reference );
  }
  ///compute v
  void updatePhaseField(const aol::Vector<RealType> &Reference,
                        const aol::Vector<RealType> &Template,
                        aol::Vector<RealType> &PhaseField,
                        const aol::MultiVector<RealType> &Phi){
    // \int (\|\nabla T \|^2 v \phi) + \int (\|\nabla R \|^2 (v \circ \Phi)(\phi \circ \Phi)
    // + \alpha \epsilon \int (\nabla v \cdot \nabla \phi)
    // + \frac{\alpha}{4 \epsilon} \int ( v \phi) = \frac{\alpha}{4 \epsilon} \int \phi
    cerr << "Assembling matrices ";
    ArrayType tmp( this->_grid );
    qc::ATMassGradU2Op<ConfiguratorType> atMassGradT2( this->_grid, Template );

    qc::ATGradSqrOp<ConfiguratorType> atGradSqrOp( this->_grid );
    atGradSqrOp.apply( Reference, this->_rhs ); //rhs is just used as temp vector here to store gradRSqrDeform

    //_mat.setZero();
    //_mass.assembleAddMatrix( _mat );
    //_solver.apply( _rhs, tmp );

    this->_lumpedMassInv.apply( this->_rhs, tmp );
/*
    typename ConfiguratorType::VecType grad;
    qc::FastILexMapper<ConfiguratorType::Dim> mapper(_grid);
    const aol::DiscreteFunctionDefault<ConfiguratorType> _discrReference(_grid, Reference);
    qc::GridDefinition::OldFullElementIterator it;
    for ( it = _grid.begin(); it != _grid.end(); ++it ) {
      _discrReference.evaluateGradientAtQuadPoint( *it, 0, grad );
      tmp[ mapper.localToGlobal( *it, 0) ] = grad.normSqr();
    }
*/
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( this->_grid );
    transform.setDeformation( Phi );
    transform.apply( tmp, this->_rhs );

    qc::ATMassGradU2DeformOp<ConfiguratorType> atMassGradR2Deform( this->_grid, this->_rhs, Phi );

    _mat.setZero();

    this->_stiff.assembleAddMatrix( _mat ); //mat = alpha*epsilon*S
    _mat *= ( this->_alpha * this->_epsilon );
    atMassGradT2.assembleAddMatrix( _mat ); //mat += M[GradT2]
    atMassGradR2Deform.assembleAddMatrix( _mat ); //mat += M[GradR2Phi^-1|det D Phi|^-1]

    _mat *= ( (4. * this->_epsilon) / this->_alpha); // mat += (alpha/(4*epsilon))*M
    this->_mass.assembleAddMatrix( _mat );
    _mat *= ( this->_alpha / ( 4. * this->_epsilon) );
    cerr << "done\n";

    cerr << "Assembling rhs ";
    tmp.setAll( this->_alpha / ( this->_epsilon * 4. ) ); // rhs = M tmp
    this->_mass.apply( tmp, this->_rhs );
    cerr << "done\n";

    _solver.apply( this->_rhs, PhaseField );
  }
};

template <typename ConfiguratorType>
class ATRegistrationOnTheFly: public ATRegistrationBase<ConfiguratorType> {
  /**Type define*/
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  aol::LinCombOp<aol::Vector<RealType> > _systemOp;

  /** The CG solver*/
  aol::CGInverse<aol::Vector<RealType> > _solver;

public:
  ATRegistrationOnTheFly(const RealType Epsilon,
                         const RealType Alpha,
                         const RealType Beta,
                         const RealType Lambda,
                         const qc::Array<RealType> &Reference0,
                         const qc::Array<RealType> &Template0)
  : ATRegistrationBase<ConfiguratorType>(Epsilon, Alpha, Beta, Lambda, Reference0, Template0),
    _solver( _systemOp, 1e-10, 1000 )
  {
    _solver.setQuietMode( false );
  }
  void updateTemplate(aol::Vector<RealType> &Template,
                      const aol::Vector<RealType> &PhaseField){
    // \beta \int (T \phi) +  \int (v^2 \nabla T \cdot \nabla \phi) = \beta \int (T_0 \phi)
    cerr << "Generating Op ";
    qc::ATStiffV2Op<ConfiguratorType> atStiffV2( this->_grid, PhaseField );

    _systemOp.clearAll(); // systemOp = 0
    _systemOp.appendReference( this->_mass, this->_beta ); // systemOp += beta * M
    _systemOp.appendReference( atStiffV2 ); // systemOp += S[v^2]
    cerr << "done\n";

    cerr << "Assembling rhs ";
    this->_mass.apply( this->_template0, this->_rhs ); //rhs = beta * M T_0
    this->_rhs *= this->_beta;
    cerr << "done\n";

    _solver.apply( this->_rhs, Template );
  }
  void updateReference(aol::Vector<RealType> &Reference,
                       const aol::Vector<RealType> &PhaseField,
                       const aol::MultiVector<RealType> &Phi){
    // \beta \int (R \phi) +  \int (v^2 \circ \Phi \nabla R \cdot \nabla \phi) = \beta \int (R_0 \phi)
    cerr << "Generating Op ";
    qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( this->_grid, PhaseField, Phi);

    _systemOp.clearAll(); // systemOp = 0
    _systemOp.appendReference( this->_mass, this->_beta ); // systemOp += beta * M
    _systemOp.appendReference( atStiffV2Deform ); // systemOp += S[v^2 \circ \Phi]
    cerr << "done\n";

    cerr << "Assembling rhs ";
    this->_mass.apply( this->_reference0, this->_rhs ); //rhs = beta * M R_0
    this->_rhs *= this->_beta;
    cerr << "done\n";

    _solver.apply( this->_rhs, Reference );
  }
  void updatePhaseField(const aol::Vector<RealType> &Reference,
                        const aol::Vector<RealType> &Template,
                        aol::Vector<RealType> &PhaseField,
                        const aol::MultiVector<RealType> &Phi){
    // \int (\|\nabla T \|^2 v \phi) + \int (\|\nabla R \|^2 (v \circ \Phi)(\phi \circ \Phi)
    // + \alpha \epsilon \int (\nabla v \cdot \nabla \phi)
    // + \frac{\alpha}{4 \epsilon} \int ( v \phi) = \frac{\alpha}{4 \epsilon} \int \phi
    cerr << "Generating Op ";
    ArrayType tmp( this->_grid );
    ArrayType tmp2( this->_grid );
    qc::ATMassGradU2Op<ConfiguratorType> atMassGradT2( this->_grid, Template );

    qc::ATGradSqrOp<ConfiguratorType> atGradSqrOp( this->_grid );
    atGradSqrOp.apply( Reference, this->_rhs ); //rhs is just used as temp vector here to store gradRSqrDeform
    this->_lumpedMassInv.apply( this->_rhs, tmp );

    qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( this->_grid );
    transform.setDeformation( Phi );
    transform.apply( tmp, tmp2 );

    qc::ATMassGradU2DeformOp<ConfiguratorType> atMassGradR2Deform( this->_grid, tmp2, Phi );

    _systemOp.clearAll(); // systemOp = 0
    _systemOp.appendReference( this->_stiff, this->_alpha * this->_epsilon ); // systemOp += alpha*epsilon*S
    _systemOp.appendReference( atMassGradT2 ); //systemOp += M[GradT2]
    _systemOp.appendReference( atMassGradR2Deform ); //systemOp += M[GradR2Phi^-1|det D Phi|^-1]
    _systemOp.appendReference( this->_mass, this->_alpha / ( 4. * this->_epsilon) ); // systemOp += (alpha/(4*epsilon))*M
    cerr << "done\n";

    cerr << "Assembling rhs ";
    tmp.setAll( this->_alpha / ( this->_epsilon * 4. ) ); // rhs = M tmp
    this->_mass.apply( tmp, this->_rhs );
    cerr << "done\n";

    _solver.apply( this->_rhs, PhaseField );

  }
};

template <typename ConfiguratorType>
class atRegMultilevelDescent : public qc::RegistrationMultilevelDescentInterface<ConfiguratorType>{
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
#ifdef GUI_PRESENT
  FXImageUpdater<RealType, ConfiguratorType::Dim> *_imageUpdater;
#endif
  qc::MultilevelArray<RealType> _template, _reference,
                                //_template0, _reference0,
                                _u;

public:
  atRegMultilevelDescent(
#ifdef GUI_PRESENT
                          FXImageUpdater<RealType, ConfiguratorType::Dim> *ImageUpdater,
#endif
                          const aol::ParameterParser &Parser ) :
    qc::RegistrationMultilevelDescentInterface<ConfiguratorType>(Parser),
#ifdef GUI_PRESENT
    _imageUpdater(ImageUpdater),
#endif
    _template( this->_grid ),
    _reference( this->_grid ),//_template0( _grid ), _reference0( _grid ),
    _u( this->_grid )
  {
    _template[ this->_maxDepth ] = this->_org_template[ this->_maxDepth ];
    _reference[ this->_maxDepth ] = this->_org_reference[ this->_maxDepth ];
    _template.levRestrict( 0, this->_curLevel );
    _reference.levRestrict( 0, this->_curLevel );
/*
    char fn[1024];
    for( int i = 0; i < ConfiguratorType::Dim; i++ ){
      sprintf( fn, "deformation_%d_08.dat.bz2", i);
      loadAndPrepareImage( fn, this->_transformation.getMultilevelArray(i), true );
    }
*/
  }

  virtual ~atRegMultilevelDescent( ) {
  }

  void setLevel( const int Level ) {
    this->_curLevel = Level;
    this->_curGrid = this->_grids[ this->_curLevel ];
    this->_transformation.setCurLevel( Level );
    _template.setCurLevel( Level );
    _reference.setCurLevel( Level );
    this->_org_template.setCurLevel( Level );
    this->_org_reference.setCurLevel( Level );
    //_template0.setCurLevel( Level );
    //_reference0.setCurLevel( Level );
    _u.setCurLevel( Level );
  }

  void prolongate( ) {
    if ( this->_curLevel < this->_grid.getGridDepth() ) {
      this->_transformation.levProlongate( );
      _u.levProlongate( );
      _template.levProlongate( );
      _reference.levProlongate( );
      setLevel( this->_curLevel + 1 );
    }
  }

  virtual void descentOnCurrentGrid();

  typename ConfiguratorType::RealType ElasticEnergy(const qc::GridDefinition &Grid, const aol::MultiVector<RealType> &phi, RealType lambda) const;
  typename ConfiguratorType::RealType DisplacementEnergy(const qc::GridDefinition &Grid, const aol::MultiVector<RealType> &phi, RealType lambda) const;
  typename ConfiguratorType::RealType FidelityEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &R, const aol::Vector<RealType> &R_0, const RealType beta) const;
  typename ConfiguratorType::RealType CoupledEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &R, const aol::Vector<RealType> &T, const aol::Vector<RealType> &v, const aol::MultiVector<RealType> &phi) const;
  typename ConfiguratorType::RealType NonzeroEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &v, const RealType alpha, const RealType epsilon) const;
  typename ConfiguratorType::RealType VSmoothEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &v, const RealType alpha, const RealType epsilon) const;
};

template <typename ConfiguratorType>
void atRegMultilevelDescent<ConfiguratorType>::descentOnCurrentGrid() {

  ArrayType tmp( *this->_curGrid );
  aol::MultiVector<RealType> mtmp(ConfiguratorType::Dim, this->_curGrid->getNumberOfNodes() );

  // make a multivector, referencing on the array in the multidimmultilevelarray _sol, which contains the deformation
  qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_transformation );

  const RealType alpha = this->getParserReference().getDouble( "alpha" );
  const RealType beta = this->getParserReference().getDouble( "beta" );
  const RealType lambda = this->getParserReference().getDouble( "lambda" );

  RealType temp = 0.;
  if( this->_curLevel < 9 )
    temp = this->H() * this->getParserReference().getDouble( "epsilonFactor" );
  else
    temp = this->_grids[ 8 ]->H() * this->getParserReference().getDouble( "epsilonFactor" );
  const RealType epsilon =  temp;

  //const RealType epsilon =  this->H() * this->getParserReference().getDouble( "epsilonFactor" );

  RealType normedDifferenceReference = 0.;
  RealType normedDifferenceTemplate = 0.;
  RealType normedDifferencePhaseField = 0.;
  RealType normedDifferenceTransformation = 0.;

  cerr << "\n>>>> epsilon = " << epsilon
       << " alpha = " << alpha
       << " beta = " << beta
       << " lambda = " << lambda
       << endl << endl;
/*
  qc::LinearSmoothOp<RealType> linSmooth0;
  linSmooth0.setCurrentGrid( *_curGrid );
  linSmooth0.setSigma( getParserReference().getDouble("presmoothFactor") *(*_curGrid).H() );
  linSmooth0.apply( _org_template[ _curLevel ], _template0[ _curLevel ] );
  linSmooth0.apply( _org_reference[ _curLevel ], _reference0[ _curLevel ] );

  writeImage<RealType> ( *_curGrid, _template0[ _curLevel ], "results/smoothed_template0" );
  writeImage<RealType> ( *_curGrid, _reference0[ _curLevel ], "results/smoothed_reference0" );
*/
  string filename;
/*
  filename = aol::strprintf ( "results/Eng_%02d.dat",_curLevel);
  ofstream EnergyOut(filename.c_str());
*/

  ///////////////////////////////////// BEGIN ITERATION ////////////////////////////////////

  ATRegistration<ConfiguratorType> atRegistration(epsilon, alpha, beta, lambda, this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ]);

  cerr << "\n----------------------------------------------------\n";
  cerr << "Registration on level " << this->_curLevel << " started\n";
  cerr << "----------------------------------------------------\n\n";

  this->writeInitialRegistration ( phi );

  filename = aol::strprintf ( "%sv_before_%02d", this->getSaveDirectory(), this->_curLevel);
  qc::writeImage<RealType> ( *this->_curGrid, _u[ this->_curLevel ], filename.c_str() );

  //_template[ _curLevel ].setAll(0.);  // initial value for (p) cg: crucial,
  //                                     // different initial values lead to different solutions!
  //_reference[ _curLevel ].setAll(0.);

  int iter = 0;
  //const RealType stopEpsilon = 0.0005*pow(_curLevel/8.,1./_dim);
  const RealType stopEpsilon = 0;
  do{
    // First, given v solve for T:
    cerr << "Solving for T at iteration " << iter << endl;
    tmp = _template[ this->_curLevel ];
    atRegistration.updateTemplate(_template[ this->_curLevel ], _u[ this->_curLevel ]);
    tmp -= _template[ this->_curLevel ];
    normedDifferenceTemplate = tmp.normSqr()/((_template[ this->_curLevel ]).normSqr()+0.0001);
    // write output
    filename = aol::strprintf ( "%sT_%02d_%03d", this->getSaveDirectory(), this->_curLevel , iter );
    qc::writeImage<RealType> ( *this->_curGrid, _template[ this->_curLevel ], filename.c_str() );

    // Then, given v solve for R:
    cerr << "Solving for R at iteration " << iter << endl;
    tmp = _reference[ this->_curLevel ];
    atRegistration.updateReference( _reference[ this->_curLevel ], _u[ this->_curLevel ], phi);
    tmp -= _reference[ this->_curLevel ];
    normedDifferenceReference = tmp.normSqr()/((_reference[ this->_curLevel ]).normSqr()+0.0001);
    // write output
    filename = aol::strprintf ( "%sR_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
    qc::writeImage<RealType> ( *this->_curGrid, _reference[ this->_curLevel ], filename.c_str() );

    // now solve for v:
    cerr << "Solving for v at iteration " << iter << endl;
    tmp = _u[ this->_curLevel ];
    atRegistration.updatePhaseField(_reference[ this->_curLevel ], _template[ this->_curLevel ], _u[ this->_curLevel ], phi);
    tmp -= _u[ this->_curLevel ];
    normedDifferencePhaseField = tmp.normSqr()/((_u[ this->_curLevel ]).normSqr()+0.0001);

    // write output...
    filename = aol::strprintf ( "%sv_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
    qc::writeImage<RealType> ( *this->_curGrid, _u[ this->_curLevel ], filename.c_str() );

#ifdef GUI_PRESENT
    _imageUpdater->updateImage( *this->_curGrid, _u[ this->_curLevel ], 0 );
#endif

    // :::::::::::::::::::::::::::: make update for phi ::::::::::::::::::::::::::::::::::::::::
    cerr << "Solving for phi at iteration " << iter << endl;
    mtmp = phi;
    atRegistration.updateTransformation(_reference[ this->_curLevel ], _u[ this->_curLevel ], phi);
    mtmp -= phi;
    normedDifferenceTransformation = tmp.normSqr()/(phi.normSqr()+0.0001);

    qc::TransformFunction<RealType, ConfiguratorType::Dim> transform( *this->_curGrid );
    transform.setDeformation( phi );
    transform.apply( this->_org_reference[ this->_curLevel ], tmp );

    this->writeCurrentRegistration ( phi, iter );

    cerr << normedDifferenceReference << " " << normedDifferenceTemplate << " " << normedDifferencePhaseField << " " << normedDifferenceTransformation << endl;

#ifdef GUI_PRESENT
    _imageUpdater->updateImage( this->_grid, match, 2 );
    //_imageUpdater->fillImageWithScalarArray( 1, this->_grid, match );
    _imageUpdater->fillImageWithTwoBlendedScalarArrays( 1, this->_grid, templ, *this->_curGrid, _u[ this->_curLevel ] );
    _imageUpdater->renderImage( 1 );
#endif
/*
    RealType EfidR, EfidT, Ecoup, Enonzero, Evsmooth, Eelast, Edisp, E;
    EfidR = FidelityEnergy(*_curGrid, _reference[ _curLevel ], _org_reference[ _curLevel ], beta);
    EfidT = FidelityEnergy(*_curGrid, _template[ _curLevel ], _org_template[ _curLevel ], beta);
    Ecoup = CoupledEnergy(*_curGrid, _reference[ _curLevel ], _template[ _curLevel ], _u[ _curLevel ], phi);
    Enonzero = NonzeroEnergy(*_curGrid, _u[ _curLevel ], alpha, epsilon);
    Evsmooth = VSmoothEnergy(*_curGrid, _u[ _curLevel ], alpha, epsilon);
    Eelast = ElasticEnergy(*_curGrid,phi,lambda);
    Edisp = DisplacementEnergy(*_curGrid,phi,lambda);
    //E = EfidR + EfidT + Ecoup + Enonzero + Evsmooth + Eelast;
    E = EfidR + EfidT + Ecoup + Enonzero + Evsmooth + Edisp;

    cerr << "Fidelity Energy for R :" << EfidR << endl;
    cerr << "Fidelity Energy for T :" << EfidT << endl;
    cerr << "Coupled Energy        :" << Ecoup << endl;
    cerr << "Nonzero Energy        :" << Enonzero << endl;
    cerr << "VSmooth Energy        :" << Evsmooth << endl;
    cerr << "Elastic Energy        :" << Eelast << endl;
    cerr << "Displacement Energy   :" << Edisp << endl;
    cerr << "Complete Energy       :" << E << endl;
    EnergyOut << iter << " " << E << endl;
*/
/*
    typedef ProcessDescentData<ConfiguratorType, aol::MultiVector<RealType> > ProcessType;
    ProcessType process( *_curGrid,
                         _template[ _curLevel ],
                         //_template0[ _curLevel ],
                         _org_template[ _curLevel ],
                         _reference[ _curLevel ],
                         //_reference0[ _curLevel ],
                         _u[ _curLevel ]  );

    process.postTimeStepCall( _curLevel, iter, phi
                              //, phi_descdir
                              );
*/
    iter++;
  } while ( ( (normedDifferenceReference > stopEpsilon)
              || (normedDifferenceTemplate > stopEpsilon)
              || (normedDifferencePhaseField > stopEpsilon)
              || (normedDifferenceTransformation > stopEpsilon) ) && (iter < this->getParserReference().getInt( "iterationsPerLevel" )) );
}

  /** @brief Print the Elastic energy E=lambda*||grad(Phi)||^2.
  */
template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
ElasticEnergy(const qc::GridDefinition &Grid,
              const aol::MultiVector<RealType> &phi,
              RealType lambda) const
{
  // stiff matrix definition
  aol::StiffOp<ConfiguratorType> stiff( Grid );
  aol::DiagonalBlockOp<RealType> blockStiff( stiff );

  aol::MultiVector<RealType> identity( this->_dim, (*this->_curGrid).getNumberOfNodes() );

  generateIdentity<RealType, ConfiguratorType::Dim>(this->_grid, identity);

  RealType Ereg = 2.;

  aol::MultiVector <RealType> mtmp( this->_dim, Grid.getNumberOfNodes() );

  blockStiff.apply( identity, mtmp);
  Ereg += 2. * (phi * mtmp);
  blockStiff.apply( phi, mtmp);
  Ereg += phi * mtmp;
  Ereg *= lambda/2.;

  return Ereg;
}

template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
DisplacementEnergy(const qc::GridDefinition &Grid,
                   const aol::MultiVector<RealType> &phi,
                   RealType lambda) const
{
  RealType E;

  // stiff matrix definition
  aol::StiffOp<ConfiguratorType> stiff( Grid );
  aol::DiagonalBlockOp<RealType> blockStiff( stiff );

  aol::MultiVector <RealType> mtmp( this->_dim, Grid.getNumberOfNodes() );

  blockStiff.apply( phi, mtmp);
  E = phi * mtmp;
  E *= lambda/2.;

  return E;
}

template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
FidelityEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &R, const aol::Vector<RealType> &R_0, const RealType beta) const
{
  RealType E;
  aol::MassOp<ConfiguratorType> mass( Grid );
  aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );
  mass.apply(R, tmp);
  E = R * tmp;
  E += (-2.)*(R_0 * tmp);
  mass.apply(R_0, tmp);
  E += R_0 * tmp;
  E *= beta/2.;
  return E; // E = beta/2 * (MR*R-2MR*R_0+MR_0*R_0)
}

template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
CoupledEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &R, const aol::Vector<RealType> &T, const aol::Vector<RealType> &v, const aol::MultiVector<RealType> &phi) const
{
  RealType E;

  qc::ATStiffV2Op<ConfiguratorType> atStiffV2( Grid, v ); // E = L[v^2]T*T
  aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );
  atStiffV2.apply(T, tmp);
  E = T * tmp;

  qc::ATStiffV2DeformOp<ConfiguratorType> atStiffV2Deform( Grid, v, phi); // E *= L[v^2\circ\phi]R*R
  atStiffV2Deform.apply(R, tmp);
  E += R * tmp;

  E /= 2.;

  return E;
}

template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
NonzeroEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &v, const RealType alpha, const RealType epsilon) const
{
  RealType E;
  aol::MassOp<ConfiguratorType> mass( Grid );
  aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );
  aol::Vector<RealType> onevector( Grid.getNumberOfNodes() );
  onevector.setAll( 1.);

  mass.apply(v, tmp);
  E = v * tmp;
  E += (-2.)*(onevector * tmp);
  mass.apply(onevector, tmp);
  E += onevector * tmp;
  E *= alpha/(8.*epsilon);
  return E; // E = 1/2 * (MV*V-2MV*1vec_0+M1vec_0*1vec  (* denotes inner product)
}

template <typename ConfiguratorType>
typename ConfiguratorType::RealType atRegMultilevelDescent<ConfiguratorType>::
VSmoothEnergy(const qc::GridDefinition &Grid, const aol::Vector<RealType> &v, const RealType alpha, const RealType epsilon) const
{
  RealType E;
  aol::StiffOp<ConfiguratorType> stiff( Grid );
  aol::Vector<RealType> tmp( Grid.getNumberOfNodes() );

  stiff.apply(v, tmp);
  E = v * tmp;
  E *= (alpha*epsilon)/2.;
  return E; // E = LV*V;
}

#endif
