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

//#define FIX_PHASE_REFERENCE
//#define INCLUDE_WEIGHTS
//#define INCLUDE_WEIGHTS2
#include <aol.h>
#include <parameterParser.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <AmbrosioTortorelli.h>
#include "ambrotelli.h"
#include "SymAtreg3d_atom.h"
#include "utilities.h"

typedef float RType;
//typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_3D, aol::GaussQuadrature<RType,qc::QC_3D,3> > ConfType;
typedef qc::QuocConfiguratorTraitMultiLin<RType, qc::QC_2D, aol::GaussQuadrature<RType,qc::QC_2D,3> > ConfType;

/**@brief The class manages the multi level scheme.*/
template <typename ConfiguratorType>
class SymATRegMultilevelDescent : public qc::RegistrationMultilevelDescentInterface<ConfiguratorType>  {
  enum MODE {
    UPDATE_TRANSFORMATIONS_ALTERNATING,
    UPDATE_TRANSFORMATIONS_SIMULTANEOUS
  };
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  typedef qc::RegistrationMultilevelDescentInterface<ConfiguratorType> SuperClass;

  typedef typename SuperClass::MultilevelArrayType MultilevelArrayType;
  typedef typename SuperClass::MultiDimMultilevelArrayType MultiDimMultilevelArrayType;
protected:
#ifdef INCLUDE_WEIGHTS
  MultilevelArrayType  _weight_reference;
  MultilevelArrayType  _weight_template;
#endif
  /// Phase fields
  MultilevelArrayType  _phase_reference;
  MultilevelArrayType  _phase_template;

  /// Smooth functions
  MultilevelArrayType  _smo_reference;
  MultilevelArrayType  _smo_template;

  /// Transformations
  // the transform such that R\circ h = T
  MultiDimMultilevelArrayType &_trans_R2T;
  // the transform such that T\circ h = R
  MultiDimMultilevelArrayType _trans_T2R;

  std::ofstream _file_RegEnergy;
  std::ofstream _file_R2T;
  std::ofstream _file_T2R;
  std::ofstream _file_tauR2T;
  std::ofstream _file_tauT2R;
  std::ofstream _file_diff_R;
  std::ofstream _file_diff_T;
  std::ofstream _file_Jacobi;
  std::ofstream _file_ATEnergy;

  //flag of output intermedinte images
  bool _flag_saveInterImage;
  int _ID;

  const MODE _mode;


public:
  SymATRegMultilevelDescent( const aol::ParameterParser &Parser ) :
    SuperClass( Parser ),
#ifdef INCLUDE_WEIGHTS
    _weight_reference(this->_grid),
    _weight_template(this->_grid),
#endif
    _phase_reference(this->_grid),
    _phase_template(this->_grid),
    _smo_reference(this->_grid),
    _smo_template(this->_grid),
    _trans_R2T(this->_transformation),
    _trans_T2R(this->_grid, ConfiguratorType::Dim),
    _flag_saveInterImage(static_cast<bool>(Parser.getInt("saveAllIntermidiateResults"))),
    _ID(0),
    _mode(UPDATE_TRANSFORMATIONS_ALTERNATING)
  {

    string filename;
    //filename = aol::strprintf ( "%sDiff_R.txt", this->getSaveDirectory());
    //_file_diff_R.open( filename.c_str() );
    //filename = aol::strprintf ( "%sDiff_T.txt", this->getSaveDirectory());
    //_file_diff_T.open( filename.c_str() );


    switch( _mode ) {
    case UPDATE_TRANSFORMATIONS_SIMULTANEOUS:
      {
        filename = aol::strprintf ( "%sEnergyRegLog.txt", this->getSaveDirectory());
        _file_RegEnergy.open( filename.c_str() );
        _file_RegEnergy << aol::strprintf ( "#ID\tE\t\tE_CC(T2R)\tE_CC(R2T)\tE_REG(T2R)\tE_REG(R2T)\tE_CON(T2R)\tE_CON(R2T)\t tau\n") << endl;
      }
      break;
    case UPDATE_TRANSFORMATIONS_ALTERNATING:
      {
        filename = aol::strprintf ( "%sEnergyLog_R2T.txt", this->getSaveDirectory());
        _file_R2T.open( filename.c_str() );
        filename = aol::strprintf ( "%sEnergyLog_T2R.txt", this->getSaveDirectory());
        _file_T2R.open( filename.c_str() );
        filename = aol::strprintf ( "%sTauLog_R2T.txt", this->getSaveDirectory());
        _file_tauR2T.open( filename.c_str() );
        filename = aol::strprintf ( "%sTauLog_T2R.txt", this->getSaveDirectory());
        _file_tauT2R.open( filename.c_str() );
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown mode", __FILE__, __LINE__ );
    }

    filename = aol::strprintf ( "%sAJacobi.txt", this->getSaveDirectory());
    _file_Jacobi.open( filename.c_str() );
    _file_Jacobi << aol::strprintf ( "#ID\tJ(T2R)1/Max\tJ(T2R)Min\tJ(R2T)1/Max\t J(R2T)Min\n") << endl;
    filename = aol::strprintf ( "%sATEnergyLog.txt", this->getSaveDirectory());
    _file_ATEnergy.open( filename.c_str() );
    _file_ATEnergy << aol::strprintf ( "#ID\tE_AT(R)\t\tE_AT(T)\n") <<endl;

    _trans_R2T.clear();
    _trans_T2R.clear();
#ifdef INCLUDE_WEIGHTS
/*
    ArrayType weight_reference_array( _weight_reference.current() );
    weight_reference_array.setAll( 1. );
    _weight_reference.levRestrict( 0, this->_curLevel );
*/
    loadAndPrepareImage( "../../../input/2D/cortex_square_pic_weight.pgm", _weight_reference);

    ArrayType weight_template_array( _weight_template.current() );
    weight_template_array.setAll( 1. );
    _weight_template.levRestrict( 0, this->_curLevel );

//    loadAndPrepareImage( "../../../input/2D/weight_test2_weight.pgm", _weight_template);
#endif
#ifdef FIX_PHASE_REFERENCE
    loadAndPrepareImage( "../../../input/2D/cortex_square_pic_phase_field.pgm", _phase_reference);
#endif
  }

  virtual ~SymATRegMultilevelDescent( ) {
  }

  /// set current level
  void setLevel( const int Level ) {
    this->_curLevel = Level;
    this->_curGrid = this->_grids[ this->_curLevel ];
    _trans_R2T.setCurLevel( Level );
    _trans_T2R.setCurLevel( Level );
    _smo_template.setCurLevel( Level );
    _smo_reference.setCurLevel( Level );
    this->_org_template.setCurLevel( Level );
    this->_org_reference.setCurLevel( Level );
    _phase_reference.setCurLevel( Level );
    _phase_template.setCurLevel( Level );
  }

  /// prolongate to next level
  void prolongate( ) {
    if ( this->_curLevel < this->_grid.getGridDepth() ) {
      _trans_R2T.levProlongate( );
      _trans_T2R.levProlongate( );
#ifndef FIX_PHASE_REFERENCE
      _phase_reference.levProlongate( );
#endif
      _phase_template.levProlongate( );
      setLevel( this->_curLevel + 1 );
    }
  }

  void writeDeformedMatch ( const aol::Vector<RealType> &Image, const char *Filename, const aol::MultiVector<RealType> &Phidofs ) const {
    ArrayType match( this->_grid );
    qc::deformImageWithCoarseDeformation<ConfiguratorType> ( Image, this->_grid, *this->_curGrid, match, Phidofs );
    qc::writeImage<RealType> ( this->_grid, match, Filename );
  }

  //! Write a checkbox image of \f$ A \f$ and \f$ B\circ\phi \f$, \f$ A \f$ and \f$ B \f$ are images on the finest level, and \f$ \phi \f$ a transformation on the current level
  void writeCheckBox ( const char *filename, const qc::Array< RealType > &ImageA, const qc::Array< RealType > &ImageB, const aol::MultiVector<RealType> &Phidofs, const int BoxWidth, const bool SliceView = false ) const {
    ArrayType checkBoxImage ( this->_grid );
    ArrayType match ( this->_grid );

    qc::deformImageWithCoarseDeformation<ConfiguratorType> ( ImageB, this->_grid, *this->_curGrid, match, Phidofs );
    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    generator.generateCheckView ( checkBoxImage, ImageA, match, BoxWidth, SliceView );

    qc::writeImage<RealType> ( this->_grid, checkBoxImage, filename );
  }

  ///loading pgm images
  void loadSliceData( );
  /// compute the transformation via gradient descent
  virtual void descentOnCurrentGrid( );
};

template <typename ConfiguratorType>
void SymATRegMultilevelDescent<ConfiguratorType>::descentOnCurrentGrid(  ) {
  // need to be implemented

  // create lumped mass inv, mass and stiff ops
  ArrayType tmp( *this->_curGrid );

  // make a multivector, referencing on the array in the multidimmultilevelarray , which contains the deformation
  aol::MultiVector<RealType> phi_T2R(0, 0),phi_R2T(0,0);
  for ( int i = 0; i < this->_dim; i++ ){
    phi_R2T.appendReference( _trans_R2T.getArray( i ) );
    phi_T2R.appendReference( _trans_T2R.getArray( i ) );
  }

  // get parameters
  const RealType alpha  = this->getParserReference().getDouble( "alpha" );
  const RealType beta   = this->getParserReference().getDouble( "beta" );
  const RealType nu     = this->getParserReference().getDouble( "nu"   );
  const RealType mu     = this->getParserReference().getDouble( "mu"   );
  const RealType kappa  = this->getParserReference().getDouble( "kappa");
  const RealType lambda = this->getParserReference().getDouble( "lambda" );
  const RealType sigma  = this->getParserReference().getDouble( "sigma" );
  //const RealType tau    = this->getParserReference().getDouble( "tau"  );
  const int numGradientSteps = this->getParserReference().getInt( "gradient_steps"  );
  //const RealType epsilon = this->H() * this->getParserReference().getDouble( "epsilonFactor");

  RealType temp = 0.;
  if( this->_curLevel < 9 )
    temp = this->H() * this->getParserReference().getDouble( "epsilonFactor" );
  else
    temp = this->_grids[ 8 ]->H() * this->getParserReference().getDouble( "epsilonFactor" );
  const RealType epsilon =  temp;


  const int checkboxWidth = this->getParserReference().getInt( "checkboxWidth");

  cerr << "\n>>>> epsilon = " << epsilon
       << " alpha = " << alpha
       << " beta = " << beta
       << " lambda = " << lambda
       << endl << endl;


  string filename;

  ///////////////////////////////////// BEGIN ITERATION ////////////////////////////////////

  SymATRegistration<ConfiguratorType> SingleLevelRegister(alpha,
                                                          beta,
                                                          nu,
                                                          mu,
                                                          kappa,
                                                          lambda,
                                                          epsilon,
                                                          sigma,
                                                          numGradientSteps,
                                                          this->_org_reference[ this->_curLevel ],
                                                          this->_org_template[ this->_curLevel ]);

  cerr << "\n\n----------------------------------------------------\n";
  cerr << "Registration on level " << this->_curLevel << " started\n";
  cerr << "----------------------------------------------------\n\n";

  qc::ATEnergyOp<ConfiguratorType> atEnergyOpReference( *this->_curGrid,
                                                    this->_org_reference[this->_curLevel],
                                                    alpha,
                                                    beta,
                                                    nu,
                                                    epsilon );

  aol::MultiVector<RealType> uvReference(0,0);
  uvReference.appendReference( _smo_reference[this->_curLevel] );
  uvReference.appendReference( _phase_reference[this->_curLevel] );

  qc::ATEnergyOp<ConfiguratorType> atEnergyOpTemplate( *this->_curGrid,
                                                   this->_org_template[this->_curLevel],
                                                   alpha,
                                                   beta,
                                                   nu,
                                                   epsilon );

  qc::ConsistencyEnergyOp<ConfiguratorType> conEnergyOpReference(*this->_curGrid,phi_R2T);
  qc::ConsistencyEnergyOp<ConfiguratorType> conEnergyOpTemplate(*this->_curGrid,phi_T2R);


  aol::MultiVector<RealType> uvTemplate(0,0);
  uvTemplate.appendReference( _smo_template[this->_curLevel] );
  uvTemplate.appendReference( _phase_template[this->_curLevel] );

  aol::Scalar<RealType> scalarTemp;

  if(this->_flag_saveInterImage)
  {
    //check box output
    filename = aol::strprintf ( "%sfine_checkview_T2R_before_%02d", this->getSaveDirectory(), this->_curLevel );
    writeCheckBox( filename.c_str(), this->_org_reference[this->_maxDepth], this->_org_template[this->_maxDepth], phi_T2R, checkboxWidth, true);

    filename = aol::strprintf ( "%sfine_checkview_R2T_before_%02d", this->getSaveDirectory(), this->_curLevel );
    writeCheckBox( filename.c_str(), this->_org_template[this->_maxDepth], this->_org_reference[this->_maxDepth], phi_R2T, checkboxWidth, true);
  }

  for ( int iter=0; iter < 10; ++iter ) {

    // First, given v solve for T:
    cerr << "Solving for R at iteration " << iter << endl;
    SingleLevelRegister.update_U(this->_org_reference[this->_curLevel],
#ifdef INCLUDE_WEIGHTS
                                 _weight_reference[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS2
                                 _weight_template[this->_curLevel],
#endif
                                 _phase_reference[this->_curLevel],
                                 _phase_template[this->_curLevel],
                                 phi_T2R,
                                 _smo_reference[this->_curLevel]);


    // Then, given v solve for T:
    cerr << "Solving for T at iteration " << iter << endl;
    SingleLevelRegister.update_U(this->_org_template[this->_curLevel],
#ifdef INCLUDE_WEIGHTS
                                 _weight_template[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS2
                                 _weight_reference[this->_curLevel],
#endif
                                 _phase_template[this->_curLevel],
                                 _phase_reference[this->_curLevel],
                                 phi_R2T,
                                 _smo_template[this->_curLevel]);

#ifndef FIX_PHASE_REFERENCE
    // now solve for Rv:
    cerr << "Solving for Rv at iteration " << iter << endl;
    SingleLevelRegister.update_V(_smo_reference[this->_curLevel],
#ifdef INCLUDE_WEIGHTS2
                                 _weight_reference[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS
                                 _weight_template[this->_curLevel],
#endif
                                 _smo_template[this->_curLevel],
                                 //_phase_template[this->_curLevel],
                                 phi_R2T,
                                 _phase_reference[this->_curLevel]);
#endif

    // now solve for Tv:
    cerr << "Solving for Tv at iteration " << iter << endl;
    SingleLevelRegister.update_V(_smo_template[this->_curLevel],
#ifdef INCLUDE_WEIGHTS2
                                 _weight_template[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS
                                 _weight_reference[this->_curLevel],
#endif
                                 _smo_reference[this->_curLevel],
                                 //_phase_reference[this->_curLevel],
                                 phi_T2R,
                                 _phase_template[this->_curLevel]);

    switch( _mode ) {
    case UPDATE_TRANSFORMATIONS_SIMULTANEOUS:
      {
        // :::::::::::::::::::::::::::: make update for phi and psi ::::::::::::::::::::::::::::::::::::::::
        cerr << "Solving for h(T2R) and g(R2T) at iteration " << iter << endl;

        SingleLevelRegister.updateBothTransformations( _smo_reference[this->_curLevel],
                                                       _smo_template[this->_curLevel],
                                                       _phase_reference[this->_curLevel],
                                                       _phase_template[this->_curLevel],
                                                       phi_T2R,
                                                       phi_R2T,
                                                       _ID,
                                                       _file_RegEnergy );

      }
      break;
    case UPDATE_TRANSFORMATIONS_ALTERNATING:
      {
        // :::::::::::::::::::::::::::: make update for phi ::::::::::::::::::::::::::::::::::::::::
        cerr << "Solving for h(T2R) at iteration " << iter << endl;

        SingleLevelRegister.update_T_AM(_smo_reference[this->_curLevel],
                                        _phase_template[this->_curLevel],
                                        phi_R2T,
#ifdef INCLUDE_WEIGHTS
                                        _weight_reference[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS2
                                        _weight_template[this->_curLevel],
#endif
                                        _ID,
                                        _file_T2R,
                                        _file_tauT2R,
                                        phi_T2R);


        // :::::::::::::::::::::::::::: make update for phi ::::::::::::::::::::::::::::::::::::::::
        cerr << "Solving for g(R2T) at iteration " << iter << endl;

        SingleLevelRegister.update_T_AM(_smo_template[this->_curLevel],
                                        _phase_reference[this->_curLevel],
                                        phi_T2R,
#ifdef INCLUDE_WEIGHTS
                                        _weight_template[this->_curLevel],
#endif
#ifdef INCLUDE_WEIGHTS2
                                        _weight_reference[this->_curLevel],
#endif
                                        _ID,
                                        _file_R2T,
                                        _file_tauR2T,
                                        phi_R2T);
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown mode", __FILE__, __LINE__ );
    }
    //evaluating the determine of jacobi
    double JacMin_T2R,JacMax_T2R,JacMin_R2T,JacMax_R2T;
    qc::JaocibiEvaluation<ConfiguratorType>(*this->_curGrid,phi_T2R, JacMin_T2R, JacMax_T2R);
    qc::JaocibiEvaluation<ConfiguratorType>(*this->_curGrid,phi_R2T, JacMin_R2T, JacMax_R2T);
    cerr<<"J(T2R)---Max: "<<JacMax_T2R<<" Min: "<<JacMin_T2R<<endl;
    cerr<<"J(R2T)---Max: "<<JacMax_R2T<<" Min: "<<JacMin_R2T<<endl;
    this->_file_Jacobi << aol::strprintf ( "%3d\t%e\t%e\t%e\t%e",_ID,1./JacMax_T2R,JacMin_T2R,1./JacMax_R2T,JacMin_R2T) <<endl;




 //   this->_file_Jacobin_R2T<<aol::strprintf ( "%3d\t%e\t%e",_ID,1./JacMax,JacMin)<<endl;
    aol::Scalar<RealType> scalarRef,scalarTem;
    atEnergyOpReference.apply( uvReference, scalarRef );
    atEnergyOpTemplate.apply( uvTemplate, scalarTem );
    this->_file_ATEnergy << aol::strprintf ( "%3d\t%e\t%e",_ID,scalarRef[0],scalarTem[0]) <<endl;

    if( _mode == UPDATE_TRANSFORMATIONS_ALTERNATING ){
      conEnergyOpReference.apply(phi_T2R,scalarTemp);
      _file_ATEnergy << "\t" << scalarTemp[0];

      conEnergyOpTemplate.apply(phi_R2T,scalarTemp);
      _file_ATEnergy << "\t" << scalarTemp[0]<<endl;
    }


    //save inter images
    if(this->_flag_saveInterImage || iter==9)
    {
      // write reconstructed images
      filename = aol::strprintf ( "%sR_%02d_%03d", this->getSaveDirectory(), this->_curLevel , iter );
      qc::writeImage<RealType> ( *this->_curGrid, _smo_reference[ this->_curLevel ], filename.c_str() );

      filename = aol::strprintf ( "%sT_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      qc::writeImage<RealType> ( *this->_curGrid, _smo_template[ this->_curLevel ], filename.c_str() );

      // write phase field
      filename = aol::strprintf ( "%sRv_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      qc::writeImage<RealType> ( *this->_curGrid, _phase_reference[ this->_curLevel ], filename.c_str() );

      filename = aol::strprintf ( "%sTv_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      qc::writeImage<RealType> ( *this->_curGrid, _phase_template[ this->_curLevel ], filename.c_str() );

      //check box output
      filename = aol::strprintf ( "%sfine_checkview_T2R_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      writeCheckBox( filename.c_str(), this->_org_reference[this->_maxDepth], this->_org_template[this->_maxDepth], phi_T2R, checkboxWidth, true);

      filename = aol::strprintf ( "%sfine_checkview_R2T_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      writeCheckBox( filename.c_str(), this->_org_template[this->_maxDepth], this->_org_reference[this->_maxDepth], phi_R2T, checkboxWidth, true);

      //deformed output
      filename = aol::strprintf ( "%sfine_deformed_T2R_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      writeDeformedMatch( this->_org_template[this->_maxDepth],filename.c_str(),phi_T2R);

      filename = aol::strprintf ( "%sfine_deformed_R2T_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      writeDeformedMatch( this->_org_reference[this->_maxDepth],filename.c_str(),phi_R2T);

      //write the difference image
      filename = aol::strprintf ( "%sfine_diff_R_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      double diffR = qc::writeDiffImage<ConfiguratorType>(this->_grid,*this->_curGrid,filename.c_str(), this->_org_reference[this->_maxDepth], this->_org_template[this->_maxDepth],phi_T2R);

      filename = aol::strprintf ( "%sfine_diff_T_%02d_%03d", this->getSaveDirectory(), this->_curLevel, iter );
      double diffT = qc::writeDiffImage<ConfiguratorType>(this->_grid,*this->_curGrid,filename.c_str(), this->_org_template[this->_maxDepth], this->_org_reference[this->_maxDepth],phi_R2T);

      this->_file_diff_R << aol::strprintf ( "%3d\t%e",_ID,diffR) << endl;

      this->_file_diff_T << aol::strprintf ( "%3d\t%e",_ID,diffT) << endl;
    }
/*
    else
    {
      qc::Array<typename ConfiguratorType::RealType> ImageTmp(this->_grid);
      double diffR = qc::ComputeDifference<ConfiguratorType>( this->_grid,*this->_curGrid, this->_org_reference[this->_maxDepth], this->_org_template[this->_maxDepth],phi_T2R,ImageTmp);
      double diffT = qc::ComputeDifference<ConfiguratorType>( this->_grid,*this->_curGrid, this->_org_template[this->_maxDepth], this->_org_reference[this->_maxDepth],phi_R2T,ImageTmp);

      this->_file_diff_R<<aol::strprintf ( "%3d\t%e",_ID,diffR)<<endl;
      this->_file_diff_T<<aol::strprintf ( "%3d\t%e",_ID,diffT)<<endl;
    }
*/
    _ID++;
  }

  filename = aol::strprintf ( "%sTransform_R2T.dat", this->getSaveDirectory() );
  qc::SaveMultiVector<ConfiguratorType>(*this->_curGrid,phi_R2T,filename.c_str() );

  filename = aol::strprintf ( "%sTransform_T2R.dat", this->getSaveDirectory() );
  qc::SaveMultiVector<ConfiguratorType>(*this->_curGrid,phi_T2R,filename.c_str() );
}

template <typename ConfiguratorType>
void SymATRegMultilevelDescent<ConfiguratorType>::loadSliceData( ) {
  char fn_reference[1024], fn_template[1024];
  this->getParserReference().getString( "reference_Fold", fn_reference );
  this->getParserReference().getString( "template_Fold", fn_template );

  /************************************ BEGIN: loading images *************************************************/
  cerr << "reading reference image from: " << fn_reference << endl;
  cerr << "reading template image from: " << fn_template << endl;

  ArrayType tem( this->_org_template.current( ) );
  ArrayType ref( this->_org_reference.current( ) );

  char fn[1024];
  sprintf(fn,"%s%%03d.pgm",fn_reference);
  ref.setQuietMode( true );
  ref.loadSlices ( fn,qc::QC_Z,0, 128);
  ref /= ref.getMaxValue();

  sprintf(fn,"%s%%03d.pgm",fn_template);
  tem.setQuietMode( true );
  tem.loadSlices ( fn,qc::QC_Z,0, 128);
  tem /= tem.getMaxValue();

  // restrice the max and min level
  this->_org_template.levRestrict( 0, this->_curLevel );
  this->_org_reference.levRestrict( 0, this->_curLevel );
}

int main( int argc, char **argv ) {

  try {
    char parameterfilename[1024];

    if ( argc > 2 ) {
      cerr << "USAGE: " << argv[0] << " <parameterfile>" << endl;
      return EXIT_FAILURE;
    }
    if ( argc == 2 ) {
      sprintf( parameterfilename, "%s",  argv[1] );
    }
    if ( argc == 1 ) {
      if( ConfType::Dim == 3)
        sprintf( parameterfilename, "ParameterFiles/SymReg_3D.par" );
      else
        sprintf( parameterfilename, "ParameterFiles/SymReg_2D.par" );
    }
    cerr << "Reading parameters from " << parameterfilename << endl;
    aol::ParameterParser parser( parameterfilename );

    aol::StopWatch watch;
    watch.start();

    SymATRegMultilevelDescent<ConfType> mld(parser);

    mld.solve();

    watch.stop();
    cerr << "duration = " << watch.elapsedCpuTimeString() << endl;
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
