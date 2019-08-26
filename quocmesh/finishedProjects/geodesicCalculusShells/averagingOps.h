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

#ifndef __AVERAGINGOPS_H
#define __AVERAGINGOPS_H

#include <aol.h>
#include <parameterParser.h>
#include <multiStreambuf.h>

#include <triangMeshConfigurators.h>
#include <triMesh.h>
#include "omp.h"

#include "deformationEnergies.h"
#include "geodCalculusOps.h"

//!===============================================================================================================================
//! ELASTIC AVERAGING
//!===============================================================================================================================

//! \brief Functional for elastic average: E[S_1, ..., S_K]( S ) := \sum_k W[ S_k, S ]
//! \author Heeren
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class ElasticAverageFunctional : public aol::Op< aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > { 
  
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::InitType MeshType;
    typedef aol::MultiVector<RealType> VectorType;
  
    const aol::ParameterParser& _pparser;
    const MeshTopologySaver<MeshType>& _topology;
    const aol::RandomAccessContainer<VectorType>& _data;
    int _numOfData;
    
  public:
    ElasticAverageFunctional( const aol::ParameterParser& pparser, 
		       const MeshTopologySaver<MeshType>& topology,
		       const aol::RandomAccessContainer<VectorType>& Data ) 
      : _pparser( pparser ), _topology( topology ), _data( Data ), _numOfData( _data.size() ){}
  
    void applyAdd( const VectorType& Arg, aol::Scalar<RealType>& Dest ) const {     
      
      aol::Vector<RealType> energies;
      evaluateEnergies( Arg, energies );
      
      for( int i = 0; i < _numOfData; i++ )
	Dest[0] += energies[i];
    }
    
    void evaluateEnergies( const VectorType& Arg, aol::Vector<RealType>& Dest ) const {
       
      Dest.reallocate( _numOfData );
      Dest.setZero();
      
#ifdef _OPENMP
#pragma omp parallel for
#endif          
      for( int i = 0; i < _numOfData; i++ ){
	aol::Scalar<RealType> temp;
        MembraneDeformationType( _topology, _pparser ).applyEnergy ( _data[i], Arg, temp );
        Dest[i] += temp[0] * _pparser.getDouble("tangWeight");
        BendingDeformationType( _topology, _pparser ).applyEnergy ( _data[i], Arg, temp );
        Dest[i] += temp[0] * _pparser.getDouble("bendWeight");      
      }
    }

};
  
//! \brief Gradient for elastic average
//! \author Heeren
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class ElasticAverageGradient : public aol::Op< aol::MultiVector<typename ConfiguratorType::RealType > >{ 
      
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::InitType MeshType;
    typedef aol::MultiVector<RealType> VectorType;
    
    const aol::ParameterParser& _pparser;
    const MeshTopologySaver<MeshType>& _topology;
    const aol::RandomAccessContainer<VectorType>& _data;
    int _numOfData;
    aol::BitVector *_bdryMask;
    bool _fixBoundary, _parallelize;
    
  public:
    ElasticAverageGradient( const aol::ParameterParser& pparser, 
			 const MeshTopologySaver<MeshType>& topology,
		     const aol::RandomAccessContainer<VectorType>& Data ) 
      : _pparser( pparser ), _topology( topology ), _data( Data ), _numOfData( _data.size() ), _bdryMask( NULL ), _fixBoundary(false), _parallelize( false ){}

    
    ~ElasticAverageGradient(){ 
      if( _fixBoundary ) 
        delete _bdryMask; 
    }
    
    void applyAdd( const VectorType& Arg, VectorType& Dest ) const {       
      
      // assemble with parallelization? (needs more storage!)
      if( _parallelize )
	assembleWithParallelization( Arg, Dest );
      else{
        for( int i = 0; i < _numOfData; i++ ){
          MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( _data[i], Arg, Dest, _pparser.getDouble("tangWeight") ) ;
          BendingDeformationType( _topology, _pparser ).applyAddDefGradient( _data[i], Arg, Dest, _pparser.getDouble("bendWeight") ) ;
        }
      }
      
      if( !_fixBoundary )
        return;
      
      Dest.setAllMasked( 0., *_bdryMask );
    }
    
    void setBoundaryMask( const aol::BitVector& mask ){
      _fixBoundary = true; 
      if ( _bdryMask )
        delete _bdryMask;
      _bdryMask = new aol::BitVector( mask );
      if( _bdryMask->size() != _topology.getNumVertices() )
        throw aol::Exception ( "AverageGradient::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
    }
    
    void useParallelization() { _parallelize = true; }
    
protected:
    void assembleWithParallelization( const VectorType& Arg, VectorType& Dest ) const {  
      
      aol::RandomAccessContainer<VectorType> gradients( _numOfData, Dest.numComponents(), _topology.getNumVertices() );
      
#ifdef _OPENMP
#pragma omp parallel for
#endif   
      for( int i = 0; i < _numOfData; i++ ){
	gradients[i].setZero();
        MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( _data[i], Arg, gradients[i], _pparser.getDouble("tangWeight") ) ;
        BendingDeformationType( _topology, _pparser ).applyAddDefGradient( _data[i], Arg, gradients[i], _pparser.getDouble("bendWeight") ) ;
      }
      
      // now collect all gradients
      for( int i = 0; i < _numOfData; i++ )
	Dest += gradients[i];      
    }
    
};
  
//! \brief Hessian for elastic average.
//! \author Heeren
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class ElasticAverageHessian : public aol::Op< aol::MultiVector<typename ConfiguratorType::RealType>, aol::SparseBlockMatrix< aol::SparseMatrix<typename ConfiguratorType::RealType> > > {
      
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::InitType MeshType;
    typedef aol::MultiVector<RealType> VectorType;
    typedef aol::SparseBlockMatrix< aol::SparseMatrix<RealType> > BlockMatrixType;
    
    const aol::ParameterParser& _pparser;
    const MeshTopologySaver<MeshType>& _topology;
    const aol::RandomAccessContainer<VectorType>& _data;
    int _numOfData;
    aol::BitVector *_bdryMask;
    bool _fixBoundary, _parallelize;
    
  public:
    ElasticAverageHessian( const aol::ParameterParser& pparser, 
			 const MeshTopologySaver<MeshType>& topology,
		     const aol::RandomAccessContainer<VectorType>& Data ) 
      : _pparser( pparser ), _topology( topology ), _data( Data ), _numOfData( _data.size() ), _bdryMask( NULL ), _fixBoundary(false), _parallelize( false ) {}

    ~ElasticAverageHessian(){ 
      if( _fixBoundary ) 
        delete _bdryMask; 
    }
    
    
    void applyAdd( const VectorType& Arg, BlockMatrixType& Dest ) const {    
      
      // assemble with parallelization? (needs more storage!)
      if( _parallelize )
	assembleWithParallelization( Arg, Dest );
      else{      
        for( int i = 0; i < _numOfData; i++ ){
          MembraneDeformationType( _topology, _pparser ).assembleAddDefHessian( _data[i], Arg, Dest, _pparser.getDouble("tangWeight") ); 
          BendingDeformationType( _topology, _pparser ).assembleAddDefHessian( _data[i], Arg, Dest, _pparser.getDouble("bendWeight") ) ;
        }      
      }

      if( !_fixBoundary )
        return;
      
      // fix boundary?
      for( int k = 0; k < _bdryMask->size(); k++ ){
        if( (*_bdryMask)[k] ){
          for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
	    Dest.getReference( i, i ).setRowColToDiagonal( k );
            Dest.getReference( i, (i+1)%3 ).setRowColToZero( k );
	    Dest.getReference( i, (i+2)%3 ).setRowColToZero( k );
	  }
        }
      }
    }
    
        
    void setBoundaryMask( const aol::BitVector& mask ){
      _fixBoundary = true; 
      if ( _bdryMask )
        delete _bdryMask;
      _bdryMask = new aol::BitVector( mask );
      if( _bdryMask->size() != _topology.getNumVertices() )
        throw aol::Exception ( "AverageHessian::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
    }
    
    void useParallelization() { _parallelize = true; }
    
protected:
    void assembleWithParallelization( const VectorType& Arg, BlockMatrixType& Dest ) const {    
      
      aol::RandomAccessContainer<BlockMatrixType> matrices( _numOfData, Dest.getNumRows(), Dest.getNumCols() );
      
#ifdef _OPENMP
#pragma omp parallel for
#endif   	  
      for( int i = 0; i < _numOfData; i++ ){
        MembraneDeformationType( _topology, _pparser ).assembleAddDefHessian( _data[i], Arg, matrices[i], _pparser.getDouble("tangWeight") ); 
        BendingDeformationType( _topology, _pparser ).assembleAddDefHessian( _data[i], Arg, matrices[i], _pparser.getDouble("bendWeight") ) ;
      }
      
      // now collect all matrices
      for( int i = 0; i < _numOfData; i++ )
	Dest.addMultiple( matrices[i], 1. );
      
    } 
};
  
  
//! \brief Operator to compute the elastic average of a set of input shapes
//! \author Heeren
//! The elastic average \bar S of input data S_1, ..., S_K is defined as
//! \bar S := argmin_S \sum_k W[ S_k, S ], where W[.,.] is an elastic deformation energy.
//! Minimization is performed by means of a second order method.
//! Boundaries can be fixed (if there are any), otherwise a Lagrange setup can be applied in the usual way.
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class ElasticAverageOp {

protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::InitType  MeshType;
  typedef typename aol::MultiVector<RealType> VectorType;

  MeshType _shell;
  const MeshTopologySaver<MeshType>* _topology;
  const aol::ParameterParser& _pparser;
  
  int _numOfInputShapes;
  aol::RandomAccessContainer<VectorType> _inputShapes; 
  
  int _solverType, _numOfDescentSteps;
  RealType _stopCriterion;
  
  char _destDirectory[1024];
  aol::DeleteFlagPointer<aol::AdditionalOutputToFile> _addOut;
  bool _quiet;
  
  aol::BitVector *_bdryMask; 
  bool _fixBoundary, _lagrangeSetup, _parallelize; 

public:
  ElasticAverageOp( const aol::ParameterParser& pparser, bool Quiet = false ) :
      _pparser( pparser ),
      _numOfInputShapes( pparser.getInt("numOfAvergageData") ),
      _inputShapes( _numOfInputShapes ),
      _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
      _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
      _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
      _quiet( Quiet ),
      _bdryMask( NULL ),
      _fixBoundary( pparser.checkAndGetBool("fixBoundary") || pparser.hasVariable("bdryMask") ),
      _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ),
      _parallelize( pparser.checkAndGetBool("parallelize") ){
	
	// read in the directory, where results are to be saved
        pparser.getString( "destDirectory", _destDirectory );
        // create the directory where we want to save the results
        aol::makeDirectory( _destDirectory );
        pparser.dumpToFile( "/parameter-dump.txt", _destDirectory );
        // if desired additionally output all console output to a log file in our results directory
        if ( pparser.checkAndGetBool( "logConsoleOutput" ) )
          _addOut.reset ( new aol::AdditionalOutputToFile ( aol::strprintf ( "%s/log.txt", _destDirectory ).c_str() ), true );

	// load data
	for( int i = 0; i < _numOfInputShapes; i++ ){
	  ostringstream dataname;
          dataname << _destDirectory << _pparser.getString("averageDataStem").c_str() << i << ".ply" << ends;
	  if(!_quiet) cerr << "Load " << dataname.str() << endl;
	  _shell.loadFromPLY( dataname.str() );
	  _shell.toVector( _inputShapes[i] );
	}
	
	// load initial shape (if available)
	if( _pparser.hasVariable("initialAverage") ){
         if(!_quiet) cerr <<"Load initial average shape " << _pparser.getString("initialAverage") << endl;
         _shell.loadFromPLY( _pparser.getString("initialAverage") );
        }

        // create topology saver
	_topology = new MeshTopologySaver<MeshType>( _shell );
	
	// boundary mask? 
	if( _fixBoundary ){
          if( pparser.hasVariable("bdryMask") ){
	    if(!_quiet) cerr << "Load boundary mask..." << endl;	
	    _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
          }
          else{
	    if(!_quiet) cerr << "Set boundary mask..." << endl;	
	    _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	    _shell.fillBoundaryMask( *_bdryMask );
	  }
	}
	
	if(!_quiet) printSettings();
      }
      
  ~ElasticAverageOp(){
    if( _fixBoundary )
      delete _bdryMask;
    delete _topology;
  }
  
  void setQuietMode( bool Quiet ){
    _quiet = Quiet;
  }

  void printSettings() const {
    cerr << endl << "===========================" << endl;
    cerr << "Solving flags:" << endl;	
    cerr << " - max number of steps = " << _numOfDescentSteps << endl;
    cerr << " - stoping criterion = " << _stopCriterion << endl;
    if( _lagrangeSetup )
      cerr << " - using Lagrange setup" << endl;
    if( _fixBoundary )
      cerr << " - fix boundary" << endl;
    cerr << "===========================" << endl << endl;
  }
  
  
  void execute( VectorType& solution ) const { 
    // initialization
    VectorType InitialGuess;
    initialize( InitialGuess );
    optimize( InitialGuess, solution );
  }
  


  // compute elastic average \bar S of input data S_1, ..., S_K via
  // \bar S := argmin_S \sum_k W[ S_k, S ], where W[.,.] is an elastic deformation energy
  void execute() const {  
    
    ElasticAverageFunctional<ConfiguratorType, MembraneDeformationType, BendingDeformationType> E( _pparser, *_topology, _inputShapes ); 
    
    // initialization
    VectorType InitialGuess;
    initialize( InitialGuess );
    
    aol::Scalar<RealType> energy;
    E.apply( InitialGuess, energy );
    if(!_quiet) cerr << "Initial energy = " << energy[0] << endl;
    
    // optimization
    VectorType solution( InitialGuess, aol::STRUCT_COPY );  
    optimize( InitialGuess, solution );
    saveSolution( solution );
    
    E.apply( solution, energy );
    if(!_quiet){
      cerr << "Final energy = " << energy[0] << endl;
    
      aol::Vector<RealType> energies;
      E.evaluateEnergies( solution, energies );
      cerr << "Single energies = ";
      for( int i = 0; i < energies.size(); i++ )
        cerr << energies[i] << ", ";
      cerr << endl;
    }

  }
  
protected:    
  // optimization
  void optimize( const VectorType& InitialGuess, VectorType& solution ) const {  
        
    ElasticAverageGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType>  dE( _pparser, *_topology, _inputShapes );
    ElasticAverageHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType>  dE2( _pparser, *_topology, _inputShapes );
    
    if( _fixBoundary ){
      dE.setBoundaryMask( *_bdryMask );
      dE2.setBoundaryMask( *_bdryMask );
    }       
    
    if( _parallelize ){
     dE.useParallelization();
     dE2.useParallelization();
    }
    
    // solving
    solution.reallocate( ConfiguratorType::Dim, _topology->getNumVertices() );
    if(!_quiet) cerr <<"Start minimizing...\n";    
    NewtonMethodSolver<ConfiguratorType>( _shell, dE, dE2, _solverType, _numOfDescentSteps, _stopCriterion, _quiet ).solve( InitialGuess, solution, _lagrangeSetup );  
  }
  
  void saveSolution( const VectorType& position ) const {
    ostringstream savename;
    savename << _destDirectory << _pparser.getString("averageName") << ends;
    cerr << "Save " << savename.str() << endl;
    MeshType mesh( _shell );
    mesh.fromVector( position );
    mesh.saveAsPLY( savename.str() );
  }
  
  void initialize( VectorType& InitialGuess ) const {
     InitialGuess.reallocate( ConfiguratorType::Dim, _topology->getNumVertices() );  
    _shell.toVector( InitialGuess );  
    
    // linear initialization ?
    if( _pparser.checkAndGetBool("initializeLinearly") ){
      if(!_quiet) cerr << "Initialize linearly..." << endl;
      InitialGuess.setZero();
      for( int n = 0; n < _numOfInputShapes; n++ )
	InitialGuess.addMultiple( _inputShapes[n], 1./_numOfInputShapes );
    }
  }

};

//!===============================================================================================================================
//! GEODESIC AVERAGING
//!===============================================================================================================================

//! \brief Functional that represents sum of geodesic distances from input shapes S_n, n = 1, ..., N to some shape S.
//! \author Heeren
//! E[S] := \sum_n E[ S_n, S ], where E[S_n, S] = \sum_k W[s_{k-1}, s_k] and s_0 = S_n and s_K = S, W[.,.] is an elastic deformation energy.
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicAverageFunctional : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> >{
  
protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename aol::MultiVector<RealType> VectorType;
  typedef typename ConfiguratorType::InitType MeshType;
  
  const aol::ParameterParser& _pparser;
  MeshTopologySaver<MeshType> _topology;
  
  aol::RandomAccessContainer<VectorType> _inputShapes; 
  int _numOfInputShapes, _numOfShells, _numOfFreeShells;
  
  RealType _tangWeight, _bendWeight;

public:
  GeodesicAverageFunctional( const aol::ParameterParser& pparser, 
			     const MeshTopologySaver<MeshType>& topology,
		             const aol::RandomAccessContainer<VectorType>& Data,
		             int numOfShells  ) 
      : _pparser( pparser ), _topology( topology ), _inputShapes( Data ), _numOfInputShapes( _inputShapes.size() ), _numOfShells( numOfShells ), _numOfFreeShells( _numOfShells - 2 ), _tangWeight( _pparser.getDouble("tangWeight")), _bendWeight(_pparser.getDouble("bendWeight")){}
    
  // Arg = {x, x_1^1, ..., x_{K-1}^1, x_1^2, ...., x_{K-1}^2, .... , x^N_{K-1} }, where N is number of input shapes    
  void applyAdd( const VectorType& Arg, aol::Scalar<RealType>& Dest ) const {
    
    int dim = ConfiguratorType::Dim;
    int size = ( 1 + _numOfInputShapes * _numOfFreeShells ) * dim;
    if( Arg.numComponents() != size )
      throw aol::Exception ( "GeodesicAverageFunctional::applyAdd(): arg has wrong size!", __FILE__, __LINE__ );
    
    // reference on average shape
    VectorType average;
    for( int j = 0; j < dim; j++ )
      average.appendReference( Arg[j] );
    
    // over all input data
#ifdef _OPENMP
#pragma omp parallel for
#endif    
    for( int n = 0; n < _numOfInputShapes; n++ ){
      //cerr << "Start for shape " << n << " of " << _numOfInputShapes << endl;
      // bring varaiables into more convenient form
      aol::RandomAccessContainer< VectorType > arg( _numOfShells );
      arg[0].appendReference( _inputShapes[n] );
      for( int i = 0; i < _numOfFreeShells; i++ )
        for( int j = 0; j < dim; j++ )
	  arg[i+1].appendReference( Arg[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
      arg[_numOfShells-1].appendReference( average );
      
      aol::Scalar<RealType> aux;
      aux.setZero();
      
      for ( int shellIdx = 1; shellIdx < _numOfShells; shellIdx++ ){  
	//cerr << shellIdx << ": " << arg[shellIdx-1].numComponents() << " and " << arg[shellIdx].numComponents() << endl;
        MembraneDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _tangWeight );
        BendingDeformationType( _topology, _pparser ).applyAddEnergy ( arg[shellIdx-1], arg[shellIdx], aux, _bendWeight );
      }        
#ifdef _OPENMP
#pragma omp critical (GeodesicAverageFunctional_applyAdd)
#endif      
      Dest += aux;
    }
  } 
  
};

//! \brief First derivative of GeodesicAverageFunctional.
//! \author Heeren
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicAverageGradient : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> >{
  
protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename aol::MultiVector<RealType> VectorType;
  typedef typename ConfiguratorType::InitType MeshType;
  
  const aol::ParameterParser& _pparser;
  MeshTopologySaver<MeshType> _topology;
  
  aol::RandomAccessContainer<VectorType> _inputShapes; 
  int _numOfInputShapes, _numOfShells, _numOfFreeShells;
  
  RealType _tangWeight, _bendWeight;
  
  aol::BitVector *_bdryMask; 
  bool _fixBoundary;

public:
  GeodesicAverageGradient( const aol::ParameterParser& pparser, 
			   const MeshTopologySaver<MeshType>& topology,
		           const aol::RandomAccessContainer<VectorType>& Data,
		           int numOfShells  ) 
      : _pparser( pparser ), _topology( topology ), _inputShapes( Data ), _numOfInputShapes( _inputShapes.size() ), _numOfShells( numOfShells ), _numOfFreeShells( _numOfShells - 2), _tangWeight( _pparser.getDouble("tangWeight")), _bendWeight(_pparser.getDouble("bendWeight")), _bdryMask(NULL), _fixBoundary(false){}
      
      
    ~GeodesicAverageGradient(){ 
      if( _fixBoundary ) 
        delete _bdryMask; 
    }
    
    void setBoundaryMask( const aol::BitVector& mask ){
      _fixBoundary = true; 
      if ( _bdryMask )
        delete _bdryMask;
      _bdryMask = new aol::BitVector( mask );
      if( _bdryMask->size() != _topology.getNumVertices() )
        throw aol::Exception ( "GeodesicAverageGradient::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
    }
      
      
  // Arg = {x, x_1^1, ..., x_{K-1}^1, x_1^2, ...., x_{K-1}^2, .... , x^N_{K-1} }, where N is number of input shapes    
  void applyAdd( const VectorType& Arg, VectorType& Dest ) const {
    
    int dim = ConfiguratorType::Dim;
    int size = ( 1 + _numOfInputShapes * _numOfFreeShells ) * dim;
    if( Arg.numComponents() != size )
      throw aol::Exception ( "GeodesicAverageGradient::applyAdd(): arg has wrong size!", __FILE__, __LINE__ );
    if( Dest.numComponents() != size )
      Dest.reallocate( size, _topology.getNumVertices() );
    
    // reference on average shape
    VectorType average, gradAverage;
    for( int j = 0; j < dim; j++ ){
      average.appendReference( Arg[j] );
      gradAverage.appendReference( Dest[j] );
    }
    
    // over all input data
    for( int n = 0; n < _numOfInputShapes; n++ ){
      // bring varaiables into more convenient form
      aol::RandomAccessContainer< VectorType > arg( _numOfShells ), dest( _numOfShells - 1 );
      arg[0].appendReference( _inputShapes[n] );      
      //cerr << "arg is number " << ( n * _numOfFreeShells + 1 )*dim << " to " << ( n * _numOfFreeShells + _numOfShells - 2 )*dim + 2 << endl;
      for( int i = 0; i < _numOfFreeShells; i++ )
        for( int j = 0; j < dim; j++ ){
	  arg[i+1].appendReference( Arg[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
	  dest[i].appendReference( Dest[ ( n * _numOfFreeShells + i + 1 )*dim + j ]);
	}
      arg[_numOfShells-1].appendReference( average );
      dest[_numOfFreeShells].appendReference( gradAverage );

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int shellIdx = 0; shellIdx < _numOfShells-1; shellIdx++ ) {
        MembraneDeformationType( _topology, _pparser ).applyAddDefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx], _tangWeight );
        BendingDeformationType( _topology, _pparser ).applyAddDefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx], _bendWeight );
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif      
      for ( int shellIdx = 1; shellIdx < _numOfShells-1; shellIdx++ ){
        MembraneDeformationType( _topology, _pparser ).applyAddUndefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx-1], _tangWeight );
        BendingDeformationType( _topology, _pparser ).applyAddUndefGradient( arg[shellIdx], arg[shellIdx+1], dest[shellIdx-1], _bendWeight );      
      }   
      
      // fix voundary ?
      if( _fixBoundary )
        for ( int shellIdx = 0; shellIdx < _numOfShells - 1; shellIdx++ )
          dest[shellIdx].setAllMasked( 0., *_bdryMask );
    }
    
  } 
  
};
      
//! \brief Second derivative of GeodesicAverageFunctional.
//! \author Heeren
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicAverageHessian : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::SparseBlockMatrix<aol::SparseMatrix<typename ConfiguratorType::RealType> > >{
  
protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename aol::MultiVector<RealType> VectorType;
  typedef typename ConfiguratorType::InitType MeshType;
  
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef aol::SparseBlockMatrix< MatrixType> BlockMatrixType;
  
  
  const aol::ParameterParser& _pparser;
  MeshTopologySaver<MeshType> _topology;
  
  aol::RandomAccessContainer<VectorType> _inputShapes; 
  int _numOfInputShapes, _numOfShells, _numOfFreeShells;
  
  RealType _tangWeight, _bendWeight;
  
  aol::BitVector *_bdryMask; 
  bool _fixBoundary;

public:
  GeodesicAverageHessian( const aol::ParameterParser& pparser, 
			   const MeshTopologySaver<MeshType>& topology,
		           const aol::RandomAccessContainer<VectorType>& Data,
		           int numOfShells  ) 
      : _pparser( pparser ), _topology( topology ), _inputShapes( Data ), _numOfInputShapes( _inputShapes.size() ), _numOfShells( numOfShells ), _numOfFreeShells( _numOfShells - 2 ), _tangWeight( _pparser.getDouble("tangWeight")), _bendWeight(_pparser.getDouble("bendWeight")), _bdryMask(NULL), _fixBoundary(false){}
      
      
    ~GeodesicAverageHessian(){ 
      if( _fixBoundary ) 
        delete _bdryMask; 
    }
    
    void setBoundaryMask( const aol::BitVector& mask ){
      _fixBoundary = true; 
      if ( _bdryMask )
        delete _bdryMask;
      _bdryMask = new aol::BitVector( mask );
      if( _bdryMask->size() != _topology.getNumVertices() )
        throw aol::Exception ( "GeodesicAverageHessian::setBoundaryMask(): mask has wrong size!", __FILE__, __LINE__ );
    }
      
      
  // Arg = {x, x_1^1, ..., x_{K-1}^1, x_1^2, ...., x_{K-1}^2, .... , x^N_{K-1} }, where N is number of input shapes    
  void applyAdd( const VectorType& Arg, BlockMatrixType& Dest ) const {
      
    
    int dim = ConfiguratorType::Dim;
    int size = ( 1 + _numOfInputShapes * _numOfFreeShells ) * dim;
    if( Arg.numComponents() != size )
      throw aol::Exception ( "GeodesicAverageHessian::applyAdd(): arg has wrong size!", __FILE__, __LINE__ );
    if( (Dest.getNumCols() != size) || (Dest.getNumRows() != size) )
      throw aol::Exception ( "GeodesicAverageHessian::applyAdd(): blockmatrix has wrong size!", __FILE__, __LINE__ );
    
    // reference on average shape
    VectorType average;
    for( int j = 0; j < dim; j++ )
      average.appendReference( Arg[j] );
    
    // over all input data
    for( int n = 0; n < _numOfInputShapes; n++ ){
      // bring varaiables into more convenient form
      aol::RandomAccessContainer< VectorType > arg( _numOfShells );
      arg[0].appendReference( _inputShapes[n] );
      for( int i = 0; i < _numOfFreeShells; i++ )
        for( int j = 0; j < dim; j++ )
	  arg[i+1].appendReference( Arg[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
      arg[_numOfShells-1].appendReference( average );

      // deformed Hessian (average DOES appear here!)
      for ( int shellIdx = 0; shellIdx < _numOfShells-1; shellIdx++ ) {      
        std::vector< std::vector< std::vector< int > > > blockPositions( dim );   
	if( shellIdx != _numOfFreeShells )
          getBlockPositions( blockPositions, n, shellIdx, 0 );
	else
	  getBlockPositionsAverage( blockPositions );
        BlockMatrixType dest(dim, dim);
        Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

        MembraneDeformationType( _topology, _pparser ).assembleAddDefHessian( arg[shellIdx], arg[shellIdx+1], dest, _tangWeight );
        BendingDeformationType( _topology, _pparser ).assembleAddDefHessian( arg[shellIdx], arg[shellIdx+1], dest, _bendWeight );
	
	// fix boundary?
	if( _fixBoundary )    
          fixBoundary( dest );
   
      }
   
      // undeformed Hessian (average does NOT appear here!)
      for ( int shellIdx = 1; shellIdx < _numOfShells-1; shellIdx++ ){      
        std::vector< std::vector< std::vector< int > > > blockPositions( dim ); 
        getBlockPositions( blockPositions, n, shellIdx-1, 0 );

        BlockMatrixType dest(dim, dim);
        Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

        MembraneDeformationType( _topology, _pparser ).assembleAddUndefHessian( arg[shellIdx], arg[shellIdx+1],  dest, _tangWeight );
        BendingDeformationType( _topology, _pparser ).assembleAddUndefHessian( arg[shellIdx], arg[shellIdx+1], dest, _bendWeight );     
	
	// fix boundary?
	if( _fixBoundary )    
          fixBoundary( dest );
      }
     
     // mixed second derivatives (average DOES appear here!)
      for ( int shellIdx = 1; shellIdx < _numOfShells-1; shellIdx++ ){
        std::vector< std::vector< std::vector< int > > > blockPositions( dim ); 
        BlockMatrixType dest(dim, dim);      
	if( shellIdx != _numOfShells-2 )
          getBlockPositions( blockPositions, n, shellIdx-1, 1 );      
	else
	  getBlockPositionsAverage( blockPositions, n, -1 );
        Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

        MembraneDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, true, _tangWeight );
        BendingDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, true, _bendWeight );   
	
	// fix boundary?
	if( _fixBoundary )    
          fixBoundary( dest );
           
        //! TODO use symmetry instead of computing again!
	if( shellIdx != _numOfShells-2 )
          getBlockPositions( blockPositions, n, shellIdx-1, -1 );     
	else
	  getBlockPositionsAverage( blockPositions, n, 1 );
        Dest.getSubBlockMatrix(dim, dim, blockPositions, dest );  

        MembraneDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, false, _tangWeight );
        BendingDeformationType( _topology, _pparser ).assembleAddMixedHessian( arg[shellIdx], arg[shellIdx+1], dest, false, _bendWeight ); 
	
	// fix boundary?
	if( _fixBoundary )    
          fixBoundary( dest );
      }
    }
  } 
  
protected:
  void getBlockPositions( std::vector< std::vector< std::vector< int > > >& blockPositions, int n, int k, int off = 0 ) const {
    const int dim =ConfiguratorType::Dim;
    int colShift = off > 0 ? off : 0;
    int rowShift = off < 0 ? off : 0;
    for ( int i = 0; i < dim; ++i ){
      blockPositions[i].resize ( dim );
      for ( int j = 0; j < dim; ++j ) {
        blockPositions[i][j].resize ( 2 );
        blockPositions[i][j][0] = (1 + n * _numOfFreeShells + k - rowShift) * dim + i;
        blockPositions[i][j][1] = (1 + n * _numOfFreeShells + k + colShift) * dim + j;
      }
    }
  }
  
  void getBlockPositionsAverage( std::vector< std::vector< std::vector< int > > >& blockPositions ) const {
    const int dim =ConfiguratorType::Dim;
    for ( int i = 0; i < dim; ++i ){
      blockPositions[i].resize ( dim );
      for ( int j = 0; j < dim; ++j ) {
        blockPositions[i][j].resize ( 2 );
        blockPositions[i][j][0] = i;
        blockPositions[i][j][1] = j;
      }
    }
  }
  
  void getBlockPositionsAverage( std::vector< std::vector< std::vector< int > > >& blockPositions, int n, int off ) const {
    const int dim =ConfiguratorType::Dim;
    int colOffset = off > 0 ? (n+1) * _numOfFreeShells * dim : 0;
    int rowOffset = off < 0 ? (n+1) * _numOfFreeShells * dim : 0;
    for ( int i = 0; i < dim; ++i ){
      blockPositions[i].resize ( dim );
      for ( int j = 0; j < dim; ++j ) {
        blockPositions[i][j].resize ( 2 );
        blockPositions[i][j][0] = rowOffset + i;
        blockPositions[i][j][1] = colOffset + j;
      }
    }
  }
  
  void fixBoundary( BlockMatrixType& Matrix ) const {
    for( int k = 0; k < _bdryMask->size(); k++ ){
      if( (*_bdryMask)[k] ){
        for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
	  Matrix.getReference( i, i ).setRowColToDiagonal( k );
          Matrix.getReference( i, (i+1)%3 ).setRowColToZero( k );
	  Matrix.getReference( i, (i+2)%3 ).setRowColToZero( k );
	}
      }
    }
  }
  
};

//! \brief operator to compute the geodesic average of a set of input shapes
//! \author Heeren
//!
//! The geodesic average \bar S of input data S_1, ..., S_K is defined as
//! \bar S := argmin_S sum_n E[ S_n, S ], where E[S_n, S] = \sum_k W[s_{k-1}, s_k] and s_0 = S_n and s_K = S, W[.,.] is an elastic deformation energy.
//! Minimization is performed by means of a second order method.
//! Boundaries can be fixed (if there are any), otherwise a Lagrange setup can be applied in the usual way. 
template < typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType >
class GeodesicAverageOp {

protected:
  typedef typename ConfiguratorType::RealType  RealType;
  typedef typename ConfiguratorType::InitType  MeshType;
  typedef typename aol::MultiVector<RealType> VectorType;
  
  typedef aol::SparseMatrix<RealType> MatrixType;
  typedef aol::SparseBlockMatrix< MatrixType > BlockMatrixType;

  mutable MeshType _shell;
  const MeshTopologySaver<MeshType>* _topology;
  const aol::ParameterParser& _pparser;
  
  int _numOfInputShapes;
  aol::RandomAccessContainer<VectorType> _inputShapes; 
  int _numOfCascadicLevels, _numOfShells, _numOfFreeShells;
  
  int _solverType, _numOfDescentSteps;
  RealType _stopCriterion;
  
  char _destDirectory[1024];
  aol::DeleteFlagPointer<aol::AdditionalOutputToFile> _addOut;
  bool _quiet;
  
  aol::BitVector *_bdryMask; 
  bool _fixBoundary;
  bool _lagrangeSetup;
  
  const static int _NumOfSingleConstraints = 6;
  
  public:
  GeodesicAverageOp( const aol::ParameterParser& pparser ) :
      _pparser( pparser ),
      _numOfInputShapes( pparser.getInt("numOfAvergageData") ),
      _inputShapes( _numOfInputShapes ),
      _numOfCascadicLevels( pparser.getIntOrDefault("cascadicLevelsInAveraging", 1) ),
      _numOfShells( pow( 2, _numOfCascadicLevels ) + 1 ), 
      _numOfFreeShells( _numOfShells - 2 ),
      _solverType( _pparser.getIntOrDefault( "solverType", LU) ),
      _numOfDescentSteps( _pparser.getIntOrDefault("numOfDescentSteps",1000) ),
      _stopCriterion( _pparser.getDoubleOrDefault("stopCriterion", 1e-8) ),
      _quiet( !pparser.checkAndGetBool("showConsoleOutput") ),
      _bdryMask( NULL ),
      _fixBoundary( pparser.checkAndGetBool("fixBoundary") || pparser.hasVariable("bdryMask") ),
      _lagrangeSetup( _fixBoundary ? false : pparser.checkAndGetBool("LagrangeSetup") ){
	
	// read in the directory, where results are to be saved
        pparser.getString( "destDirectory", _destDirectory );
        // create the directory where we want to save the results
        aol::makeDirectory( _destDirectory );
        pparser.dumpToFile( "/parameter-dump.txt", _destDirectory );
        // if desired additionally output all console output to a log file in our results directory
        if ( pparser.checkAndGetBool( "logConsoleOutput" ) )
          _addOut.reset ( new aol::AdditionalOutputToFile ( aol::strprintf ( "%s/log.txt", _destDirectory ).c_str() ), true );
	
	// load data
	for( int i = 0; i < _numOfInputShapes; i++ ){
	  ostringstream dataname;
          dataname << _destDirectory << _pparser.getString("averageDataStem").c_str() << i << ".ply" << ends;
	  cerr << "Load " << dataname.str() << endl;
	  _shell.loadFromPLY( dataname.str() );
	  _shell.toVector( _inputShapes[i] );
	}
	
	// create topology saver
	_topology = new MeshTopologySaver<MeshType>( _shell );
	
	// boundary mask? 
	if( _fixBoundary ){
          if( pparser.hasVariable("bdryMask") ){
	    cerr << "Load boundary mask..." << endl;	
	    _bdryMask = new aol::BitVector( pparser.getString("bdryMask").c_str() );	
          }
          else{
	    cerr << "Set boundary mask..." << endl;	
	    _bdryMask = new aol::BitVector( _topology->getNumVertices() );
	    _shell.fillBoundaryMask( *_bdryMask );
	  }
	}
      }
      
  ~GeodesicAverageOp(){
    if( _fixBoundary )
      delete _bdryMask;
    delete _topology;
  }
  
  void setQuietMode( bool Quiet ){
    _quiet = Quiet;
  }


  void execute() const {  
    
    int dim = ConfiguratorType::Dim;
    int size = (1 + _numOfFreeShells * _numOfInputShapes ) * dim;
    VectorType Solution( size, _topology->getNumVertices() );
   
    // initialization
    if( !_quiet ) cerr << "\n=============================================" << endl;
    if( !_quiet ) cerr << "Compute elastic average as initialization..." << endl;
    VectorType elasticAverage;
    ElasticAverageOp<ConfiguratorType, MembraneDeformationType, BendingDeformationType>( _pparser, true ).execute( elasticAverage );
    
    // saving elastic average
    ostringstream AverageName;
    AverageName << _destDirectory << _pparser.getString("averageName");
    ostringstream elasticAverageName;
    string tempStr( AverageName.str() );
    elasticAverageName << tempStr.substr( 0, tempStr.rfind(".") ) << "_elastic.ply";
    _shell.fromVector( elasticAverage );
    _shell.saveAsPLY( elasticAverageName.str() );
    
    // average is first DOF in vector
    for( int j = 0; j < dim; j++ )
      Solution[j] = elasticAverage[j];   
    
    // initialize 
    bool hasToBeInitialized = true;
    if( _pparser.hasVariable( "averagePathInitialStem" ) ){
      if( !_quiet ) cerr << endl << "Start initialization for averaging... "  << endl;
      for( int n = 0; n < _numOfInputShapes; n++ ){      
        for( int i = 0; i < _numOfFreeShells; i++ ){
	  VectorType shape;
          for( int j = 0; j < dim; j++ )
	    shape.appendReference( Solution[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
	  ostringstream loadname;
	  loadname << _destDirectory << _pparser.getString( "averagePathInitialStem" ) << n << "_" << i+1 << ".ply" << ends;
	  if( !_quiet ) cerr << "Load " << loadname.str() << endl;
	  _shell.loadFromPLY( loadname.str() );
	  _shell.toVector( shape );
        }       
      }  
      hasToBeInitialized  = false;
    }
   
    // multilevel
    int startLevel = hasToBeInitialized ? 1 : _numOfCascadicLevels;
    if( !_quiet ) cerr << endl << "=============================================" << endl;
    if( !_quiet ) cerr << "Start cascadic optimization on " << _numOfCascadicLevels << " levels :" << endl;  
    if( !_quiet ) cerr << "=============================================" << endl;
    for( int l = startLevel; l <= _numOfCascadicLevels; l++ ){
      
      if( !_quiet ) cerr << endl << "=============================================" << endl;
      if( !_quiet ) cerr << "Level " << l << ": " << endl;
      
      int numOfShells = pow(2, l) + 1;
      int numOfFreeShells = numOfShells - 2;
      
      // get average
      VectorType SolutionOnThisLevel;
      for( int j = 0; j < dim; j++ )
	SolutionOnThisLevel.appendReference( Solution[j] );

      // get geodesic paths
      for( int n = 0; n < _numOfInputShapes; n++ )
        for( int i = 0; i < numOfFreeShells; i++ )
          for( int j = 0; j < dim; j++ )
	    SolutionOnThisLevel.appendReference( Solution[ ( n * _numOfFreeShells + i * pow(2,_numOfCascadicLevels - l) + ( _numOfCascadicLevels > l ? pow(2,_numOfCascadicLevels - l - 1) : 0) + 1 ) * dim + j  ] );
      
      // compute separate geodesics as initializations  
      if( !_quiet ) cerr << endl << "-----------------------------------------" << endl;
      if( !_quiet )  cerr << "Compute " << _numOfInputShapes << " geodesics of length " << numOfShells << " as initializations :" << endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for( int n = 0; n < _numOfInputShapes; n++ ){
        //if( !_quiet ) cerr << "Compute geodesic of length " << numOfShells << " as initialization of path "<< n << ":" << endl;
        VectorType Path;
        for( int i = 0; i < numOfFreeShells; i++ )
          for( int j = 0; j < dim; j++ )
	    Path.appendReference( SolutionOnThisLevel[ ( n * numOfFreeShells + i + 1 )*dim + j ] );
        GeodesicOp<ConfiguratorType, MembraneDeformationType, BendingDeformationType>( _pparser, *_topology, true ).computeGeodesicRedBlack( _inputShapes[n], elasticAverage, numOfShells, Path, hasToBeInitialized ); 
      }    
      
      // now optimize over all shapes simultaneously
      if( !_quiet ) cerr << endl << "------------------------------------------" << endl;
      if( !_quiet ) cerr << "Start computation for geodesic average..." << endl;
      optimize( elasticAverage, SolutionOnThisLevel, numOfShells );
      
      // no need for prolongatin at end
      if( l == _numOfCascadicLevels )
	continue;
      
      // prolongation to finer level
      if( !_quiet ) cerr << "Prolongation..." << endl;
      for( int n = 0; n < _numOfInputShapes; n++ ){
	// full path on this level
	VectorType Path;
	for( int j = 0; j < dim; j++ )
	  Path.appendReference( _inputShapes[n][j] );
	for( int i = 0; i < numOfFreeShells; i++ )
          for( int j = 0; j < dim; j++ )
	    Path.appendReference(  SolutionOnThisLevel[ ( n * numOfFreeShells + i + 1 )*dim + j ] );
	
	// prolongate
        for( int i = 0; i < numOfShells - 1; i++ )
	  for( int k = 1; k < pow(2,_numOfCascadicLevels - l); k++ )
            for( int j = 0; j < dim; j++ )
	      Solution[ ( n * _numOfFreeShells + i * pow(2,_numOfCascadicLevels - l) + k ) * dim + j  ] = Path[i*dim + j];
      }
      
      // now we have an initialization!
      hasToBeInitialized = false;
    }

    if( !_quiet ) cerr << endl << "=============================================" << endl;
    if( !_quiet ) cerr << endl << "Saving..." << endl;
    // saving geodesic average
    VectorType average;    
    for( int j = 0; j < dim; j++ )
      average.appendReference( Solution[j] );
    _shell.fromVector( average );
    _shell.saveAsPLY( elasticAverageName.str() );
    
    // saving paths
    savePaths( average, Solution );

  }
  
protected:
  void optimize( const VectorType& elasticAverage, VectorType& Solution, int numOfShells ) const {
    
    int numOfFreeShells = numOfShells - 2;
        
    GeodesicAverageFunctional<ConfiguratorType, MembraneDeformationType, BendingDeformationType> E( _pparser, *_topology, _inputShapes, numOfShells );
    GeodesicAverageGradient<ConfiguratorType, MembraneDeformationType, BendingDeformationType> dE( _pparser, *_topology, _inputShapes, numOfShells );
    GeodesicAverageHessian<ConfiguratorType, MembraneDeformationType, BendingDeformationType> dE2( _pparser, *_topology, _inputShapes, numOfShells );    
    
    // fix boundary?
    if( _fixBoundary ){
      dE.setBoundaryMask( *_bdryMask );
      dE2.setBoundaryMask( *_bdryMask );
    }
    
    aol::Scalar<RealType> energy;
    E.apply( Solution, energy );
    if( !_quiet ) cerr << "Initial energy = " <<  energy[0] << endl;
 
    // rigid body motion handler
    int NumOfSingleConstraints = _lagrangeSetup ? _NumOfSingleConstraints : 0;
    int NumOfConstraints = _lagrangeSetup ? (1 + numOfFreeShells * _numOfInputShapes) * NumOfSingleConstraints : 0;
    typedef MultipleRigidBodyMotionsConstraintHandler<ConfiguratorType> RBMHandlerType;
    RBMHandlerType CHandler( *_topology, elasticAverage, NumOfSingleConstraints, 1 + numOfFreeShells * _numOfInputShapes );

    GenericLagrangeGradient< ConfiguratorType, RBMHandlerType > dL( dE, CHandler );
    GenericLagrangeHessian< ConfiguratorType, RBMHandlerType > dL2( dE2, CHandler );
   
    // append Lagrange multipliers
    VectorType Destination;           
    VectorType lagrangeMult( NumOfConstraints, 1 );
    Destination.appendReference( Solution );
    if( _lagrangeSetup ) Destination.appendReference( lagrangeMult );
    VectorType Argument( Destination );
    
    // initialize solver and solve
    RealType stopCriterion = max( _stopCriterion, 1e-8 * sqrt( _topology->getNumVertices() * (1 + numOfFreeShells * _numOfInputShapes) ) );
    if( !_quiet ) cerr <<"Start solving with stopping criterion " << stopCriterion << " ... " << endl;
    NewtonMethod<ConfiguratorType, VectorType, BlockMatrixType> zeroFinder( *_topology, dL, dL2, ConfiguratorType::Dim * (1 + numOfFreeShells * _numOfInputShapes), NumOfConstraints, _solverType, _numOfDescentSteps, stopCriterion );
    zeroFinder.setTimestepController( aol::NewtonInfo<RealType>::NEWTON_OPTIMAL );
    zeroFinder.setSigma( _pparser.getDoubleOrDefault("NewtonSigma", 0.1) ); 
    zeroFinder.setQuietMode( _quiet );      
    zeroFinder.apply( Argument, Destination ); 
    
    if( !_quiet && _lagrangeSetup ) cerr << "Constraint check: " <<  CHandler.checkConstraints( Solution ) << endl;   
    
    E.apply( Solution, energy );
    if( !_quiet ) cerr << "Final energy = " <<  energy[0] << endl;
  }
  
  void savePaths( const VectorType& Average, const VectorType& Solution ) const {
    int dim = ConfiguratorType::Dim;
    for( int n = 0; n < _numOfInputShapes; n++ ){      
      // first shape is input data
      _shell.fromVector( _inputShapes[n] );
      ostringstream firstname;
      firstname << _destDirectory << _pparser.getString("averagePathStem") << n << "_" << 0 << ".ply"<< ends;
      _shell.saveAsPLY( firstname.str() );
      
      VectorType fullPath;
      for( int i = 0; i < _numOfFreeShells; i++ ){	
	VectorType temp;
        for( int j = 0; j < dim; j++ ){
	  temp.appendReference( Solution[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
	  fullPath.appendReference( Solution[ ( n * _numOfFreeShells + i + 1 )*dim + j ] );
	}
	_shell.fromVector( temp );
	ostringstream savename;
	savename << _destDirectory << _pparser.getString("averagePathStem") << n << "_" << i+1 << ".ply" << ends;
	_shell.saveAsPLY( savename.str() );
      }
      // last shell is average
      _shell.fromVector( Average );
      ostringstream lastname;
      lastname << _destDirectory << _pparser.getString("averagePathStem") << n << "_" << _numOfFreeShells+1 << ".ply"<< ends;
      _shell.saveAsPLY( lastname.str() );
      
      aol::Scalar<RealType> aux;
      GeodesicEnergy<ConfiguratorType, MembraneDeformationType, BendingDeformationType>( *_topology, _pparser, _inputShapes[n], Average, _numOfShells ).apply( fullPath, aux );
      cerr << "Path energy of " << n << "th path is " << aux[0] * (_numOfShells-1) << endl;
    } 
  }
};

#endif