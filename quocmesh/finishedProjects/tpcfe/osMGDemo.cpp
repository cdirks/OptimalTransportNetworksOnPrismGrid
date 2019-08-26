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

#if 0

// multigrid demo that outputs various intermediate results to show how smoothing and coarse grid correction perform

#define USE_ONTHEFLY_COARSENING 1

// #define COARSENING_VERBOSE 1

#include <tpCFEStandardOp.h>
#include <tpCFEMultigrid.h>

#include "affine.h"
#include "levelsets.h"


namespace tpcfe{

template< class CFEOpType >
class CFEDemoMultigrid : public CFEMultigrid<CFEOpType> {
  typedef typename CFEMultigrid<CFEOpType>::GridType GridType;
  typedef typename CFEMultigrid<CFEOpType>::RealType RealType;
  typedef typename CFEMultigrid<CFEOpType>::VectorType VectorType;
  typedef typename CFEMultigrid<CFEOpType>::MatrixType MatrixType;

public:
  CFEDemoMultigrid ( const GridType &Grid, MatrixType &fineMatrix,
                     int ExplicitLevel = 2,
                     int PreSmooth = 3, int PostSmooth = 3,
                     aol::GaussSeidelSweepingMode GSmode = aol::GAUSS_SEIDEL_FORWARD,
                     RealType Relax = 1.0,
                     mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                     RealType CycleEps = 1.0e-16,
                     int MaxCycle = 100,
                     aol::StoppingMode stop = aol::STOPPING_UNSET ) : CFEMultigrid<CFEOpType> ( Grid, fineMatrix, ExplicitLevel, PreSmooth, PostSmooth, GSmode, Relax, CycleMode, CycleEps, MaxCycle, stop ){
    // do nothing else
  }

protected:
    //! pre-smoothing
  virtual void preSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    {
      VectorType dummy ( X, aol::STRUCT_COPY ), residuum ( X, aol::STRUCT_COPY );
      this->applyOperator ( X, dummy, Level );
      residuum = Rhs;
      residuum -= dummy;
      const int
        ctr = ( 1 << ( Level - 1 ) ) + 1,
        wi  = ( 1 << ( Level     ) ) + 1;
      for ( int i = 0; i < wi ; ++i ) {
        cout << residuum.get( wi * wi * i + wi * ctr + ctr  ) << " ";
      }
      cout << "DEMO_before_presmooth_L" << Level << endl;
    }

    if ( this->_verbose < 5 ) {
      this->smooth ( Rhs, X, Level, this->_preSmoothSteps );
    } else {
      for( int smo = 0; smo < this->_preSmoothSteps; ++smo ) {
        this->smooth ( Rhs, X, Level, 1 );
        cerr << "Pre-smoothing step " << smo+1 << ": L_2 norm of residuum: " << aol::mixedFormat ( sqrt ( this->getResiduum( Rhs, X, Level ) / this->getVectorSize ( Level ) ) ) << endl;
      }
    }

    {
      VectorType dummy ( X, aol::STRUCT_COPY ), residuum ( X, aol::STRUCT_COPY );
      this->applyOperator ( X, dummy, Level );
      residuum = Rhs;
      residuum -= dummy;
      const int
        ctr = ( 1 << ( Level - 1 ) ) + 1,
        wi  = ( 1 << ( Level     ) ) + 1;
      for ( int i = 0; i < wi ; ++i ) {
        cout << residuum.get( wi * wi * i + wi * ctr + ctr  ) << " ";
      }
      cout << "DEMO_after_presmooth_L" << Level << endl;
    }

  }

  //! post-smoothing
  virtual void postSmooth ( const VectorType &Rhs, VectorType &X, const int Level ) const {
    {
      VectorType dummy ( X, aol::STRUCT_COPY ), residuum ( X, aol::STRUCT_COPY );
      this->applyOperator ( X, dummy, Level );
      residuum = Rhs;
      residuum -= dummy;
      const int
        ctr = ( 1 << ( Level - 1 ) ) + 1,
        wi  = ( 1 << ( Level     ) ) + 1;
      for ( int i = 0; i < wi ; ++i ) {
        cout << residuum.get( wi * wi * i + wi * ctr + ctr  ) << " ";
      }
      cout << "DEMO_before_postsmooth_L" << Level << endl;
    }

    if ( this->_verbose < 5 ) {
      smooth ( Rhs, X, Level, this->_postSmoothSteps );
    } else {
      for( int smo = 0; smo < this->_postSmoothSteps; ++smo ) {
        this->smooth ( Rhs, X, Level, 1 );
        cerr << "Post-smoothing step " << smo+1 << ": L_2 norm of residuum: " << aol::mixedFormat ( sqrt ( this->getResiduum( Rhs, X, Level ) / this->getVectorSize ( Level ) ) ) << endl;
      }
    }

    {
      VectorType dummy ( X, aol::STRUCT_COPY ), residuum ( X, aol::STRUCT_COPY );
      this->applyOperator ( X, dummy, Level );
      residuum = Rhs;
      residuum -= dummy;
      const int
        ctr = ( 1 << ( Level - 1 ) ) + 1,
        wi  = ( 1 << ( Level     ) ) + 1;
      for ( int i = 0; i < wi ; ++i ) {
        cout << residuum.get( wi * wi * i + wi * ctr + ctr  ) << " ";
      }
      cout << "DEMO_after_postsmooth_L" << Level << endl;
    }

  }


}; // end class


} // end namespace




typedef tpcfe::CFEGrid< double, tpcfe::CFE_CD >                              GridType;

typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEBandMatrix<GridType> >  ConfiguratorType;
typedef ConfiguratorType::MatrixType                                         MatrixType;

typedef tpcfe::CFEMassOp  < ConfiguratorType >                               MassOpType;

static const int VMODE = 0;

static const int NUM_SMOOTH = 2;
static const aol::GaussSeidelSweepingMode GSMODE = aol::GAUSS_SEIDEL_FORWARD;
static const double RELAX = 1.0;

static const int depth     = 6;
static const int pcgsteps  = 2000;
static const int mgcycles  = 1;

static const double eps  = 1e-16;


int main ( int, char** ) {
  try {

    GridType grid ( depth, qc::QC_3D );

    qc::ScalarArray<double, qc::QC_3D>
      levelset ( grid ),
      orig ( grid ),
      rhs ( grid ),
      soln ( grid );

    levelset.setAll ( -1.0 );

    // Initialize the grid
    grid.setDomainFrom ( levelset );
    grid.detectVirtualNodes();

    qc::ScalarArray<double, qc::QC_3D> coeff ( grid );
    coeff.setAll ( +1.0 );
    tpcfe::CFEWeightInterface<double> weightInterface ( coeff );
    grid.initVirtualNodes ( weightInterface );

    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    MassOpType massOp ( grid, aol::ASSEMBLED );
    MatrixType massMat ( grid );
    massMat = massOp.getMatrixRef();

    restrictNonDomainEntries ( grid, massMat , 1.0 );

    for ( int i = 0; i < grid.getWidth(); ++i ) {
      for ( int j = 0; j < grid.getWidth(); ++j ) {
        for ( int k = 0; k < grid.getWidth(); ++k ) {
          double value = ( sin ( ( 1 << (depth - 3) ) * M_PI * i / static_cast<double>( grid.getWidth() - 1 ) ) *
                           sin ( ( 1 << (depth - 3) ) * M_PI * j / static_cast<double>( grid.getWidth() - 1 ) ) *
                           sin ( ( 1 << (depth - 3) ) * M_PI * k / static_cast<double>( grid.getWidth() - 1 ) ) +
                           0.4 * ( sin ( ( 1 << (depth - 1) ) * M_PI * i / static_cast<double>( grid.getWidth() - 1 ) ) *
                                   sin ( ( 1 << (depth - 1) ) * M_PI * j / static_cast<double>( grid.getWidth() - 1 ) ) *
                                   sin ( ( 1 << (depth - 1) ) * M_PI * k / static_cast<double>( grid.getWidth() - 1 ) )   ) );
          orig.set ( i, j, k, value );
        }
      }
    }

    aol::RandomGenerator<double> rg;
    for ( int i = 0; i < soln.size(); ++i )
      soln[i] = rg.rReal();

    grid.restrictToDomain ( orig );

    massOp.apply ( orig, rhs );

#if 0
    tpcfe::CFEDemoMultigrid<MassOpType> mg_solver ( grid, massMat, 2, NUM_SMOOTH, NUM_SMOOTH, GSMODE, RELAX, mg::MULTIGRID_V_CYCLE, eps, mgcycles );
    //    demomg::CFEMultigrid<MassOpType> mg_solver ( grid, massMat, 2, NUM_SMOOTH, NUM_SMOOTH, GSMODE, RELAX, demomg::MULTIGRID_V_CYCLE, eps, mgcycles );
    //    demomg::CFEMultigrid<MassOpType> mg_solver ( grid, massMat, 2, NUM_SMOOTH, NUM_SMOOTH, RELAX, demomg::MULTIGRID_V_CYCLE, eps, mgcycles );
    mg_solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    mg_solver.setVerboseMode ( VMODE );

    mg_solver.apply ( rhs, soln );
#else

    soln.setZero();
    mg::GaussSeidelSmoother< aol::Vector<double>, MatrixType > smo( massMat, 1 );
    aol::MultiVector<double> residuae( 10, grid.getWidth() );

    for ( int step = 0; step < 10; ++step ) {
      qc::ScalarArray<double, qc::QC_3D> dummy ( grid ), residuum ( grid );
      massMat.apply ( soln, dummy );
      residuum = rhs;
      residuum -= dummy;
      const int
        ctr = grid.getWidth() / 2,
        wi  = grid.getWidth();
      for ( int i = 0; i < wi ; ++i ) {
        residuae[step][i] = residuum.get( wi * wi * i + wi * ctr + ctr  );
      }
      smo.apply ( rhs, soln );
    }

    for ( int i = 0; i < grid.getWidth(); ++i ) {
      for ( int j = 0; j < 10; ++j ) {
        cout << aol::detailedFormat ( residuae[j][i] ) << " ";
      }
      cout << endl;
    }
#endif

  } catch ( aol::Exception &ex ) {
    ex.dump();
  }

}


#else


#include<multigrid.h>
#include<smoother.h>
#include<bandMatrix.h>
#include<scalarArray.h>
#include<configurators.h>
#include<quocMatrices.h>

void gnuplotDump ( qc::ScalarArray<double, qc::QC_2D> &arr, const char* filename ) {
  ofstream out ( filename );
  const int wi = arr.getNumX();
  for ( int i = 0; i < wi-1; ++i ) {
    for ( int j = 0; j < wi-1; ++j ) {
      out << i+0 << " " << j+0 << " " << arr.get ( i+0, j+0 ) << endl;
      out << i+0 << " " << j+1 << " " << arr.get ( i+0, j+1 ) << endl;
      out << i+1 << " " << j+1 << " " << arr.get ( i+1, j+1 ) << endl;
      out << i+1 << " " << j+0 << " " << arr.get ( i+1, j+0 ) << endl;
      out << i+0 << " " << j+0 << " " << arr.get ( i+0, j+0 ) << endl << endl << endl;
    }
  }
}


typedef aol::StiffOp< qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> > > StiffOpType;
typedef aol::MassOp<  qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> > >  MassOpType;


int main ( int, char** ) {

  qc::GridDefinition grid ( 5, qc::QC_2D );
  qc::ScalarArray<double, qc::QC_2D> orig ( grid ), rhs ( grid ), soln ( grid ), iter ( grid ), resid ( grid ), error ( grid );

  StiffOpType stiffOp ( grid, aol::ONTHEFLY );
  MassOpType   massOp ( grid, aol::ONTHEFLY );

  const int wi = grid.getWidth();
  const double h = grid.H();
  for ( int i = 0; i < wi; ++i ) {
    for ( int j = 0; j < wi; ++j ) {
      //       double value = sin ( aol::NumberTrait<double>::pi * i * h ) * sin ( aol::NumberTrait<double>::pi * j * h);
      //       value += 0.4 * sin ( aol::NumberTrait<double>::pi * aol::Sqr ( 4 * i * h ) ) * sin ( aol::NumberTrait<double>::pi * aol::Sqr ( 4 * j * h) ) ;
      //       orig.set ( i, j, value );
      const double value = 1.0 + 32 * sin ( aol::NumberTrait<double>::pi * aol::Sqr ( 4 * i * h ) ) * sin ( aol::NumberTrait<double>::pi * aol::Sqr ( 4 * j * h ) );
      rhs.set( i, j, value );
    }
  }

  qc::MultilinFEBandMatrix<double,qc::QC_2D> stiffMat ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
  stiffOp.assembleAddMatrix ( stiffMat );

  qc::OTFILexMapper< qc::QC_2D > IMapper ( grid );
  for ( int i = 0; i < wi; ++i ) {
    int ind = 0;
    ind = IMapper.getGlobalIndex (    i,    0 );    stiffMat.setRowColToZero ( ind  );    stiffMat.set( ind, ind, 1.0 );
    ind = IMapper.getGlobalIndex (    i, wi-1 );    stiffMat.setRowColToZero ( ind  );    stiffMat.set( ind, ind, 1.0 );
    ind = IMapper.getGlobalIndex (    0,    i );    stiffMat.setRowColToZero ( ind  );    stiffMat.set( ind, ind, 1.0 );
    ind = IMapper.getGlobalIndex ( wi-1,    i );    stiffMat.setRowColToZero ( ind  );    stiffMat.set( ind, ind, 1.0 );
  }

  mg::GaussSeidelSmoother< aol::Vector<double>, qc::MultilinFEBandMatrix<double,qc::QC_2D> > smoother ( stiffMat, 2 );

  aol::CGInverse< aol::Vector<double> > solver ( stiffMat );
  solver.apply ( rhs, soln );
  gnuplotDump ( rhs, "./out/rhs.dat" );
  gnuplotDump ( soln, "./out/soln.dat" );


  gnuplotDump ( iter, "./out/iterate_before_smooth.dat" );

  error = soln;
  error -= iter;
  gnuplotDump ( error, "./out/error_before_smooth.dat" );

  stiffMat.apply ( iter, resid );
  resid -= rhs;

  gnuplotDump ( resid, "./out/resid_before_smooth.dat" );

  // NOW APPLY SMOOTHER
  smoother.apply ( rhs, iter );

  gnuplotDump ( iter, "./out/iterate_after_smooth.dat" );

  error = soln;
  error -= iter;
  gnuplotDump ( error, "./out/error_after_smooth.dat" );

  stiffMat.apply ( iter, resid );
  resid -= rhs;
  gnuplotDump ( resid, "./out/resid_after_smooth.dat" );


  return( EXIT_SUCCESS );
}


#endif
