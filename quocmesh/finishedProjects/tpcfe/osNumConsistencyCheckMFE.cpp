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

#include <scalarArray.h>
#include <indexMapper.h>
#include <FEOpInterface.h>
#include <preconditioner.h>
#include <quocMatrices.h>
#include <configurators.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


#include "tpcfe_utils.h"
#include "osTestfunction.h"


typedef double RealType;
static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
typedef tpcfe::CFEGrid < RealType, AT > CFEGridType;

typedef tpcfe::QuocCompatGridDefinition< RealType > QuocGridType;


static const RealType DC_PLUS = 1.0, DC_MINUS = 100.0;
static const int numSmooth = 3, explicitLevel = 3;
static const RealType eps = 1.0e-16;
static const RealType relax = 1.0;

static const RealType RtoLratio = 0.201;

const int GStdLevel = 8;


static const short nRods = 1;

const char* sfnmask  = "out/osNCC1Rod_L%d_MFEsoln.dat.bz2";
const char* Gsfnmask = "out/osNCC1Rod_L%d_soln.dat.bz2";


template< typename ArrayType >
void NCCSetLevelset ( ArrayType &levelset ) {
#if 0
  tpcfe::SwissCheeseInterfaceTestFunction<RealType> testFct ( 1.0f * nRods );
  testFct.setWidth ( levelset.getNumXYZ() );
  testFct.createInterface ( levelset );
#else
  qc::ShapeLevelsetGenerator<RealType>::generateRandomRodsRemoved3DRodsLevelset ( levelset, nRods, RtoLratio/nRods, 0.0 );
#endif
}

int main ( int argc, char** argv ) {
  try {

    if ( argc != 3 || ( argv[1][0] != 'S' && argv[1][0] != 'C' ) )
      cerr << "usage: " << argv[0] << " {S,C} level" << endl;

    const char doWhat = argv[1][0];

    const int level = atoi ( argv[2] );

    char solnfilename[1024];
    sprintf ( solnfilename, sfnmask, level );

    QuocGridType grid ( level );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

    NCCSetLevelset ( levelset );

    levelset.save ( "out/NCClevelset.dat.bz2", qc::PGM_DOUBLE_BINARY );

    if ( doWhat == 'S' ) {
      // do simulation

      qc::ScalarArray<RealType, qc::QC_3D> coeff ( grid );
      for ( int i = 0; i < levelset.size(); ++i )
        coeff[i] = ( levelset[i] < 0 ? DC_MINUS : DC_PLUS );

      typedef qc::MultilinFEBandMatrix<RealType,qc::QC_3D> MatrixType;
      typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 23> > ConfiguratorType;

      typedef aol::MyMFEWeightedStiffOp < ConfiguratorType >  WStiffOpType;

      qc::ScalarArray<RealType, qc::QC_3D> dwarf ( grid ), rhs ( grid ), bcond ( grid ), soln ( grid ), Lbcond ( grid );
      qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );

      for ( qc::RectangularBoundaryIterator<qc::QC_3D> bnit ( grid ); bnit.notAtEnd(); ++bnit ) {
        if ( ( *bnit ) [1] == 0 ) {
          DirichletMask.set ( *bnit, true );
          bcond.set ( *bnit, 0.0 );
        }

        if ( ( *bnit ) [1] == grid.getNumY() - 1 ) {
          DirichletMask.set ( *bnit, true );
          bcond.set ( *bnit, 1.0 );
        }
      }

      cerr << "Setting up stiffOp" << endl;
      tpcfe::FrickelInterfaceTestFunction<RealType> fITF ( DC_PLUS, DC_MINUS, levelset );
      WStiffOpType stiffOp ( grid, fITF, aol::ONTHEFLY );

      cerr << "Assembling matrix" << endl;
      MatrixType sysMat ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
      stiffOp.assembleAddMatrix ( sysMat );

      sysMat.apply ( bcond, Lbcond );
      rhs -= Lbcond;

      cerr << "Enforcing Dirichlet BCs" << endl;
      for ( int i = 0; i < DirichletMask.size(); ++i ) {
        if ( DirichletMask[i] == true ) {
          sysMat.setRowColToZero ( i );
          sysMat.set ( i, i, 1.0 );
          rhs[i] = 0;
        }
      }

      aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( sysMat );
      aol::PCGInverse< aol::Vector<RealType> > pcgSolver ( sysMat, prec, 1e-16, 2000 );
      pcgSolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
      pcgSolver.apply ( rhs, soln );

      cerr << "Total memusage: " << aol::memusage() / 1048576.0 << " MB" << endl;

      soln += bcond;

      soln.save ( solnfilename, qc::PGM_DOUBLE_BINARY );

    } else if ( doWhat == 'C' ) {

      CFEGridType GStdGrid ( GStdLevel );
      qc::ScalarArray<RealType, qc::QC_3D> GStdLevelset ( GStdGrid );

      NCCSetLevelset ( GStdLevelset );

      GStdGrid.addStructureFrom ( GStdLevelset );

      qc::AArray<RealType, qc::QC_3D> GStdCoeff ( GStdGrid );
      for ( int i = 0; i < GStdLevelset.size(); ++i )
        GStdCoeff[i] = ( GStdLevelset[i] < 0 ? DC_MINUS : DC_PLUS );

      GStdGrid.relaxedDetectAndInitVirtualNodes ( GStdCoeff, 1.0e-15, 1.0e-15 );

      qc::ScalarArray<RealType, qc::QC_3D> GStdSoln ( GStdGrid ), compSoln ( grid );
      char filename[1024];
      sprintf ( filename, Gsfnmask, GStdLevel );
      GStdSoln.load ( filename );
      compSoln.load ( solnfilename );

      RealType L1sum = 0, L2sum = 0, Linfb = 0, Linfv = 0, volSum = 0;
      aol::Vec3<RealType> H1Sum, wH1Sum;

      RealType GL1S = 0, GL2S = 0, GLinf = 0;
      aol::Vec3<RealType> GH1S, GwH1S;

      cerr << "Total memusage: " << aol::memusage() / 1048576.0 << " MB" << endl;

      aol::ProgressBar<> pb ( "Computing Difference" );
      pb.start ( GStdGrid.getNumberOfElements() );

      const RealType h = GStdGrid.H();

      const qc::GridSize<qc::QC_3D> gridSize ( GStdGrid );
      for ( CFEGridType::FullElementIterator it ( GStdGrid ); it.notAtEnd(); ++it ) {
        tpcfe::CFEElement<RealType> el ( *it, gridSize, GStdGrid.getElType ( *it ) );

        el.computeAssembleData ( GStdGrid );

        const int BO = 0; // offset from boundary (grid cells on GStdGrid, if we loop over this) for computation of errors

        if ( ! ( ( el[0] < BO ) || ( el[1] < BO ) || ( el[2] < BO ) || ( el[0] > GStdGrid.getNumX() - 1 - BO ) || ( el[1] > GStdGrid.getNumY() - 1 - BO ) || ( el[2] > GStdGrid.getNumZ() - 1 - BO ) ) ) {

          for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
            const tpcfe::CFETetra<RealType> &tetra = *tit;

            aol::Vec3<RealType> barycenterGlobalCO;
            tetra.computeGlobalCoordinateOfBarycenter ( barycenterGlobalCO, el );

            aol::Vec<4, RealType> bary;  bary.setAll( 0.25 );

            const RealType
              GStdValue = tpcfe::interpolateDataAtPositionBary ( GStdGrid, GStdSoln, el, tetra, bary ),
              compValue = compSoln.interpolate ( ( h / grid.H() ) * barycenterGlobalCO ),
              vol = tetra.getVolume() * aol::Cub ( h ),
              diffNorm = fabs( GStdValue - compValue );

            L1sum += vol * ( diffNorm );
            L2sum += vol * aol::Sqr ( diffNorm );
            Linfb = aol::Max ( Linfb, diffNorm );

            GL1S += vol * fabs ( GStdValue );
            GL2S += vol * aol::Sqr ( GStdValue );
            GLinf = aol::Max ( GLinf, fabs ( GStdValue) );

            volSum += vol;

            // more work for approximate gradients:
            aol::Vec3<RealType> verticesGlobalCO[4], vertices[4];
            RealType GStdValues[4], compValues[4];

            for ( short v = 0; v < 4; ++v ) {
              tetra.computeGlobalCoordinate ( verticesGlobalCO[v], el, v );
              vertices[v] = h * verticesGlobalCO[v];
              aol::Vec<4, RealType> bary;  bary[v] = 1.0;
              GStdValues[v] = tpcfe::interpolateDataAtPositionBary ( GStdGrid, GStdSoln, el, tetra, bary );
              compValues[v] = compSoln.interpolate ( ( h / grid.H() ) * verticesGlobalCO[v] );

              Linfv = aol::Max ( Linfv, fabs ( GStdValues[v] - compValues[v] ) ); // these are more points for evaluation of max absolute difference
              GLinf = aol::Max ( GLinf, fabs ( GStdValues[v] ) );
            }

            aol::Vec3<RealType> GStdDifference, compDifference, GStdApproxGrad, compApproxGrad;
            aol::Matrix33<RealType> dirMat, invDirMat;

            for ( short v = 0; v < 3; ++v ) {
              dirMat.setRow ( v, vertices[v] - vertices[3] );
              GStdDifference[v] = GStdValues[v] - GStdValues[3];
              compDifference[v] = compValues[v] - compValues[3];
            }
            invDirMat = dirMat.inverse();

            GStdApproxGrad = ( invDirMat * GStdDifference );
            compApproxGrad = ( invDirMat * compDifference );

            H1Sum[0] += vol * aol::Sqr ( GStdApproxGrad[0] - compApproxGrad[0] );
            H1Sum[1] += vol * aol::Sqr ( GStdApproxGrad[1] - compApproxGrad[1] );
            H1Sum[2] += vol * aol::Sqr ( GStdApproxGrad[2] - compApproxGrad[2] );

            const RealType currentCoeff = ( tetra.getSign() < 0 ? DC_MINUS : DC_PLUS );
            wH1Sum[0] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[0] - compApproxGrad[0] );
            wH1Sum[1] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[1] - compApproxGrad[1] );
            wH1Sum[2] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[2] - compApproxGrad[2] );

            GH1S[0] += vol * aol::Sqr ( GStdApproxGrad[0] );
            GH1S[1] += vol * aol::Sqr ( GStdApproxGrad[1] );
            GH1S[2] += vol * aol::Sqr ( GStdApproxGrad[2] );

            GwH1S[0] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[0] );
            GwH1S[1] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[1] );
            GwH1S[2] += currentCoeff * vol * aol::Sqr ( GStdApproxGrad[2] );

          }
        }
      }
      pb.finish();

      const RealType
        L1Norm = L1sum,
        L2Norm = sqrt ( L2sum ),
        H1Norm = sqrt ( L2sum + H1Sum[0] + H1Sum[1] + H1Sum[2] ),
        wH1Norm = sqrt ( L2sum + wH1Sum[0] + wH1Sum[1] + wH1Sum[2] ),
        GL1N = GL1S,
        GL2N = sqrt ( GL2S ),
        GH1N = sqrt ( GL2S + GH1S[0] + GH1S[1] + GH1S[2] ),
        GwH1N = sqrt ( GL2S + GwH1S[0] + GwH1S[1] + GwH1S[2] );


      cerr << "Volume " << aol::detailedFormat ( volSum )
        //            << ", Linfb difference: " << aol::detailedFormat ( Linfb )
        //            << ", Linfv difference: " << aol::detailedFormat ( Linfv )
           << ", Linf difference: " << aol::detailedFormat ( aol::Max ( Linfb, Linfv ) )
           << " L1 difference: " << aol::detailedFormat ( L1Norm )
           << " L2 difference: " << aol::detailedFormat ( L2Norm )
           << " H1 difference: " << aol::detailedFormat ( H1Norm )
           << " wH1 difference: " << aol::detailedFormat ( wH1Norm ) << endl;

      cerr << "Relative differences: " << aol::detailedFormat ( aol::Max ( Linfb, Linfv ) / GLinf )
           << " " << aol::detailedFormat ( L1Norm / GL1N )
           << " " << aol::detailedFormat ( L2Norm / GL2N )
           << " " << aol::detailedFormat ( H1Norm / GH1N )
           << " " << aol::detailedFormat ( wH1Norm / GwH1N ) << endl;
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
