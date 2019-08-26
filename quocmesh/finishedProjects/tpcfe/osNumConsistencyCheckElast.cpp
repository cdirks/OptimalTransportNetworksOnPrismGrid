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
#include <preconditioner.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEUtils.h>

#include "tpcfe_utils.h"
#include "osTestfunction.h"


typedef double RealType;
static const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;
typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;


static const int numSmooth = 3, explicitLevel = 3;
static const RealType eps = 1.0e-16;
static const RealType relax = 1.0;

const int GStdLevel = 8;

const RealType
  EMinus  = 5.0,
  nuMinus = 0.2,
  EPlus   = 1.0,
  nuPlus  = 0.2;


static const short nRods = 3;

const char* sfnmask = "out/osNCCE3Rod_L%d_soln%%d.dat.bz2";


template< typename ArrayType >
void NCCSetLevelset ( ArrayType &levelset ) {
  qc::ShapeLevelsetGenerator<RealType>::generate3DRodsLevelset ( levelset, nRods );
}


int main ( int argc, char** argv ) {
  try {

    if ( argc != 3 || ( argv[1][0] != 'S' && argv[1][0] != 'C' ) )
      cerr << "usage: " << argv[0] << " {S,C} level" << endl;

    const char doWhat = argv[1][0];

    const int level = atoi ( argv[2] );
    GridType grid ( level );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

    NCCSetLevelset ( levelset );

    levelset.save ( "out/NCClevelset.dat.bz2", qc::PGM_DOUBLE_BINARY );

    grid.addStructureFrom ( levelset );

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
    cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;
    const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );
    tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
    grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-14, 2.0e-15 );

    char solnfilename[1024];
    sprintf ( solnfilename, sfnmask, level );

    if ( doWhat == 'S' ) {
      // do simulation

      typedef tpcfe::CFEHybridMatrix<GridType>                 MatrixType;
      typedef tpcfe::CFEConfigurator < GridType, MatrixType >  ConfiguratorType;
      typedef tpcfe::CFEMassOp < ConfiguratorType >            MassOpType;
      typedef tpcfe::CFEStiffOpWI < ConfiguratorType >         WStiffOpType;

      qc::MultiArray< RealType, qc::QC_3D > DirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );
      qc::BitArray<qc::QC_3D> DirichletMask ( grid );

      tpcfe::setGeneralShearing ( grid, DirichletBCs, DirichletMask, 2, 2, 1.0 );

      grid.setDirichletMask ( DirichletMask );

      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      ElastOpType elastOp ( grid, coeff );

      dwarf -= DirichletBCs;

      elastOp.apply ( dwarf, rhs );

      elastOp.restrictDirichletEntries();
      grid.restrictDirichletNodes ( rhs );

      aol::BlockGaussSeidelPreconditioner<aol::MultiVector<RealType>, ElastOpType::BlockMatrixType> precond ( elastOp.getBlockMatrixRef() );
      aol::PCGInverse< aol::MultiVector<RealType> > pcgsolver ( elastOp.getBlockMatrixRef(), precond, eps, 100000 );
      pcgsolver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

      cerr << "Total memusage: " << ( aol::memusage() >> 20 ) << " MiB" << endl;

      pcgsolver.apply ( rhs, soln );

      soln += DirichletBCs;

      soln.save ( solnfilename, qc::PGM_DOUBLE_BINARY );

    } else if ( doWhat == 'C' ) {

      GridType GStdGrid ( GStdLevel );
      qc::ScalarArray<RealType, qc::QC_3D> GStdLevelset ( GStdGrid );

      NCCSetLevelset ( GStdLevelset );

      GStdGrid.addStructureFrom ( GStdLevelset );

      qc::AArray< GridType::NodalCoeffType, qc::QC_3D > GStdCoeff ( GStdGrid );

      tpcfe::setCoeffForLevelset ( GStdCoeff, GStdLevelset, ENuMinus, ENuPlus );
      GStdGrid.relaxedDetectAndInitVirtualNodes ( GStdCoeff, 1.0e-14, 2.0e-15 );

      qc::MultiArray< RealType, qc::QC_3D > GStdSoln ( GStdGrid ), compSoln ( grid );
      char filename[1024];
      sprintf ( filename, sfnmask, GStdLevel );
      GStdSoln.load ( filename );
      compSoln.load ( solnfilename );

      RealType L1sum = 0, L2sum = 0, Linfb = 0, Linfv = 0, H1sum = 0, volSum = 0, GL1S = 0, GL2S = 0, GLinf = 0, GH1S = 0;

      cerr << "Total memusage: " << aol::memusage() / 1048576.0 << " MB" << endl;

      aol::ProgressBar<> pb ( "Computing Difference" );
      pb.start ( GStdGrid.getNumberOfElements() );

      const RealType h = GStdGrid.H();

      const qc::GridSize<qc::QC_3D> gridSize ( GStdGrid );
      for ( GridType::FullElementIterator it ( GStdGrid ); it.notAtEnd(); ++it ) {
        tpcfe::CFEElement<RealType> el ( *it, gridSize, GStdGrid.getElType ( *it ) );

        el.computeAssembleData ( GStdGrid );

        const int BO = -1; // offset from boundary (grid cells on GStdGrid, if we loop over this) for computation of errors

        if ( ! ( ( el[0] < BO ) || ( el[1] < BO ) || ( el[2] < BO ) || ( el[0] > GStdGrid.getNumX() - 1 - BO ) || ( el[1] > GStdGrid.getNumY() - 1 - BO ) || ( el[2] > GStdGrid.getNumZ() - 1 - BO ) ) ) {

          for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
            const tpcfe::CFETetra<RealType> &tetra = *tit;
            const RealType vol = tetra.getVolume() * aol::Cub ( h );

            {
              aol::Vec3<RealType> barycenter;
              tetra.computeGlobalCoordinateOfBarycenter ( barycenter, el );
              barycenter *= h; // global -> world coordinates; barycenter of current tetrahedron     // use better quadrature??

              aol::Vec<4, RealType> bary;  bary.setAll( 0.25 );

              const aol::Vec3<RealType>
                GStdValue = tpcfe::interpolateDataAtPositionBary ( GStdGrid, GStdSoln, el, tetra, bary ),
                compValue = tpcfe::interpolateDataAtPositionWC ( grid, compSoln, barycenter );

              const RealType diffNorm = ( GStdValue - compValue ).norm(), GStdNorm = GStdValue.norm();

              L1sum += vol * diffNorm;
              L2sum += vol * aol::Sqr ( diffNorm );
              Linfb = aol::Max ( Linfb, diffNorm );

              GL1S += vol * fabs ( GStdNorm );
              GL2S += vol * aol::Sqr ( GStdNorm );
              GLinf = aol::Max ( GLinf, GStdNorm );

              volSum += vol;
            }

            {
              // more work for approximate gradients:
              aol::Vec3<RealType> vertices[4], GStdValues[4], compValues[4];

              for ( short v = 0; v < 4; ++v ) {
                tetra.computeGlobalCoordinate ( vertices[v], el, v );
                vertices[v] *= GStdGrid.H(); // convert to world coordinate
                aol::Vec<4, RealType> bary;  bary[v] = 1.0;
                GStdValues[v] = tpcfe::interpolateDataAtPositionBary ( GStdGrid, GStdSoln, el, tetra, bary );
                compValues[v] = tpcfe::interpolateDataAtPositionWC ( grid, compSoln, vertices[v], 1.0e-10 );

                Linfv = aol::Max ( Linfv, ( GStdValues[v] - compValues[v] ).norm() );
                GLinf = aol::Max ( GLinf, ( GStdValues[v] ).norm() );
              }

              aol::Matrix33<RealType> GStdDifference, compDifference, GStdApproxGrad, compApproxGrad, dirMat, invDirMat;

              for ( short v = 0; v < 3; ++v ) {
                dirMat.setRow ( v, vertices[v] - vertices[3] );
                GStdDifference[v] = GStdValues[v] - GStdValues[3];
                compDifference[v] = compValues[v] - compValues[3];
              }
              invDirMat = dirMat.inverse();

              GStdApproxGrad.makeProduct ( invDirMat,  GStdDifference ); // check order!
              compApproxGrad.makeProduct ( invDirMat,  compDifference );

              aol::Matrix33<RealType> diffMat = GStdApproxGrad; diffMat -= compApproxGrad;
              const RealType
                diffGradNorm = ( diffMat ).norm(),
                GStdApproxGradNorm = GStdApproxGrad.norm();

              H1sum += vol * aol::Sqr ( diffGradNorm );
              GH1S += vol * aol::Sqr ( GStdApproxGradNorm );
            }
          }
        }
      }

      pb.finish();

      cout << "Volume " << aol::detailedFormat ( volSum ) << endl;
      cout << "Linf difference: " << aol::detailedFormat ( aol::Max ( Linfb, Linfv ) ) << endl;
      cout << "L1   difference: " << aol::detailedFormat ( L1sum ) << endl;
      cout << "L2   difference: " << aol::detailedFormat ( sqrt ( L2sum ) ) << endl;
      cout << "H1   difference: " << aol::detailedFormat ( sqrt ( H1sum ) ) << endl;


      const RealType
        L1Norm = L1sum,
        L2Norm = sqrt ( L2sum ),
        H1Norm = sqrt ( L2sum + H1sum ),
        GL1N = GL1S,
        GL2N = sqrt ( GL2S ),
        GH1N = sqrt ( GL2S + GH1S );

      cerr << "Volume " << aol::detailedFormat ( volSum ) << endl
           << "Linf difference: " << aol::detailedFormat ( aol::Max ( Linfb, Linfv ) ) << endl
           << "L1 difference: " << aol::detailedFormat ( L1Norm ) << endl
           << "L2 difference: " << aol::detailedFormat ( L2Norm ) << endl
           << "H1 difference: " << aol::detailedFormat ( H1Norm ) << endl;

      cerr << "Relative differences:"
           << " " << aol::detailedFormat ( aol::Max ( Linfb, Linfv ) / GLinf )
           << " " << aol::detailedFormat ( L1Norm / GL1N )
           << " " << aol::detailedFormat ( L2Norm / GL2N )
           << " " << aol::detailedFormat ( H1Norm / GH1N ) << endl;
    }

  } catch ( aol::Exception &ex ) {
    ex.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
