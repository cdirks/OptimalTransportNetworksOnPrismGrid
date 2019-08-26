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

#include <preconditioner.h>
#include <eigenvectors.h>
#include <shapeLevelsetGenerator.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMatrices.h>
#include <tpCFEStandardOp.h>
#include <tpCFEUtils.h>


template< typename RealType >
void computeConditionNumbers ( const aol::Matrix<RealType> &matrix, const aol::BitVector &useEntry ) {

  aol::SparseMatrix<RealType> mat ( useEntry.numOccurence(true), useEntry.numOccurence(true) );
  aol::Vector<RealType> indMap ( useEntry.numOccurence(true) );
  for ( int i = 0, oi= 0 ; i < useEntry.size(); ++i ) {
    if ( useEntry[i] == true ) {
      indMap[oi] = i;
      ++oi;
    }
  }


  for ( int i = 0; i < indMap.size(); ++i ) {
    for ( int j = 0; j < indMap.size(); ++j ) {
      mat.set ( i, j, matrix.get( indMap[i], indMap[j] ) );
    }
  }

  cerr << "condition number " << aol::longScientificFormat ( aol::computeConditionNumberByVectorIteration ( mat, 1.0e-14 ) ) << endl;

  {
    aol::DiagonalPreconditioner< aol::Vector<RealType> > DiagPrec ( mat );
    aol::DiagonalMatrix<RealType> DiagPrecMat ( mat.getNumRows() );
    DiagPrec.getPreconditionerMatrix ( DiagPrecMat );
    for ( int i = 0; i < DiagPrecMat.getNumRows(); ++i )
      DiagPrecMat.set ( i, i, sqrt ( DiagPrecMat.get(i,i) ) );
    aol::SparseMatrix<RealType> DiagPrecondMat ( mat.getNumRows(), mat.getNumCols() ), dummyMat ( mat.getNumRows(), mat.getNumCols() );
    dummyMat.addMatrixProduct ( DiagPrecMat, mat );
    DiagPrecondMat.addMatrixProduct ( dummyMat, DiagPrecMat );

    cerr << "condition number of preconditioned matrix (diagonal) " << aol::longScientificFormat ( aol::computeConditionNumberByVectorIteration ( DiagPrecondMat, 1.0e-14 ) ) << endl;
  }

  {
    aol::GeometricScalingPreconditioner< aol::Vector<RealType> > GScaPrec ( mat );
    aol::DiagonalMatrix<RealType> GScaPrecMat ( mat.getNumRows() );
    GScaPrec.getPreconditionerMatrix ( GScaPrecMat );
    for ( int i = 0; i < GScaPrecMat.getNumRows(); ++i )
      GScaPrecMat.set ( i, i, sqrt ( GScaPrecMat.get(i,i) ) );
    aol::SparseMatrix<RealType> GScaPrecondMat ( mat.getNumRows(), mat.getNumCols() ), dummyMat ( mat.getNumRows(), mat.getNumCols() );
    dummyMat.addMatrixProduct ( GScaPrecMat, mat );
    GScaPrecondMat.addMatrixProduct ( dummyMat, GScaPrecMat );

    cerr << "condition number of preconditioned matrix (geometric scaling) " << aol::longScientificFormat ( aol::computeConditionNumberByVectorIteration ( GScaPrecondMat, 1.0e-14) ) << endl;
  }

}


typedef double RealType;

int main ( int, char** ) {
  try {

    const int START_LEVEL = 2;
    const int STOP_LEVEL = 3;

    for ( int lev = START_LEVEL; lev < STOP_LEVEL; ++lev ) {
      cerr << "Level: " << lev << endl;

      const int INTERFACESTEPS = 8;

      for ( int updown = 0; updown < 2; ++updown ) {
        for ( int exponent = lev + 1; exponent < lev + 1 + INTERFACESTEPS; ++exponent ) {
          const RealType position = ( updown == 0 ? 0.5 + 1.0 / ( 1 << exponent ) : 0.5 + 1.0 / ( 1 << lev ) - 1.0 / ( 1 << exponent ) );
          cerr << "Position: " << position << endl;

          {
            const tpcfe::ConstraintType CT = tpcfe::CFE_CD;

            typedef tpcfe::CFEGrid< RealType, CT >                                               GridType;
            typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEBandMatrix<GridType> >  ConfiguratorType;
            typedef ConfiguratorType::MatrixType                                               MatrixType;

            typedef tpcfe::CFEMassOp  < ConfiguratorType >                                     MassOpType;
            typedef tpcfe::CFEStiffOp  < ConfiguratorType >                                   StiffOpType;

            GridType grid ( lev );
            qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
            qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, position );

            grid.setDomainFrom ( levelset );
            grid.detectAndInitVirtualNodes();
            grid.setDOFMaskFromDirichletAndDomainNodeMask();

            {
              MatrixType sysMat ( grid );
              {
                StiffOpType StiffOp ( grid, aol::ONTHEFLY );
                StiffOp._quietMode = true;
                MassOpType MassOp ( grid, aol::ONTHEFLY );
                MassOp._quietMode = true;
                StiffOp.assembleAddMatrix ( sysMat );
                sysMat *= grid.H();
                MassOp.assembleAddMatrix ( sysMat );
              }
              restrictNonDomainEntries ( grid,  sysMat, 1.0 );

              cerr << "condition numbers of M + h L matrix (Complicated Domain): " << endl;
              computeConditionNumbers ( sysMat, grid.getDOFMaskRef() );

            }

            {
              qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
              DirichletMask.setAll ( false );

              for ( int i = 0; i < grid.getWidth(); ++i ) {
                for ( int j = 0; j < grid.getWidth(); ++j ) {
                  DirichletMask.set ( i, j, 0, true );
                  DirichletMask.set ( i, j, grid.getWidth() - 1, true );
                }
              }

              grid.setDirichletMask ( DirichletMask );
              grid.setDOFMaskFromDirichletAndDomainNodeMask();

              typedef tpcfe::CFEElastOp< ConfiguratorType > ElastOpType;
              ElastOpType ElastOp ( grid, tpcfe::computeLambda ( 1.0, 0.33 ), tpcfe::computeMu ( 1.0, 0.33 ) );

              ElastOp.restrictNonDomainEntries();
              ElastOp.restrictDirichletEntries();

              aol::SparseMatrix<RealType> ElastMatUnblocked ( 0, 0 );
              ElastOp.getBlockMatrixRef().addUnblockedMatrixTo( ElastMatUnblocked );

              aol::BitVector DOFMaskUnblocked ( 3 * grid.getNumberOfNodes() );
              for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
                DOFMaskUnblocked.set ( i, grid.isDOFNode(i) );
                DOFMaskUnblocked.set ( grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
                DOFMaskUnblocked.set ( 2 * grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
              }

              cerr << "condition numbers of elasticity matrix (Complicated Domain):" << endl;
              computeConditionNumbers ( ElastMatUnblocked, DOFMaskUnblocked );

            }
          } // end CFE_CD

          {
            const tpcfe::ConstraintType CT = tpcfe::CFE_TPOS;

            typedef tpcfe::CFEGrid< RealType, CT >                                                        GridType;
            typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >  ConfiguratorType;
            typedef ConfiguratorType::MatrixType                                                          MatrixType;

            typedef tpcfe::CFEMassOp  < ConfiguratorType >                                                MassOpType;
            typedef tpcfe::CFEStiffOp  < ConfiguratorType >                                               StiffOpType;

            GridType grid ( lev );
            qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
            qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, position );

            grid.addStructureFrom ( levelset );
            qc::AArray< RealType, qc::QC_3D > coeff ( grid );
            tpcfe::setCoeffForLevelset ( coeff, levelset, 42.0, 1.0 );
            grid.relaxedDetectAndInitVirtualNodes( coeff, 1.0e-14, 1.0e-14 );
            grid.setDOFMaskFromDirichletAndDomainNodeMask();

            MatrixType sysMat ( grid );
            {
              StiffOpType StiffOp ( grid, aol::ONTHEFLY );
              StiffOp._quietMode = true;
              MassOpType MassOp ( grid, aol::ONTHEFLY );
              MassOp._quietMode = true;
              StiffOp.assembleAddMatrix ( sysMat );
              sysMat *= grid.H();
              MassOp.assembleAddMatrix ( sysMat );
            }


            cerr << "condition number of M + h L matrix (Jump Coeff)" << endl;
            computeConditionNumbers ( sysMat, grid.getDOFMaskRef() );

          } // end CFE_TPOS

          {
            const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;

            typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> >  GridType;
            typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >             ConfiguratorType;
            typedef ConfiguratorType::MatrixType                                                      MatrixType;

            GridType grid ( lev );
            qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
            qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, position );

            grid.addStructureFrom ( levelset );
            qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
            const GridType::NodalCoeffType ENuMinus ( 70.0, 0.35 ), ENuPlus ( 3.0, 0.38 );
            tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
            grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

            qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
            DirichletMask.setAll ( false );

            for ( int i = 0; i < grid.getWidth(); ++i ) {
              for ( int j = 0; j < grid.getWidth(); ++j ) {
                DirichletMask.set ( i, j, 0, true );
                DirichletMask.set ( i, j, grid.getWidth() - 1, true );
              }
            }

            grid.setDirichletMask ( DirichletMask );
            grid.setDOFMaskFromDirichletAndDomainNodeMask();

            typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;
            ElastOpType ElastOp ( grid, coeff );

            aol::SparseMatrix<RealType> ElastMatUnblocked ( 0, 0 );
            ElastOp.getBlockMatrixRef().addUnblockedMatrixTo( ElastMatUnblocked );

            aol::BitVector DOFMaskUnblocked ( 3 * grid.getNumberOfNodes() );
            for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
              DOFMaskUnblocked.set ( i, grid.isDOFNode(i) );
              DOFMaskUnblocked.set ( grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
              DOFMaskUnblocked.set ( 2 * grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
            }

            cerr << "condition numbers of elasticity matrix (Jump Coeff):" << endl;
            computeConditionNumbers ( ElastMatUnblocked, DOFMaskUnblocked );

          } // end CFE_TPOSELAST

        }
      }

      const int KAPPASTEPS = 25;
      RealType kappa = 1.0e-06;
      for ( int i = 0; i < KAPPASTEPS; ++i, kappa *= sqrt ( 10.0 ) ) {
        const RealType position = 0.5 + 1.0 / ( 1 << ( lev + 1) );
        cerr << "Kappa = " << kappa << endl;

        {
          const tpcfe::ConstraintType CT = tpcfe::CFE_TPOS;

          typedef tpcfe::CFEGrid< RealType, CT >                                                 GridType;
          typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >  ConfiguratorType;
          typedef ConfiguratorType::MatrixType                                                 MatrixType;

          typedef tpcfe::CFEMassOp  < ConfiguratorType >                                       MassOpType;
          typedef tpcfe::CFEStiffOp  < ConfiguratorType >                                     StiffOpType;

          GridType grid ( lev );
          qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
          qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, position );

          grid.addStructureFrom ( levelset );
          qc::AArray< RealType, qc::QC_3D > coeff ( grid );
          tpcfe::setCoeffForLevelset ( coeff, levelset, kappa, 1.0 );
          grid.relaxedDetectAndInitVirtualNodes( coeff, 1.0e-14, 1.0e-14 );
          grid.setDOFMaskFromDirichletAndDomainNodeMask();

          MatrixType sysMat ( grid );
          {
            StiffOpType StiffOp ( grid, aol::ONTHEFLY );
            StiffOp._quietMode = true;
            MassOpType MassOp ( grid, aol::ONTHEFLY );
            MassOp._quietMode = true;
            StiffOp.assembleAddMatrix ( sysMat );
            sysMat *= grid.H();
            MassOp.assembleAddMatrix ( sysMat );
          }

          cerr << "condition number of M + h L matrix (Jump Coeff)" << endl;
          computeConditionNumbers ( sysMat, grid.getDOFMaskRef() );
        } // end CFE_TPOS

        {
          const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;

          typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> >  GridType;
          typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> >             ConfiguratorType;
          typedef ConfiguratorType::MatrixType                                                      MatrixType;

          GridType grid ( lev );
          qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );
          qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, position );

          grid.addStructureFrom ( levelset );
          qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );
          const GridType::NodalCoeffType ENuMinus ( kappa, 0.1 ), ENuPlus ( 1, 0.3 );
          tpcfe::setCoeffForLevelset ( coeff, levelset, ENuMinus, ENuPlus );
          grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-15, 1.0e-15 );

          qc::BitArray<qc::QC_3D> DirichletMask ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
          DirichletMask.setAll ( false );

          for ( int i = 0; i < grid.getWidth(); ++i ) {
            for ( int j = 0; j < grid.getWidth(); ++j ) {
              DirichletMask.set ( i, j, 0, true );
              DirichletMask.set ( i, j, grid.getWidth() - 1, true );
            }
          }

          grid.setDirichletMask ( DirichletMask );
          grid.setDOFMaskFromDirichletAndDomainNodeMask();

          typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;
          ElastOpType ElastOp ( grid, coeff );

          aol::SparseMatrix<RealType> ElastMatUnblocked ( 0, 0 );
          ElastOp.getBlockMatrixRef().addUnblockedMatrixTo( ElastMatUnblocked );

          aol::BitVector DOFMaskUnblocked ( 3 * grid.getNumberOfNodes() );
          for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
            DOFMaskUnblocked.set ( i, grid.isDOFNode(i) );
            DOFMaskUnblocked.set ( grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
            DOFMaskUnblocked.set ( 2 * grid.getNumberOfNodes() + i, grid.isDOFNode(i) );
          }

          cerr << "condition numbers of elasticity matrix (Jump Coeff):" << endl;
          computeConditionNumbers ( ElastMatUnblocked, DOFMaskUnblocked );

        } // end CFE_TPOSELAST

      }
    }
  } catch ( aol::Exception &exc ) {
    exc.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
