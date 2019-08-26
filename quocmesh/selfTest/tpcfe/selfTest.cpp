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

// include all tpcfe header files
#include <FELevelsetOpInterface.h>
#include <integrateLevelSet.h>
#include <signedDistanceSweeping.h>
#include <tfeBaseFunctionSet.h>
#include <tpCFEBasics.h>
#include <tpCFEElastOp.h>
#include <tpCFEElement.h>
#include <tpCFEGrid.h>
#include <tpCFELevelsets.h>
#include <tpCFELookup.h>
#include <tpCFEMatrices.h>
#include <tpCFEMultigrid.h>
#include <tpCFEPeriodicBC.h>
#include <tpCFEStandardOp.h>
#include <tpCFETester.h>
#include <tpCFETetrahedron.h>
#include <tpCFEUtils.h>
#include <tpCFEVirtualNode.h>

// include additional quoc and aol headers necessary for generating test problems
#include <scalarArray.h>
#include <shapeLevelsetGenerator.h>
#include <generator.h>



void testFailure ( const bool failed, const char *errorText ) {
  if ( !failed ) {
    cerr << aol::color::green << "OK." << aol::color::reset << endl;
  } else {
    cerr << aol::color::red << errorText << endl << aol::color::reset;
  }
}


/**
 * \brief This class computes via "apply" the integral of the levelset passed to the constructor.
 */
template <typename ConfiguratorType>
class TestLevelsetIntegration : public qc::FENonlinLevelsetIntegrationVectorInterface<ConfiguratorType, ConfiguratorType::Dim, 1, TestLevelsetIntegration<ConfiguratorType> > {
  typedef typename ConfiguratorType::RealType RealType;

public:
  TestLevelsetIntegration ( const typename ConfiguratorType::InitType &Grid, const typename ConfiguratorType::ArrayType &LevelSetFunction ) :
      qc::FENonlinLevelsetIntegrationVectorInterface<ConfiguratorType, ConfiguratorType::Dim, 1, TestLevelsetIntegration<ConfiguratorType> > ( Grid, LevelSetFunction ) {}

  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType, 1> &/*DiscFuncs*/,
                               const typename ConfiguratorType::ElementType &/*El*/,
                               const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return static_cast<RealType> ( 1. );
  }
};


int main ( int, char** ) {

  try {

    bool failed = false;

    {
      cerr << "--- Testing whether indexing is correct ... ";
      aol::BaseFunctionSetMultiLin<double, qc::QC_3D, aol::GaussQuadrature< double, qc::QC_3D, 3 > > bfcts ( aol::NumberTrait<double>::NaN ); // value for h; should not be used!
      for ( unsigned short v = 0; v < 8; ++v ) {
        failed |= ( fabs ( bfcts.evaluate ( v, tpcfe::CFELookup<double>::_hexNbOffs[v] ) - 1.0 ) > 1e-10 );
      }
      testFailure ( failed, "<---------------- indexing test failed!\n" );
    }


    {
      cerr << "--- Testing whether volume of known object is computed correctly ... ";
      tpcfe::CFEGrid< double, tpcfe::CFE_CD> grid ( 4 );
      qc::ScalarArray< double, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<double>::generateColumnLevelset ( levelset, 0.3 );
      grid.setDomainFrom ( levelset );
      grid.detectAndInitVirtualNodes();
      failed |= ( fabs ( grid.getTotalInnerVolume() - 0.2783911111111 ) > 1e-10 );
      testFailure ( failed, "<---------------- volume test failed!\n" );
    }

    {
      cerr << "--- Testing whether constraints are computed correctly ... ";
      typedef tpcfe::CFEGrid< double, tpcfe::CFE_TPOS> GridType;
      GridType grid ( 4 );
      qc::ScalarArray< double, qc::QC_3D> levelset ( grid );
      qc::ShapeLevelsetGenerator<double>::generateBallLevelset ( levelset, 0.3 );
      grid.addStructureFrom ( levelset );

      qc::AArray<double, qc::QC_3D> coeff ( grid );
      for ( int i = 0; i < levelset.size(); ++i )
        coeff[i] = ( levelset[i] < 0 ? 1.0 : 2.0 );

      grid.detectAndInitVirtualNodes ( coeff );

      aol::Vector< GridType::VNIndexType > vnIndices ( grid.virtualNodeMapRef().size() );
      aol::Vector< int > constraints ( 8729 ); // reference length

      int i = 0, k = 0;
      for ( GridType::VNMapTypeConstIt it = grid.virtualNodeMapRef().begin(); it != grid.virtualNodeMapRef().end(); ++it, ++i ) {
        vnIndices[i] = it->first;
        for ( int j = 0; j < it->second->numConstraints(); ++j, ++k ) {
          constraints[k] = it->second->constraint ( j );
        }
      }

      failed |= ( aol::crc32 ( constraints.getData(), constraints.size() * sizeof ( int ) ) != 1496955181UL ); // reference value
      failed |= ( aol::crc32 ( vnIndices.getData(), vnIndices.size() * sizeof ( GridType::VNIndexType ) ) != 1181797016UL ); // reference value

      testFailure ( failed, "<---------------- constraint test failed!\n" );


      cerr << "--- Testing whether interface triangulation is determined correctly ... ";
      typedef     tpcfe::CFEInterfaceTriangulationGenerator<GridType> ITGType;

      ITGType itg( grid ); // not necessary to load any data
      itg.determineInterfaceTriangulation();
      ITGType::VertexVectorType vertices;
      ITGType::TriangVectorType triangles;
      itg.getVertexVector( vertices );
      itg.getTriangVector( triangles );

      double checkSum = 0.;
      for ( unsigned int i = 0; i < vertices.size(); ++i ) {
        for ( short j = 0; j < 3; ++j ) {
          checkSum += vertices[i][j] * sin ( static_cast<double> ( i ) ) * cos ( static_cast<double> ( j ) ); // this is probably also instable!
        }
      }

      failed |= ( fabs ( checkSum - ( - 2.82424662587386033e+00)  ) > 1e-6 );

      aol::Vector<int> triangInd ( 3 * triangles.size() );
      for ( unsigned int i = 0; i < triangles.size(); ++i ) {
        for ( short j = 0; j < 3; ++j ) {
          triangInd [ 3 * i + j ] = triangles[i][j];
        }
      }

      failed |= ( aol::crc32( triangInd.getData(), triangInd.size() * sizeof ( int ) ) != 1004437188UL );
      testFailure ( failed, "<---------------- interface triangulation test failed!\n" );
    }

    {
      cerr << "--- Testing whether CFEMassOp is assembled correctly ... ";

      typedef double RealType;
      static const tpcfe::ConstraintType AT = tpcfe::CFE_TPOS;
      typedef tpcfe::CFEGrid < RealType, AT >                                    GridType;
      typedef tpcfe::CFEConfigurator < GridType, aol::SparseMatrix<RealType> >   ConfiguratorType;
      typedef tpcfe::CFEMassOp  < ConfiguratorType >                             MassOpType;
      typedef aol::SparseMatrix<RealType>                                        MatrixType;

      const int depth = 1;

      const RealType DC_PLUS = 1.0, DC_MINUS = 2.0, threshold = 1. / 3.;

      GridType grid ( depth );

      qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid ); qc::ShapeLevelsetGenerator<RealType>::generateHalfcubeLevelset ( levelset, threshold );
      grid.addStructureFrom ( levelset );

      qc::AArray<RealType, qc::QC_3D> coeff ( grid );
      for ( int i = 0; i < levelset.size(); ++i )
        coeff[i] = ( levelset[i] < 0 ? DC_MINUS : DC_PLUS );

      grid.detectAndInitVirtualNodes ( coeff );



      qc::BitArray<qc::QC_3D> AF ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );                           AF.setAll ( false );
      grid.setDirichletMask ( AF );
      grid.setDOFMaskFromDirichletAndDomainNodeMask();

      MassOpType massOp ( grid );
      massOp._quietMode = true;

      {
        // check whether manually assembled mass matrix coincides with CFEMassOp
        MatrixType ManualMassMat ( grid );

        const qc::GridSize<qc::QC_3D> gridSize ( grid );
        for ( GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
          tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

          el.computeAssembleData ( grid );

          for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
            const tpcfe::CFETetra<RealType>& tetra = *tit;

            aol::Vec3<RealType> vertex[4]; // in world coordinates
            for ( int v = 0; v < 4; ++v ) {
              tetra.computeGlobalCoordinate ( vertex[v], el, v );
              vertex[v] *= grid.H();
            }

            aol::Matrix33<RealType> edge;
            for ( int e = 0; e < 3; ++e )
              edge[e] = vertex[e] - vertex[3];


            RealType volume = fabs ( edge.det() ) / 6; // this is the volume of the current tetrahedron: determinant of matrix of vectors spanning it divided by 6.

            const RealType referenceVolume = 1. / 6.;

            ConfiguratorType::MatrixType localMassMat ( 4, 4 );

            for ( int ilocal = 0; ilocal < 4; ++ilocal ) {
              for ( int jlocal = 0; jlocal < 4; ++jlocal ) {
                localMassMat.set ( ilocal, jlocal, volume / referenceVolume * ( ilocal == jlocal ? 1. / 60. : 1. / 120. ) ); // 1./60. / 1./120. are the local mass matrix entries for a reference tetrahedron spanned by (0,0,0), (0,0,1), (0,1,0), (1,0,0)
              }
            }

            for ( int ilocal = 0; ilocal < 4; ++ilocal ) {
              aol::Vector<int> iconstraints, jconstraints;
              aol::RandomAccessContainer<RealType> iweights, jweights;
              if ( tetra ( ilocal, 1 ) == tpcfe::NODE_NOT_VIRTUAL ) {

                iconstraints.reallocate ( 1 );
                iconstraints[0] = el.globalIndex ( tetra ( ilocal, 0 ) );
                iweights.reallocate ( 1 );
                iweights[0] = 1.0;

              } else {

                const int lIdx0i = tetra ( ilocal, 0 );
                const int lIdx1i = tetra ( ilocal, 1 );
                const int gIdx0i = el.globalIndex ( lIdx0i );
                const int gIdx1i = ( lIdx1i == tpcfe::NODE_NOT_VIRTUAL ? 0 : el.globalIndex ( lIdx1i ) );

                const GridType::VNType& VNi = grid.getVirtualNodeRef ( gIdx0i, gIdx1i );

                iconstraints.reallocate ( VNi._numOfConstraints );
                iconstraints = VNi._constraints;
                iweights.reallocate ( VNi._numOfConstraints );
                iweights = VNi._weights;
              }

              for ( int jlocal = 0; jlocal < 4; ++jlocal ) {
                if ( tetra ( jlocal, 1 ) == tpcfe::NODE_NOT_VIRTUAL ) {

                  jconstraints.reallocate ( 1 );
                  jconstraints[0] = el.globalIndex ( tetra ( jlocal, 0 ) );
                  jweights.reallocate ( 1 );
                  jweights[0] = 1.0;

                } else {

                  const int lIdx0j = tetra ( jlocal, 0 );
                  const int lIdx1j = tetra ( jlocal, 1 );
                  const int gIdx0j = el.globalIndex ( lIdx0j );
                  const int gIdx1j = ( lIdx1j == tpcfe::NODE_NOT_VIRTUAL ? 0 : el.globalIndex ( lIdx1j ) );

                  const GridType::VNType& VNj = grid.getVirtualNodeRef ( gIdx0j, gIdx1j );

                  jconstraints.reallocate ( VNj._numOfConstraints );
                  jconstraints = VNj._constraints;
                  jweights.reallocate ( VNj._numOfConstraints );
                  jweights = VNj._weights;
                }

                for ( int ii = 0; ii < iweights.size(); ++ii ) {
                  for ( int jj = 0; jj < jweights.size(); ++jj ) {

                    ManualMassMat.add ( iconstraints[ii], jconstraints[jj], iweights[ii] * jweights[jj] * localMassMat.get ( ilocal, jlocal ) );

                  }
                }

              }
            }

          }
        }

        failed |= !compareOps ( massOp, ManualMassMat, grid.getNumberOfNodes(), grid.getNumberOfNodes(), 1.0e-16 );

      }

      {
        // check whether manually assembled and MG-restricted mass matrix coincides with MassOp

        std::map< GridType::VNIndexType, int > vnodeIM;
        int counter = grid.getNumberOfNodes();
        for ( GridType::VNMapTypeConstIt vnit = grid.virtualNodeMapRef().begin(); vnit != grid.virtualNodeMapRef().end(); ++vnit ) {
          vnodeIM[ vnit->first ] = counter;
          ++counter;
        }

        MatrixType ManualBigMassMat ( counter, counter );

        const qc::GridSize<qc::QC_3D> gridSize ( grid );
        for ( GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
          tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

          el.computeAssembleData ( grid );

          for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
            const tpcfe::CFETetra<RealType>& tetra = ( *tit );

            aol::Vec3<RealType> vertex[4]; // in world coordinates
            for ( int v = 0; v < 4; ++v ) {
              tetra.computeGlobalCoordinate ( vertex[v], el, v );
              vertex[v] *= grid.H();
            }

            aol::Matrix33<RealType> edge;
            for ( int e = 0; e < 3; ++e )
              edge[e] = vertex[e] - vertex[3];


            RealType volume = fabs ( edge.det() ) / 6; // this is the volume of the current tetrahedron: determinant of matrix of vectors spanning it divided by 6.

            const RealType referenceVolume = 1. / 6.;

            MatrixType localMassMat ( 4, 4 );

            for ( int ilocal = 0; ilocal < 4; ++ilocal ) {
              for ( int jlocal = 0; jlocal < 4; ++jlocal ) {
                localMassMat.set ( ilocal, jlocal, volume / referenceVolume * ( ilocal == jlocal ? 1. / 60. : 1. / 120. ) ); // 1./60. / 1./120. are the local mass matrix entries for a reference tetrahedron spanned by (0,0,0), (0,0,1), (0,1,0), (1,0,0)
              }
            }

            for ( int ilocal = 0; ilocal < 4; ++ilocal ) {
              int iposition = -1;

              if ( tetra ( ilocal, 1 ) == tpcfe::NODE_NOT_VIRTUAL ) {
                iposition = el.globalIndex ( tetra ( ilocal, 0 ) );
              } else {
                const int lIdx0i = tetra ( ilocal, 0 );
                const int lIdx1i = tetra ( ilocal, 1 );
                const int gIdx0i = el.globalIndex ( lIdx0i );
                const int gIdx1i = ( lIdx1i == tpcfe::NODE_NOT_VIRTUAL ? 0 : el.globalIndex ( lIdx1i ) );

                iposition = vnodeIM[ GridType::VNType::mapIndex ( gIdx0i, gIdx1i ) ];
              }

              for ( int jlocal = 0; jlocal < 4; ++jlocal ) {
                int jposition = -1;
                if ( tetra ( jlocal, 1 ) == tpcfe::NODE_NOT_VIRTUAL ) {
                  jposition = el.globalIndex ( tetra ( jlocal, 0 ) );
                } else {

                  const int lIdx0j = tetra ( jlocal, 0 );
                  const int lIdx1j = tetra ( jlocal, 1 );
                  const int gIdx0j = el.globalIndex ( lIdx0j );
                  const int gIdx1j = ( lIdx1j == tpcfe::NODE_NOT_VIRTUAL ? 0 : el.globalIndex ( lIdx1j ) );

                  jposition = vnodeIM[ GridType::VNType::mapIndex ( gIdx0j, gIdx1j ) ];
                }

                ManualBigMassMat.add ( iposition, jposition, localMassMat.get ( ilocal, jlocal ) );


              }
            }

          }
        }

        // now set up restriction/prolongation between virtual and regular grid and coarsen virtual mass mat to regular mass mat.
        MatrixType restrMat ( grid.getNumberOfNodes(), counter ), prolongMat ( counter, grid.getNumberOfNodes() );
        for ( int i = 0; i < grid.getNumberOfNodes(); ++i ) {
          restrMat.set ( i, i, 1.0 );
          prolongMat.set ( i, i, 1.0 );
        }

        for ( GridType::VNMapTypeConstIt vnit = grid.virtualNodeMapRef().begin(); vnit != grid.virtualNodeMapRef().end(); ++vnit ) {
          for ( int i = 0; i < vnit->second->_numOfConstraints; ++i ) {
            prolongMat.set ( vnodeIM[ vnit->first ], vnit->second->constraint ( i ), vnit->second->weight ( i ) );
            restrMat.set ( vnit->second->constraint ( i ), vnodeIM[ vnit->first ], vnit->second->weight ( i ) );
          }
        }

        MatrixType tmpMat ( counter, grid.getNumberOfNodes() );
        tmpMat.addMatrixProduct ( ManualBigMassMat, prolongMat );
        MatrixType ManualMassMat ( grid );
        ManualMassMat.addMatrixProduct ( restrMat, tmpMat );

        failed |= !compareOps ( massOp, ManualMassMat, grid.getNumberOfNodes(), grid.getNumberOfNodes(), 1.0e-16 );
      }

      testFailure ( failed, "<---------------- CFEMassOp assembly test failed!\n" );

      cerr << "--- Testing whether CFEStiffOpWI coincides with reference matrix ... ";
      tpcfe::CFEStiffOpWI  < ConfiguratorType > stiffOp ( coeff, grid );
      stiffOp._quietMode = true;
      MatrixType stiffMat ( grid );
      stiffOp.assembleAddMatrix ( stiffMat );
      RealType checkSum = 0;
      for ( int i = 0; i < stiffMat.getNumRows(); ++i ) {
        for ( int j = 0; j < stiffMat.getNumCols(); ++j ) {
          checkSum += aol::Sqr ( stiffMat.get( i, j ) * sin ( static_cast<double> ( i ) ) * cos ( static_cast<double> ( j ) ) );
        }
      }
      failed |= ( fabs ( checkSum - 1.36676154786314843e+01 ) > 1e-10 );
      testFailure ( failed, "<---------------- CFEStiffOpWI test failed!\n" );
    }

    {
      cerr << "--- Testing whether CFEElastOp and CFEFullyAnisotropicElastOp with isotropic tensor coincide ... ";
      typedef tpcfe::CFEGrid<double, tpcfe::CFE_CD> GridType;
      GridType cfeGrid ( 3 );
      qc::ScalarArray<double, qc::QC_3D> levelset ( cfeGrid );
      levelset.setAll ( -1.0 );
      cfeGrid.setDomainFrom ( levelset );

      cfeGrid.detectAndInitVirtualNodes();

      qc::BitArray<qc::QC_3D> allfalse ( qc::GridSize<qc::QC_3D>::createFrom ( cfeGrid ) );
      allfalse.setAll ( false );
      cfeGrid.setDirichletMask ( allfalse );
      cfeGrid.setDOFMaskFromDirichletAndDomainNodeMask();

      const double E = 42, nu = 0.23;
      const double lambda = tpcfe::computeLambda ( E, nu ), mu = tpcfe::computeMu ( E, nu );

      tpcfe::CFEElastOp< tpcfe::CFEConfigurator< GridType, tpcfe::CFEBandMatrix< GridType > > >  ElastOp ( cfeGrid, lambda, mu );

      aol::Vec3<double> tensor[3][3][3];
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          for ( int k = 0; k < 3; ++k )
            for ( int l = 0; l < 3; ++l )
              tensor[i][j][k][l] = lambda * ( i == j ) * ( k == l ) + mu * ( ( i == k ) * ( j == l ) + ( i == l ) * ( j == k ) );

      tpcfe::CFEFullyAnisotropicElastOp< tpcfe::CFEConfigurator< GridType, tpcfe::CFEHybridMatrix< GridType > > >  AElastOp ( cfeGrid, tensor );

      bool okay = true;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
          okay &= aol::compareOps ( ElastOp.getBlockMatrixRef().getReference ( 0, 0 ), AElastOp.getBlockMatrixRef().getReference ( 0, 0 ), cfeGrid.getNumberOfNodes(), cfeGrid.getNumberOfNodes(), 1.E-14 );

      failed |= !okay;

      testFailure ( failed, "<---------------- CFEElastOpTest failed!\n" );
    }


    {
      cerr << "--- Testing whether CFEElastOp (no complicated domain) and CFEJCEElastOp (interface with trivial jump) coincide ... ";

      aol::SparseMatrix<double> diffMassMat, diffElastMat;
      const int level = 2;

      {
        const tpcfe::ConstraintType CT = tpcfe::CFE_TPOSELAST;
        typedef tpcfe::CFEGrid< double, CT, tpcfe::IsotropicElasticityCoefficient<double> > GridType;
        typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
        typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

        GridType grid ( level );

        qc::ScalarArray<double, qc::QC_3D> levelset ( grid );
        for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
          levelset.set ( *bit, grid.H() * ( *bit ).norm() - 0.85 * sqrt ( 3.0 ) );
        }

        const tpcfe::IsotropicElasticityCoefficient<double> ENuMinus ( 1.0, 0.1 ), ENuPlus ( 1.0, 0.1 );

        qc::AArray< tpcfe::IsotropicElasticityCoefficient<double>, qc::QC_3D > coeff ( grid );
        for ( int i = 0; i < coeff.size(); ++i ) {
          coeff[i] = ( levelset[i] < 0 ? ENuMinus : ENuPlus );
        }

        grid.addStructureFrom ( levelset );
        grid.relaxedDetectAndInitVirtualNodes ( coeff, 5.0e-16, 1.0e-16 );

        tpcfe::CFEJCEMassOp< ConfiguratorType > massOp ( grid, true /*quiet*/ );
        tpcfe::CFEJCElastOp< ConfiguratorType > elastOp ( grid, coeff, true /*quiet*/ );

        double checkSum = 0;
        for ( int i = 0; i < massOp.getBlockMatrixRef().getReference(1,1).getNumRows(); ++i ) {
          for ( int j = 0; j < massOp.getBlockMatrixRef().getReference(1,1).getNumCols(); ++j ) {
            checkSum += aol::Sqr ( massOp.getBlockMatrixRef().getReference(1,1).get( i, j ) * sin ( static_cast<double> ( i ) ) * cos ( static_cast<double> ( j ) ) );
          }
        }
        failed |= ( fabs ( checkSum - 2.99435986406196709e-04 ) > 1e-13 );
        checkSum = 0;
        for ( int i = 0; i < elastOp.getBlockMatrixRef().getReference(0,1).getNumRows(); ++i ) {
          for ( int j = 0; j < elastOp.getBlockMatrixRef().getReference(0,1).getNumCols(); ++j ) {
            checkSum += aol::Sqr ( elastOp.getBlockMatrixRef().getReference(0,1).get( i, j ) * sin ( static_cast<double> ( i ) ) * cos ( static_cast<double> ( j ) ) );
          }
        }
        failed |= ( fabs ( checkSum - 6.80178709729035935e-01 ) > 1e-10 );

        diffMassMat.resize ( grid.getNumberOfNodes(), grid.getNumberOfNodes() );
        diffMassMat += massOp.getBlockMatrixRef().getReference ( 0, 0 );
        diffElastMat.resize ( 3 * grid.getNumberOfNodes(), 3 * grid.getNumberOfNodes() );
        elastOp.getBlockMatrixRef().addUnblockedMatrixTo ( diffElastMat );
      }

      diffMassMat *= -1.0;
      diffElastMat *= -1.0;

      {
        const tpcfe::ConstraintType CT = tpcfe::CFE_CD;
        typedef tpcfe::CFEGrid< double, CT > GridType;
        typedef aol::SparseMatrix<double> MatrixType;
        typedef tpcfe::CFEConfigurator < GridType, MatrixType  > ConfiguratorType;

        GridType grid ( level );

        qc::ScalarArray<double, qc::QC_3D> levelset ( grid );
        levelset.setAll ( -1 );

        grid.setDomainFrom ( levelset );
        grid.detectAndInitVirtualNodes();

        tpcfe::CFEMassOp< ConfiguratorType > massOp ( grid );
        massOp._quietMode = true;
        massOp.assembleAddMatrix ( diffMassMat );

        const double E = 1.0, nu = 0.1;
        tpcfe::CFEElastOp< ConfiguratorType > elastOp ( grid, tpcfe::computeLambda ( E, nu ), tpcfe::computeMu ( E, nu ) );

        elastOp.getBlockMatrixRef().addUnblockedMatrixTo ( diffElastMat );
      }

      failed |= ( diffMassMat.getFrobeniusNormSqr() > 1.0e-28 ) || ( diffElastMat.getFrobeniusNormSqr() > 1.0e-28 );

      testFailure ( failed, "<---------------- CFEJCEElastOpTest failed!\n" );
    }

    {
      cerr << "--- Testing integrateLevelSet in 3d on a plane ... ";
      qc::GridDefinition grid ( 4, qc::QC_3D );
      const double h = grid.H();
      qc::ScalarArray<double, qc::QC_3D> function ( grid );
      qc::ScalarArray<double, qc::QC_3D> levelSet ( grid );

      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          for ( int k = 0; k < grid.getWidth(); ++k ) {
            function.set ( i, j, k, k*h );
            levelSet.set ( i, j, k, k*h - 0.27 );
          }
        }
      }

      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> testPlaneIntegrator ( grid, levelSet, function );
      failed = failed || ( fabs ( testPlaneIntegrator.integrate() - 0.27 ) > 1.e-6 );

      testFailure ( failed, "<---------------- IntegrateFunctionOverLevelSet failed!\n" );

      // test the integrateSaveEachValue-method
      cerr << "--- Testing integrateSaveEachValue in 3d on a plane... ";
      qc::ScalarArray<double, qc::QC_3D> integrals ( grid );
      testPlaneIntegrator.integrateSaveEachValue ( integrals );
      failed = failed || ( fabs ( integrals.sum() - 0.27 ) > 1.e-6 );
      testFailure ( failed, "<---------------- integrateSaveEachValue on a plane" );

      // test the integrateBall-method
      cerr << "--- Testing integrateBall in 3d ... ";
      aol::Vec3<double> coords ( 0.5, 0.5, 0.27 );
      function.setAll ( 1. );
      failed = failed || ( fabs ( testPlaneIntegrator.integrateBall ( 0.3, coords ) - 0.34375 ) > 1.e-4 );
      testFailure ( failed, "<---------------- IntegrateBall failed!\n" );

      // generate a sphere as leveset-surface
      for ( int i = 0; i < grid.getWidth(); ++i ) {
        for ( int j = 0; j < grid.getWidth(); ++j ) {
          for ( int k = 0; k < grid.getWidth(); ++k ) {
            const double x = ( i * h ) - 0.5;
            const double y = ( j * h ) - 0.5;
            const double z = ( k * h ) - 0.5;
            const double value = sqrt ( x * x + y * y + z * z );
            levelSet.set ( i, j, k, value - 0.25 );
            function.set ( i, j, k, 1. );
          }
        }
      }

      cerr << "--- Testing integrate in 3d on a sphere... ";
      qc::IntegrateFunctionOverLevelSet<double, qc::QC_3D> testBallIntegrator ( grid, levelSet, function );
      failed = failed || ( fabs ( testBallIntegrator.integrate() - 0.77318 ) > 1.e-4 );
      testFailure ( failed, "<---------------- integrate on a sphere" );

      cerr << "--- Testing integrateSaveEachValue in 3d on a sphere... ";
      integrals.setZero();
      testBallIntegrator.integrateSaveEachValue ( integrals );
      failed = failed || ( fabs ( integrals.sum() - 0.77318 ) > 1.e-4 );

      testFailure ( failed, "<---------------- integrateSaveEachValue in 3d on a sphere failed" );
    }


    {
      cerr << "--- Testing FEOps operating on levelsets ...";

      { // test in 2D
        // define a circle or sphere as levelset
        typedef double RealType;
        typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
        ConfiguratorType::InitType grid ( 4, ConfiguratorType::Dim );
        ConfiguratorType::ArrayType levelset ( grid );
        ConfiguratorType::VecType center;
        center.setAll ( .5 );
        RealType radius = .4;
        qc::DataGenerator<ConfiguratorType> ( grid ).generateSphereLevelset ( levelset, center, radius );
        aol::MultiVector<RealType> levelsetMultiVector;
        levelsetMultiVector.appendReference ( levelset );
        aol::MultiVector<RealType> aux ( levelsetMultiVector, aol::STRUCT_COPY );
        RealType tol = .01;

        { // test FENonlinLevelsetIntegrationVectorInterface
          TestLevelsetIntegration<ConfiguratorType> integrator ( grid, levelset );
          aol::Scalar<RealType> integral;
          integrator.apply ( levelsetMultiVector, integral );
          RealType diff = aol::Abs ( integral - 2 * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetIntegrationVectorInterface accuracy on level 4: " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FENonlinLevelsetVectorDiffOpInterface
          qc::LevelsetVectorStiffOp<ConfiguratorType, 1, 1> stiffOp ( grid, levelset );
          stiffOp.apply ( levelsetMultiVector, aux );
          RealType diff = aol::Abs ( aux * levelsetMultiVector - 2 * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetVectorDiffOpInterface accuracy on level 4:      " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FENonlinLevelsetVectorOpInterface
          qc::LevelsetVectorMassOp<ConfiguratorType, 1> massOp ( grid, levelset );
          aol::MultiVector<RealType> ones ( aux, aol::STRUCT_COPY );
          ones.setAll ( 1. );
          massOp.apply ( ones, aux );
          RealType diff = aol::Abs ( aux * ones - 2 * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetVectorOpInterface accuracy on level 4:          " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FELinLevelsetOpInterface
          qc::LevelsetMassOp<ConfiguratorType> massOp ( grid, levelset );
          aol::Vector<RealType> ones ( levelset, aol::STRUCT_COPY );
          ones.setAll ( 1. );
          massOp.apply ( ones, aux[0] );
          RealType diff = aol::Abs ( aux[0] * ones - 2 * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FELinLevelsetOpInterface accuracy on level 4:                   " << diff << endl;
          failed = failed || ( diff > tol );
        }
      }

      { // test in 3D (exactly the same code as in 2D, only typedef ConfiguratorType and sphere surface changes)
        // define a circle or sphere as levelset
        typedef double RealType;
        typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_3D, aol::GaussQuadrature<RealType, qc::QC_3D, 3> > ConfiguratorType;
        ConfiguratorType::InitType grid ( 4, ConfiguratorType::Dim );
        ConfiguratorType::ArrayType levelset ( grid );
        ConfiguratorType::VecType center;
        center.setAll ( .5 );
        RealType radius = .4;
        qc::DataGenerator<ConfiguratorType> ( grid ).generateSphereLevelset ( levelset, center, radius );
        aol::MultiVector<RealType> levelsetMultiVector;
        levelsetMultiVector.appendReference ( levelset );
        aol::MultiVector<RealType> aux ( levelsetMultiVector, aol::STRUCT_COPY );
        RealType tol = .04;

        { // test FENonlinLevelsetIntegrationVectorInterface
          TestLevelsetIntegration<ConfiguratorType> integrator ( grid, levelset );
          aol::Scalar<RealType> integral;
          integrator.apply ( levelsetMultiVector, integral );
          RealType diff = aol::Abs ( integral - 4 * radius * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetIntegrationVectorInterface accuracy on level 4: " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FENonlinLevelsetVectorDiffOpInterface
          qc::LevelsetVectorStiffOp<ConfiguratorType, 1, 1> stiffOp ( grid, levelset );
          stiffOp.apply ( levelsetMultiVector, aux );
          RealType diff = aol::Abs ( aux * levelsetMultiVector - 4 * radius * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetVectorDiffOpInterface accuracy on level 4:      " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FENonlinLevelsetVectorOpInterface
          qc::LevelsetVectorMassOp<ConfiguratorType, 1> massOp ( grid, levelset );
          aol::MultiVector<RealType> ones ( aux, aol::STRUCT_COPY );
          ones.setAll ( 1. );
          massOp.apply ( ones, aux );
          RealType diff = aol::Abs ( aux * ones - 4 * radius * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FENonlinLevelsetVectorOpInterface accuracy on level 4:          " << diff << endl;
          failed = failed || ( diff > tol );
        }

        { // test FELinLevelsetOpInterface
          qc::LevelsetMassOp<ConfiguratorType> massOp ( grid, levelset );
          aol::Vector<RealType> ones ( levelset, aol::STRUCT_COPY );
          ones.setAll ( 1. );
          massOp.apply ( ones, aux[0] );
          RealType diff = aol::Abs ( aux[0] * ones - 4 * radius * radius * aol::NumberTrait<RealType>::pi );
          // cerr << "    " << ConfiguratorType::Dim << "D FELinLevelsetOpInterface accuracy on level 4:                   " << diff << endl;
          failed = failed || ( diff > tol );
        }
      }

      testFailure ( failed, "<---------------- FEOp on levelsets failed" );
    }

    {
      cerr << "--- Testing SweepingSignedDistanceOp3D  ...";
      typedef double RealType;

      { // planar interface
        typedef tpcfe::CFEGrid < RealType, tpcfe::CFE_NONE > GridType;
        GridType grid ( 4 );
        qc::ScalarArray< RealType, qc::QC_3D > soln ( grid ), SDF ( grid ), difference ( grid );
        qc::ShapeLevelsetGenerator<RealType>::generateRotatedHalfcubeLevelset ( soln, 1./3., 42 );
        qc::SignedDistanceSweepingOp3D<RealType> SDsweeper;
        SDsweeper.apply ( soln, SDF );

        // soln *= grid.H();
        difference = SDF;   difference -= soln;
        cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;
        failed &= ( difference.norm() / difference.size() < 1.0e-4 ) && ( difference.getMaxAbsValue() < 0.05 );
      }

      { // spherical interface
        qc::ScalarArray< RealType, qc::QC_3D > soln ( 14, 17, 21 );
        qc::ScalarArray< RealType, qc::QC_3D > SDF ( soln, aol::STRUCT_COPY ), difference ( soln, aol::STRUCT_COPY );

        const RealType h = 1.0 / aol::Max ( soln.getNumX() - 1, soln.getNumY() - 1 , soln.getNumZ() - 1 );
        for ( qc::RectangularIterator<qc::QC_3D> rit ( soln ); rit.notAtEnd(); ++rit ) {
          const int i = ( *rit ) [0], j = ( *rit ) [1], k = ( *rit ) [2];
          const RealType x = h * i, y = h * j, z = h * k;
          const RealType value = sqrt ( aol::Sqr ( x - 0.5 ) + aol::Sqr ( y - 0.5 ) + aol::Sqr ( z - 0.5 ) ) - 1./3.;
          soln.set ( *rit, value );
        }

        qc::SignedDistanceSweepingOp3D<RealType> SDsweeper;
        SDsweeper.apply ( soln, SDF );

        difference = SDF;   difference -= soln;
        cerr << "approximate L2 norm of difference = " << difference.norm() / difference.size() << ", Linf norm = " << difference.getMaxAbsValue() << endl;
        failed &= ( difference.norm() / difference.size() < 2.0e-4 ) && ( difference.getMaxAbsValue() < 0.1 );
      }

      testFailure ( failed, "<---------------- SweepingSignedDistanceOp3D failed" );

    }

    if ( !failed ) {
      aol::printSelfTestSuccessMessage ( "--                      TPCFE SelfTest Successful                             --" );
      aol::callSystemPauseIfNecessaryOnPlatform();
      return ( EXIT_SUCCESS );
    }

  } catch ( aol::Exception &ex ) {
    cerr << aol::color::error << "\n\naol::Exception caught:\n" << aol::color::reset;
    ex.dump ();
  }

  aol::printSelfTestFailureMessage ( "!!                      TPCFE SelfTest FAILED                                 !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return ( EXIT_FAILURE );
}
