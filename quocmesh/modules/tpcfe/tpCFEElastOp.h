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

#ifndef __TPCFEELASTOP_H
#define __TPCFEELASTOP_H

#include <tpCFEUtils.h>
#include <tpCFEStandardOp.h>
#include <quoc.h>

namespace tpcfe {

/** Basis class for CFEMixedDerivativeOp and CFEMixedDerivativeOpWI. Instances of this class are probably not useful.
 * \author Schwen
 */
template < typename ConfiguratorType >
class CFEMixedDerivativeOpBase {
protected:
  typedef typename ConfiguratorType::TetraMatrixType TetraMatrixType;
  typedef typename ConfiguratorType::RealType RealType;

  int _i, _j;
  const typename ConfiguratorType::GridType &_grid;

  /** This constructor is protected so that it cannot create instances */
  CFEMixedDerivativeOpBase ( const int i, const int j,
                             const typename ConfiguratorType::GridType & Grid ) : _i ( i ), _j ( j ), _grid ( Grid ) {}

  void createMatrixForTetra ( TetraMatrixType &matrix, const CFETetra<RealType> & t ) const {
    // Prepare barycentric coordinates
    t.computeInverseTransformation();

    // Compute scaling factor
    const RealType det = t.determinant();
    const RealType scale = this->_grid.H() / ( 6 * det ); // note that the computation via barycentric coordinates has a factor det^2 which needs to be divided out here (see CFEStiffOp).
    RealType DiL0 = 0, DjL0 = 0;

    for ( int a = 0; a < 3; ++a ) {
      DiL0 -= t._barycentricCoordinates[a][_i];
      DjL0 -= t._barycentricCoordinates[a][_j];
      for ( int b = 0; b < 3; ++b ) {
        matrix[a+1][b+1] = t._barycentricCoordinates[a][_i] * t._barycentricCoordinates[b][_j] * scale;
      }
    }

    matrix[0][0] = DiL0 * DjL0 * scale;
    for ( int a = 0; a < 3; ++a ) {
      matrix[a+1][0] = t._barycentricCoordinates[a][_i] * DjL0 * scale;
      matrix[0][a+1] = DiL0 * t._barycentricCoordinates[a][_j] * scale;
    }
  }

  // end class CFEMixedDerivativeOpBase
};


/** Composite Finite Element Matrix with mixed derivatives, i. e.
 *  \f$ \int_\Omega \partial_I\phi \,\partial_J\psi \f$
 */
template < typename ConfiguratorType >
class CFEMixedDerivativeOp : public CFEStandardOp < ConfiguratorType, CFEMixedDerivativeOp < ConfiguratorType > > , public CFEMixedDerivativeOpBase< ConfiguratorType > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

public:
  CFEMixedDerivativeOp ( const int i, const int j,
                         const typename ConfiguratorType::GridType & Grid,
                         const aol::OperatorType OpType = aol::ONTHEFLY ) :
      CFEStandardOp < ConfiguratorType,  CFEMixedDerivativeOp< ConfiguratorType > > ( Grid, OpType ), CFEMixedDerivativeOpBase < ConfiguratorType > ( i, j, Grid ) {
    this->init();
  }

  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    std::vector< tpcfe::CFETetra<RealType> > tv;
    tpcfe::CFELookup<RealType>::writeRefTetraVecTo ( tv );

    // loop over the tetrahedra of a standard (non interfaced) element
    for ( typename std::vector< tpcfe::CFETetra<RealType> >::const_iterator it = tv.begin(); it != tv.end(); ++it ) {
      this->createMatrixForTetra ( this->_defaultLocalTetraMatrix[ it->getParent() ], *it );
    }
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType >& /*e*/, const CFETetra<RealType> & t ) const {
    this->createMatrixForTetra ( this->_localTetraMatrix, t );
  }

  // end class CFEMixedDerivativeOp
};


template < typename ConfiguratorType >
class CFEAssembledBlockOp : public aol::Op < aol::MultiVector< typename ConfiguratorType::RealType > > {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::GridType   GridType;
  typedef typename ConfiguratorType::MatrixType MatrixType;
  typedef aol::BlockMatrix< MatrixType >        BlockMatrixType;

protected:
  const GridType &_grid;
  aol::BlockMatrix< MatrixType > _blockMatrix;

protected:
  CFEAssembledBlockOp ( const GridType &Grid ) : _grid ( Grid ), _blockMatrix ( Grid ) {
  }

public:
  void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _blockMatrix.apply ( Arg, Dest );
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _blockMatrix.applyAdd ( Arg, Dest );
  }

  BlockMatrixType& getBlockMatrixRef() {
    return ( _blockMatrix );
  }

  const BlockMatrixType& getBlockMatrixRef() const {
    return ( _blockMatrix );
  }

  void assembleAddBlockMatrix ( aol::BlockMatrix< MatrixType > &dest ) {
    for ( short i = 0; i < 3; ++i ) {
      for ( short j = 0; j < 3; ++j ) {
        dest._blockMatrix.getReference ( i, j ) += this->_blockMatrix.getReference ( i, j ) ;
      }
    }
  }

  void restrictNonDomainEntries() {
    for ( short i = 0; i < 3; ++i ) {
      for ( short j = 0; j < 3; ++j ) {
        tpcfe::restrictNonDomainEntries ( _grid, this->_blockMatrix.getReference ( i, j ), static_cast<RealType> ( i == j ) );
      }
    }
  }

  void restrictDirichletEntries() {
    for ( short i = 0; i < 3; ++i ) {
      for ( short j = 0; j < 3; ++j ) {
        tpcfe::restrictDirichletEntries ( _grid, this->_blockMatrix.getReference ( i, j ), static_cast<RealType> ( i == j ) );
      }
    }
  }

  void restrictNonDOFEntries() {
    for ( short i = 0; i < 3; ++i ) {
      for ( short j = 0; j < 3; ++j ) {
        tpcfe::restrictNonDOFEntries ( _grid, this->_blockMatrix.getReference ( i, j ), static_cast<RealType> ( i == j ) );
      }
    }
  }

private:
  CFEAssembledBlockOp ( );
  CFEAssembledBlockOp ( const CFEAssembledBlockOp<ConfiguratorType>& );
  CFEAssembledBlockOp<ConfiguratorType>& operator= ( const CFEAssembledBlockOp<ConfiguratorType>& );
  // end class CFEAssembledBlockOp
};


/** Composite Finite Element Elasticity Operator for complicated domains based on CFEMixedDerivativeOp, only for constant isotropic elasticity given by Lame numbers
 */
template < typename _ConfiguratorType >
class CFEElastOp : public CFEAssembledBlockOp < _ConfiguratorType > {
public:
  typedef _ConfiguratorType                                      ConfiguratorType;
  typedef aol::MultiVector<typename ConfiguratorType::RealType>  VectorType;

  typedef typename CFEAssembledBlockOp<ConfiguratorType>::GridType    GridType;
  typedef typename CFEAssembledBlockOp<ConfiguratorType>::RealType    RealType;
  typedef typename CFEAssembledBlockOp<ConfiguratorType>::MatrixType  MatrixType;

protected:
  const RealType _lambda, _mu;

public:
  CFEElastOp ( const GridType &Grid, const RealType lambda, const RealType mu ) : CFEAssembledBlockOp < ConfiguratorType > ( Grid ), _lambda ( lambda ), _mu ( mu ) {
    if ( GridType::CT != CFE_CD )
      throw aol::Exception ( "CFEElastOp may only be used with CFE_CD", __FILE__, __LINE__ );

    for ( int a = 0; a < 3; ++a ) {
      for ( int b = 0; b < 3; ++b ) {
        MatrixType mdMatrix ( this->_grid );
        {
          tpcfe::CFEMixedDerivativeOp< ConfiguratorType > MDOp ( a, b, this->_grid );
#ifndef VERBOSE
          MDOp._quietMode = true;
#endif
          MDOp.assembleAddMatrix ( mdMatrix );
        }

        this->_blockMatrix.getReference ( a, b ).addMultiple ( mdMatrix, _lambda );
        this->_blockMatrix.getReference ( b, a ).addMultiple ( mdMatrix,     _mu );

        if ( a == b ) {
          for ( int c = 0; c < 3; ++c ) {
            this->_blockMatrix.getReference ( c, c ).addMultiple ( mdMatrix, _mu );
          }
        }
      }
    }
  }

  // end class CFEElastOp
};


/** Composite Finite Element Elasticity Operator for fully anisotropic material (given by fourth-order tensor)
 *  Only for complicated domain!
 */
template < typename _ConfiguratorType >
class CFEFullyAnisotropicElastOp : public CFEAssembledBlockOp < _ConfiguratorType > {
public:
  typedef _ConfiguratorType                                      ConfiguratorType;
  typedef aol::MultiVector<typename ConfiguratorType::RealType>  VectorType;

  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::GridType   GridType;
  typedef typename ConfiguratorType::MatrixType MatrixType;

public:
  //! constructor for constant coefficients throughout interior
  CFEFullyAnisotropicElastOp ( const GridType &Grid, const aol::Vec3<RealType> Tensor[3][3][3] ) : CFEAssembledBlockOp<ConfiguratorType> ( Grid ) {

    if ( this->_grid.CT != CFE_CD ) {
      throw aol::Exception ( "tpcfe::CFEFullyAnisotropicElastOp is currently only for CFE_CD.", __FILE__, __LINE__ );
    }

    for ( int r = 0; r < 3; ++r ) {
      for ( int s = 0; s < 3; ++s ) {
        MatrixType md_matrix ( this->_grid );
        {
          tpcfe::CFEMixedDerivativeOp< ConfiguratorType > md_op ( r, s, this->_grid );
#ifndef VERBOSE
          md_op._quietMode = true;
#endif
          md_op.assembleAddMatrix ( md_matrix );
        }

        for ( int a = 0; a < 3; ++a )
          for ( int b = 0; b < 3; ++b )
            this->_blockMatrix.getReference ( a, b ).addMultiple ( md_matrix, ( Tensor[r][a][b][s] + Tensor[r][a][s][b] + Tensor[a][r][b][s] + Tensor[a][r][s][b] ) / 4 );
      }
    }

  }

  //! constructor for varying coefficients throughout interior (make sure to use appropriate matrix structure!)
  //! copy&paste&modify from CFEJCElastOp, think about common base class or something different ...
  CFEFullyAnisotropicElastOp ( const GridType &Grid, const qc::AArray<typename GridType::NodalCoeffType, qc::QC_3D> &Coeff ) : CFEAssembledBlockOp<ConfiguratorType> ( Grid ) {

    if ( this->_grid.CT != CFE_CD ) {
      throw aol::Exception ( "tpcfe::CFEFullyAnisotropicElastOp is currently only for CFE_CD.", __FILE__, __LINE__ );
    }

    // buildDefaultTetraMatrices
    aol::Mat<4,4,RealType> defaultLocalTetraMatrix[6][3][3]; // standard tetra with edge length h; first partial derivative; second partial derivative

    {
      std::vector< tpcfe::CFETetra<RealType> > tv;
      tpcfe::CFELookup<RealType>::writeRefTetraVecTo ( tv );

      for ( typename std::vector< tpcfe::CFETetra<RealType> >::const_iterator it = tv.begin(); it != tv.end(); ++it ) {
        computeLocalTetraMatrix ( *it, defaultLocalTetraMatrix[ it->getParent() ] );
      }
    }

    // objects that we do not want to construct and destroy in the inner loops ...
    aol::Vector<int> iConstrNonvirt ( 1 ), jConstrNonvirt ( 1 );
    aol::RandomAccessContainer< typename GridType::VNType::ExtrapolationWeightType > weiNonvirt ( 1 );
    weiNonvirt[ 0 ] = aol::ZOTrait< typename GridType::VNType::ExtrapolationWeightType >::one;


    aol::ProgressBar<> pb ( "Assembling anisotropic ElastOp" );
    pb.start ( this->_grid.getNumberOfElements() );

    const qc::GridSize<qc::QC_3D> gridSize ( this->_grid );
    for ( typename GridType::FullElementIterator elit ( this->_grid ); elit.notAtEnd(); ++elit, pb++ ) {
      tpcfe::CFEElement<RealType> el ( *elit, gridSize, this->_grid.getElType ( *elit ) );

      el.computeAssembleData ( this->_grid );
      CFEWeightProvider< RealType, typename GridType::NodalCoeffType > weightProvider ( Coeff, el );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); tit.notAtEnd(); ++tit ) {
        const CFETetra<RealType> &tet = *tit;

        RealType localTensor[3][3][3][3];
        weightProvider.meanWeight ( tet.getSign() ).getAnisotropicTensor ( localTensor );

        aol::Mat<4,4,RealType> localTetraMatrix[3][3];

        if ( tet.isStandard() ) {

          const RealType volumeFactor = tet.getVolume() / tpcfe::CFELookup<RealType>::_stdTetraVolume;
          for ( short a = 0; a < 3; ++a )
            for ( short b = 0; b < 3; ++b )
              for ( short i = 0; i < 4; ++i )
                for ( short j = 0; j < 4; ++j )
                  localTetraMatrix[a][b][i][j] = defaultLocalTetraMatrix[ tet.getParent() ][a][b][i][j] * volumeFactor; // tet is its own parent if standard

        } else {

          CFETetra<RealType> tetCopy ( tet ); // need to modify
          el.getTetraVertices ( tetCopy );

          computeLocalTetraMatrix ( tetCopy, localTetraMatrix );

        }

        // think about optimization of these loops (at most two weights! scalar weights!)
        for ( short i = 0; i < 4; ++i ) {
          const bool iVirt = tet.isVirtualNode ( i );
          iConstrNonvirt[0] = el.globalIndex ( tet ( i, 0 ) );
          const aol::Vector<int> &iConstr = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._constraints : iConstrNonvirt );
          const aol::RandomAccessContainer< typename GridType::VNType::ExtrapolationWeightType > &iWei = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._weights : weiNonvirt );

          for ( short j = 0; j < 4; ++j ) {
            const bool jVirt = tet.isVirtualNode ( j );
            jConstrNonvirt[0] = el.globalIndex ( tet ( j, 0 ) );
            const aol::Vector<int> &jConstr = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._constraints : jConstrNonvirt );
            const aol::RandomAccessContainer< typename GridType::VNType::ExtrapolationWeightType > &jWei = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._weights : weiNonvirt );

            for ( short i1 = 0; i1 < iConstr.size(); ++i1 ) {
              for ( short j1 = 0; j1 < jConstr.size(); ++j1 ) {
                aol::Matrix33<RealType>
                  wI ( iWei[i1], 0, 0, 0, iWei[i1], 0, 0, 0, iWei[i1] ),
                  wJ ( jWei[j1], 0, 0, 0, jWei[j1], 0, 0, 0, jWei[j1] );
                const int cI = iConstr[i1], cJ = jConstr[j1];

                for ( short A = 0; A < 3; ++A ) {
                  for ( short B = 0; B < 3; ++B ) {
                    RealType dwarf = 0;
                    for ( short k = 0; k < 3; ++k )
                      for ( short l = 0; l < 3; ++l )
                        for ( short m = 0; m < 3; ++m )
                          for ( short n = 0; n < 3; ++n )
                            dwarf += localTensor[k][l][m][n] / 4 * ( wI[n][A] * wJ[l][B] * localTetraMatrix[m][k][i][j] +
                                                                     wI[n][A] * wJ[k][B] * localTetraMatrix[m][l][i][j] +
                                                                     wI[m][A] * wJ[l][B] * localTetraMatrix[n][k][i][j] +
                                                                     wI[m][A] * wJ[k][B] * localTetraMatrix[n][l][i][j]   );

                    this->_blockMatrix.getReference ( A, B ).add ( cI, cJ, dwarf );
                  }
                }
              }
            }
          }
        }

      }
    }

    pb.finish();

  }

protected:
  //! copy&paste from CFEJCElastOp, think about common base class or something different ...
  void computeLocalTetraMatrix ( const CFETetra<RealType> &tet, aol::Mat<4,4,RealType> localTetraMat[3][3] /*automatically call-by-reference*/ ) const {
    // due to the non-obvious construction of the "barycentricCoordinates" of CFETetrahedron, we do not use the method used in the CFEMixedDerivativeOp

    const RealType scale = this->_grid.H() * tet.getVolume(); // h^2 cancels with the two factors h due to barycentric coordinates wrt unit cube

    aol::Matrix33<RealType> dirMat, dirMatInv;
    for ( short i = 0; i < 3; ++i )
      for ( short j = 0; j < 3; ++j )
        dirMat.set ( i, j, tet._edge[i][j] );
    dirMatInv = dirMat.inverse();

    aol::Vec3<RealType> partDerivBFct[4];
    for ( short a = 0; a < 3; ++a ) {
      for ( short i = 0; i < 3; ++i ) {
        partDerivBFct[ 0 ][a] += dirMatInv[a][i];
        partDerivBFct[i+1][a] -= dirMatInv[a][i];
      }
    }

    for ( short a = 0; a < 3; ++a ) {
      for ( short b = 0; b < 3; ++b ) {
        for ( short i = 0; i < 4; ++i ) {
          for ( short j = 0; j < 4; ++j ) {
            localTetraMat[a][b][i][j] = scale * partDerivBFct[i][a] * partDerivBFct[j][b];
          }
        }
      }
    }
  }


private:
  CFEFullyAnisotropicElastOp ( );
  CFEFullyAnisotropicElastOp ( const CFEFullyAnisotropicElastOp<_ConfiguratorType>& );
  CFEFullyAnisotropicElastOp<_ConfiguratorType>& operator= ( const CFEFullyAnisotropicElastOp<_ConfiguratorType>& );

  // end class CFEFullyAnisotropicElastOp
};

/** MassOp for CFE_TPOSELAST
 *  Note that it has a different sparsity structure and assembling it is currently incompatible to the usual CFE Operators.
 *  \todo adapt to CFE_CDWI_ELAST case
 */
template < typename _ConfiguratorType >
class CFEJCEMassOp : public CFEAssembledBlockOp< _ConfiguratorType > {
public:
  typedef _ConfiguratorType                                       ConfiguratorType;
  typedef aol::MultiVector<typename ConfiguratorType::RealType>   VectorType;

protected:
  typedef typename ConfiguratorType::RealType                     RealType;
  typedef typename ConfiguratorType::GridType                     GridType;
  typedef typename GridType::VNType::ExtrapolationWeightType      ExtrapolationWeightType;

  aol::Mat< 8, 8, RealType > _defaultLocalHexaMatrix;
  aol::Mat< 4, 4, RealType > _defaultLocalTetraMatrix; // these are the same for all six standard tetra

public:
  CFEJCEMassOp ( const GridType &Grid, const bool QuietMode = false ) : CFEAssembledBlockOp< _ConfiguratorType > ( Grid ) {
    if ( GridType::CT != CFE_TPOSELAST ) {
      throw ( aol::Exception ( "Illegal constraint type.", __FILE__, __LINE__ ) );
      // not implemented yet for CFE_TPOSELAST, but should work similar to the CFEJCElastOp
    }

    buildDefaultTetraMatrix();
    buildHexaMatrix();

    // objects that we do not want to construct and destroy in the inner loops ...
    aol::Vector<int> iConstrNonvirt ( 1 ), jConstrNonvirt ( 1 );
    aol::RandomAccessContainer< ExtrapolationWeightType > iWeiNonvirt ( 1 ), jWeiNonvirt ( 1 ); // one is sufficient here!
    iWeiNonvirt[ 0 ] = aol::ZOTrait< ExtrapolationWeightType >::one;
    jWeiNonvirt[ 0 ] = aol::ZOTrait< ExtrapolationWeightType >::one;

    aol::ProgressBar<> pb ( "Assembling JCEMassOp" );
    pb.start ( this->_grid.getNumberOfElements() );

    const qc::GridSize<qc::QC_3D> gridSize ( this->_grid );
    for ( typename GridType::FullElementIterator elit ( this->_grid ); elit.notAtEnd(); ++elit ) {
      if ( !QuietMode )
        pb++;

      tpcfe::CFEElement<RealType> el ( *elit, gridSize, this->_grid.getElType ( *elit ) );

      if ( ! ( el.cfeType().representsInterfaced() ) ) {

        // add up contributions to global matrix for non-interfaced hexahedron
        for ( short i = 0; i < 8; ++i ) {
          const int i1 = el.globalIndex ( i );
          for ( short j = 0; j < 8; ++j ) {
            const int j1 = el.globalIndex ( j );
            for ( short A = 0; A < 3; ++A ) {
              this->_blockMatrix.getReference ( A, A ).add ( i1, j1, _defaultLocalHexaMatrix[i][j] ); // same for all three
            }
          }
        }

      } else {

        el.computeAssembleData ( this->_grid );
        for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {
          const CFETetra<RealType> &tet = *tit;

          aol::Mat<4,4,RealType> localTetraMatrix;
          const RealType volumeFactor = tet.getVolume() / tpcfe::CFELookup<RealType>::_stdTetraVolume;
          for ( short i = 0; i < 4; ++i )
            for ( short j = 0; j < 4; ++j )
              localTetraMatrix[i][j] = _defaultLocalTetraMatrix[i][j] * volumeFactor;

          for ( short i = 0; i < 4; ++i ) {
            const bool iVirt = tet.isVirtualNode ( i );
            iConstrNonvirt[0] = el.globalIndex ( tet ( i, 0 ) );
            const aol::Vector<int> &iConstr = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._constraints : iConstrNonvirt );
            const aol::RandomAccessContainer< ExtrapolationWeightType > &iWei = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._weights : iWeiNonvirt );

            for ( short j = 0; j < 4; ++j ) {
              const bool jVirt = tet.isVirtualNode ( j );
              jConstrNonvirt[0] = el.globalIndex ( tet ( j, 0 ) );
              const aol::Vector<int> &jConstr = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._constraints : jConstrNonvirt );
              const aol::RandomAccessContainer< ExtrapolationWeightType > &jWei = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._weights : jWeiNonvirt );

              for ( short i1 = 0; i1 < iConstr.size(); ++i1 ) {
                for ( short j1 = 0; j1 < jConstr.size(); ++j1 ) {
                  aol::Matrix33<RealType> coeff;
                  coeff.makeProductAtransposedB ( iWei[i1], jWei[j1] );
                  for ( short A = 0; A < 3; ++A ) {
                    for ( short B = 0; B < 3; ++B ) {
                      this->_blockMatrix.getReference ( A, B ).add ( iConstr[i1], jConstr[j1], coeff.get(A,B) * localTetraMatrix[i][j] );
                    }
                  }
                }
              }
            }
          }
        }

      }
    }

    if ( !QuietMode )
      pb.finish();

  }

protected:
  void buildDefaultTetraMatrix ( ) {
    const RealType stdTetraVolFactor = aol::Cub ( this->_grid.H() );
    for ( short i = 0; i < 4; ++i ) {
      for ( short j = 0; j < 4; ++j ) {
        _defaultLocalTetraMatrix[i][j] = tpcfe::CFELookup<RealType>::_localMassMatrixRefTetra[i][j] * stdTetraVolFactor;
      }
    }
  }

  void buildHexaMatrix ( ) {
    for ( CFETopoTetraIterator tit ( CFEType ( -1, 0 ) ); tit.notAtEnd(); ++tit ) {
      for ( short i = 0; i < 4; ++i ) {
        const short i1 = ( *tit ) ( i, 0 );
        for ( short j = 0; j < 4; ++j ) {
          const short j1 = ( *tit ) ( j, 0 );
          _defaultLocalHexaMatrix[i1][j1] += _defaultLocalTetraMatrix[i][j];
        }
      }
    }

  }

private:
  CFEJCEMassOp ( );
  CFEJCEMassOp ( const CFEJCEMassOp<_ConfiguratorType>& );
  CFEJCEMassOp<_ConfiguratorType>& operator= ( const CFEJCEMassOp<_ConfiguratorType>& );

  // end class CFEJCEMassOp
};


/** ElastOp for CFE_TPOSELAST
 *  Note that it has a different sparsity structure and assembling it is currently incompatible to the usual CFE Operators.
 *  Note also that currently the WeightInterface contains local values of E and nu (isotropic linear elasticity). If used with a CFEHybridMatrix, these must be constant throughout the plus/minus domains.
 *  \todo adapt to CFE_CDWI_ELAST case
 */
template < typename _ConfiguratorType >
class CFEJCElastOp : public CFEAssembledBlockOp< _ConfiguratorType > {
public:
  typedef _ConfiguratorType                                       ConfiguratorType;
  typedef aol::MultiVector<typename ConfiguratorType::RealType>   VectorType;

  typedef typename ConfiguratorType::RealType                     RealType;
  typedef typename ConfiguratorType::GridType                     GridType;
  typedef typename GridType::VNType::ExtrapolationWeightType      ExtrapolationWeightType;
  typedef typename GridType::NodalCoeffType                       NodalCoeffType;

protected:
  const qc::AArray<NodalCoeffType, qc::QC_3D> &_nodalCoeff;
  aol::Mat<4,4,RealType> _defaultLocalTetraMatrix[6][3][3]; // last indices are first partial derivative, second partial derivative

public:
  CFEJCElastOp ( const GridType &Grid, const qc::AArray<NodalCoeffType, qc::QC_3D> &NodalCoeff, const bool QuietMode = false )
    : CFEAssembledBlockOp< _ConfiguratorType > ( Grid ),
      _nodalCoeff ( NodalCoeff ) {

    if ( ( GridType::CT != CFE_TPOSELAST ) && ( GridType::CT != CFE_CDWI_ELAST ) ) {
      throw ( aol::Exception ( "Illegal constraint type.", __FILE__, __LINE__ ) );
    }

    buildDefaultTetraMatrices();

    // objects that we do not want to construct and destroy in the inner loops ...
    aol::Vector<int> iConstrNonvirt ( 1 ), jConstrNonvirt ( 1 );
    aol::RandomAccessContainer< ExtrapolationWeightType > weiNonvirt ( 1 );
    weiNonvirt[ 0 ] = aol::ZOTrait< typename GridType::VNType::ExtrapolationWeightType >::one;

    aol::ProgressBar<> pb ( "Assembling JCElastOp" );
    pb.start ( this->_grid.getNumberOfElements() );

    const qc::GridSize<qc::QC_3D> gridSize ( this->_grid );
    for ( typename GridType::FullElementIterator elit ( this->_grid ); elit.notAtEnd(); ++elit ) {
      if ( !QuietMode )
        pb++;

      tpcfe::CFEElement<RealType> el ( *elit, gridSize, this->_grid.getElType ( *elit ) );

      el.computeAssembleData ( this->_grid );

      signed char titSign = -2;

      if ( this->_grid.hasDomain() ) {
        const int structureNo =  el.cfeType()._structureNo;
        if ( el.cfeType().representsInterfaced() ) {
          if ( structureNo == tpcfe::MAX_STRUCT_ID ) {
            // at domain boundary, only loop over inner tet
            titSign = -1;
          } else {
            // at interface (between coefficient domains) boundary
            titSign = 0;
          }
        } else {
          // not interfaced: pureType refers to domain, so all inner will use whichever coefficient domain
          titSign = -1;
        }
      } else {
        titSign = 0; // loop over all tet if there is no domain
      }

#if 0
      // if not interfaced
      //   if not outside
      //     "add up contributions"
      //   else do nothing
      // else
      //   if interior
      //     loop over all subtetra
      //   else
      //     loop over inner subtetra only
      //   endif
      // endif

      // const int structureNo =  el.cfeType()._structureNo;
      // if(el.pureCFEType() == 0){ continue; } // bit set ^= minus sign ^= interior, i.e. 0 is all outside
      // if(structureNo == tpcfe::MAX_STRUCT_ID && (tetra.getSign() != -1)){ continue; }
#endif
      CFEWeightProvider< RealType, NodalCoeffType > weightProvider ( NodalCoeff, el );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, titSign ); tit.notAtEnd(); ++tit ) {
        const CFETetra<RealType> &tet = *tit;

        RealType localTensor[3][3][3][3];
        weightProvider.meanWeight ( tet.getSign() ).getAnisotropicTensor ( localTensor );

        aol::Mat<4,4,RealType> localTetraMatrix[3][3];

        if ( tet.isStandard() ) {

          const RealType volumeFactor = tet.getVolume() / tpcfe::CFELookup<RealType>::_stdTetraVolume;
          for ( short a = 0; a < 3; ++a )
            for ( short b = 0; b < 3; ++b )
              for ( short i = 0; i < 4; ++i )
                for ( short j = 0; j < 4; ++j )
                  localTetraMatrix[a][b][i][j] = _defaultLocalTetraMatrix[ tet.getParent() ][a][b][i][j] * volumeFactor; // tet is its own parent if standard

        } else {

          computeLocalTetraMatrix ( tet, localTetraMatrix );

        }

        // this loop structure is similar to the CFEJCEMassOp, but the overall assembling routines are sufficiently different to justify similar-duplicate code
        for ( short i = 0; i < 4; ++i ) {
          const bool iVirt = tet.isVirtualNode ( i );
          iConstrNonvirt[0] = el.globalIndex ( tet ( i, 0 ) );
          const aol::Vector<int> &iConstr = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._constraints : iConstrNonvirt );
          const aol::RandomAccessContainer< ExtrapolationWeightType > &iWei = ( iVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( i, 0 ) ), el.globalIndex ( tet ( i, 1 ) ) )._weights : weiNonvirt );

          for ( short j = 0; j < 4; ++j ) {
            const bool jVirt = tet.isVirtualNode ( j );
            jConstrNonvirt[0] = el.globalIndex ( tet ( j, 0 ) );
            const aol::Vector<int> &jConstr = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._constraints : jConstrNonvirt );
            const aol::RandomAccessContainer< ExtrapolationWeightType > &jWei = ( jVirt ? this->_grid.getVirtualNodeRef ( el.globalIndex ( tet ( j, 0 ) ), el.globalIndex ( tet ( j, 1 ) ) )._weights : weiNonvirt );

            for ( short i1 = 0; i1 < iConstr.size(); ++i1 ) {
              for ( short j1 = 0; j1 < jConstr.size(); ++j1 ) {
                const aol::Matrix33<RealType> &wI = iWei[i1], &wJ = jWei[j1];
                const int cI = iConstr[i1], cJ = jConstr[j1];

                for ( short A = 0; A < 3; ++A ) {
                  for ( short B = 0; B < 3; ++B ) {
                    RealType dwarf = 0;
                    for ( short k = 0; k < 3; ++k )
                      for ( short l = 0; l < 3; ++l )
                        for ( short m = 0; m < 3; ++m )
                          for ( short n = 0; n < 3; ++n )
                            dwarf += localTensor[k][l][m][n] / 4 * ( wI[n][A] * wJ[l][B] * localTetraMatrix[m][k][i][j] +
                                                                     wI[n][A] * wJ[k][B] * localTetraMatrix[m][l][i][j] +
                                                                     wI[m][A] * wJ[l][B] * localTetraMatrix[n][k][i][j] +
                                                                     wI[m][A] * wJ[k][B] * localTetraMatrix[n][l][i][j]   );

                    this->getBlockMatrixRef().getReference ( A, B ).add ( cI, cJ, dwarf );
                  }
                }
              }
            }
          }
        }

      }
    }

    if ( !QuietMode )
      pb.finish();

  }

private:
  CFEJCElastOp ( );
  CFEJCElastOp ( const CFEJCElastOp<_ConfiguratorType>& );
  CFEJCElastOp<_ConfiguratorType>& operator= ( const CFEJCElastOp<_ConfiguratorType>& );

  void buildDefaultTetraMatrices ( ) {
    std::vector< tpcfe::CFETetra<RealType> > tv;
    tpcfe::CFELookup<RealType>::writeRefTetraVecTo ( tv );

    // loop over the tetrahedra of a standard (non interfaced) element
    for ( typename std::vector< tpcfe::CFETetra<RealType> >::const_iterator it = tv.begin(); it != tv.end(); ++it ) {
      computeLocalTetraMatrix ( *it, _defaultLocalTetraMatrix[ it->getParent() ] );
    }
  }

  void computeLocalTetraMatrix ( const CFETetra<RealType> &tet, aol::Mat<4,4,RealType> localTetraMat[3][3] /*automatically call-by-reference*/ ) const {
    // due to the non-obvious construction of the "barycentricCoordinates" of CFETetrahedron, we do not use the method used in the CFEMixedDerivativeOp

    const RealType scale = this->_grid.H() * tet.getVolume(); // h^2 cancels with the two factors h due to barycentric coordinates wrt unit cube

    aol::Matrix33<RealType> dirMat, dirMatInv;
    for ( short i = 0; i < 3; ++i )
      for ( short j = 0; j < 3; ++j )
        dirMat.set ( i, j, tet._edge[i][j] );
    dirMatInv = dirMat.inverse();

    aol::Vec3<RealType> partDerivBFct[4];
    for ( short a = 0; a < 3; ++a ) {
      for ( short i = 0; i < 3; ++i ) {
        partDerivBFct[ 0 ][a] += dirMatInv[a][i];
        partDerivBFct[i+1][a] -= dirMatInv[a][i];
      }
    }

    for ( short a = 0; a < 3; ++a ) {
      for ( short b = 0; b < 3; ++b ) {
        for ( short i = 0; i < 4; ++i ) {
          for ( short j = 0; j < 4; ++j ) {
            localTetraMat[a][b][i][j] = scale * partDerivBFct[i][a] * partDerivBFct[j][b];
          }
        }
      }
    }
  }

  // end class CFEJCElastOp
};

} // end namespace tpcfe

#endif
