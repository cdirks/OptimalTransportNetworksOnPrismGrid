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

#ifndef __QUOCMULTIGRID_H
#define __QUOCMULTIGRID_H

#include <restriction.h>
#include <prolongation.h>
#include <multigrid.h>
#include <smoother.h>
#include <quocMatrices.h>

namespace mg {

inline bool isInside ( const int _x, const int _y, const int _z, const int _n ) {
  return ( ( _x >= 0 ) && ( _x < _n ) && ( _y >= 0 ) && ( _y < _n ) && ( _z >= 0 ) && ( _z < _n ) );
}


/** Onthefly-coarsening of 3D quoc matrices */
template< class MatrixType, class GridType >
void coarsenQuocMatrix3D ( const MatrixType &fineMat, MatrixType &coarseMat,
                           const GridType &fineGrid, const GridType &coarseGrid ) {

  typedef typename MatrixType::DataType RealType;

  // this may be redundant
  static const int _quocNeighborOffsets[27][3] = { { -1, -1, -1 },
                                                   {  0, -1, -1 },
                                                   {  1, -1, -1 },
                                                   { -1,  0, -1 },
                                                   {  0,  0, -1 },
                                                   {  1,  0, -1 },
                                                   { -1,  1, -1 },
                                                   {  0,  1, -1 },
                                                   {  1,  1, -1 },
                                                   { -1, -1,  0 },
                                                   {  0, -1,  0 },
                                                   {  1, -1,  0 },
                                                   { -1,  0,  0 },
                                                   {  0,  0,  0 },
                                                   {  1,  0,  0 },
                                                   { -1,  1,  0 },
                                                   {  0,  1,  0 },
                                                   {  1,  1,  0 },
                                                   { -1, -1,  1 },
                                                   {  0, -1,  1 },
                                                   {  1, -1,  1 },
                                                   { -1,  0,  1 },
                                                   {  0,  0,  1 },
                                                   {  1,  0,  1 },
                                                   { -1,  1,  1 },
                                                   {  0,  1,  1 },
                                                   {  1,  1,  1 } } ;

  // this may be redundant
  static const RealType _quocWeight[27] = { 0.125 ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.5   ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.5   ,
                                            0.25  ,
                                            0.5   ,
                                            1.0   ,
                                            0.5   ,
                                            0.25  ,
                                            0.5   ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.5   ,
                                            0.25  ,
                                            0.125 ,
                                            0.25  ,
                                            0.125  };



  const int nfi = fineGrid.getWidth(), nco = coarseGrid.getWidth();

  aol::ProgressBar<> pb ( "Coarsening matrix" );
  pb.start ( aol::Cub ( nco ) );

  qc::OTFILexMapper<qc::QC_3D>  coarseMapper ( nco ), fineMapper( nfi );

  // todo: optimize these nested loops.
  // todo: implement case of Dirichlet boundary conditions

  // \todo use brick iterator here

  for ( int iz = 0; iz < nco; ++iz ) {
    for ( int iy = 0; iy < nco; ++iy ) {
      for ( int ix = 0; ix < nco; ++ix ) {
        pb++;
        const int i = coarseMapper.getGlobalIndex ( ix, iy, iz ); // global node index in coarse grid

        if ( isInside ( ix, iy, iz, nco ) ) {

          for ( int _l = 0; _l < 27; ++_l ) {
            const int
              lx = ix + _quocNeighborOffsets[_l][0],
              ly = iy + _quocNeighborOffsets[_l][1],
              lz = iz + _quocNeighborOffsets[_l][2],
              l  = coarseMapper.getGlobalIndex ( lx, ly, lz );

            if ( isInside ( lx, ly, lz, nco ) ) {

              for ( int _j = 0; _j < 27; ++_j ) {
                const int
                  jx = 2 * ix + _quocNeighborOffsets[_j][0] ,
                  jy = 2 * iy + _quocNeighborOffsets[_j][1] ,
                  jz = 2 * iz + _quocNeighborOffsets[_j][2] ,
                  j  = fineMapper.getGlobalIndex ( jx, jy, jz );

                if ( isInside ( jx, jy, jz, nfi ) ) {

                  for ( int _k = 0; _k < 27; ++_k ) {
                    const int
                      kx = 2 * lx + _quocNeighborOffsets[_k][0] ,
                      ky = 2 * ly + _quocNeighborOffsets[_k][1] ,
                      kz = 2 * lz + _quocNeighborOffsets[_k][2] ,
                      k  = fineMapper.getGlobalIndex ( kx, ky, kz );

                    if ( isInside ( kx, ky, kz, nfi ) ) {

                      const RealType
                        re_ij = _quocWeight[ _j ] ,    // 1, if coarse and fine point coincide, 1/2, 1/4, 1/8  if interpolation point nearby
                        fi_jk = fineMat.get ( j, k ) ,
                        pr_kl = _quocWeight[ _k ] ;

                      coarseMat.add ( i, l, re_ij * fi_jk * pr_kl );

                    } // k dof
                  } // _k
                } // j dof
              } // _j
            } // l dof
          } // _l
        } // i dof
      }  // ix
    } // iy
  } // iz

  pb.finish();

}


// this seems to be inefficient compared to assembling restriction / prolongation operators and matrix-multiplication for coarsening operators


/** A general multigrid solver for quoc computations (multilinear finite elements on quoc grids) with Neumann (and without Dirichlet) boundary conditions.
 *  The matrix on the finest grid must be given, it is coarsened automatically.
 *  This solver is useful (i. e. efficient) only if the coarse matrices cannot be generated directly.
 *  \author Schwen
 *  \todo allow using other MatrixType than UniformGridSparseMatrix
 *  \todo do not specify CMode as template parameter, but pass it in the constructor
 */
template< typename RealType, typename MatrixType, CoarseningMode CMode = MATRIXMULT_COARSENING >
class QuocMultigrid : public mg::ExplicitOperatorHierarchyMultigrid < aol::Vector<RealType>,
                                                                      MatrixType,
                                                                      mg::GaussSeidelSmoother < aol::Vector<RealType>, MatrixType > ,
                                                                      qc::RestrictOp< RealType, qc::STD_MG_RESTRICT >,
                                                                      qc::ProlongOp< RealType >,
                                                                      qc::GridDefinition > {

public:
  typedef qc::GridDefinition      GridType;
  typedef aol::Vector<RealType>   VectorType;

public:
  QuocMultigrid ( const MatrixType &fineMatrix,
                  const GridType &fineGrid,
                  const int ExplicitLevel = 2,
                  const int PreSmooth = 3, const int PostSmooth = 3,
                  const RealType Relax = 1.0,
                  const mg::MGCycleMode CycleMode = mg::MULTIGRID_V_CYCLE,
                  const RealType CycleEps = 1.0e-16,
                  const int MaxCycle = 100,
                  const aol::StoppingMode stop = aol::STOPPING_UNSET  ) : mg::ExplicitOperatorHierarchyMultigrid < aol::Vector<RealType>,
                                                                                                                   MatrixType,
                                                                                                                   mg::GaussSeidelSmoother < aol::Vector<RealType>, MatrixType > ,
                                                                                                                   qc::RestrictOp< RealType, qc::STD_MG_RESTRICT >,
                                                                                                                   qc::ProlongOp< RealType >,
                                                                                                                   qc::GridDefinition > ( fineMatrix, fineGrid, ExplicitLevel, PreSmooth, PostSmooth, Relax, CycleMode, CycleEps, MaxCycle, stop ) {

    if ( fineGrid.getDimOfWorld() != 3 ) {
      throw aol::UnimplementedCodeException ( "mg::QuocMultigrid: not implemented for dimension != 3", __FILE__, __LINE__ );
    }

    this->init();

  }

public:
  virtual void coarsenGrid ( const int coarseLevel ) {
    GridType *pCoarseGrid = new GridType ( coarseLevel, this->getGridRef ( coarseLevel + 1 ).getDimOfWorld() );
    this->setpGrid ( coarseLevel, pCoarseGrid );
  }

  virtual void coarsenOperator ( const int coarseLevel ) {
    typedef aol::SparseMatrix<RealType> SparseMatrixType;

    // this will become the operator on coarseLevel
    MatrixType *re_X_op_X_pr = new MatrixType ( this->getGridRef ( coarseLevel ) );

    switch ( CMode ) {
    case ONTHEFLY_COARSENING: {

      mg::coarsenQuocMatrix3D ( this->getOperatorRef ( coarseLevel + 1 ), ( *re_X_op_X_pr ),
                                  this->getGridRef ( coarseLevel + 1 ), this->getGridRef ( coarseLevel ) );

    } break;

    case MATRIXMULT_COARSENING: {

      const int
        fineLevel   = coarseLevel + 1,
        coarseSize  = this->getGridRef ( coarseLevel ).getNumberOfNodes(),
        fineSize    = this->getGridRef ( fineLevel ).getNumberOfNodes();

      cerr << " left-multiplying by restriction matrix ... ";

      SparseMatrixType *re_mat = new SparseMatrixType ( coarseSize, fineSize );
      this->getRestrictionOpRef ( coarseLevel ).assembleAddMatrix ( ( *re_mat ) );

      SparseMatrixType *re_X_fineop = new SparseMatrixType ( coarseSize, fineSize );

      ( *re_X_fineop ).addMatrixProduct ( ( *re_mat ), this->getOperatorRef ( fineLevel ) ); // i. e. get matrix

      delete ( re_mat );

      cerr << " right-multiplying by prolongation matrix ... ";

      SparseMatrixType *pr_mat  = new SparseMatrixType ( fineSize, coarseSize );
      this->getProlongationOpRef ( coarseLevel ).assembleAddMatrix ( ( *pr_mat ) );

      ( *re_X_op_X_pr ).addMatrixProduct ( ( *re_X_fineop ), ( *pr_mat ) );

      delete ( re_X_fineop );
      delete ( pr_mat );
    } break;

    default:
      throw aol::Exception ( "Illegal coarsening mode", __FILE__, __LINE__ );

    }

    this->setpOperator ( coarseLevel, re_X_op_X_pr );

    // not deleted here (on purpose): re_X_op_X_pr
    // this one is deleted in the destructor

  }

  virtual void setRestrictionOperator ( const int level ) {
    this->setpRestrictionOperator ( level, new qc::RestrictOp<RealType, qc::STD_MG_RESTRICT>( this->getGridRef ( level ), this->getGridRef ( level+1 ), aol::ONTHEFLY) );
    // ONTHEFLY_COARSENING only applies to operator coarsening. The restriction operator should be onthefly.
  }

  virtual void setProlongationOperator( const int level ) {
    this->setpProlongationOperator( level, new qc::ProlongOp<RealType>( this->getGridRef ( level ), this->getGridRef ( level + 1 ), aol::ONTHEFLY ) );
  }

}; // end of class mg::QuocMultigrid


}

#endif
