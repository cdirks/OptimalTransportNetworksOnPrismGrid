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

#ifndef __TPCFE_ELASTOP_H
#define __TPCFE_ELASTOP_H

#include "tpcfe_utils.h"

// here, we try to implement an elasticity ooperator for the tpcfe framework.

// one based on TP's mixed derivative ops has been moved to the modules directory, this one is for testing.

namespace tpcfe {

template <typename RealType>
class CFE_Elast_Lookup {
public:
  CFE_Elast_Lookup() {

    createElastLookups();

  };

public:
  RealType
  _LmixedM_StdTetra[6][3][3][4][4],
  _LEM_StdTetra[6][3][3][4][4],
  _LmixedM_StdCube[3][3][8][8],
  _LEM_StdCube[3][3][8][8];

  void createElastLookups() {

    // _LmixedM_StdTera and _LEM_StdTetra
    for ( int t = 0; t < 6; ++t ) {
      // D Phi:
      aol::Matrix33<double> dphi, dphiinv;
      for ( int m = 0; m < 3; ++m ) {
        for ( int n = 0; n < 3; ++n ) {
          dphi.set ( m, n, tpcfe::CFELookup<RealType>::_hexNbOffs[ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][m+1] ][ n ] -  tpcfe::CFELookup<RealType>::_hexNbOffs[ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][0] ][ n ] );
        };
      };

      double det_dphi = aol::Abs ( dphi.det() );

      dphiinv = dphi.inverse();

      aol::Vec3<double> grad_phi[4], def_grad_phi[4];
      for ( int i = 0; i < 4; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          grad_phi[i][j] = tpcfe::CFELookup<RealType>::_gradientStdTetraBasisFct[i][j];
        };
        def_grad_phi[i] = dphiinv * grad_phi[i];
      };

      for ( int a = 0; a < 3; ++a ) {
        for ( int b = 0; b < 3; ++b ) {
          for ( int i = 0; i < 4; ++i ) {
            for ( int j = 0; j < 4; ++j ) {

              _LmixedM_StdTetra[t][a][b][i][j]  = det_dphi * 1.0 * ( def_grad_phi[i][a] * def_grad_phi[j][b] );

              _LEM_StdTetra    [t][a][b][i][j]  = det_dphi * 0.5 * ( def_grad_phi[i][b] * def_grad_phi[j][a] );
              if ( a == b ) {
                _LEM_StdTetra  [t][a][b][i][j] += det_dphi * 0.5 * ( def_grad_phi[i][0] * def_grad_phi[j][0] );
                _LEM_StdTetra  [t][a][b][i][j] += det_dphi * 0.5 * ( def_grad_phi[i][1] * def_grad_phi[j][1] );
                _LEM_StdTetra  [t][a][b][i][j] += det_dphi * 0.5 * ( def_grad_phi[i][2] * def_grad_phi[j][2] );
              };

            };
          };
        };
      };
    };

    for ( int a = 0; a < 3; ++a ) {
      for ( int b = 0; b < 3; ++b ) {
        for ( int c = 0; c < 8; ++c ) {
          for ( int d = 0; d < 8; ++d ) {
            _LmixedM_StdCube[a][b][c][d] = 0.0;
            _LEM_StdCube[a][b][c][d] = 0.0;
          };
        };
        for ( int i = 0; i < 4; ++i ) {
          for ( int j = 0; j < 4; ++j ) {
            for ( int t = 0; t < 6; ++t ) {
              _LmixedM_StdCube[a][b][ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][i] ][ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][j] ] += _LmixedM_StdTetra[t][a][b][i][j];
              _LEM_StdCube    [a][b][ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][i] ][ tpcfe::CFELookup<RealType>::_stdTetraVertex[t][j] ] += _LEM_StdTetra    [t][a][b][i][j];
            };
          };
        };
      };

    };
  }

};



enum BlockEntry { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };

template < typename ConfiguratorType >
class CFEElastOp_Single : public CFEStandardOp < ConfiguratorType, CFEElastOp_Single < ConfiguratorType > > {

public:
  typedef typename ConfiguratorType::RealType RealType;

public:
  // for now, we explicitely modify matrices for the nonDOF nodes, thus operators are assembled in constructor.

  CFEElastOp_Single ( const typename ConfiguratorType::GridType &Grid, const BlockEntry block, const RealType lambda, const RealType mu ) :
      CFEStandardOp < ConfiguratorType,  CFEElastOp_Single< ConfiguratorType > > ( Grid, aol::ASSEMBLED ),
      _block ( block ),
      _lambda ( lambda ),
      _mu ( mu ),
      _a ( 3 ), _b ( 3 ),
      _grid ( Grid ) {

    switch ( _block ) {
    case XX:
      _a = 0; _b = 0; break;
    case XY:
      _a = 0; _b = 1; break;
    case XZ:
      _a = 0; _b = 2; break;
    case YX:
      _a = 1; _b = 0; break;
    case YY:
      _a = 1; _b = 1; break;
    case YZ:
      _a = 1; _b = 2; break;
    case ZX:
      _a = 2; _b = 0; break;
    case ZY:
      _a = 2; _b = 1; break;
    case ZZ:
      _a = 2; _b = 2; break;
    default:
      cerr << "This should not happen" << __FILE__ << " " << __LINE__ << endl;
    };

    this->init();

    this->assembleMatrix();

  }

public:
  // TODO: make this nicer.
  CFEElastOp_Single ( const typename ConfiguratorType::GridType &Grid, const BlockEntry block ) {
    _block = block;

    switch ( _block ) {
    case XX:
      _a = 0; _b = 0; break;
    case XY:
      _a = 0; _b = 1; break;
    case XZ:
      _a = 0; _b = 2; break;
    case YX:
      _a = 1; _b = 0; break;
    case YY:
      _a = 1; _b = 1; break;
    case YZ:
      _a = 1; _b = 2; break;
    case ZX:
      _a = 2; _b = 0; break;
    case ZY:
      _a = 2; _b = 1; break;
    case ZZ:
      _a = 2; _b = 2; break;
    default:
      cerr << "This should not happen" << __FILE__ << " " << __LINE__ << endl;
    };

  };

public:
  /** Create the standard matrices of the non-interfaced standard tetrahedra
   */
  void createDefaultTetraMatrices() {
    CFE_Elast_Lookup<RealType> lookup;

    const RealType scale = this->_grid.H() ; // sure?

    for ( int l = 0; l < 6; l++ ) {
      for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
          this->_defaultLocalTetraMatrix[l][i][j] = ( _lambda * lookup._LmixedM_StdTetra[l][_a][_b][i][j] + 2.0 * _mu * lookup._LEM_StdTetra[l][_a][_b][i][j] ) * scale ; // scale? // same as below!

        }
      }
    }
  }

  /** fill the matrix for a sub-tetrahedron
   */
  void preprocessLocalTetraMatrix ( const CFEElement< RealType > &/*e*/, const CFETetra<RealType> & t ) const {
    RealType local_mixed_m[4][4], local_etens_m[4][4];

    for ( int i = 0; i < 4; ++i ) {
      for ( int j = 0; j < 4; ++j ) {
        local_mixed_m[i][j] = 0.0;
        local_etens_m[i][j] = 0.0;
      };
    };

    // Prepare stiffness matrix
    t.computeInverseTransformation();

    aol::Vec3<RealType> grad_ref_bf[4], def_grad_ref_bf[4]; // gradients of basis functions on reference tetrahedron
    for ( int i = 0; i < 4; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        grad_ref_bf[i][j] = tpcfe::CFELookup<RealType>::_gradientStdTetraBasisFct[i][j];
      };
      def_grad_ref_bf[i] = t.barycentricCoordinates * grad_ref_bf[i];
    };

    for ( int i = 0; i < 4; ++i ) {
      for ( int j = 0; j < 4; ++j ) {
        local_mixed_m[i][j] +=     1.0 * def_grad_ref_bf[i][_a] * def_grad_ref_bf[j][_b];

        local_etens_m[i][j] +=     0.5 * def_grad_ref_bf[i][_b] * def_grad_ref_bf[j][_a];
        if ( _a == _b ) {
          local_etens_m[i][j]  +=  0.5 * def_grad_ref_bf[i][ 0] * def_grad_ref_bf[j][ 0] ;
          local_etens_m[i][j]  +=  0.5 * def_grad_ref_bf[i][ 1] * def_grad_ref_bf[j][ 1] ;
          local_etens_m[i][j]  +=  0.5 * def_grad_ref_bf[i][ 2] * def_grad_ref_bf[j][ 2] ;
        };
      };
    };

    // Compute scaling factor
#if 0
    const RealType det = t.determinant();
    const RealType volume = t.getVolume();
    const RealType scale = ( this->_grid.H() / aol::Sqr ( det ) ) * volume; // explanation: barycentricCoordinates is determinant times inverse, another factor det comes from the transformation theorem
#endif
    // this is what Tobias uses.
    const RealType det = t.determinant();
    const RealType scale = this->_grid.H() / ( 6.0 * det );

    // OS TODO: think about whether this is the correct thing to do. Why ignore degenerate tetrahedra of nonzero volume? Which tetrahedra are these? Maybe those outside? But why does the error happen at computational nodes? At boundary?

    if ( det != 0.0 ) { // !qc::isNan( scale ) ){ // is this correct? if so, why?
      for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
          this->_localTetraMatrix[i][j] = ( _lambda * local_mixed_m[i][j] + 2.0 * _mu * local_etens_m[i][j] ) * scale;; // scale ??
        };
      };
    } else {
      cerr << "Problem: detected degenerate tetrahedron: " << t.edge[0] << " " << t.edge[1] << " " << t.edge[2] << endl;
    };

  }

  // TODO: encapsulate this operator in CFEDDRestrOp instead of modifying matrices for Dirichlet/non-domain entries

public:

  RealType getMatrixEntry ( const int i, const int j ) const {
    return ( this->_mat->get ( i, j ) );
  }

  void restrict_NonDOFEntries() {
    //      cerr << _a << " " << _b << " " << static_cast<RealType>( _a == _b ) << endl;
    restrictNonDOFEntries ( this->_grid, * ( this->_mat ), static_cast<RealType> ( _a == _b ) ); // diagonal blocks: 1, off-diagonal blocks: 0
  }

  void restrictNonDomainEntries() {
    //      cerr << _a << " " << _b << " " << static_cast<RealType>( _a == _b ) << endl;
    restrictNonDomainEntries ( this->_grid, * ( this->_mat ), static_cast<RealType> ( _a == _b ) );
  }

  void restrictDirichletEntries() {
    //      cerr << _a << " " << _b << " " << static_cast<RealType>( _a == _b ) << endl;
    restrictDirichletEntries ( this->_grid, * ( this->_mat ), static_cast<RealType> ( _a == _b ) );
  }


protected:
  BlockEntry _block;
  RealType _lambda, _mu;
  int _a, _b;
  const typename ConfiguratorType::GridType &_grid;
};


template < typename T_ConfiguratorType >
class MyCFEElastOp : public aol::BlockOp< typename T_ConfiguratorType::RealType, CFEElastOp_Single< T_ConfiguratorType > > {

public:
  typedef T_ConfiguratorType                    ConfiguratorType;
  typedef typename ConfiguratorType::RealType   RealType;
  typedef qc::MultiArray<RealType, 3>          VectorType;
  typedef typename ConfiguratorType::GridType   GridType;

public:
  MyCFEElastOp ( const typename ConfiguratorType::GridType &Grid, const RealType lambda, const RealType mu ) :
      aol::BlockOp< typename ConfiguratorType::RealType, CFEElastOp_Single< ConfiguratorType > > ( 3, 3 ),
      _lambda ( lambda ),
      _mu ( mu ) {

    CFEElastOp_Single< ConfiguratorType >
      *op_xx = new CFEElastOp_Single< ConfiguratorType > ( Grid, XX, _mu, _lambda ),
      *op_xy = new CFEElastOp_Single< ConfiguratorType > ( Grid, XY, _mu, _lambda ),
      *op_xz = new CFEElastOp_Single< ConfiguratorType > ( Grid, XZ, _mu, _lambda ),
      *op_yx = new CFEElastOp_Single< ConfiguratorType > ( Grid, YX, _mu, _lambda ),
      *op_yy = new CFEElastOp_Single< ConfiguratorType > ( Grid, YY, _mu, _lambda ),
      *op_yz = new CFEElastOp_Single< ConfiguratorType > ( Grid, YZ, _mu, _lambda ),
      *op_zx = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZX, _mu, _lambda ),
      *op_zy = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZY, _mu, _lambda ),
      *op_zz = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZZ, _mu, _lambda );

    this->setReference ( 0, 0, *op_xx );   this->setReference ( 0, 1, *op_xy );   this->setReference ( 0, 2, *op_xz );
    this->setReference ( 1, 0, *op_yx );   this->setReference ( 1, 1, *op_yy );   this->setReference ( 1, 2, *op_yz );
    this->setReference ( 2, 0, *op_zx );   this->setReference ( 2, 1, *op_zy );   this->setReference ( 2, 2, *op_zz );

  };

public:
  MyCFEElastOp ( const typename ConfiguratorType::GridType &Grid, aol::OperatorType OpType ) : aol::BlockOp< typename ConfiguratorType::RealType, CFEElastOp_Single< ConfiguratorType > > ( 3, 3 ) {
    CFEElastOp_Single< ConfiguratorType >
      *op_xx = new CFEElastOp_Single< ConfiguratorType > ( Grid, XX ),
      *op_xy = new CFEElastOp_Single< ConfiguratorType > ( Grid, XY ),
      *op_xz = new CFEElastOp_Single< ConfiguratorType > ( Grid, XZ ),
      *op_yx = new CFEElastOp_Single< ConfiguratorType > ( Grid, YX ),
      *op_yy = new CFEElastOp_Single< ConfiguratorType > ( Grid, YY ),
      *op_yz = new CFEElastOp_Single< ConfiguratorType > ( Grid, YZ ),
      *op_zx = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZX ),
      *op_zy = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZY ),
      *op_zz = new CFEElastOp_Single< ConfiguratorType > ( Grid, ZZ );

    this->setReference ( 0, 0, *op_xx );   this->setReference ( 0, 1, *op_xy );   this->setReference ( 0, 2, *op_xz );
    this->setReference ( 1, 0, *op_yx );   this->setReference ( 1, 1, *op_yy );   this->setReference ( 1, 2, *op_yz );
    this->setReference ( 2, 0, *op_zx );   this->setReference ( 2, 1, *op_zy );   this->setReference ( 2, 2, *op_zz );

  };

public:
  ~MyCFEElastOp() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        if ( this->getPointer ( i, j ) ) {
          delete ( this->getPointer ( i, j ) );
        };
      };
    };
  };

public:

  void allocateMatrix() {
    this->get ( 0, 0 )->allocateMatrix(); // TODO: continue!
  };

  void restrict_NonDOFEntries() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        if ( this->getPointer ( i, j ) ) {
          this->get ( i, j )->restrict_NonDOFEntries();
        };
      };
    };
  }

  void restrictNonDomainEntries() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        if ( this->getPointer ( i, j ) ) {
          this->getPointer ( i, j )->restrictNonDomainEntries();
        };
      };
    };
  }

  void restrictDirichletEntries() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        if ( this->getPointer ( i, j ) ) {
          this->getPointer ( i, j )->restrictDirichletEntries();
        };
      };
    };
  }

protected:
  RealType _lambda, _mu;

}; // end class


} // end namespace

#endif
