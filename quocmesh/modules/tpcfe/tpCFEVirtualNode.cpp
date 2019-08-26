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

#include <tpCFEVirtualNode.h>
#include <tpCFEGrid.h>
#include <tpCFEUtils.h>
#include <matrixInverse.h>

namespace tpcfe {

template < typename RealType >
CFEGeomVirtualNode< RealType >::CFEGeomVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) {
  for ( short i = 0; i < 3; ++i ) {
    _coord[i] = tn[i];
    _normal[i] = no[i];
  }
}


template < typename RealType >
typename CFEGeomVirtualNode< RealType >::IndexType CFEGeomVirtualNode< RealType >::mapIndex ( const int node0,  const int node1 ) {
  const int shift = sizeof ( int ) * 8;
#ifdef BOUNDS_CHECK
  if ( node0 < -1 || node1 < -1 ) {
    cerr << node0 << " " << node1 << endl;
    throw aol::Exception ( "tpcfe::CFEGeomVirtualNode<RealType>::mapIndex: node indices must be >= -1 (corresponding to other node regular)", __FILE__, __LINE__ );
  }
#endif
  if ( node0 == -1 ) return ( static_cast<IndexType> ( node1 ) );
  if ( node1 == -1 ) return ( static_cast<IndexType> ( node0 ) );

  IndexType Node0 = node0, Node1 = node1;
  if ( node0 > node1 ) { // this is different from mapIndex
    return ( ( Node0 << shift ) | Node1 );
  } else {
    return ( ( Node1 << shift ) | Node0 );
  }
}


template < typename RealType >
void CFEGeomVirtualNode< RealType >::splitIndex ( const IndexType index, int &node0, int &node1 ) {
  const int shift = sizeof ( node0 ) * 8;
  node0 = static_cast< int > ( index );          // lower shift bits
  node1 = static_cast< int > ( index >> shift ); // upper shift bits
}


template < typename RealType >
void CFEGeomVirtualNode< RealType >::addInElement ( const CFEElement<RealType> &inEl, const int ParentTetra ) {
  typename PTMapType::iterator vnIt = _parentTetra.find ( CFETopoElement ( inEl ) );

  if ( vnIt == _parentTetra.end() )  {   // This hexahedron has not been added before
    ParentTetraList pt;
    pt.addPT ( ParentTetra );
    _parentTetra.insert ( pair<CFETopoElement, ParentTetraList> ( CFETopoElement ( inEl ), pt ) );
  } else {                            // It has been added before
    for ( int i = 0; i < vnIt->second.getNumPT(); ++i ) {     // Check whether parent tetra has been added before
      if ( ( vnIt->second ) [i] == ParentTetra ) return;
    }
    vnIt->second.addPT ( ParentTetra );
  }
}


template < typename RealType >
void CFEGeomVirtualNode< RealType >::determineGlobalCoordinates ( const CFEElement<RealType> &el ) {
  aol::Vec3<int> _gPos0, _gPos1;
  for ( short l = 0; l < 3; ++l ) {
    _gPos0[l] = el[l] + static_cast<int> ( CFELookup<RealType>::_hexNbOffs[_local0][l] );
  }

  if ( _local1 != NODE_NOT_VIRTUAL ) {
    for ( short l = 0; l < 3; ++l ) {
      _gPos1[l] = el[l] + static_cast<int> ( CFELookup<RealType>::_hexNbOffs[_local1][l] );
    }
    for ( short l = 0; l < 3; ++l ) {
      this->_coord[l] = ( el.getCutRelation ( _local0, _local1 ) * _gPos1[l] + el.getCutRelation ( _local1, _local0 ) * _gPos0[l] );
    }
  }
}


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
void CFEVirtualNodeBase< RealType, CT, NodalCoeffType >::determineApproximateNormal ( const CFEStructure<RealType> &s ) {

  this->_normal.setZero();

  // Run through all elements which contain this virtual node
  for ( typename PTMapType::iterator it = this->_parentTetra.begin(); it != this->_parentTetra.end(); ++it ) {
    const CFETopoElement &sparseEl = it->first;
    CFEElement<RealType> fullEl ( sparseEl );

    const int numConstrainingTetra = it->second.getNumPT();

    for ( int i = 0; i < 8; ++i )
      fullEl.setStructureValue ( i, s.getValue ( fullEl.getGlobalPos ( i ) ) );

    // For each parent tetrahedron constraining this node
    for ( int i = 0; i < numConstrainingTetra; ++i ) {
      const int parentTetra = ( it->second ) [i];

      // For each vertex of this tetra
      for ( int j = 0; j < 4; ++j ) {
        const int vertexIndex = CFELookup<RealType>::_stdTetraVertex[parentTetra][j];
        const RealType value = fullEl.getStructureValue ( vertexIndex );
        aol::Vec3<RealType> localNormal = CFELookup<RealType>::_stdTetraBasisGradients[parentTetra][j];
        localNormal *= value;
        this->_normal += localNormal;
      }
    }
  }

  // above, we have added all tetra-wise gradients (non-normalized normals),
  // averaging is done implicitly when normalizing now
  this->_normal.normalize();

}


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
void CFEVirtualNodeBase< RealType, CT, NodalCoeffType >::determineConstraintsCD ( const CFEGridBase< RealType, CT, NodalCoeffType > &grid ) {
  const CFETopoElement &sparseEl = this->getInElement();
  CFEElement<RealType> fullEl ( sparseEl );

  this->_numOfConstraints = 2;
  this->_weights.reallocate ( 2 );
  this->_constraints.reallocate ( 2 );

  this->_constraints[0] = this->_global0;
  this->_constraints[1] = this->_global1;

  grid.getCutRelations ( fullEl );

  const RealType cut0 = fullEl.getCutRelation ( this->_local0, this->_local1 ), cut1 = fullEl.getCutRelation ( this->_local1, this->_local0 );

  this->_weights[0] = aol::ZOTrait<ExtrapolationWeightType>::one; // opposite way around
  this->_weights[0] *= cut1;
  this->_weights[1] = aol::ZOTrait<ExtrapolationWeightType>::one;
  this->_weights[1] *= cut0;
}


template < typename RealType, tpcfe::ConstraintType CT >
void CFEVirtualNodeScalar < RealType, CT >::determineConstraintsLIEHR ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs,
                                                                        const CFEGridBase< RealType, CT, RealType > &grid ) {

  const CFETopoElement &sparseEl = this->getInElement();
  CFEElement<RealType> fullEl ( sparseEl );

  this->_numOfConstraints = 2;
  this->_weights.reallocate ( 2 );
  this->_constraints.reallocate ( 2 );

  this->_constraints[0] = this->_global0;
  this->_constraints[1] = this->_global1;

  grid.getCutRelations ( fullEl );

  const RealType cut0    = fullEl.getCutRelation ( this->_local0, this->_local1 );
  const RealType cut1    = fullEl.getCutRelation ( this->_local1, this->_local0 );
  const RealType coeff0  = nodalCoeffs[ this->_global0 ];
  const RealType coeff1  = nodalCoeffs[ this->_global1 ];
  const RealType weight1 = cut0 / ( cut0 + cut1 * coeff0 / coeff1 );
  const RealType weight0 = 1 - weight1;

#ifdef VERBOSE
  cerr << "cut0 = " << cut0 << ", coeff0 = " << coeff0 << ", _global0 = " << this->_global0 << ", weight0 = " << weight0 << ";  " << "cut1 = " << cut1 << ", coeff1 = " << coeff1 << ", _global1 = " << this->_global1 << ", weight1 = " << weight1 << endl;
#endif

  this->_weights[0] = weight0;
  this->_weights[1] = weight1;
}



template < typename RealType, tpcfe::ConstraintType CT >
void CFEVirtualNodeScalar < RealType, CT >::determineConstraintsTPOS ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs,
                                                                       const CFEGridBase< RealType, CT, RealType > &grid ) {
  map<int, RealType>   weightmap;
  int                  globalNumConstrainingTetra = 0;

  aol::Vec<4, RealType> localWeight;
  aol::Vec<4, int>      localConstraint;

  // For each hexahedral element containing this node
  for ( typename PTMapType::iterator it = this->_parentTetra.begin(); it != this->_parentTetra.end(); ++it ) {
    const int numConstrainingTetra = it->second.getNumPT();

    // For each parent tetrahedron constraining this node
    for ( int i = 0; i < numConstrainingTetra; ++i ) {
      const int parentTetra = ( it->second ) [i];
      const CFETopoElement &sparseEl = it->first;

      // Do the local construction
      bool inversionReliable = false;
      this->computeTPOSLocalWeightsAndConstraints ( nodalCoeffs, sparseEl, parentTetra, localWeight, localConstraint, inversionReliable );

      if ( inversionReliable ) {
        ++globalNumConstrainingTetra;
        // Add the weights of all local constructions up
        for ( int j = 0; j < 4; ++j ) {
          weightmap[ localConstraint[j] ] += localWeight[j];
        }
      } else {
#ifdef VERBOSE
        qc::CoordType pos0, pos1;
        grid.getIndexMapperRef().splitGlobalIndex ( this->_global0, pos0 );
        grid.getIndexMapperRef().splitGlobalIndex ( this->_global1, pos1 );
        cerr << "VN " << pos0 << " - " << pos1
        << ", parent tetra " << parentTetra
        << " in element " << sparseEl
        << " not considered reliable for weight construction" << endl;
#else
        aol::doNothingWithArgumentToPreventUnusedParameterWarning ( grid );
#endif
      }
    }
  }

  // We do not want to store constraints with weights exactly equal to zero. Thresholding better?
  for ( typename map<int, RealType>::iterator mit = weightmap.begin(); mit != weightmap.end(); /*increment must be done below*/ ) {
#ifdef VERBOSE
    cerr << "Dropping constraint " << mit->first << " with zero extrapolation weight" << endl;
#endif
    if ( mit->second == aol::NumberTrait<RealType>::zero ) {
      weightmap.erase ( mit++ );
    } else {
      ++mit;
    }
  }

  // Divide by the number of constraining parent tetrahedra and store the result

  this->_numOfConstraints = weightmap.size();
  this->_weights.reallocate ( this->_numOfConstraints );
  this->_constraints.reallocate ( this->_numOfConstraints );

  typename map<int, RealType>::iterator mit = weightmap.begin();
  for ( int i = 0; mit != weightmap.end(); ++mit, ++i ) {
    this->_weights[i] = mit->second / globalNumConstrainingTetra;
    this->_constraints[i] = mit->first;
  }
}


template < typename RealType, tpcfe::ConstraintType CT  >
void CFEVirtualNodeScalar < RealType, CT >::getTPOSLocalMatrix ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs,
                                                                 const CFETopoElement &sparseEl,
                                                                 const int parentTetra,
                                                                 aol::Vec<4, int> &localConstraint,
                                                                 aol::Matrix44<RealType> &locMat ) {

  if ( CT != CFE_TPOS && CT != CFE_CDWI_TPOS ) {
    throw aol::Exception ( "CFEVirtualNode::getTPOSLocalMatrix may not be called for CT != CFE_TPOS", __FILE__, __LINE__ );
  }

  aol::Vec<4, int>          sign;
  aol::Vec3<RealType>       vnToVertex[4];
  aol::Vec3<RealType>       tangent[2];

  CFEElement<RealType>      fullEl ( sparseEl );

  tpcfe::CFEWeightProvider<RealType, RealType> weightProvider ( nodalCoeffs, fullEl );

  RealType diffCoeffPlus = weightProvider.meanPlusWeight();
  RealType diffCoeffMinus = weightProvider.meanMinusWeight();
  RealType quotient = diffCoeffMinus / diffCoeffPlus;

  // Compute the vertices relative to the virtual node
  for ( short i = 0; i < 4; ++i ) {
    const int vertexIndex = CFELookup<RealType>::_stdTetraVertex[parentTetra][i];
    const int globalIndex = fullEl.globalIndex ( vertexIndex );

    vnToVertex[i] =  aol::Vec3<RealType> ( sparseEl );
    vnToVertex[i] += aol::Vec3<RealType> ( CFELookup<RealType>::_hexNbOffs[vertexIndex] );
    vnToVertex[i] -= this->_coord;

    sign[i] = weightProvider.getSign ( vertexIndex );
    localConstraint[i] = globalIndex;
  }

  aol::Vec3<RealType> normalUsed;
  if ( CT == CFE_TPOS || CT == CFE_CDWI_TPOS ) {
    normalUsed = this->_normal;
  } else {
    throw aol::Exception ( "tpcfe::CFEVirtualNode<RealType, CT>::getTPOSLocalMatrix should not be called for CT != TPOS, CDWI_TPOS", __FILE__, __LINE__ );
  }

  // Compute the tangent vectors
  {
    short j0 = 0;
    for ( short i = 1; i < 3; ++i ) {
      if ( fabs ( normalUsed[i] ) > fabs ( normalUsed[j0] ) ) {
        j0 = i;
      }
    }

    tangent[0][ ( j0+1 ) %3 ] = -normalUsed[j0];
    tangent[0][   j0        ] =  normalUsed[ ( j0+1 ) %3];
    tangent[0][ ( j0+2 ) %3 ] =  aol::NumberTrait<RealType>::zero;

    tangent[0].normalize();

    tangent[1] = normalUsed.crossProduct ( tangent[0] );

    QUOC_ASSERT ( ( fabs ( normalUsed.norm() - 1 ) < 1.0e-6 ) && ( fabs ( tangent[0].norm() - 1 ) < 1.0e-6 ) && ( fabs ( tangent[1].norm() - 1 ) < 1.0e-6 ) );
  }

#ifdef VERBOSE
  {
    cerr << "Virtual Node at " << this->_coord << ", dc-/dc+ = " << quotient << endl
         << "At position: " << sparseEl[0] << " " << sparseEl[1] << " " << sparseEl[2] << endl
         << "Edges (VN to vertex): " << vnToVertex[0] << ", " << vnToVertex[1] << ", " << vnToVertex[2] << ", " << vnToVertex[3] << endl
         << "NormalUsed: " << normalUsed << " first tangent: " << tangent[0] << " second tangent: " << tangent[1] << endl
         << normalUsed * tangent[0] << endl
         << normalUsed * tangent[1] << endl
         << tangent[0] * tangent[1] << endl << endl
         << aol::detailedFormat ( normalUsed * vnToVertex[0] ) << ", " << aol::detailedFormat ( normalUsed * vnToVertex[1] ) << ", " << aol::detailedFormat ( normalUsed * vnToVertex[2] ) << endl
         << aol::detailedFormat ( tangent[0] * vnToVertex[0] ) << ", " << aol::detailedFormat ( tangent[0] * vnToVertex[1] ) << ", " << aol::detailedFormat ( tangent[0] * vnToVertex[2] ) << endl
         << aol::detailedFormat ( tangent[1] * vnToVertex[0] ) << ", " << aol::detailedFormat ( tangent[1] * vnToVertex[1] ) << ", " << aol::detailedFormat ( tangent[1] * vnToVertex[2] ) << endl;
  }
#endif

  // Now set the matrix entries
  for ( int i = 0; i < 4; ++i ) {
    locMat.set ( i, 0, aol::NumberTrait<RealType>::one );
    locMat.set ( i, 1, ( sign[i] == + 1 ? quotient : aol::NumberTrait<RealType>::one ) *  normalUsed * vnToVertex[i] );
    locMat.set ( i, 2, tangent[0] * vnToVertex[i] );
    locMat.set ( i, 3, tangent[1] * vnToVertex[i] );
  }

}


template < typename RealType, tpcfe::ConstraintType CT >
void CFEVirtualNodeScalar < RealType, CT >::computeTPOSLocalWeightsAndConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs,
                                                                                    const CFETopoElement &sparseEl,
                                                                                    const int parentTetra,
                                                                                    aol::Vec<4, RealType> &localWeight, aol::Vec<4, int> &localConstraint,
                                                                                    bool &inversionReliable ) {

  aol::Matrix44<RealType>   matrix, coeffMatrix;

  getTPOSLocalMatrix ( nodalCoeffs, sparseEl, parentTetra, localConstraint, matrix );

  try {
    coeffMatrix.makeInverse ( matrix );
#ifdef VERBOSE
    {
      aol::Matrix44<RealType> test = matrix;
      test *= coeffMatrix;
      cerr << "Matrix " << endl << matrix << endl;
      cerr << "coeffMatrix " << endl << coeffMatrix << endl;
      cerr << "Testing for Identity: " << endl << test << endl;
      aol::Matrix44<RealType> ident;
      ident.setIdentity();
      test -= ident;
      cerr << "Frobenius norm of difference to identity = " << test.norm() << endl;

      if ( false ) {
        aol::FullMatrix<RealType> fmat ( 4, 4 );
        for ( int i = 0; i < 4; ++i )
          for ( int j = 0; j < 4; ++j )
            fmat.set ( i, j, matrix.get ( i, j ) );
        aol::computeConditionNumberViaOctave<aol::FullMatrix<RealType>, RealType > ( fmat );
      }
    }
#endif
  } catch ( aol::Exception &e ) {
    e.dump();
    matrix.print ( cerr );
  }

  for ( short i = 0; i < 4; ++i ) {
    localWeight[i] = coeffMatrix.get ( 0, i );
  }

  // check whether Frobenius norm of matrix times its inverse is sufficiently small
  aol::Matrix44<RealType> test = matrix;
  aol::Matrix44<RealType> ident;
  test *= coeffMatrix;
  ident.setIdentity();
  test -= ident;
  inversionReliable = ( test.norm() < this->_maximalInversionError );
}



template < typename RealType >
void CFEVirtualNode<RealType, CFE_NONE, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &/*nodalCoeffs*/, const CFEGridBase< RealType, CFE_NONE, RealType > &/*grid*/ ) {
  throw aol::Exception ( "Illegal CFE constraint type", __FILE__, __LINE__ );
}


template < typename RealType, typename NodalCoeffType >
void CFEVirtualNode<RealType, CFE_CD, NodalCoeffType >::determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &/*nodalCoeffs*/, const CFEGridBase< RealType, CFE_CD, NodalCoeffType > &grid ) {
  this->determineConstraintsCD ( grid );
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_LIEHR, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_LIEHR, RealType > &grid ) {
  this->determineConstraintsLIEHR ( nodalCoeffs, grid );
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_TPOS, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_TPOS, RealType > &grid ) {
  this->determineConstraintsTPOS ( nodalCoeffs, grid );
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_TP, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_TP, RealType > &/*grid*/ ) {
  this->_numOfConstraints = 4;
  this->_weights.reallocate ( 4 );
  this->_constraints.reallocate ( 4 );

  aol::Vec<4, RealType> localWeight;
  aol::Vec<4, int>      localConstraint;

  const int parentTetra = this->getParentTetra();

  bool iR;
  this->computeTPOSLocalWeightsAndConstraints ( nodalCoeffs, this->getInElement(), parentTetra, localWeight, localConstraint, iR );

  for ( int i = 0; i < 4; ++i ) {
    this->_weights[i] = localWeight[i];
    this->_constraints[i] = localConstraint[i];
  }
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_DOMAIN, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &/*nodalCoeffs*/, const CFEGridBase< RealType, CFE_DOMAIN, RealType > &grid ) {
  this->_weights.reallocate ( 1 );
  this->_constraints.reallocate ( 1 );
  this->_weights[0] =  1.0;
  // We set the constraint to a negavtive value, since this is
  // the way how the matrix assembler distinguishes regular from virtual nodes
  this->_constraints[0] = -1 * dynamic_cast< const CFEGrid<RealType,CFE_DOMAIN>& >( grid ).getFreeGlobalDof ( this ); // todo: change this ...
  this->_numOfConstraints = 1;
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_CDWI_LIEHR, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_CDWI_LIEHR, RealType > &grid ) {
  const CFETopoElement &sparseEl = this->getInElement();
  const CFEType type = sparseEl.cfeType();
  const int structureNo = type._structureNo;
  if ( structureNo == MAX_STRUCT_ID ) {
    // now we are in the CFE_CD case
    this->determineConstraintsCD ( grid );
  } else if ( structureNo < MAX_STRUCT_ID ) {
    this->determineConstraintsLIEHR ( nodalCoeffs, grid );
  } else {
    throw aol::Exception ( "Cannot find this structure", __FILE__, __LINE__ );
  }
  this->postprocessDetermineConstraints();
}


template < typename RealType >
void CFEVirtualNode<RealType, CFE_CDWI_TPOS, RealType >::determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_CDWI_TPOS, RealType > &grid ) {
  const CFETopoElement &sparseEl = this->getInElement();
  const CFEType type = sparseEl.cfeType();
  const int structureNo = type._structureNo;
  if ( structureNo == MAX_STRUCT_ID ) {
    // now we are in the CFE_CD case
    this->determineConstraintsCD ( grid );
  } else if ( structureNo < MAX_STRUCT_ID ) {
    this->determineConstraintsTPOS ( nodalCoeffs, grid );
  } else {
    throw aol::Exception ( "Cannot find this structure", __FILE__, __LINE__ );
  }
  this->postprocessDetermineConstraints();
}


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
void CFEVirtualNodeVector < RealType, CT, NodalCoeffType >::extrapolate ( const qc::MultiArray< RealType, qc::QC_3D > &data, aol::Vec3< RealType > &dest ) const {
  dest.setZero();

  for ( int i = 0; i < this->_numOfConstraints; ++i ) {
    dest += this->weight ( i ) * data.get ( this->constraint ( i ) ) ; // Matrix33 * Vec3
  }
}


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
void CFEVirtualNodeVector < RealType, CT, NodalCoeffType >::determineConstraintsTPOSELAST ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CT, NodalCoeffType > &grid ) {
  if ( CT != CFE_TPOSELAST && CT != CFE_CDWI_ELAST ) {
    throw aol::Exception ( "Illegal constraint type", __FILE__, __LINE__ );
  }

  map<int, ExtrapolationWeightType> weightmap;
  int globalNumConstrainingTetra = 0;

  // For each hexahedral element containing this node
  for ( typename PTMapType::iterator it = this->_parentTetra.begin(); it != this->_parentTetra.end(); ++it ) {
    const int numConstrainingTetra = it->second.getNumPT();

    // For each parent tetrahedron constraining this node
    for ( int cT = 0; cT < numConstrainingTetra; ++cT ) {
      const int parentTetra = ( it->second ) [cT];
      const CFETopoElement &sparseEl = it->first;

      // Do the local construction
      bool inversionReliable = true;

      aol::RandomAccessContainer< ExtrapolationWeightType > localWeight ( 4 );
      aol::RandomAccessContainer< int > localConstraint ( 4 );

      const aol::Vec3<RealType> N = this->_normal;
      aol::Vec3<RealType> S, T;
      aol::Mat<4,3,RealType> vnToVertexTimesNST;
      aol::Vec<4, signed char> sign;

      // compute difference between regular nodes (of tetrahedron) and virtual node
      const CFEElement<RealType> fullEl ( sparseEl );
      CFEWeightProvider< RealType, NodalCoeffType> weightProvider ( nodalCoeffs, fullEl );

      { // geometric information [OK]
        aol::Mat<4,3,RealType> vnToVertex;

        for ( short v = 0; v < 4; ++v ) {
          const int vertexIndex = CFELookup<RealType>::_stdTetraVertex[parentTetra][v];
          const int globalIndex = fullEl.globalIndex ( vertexIndex );

          vnToVertex[v] =  aol::Vec3<RealType> ( sparseEl );
          vnToVertex[v] += aol::Vec3<RealType> ( CFELookup<RealType>::_hexNbOffs[vertexIndex] );
          vnToVertex[v] -= this->_coord;

          sign[v] = weightProvider.getSign ( vertexIndex );
          localConstraint[v] = globalIndex;
        }

        // Compute the tangent vectors
        short j0 = 0;
        for ( short i = 1; i < 3; ++i ) {
          if ( fabs ( N[i] ) > fabs ( N[j0] ) ) {
            j0 = i;
          }
        }

        S [ (j0+1)%3 ] = -N[j0];
        S [   j0     ] =  N[ (j0+1)%3 ];
        S [ (j0+2)%3 ] =  aol::NumberTrait<RealType>::zero;

        S.normalize();

        T = N.crossProduct ( S );

        QUOC_ASSERT ( ( fabs ( N.norm() - 1 ) < 1.0e-6 ) && ( fabs ( S.norm() - 1 ) < 1.0e-6 ) && ( fabs ( T.norm() - 1 ) < 1.0e-6 ) && ( fabs ( N * S ) < 1.0e-6 ) && ( fabs ( N * T ) < 1.0e-6 ) && ( fabs ( S * T ) < 1.0e-6 ));

        for ( short v = 0; v < 4; ++v ) {
          vnToVertexTimesNST[v] = aol::Vec3<RealType> ( vnToVertex[v] * N, vnToVertex[v] * S, vnToVertex[v] * T );
        }
      }

      aol::FullMatrix<RealType> couplingMat ( 9, 9 );

      { // computation of the coupling matrix [OK]
        aol::Matrix33<RealType> Amat;
        Amat[0][0] = N[0];   Amat[0][1] = N[1];   Amat[0][2] = N[2];
        Amat[1][0] = S[0];   Amat[1][1] = S[1];   Amat[1][2] = S[2];
        Amat[2][0] = T[0];   Amat[2][1] = T[1];   Amat[2][2] = T[2];

        aol::Matrix33<RealType> KN[2], KS[2], KT[2]; // 0: plus, 1: minus
        RealType CTensor[2][3][3][3][3];
        weightProvider.meanPlusWeight().getAnisotropicTensor( CTensor[0] );
        weightProvider.meanMinusWeight().getAnisotropicTensor( CTensor[1] );

        for ( short sign = 0; sign < 2; ++sign ) {
          for ( short i = 0; i < 3; ++i ) {
            for ( short l = 0; l < 3; ++l ) {
              for ( short j = 0; j < 3; ++j ) {
                for ( short k = 0; k < 3; ++k ) {
                  KN[sign][i][l] += CTensor[sign][i][j][k][l] * Amat[0][k] * N[j] ;
                  KS[sign][i][l] += CTensor[sign][i][j][k][l] * Amat[1][k] * N[j] ;
                  KT[sign][i][l] += CTensor[sign][i][j][k][l] * Amat[2][k] * N[j] ;
                }
              }
            }
          }
        }

        aol::FullMatrix<RealType> KNSTPlus ( 9, 9 ), KNSTMinus ( 9, 9 );
        for ( short i = 0; i < 3; ++i ) {
          for ( short j = 0; j < 3; ++j ) {
            KNSTPlus.set ( i,   j, KN[0][i][j] );    KNSTMinus.set ( i,   j, KN[1][i][j] );
            KNSTPlus.set ( i, 3+j, KS[0][i][j] );    KNSTMinus.set ( i, 3+j, KS[1][i][j] );
            KNSTPlus.set ( i, 6+j, KT[0][i][j] );    KNSTMinus.set ( i, 6+j, KT[1][i][j] );
          }
        }
        for ( short i = 3; i < 9; ++i ) {
          KNSTPlus.set ( i, i, 1 );     KNSTMinus.set ( i, i, 1 );
        }

        aol::FullMatrix<RealType> KNSTPlusInv ( 9, 9 );
        KNSTPlusInv.makeInverse ( KNSTPlus );
        couplingMat.makeProduct ( KNSTPlusInv, KNSTMinus );

        // compute reliability coefficient
        aol::FullMatrix<RealType> test ( 9, 9 ), ident ( 9, 9 );
        test = KNSTPlus;
        test *= KNSTPlusInv;
        ident.setIdentity();
        test -= ident;
        inversionReliable &= ( test.getFrobeniusNormSqr() < this->_maximalInversionError );
      }

#ifdef DEBUG
      { // check against alternative implementation
        aol::FullMatrix<RealType> pmMat[2];  pmMat[0].reallocate ( 9, 9 ); pmMat[1].reallocate ( 9, 9 );

        const aol::Vec3<RealType> RMat[3] = { N, S, T }; // i. e. row-wise N, S, T
        RealType GTensor[3][3][3][3], CTensor[2][3][3][3][3];

        for ( short i = 0; i < 3; ++i )
          for ( short j = 0; j < 3; ++j )
            for ( short k = 0; k < 3; ++k )
              for ( short l = 0; l < 3; ++l )
                GTensor[i][j][k][l] = 0.5 * ( ( i == k ) * RMat[l][j] + ( j == k ) * RMat[l][i] );

        weightProvider.meanPlusWeight().getAnisotropicTensor( CTensor[0] );
        weightProvider.meanMinusWeight().getAnisotropicTensor( CTensor[1] );

        for ( short pm = 0; pm < 2; ++pm ) {
          for ( short i = 0; i < 3; ++i ) {
            for ( short m = 0; m < 3; ++m ) {
              for ( short n = 0; n < 3; ++n ) {
                RealType dwarf = 0;
                for ( short j = 0; j < 3; ++j ) {
                  for ( short k = 0; k < 3; ++k ) {
                    for ( short l = 0; l < 3; ++l ) {
                      dwarf += CTensor[pm][i][j][k][l] * GTensor[k][l][m][n] * N[j];
                    }
                  }
                }
                pmMat[pm].add ( i, m+3*n, dwarf ); // partial derivatives in col-wise order
              }
            }
          }

          for ( short i = 3; i < 9; ++i )
            pmMat[pm].set ( i, i, 1 );
        }

        aol::FullMatrix<RealType> plusMatInv ( 9, 9 ), checkCouplingMat ( 9, 9 );
        plusMatInv.makeInverse ( pmMat[0] );
        checkCouplingMat.makeProduct ( plusMatInv, pmMat[1] );

        checkCouplingMat -= couplingMat;
        if ( checkCouplingMat.getFrobeniusNormSqr() > 1.0e-24 )
          cerr << "Check of computing coupling matrix vs. other implementation failed! " << checkCouplingMat.getFrobeniusNormSqr() << endl;


        // 2nd check: check whether coupling mat describes the right thing
        for ( short eta = 0; eta < 9; ++eta ) {
          aol::Matrix33<RealType> epsilonPlus, epsilonMinus, NST;
          NST.setRow ( 0, N );      NST.setRow ( 1, S );      NST.setRow ( 2, T );

          {
            aol::Matrix33<RealType> tmp;
            switch ( eta ) {
            case 0:
            case 1:
            case 2:
              tmp.setRow ( eta, N );
              break;
            case 3:
            case 4:
            case 5:
              tmp.setRow ( eta-3, S );
              break;
            case 6:
            case 7:
            case 8:
              tmp.setRow ( eta-6, T );
              break;
            default: // this cannot happen as eta ranges from 0 to 8
              throw aol::UnimplementedCodeException ( "eta out of bounds", __FILE__, __LINE__ );
            }

            epsilonMinus += tmp;
            epsilonMinus += tmp.transposed();
            epsilonMinus *= 0.5;
          }

          {
            aol::Matrix33<RealType> tmp, tmp2;
            tmp.set ( 0, 0, couplingMat.get(0,eta) );      tmp.set ( 0, 1, couplingMat.get(3,eta) );      tmp.set ( 0, 2, couplingMat.get(6,eta) );
            tmp.set ( 1, 0, couplingMat.get(1,eta) );      tmp.set ( 1, 1, couplingMat.get(4,eta) );      tmp.set ( 1, 2, couplingMat.get(7,eta) );
            tmp.set ( 2, 0, couplingMat.get(2,eta) );      tmp.set ( 2, 1, couplingMat.get(5,eta) );      tmp.set ( 2, 2, couplingMat.get(8,eta) );
            tmp2.makeProduct ( tmp, NST );
            epsilonPlus += tmp2;
            epsilonPlus += tmp2.transposed();
            epsilonPlus *= 0.5;
          }

          aol::Matrix33<RealType> sigmaMinus, sigmaPlus;
          for ( short i = 0; i < 3; ++i ) {
            for ( short j = 0; j < 3; ++j ) {
              for ( short k = 0; k < 3; ++k ) {
                for ( short l = 0; l < 3; ++l ) {
                  sigmaMinus.add ( i, j, CTensor[1][i][j][k][l] * epsilonMinus.get ( k, l ) );
                  sigmaPlus.add ( i, j, CTensor[0][i][j][k][l] * epsilonPlus.get ( k, l ) );
                }
              }
            }
          }

          if ( ( sigmaMinus * N - sigmaPlus * N ).norm() > 1.0e-14 * aol::Max ( ( sigmaMinus * N ).norm(), aol::NumberTrait<RealType>::one ) )
            cerr << "Problems determining coupling matrix: " << endl << sigmaMinus * N << endl << " and " << endl << sigmaPlus * N << endl << " differ too much." << endl;
        }
      }
#endif

      if ( inversionReliable ) { // use coupling conditions to determine the extrapolation weights [OK]; else don't waste time doing the work ...

        // for row and column indices in the matrix rcI[l][n][k] and vectors rcI[i][j][n]
        static const unsigned char rcI[4][3][3] = { { {  0,  1,  2 }, {  3,  4,  5 }, {  6,  7,  8 } },
                                                    { {  9, 10, 11 }, { 12, 13, 14 }, { 15, 16, 17 } },
                                                    { { 18, 19, 20 }, { 21, 22, 23 }, { 24, 25, 26 } },
                                                    { { 27, 28, 29 }, { 30, 31, 32 }, { 33, 34, 35 } } };
        // for entries from the constraint matrix
        static const unsigned char cMatI[3][3] = { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 } };
        aol::FullMatrix<RealType> conMat ( 36, 36 );
        aol::Vector<RealType> conSoln ( 36 ), conRhs ( 36 );

        short i = -1;
        for ( short j = 0; j < 3; ++j ) {
          for ( short l = 0; l < 4; ++l ) {
            for ( short n = 0; n < 3; ++n ) {
              for ( short k = 0; k < 3; ++k ) {
                i = 0; {
                  switch ( sign[l] ) {
                  case -1: conMat.set ( rcI[i][j][n], rcI[l][n][k], ( j == k ? 1 : 0 ) * vnToVertexTimesNST[l][i] ); break;
                  case +1: conMat.set ( rcI[i][j][n], rcI[l][n][k], couplingMat.get ( cMatI[i][k], j ) * vnToVertexTimesNST[l][i] ); break;
                  default: throw aol::Exception ( "Illegal sign", __FILE__, __LINE__ );
                  }
                }
                for ( i = 1; i < 3; ++i ) {
                  switch ( sign[l] ) {
                  case -1: conMat.set ( rcI[i][j][n], rcI[l][n][k], ( j == k ? 1 : 0 ) * vnToVertexTimesNST[l][i] ); break;
                  case +1: conMat.set ( rcI[i][j][n], rcI[l][n][k], ( j == k ? 1 : 0 ) * vnToVertexTimesNST[l][i] + couplingMat.get ( k, cMatI[i][j] ) * vnToVertexTimesNST[l][0] ); break;
                  default: throw aol::Exception ( "Illegal sign", __FILE__, __LINE__ );
                  }
                }
                i = 3; {
                  conMat.set ( rcI[i][j][n], rcI[l][n][k], ( j == k ? aol::ZOTrait<RealType>::one : aol::ZOTrait<RealType>::zero ) );
                  conRhs.set( rcI[i][j][n], ( j == n ? aol::ZOTrait<RealType>::one : aol::ZOTrait<RealType>::zero ) );
                }
              }
            }
          }
        }

        aol::LUInverse<RealType> conSolver ( conMat );
        conSolver.apply ( conRhs, conSoln );

        { // check reliability
          aol::Vector<RealType> dummy (36);
          conMat.apply ( conSoln, dummy );
          dummy -= conRhs;

          inversionReliable &= dummy.norm() < this->_maximalInversionError;
        }

        for ( short l = 0; l < 4; ++l )
          for ( short n = 0; n < 3; ++n )
            for ( short k = 0; k < 3; ++k )
              localWeight[l].set ( n, k, conSoln[ rcI[l][n][k] ] );


#ifdef DEBUG
        {
          aol::FullMatrix<RealType> checkConMat ( 36, 36 );
          aol::Vector<RealType> checkConSoln ( 36 ), checkConRhs ( 36 );

          // this is how it is written out
          aol::FullMatrix<RealType> bMC ( 9, 9 );
          bMC.set ( 0, 0, 1.0 );        bMC.set ( 1, 3, 1.0 );        bMC.set ( 2, 6, 1.0 );
          bMC.set ( 3, 1, 1.0 );        bMC.set ( 4, 4, 1.0 );        bMC.set ( 5, 7, 1.0 );
          bMC.set ( 6, 2, 1.0 );        bMC.set ( 7, 5, 1.0 );        bMC.set ( 8, 8, 1.0 );

          // cerr << "bMC: " << endl << bMC << endl << endl;

          aol::FullMatrix<RealType> bMK[3];
          for ( short i = 0; i < 3; ++i )
            bMK[i].reallocate ( 9, 9 );

          for ( short i = 0; i < 3; ++i ) {
            bMK[0].set ( 0+i, 3*i+0, couplingMat.get( 0, 0 ) );          bMK[0].set ( 0+i, 3*i+1, couplingMat.get( 1, 0 ) );          bMK[0].set ( 0+i, 3*i+2, couplingMat.get( 2, 0 ) );
            bMK[0].set ( 3+i, 3*i+0, couplingMat.get( 0, 1 ) );          bMK[0].set ( 3+i, 3*i+1, couplingMat.get( 1, 1 ) );          bMK[0].set ( 3+i, 3*i+2, couplingMat.get( 2, 1 ) );
            bMK[0].set ( 6+i, 3*i+0, couplingMat.get( 0, 2 ) );          bMK[0].set ( 6+i, 3*i+1, couplingMat.get( 1, 2 ) );          bMK[0].set ( 6+i, 3*i+2, couplingMat.get( 2, 2 ) );

            bMK[1].set ( 0+i, 3*i+0, couplingMat.get( 0, 3 ) );          bMK[1].set ( 0+i, 3*i+1, couplingMat.get( 1, 3 ) );          bMK[1].set ( 0+i, 3*i+2, couplingMat.get( 2, 3 ) );
            bMK[1].set ( 3+i, 3*i+0, couplingMat.get( 0, 4 ) );          bMK[1].set ( 3+i, 3*i+1, couplingMat.get( 1, 4 ) );          bMK[1].set ( 3+i, 3*i+2, couplingMat.get( 2, 4 ) );
            bMK[1].set ( 6+i, 3*i+0, couplingMat.get( 0, 5 ) );          bMK[1].set ( 6+i, 3*i+1, couplingMat.get( 1, 5 ) );          bMK[1].set ( 6+i, 3*i+2, couplingMat.get( 2, 5 ) );

            bMK[2].set ( 0+i, 3*i+0, couplingMat.get( 0, 6 ) );          bMK[2].set ( 0+i, 3*i+1, couplingMat.get( 1, 6 ) );          bMK[2].set ( 0+i, 3*i+2, couplingMat.get( 2, 6 ) );
            bMK[2].set ( 3+i, 3*i+0, couplingMat.get( 0, 7 ) );          bMK[2].set ( 3+i, 3*i+1, couplingMat.get( 1, 7 ) );          bMK[2].set ( 3+i, 3*i+2, couplingMat.get( 2, 7 ) );
            bMK[2].set ( 6+i, 3*i+0, couplingMat.get( 0, 8 ) );          bMK[2].set ( 6+i, 3*i+1, couplingMat.get( 1, 8 ) );          bMK[2].set ( 6+i, 3*i+2, couplingMat.get( 2, 8 ) );
          }


          for ( short j = 0; j < 4; ++j ) {
            for ( int i = 0; i < 1; ++i ) {
              aol::FullMatrix<RealType> currentBlock ( 9, 9 );
              switch ( sign[j] ) {
              case +1:
                currentBlock = bMK[i];
                break;
              case -1:
                currentBlock = bMC;
                break;
              default:
                throw aol::Exception ("Illegal sign", __FILE__, __LINE__ );
              }
              currentBlock *= vnToVertexTimesNST[j][i];

              checkConMat.setBlock ( i*9, j*9, currentBlock );
            }
            for ( int i = 1; i < 3; ++i ) {
              aol::FullMatrix<RealType> currentBlock ( 9, 9 );
              // in all cases:
              currentBlock = bMC;
              currentBlock *= vnToVertexTimesNST[j][i];
              switch ( sign[j] ) {
              case +1:
                {
                  aol::FullMatrix<RealType> dwarf = bMK[i];
                  dwarf *= vnToVertexTimesNST[j][0];
                  currentBlock += dwarf;
                }
                break;
              case -1:
                // no additional term
                break;
              default:
                throw aol::Exception ("Illegal sign", __FILE__, __LINE__ );
              }

              checkConMat.setBlock ( i*9, j*9, currentBlock );
            }
            for ( int i = 3; i < 4; ++i ) {
              checkConMat.setBlock ( i*9, j*9, bMC );
            }
          }
          checkConRhs[27] = 1;        checkConRhs[31] = 1;        checkConRhs[35] = 1;

          aol::LUInverse<RealType> checkConSolver ( checkConMat );
          checkConSolver.apply ( checkConRhs, checkConSoln );

          aol::RandomAccessContainer< ExtrapolationWeightType > checkLocalWeight ( 4 );
          for ( int v = 0, c = 0; v < 4; ++v )
            for ( int i = 0; i < 3; ++i )
              for ( int j = 0; j < 3; ++j, ++c )
                checkLocalWeight[v].set ( i, j, conSoln[c] );

          for ( short l = 0; l < 4; ++l ) {
            checkLocalWeight[l] -= localWeight[l];
            if ( checkLocalWeight[l].norm() > 1.0e-14 ) {
              cerr << "Problems determining construction weights: " << endl << checkLocalWeight[l] << endl << " and " << endl << localWeight[l] << endl << " differ too much." << endl;
            }
          }
        }
#endif

#ifdef DEBUG
        {
          aol::FullMatrix<RealType> checkConMat ( 36, 36 );
          aol::Vector<RealType> checkConSoln ( 36 ), checkConRhs ( 36 );

          aol::Vec<3,RealType> lr[4][3];       // indices: [b][d] [c]       type of spanning function b, spatial direction of spanning function d, spatial component c
          aol::Vec<3,RealType> lv[4][3][4]; // indices: [b][d][v] [e] type of spanning function b, spatial direction of spanning function d, vertex (regular node) v; spatial component e
          // note: c is not necessary because the entries are the same for all c(!!)
          for ( short b = 0; b < 3; ++b ) {
            // b: eta in directions N, S, T
            for ( short d = 0; d < 3; ++d ) {
              aol::Mat<3,3,RealType> kLocMat;
              for ( short kli = 0; kli < 3; ++kli ) {
                for ( short klj = 0; klj < 3; ++klj ) {
                  kLocMat.set ( kli, klj, couplingMat.get ( 3 * klj + kli, 3 * b + d ) );
                }
              }

              for ( short v = 0; v < 4; ++v ) {
                if ( sign[v] == +1 ) {
                  const aol::Vec<3,RealType> &dummy = vnToVertexTimesNST[v];
                  lv[b][d][v] = kLocMat * dummy;
                } else if ( sign[v] == -1 ) {
                  lv[b][d][v][d] = vnToVertexTimesNST[v][b];
                } else {
                  throw aol::Exception ("Illegal sign", __FILE__, __LINE__ );
                }
                // lr[b][d][c(=d)] remains zero
              }

            }
          }

          { // constant eta;
            const short b = 3;
            for ( short d = 0; d < 3; ++d ) {
              for ( short v = 0; v < 4; ++v ) {
                const short e = d;
                lv[b][d][v][e] = aol::NumberTrait<RealType>::one;
              }
              const short c = d;
              lr[b][d][c] = aol::NumberTrait<RealType>::one;
            }
          }

          for ( short b = 0; b < 4; ++b ) {
            for ( short d = 0; d < 3; ++d ) {
              for ( short c = 0; c < 3; ++c ) {
                const short cI = 9 * b + 3 * d + c;

                // entry in vector
                checkConRhs.set ( cI, lr[b][d][c] );

                // row in matrix
                for ( short v = 0; v < 4; ++v ) {
                  for ( short e = 0; e < 3; ++e ) {
                    const short cJ = 9 * v + 3 * c + e;
                    checkConMat.set ( cI, cJ, lv[b][d][v][e] );
                  }
                }
              }
            }
          }

          aol::LUInverse<RealType> checkConSolver ( checkConMat );
          checkConSolver.apply ( checkConRhs, checkConSoln );

          aol::RandomAccessContainer< ExtrapolationWeightType > checkLocalWeight ( 4 );
          for ( int v = 0, c = 0; v < 4; ++v )
            for ( int i = 0; i < 3; ++i )
              for ( int j = 0; j < 3; ++j, ++c )
                checkLocalWeight[v].set ( i, j, checkConSoln[c] );

          for ( short l = 0; l < 4; ++l ) {
            checkLocalWeight[l] -= localWeight[l];
            if ( checkLocalWeight[l].norm() > 1.0e-14 ) {
              cerr << "Problems determining construction weights: " << endl << checkLocalWeight[l] << endl << " and " << endl << localWeight[l] << endl << " differ too much." << endl;
            }
          }
        }
#endif
      }


      // Store local weights if tetra considered reliable

      if ( inversionReliable ) {
        ++globalNumConstrainingTetra;
        // Add the weights of all local constructions up
        for ( int j = 0; j < 4; ++j ) {
          weightmap[ localConstraint[j] ] += localWeight[j];
        }
      } else {
#ifdef VERBOSE
        qc::CoordType pos0, pos1;
        grid.getIndexMapperRef().splitGlobalIndex ( this->_global0, pos0 );
        grid.getIndexMapperRef().splitGlobalIndex ( this->_global1, pos1 );
        cerr << "VN " << pos0 << " - " << pos1
             << ", parent tetra " << parentTetra
             << " in element " << sparseEl
             << " not considered reliable for weight construction" << endl;
#else
        aol::doNothingWithArgumentToPreventUnusedParameterWarning ( grid );
#endif
      }
    }
  }

  // We do not want to store constraints with weights exactly equal to zero. Maybe could also use threshold?
  for ( typename map<int, ExtrapolationWeightType>::iterator mit = weightmap.begin(); mit != weightmap.end(); /*increment must be done below*/ ) {
    if ( mit->second == aol::ZOTrait<ExtrapolationWeightType>::zero ) {
#ifdef VERBOSE
      cerr << "Dropping constraint " << mit->first << " with extrapolation weight " << mit->second << endl;
#endif
      weightmap.erase ( mit++ );
    } else {
      ++mit;
    }
  }

  // Divide by the number of constraining parent tetrahedra and store the result

  this->_numOfConstraints = weightmap.size();
  this->_weights.reallocate ( this->_numOfConstraints );
  this->_constraints.reallocate ( this->_numOfConstraints );

  typename map<int, ExtrapolationWeightType>::iterator mit = weightmap.begin();
  for ( int i = 0; mit != weightmap.end(); ++mit, ++i ) {
    ExtrapolationWeightType tmpW = mit->second;
    tmpW /= static_cast<RealType>( globalNumConstrainingTetra );
    this->_weights[i] = tmpW;
    this->_constraints[i] = mit->first;
  }

  this->postprocessDetermineConstraints();
}


template < typename RealType, typename NodalCoeffType >
void CFEVirtualNode<RealType, CFE_TPOSELAST, NodalCoeffType >::determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_TPOSELAST, NodalCoeffType > &grid ) {
  CFEVirtualNodeVector<RealType, CFE_TPOSELAST, NodalCoeffType >::determineConstraintsTPOSELAST ( nodalCoeffs,  grid );
}


template < typename RealType, typename NodalCoeffType >
void CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType >::determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase<RealType, CFE_CDWI_ELAST, NodalCoeffType> &grid ) {
  const CFETopoElement &sparseEl = this->getInElement();
  const CFEType type = sparseEl.cfeType();
  const int structureNo = type._structureNo;
  if ( structureNo == MAX_STRUCT_ID ) {
    this->determineConstraintsCD ( grid );
  } else if ( structureNo < MAX_STRUCT_ID ) {
    this->determineConstraintsTPOSELAST ( nodalCoeffs, grid );
  } else {
    throw aol::Exception ( "Cannot find this structure", __FILE__, __LINE__ );
  }
  this->postprocessDetermineConstraints();
}


template class CFEGeomVirtualNode < float >;
template class CFEGeomVirtualNode < double  >;
template class CFEGeomVirtualNode < long double >;


template<> float CFEGeomVirtualNode< float >::_maximalInversionError = aol::NumberTrait<float>::getInf();
template<> double CFEGeomVirtualNode< double >::_maximalInversionError = aol::NumberTrait<double>::getInf();
template<> long double CFEGeomVirtualNode< long double >::_maximalInversionError = aol::NumberTrait<long double>::getInf();


#define INSTANTIATE_CFEClass2(ClassName, CT)\
template class ClassName < float,       CT >;\
template class ClassName < double,      CT >;\
template class ClassName < long double, CT >;

#define INSTANTIATE_CFEClass3(ClassName, CT)\
template class ClassName < float,       CT, float       >;\
template class ClassName < double,      CT, double      >;\
template class ClassName < long double, CT, long double >;

#define INSTANTIATE_CFEClassElast32Arg(ClassName, CT, ElastCoeffT)\
template class ClassName < float,       CT, ElastCoeffT<float>       >;\
template class ClassName < double,      CT, ElastCoeffT<double>      >;\
template class ClassName < long double, CT, ElastCoeffT<long double> >;

#define INSTANTIATE_CFEClassElast3(ClassName, CT)\
INSTANTIATE_CFEClassElast32Arg(ClassName, CT, tpcfe::IsotropicElasticityCoefficient)\
INSTANTIATE_CFEClassElast32Arg(ClassName, CT, tpcfe::VoigtElasticityCoefficient)


INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_NONE )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_TP )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_TPOS )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_LIEHR )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_DOMAIN )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_CD )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_CDWI_TPOS )
INSTANTIATE_CFEClass3 ( CFEVirtualNodeBase, CFE_CDWI_LIEHR )

INSTANTIATE_CFEClassElast3 ( CFEVirtualNodeBase, CFE_CD )
INSTANTIATE_CFEClassElast3 ( CFEVirtualNodeBase, CFE_TPOSELAST )
INSTANTIATE_CFEClassElast3 ( CFEVirtualNodeBase, CFE_CDWI_ELAST )

INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_NONE )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_TP )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_TPOS )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_LIEHR )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_DOMAIN )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_CD )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_CDWI_TPOS )
INSTANTIATE_CFEClass2( CFEVirtualNodeScalar, CFE_CDWI_LIEHR )

// no CFE_CD here!
INSTANTIATE_CFEClassElast3 ( CFEVirtualNodeVector, CFE_TPOSELAST )
INSTANTIATE_CFEClassElast3 ( CFEVirtualNodeVector, CFE_CDWI_ELAST )

INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_NONE )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_TP )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_TPOS )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_LIEHR )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_DOMAIN )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_CD )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_CDWI_TPOS )
INSTANTIATE_CFEClass3 ( CFEVirtualNode, CFE_CDWI_LIEHR )

INSTANTIATE_CFEClassElast3 ( CFEVirtualNode, CFE_CD )
INSTANTIATE_CFEClassElast3 ( CFEVirtualNode, CFE_TPOSELAST )
INSTANTIATE_CFEClassElast3 ( CFEVirtualNode, CFE_CDWI_ELAST )

#undef INSTANTIATE_CFEClass3
#undef INSTANTIATE_CFEClassElast32Arg
#undef INSTANTIATE_CFEClassElast3

}
