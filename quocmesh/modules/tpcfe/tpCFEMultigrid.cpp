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

#include <tpCFEMultigrid.h>

namespace tpcfe {

template< typename GridType >
void SlopeInterfaceBase<GridType>::dumpSlopes ( ostream &os ) const {
  os << "Dumping slopes:" << endl;
  for ( qc::RectangularIterator<qc::QC_3D> pit ( _slopeStorage ); pit.notAtEnd(); ++pit ) {
    if ( _slopeStorage.getRef ( *pit ).size() > 0 )
      os << "Position " << *pit << endl;
    for ( typename std::map< EdgeIndexType, SlopesType >::const_iterator mappos = _slopeStorage.getRef ( *pit ).begin(); mappos != _slopeStorage.getRef ( *pit ).end(); ++mappos ) {
      int n1, n2;
      splitEdgeIndex ( mappos->first, n1, n2 );
      qc::CoordType pn1, pn2;
      _iMap.splitGlobalIndex ( n1, pn1 );
      _iMap.splitGlobalIndex ( n2, pn2 );

      const RealType outSl = getOutSlope ( *pit, pn1, pn2 ), arrSl = getArrSlope ( *pit, pn1, pn2 );
      os << pn1 << " -- " << pn2 << ": " << aol::detailedFormat ( outSl  ) << " " << aol::detailedFormat ( arrSl );
      if ( ( fabs ( outSl - arrSl ) < 1.0e-6 ) && ( fabs ( fabs ( outSl ) - aol::NumberTrait<RealType>::one ) < 1.0e-6 ) ) { // standard radial edge
        os << " (standard)";
#ifdef DEBUG
        if ( ( ( fabs ( outSl - aol::NumberTrait<RealType>::one ) < 1.0e-6 ) &&  ( *pit != pn2 ) ) || ( ( fabs ( outSl + aol::NumberTrait<RealType>::one ) < 1.0e-6 ) &&  ( *pit != pn1 ) ) )
          os << " but incorrect";
#endif
      } else if ( outSl * arrSl < 0 ) { // interfaced radial edge
        os << " (with sign change)" << aol::detailedFormat ( outSl / arrSl );
#ifdef DEBUG
        if ( ( *pit == pn1 ) || ( *pit == pn2 ) ) {
          os << " but radial: " << *pit << " " << pn1 << " " << pn2;
        }
#endif
      } else { // interfaced circumferential edge
        os << " (with kink) " << aol::detailedFormat ( outSl / arrSl );
#ifdef DEBUG
        if ( ( *pit != pn1 ) && ( *pit != pn2 ) ) {
          os << " but not radial";
        }
#endif
      }

      os << endl;

    }
  }
}

template< typename GridType >
void SlopeInterfaceBase<GridType>::dumpCoarseningWeights ( ostream &os ) const {
  os << "Dumping coarsening weights:" << endl;
  const int cWidth = _iMap.getGlobalIndex ( qc::CoordType ( 0, 1, 0 ) ) / 2 + 1;
  for ( qc::RectangularIterator<qc::QC_3D> it ( aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( cWidth, cWidth, cWidth ) ); it.notAtEnd(); ++it ) {
    for ( typename std::map < qc::CoordType, RealType >::const_iterator mit = _coarseningWeights.getRef ( *it ).begin(); mit != _coarseningWeights.getRef ( *it ).end(); ++mit ) {
      os << *it << " " << mit->first << ": " << aol::detailedFormat ( mit->second ) << endl;
    }
  }
}

template< typename GridType >
void SlopeInterfaceBase<GridType>::checkSumOfCoarseningWeights ( ostream &os ) const {
  std::map< qc::CoordType, RealType > sums;
  const int cWidth = _iMap.getGlobalIndex ( qc::CoordType ( 0, 1, 0 ) ) / 2 + 1;
  for ( qc::RectangularIterator<qc::QC_3D> it ( aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( cWidth, cWidth, cWidth ) ); it.notAtEnd(); ++it ) {
    for ( typename std::map < qc::CoordType, RealType >::const_iterator mit = _coarseningWeights.getRef ( *it ).begin(); mit != _coarseningWeights.getRef ( *it ).end(); ++mit ) {
      sums[ mit->first ] += mit->second;
#ifdef DEBUG
      if ( mit->second < -0.05 || mit->second > 1 )
        os << *it << " " << mit->first << ": restriction weight " << aol::detailedFormat ( mit->second ) << " out of range [-0.05, 1]" << endl;
#endif
    }
  }

  for ( typename std::map< qc::CoordType, RealType >::iterator sumspos = sums.begin(); sumspos != sums.end(); ++sumspos ) {
    if ( fabs ( sumspos->second - 1 ) > 1e-6 ) {
      os << sumspos->first << " " << aol::detailedFormat ( sumspos->second ) << endl;
    }
  }
}



template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_LIEHR > >::determineSlopes ( const GridType &Grid ) {
  // regular slopes for basis functions on non-interfaced edges
  for ( qc::RectangularIterator<qc::QC_3D> bit ( Grid ); bit.notAtEnd(); ++bit ) {
    const qc::CoordType pos = *bit;
    for ( short nbi = 0; nbi < NUM_NEIGHBORS; ++nbi ) {
      if ( nbi != NEIGHBOR_SELF ) { // node to node itself is not an edge
        const qc::CoordType neigh = pos + CFETopoLookup::_tetNbOffsets[ nbi ];
        if ( Grid.isAdmissibleNode ( neigh ) && !Grid.isEdgeInterfaced ( pos, neigh ) ) { // edge present (node is present and neighbor needs to be) and not interfaced

          this->setSlopes ( pos, pos, neigh, - aol::NumberTrait<RealType>::one, - aol::NumberTrait<RealType>::one );

        }
      }
    }
  }

  for ( typename GridType::VNMapTypeConstIt it = Grid.virtualNodeMapRef().begin(); it != Grid.virtualNodeMapRef().end(); ++it ) {
    aol::Vec2<int> node;
    it->second->splitIndex ( it->first, node[0], node[1] );

    const aol::Vector<int>& constraints = it->second->_constraints;
    const aol::RandomAccessContainer<RealType>& weights = it->second->_weights;
    const int numConstraints = it->second->numConstraints();

    aol::Vec2<int> ci;
    ci[0] = ci[1] = -1;

    // determine which two constraints are those corresponding to the edge on which the virtual node lies
    for ( int i = 0; i < 2; ++i ) {
      for ( int s = 0; s < numConstraints; ++s ) {
        if ( constraints[s] == node[i] ) {
          ci[i] = s;
        }
      }
    }

    if ( ci[0] == -1 || ci[1] == -1 ) throw aol::Exception ( "tpcfe::SlopeInterface: Constraint not found", __FILE__, __LINE__ );

    qc::CoordType pos[2];
    for ( short i = 0; i < 2; ++i ) {
      this->_iMap.splitGlobalIndex ( node[i], pos[i] );
    }

    const RealType CR = Grid.getCutRelationBetween ( pos[0], pos[1], 0 ); // only one interface
    const aol::Vec2<RealType> W ( weights[ci[0]], weights[ci[1]] );   // weight of virtual node wrt node[0], node[1]

    const RealType slope00 = - ( 1 - W[0] ) / (     CR );
    const RealType slope01 = - (     W[0] ) / ( 1 - CR );
    const RealType slope10 = - ( 1 - W[1] ) / ( 1 - CR );
    const RealType slope11 = - (     W[1] ) / (     CR );

    this->setSlopes ( pos[0], pos[0], pos[1], slope00, slope01 );
    this->setSlopes ( pos[1], pos[1], pos[0], slope10, slope11 );

  }
}

template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_LIEHR > >::determineCoarseningWeights ( const GridType &Grid ) {
  const int width = Grid.getNumXYZ();

  const int coarseWidth = width / 2 + 1;
  qc::OTFILexMapper< qc::QC_3D > coarseMapper ( coarseWidth );

  for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( coarseWidth, coarseWidth, coarseWidth ) ); cPosIt.notAtEnd(); ++cPosIt ) {
    const qc::CoordType cPos = *cPosIt;
    for ( int nbn = 0; nbn < NUM_NEIGHBORS; ++nbn ) {
      if ( nbn != NEIGHBOR_SELF ) {
        const qc::CoordType nbOffset = CFETopoLookup::_tetNbOffsets[nbn];  // offset to current neighbor
        const qc::CoordType cNPos  = cPos + nbOffset;                               // neighbor on coarse grid
        const qc::CoordType fPos   = static_cast<short> ( 2 ) * cPos;               // point itself on fine grid
        const qc::CoordType fNPos  = fPos + nbOffset;                               // intermediate node that is not present on coarse grid
        const qc::CoordType fNNPos = fPos + static_cast<short> ( 2 ) * nbOffset;    // next neighbor that is again present on coarse grid

        if ( CFEMG_is_inside ( cNPos, coarseWidth ) ) { // this is equivalent to all being inside

          // Attention: duplicate code!!! Cf. coarsening!
          const RealType restrictionWeight = ( this->getArrSlope ( fPos, fPos, fNPos ) ) / ( this->getOutSlope ( fNPos, fNPos, fNNPos ) + this->getOutSlope ( fNPos, fNPos, fPos ) );

          this->setCoarseningWeight ( fNPos, cPos, restrictionWeight );

        }
      }
    }
  }
}


template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_LIEHR > >::coarsenSlopeInterface ( SlopeInterfaceBase<GridType> &coarseSlI, const GridType &fineGrid, const GridType &coarseGrid ) {
  if ( fineGrid.getGridDepth() != coarseGrid.getGridDepth() + 1 )
    throw aol::Exception ( "tpcfe::SlopeInterface::coarsenSlopeInterface: incorrect grid depths", __FILE__, __LINE__ );

  for ( qc::RectangularIterator<qc::QC_3D> cPosIt ( coarseGrid ); cPosIt.notAtEnd(); ++cPosIt ) {
    const qc::CoordType cPos = *cPosIt;
    for ( int nbn = 0; nbn < NUM_NEIGHBORS; ++nbn ) {
      if ( nbn != NEIGHBOR_SELF ) {
        const qc::CoordType nbOffset = CFETopoLookup::_tetNbOffsets[nbn];  // offset to current neighbor
        const qc::CoordType cNPos  = cPos + nbOffset;                               // neighbor on coarse grid
        const qc::CoordType fPos   = static_cast<short> ( 2 ) * cPos;               // point itself on fine grid
        const qc::CoordType fNPos  = fPos + nbOffset;                               // intermediate node that is not present on coarse grid
        const qc::CoordType fNNPos = fPos + static_cast<short> ( 2 ) * nbOffset;    // next neighbor that is again present on coarse grid

        if ( coarseGrid.isAdmissibleNode ( cNPos ) ) { // this is equivalent to all being inside

          const RealType newCoarseningWeight = ( this->getArrSlope ( fPos, fPos, fNPos ) ) / ( this->getOutSlope ( fNPos, fNPos, fNNPos ) + this->getOutSlope ( fNPos, fNPos, fPos ) );
          const RealType newOutSlope = 2 * ( this->getOutSlope ( fPos, fPos, fNPos ) - newCoarseningWeight * this->getArrSlope ( fNPos, fNPos, fPos ) );
          const RealType newArrSlope = 2 * ( newCoarseningWeight * this->getArrSlope ( fNPos, fNPos, fNNPos ) );

          coarseSlI.setSlopes ( cPos, cPos, cNPos, newOutSlope, newArrSlope );

          this->setCoarseningWeight ( fNPos, cPos, newCoarseningWeight );

        }
      }
    }
  }
}



template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_TPOS > >::determineSlopes ( const GridType &Grid ) {

  // regular slopes for basis functions on non-interfaced edges
  for ( qc::RectangularIterator<qc::QC_3D> bit ( Grid ); bit.notAtEnd(); ++bit ) {
    const qc::CoordType pos = *bit;
    for ( short nbi = 0; nbi < NUM_NEIGHBORS; ++nbi ) {
      if ( nbi != NEIGHBOR_SELF ) { // node to node itself is not an edge
        const qc::CoordType neigh = pos + CFETopoLookup::_tetNbOffsets[ nbi ];
        if ( Grid.isAdmissibleNode ( neigh ) && !Grid.isEdgeInterfaced ( pos, neigh ) ) { // edge present (node is present and neighbor needs to be) and not interfaced
          this->storeRelevantEdge ( pos, NEIGHBOR_SELF, nbi );
          this->setSlopes ( pos, pos, neigh, - aol::NumberTrait<RealType>::one, - aol::NumberTrait<RealType>::one );
        }
      }
    }
  }


  // now the slopes for interfaced edges.
  for ( typename GridType::VNMapTypeConstIt vnIt = Grid.virtualNodeMapRef().begin(); vnIt != Grid.virtualNodeMapRef().end(); ++vnIt ) {
    aol::Vec2<int> edgeNodeIndex;
    vnIt->second->splitIndex ( vnIt->first, edgeNodeIndex[1], edgeNodeIndex[0] ); // this way, the smaller one is the first
    const int numConstraints = vnIt->second->numConstraints();
    const aol::Vector<int>& constraints = vnIt->second->_constraints;

    qc::CoordType edgepos[2];
    for ( short i = 0; i < 2; ++i ) {
      Grid.getIndexMapperRef().splitGlobalIndex ( edgeNodeIndex[i], edgepos[i] );
    }

    for ( int i = 0; i < numConstraints; ++i ) {
      qc::CoordType constrIPos;
      Grid.getIndexMapperRef().splitGlobalIndex ( constraints[i], constrIPos );

      const RealType bfctAtEN0 = ( edgeNodeIndex[0] == constraints[i] ? aol::NumberTrait<RealType>::one : aol::NumberTrait<RealType>::zero );
      const RealType bfctAtEN1 = ( edgeNodeIndex[1] == constraints[i] ? aol::NumberTrait<RealType>::one : aol::NumberTrait<RealType>::zero );
      const RealType CR        = Grid.getCutRelationBetween ( edgepos[0], edgepos[1], 0 );
      const RealType bfctAtVN  = vnIt->second->_weights[i];

      const RealType slopes0 = ( bfctAtVN  - bfctAtEN0 ) / (   CR   );
      const RealType slopes1 = ( bfctAtEN1 - bfctAtVN  ) / ( 1 - CR );

      if ( Grid.isEdgeInterfaced ( constrIPos, edgepos[0] ) ) {
        this->storeRelevantEdge ( constrIPos, edgepos[1], edgepos[0] ); // first and second must be on same side of interface
      } else {
        this->storeRelevantEdge ( constrIPos, edgepos[0], edgepos[1] );
      }

      this->setSlopes ( constrIPos, edgepos[0], edgepos[1], slopes0, slopes1 );

    }
  }
}


template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_TPOS > >::determineCoarseningWeights ( const GridType &Grid ) {

  if ( _useIncorrectStandardCoarsening ) {

    for ( qc::RectangularIterator<qc::QC_3D> posIt ( Grid ); posIt.notAtEnd(); ++posIt ) {
      if ( ( ( *posIt ) [0] % 2 == 0 ) && ( ( *posIt ) [1] % 2 == 0 ) && ( ( *posIt ) [2] % 2 == 0 ) ) { // we're only interested in those nodes present on the coarsened grid
        const qc::CoordType  fPos = *posIt, cPos = fPos / 2;
        for ( short i = 0; i < NUM_NEIGHBORS; ++i ) {
          const qc::CoordType NBPosf = fPos + CFELookup<RealType>::_tetNbOffsets[i];
          if ( Grid.isAdmissibleNode ( NBPosf ) && i != NEIGHBOR_SELF ) { // trivial coarsening weights are not stored
            this->setCoarseningWeight ( NBPosf, cPos, 0.5 );
          }
        }

        for ( unsigned int n = 0; n < _relevantEdges.getRef ( fPos ).size(); ++n ) {
          const qc::CoordType fNBPos = _relevantEdges.getRef ( fPos ) [ n ].getCOAneighbor ( fPos );
          if ( this->_coarseningWeights.getRef ( cPos ).find ( fNBPos ) == this->_coarseningWeights.getRef ( cPos ).end() ) {
            this->setCoarseningWeight ( fNBPos, cPos, 0.0 );
          }
        }
      }
    }

  } else {

    qc::OTFILexMapper< qc::QC_3D > coarseMapper ( Grid.getNumX() / 2 + 1, Grid.getNumY() / 2 + 1, Grid.getNumZ() / 2 + 1 );

    std::map < qc::CoordType, RealType > circumferentialCorrectionNumerators, circumferentialCorrectionDenominators;
    std::vector < std::pair < qc::CoordType, qc::CoordType > > circumferentialRestrictionWeights;

    for ( qc::RectangularIterator<qc::QC_3D> posIt ( Grid ); posIt.notAtEnd(); ++posIt ) {
      if ( ( ( *posIt ) [0] % 2 == 0 ) && ( ( *posIt ) [1] % 2 == 0 ) && ( ( *posIt ) [2] % 2 == 0 ) ) { // we're only interested in those nodes present on the coarsened grid
        const qc::CoordType  fPos = *posIt, cPos = fPos / 2;
        const int numRelEdges = _relevantEdges.getRef ( fPos ).size();

        // todo: find out wheter sparse or full matrix is faster here!
        aol::SparseMatrix<RealType> currentSystem ( numRelEdges, numRelEdges );
        aol::Vector<RealType> currentRHS ( numRelEdges ), currentSoln ( numRelEdges );

        for ( int m = 0; m < numRelEdges; ++m ) { // indexes equation, i. e. edge (through S neighbor)
          qc::CoordType E0, E1, E2;
          _relevantEdges.getRef ( fPos ) [ m ].getNAKedges ( fPos, E0, E1, E2 );
          for ( int n = 0; n < numRelEdges; ++n ) { // indexes (unknown) weights, i. e. coarsening neighbor considered
            const qc::CoordType fNBPos = _relevantEdges.getRef ( fPos ) [ n ].getCOAneighbor ( fPos );
            currentSystem.set ( m, n,  this->getArrSlope ( fNBPos, E0, E1 ) - this->getOutSlope ( fNBPos, E1, E2 ) );
          }
          currentRHS.set ( m, this->getOutSlope ( fPos, E1, E2 ) - this->getArrSlope ( fPos, E0, E1 ) );
        }

        aol::LUInverse<RealType> solver ( currentSystem );
        solver.apply ( currentRHS, currentSoln );

#ifdef DEBUG
        {
          aol::Vector<RealType> dummy ( currentRHS, aol::DEEP_COPY );
          dummy *= -1;
          currentSystem.applyAdd ( currentSoln, dummy );
          if ( dummy.norm() > ( 1.e-12 * currentSoln.norm() ) )
            cerr << "LU resid " << dummy.norm() << endl;
        }
#endif

        // now store restriction weights:
        for ( int n = 0; n < numRelEdges; ++n ) {
          const qc::CoordType fNBPos = _relevantEdges.getRef ( fPos ) [ n ].getCOAneighbor ( fPos );
          this->setCoarseningWeight ( fNBPos, cPos, currentSoln[n] );
          if ( this->_neighborOffsetToIndex.find ( fNBPos - fPos ) != this->_neighborOffsetToIndex.end() ) {
            // standard, i.e. radial neighbor
            circumferentialCorrectionNumerators[ fNBPos ] -= currentSoln[n];
          } else {
            // nonstandard, i.e. circumferential neighbor
            circumferentialCorrectionDenominators[ fNBPos ] += currentSoln[n];
            circumferentialRestrictionWeights.push_back ( std::pair< qc::CoordType, qc::CoordType > ( fNBPos, cPos ) );
          }
        }

      }
    }

    // numerators are 1 - ( sum of coarsening weights w.r.t. nodes forming coarse grid edge on which fine grid node lies.
    for ( typename std::map< qc::CoordType, RealType >::iterator it = circumferentialCorrectionNumerators.begin(); it != circumferentialCorrectionNumerators.end(); ++it ) {
      it->second += aol::NumberTrait<RealType>::one;
    }

    if ( _performNAKRelaxationAtCircumferentialEdges ) {
      // ==> first correction: allow for additional kinks along circumferential edges to enforce partition of unity condition
      for ( typename std::vector< std::pair< qc::CoordType, qc::CoordType > >::const_iterator it = circumferentialRestrictionWeights.begin(); it != circumferentialRestrictionWeights.end(); ++it ) {
        const qc::CoordType &fNPos = it->first;
        const qc::CoordType &cPos  = it->second;
        const RealType correctionFactor = circumferentialCorrectionNumerators[ fNPos ] / circumferentialCorrectionDenominators[ fNPos ]; // those two entries should exist by construction
        if ( fabs ( circumferentialCorrectionDenominators[ fNPos ] ) > 1e-10 /* some condition telling us we want to correct on circumferential edges only */ ) {
          this->_coarseningWeights.getRef ( cPos ) [ fNPos ] *= correctionFactor;
        } else {
#ifdef VERBOSE
          cerr << "Avoiding NAN (modifying circumferential restriction weights to preserve partition of unity) for " << it->first << " " << it->second << endl;
#endif
          // do nothing
        }
      }
    }

    if ( _performNAKRelaxationAtAllRemainingEdges ) {
      // finally, we may have the case that radial coarsening weights do not sum up to one (for a given fine grid node) but there is no circumferential weight via which to correct this.
      // ==> second correction: if the above fails, allow for additional kinks also on radial edges to enforce partition of unity.
      const int cWidth = coarseMapper ( qc::CoordType ( 0, 1, 0 ) );
      std::map< qc::CoordType, RealType > sums;
      for ( qc::RectangularIterator<qc::QC_3D> it ( aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( cWidth, cWidth, cWidth ) ); it.notAtEnd(); ++it ) {
        for ( typename std::map < qc::CoordType, RealType >::const_iterator mit = this->_coarseningWeights.getRef ( *it ).begin(); mit != this->_coarseningWeights.getRef ( *it ).end(); ++mit ) {
          sums[ mit->first ] += mit->second;
        }
      }

      for ( qc::RectangularIterator<qc::QC_3D> it ( aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( cWidth, cWidth, cWidth ) ); it.notAtEnd(); ++it ) {
        for ( typename std::map < qc::CoordType, RealType >::iterator mit = this->_coarseningWeights.getRef ( *it ).begin(); mit != this->_coarseningWeights.getRef ( *it ).end(); ++mit ) {
          const RealType correctionFactor = sums[ mit->first ];
          if ( fabs ( correctionFactor - 1 ) > 1.0e-6 ) {
            mit->second /= correctionFactor;
          }
        }
      }
    }

  }

  if ( _checkNAKQualityAfterCorrection != NULL ) {
    char filename[1024];
    sprintf ( filename, _checkNAKQualityAfterCorrection, Grid.getGridDepth() );
    ofstream nakQout ( filename );

    for ( qc::RectangularIterator<qc::QC_3D> posIt ( Grid ); posIt.notAtEnd(); ++posIt ) {
      if ( ( ( *posIt ) [0] % 2 == 0 ) && ( ( *posIt ) [1] % 2 == 0 ) && ( ( *posIt ) [2] % 2 == 0 ) ) { // we're only interested in those nodes present on the coarsened grid
        const qc::CoordType  fPos = *posIt, cPos = fPos / 2;
        const int numRelEdges = _relevantEdges.getRef ( fPos ).size();

        aol::Vector<RealType> arrSlopes ( numRelEdges ), outSlopes ( numRelEdges );
        std::vector< qc::CoordType > pts ( numRelEdges );

        for ( int m = 0; m < numRelEdges; ++m ) { // indexes position of kink
          qc::CoordType E0, E1, E2;
          _relevantEdges.getRef ( fPos ) [ m ].getNAKedges ( fPos, E0, E1, E2 );
          pts[m] = E1;

          for ( int n = 0; n < numRelEdges; ++n ) { // indexes fine grid neighbor considered for coarsening at fixed position
            const qc::CoordType fNBPos = _relevantEdges.getRef ( fPos ) [ n ].getCOAneighbor ( fPos );
            if ( this->_coarseningWeights.getRef ( cPos ).find ( fNBPos ) == this->_coarseningWeights.getRef ( cPos ).end() ) {
              throw ( aol::Exception ( "coarsening weight not found", __FILE__, __LINE__ ) );
            }
            const RealType coarseningWeight = this->_coarseningWeights.getRef ( cPos ).find ( fNBPos )->second;
            arrSlopes[m] += coarseningWeight * this->getArrSlope ( fNBPos, E0, E1 );
            outSlopes[m] += coarseningWeight * this->getOutSlope ( fNBPos, E1, E2 );
          }
          // and slopes of trivial neighbor, ie at bfct itself
          arrSlopes[m] += 1.0 * this->getArrSlope ( fPos, E0, E1 );
          outSlopes[m] += 1.0 * this->getOutSlope ( fPos, E1, E2 );
        }

        for ( int m = 0; m < numRelEdges; ++m ) {
          if ( fabs ( arrSlopes[m] / outSlopes[m] - 1.0 ) > 1.0e-10 )
            nakQout << arrSlopes[m] << " " << outSlopes[m] << " " << arrSlopes[m] / outSlopes[m] << " " << pts[m] << endl;
        }
      }
    }
  }

}


template< typename RealType >
void SlopeInterface< CFEGrid< RealType, CFE_TPOS > >::coarsenSlopeInterface ( SlopeInterfaceBase<GridType> &coarseSlIntf, const GridType &fineGrid, const GridType &coarseGrid ) {

  determineCoarseningWeights ( fineGrid );

  SlopeInterface<GridType> &coarseSlI = dynamic_cast< SlopeInterface<GridType>& > ( coarseSlIntf );

  for ( qc::RectangularIterator<qc::QC_3D> posIt ( coarseGrid ); posIt.notAtEnd(); ++posIt ) {
    const qc::CoordType cPos = *posIt, fPos = static_cast<short> ( 2 ) * cPos;

    const int numRelEdges = _relevantEdges.getRef ( fPos ).size();

    for ( int n = 0; n < numRelEdges; ++n ) {
      const qc::CoordType offset = _relevantEdges.getRef ( fPos ) [ n ].getOffset();
      const qc::CoordType direction = _relevantEdges.getRef ( fPos ) [ n ].getDirection();
      const qc::CoordType fInit = fPos  + static_cast<short> ( 2 ) * offset;
      const qc::CoordType fCtr  = fInit + direction;
      const qc::CoordType fTerm = fCtr  + direction;
      const qc::CoordType cInit = cPos  + offset;
      const qc::CoordType cTerm = cInit + direction;

      const short offsI = _relevantEdges.getRef ( fPos ) [ n ].getOffsetI(), direI = _relevantEdges.getRef ( fPos ) [ n ].getDirectionI();

      // collect terms for new slopes: start with fPos as trivial neighbor, then consider all relevant neighbors.
      RealType newOutSlope = this->getOutSlope ( fPos, fInit, fCtr ), newArrSlope = this->getArrSlope ( fPos, fCtr, fTerm ); // second should be zero.

      for ( int nn = 0; nn < numRelEdges; ++nn ) {
        const qc::CoordType fNPos = _relevantEdges.getRef ( fPos ) [ nn ].getCOAneighbor ( fPos );
        const RealType coaWeight = this->getCoarseningWeight ( fNPos, cPos );

        newOutSlope += coaWeight * this->getOutSlope ( fNPos, fInit, fCtr );
        newArrSlope += coaWeight * this->getArrSlope ( fNPos, fCtr, fTerm );
      }

      newOutSlope *= 2; newArrSlope *= 2; // rescale because slopes are relative to current grid spacing

      if ( ( aol::Sqr ( newOutSlope ) + aol::Sqr ( newArrSlope ) ) == aol::NumberTrait<RealType>::zero ) {
#ifdef VERBOSE
        // todo: think about whether this requires additional treatment
        cerr << "Fine level " << fineGrid.getGridDepth() << ", on edge " << fInit << " " << fCtr << " " << fTerm << "for " << cPos << ": zero slopes; " << numRelEdges << endl;
#endif
      }

      // store data:
      coarseSlI.storeRelevantEdge ( cPos, offsI, direI );
      coarseSlI.setSlopes ( cPos, cInit, cTerm, newOutSlope, newArrSlope );

    }
  }
}


template< typename GridType >
void CFEProlongOp<GridType>::apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
  if ( _opType == aol::ASSEMBLED ) {

    _pmat->apply ( Arg, Dest );

  } else {

    if ( ( Arg.size() != _coarseGrid.getNumberOfNodes() ) || Dest.size() != _fineGrid.getNumberOfNodes() ) {
      throw aol::Exception ( "Cannot apply CFE Prolongation operator: incorrect vector size.", __FILE__, __LINE__ );
    }

    Dest.setZero();

    switch ( _coarseGrid.CT ) {
      case CFE_CD : // intended jump across labels
      case CFE_LIEHR:
      case CFE_TPOSELAST:
        for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
          if ( _coarseGrid.isDOFNode ( *bit ) ) { // coarse is arg
            for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
              const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
              if ( _fineGrid.isAdmissibleNode ( NBPosf ) ) {
                const int fine_co   = _fineGrid.getIndexMapperRef().getGlobalIndex ( NBPosf );
                const int coarse_co = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
                if ( _fineGrid.isDOFNode ( fine_co ) ) { // fine is dest
                  Dest.add ( fine_co, ( ( i == NEIGHBOR_SELF ? 1.0 : _mgWeightProvider.getProperProlongationWeight ( NBPosf, *bit ) ) *  Arg.get ( coarse_co ) ) );
                }
              }
            }
          }
        }
        break;

      case CFE_TPOS:
        for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
          if ( _coarseGrid.isDOFNode ( *bit ) ) {
            for ( typename std::map< qc::CoordType, RealType >::const_iterator mit = _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).begin(); mit != _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).end(); ++mit ) {
              if ( _fineGrid.isDOFNode ( mit->first ) ) {
                const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
                const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( mit->first );
                Dest.add ( fineInd, mit->second * Arg.get ( coarseInd ) );
              }
            }

            // the following ones are not stored in the map
            const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
            const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( *bit + *bit );
            if ( _fineGrid.isDOFNode ( fineInd ) ) {
              Dest.add ( fineInd, /* 1 * */ Arg.get ( coarseInd ) );
            }
          }
        }
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFEProlongOp::apply is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

    }
  }
}


template< typename GridType >
void CFEProlongOp<GridType>::assembleAddMatrix ( aol::SparseMatrix<RealType> &mat ) const {
  switch ( _coarseGrid.CT ) {
    case CFE_CD : // intended jump across labels
    case CFE_LIEHR:
    case CFE_TPOSELAST:
      for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
        if ( _coarseGrid.isDOFNode ( *bit ) ) { // coarse is arg
          for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
            const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
            if ( _fineGrid.isAdmissibleNode ( NBPosf ) ) {
              const int coarse_co = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
              const int fine_co   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( NBPosf );

              if ( _fineGrid.isDOFNode ( fine_co ) ) {
                mat.add ( fine_co, coarse_co,  ( i == NEIGHBOR_SELF ? aol::NumberTrait<RealType>::one : _mgWeightProvider.getProperProlongationWeight ( NBPosf, *bit ) ) );
              }
            }
          }
        }
      }
      break;

    case CFE_TPOS:
      for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
        if ( _coarseGrid.isDOFNode ( *bit ) ) {
          for ( typename std::map< qc::CoordType, RealType >::const_iterator mit = _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).begin(); mit != _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).end(); ++mit ) {
            if ( _fineGrid.isDOFNode ( mit->first ) ) {
              const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
              const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( mit->first );
              mat.add ( fineInd, coarseInd, mit->second );
            }
          }

          // the following ones are not stored in the map
          const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
          const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( *bit + *bit );
          if ( _fineGrid.isDOFNode ( fineInd ) ) {
            mat.add ( fineInd, coarseInd, aol::NumberTrait<RealType>::one );
          }
        }
      }
      break;


    default:
      throw aol::UnimplementedCodeException ( "CFEProlongOp::assembleAddMatrix is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

  }
}


template< typename GridType >
void CFERestrictOp<GridType>::apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const  {
  if ( _opType == aol::ASSEMBLED ) {

    _pmat->apply ( Arg, Dest );

  } else {

    if ( ( Arg.size() != _fineGrid.getNumberOfNodes() ) || Dest.size() != _coarseGrid.getNumberOfNodes() ) {
      throw aol::Exception ( "Cannot apply CFE Restriction operator: incorrect vector size.", __FILE__, __LINE__ );
    }

    Dest.setZero();

    switch ( _coarseGrid.CT ) {
      case CFE_CD: // intended jump across labels
      case CFE_LIEHR:
      case CFE_TPOSELAST:
        for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
          if ( _coarseGrid.isDOFNode ( *bit ) ) { // coarse is dest
            RealType res = aol::NumberTrait<RealType>::zero;
            const int coarse_co = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
            for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
              const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
              if ( _fineGrid.isAdmissibleNode ( NBPosf ) ) {
                const int fine_co   = _fineGrid.getIndexMapperRef().getGlobalIndex ( NBPosf );

                res += ( _fineGrid.isDOFNode ( fine_co ) *    // if use fine value, fine is arg
                         ( i == NEIGHBOR_SELF ? 1.0 : _mgWeightProvider.getProperRestrictionWeight ( NBPosf, *bit ) ) *
                         Arg.get ( fine_co )               ); // corresponding arg value

              }
            }
            Dest.add ( coarse_co, res );
          }
        }
        break;

      case CFE_TPOS:
        for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
          if ( _coarseGrid.isDOFNode ( *bit ) ) {
            RealType res = aol::NumberTrait<RealType>::zero;
            for ( typename std::map< qc::CoordType, RealType >::const_iterator mit = _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).begin(); mit != _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).end(); ++mit ) {
              const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( mit->first );
              res += _fineGrid.isDOFNode ( mit->first ) * mit->second * Arg.get ( fineInd );
            }
            // the following ones are not stored:
            const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
            const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( *bit + *bit );
            res += _fineGrid.isDOFNode ( fineInd ) * Arg.get ( fineInd );
            Dest.add ( coarseInd, res );
          }
        }
        break;

      default:
        throw aol::UnimplementedCodeException ( "CFERestrictOp::apply is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

    }

  }

}



template< typename GridType >
void CFERestrictOp<GridType>::assembleAddMatrix ( aol::SparseMatrix<RealType> &mat ) const {

  switch ( _coarseGrid.CT ) {
    case CFE_CD: // intended jump across labels
    case CFE_LIEHR:
    case CFE_TPOSELAST:
      for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
        if ( _coarseGrid.isDOFNode ( *bit ) ) { // coarse is dest
          for ( int i = 0; i < NUM_NEIGHBORS; ++i ) {
            const qc::CoordType NBPosf = static_cast<short> ( 2 ) * ( *bit ) + CFELookup<RealType>::_tetNbOffsets[i];
            if ( _fineGrid.isAdmissibleNode ( NBPosf ) ) {
              const int coarse_co = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
              const int fine_co   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( NBPosf );

              if ( _fineGrid.isDOFNode ( fine_co ) ) {
                mat.add ( coarse_co, fine_co,  ( i == NEIGHBOR_SELF ? 1.0 : _mgWeightProvider.getProperRestrictionWeight ( NBPosf, *bit ) ) );
              }
            }
          }
        }
      }
      break;

    case CFE_TPOS:
      for ( qc::RectangularIterator<qc::QC_3D> bit ( _coarseGrid ); bit.notAtEnd(); ++bit ) {
        if ( _coarseGrid.isDOFNode ( *bit ) ) {
          for ( typename std::map< qc::CoordType, RealType >::const_iterator mit = _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).begin(); mit != _mgWeightProvider.getCoarseningWeightMapRef ( *bit ).end(); ++mit ) {
            const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
            const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( mit->first );
            if ( _fineGrid.isDOFNode ( fineInd ) ) {
              mat.add ( coarseInd, fineInd, mit->second );
            }
          }
          // the following ones are not stored:
          const int coarseInd = _coarseGrid.getIndexMapperRef().getGlobalIndex ( *bit );
          const int fineInd   =   _fineGrid.getIndexMapperRef().getGlobalIndex ( *bit + *bit );
          if ( _fineGrid.isDOFNode ( fineInd ) ) {
            mat.add ( coarseInd, fineInd, aol::NumberTrait<RealType>::one );
          }
        }
      }
      break;

    default:
      throw aol::UnimplementedCodeException ( "CFERestrictOp::assembleAddMatrix is not implemented yet for this ConstraintType.", __FILE__, __LINE__ );

  }
}


#define INSTANTIATE_CFEClassG2(ClassName, CT)\
template class ClassName < CFEGrid < float,       CT > >;\
template class ClassName < CFEGrid < double,      CT > >;\
template class ClassName < CFEGrid < long double, CT > >;

#define INSTANTIATE_CFEClassElastG32Arg(ClassName, CT, ElastCoeffT)\
template class ClassName < CFEGrid < float,       CT, ElastCoeffT<float>       > >;\
template class ClassName < CFEGrid < double,      CT, ElastCoeffT<double>      > >;\
template class ClassName < CFEGrid < long double, CT, ElastCoeffT<long double> > >;

#define INSTANTIATE_CFEClassElastG3(ClassName, CT)\
INSTANTIATE_CFEClassElastG32Arg(ClassName, CT, tpcfe::IsotropicElasticityCoefficient)\
INSTANTIATE_CFEClassElastG32Arg(ClassName, CT, tpcfe::VoigtElasticityCoefficient)

INSTANTIATE_CFEClassG2 ( SlopeInterfaceBase, CFE_CD )
INSTANTIATE_CFEClassG2 ( SlopeInterfaceBase, CFE_LIEHR )
INSTANTIATE_CFEClassG2 ( SlopeInterfaceBase, CFE_TPOS )

INSTANTIATE_CFEClassElastG3 ( SlopeInterfaceBase, CFE_CD )
INSTANTIATE_CFEClassElastG3 ( SlopeInterfaceBase, CFE_TPOSELAST )


template<> bool SlopeInterface < CFEGrid<       float, CFE_TPOS > >::_performNAKRelaxationAtCircumferentialEdges = false;
template<> bool SlopeInterface < CFEGrid<      double, CFE_TPOS > >::_performNAKRelaxationAtCircumferentialEdges = false;
template<> bool SlopeInterface < CFEGrid< long double, CFE_TPOS > >::_performNAKRelaxationAtCircumferentialEdges = false;

template<> bool SlopeInterface < CFEGrid<       float, CFE_TPOS > >::_performNAKRelaxationAtAllRemainingEdges = true;
template<> bool SlopeInterface < CFEGrid<      double, CFE_TPOS > >::_performNAKRelaxationAtAllRemainingEdges = true;
template<> bool SlopeInterface < CFEGrid< long double, CFE_TPOS > >::_performNAKRelaxationAtAllRemainingEdges = true;

template<> const char* SlopeInterface < CFEGrid<       float, CFE_TPOS > >::_checkNAKQualityAfterCorrection = NULL;
template<> const char* SlopeInterface < CFEGrid<      double, CFE_TPOS > >::_checkNAKQualityAfterCorrection = NULL;
template<> const char* SlopeInterface < CFEGrid< long double, CFE_TPOS > >::_checkNAKQualityAfterCorrection = NULL;

template<> bool SlopeInterface < CFEGrid<       float, CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
template<> bool SlopeInterface < CFEGrid<      double, CFE_TPOS > >::_useIncorrectStandardCoarsening = false;
template<> bool SlopeInterface < CFEGrid< long double, CFE_TPOS > >::_useIncorrectStandardCoarsening = false;


INSTANTIATE_CFEClassG2 ( SlopeInterface, CFE_LIEHR )
INSTANTIATE_CFEClassG2 ( SlopeInterface, CFE_TPOS )

INSTANTIATE_CFEClassG2 ( CFEProlongOp, CFE_CD )
INSTANTIATE_CFEClassG2 ( CFEProlongOp, CFE_LIEHR )
INSTANTIATE_CFEClassG2 ( CFEProlongOp, CFE_TPOS )

INSTANTIATE_CFEClassElastG3 ( CFEProlongOp, CFE_CD )
INSTANTIATE_CFEClassElastG3 ( CFEProlongOp, CFE_TPOSELAST )

INSTANTIATE_CFEClassG2 ( CFERestrictOp, CFE_CD )
INSTANTIATE_CFEClassG2 ( CFERestrictOp, CFE_LIEHR )
INSTANTIATE_CFEClassG2 ( CFERestrictOp, CFE_TPOS )

INSTANTIATE_CFEClassElastG3 ( CFERestrictOp, CFE_CD )
INSTANTIATE_CFEClassElastG3 ( CFERestrictOp, CFE_TPOSELAST )

#undef INSTANTIATE_CFEClassG2
#undef INSTANTIATE_CFEClassElastG32Arg
#undef INSTANTIATE_CFEClassElastG3

}
