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

#include <tpCFEGrid.h>
#include <progressBar.h>

namespace tpcfe {

template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
CFEGridBase < _RealType, _CT, NodalCoeffType >::CFEGridBase ( const int Depth ) :
    CFEGridDefinition ( Depth ),
    _inadmissibleSignatureWarningPrinted ( false ),
    _structure ( 0 ),
    _deleteStructures ( MAX_STRUCT_ID + 1 ),
    _virtualNodesInitialized ( false ),
    _verboseMode ( false ),
    _pDomain ( NULL ),
    _pDomainNodeMask ( NULL ), _pDirichletMask ( NULL ), _pDOFMask ( NULL ),
    _adjustLevelset ( 1.0e-6 ), _adjustVirtualNodes ( 0.0 ),
    _isCoarsened ( false ) {

  int mem = CFELookup<RealType>::construct();
  if ( _verboseMode ) {
    cerr << "Lookup table is of size " << mem / 1024 << "KiB " << endl;
  }

  _deleteStructures.setAll ( false );
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
CFEGridBase < _RealType, _CT, NodalCoeffType >::CFEGridBase ( const qc::GridSize<qc::QC_3D> &Sizes ) :
    CFEGridDefinition ( Sizes ),
    _inadmissibleSignatureWarningPrinted ( false ),
    _structure ( 0 ),
    _deleteStructures ( MAX_STRUCT_ID + 1 ),
    _virtualNodesInitialized ( false ),
    _verboseMode ( false ),
    _pDomain ( NULL ),
    _pDomainNodeMask ( NULL ), _pDirichletMask ( NULL ), _pDOFMask ( NULL ),
    _adjustLevelset ( 1.0e-6 ), _adjustVirtualNodes ( 0.0 ),
    _isCoarsened ( false ) {

  int mem = CFELookup<RealType>::construct();
  if ( _verboseMode ) {
    cerr << "Lookup table is of size " << mem / 1024 << "KiB " << endl;
  }

  _deleteStructures.setAll ( false );
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
CFEGridBase < _RealType, _CT, NodalCoeffType >::~CFEGridBase ( ) {
  removeVirtualNodes();

  for ( unsigned int i = 0; i < _structure.size(); ++i ) {
    if ( _deleteStructures[i] ) {
      delete ( _structure[i] );
    }
  }

  if ( _deleteStructures[ MAX_STRUCT_ID ] && _pDomain )  delete _pDomain;
  if ( _pDomainNodeMask )  delete _pDomainNodeMask;
  if ( _pDirichletMask )   delete _pDirichletMask;
  if ( _pDOFMask )         delete _pDOFMask;

  CFELookup<RealType>::destruct();
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
_RealType CFEGridBase < _RealType, _CT, NodalCoeffType >::relaxedDetectAndInitVirtualNodes ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeff, const RealType initThreshold, const RealType thresholdStep ) {
  VNType::_maximalInversionError = initThreshold;
  bool initialized = false;
  while ( !initialized ) {
    try {
      detectAndInitVirtualNodes ( nodalCoeff );
      initialized = true;
    } catch ( aol::Exception &ex ) {
      ex.dump();
      removeVirtualNodes();
      VNType::_maximalInversionError += thresholdStep;
      cerr << "retrying with VNType::_maximalInversionError = " << VNType::_maximalInversionError << endl;
    }
  }
  return ( VNType::_maximalInversionError );
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
bool CFEGridBase < _RealType, _CT, NodalCoeffType >::isInterfaced ( const CFEElement < RealType > &e ) const {
  const int index = e.globalIndex ();
  return ( _elType[index].representsInterfaced() );
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase < _RealType, _CT, NodalCoeffType >::setElTypes ( int structureNo ) {
  CFEStructure<RealType> *s = ( structureNo == MAX_STRUCT_ID ? _pDomain : _structure[structureNo] );

  const qc::GridSize<qc::QC_3D> gridSize ( *this );
  for ( FullElementIterator it ( *this ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> curEl ( *it, gridSize, this->getElType ( *it ) );

    CFEType tempType;
    unsigned char  bitMask = 1;
    const int globalIndex = curEl.globalIndex();

    for ( int i = 0; i < 8; ++i ) {                   // Run through all nodes of element
      const qc::CoordType epos = *it + CFETopoLookup::_hexNbOffsets[i];
      const RealType p = s->getValue ( epos );
      if ( p < 0 ) tempType._pureType |= bitMask;
      bitMask <<= 1;
    }

    if ( tempType.representsInterfaced() ) { // this is necessary because an element can be non-interfaced by more than one structure
      tempType._structureNo = structureNo;
    }

    // An element is interfaced by more than one structure if the element is interfaced by another structure
    // and the element is interfaced by the new structure.
    if ( _elType[globalIndex].representsInterfaced() && tempType.representsInterfaced() ) {
      throw aol::Exception ( "CFEGridBase::setElTypes: This element is interfaced by more than one structure!", __FILE__, __LINE__ );
    } else {
      //     } else if ( CFETopoLookup::_admissibleSignature [ tempType._pureType ] == false ) {
      if ( ( CFETopoLookup::_admissibleSignature [ tempType._pureType ] == false ) && ( _inadmissibleSignatureWarningPrinted == false ) ) {
        cerr << "CFEGridBase::setElTypes: Element " << *it << " is interfaced more than once by the interface. Signature is " << static_cast<short>( tempType._pureType ) << endl;
        cerr << "Currently not throwing exception. This warning will not be printed again." << endl;
        _inadmissibleSignatureWarningPrinted = true;
      // throw aol::Exception ( "CFEGridBase::setElTypes: This element is interfaced more than once by the interface!", __FILE__, __LINE__ );
      }

      // The type of the element needs an update if the element is interfaced by the structure or the structure is a domain.
      // (later the type is used to specify whether an element is in a domain or not)
      if ( ( tempType._pureType != 0 ) || structureNo == MAX_STRUCT_ID ) {
        _elType[ globalIndex ] = tempType;
      }
    }

  }

}

template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::detectVirtualNodes() {

  aol::ProgressBar<> pb ( "Detecting virtual nodes" );
  pb.start ( this->getNumberOfElements() );

  // Check all elements
  const qc::GridSize<qc::QC_3D> gridSize ( *this );
  for ( FullElementIterator it ( *this ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, this->getElType ( *it ) );
    const CFEType type = el.cfeType();

    if ( _verboseMode ) pb++;

    if ( type.representsInterfaced() ) {

      el.computeAssembleData ( *this );

      // Run through tetrahedra and get virtual nodes
      for ( CFETopoTetraIterator it ( type ); it.notAtEnd(); ++it ) {
        for ( short i = 0; i < 4; ++i ) { // Check whether any node of this tetra is virtual
          const short li1 = (*it) ( i, 1 );
          if ( li1 != NODE_NOT_VIRTUAL ) { // This is a virtual node
            // Create virtual node object
            const short li0 = (*it) ( i, 0 );
            const VNIndexType mapIndex = VNType::mapIndex ( el.globalIndex ( li0 ), el.globalIndex ( li1 ) );

            VNMapTypeConstIt vnIt = _vNode.find ( mapIndex );
            if ( vnIt == _vNode.end() ) { // virtual node is new
              VNType *vnode = new VNType ( li0, li1,
                                           el.globalIndex ( li0 ),
                                           el.globalIndex ( li1 ),
                                           el, it->getParent() );
              _vNode[mapIndex] = vnode;

            } else { // virtual node has been found before
              vnIt->second->addInElement ( el, it->getParent() );
            }
          }
        }
      }
    }
  }

  if ( _verboseMode ) {
    pb.finish();
    cerr << "Detected " << _vNode.size() << " virtual nodes " << endl;
  }

}

template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::initVirtualNodes ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeff ) {

  aol::ProgressBar<> pb ( "Computing weights" );
  pb.start ( _vNode.size() );

  uint64_t mem = 0;
  for ( VNMapTypeConstIt vnIt = _vNode.begin(); vnIt != _vNode.end(); ++vnIt ) {
    VNType                       &vn = * ( vnIt->second );
    const CFETopoElement         &el = vn.getInElement();
    const CFEStructure<RealType> &s = getStructureRefForType ( el.cfeType() );

    if ( CT != CFE_DOMAIN ) vn.determineApproximateNormal ( s );
#ifdef VERBOSE
    int gIdx0, gIdx1;
    vn.splitIndex ( vnIt->first, gIdx0, gIdx1 );
    cerr << "determining constraints for " << vnIt->first << ": " << gIdx0 << " " << gIdx1 << " in element " << el << endl;
#endif
    vn.determineConstraints ( nodalCoeff, *this );

    mem += vn.memoryUsed();
    mem += sizeof ( VNIndexType );
    mem += sizeof ( VNType* );

    if ( _verboseMode ) pb++;
  }

  if ( _verboseMode ) {
    pb.finish();
    cerr << "   Memory used is " << mem << " bytes = " << ( mem*1.0 / 1048576 ) << " MiB " << endl;
    if ( _vNode.size() ) {
      cerr << "   That is " << mem / _vNode.size() << " bytes = " << ( mem*1.0 / ( 1024*_vNode.size() ) ) << " KiB per node" << endl;
    }
  }
  _virtualNodesInitialized = true;
}

template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::setDomainRef ( tpcfe::CFEStructure<RealType>& Domain, const bool deleteDomain ) {
  if ( CT != CFE_CD && CT != CFE_CDWI_TPOS && CT != CFE_CDWI_LIEHR && CT != CFE_CDWI_ELAST )
    throw aol::Exception ( "CFEGridBase::setDomainRef called for illegal constraint type", __FILE__, __LINE__ );

  this->removeDomain();

  _pDomain = &Domain;
  _deleteStructures.set ( MAX_STRUCT_ID, deleteDomain );

  setElTypes ( MAX_STRUCT_ID );

  if ( _pDomainNodeMask )
    delete ( _pDomainNodeMask );

  _pDomainNodeMask = new BitArrayType ( *this );
  _pDomainNodeMask->setAll ( false );

  for ( qc::RectangularIterator<qc::QC_3D> rit ( *this ); rit.notAtEnd(); ++rit ) {
    if ( getDomainRef().getValue ( *rit ) < 0 ) {
      for ( int l = 0; l < NUM_NEIGHBORS; ++l ) {
        const qc::CoordType nbpos = ( *rit ) + CFETopoLookup::_tetNbOffsets[l];
        if ( this->isAdmissibleNode ( nbpos ) ) {
          _pDomainNodeMask->set ( nbpos, true );
        }
      }
    }
  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::doRestrictToDomain ( aol::Vector<RealType> &arg ) const {
  if ( arg.size() != this->getNumberOfNodes() )
    throw aol::Exception ( "CFEGridBase::restrictToDomain: size mismatch", __FILE__, __LINE__ );

  if ( CT != CFE_CD && CT != CFE_CDWI_TPOS && CT != CFE_CDWI_LIEHR )
    throw aol::Exception ( "CFEGridBase::restrictToDomain: illegal constraint type", __FILE__, __LINE__ );

  if ( _pDomainNodeMask == NULL ) {
    // do nothing
  } else {
    for ( int i = 0; i < arg.size(); ++i )
      if ( !isDomainNode ( i ) )
        arg[i] = aol::NumberTrait<RealType>::zero;

  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::restrictDirichletNodes ( aol::Vector<_RealType> &arg ) const {
  if ( ( CT == CFE_CD || CT == CFE_TPOS || CT == CFE_TPOSELAST || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) && ( _pDirichletMask != NULL ) ) {
    for ( int i = 0; i < arg.size(); ++i )
      if ( isDirichletNode ( i ) )
        arg[i] = aol::NumberTrait<RealType>::zero;

  } else {
    throw aol::Exception ( "CFEGridBase::restrictDirichletNodes not implemented for CT != CFE_CD, CFE_TPOS, CFE_TPOSELAST, CFE_CDWI_TPOS, CFE_CDWI_LIEHR, CFE_CDWI_ELAST or domain not set.", __FILE__, __LINE__ );
  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::restrictToDofs ( aol::Vector<_RealType> &arg ) const {
  if ( CT == CFE_CD || CT == CFE_TPOS || CT == CFE_TPOSELAST || CT == CFE_LIEHR || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
    if ( _pDOFMask == NULL ) {
      // do nothing
    } else {
      for ( int i = 0; i < arg.size(); ++i )
        if ( !isDOFNode ( i ) )
          arg[i] = aol::NumberTrait<RealType>::zero;
    }
  } else {
    throw aol::Exception ( "CFEGridBase::restrictToDOFs not implemented for CT != CFE_CD, CFE_TPOS, CFE_TPOSELAST, CFE_LIEHR, CFE_CDWI_TPOS, CFE_CDWI_LIEHR, CFE_CDWI_ELAST", __FILE__, __LINE__ );
  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::restrictToInnerNodes ( qc::ScalarArray<_RealType, qc::QC_3D> &arg ) const {
  if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR ) {
    if ( _pDomainNodeMask == NULL ) {
      // do nothing
    } else {
      for ( int i = 0; i < arg.size(); ++i ) {
        if ( !isInnerNode ( i ) ) {
          arg.set ( i, aol::NumberTrait<RealType>::zero );
        }
      }

    }
  } else {
    throw aol::Exception ( "CFEGridBase::restrictToInnerNodes not implemented for CT != CFE_CD, CFE_CDWI_TPOS, CFE_CDWI_LIEHR", __FILE__, __LINE__ );
  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
void CFEGridBase< _RealType, _CT, NodalCoeffType >::getCutRelations ( CFEElement<_RealType> &el ) const {
  for ( int edg = 0; edg < tpcfe::NUM_TETRA_EDGES_IN_CUBE; ++edg ) {
    const int i = CFETopoLookup::_tetraEdges[ edg ][0], j = CFETopoLookup::_tetraEdges[ edg ][1];
    const qc::CoordType posi = el + CFETopoLookup::_hexNbOffsets[ i ], posj = el + CFETopoLookup::_hexNbOffsets[ j ];

    RealType CR = this->getCutRelationBetween ( posi, posj, el.cfeType()._structureNo );

    // adjustment:
    const RealType adjustVN = this->getAdjustVirtualNodes();
    if ( CR >= 0.0 && CR < 0.0 + adjustVN ) {
#ifdef VERBOSE
      cerr << "Adjusting virtual node " << posi << " " << posj << ", CR = " << CR << " adjusted to " << 0.0 + adjustVN << endl;
#endif
      CR = 0.0 + adjustVN;
    }
    if ( CR <= 1.0 && CR > 1.0 - adjustVN ) {
#ifdef VERBOSE
      cerr << "Adjusting virtual node " << posi << " " << posj << ", CR = " << CR << " adjusted to " << 1.0 - adjustVN << endl;
#endif
      CR = 1.0 - adjustVN;
    }

    el.setCutRelation ( i, j, CR );
    el.setCutRelation ( j, i, aol::NumberTrait<RealType>::one - CR );
  }
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
bool CFEGridBase< _RealType, _CT, NodalCoeffType >::getStructureValues ( CFEElement<_RealType> &el ) const {
  const int structureNo = el.cfeType()._structureNo;

  // Now fill the given array with values
  if ( ( structureNo >= 0 ) && ( structureNo < static_cast<int> ( _structure.size() ) ) ) {
    CFEStructure<RealType> &s = * ( _structure[structureNo] );
    getStructureValues ( s, el );
  } else if ( structureNo == MAX_STRUCT_ID ) {
    CFEStructure<RealType> &s = *_pDomain;
    getStructureValues ( s, el );
  } else {
    cerr << "structure " << structureNo << " does not exist." << endl;
    return false;
  }
  return true;
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
_RealType CFEGridBase< _RealType, _CT, NodalCoeffType >::getTotalPartialVolume ( const signed char which ) const {
  if ( ( which != 1 ) && ( which != -1 ) )
    throw aol::Exception ( "CFEGridBase::getTotalPartialVolume only works for inner or outer volume", __FILE__, __LINE__ );

  RealType totalVolume = 0;

  const qc::GridSize<qc::QC_3D> gridSize ( *this );
  for ( FullElementIterator it ( *this ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> curEl ( *it, gridSize, this->getElType ( *it ) );
    const CFEType cfeType = curEl.cfeType();

    if ( ( ( cfeType._pureType == 0x00 ) && ( which == +1 ) ) ||  ( ( cfeType._pureType == 0xFF ) && ( which == -1 ) ) ) {
      totalVolume += 1;
    } else if ( ( ( cfeType._pureType == 0x00 ) && ( which == -1 ) ) ||  ( ( cfeType._pureType == 0xFF ) && ( which == +1 ) ) ) {
      // no volume to add
    } else {
      curEl.computeAssembleData ( *this );

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( curEl, which ); tit.notAtEnd(); ++tit ) {
        totalVolume += tit->getVolume();
      }
    }

  }

  totalVolume *= aol::Cub ( this->H() );

  return ( totalVolume );
}


template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
bool CFEGridBase< _RealType, _CT, NodalCoeffType >::interfaced ( const CFEElement<_RealType> &el, int parentTetra ) const {
  int numPlusNodes = 0;
  int numMinusNodes = 0;

  const CFEStructure<RealType>& cfeInterface = getStructureRef ( 0 );

  // Compute the vertices relative to the virtual node
  for ( int i = 0; i < 4; ++i ) {
    const int vertexIndex = CFELookup<RealType>::_stdTetraVertex[parentTetra][i];
    const qc::CoordType epos = el + CFETopoLookup::_hexNbOffsets[ vertexIndex ];

    if ( cfeInterface.getValue ( epos ) > 0 ) {
      ++numPlusNodes;
    } else {
      ++numMinusNodes;
    }
  }
  return ( numMinusNodes != 0 && numPlusNodes != 0 );
}


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

INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_TP )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_TPOS )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_LIEHR )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_DOMAIN )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_CD )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_CDWI_TPOS )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_CDWI_LIEHR )
INSTANTIATE_CFEClass3 ( CFEGridBase, CFE_NONE )

INSTANTIATE_CFEClassElast3 ( CFEGridBase, CFE_CD )
INSTANTIATE_CFEClassElast3 ( CFEGridBase, CFE_TPOSELAST )
INSTANTIATE_CFEClassElast3 ( CFEGridBase, CFE_CDWI_ELAST )

#undef INSTANTIATE_CFEClass3
#undef INSTANTIATE_CFEClassElast32Arg
#undef INSTANTIATE_CFEClassElast3

}
