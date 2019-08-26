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

#ifndef __TPCFEGRID_H
#define __TPCFEGRID_H

#include <rectangularGrid.h>

#include <bitVector.h>
#include <scalarArray.h>
#include <multiArray.h>

#include <gridSize.h>
#include <indexMapper.h>
#include <iterators.h>

#include <tpCFEBasics.h>
#include <tpCFEElement.h>
#include <tpCFELookup.h>
#include <tpCFEVirtualNode.h>


namespace tpcfe {

/** Non-templated basis class for CFEGrid
 */
class CFEGridDefinition: public qc::RectangularGrid<qc::QC_3D> {
protected:
  std::vector< CFEType > _elType;            //!< Stores the types of the elements
  int _cfeGridIndexOffset[NUM_NEIGHBORS];    //!< This is the analogy to uniformGridIndexOffset.
  const bool _isCubic;
  const int _depth;

public:

  //! constructor from grid depth
  explicit CFEGridDefinition ( const int Depth ) : qc::RectangularGrid<qc::QC_3D> ( Depth, qc::QC_3D ), _isCubic ( true ), _depth ( Depth ) {
    initialize();
  }

  //! constructor from size in all space directions
  explicit CFEGridDefinition ( const qc::GridSize<qc::QC_3D> &Sizes ) : qc::RectangularGrid<qc::QC_3D> ( Sizes ), _isCubic ( ( Sizes[0] == Sizes[1] && Sizes[1] == Sizes[2] ) ), _depth ( -1 ) {
    initialize();
  }

  ~CFEGridDefinition () {}

  /** Return the vector of element types
   *  \warning provides direct access to protected member
   */
  std::vector< CFEType >::const_iterator elTypeData() const {
    return ( _elType.begin() );
  }

  CFEType getElType ( const int x, const int y, const int z ) const {
    return ( _elType[ getIndexMapperRef().getGlobalIndex ( x, y, z ) ] );
  }

  CFEType getElType ( const qc::CoordType &pos ) const {
    return ( _elType[ getIndexMapperRef().getGlobalIndex ( pos ) ] );
  }

  //! numX, numY and numZ may be the same or different
  inline int getNumXYZ ( ) const {
    if ( ! _isCubic ) {
      throw aol::OperationNotPermittedException ( "tpcfe::CFEGridDefinition::getNumXYZ() may not be called if grid does not represent cubic domain and was created as such.", __FILE__, __LINE__ );
    }

    return ( this->getNumX() );
  }

  int getGridDepth() const {
    if ( ! _isCubic ) {
      throw aol::OperationNotPermittedException ( "tpcfe::CFEGridDefinition::getGridDepth() may not be called if grid does not represent cubic domain and was created as such.", __FILE__, __LINE__ );
    }

    return ( _depth );
  }

  virtual short getElementGridDepth ( ) const {
    return ( _depth );
  }

protected:
  void initialize () {
    _elType.resize ( this->getNumberOfNodes() );

    for ( short i = 0; i < NUM_NEIGHBORS; ++i ) {
      _cfeGridIndexOffset[i] = qc::ILexCombine3 ( CFETopoLookup::_tetNbOffsets[i][0], CFETopoLookup::_tetNbOffsets[i][1], CFETopoLookup::_tetNbOffsets[i][2], this->getNumX(), this->getNumY() );
    }
  }

};


/** Class for grids with interfaces
 *  A Grid contains a domain and it contains Structures
 */
template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType >
class CFEGridBase: public CFEGridDefinition {
public:
  // re-export template parameters:
  typedef _RealType RealType;
  static const tpcfe::ConstraintType CT = _CT;

  // typedefs:
  typedef qc::BitArray<qc::QC_3D>                                                 BitArrayType;
  typedef tpcfe::CFEVirtualNode< RealType, _CT, NodalCoeffType >                  VNType;
  typedef typename VNType::IndexType                                              VNIndexType;
  typedef std::map<typename VNType::IndexType, VNType*>                           VNMapType;
  typedef typename std::map<typename VNType::IndexType, VNType*>::iterator        VNMapTypeIt;
  typedef typename std::map<typename VNType::IndexType, VNType*>::const_iterator  VNMapTypeConstIt;

  bool _inadmissibleSignatureWarningPrinted;

protected:
  // members
  std::vector < CFEStructure<RealType>* > _structure;         //! here, we may store multiple level set functions describing regions with different coefficients
  aol::BitVector                          _deleteStructures;  //!< bookkeeping whether structure should be deleted in destructor or belongs to someone else, MAX_STRUCT_ID refers to domain
  bool                                    _virtualNodesInitialized;
  bool                                    _verboseMode;
  VNMapType                               _vNode;

  // use virtual index MAX_STRUCT_ID to refer to domain.
  CFEStructure<RealType>                  *_pDomain;          //!< this is used for describing the domain of computation (e.g. based on a level set).

  // Here come the CFE_CD member variables

  BitArrayType                            *_pDomainNodeMask;  //!< nodes on which computation is performed in case of complicated domains (interior and neighboring), NULL: all inside, true <=> node inside
  BitArrayType                            *_pDirichletMask;   //!< dirichlet nodes, NULL: no dirichlet nodes, true <=> dirichlet node
  mutable BitArrayType                    *_pDOFMask;         //!< degrees of freedom. note: except for finest grid, *not* all non-Dirichlet domainNodes are DOFs.

  /*  Only domainNodes can be dirichletNodes.
   *  A domainNode on a coarse grid need not be domain node on finer grids.
   *  A Dirichlet node on a coarse grid must be dirichlet node on finer grids. I. e. coarsening does not introduce new Dirichlet nodes.
   *  Maybe the name domainNode is misleading - one layer of nodes around the interior of the domain of computation consists of domainNodes.
   *  We could call them computationalNodes, but then Dirichlet nodes would be computationalNodes. Not significantly less confusing ...
   */

  //! If levelset function is closer to zero than this value (in 1), it is adjusted by this value.
  RealType _adjustLevelset;

  //! If virtual node is closer to regular grid point than this value times the grid width, it is moved away from the regular grid point.
  RealType _adjustVirtualNodes;

  //! Store whether grid was created by multigrid solver and represents coarsened level (i.e. no level set function describing geometry present). This is important for creating appropriate matrices.
  bool _isCoarsened;

public:

  // --- con/destructors ---

  /** Construct a CFE grid which does not contain any structure
   */
  explicit CFEGridBase ( const int Depth );

  explicit CFEGridBase ( const qc::GridSize<qc::QC_3D> &Sizes );

  ~CFEGridBase ();

  // --- "small" access functions ---

  void setVerboseMode ( const bool whichMode ) {
    _verboseMode = whichMode;
  }

  /** Return reference to the map of virtual nodes
   */
  const VNMapType& virtualNodeMapRef() const {
    return ( _vNode );
  }

  /** Return Reference to the virtual node from the map which is
   *  interpolated by the given global dofs
   */
  const VNType& getVirtualNodeRef ( const int gIdx0, const int gIdx1 ) const {
    const typename VNType::IndexType hashIdx = VNType::mapIndex ( gIdx0, gIdx1 );

    typename VNMapType::const_iterator vnIt = _vNode.find ( hashIdx );

    if ( vnIt == _vNode.end() ) {
      throw aol::Exception ( "Virtual node not found", __FILE__, __LINE__ );
    }

    return ( * ( vnIt->second ) );
  }


  VNType& getVirtualNodeRef ( const int gIdx0, const int gIdx1 ) {
    typename VNType::IndexType hashIdx = VNType::mapIndex ( gIdx0, gIdx1 );

    typename VNMapType::const_iterator vnIt = _vNode.find ( hashIdx );

    if ( vnIt == _vNode.end() ) {
      throw aol::Exception ( "Virtual node not found", __FILE__, __LINE__ );
    }

    return ( * ( vnIt->second ) );
  }


  void setAdjustVirtualNodes ( const RealType AdjustVirtualNodes ) {
    _adjustVirtualNodes = AdjustVirtualNodes;
  }

  RealType getAdjustVirtualNodes ( ) const {
    return ( _adjustVirtualNodes );
  }

  //! must be called before setting domain or adding structure
  void setAdjustLevelset ( const RealType AdjustLevelset ) {
    _adjustLevelset = AdjustLevelset;
  }

  inline ConstraintType getCT ( ) const {
    return ( _CT );
  }

  /** Set reference to the computational domain (not copied)
   *  For the CFE_DOMAIN mode this method also computes the remapping of standard degrees of freedom
   *  \attention Changing the Structure afterwards may not have the desired effects on the grid.
   */
  void setDomainRef ( tpcfe::CFEStructure<RealType>& Domain ) {
    setDomainRef ( Domain, false );
  }

  /** Add a structure constructed from a ScalarArray<QC_3D> to this grid (level set data will be deep copied)
   */
  void addStructureFrom ( const qc::ScalarArray<RealType, qc::QC_3D> &Strct ) {
    CFEStructureAffine<RealType> *Structure = new CFEStructureAffine<RealType>;
    Structure->setImageFrom ( Strct, _adjustLevelset );
    _deleteStructures.set ( _structure.size(), true );
    addStructureRef ( *Structure );
  }

  /** Add reference to a CFEStructure
  */
  void addStructureRef ( CFEStructure<RealType> &Strct ) {
    if ( _structure.size() > ( tpcfe::MAX_STRUCT_ID - 2 ) ) {
      throw aol::Exception ( "Cannot add another structure.", __FILE__, __LINE__ );
    }
    if ( _structure.size() == 1 ) {
      throw aol::Exception ( "tpcfe::CFEGridBase: trying to add more than one structure. This is completely untested.", __FILE__, __LINE__ );
    }
    _structure.push_back ( &Strct );
    setElTypes ( _structure.size () - 1 );
  }

  /** Get Reference to the computational domain
   */
  const CFEStructure<RealType>& getDomainRef() const {
    if ( _pDomain == NULL ) {
      throw aol::Exception ( "No domain set.", __FILE__, __LINE__ );
    }
    return ( *_pDomain );
  }

  /** Return reference to the structure for the given type (not for index of structure!)
   */
  const CFEStructure<RealType>& getStructureRefForType ( const CFEType type ) const {
    const short structureNo = type._structureNo;
    return ( getStructureRef ( structureNo ) );
  }

  /** Return reference to structure identified by number
   */
  const CFEStructure<RealType>& getStructureRef ( const short num ) const {
    if ( num == MAX_STRUCT_ID ) {
      if ( _pDomain == NULL ) {
        throw aol::Exception ( "No domain set.", __FILE__, __LINE__ );
      }
      return ( *_pDomain );
    } else {
      if ( num < 0 || num >= MAX_STRUCT_ID ) {
        throw aol::OutOfBoundsException ( "Cannot return this structure", __FILE__, __LINE__ );
      }
      if ( _structure[num] == NULL ) {
        throw aol::Exception ( "Structure not set.", __FILE__, __LINE__ );
      }
      return ( * ( _structure[num] ) );
    }
  }

  /** removes domain from the grid
  */
  void removeDomain () {
    removeStructureFromElTypes ( MAX_STRUCT_ID );
    if ( _deleteStructures[ MAX_STRUCT_ID ] )
      delete ( _pDomain );
    _pDomain = NULL;
  }

  /** Remove a structure from the grid. This modifies the indices of other structures. Untested.
   */
  void removeStructure ( const int which ) {
    cerr << "CFEGridBase::removeStructure has not been tested (remove this output if it has been tested)." << endl;
    if ( which >= static_cast<int> ( _structure.size() ) ) {
      throw aol::Exception ( "Cannot remove structure. Structure does not exist." );
    }

    removeStructureFromElTypes ( which );

    // Remove pointer to structure object from _structure vector
    _structure.erase ( _structure.begin() + which );
  }

  /** Check whether a domain is in use or not
   */
  bool hasDomain() const { return ( _pDomain != NULL ); }

  /** Return number of structures (i. e. interfaces for jumping coefficients)
   */
  int getNumStructures ( ) const {
    return ( _structure.size() );
  }

  /** Detection of slave and virtual nodes
   */
  void detectVirtualNodes();

  /** Performs detection and initialization of weights and constraints
   *  in the jumping coefficient case for virtual nodes based on an
   *  AArray of coefficients at the grid nodes
   */
  void detectAndInitVirtualNodes ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeff ) {
    if ( ( CT == CFE_TPOS ) || ( CT == CFE_TPOSELAST ) || ( CT == CFE_LIEHR ) || ( CT == CFE_CDWI_TPOS ) || ( CT == CFE_CDWI_LIEHR ) || ( CT == CFE_CDWI_ELAST ) ) {
      this->detectVirtualNodes();
      this->initVirtualNodes ( nodalCoeff );
    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::detectAndInitVirtualNodes ( nodalCoeff ): illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  /** Try to detect and initialize virtual nodes with initial maximal
      inversion error threshold in local construction, increasing the
      threshold until initialization was successful.
   */
  RealType relaxedDetectAndInitVirtualNodes ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeff, const RealType initThreshold, const RealType thresholdStep );


  // --- CFE_CD specific methods ---

  /** Get reference to the domain node mask
   */
  const BitArrayType& getDomainNodeMaskRef() const {
    if ( !_pDomainNodeMask ) {
      throw aol::Exception ( "No domain node mask set", __FILE__, __LINE__ );
    }
    return ( ( *_pDomainNodeMask ) );
  }

  /** Returns whether the given node lies inside a domain
   */
  inline bool isInsideDomain ( const qc::CoordType &pos ) const {
    return ( getDomainRef().getValue ( pos ) < 0 );
  }

  //! this function only sets the domainNode mask but does not initialize any virtual nodes or things like that.
  void setDomainNodeMask ( BitArrayType &dmask ) {
    if ( &dmask == _pDomainNodeMask ) {
      throw aol::Exception ( "tpcfe::CFEGridBase::setDomainNodeMask: cannot set _pDomainNodeMask to itself", __FILE__, __LINE__ );
      // maybe just do nothing - but I'm not sure whether that's what we want.
    }

    if ( _pDomainNodeMask ) {
      delete ( _pDomainNodeMask );
    }
    _pDomainNodeMask = new BitArrayType ( dmask ); // deep copy
  }

  void setDirichletMask ( BitArrayType &dmask ) {
    if ( &dmask == _pDirichletMask ) {
      throw aol::Exception ( "tpcfe::CFEGridBase::setDirichletMask: cannot set _pDirichletMask to itself", __FILE__, __LINE__ );
      // maybe just do nothing - but I'm not sure whether that's what we want.
    }

    if ( _pDirichletMask ) {
      delete ( _pDirichletMask );
    }
    _pDirichletMask = new BitArrayType ( dmask ); // deep copy
  }

  void setDOFMask ( BitArrayType &dmask ) {
    if ( &dmask == _pDOFMask ) {
      throw aol::Exception ( "tpcfe::CFEGridBase::setDOFMask: cannot set _pDOFMask to itself", __FILE__, __LINE__ );
      // maybe just do nothing - but I'm not sure whether that's what we want.
    }

    if ( _pDOFMask ) {
      delete ( _pDOFMask );
    }
    _pDOFMask = new BitArrayType ( dmask ); // deep copy
  }

  //! this function is not called automatically after domainNodeMask or DirichletMask is set, but if isDOFNode requires DOF mask to be set.
  //! all non-Dirichlet domain nodes are set to be DOFs
  void setDOFMaskFromDirichletAndDomainNodeMask( ) const { // possibly called from isDOFNode -> must be const -> _pDOFMask must be mutable
    if ( _pDOFMask ) {
      delete ( _pDOFMask );
    }

    _pDOFMask = new BitArrayType ( *this );
    for ( int i = 0; i < this->getNumberOfNodes(); ++i ) {
      _pDOFMask->set ( i, ( this->isDomainNode ( i ) && ! ( this->isDirichletNode ( i ) ) ) );
    }
  }

  void setDOFMaskFromDirichletMask ( ) const {
    if ( _pDOFMask ) {
      delete ( _pDOFMask );
    }

    _pDOFMask = new BitArrayType ( *this );
    for ( int i = 0; i < this->getNumberOfNodes(); ++i ) {
      _pDOFMask->set ( i, ! ( this->isDirichletNode ( i ) ) );
    }

  }

  bool hasDirichletMask() const { return ( _pDirichletMask != NULL ); }

  bool hasDomainNodeMask() const { return ( _pDomainNodeMask != NULL ); }

  bool hasDOFMask() const { return ( _pDOFMask != NULL ); }

  const BitArrayType& getDirichletMaskRef() const {
    if ( !_pDirichletMask ) {
      throw aol::Exception ( "No Dirichlet mask set", __FILE__, __LINE__ );
    }
    return ( ( *_pDirichletMask ) );
  }

  const BitArrayType& getDOFMaskRef() const {
    if ( !_pDOFMask ) {
      throw aol::Exception ( "No DOFMask set", __FILE__, __LINE__ );
    }
    return ( ( *_pDOFMask ) );
  }

  //! set all entries of the given BitArray<qc::QC_3D> to true where levelset function has negative values
  void setInnerNodeMask ( qc::BitArray<qc::QC_3D> &innerNodes ) const {
    if ( innerNodes.size() != this->getNumberOfNodes() ) {
      throw aol::Exception ( "tpcfe::CFEGridBase::setInnerNodeMask: incorrect size of BitArray", __FILE__, __LINE__ );
    }
    if ( _pDomain == NULL ) {
      throw aol::Exception ( "No domain set.", __FILE__, __LINE__ );
    }

    for ( qc::RectangularIterator<qc::QC_3D> bit ( *this ); bit.notAtEnd(); ++bit ) {
      innerNodes.set ( *bit, ( _pDomain->getValue ( *bit ) < 0 ) ); // check for non-NULLness
    }

  }

  inline bool isDomainNode ( const int globalIndex ) const {
    if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
      if ( _pDomainNodeMask == NULL ) {
        throw aol::Exception ( "tpcfe::CFEGridBase::isDomainNode: no domain node mask set!", __FILE__, __LINE__ );
      }

      return ( _pDomainNodeMask->get ( globalIndex ) );

    } else if ( CT == CFE_LIEHR || CT == CFE_TPOS || CT == CFE_TPOSELAST ) {
      return ( true );
    } else {
      // maybe just return true?
      throw aol::Exception ( "tpcfe::CFEGridBase::isDomainNode: illegal constraint type", __FILE__, __LINE__ );
    }

  }

  inline bool isDomainNode ( const qc::CoordType &pos ) const {
    if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
      if ( _pDomainNodeMask == NULL ) {
        throw aol::Exception ( "tpcfe::CFEGridBase::isDomainNode: no domain node mask set!", __FILE__, __LINE__ );
      }

      return ( _pDomainNodeMask->get ( pos ) );

    } else if ( CT == CFE_LIEHR || CT == CFE_TPOS || CT == CFE_TPOSELAST ) {
      return ( true );
    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isDomainNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isDirichletNode ( const int globalIndex ) const {
    if ( CT == CFE_CD || CT == CFE_TPOS || CT == CFE_LIEHR || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_TPOSELAST || CT == CFE_CDWI_ELAST ) {
      if ( _pDirichletMask == NULL ) {
        return ( false );
      } else {
        return ( ( _pDirichletMask->get ( globalIndex ) ) );
      }
    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isDirichletNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isDirichletNode ( const qc::CoordType &pos ) const {
    if ( CT == CFE_CD || CT == CFE_TPOS || CT == CFE_LIEHR || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_TPOSELAST || CT == CFE_CDWI_ELAST ) {
      if ( _pDirichletMask == NULL ) {
        return ( false );
      } else {
        return ( ( _pDirichletMask->get ( pos ) ) );
      }
    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isDirichletNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isDOFNode ( const int globalIndex ) const {
    if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
      if ( _pDOFMask == NULL ) {
        this->setDOFMaskFromDirichletAndDomainNodeMask();
      }

      return ( ( _pDOFMask->get ( globalIndex ) ) );

    } else if ( CT == CFE_TPOS || CT == CFE_LIEHR || CT == CFE_TPOSELAST ) {

      return ( !isDirichletNode ( globalIndex ) );

    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isDOFNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isDOFNode ( const qc::CoordType &pos ) const {
    if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
      if ( _pDOFMask == NULL ) {
        this->setDOFMaskFromDirichletAndDomainNodeMask();
      }

      return ( ( _pDOFMask->get ( pos ) ) );

    } else if ( CT == CFE_TPOS || CT == CFE_LIEHR || CT == CFE_TPOSELAST ) {

      return ( !isDirichletNode ( pos ) );

    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isDOFNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isInnerNode ( const qc::CoordType &pos ) const {
    if ( CT == CFE_CD || CT == CFE_CDWI_TPOS || CT == CFE_CDWI_LIEHR || CT == CFE_CDWI_ELAST ) {
      if ( _pDomain == NULL ) {
        throw aol::Exception ( "tpcfe::CFEGridBase::isInnerNode: no domain set!", __FILE__, __LINE__ );
      }

      return ( _pDomain->getValue ( pos ) < 0.0 );

    } else if ( CT == CFE_TPOS ) {
      return ( true );
    } else {
      throw aol::Exception ( "tpcfe::CFEGridBase::isInnerNode: illegal constraint type.", __FILE__, __LINE__ );
    }
  }

  inline bool isInnerNode ( const int globalIndex ) const {
    qc::CoordType pos;
    getIndexMapperRef().splitGlobalIndex ( globalIndex, pos );
    return ( isInnerNode ( pos ) );
  }

  /** for CFE_TPOS, determine whether any of the 15 neighbors lies on the other side of interface struc (maybe check all?)
   */
  bool isNearInterface ( const qc::CoordType &pos, const short struc = 0 ) const {
    if ( CT != CFE_TPOS && CT != CFE_TPOSELAST && CT != CFE_CDWI_ELAST ) {
      throw aol::Exception ( "tpcfe::CFEGridBase::isNearInterface: illegal constraint type.", __FILE__, __LINE__ );
    }
    const RealType hereVal = getStructureRef ( struc ).getValue ( pos );
    for ( short nb = 0; nb < NUM_NEIGHBORS /* sic: CFE_CD neighbors only! */; ++nb ) {
      const qc::CoordType nbpos = pos + CFETopoLookup::_tetNbOffsets[nb];
      if ( this->isAdmissibleNode ( nbpos ) && getStructureRef ( struc ).getValue ( nbpos  ) * hereVal < 0 ) { // has opposite sign
        return ( true );
      }
    }
    return ( false );
  }

protected:
  //! Sets all entries of arg outside of domain of computation to zero
  void doRestrictToDomain ( aol::Vector<RealType> &arg ) const;

  void doRestrictToDomain ( aol::MultiVector<RealType> &arg ) const {
    for ( int i = 0; i < arg.numComponents(); ++i )
      doRestrictToDomain ( arg[i] );
  }

  void setDomainRef ( tpcfe::CFEStructure<RealType>& Domain, const bool deleteDomain );

public:
  //! Sets all entries of arg at Dirichlet nodes to zero
  // sure we want to use this?
  void restrictDirichletNodes ( aol::Vector<RealType> &arg ) const;

  void restrictDirichletNodes ( aol::MultiVector<RealType> &arg ) const {
    for ( int i = 0; i < arg.numComponents(); ++i )
      restrictDirichletNodes ( arg[i] );
  }

  //! Sets all entries of arg at non-dofs to zero (outside of domain of computation and Dirichlet nodes)
  // sure we want to use this?
  void restrictToDofs ( aol::Vector<RealType> &arg ) const;

  // sure we want to use this?
  void restrictToDofs ( aol::MultiVector<RealType> &arg ) const {
    for ( int i = 0; i < arg.numComponents(); ++i )
      restrictToDofs ( arg[i] );
  }

  void restrictToInnerNodes ( qc::ScalarArray<RealType, qc::QC_3D> &arg ) const;


  // --- CFE_jumping coefficient specific methods ---

  /** Fill the element with appropriate structure values at its vertices
   */
  bool getStructureValues ( CFEElement<RealType> &el ) const;


  // --- general stuff ---

  void getCutRelations ( CFEElement<RealType> &el ) const;

  /** Print out the number of structures
   */
  void dump ( ostream & out ) const {
    const int num = _structure.size ();
    out << "This grid has " << num << " structure" << ( ( num > 1 ) ? "s." : "." ) << endl;
  }

  /** Return whether an elemnt has type > 0 or not
   */
  bool isInterfaced ( const CFEElement < RealType > &e ) const;

  //! Return whether the edge between nodes is interfaced by structure.
  bool isEdgeInterfaced ( const qc::CoordType pos0, const qc::CoordType pos1, const short stru = 0 ) const {
    return ( getStructureRef ( stru ).getValue ( pos0 ) * getStructureRef ( stru ).getValue ( pos1 ) < 0 );
  }

  RealType getTotalInnerVolume ( ) const {
    return ( getTotalPartialVolume ( -1 ) );
  }

  RealType getTotalPartialVolume ( const signed char which ) const;

  int getNumDOFs ( ) const {
    // print warning if _pDOFMask not set but _DirichletMask set??
    if ( _pDOFMask == NULL )
      return ( this->getNumberOfNodes() );
    else
      return ( _pDOFMask->numTrue() );
  }

  inline int getCfeGridIndexOffset ( const int i ) const {
    return ( _cfeGridIndexOffset[i] );
  }

  //! Returns cut relation of structure stru between two nodes (given by global index). This only makes sense if the two nodes are connected by an edge (which would be expensive to check).
  RealType getCutRelationBetween ( const qc::CoordType& pos0, const qc::CoordType& pos1, const signed char stru ) const {
    return ( getStructureRef ( stru ).getCutRelationBetween ( pos0, pos1 ) );
  }

  void removeVirtualNodes ( ) {
    for ( VNMapTypeIt it = _vNode.begin(); it != _vNode.end(); ++it ) {
      delete ( ( *it ).second );
    }
    _vNode.clear();
  }

  //! For grids not on the finest level, a different matrix structure may be necessary.
  void setCoarsenedFlag ( const bool isCo ) {
    _isCoarsened = isCo;
  }

  bool isCoarsened ( ) const {
    return ( _isCoarsened );
  }

private:
  CFEGridBase<_RealType, _CT, NodalCoeffType > ( const CFEGridBase<_RealType, _CT, NodalCoeffType > &other );
  CFEGridBase<_RealType, _CT, NodalCoeffType >& operator= ( const CFEGridBase<_RealType, _CT, NodalCoeffType >& );

  int getUniformGridIndexOffset ( int /*i*/ ) const {
    throw aol::Exception ( "getUniformGridIndexOffset: This function must not be used for tpCFE grids" , __FILE__, __LINE__ );
  }

protected:
  /** Performs the initialization of weights and constrains of virtual nodes
   */
  void initVirtualNodes ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeff );

  /** Must be called whenever a structure has been added
   */
  void setElTypes ( int structureNo );

  /** adjust the element types if a struture has been removed
   */
  void removeStructureFromElTypes ( const signed char which ) {
    /*  int          adjustCounter = 0; */
    /*  CFEStructure<RealType> *s = ( which == MAX_STRUCT_ID ? _pDomain : _structure[which] ); */

    // Remove all element type <> 0 that have been caused by this structure
    for ( unsigned int i = 0; i < _elType.size(); ++i ) {
      const signed char sno = _elType[i]._structureNo;
      if ( sno == which ) {
        _elType[i] = CFEType();
      }
    }
  }

  /** Fill the element with structure values of the structure s at its vertices
   */
  void getStructureValues ( const CFEStructure<RealType> &s, CFEElement<RealType> &el ) const {
    for ( int i = 0; i < 8; ++i )
      el.setStructureValue ( i, s.getValue ( el + CFETopoLookup::_hexNbOffsets[i] ) );
  }


  /** Tells whether an element is interfaced or not.
   *  Please note that it is much faster to query the cfeType
   *  of the element than calling this function. For any cfeType
   *  other than 0x00 and 0xFF the element is interfaced!
   */
  bool interfaced ( const CFEElement<RealType> &el, int parentTetra ) const;

};


//! class to collect methods necessary if a complicated domain boundary is present (CFE_CD and CFE_CDWI_*)
template < typename _RealType, tpcfe::ConstraintType _CT, typename _NodalCoeffType >
class CFEGridBase_CD : public CFEGridBase < _RealType, _CT, _NodalCoeffType > {
public:
  typedef _RealType RealType;
  typedef _NodalCoeffType NodalCoeffType;

  explicit CFEGridBase_CD ( const int Depth ) : CFEGridBase< _RealType, _CT, _NodalCoeffType > ( Depth ) {
  }

  explicit CFEGridBase_CD ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, _CT, _NodalCoeffType > ( Sizes ) {
  }

  void setDomainFrom ( const qc::ScalarArray<RealType, qc::QC_3D> &Domain ) {
    CFEStructureAffine<RealType>* pStrAff = new CFEStructureAffine<RealType>;
    pStrAff->setImageFrom ( Domain, this->_adjustLevelset );
    this->setDomainRef ( *pStrAff, true );
  }

  void restrictToDomain ( aol::Vector<RealType> &arg ) const {
    this->doRestrictToDomain ( arg );
  }

  void restrictToDomain ( aol::MultiVector<RealType> &arg ) const {
    this->doRestrictToDomain ( arg );
  }

};


//! For template specialization
template < typename _RealType, tpcfe::ConstraintType _CT, typename NodalCoeffType = _RealType >
class CFEGrid : public CFEGridBase < _RealType, _CT, NodalCoeffType > { };

template < typename _RealType >
class CFEGrid <_RealType, CFE_NONE, _RealType > : public CFEGridBase < _RealType, CFE_NONE, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_NONE, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_NONE, _RealType > ( Sizes ) {
  }

  void setDomainFrom ( const qc::ScalarArray<RealType, qc::QC_3D> & Domain ) {
    this->removeDomain();

    CFEStructureAffine<RealType>* tmp = new CFEStructureAffine<RealType>;
    tmp->setImageFrom ( Domain, this->_adjustLevelset );
    this->_pDomain = tmp;

    this->setElTypes ( MAX_STRUCT_ID );
    this->_deleteStructures.set ( MAX_STRUCT_ID, true );
  }
};


template < typename _RealType, typename _NodalCoeffType >
class CFEGrid <_RealType, CFE_CD, _NodalCoeffType > : public CFEGridBase_CD < _RealType, CFE_CD, _NodalCoeffType > {
public:
  typedef _RealType RealType;
  typedef _NodalCoeffType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase_CD< _RealType, CFE_CD, _NodalCoeffType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase_CD< _RealType, CFE_CD, _NodalCoeffType > ( Sizes ) {
  }

  /** Performs the detection and initialization of weights and
   *  constraints of virtual nodes in the CFE_CD case. This requires
   *  dummy coefficients and is encapsulated here.
   */
  void detectAndInitVirtualNodes ( ) {
    this->detectVirtualNodes();

    qc::AArray< NodalCoeffType, qc::QC_3D > dummyCoeff ( *this );
    dummyCoeff.setAll ( aol::ZOTrait<NodalCoeffType>::zero  );
    this->initVirtualNodes ( dummyCoeff );
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_LIEHR, _RealType > : public CFEGridBase < _RealType, CFE_LIEHR, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_LIEHR, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_LIEHR, _RealType > ( Sizes ) {
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_TPOS, _RealType > : public CFEGridBase < _RealType, CFE_TPOS, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_TPOS, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_TPOS, _RealType > ( Sizes ) {
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_TP, _RealType > : public CFEGridBase < _RealType, CFE_TP, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_TP, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_TP, _RealType > ( Sizes ) {
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_CDWI_LIEHR, _RealType > : public CFEGridBase_CD < _RealType, CFE_CDWI_LIEHR, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase_CD < _RealType, CFE_CDWI_LIEHR, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase_CD < _RealType, CFE_CDWI_LIEHR, _RealType > ( Sizes ) {
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_CDWI_TPOS, _RealType > : public CFEGridBase_CD < _RealType, CFE_CDWI_TPOS, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase_CD < _RealType, CFE_CDWI_TPOS, _RealType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase_CD < _RealType, CFE_CDWI_TPOS, _RealType > ( Sizes ) {
  }

};


template < typename _RealType, typename _NodalCoeffType >
class CFEGrid <_RealType, CFE_TPOSELAST, _NodalCoeffType > : public CFEGridBase < _RealType, CFE_TPOSELAST, _NodalCoeffType > {
public:
  typedef _RealType RealType;
  typedef _NodalCoeffType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_TPOSELAST, NodalCoeffType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_TPOSELAST, NodalCoeffType > ( Sizes ) {
  }
};


template < typename _RealType, typename _NodalCoeffType >
class CFEGrid <_RealType, CFE_CDWI_ELAST, _NodalCoeffType > : public CFEGridBase_CD < _RealType, CFE_CDWI_ELAST, _NodalCoeffType > {
public:
  typedef _RealType RealType;
  typedef _NodalCoeffType NodalCoeffType;

  explicit CFEGrid ( const int Depth ) : CFEGridBase_CD < _RealType, CFE_CDWI_ELAST, NodalCoeffType > ( Depth ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase_CD < _RealType, CFE_CDWI_ELAST, NodalCoeffType > ( Sizes ) {
  }
};


template < typename _RealType >
class CFEGrid <_RealType, CFE_DOMAIN, _RealType > : public CFEGridBase < _RealType, CFE_DOMAIN, _RealType > {
public:
  typedef _RealType RealType;
  typedef _RealType NodalCoeffType;

  typedef tpcfe::CFEVirtualNode< RealType, CFE_DOMAIN, RealType >     VNType;
  typedef std::map<typename VNType::IndexType, VNType*>            VNMapType;

protected:
  aol::Vector<int>                        *_pDOFMap;  //!< Remapping of the degrees of freedom in case of domain mode
  int                                     _domainSize;
  mutable int                             _globalDof;
  mutable VNMapType                       _gDofVNode;

public:
  explicit CFEGrid ( const int Depth ) : CFEGridBase< _RealType, CFE_DOMAIN, RealType > ( Depth ),
      _pDOFMap ( NULL ), _domainSize ( 0 ), _globalDof ( -1 ) {
  }

  explicit CFEGrid ( const qc::GridSize<qc::QC_3D> &Sizes ) : CFEGridBase< _RealType, CFE_DOMAIN, RealType > ( Sizes ),
      _pDOFMap ( NULL ), _domainSize ( 0 ), _globalDof ( -1 ) {
  }

  ~CFEGrid ( ) {
    if ( _pDOFMap )
      delete _pDOFMap;
  }

  /** Return the new (remapped) global index for this old index.
   *  Makes sense only in CFE_DOMAIN mode where the dofs have been remapped
   */
  int remappedDof ( const int oldIndex ) const {
    return this->_pDOFMap->get ( oldIndex );
  }

  /** fill a vector with original dof numbering with the values of a vector with
   *  renumbered dofs. Only the standard nodex of the grid can be assigned.
   */
  void mapVectorBack ( aol::Vector<RealType> &img, aol::Vector<RealType> &gImg ) {
    const int gSize = gImg.size();
    qc::CoordType pos;
    for ( int i = 0; i < gSize; ++i ) {
      this->getIndexMapperRef().splitGlobalIndex ( i, pos );
      if ( this->isInsideDomain ( pos ) ) {
        gImg.set ( i, img.get ( remappedDof ( i ) ) );
      }
    }
  }

  /** This method is called successively to distribute the free global dofs among the
   *  virtual modes if the grid runs in the CFE_DOMAIN mode
   */
  int getFreeGlobalDof ( VNType* vNode ) const {
    int tmp = this->_globalDof;
    _gDofVNode[tmp] = vNode;
    ++_globalDof;
    return tmp;
  }

  /** Return the virtual node object that has the given global index
   *  Makes sense only in CFE_DOMAIN mode
   *  \todo name misleading
   */
  VNType *virtualNodeForGlobalDof ( const int index ) const {
    typename VNMapType::const_iterator vnIt = _gDofVNode.find ( index );

    if ( vnIt == _gDofVNode.end() )
      return NULL;
    else
      return vnIt->second;
  }

  /** Return the size of the domain vector
   */
  int domainSize() const { return _domainSize; }

  /** Return largest global index (makes sense in domain mode CFE_DOMAIN only)
   */
  int largestDofIndex() const {
    return ( this->getNumberOfNodes() );
  }

  void setDomainFrom ( const qc::ScalarArray<RealType, qc::QC_3D>& Domain ) {
    this->removeDomain();
    {
      CFEStructureAffine<RealType>* tmp = new CFEStructureAffine<RealType>;
      tmp->setImageFrom ( Domain, this->_adjustLevelset );
      this->_pDomain = tmp;
    }

    this->setElTypes ( MAX_STRUCT_ID );
    this->_deleteStructures.set ( MAX_STRUCT_ID, true );

    _globalDof = 0;
    _domainSize = Domain.size();
    _pDOFMap = new aol::Vector<int> ( _domainSize );

    for ( qc::RectangularIterator<qc::QC_3D> bit ( *this ); bit.notAtEnd(); ++bit ) {
      if ( this->getDomainRef().getValue ( *bit ) < 0 ) {
        _pDOFMap->set ( this->getIndexMapperRef() ( *bit ), _globalDof );
        /*cerr << _globalDof << " "; */
        ++_globalDof;
      }
    }
    if ( _globalDof == 0 )
      throw aol::Exception ( "Domain does not contain any dofs", __FILE__, __LINE__ );
  }

  /** This iterator runs over the nodes of the given grid, but returns only those, which are inside
   *  the domain.
   *  In CFE_DOMAIN mode it additionally returns all the virtual nodes (after the regular nodes have
   *  been traversed)
   */
  class DomainNodeIterator {
  protected:
    typedef CFEGrid <_RealType, CFE_DOMAIN, _RealType > ParentGridType;
    aol::Vec3<RealType>          _coord;
    const ParentGridType        &_grid;
    int                          _id;
    int                          _index;
    typename VNMapType::iterator _vnIt;
    bool                         _hasVirtualNodes;
    bool                         _isVirtualNode;

  public:
    DomainNodeIterator ( const ParentGridType &grid, const int id = 0 ) :
        _coord(), _grid ( grid ), _id ( id ),
        _index( ), _vnIt ( ), _hasVirtualNodes ( _grid._gDofVNode.size() != 0 ),
        _isVirtualNode ( false ) {
    }

    /** Returns the current (remapped) index.
     */
    int operator*() const { return _index; }

    /** Increments the iterator and returns the remapped index
     */
    int operator++() {
      ++_id;
      // We are in the range of regular nodes
      if ( _id < _grid._domainSize ) {
        qc::CoordType idpos;
        _grid.getIndexMapperRef().splitGlobalIndex ( _id, idpos );
        while ( _id < _grid._domainSize && !_grid.isInsideDomain ( idpos ) ) ++_id;
        if ( _id < _grid._domainSize ) _index = _grid.remappedDof ( _id );
      }
      // Now comes the first virtual node
      if ( _id == _grid._domainSize && _hasVirtualNodes ) {
        VNType *vNode = _vnIt->second;
        _index = -1 * vNode->constraint ( 0 );
        return _index;
      }
      // Here come the rest of the virtual nodes.
      if ( _id > _grid._domainSize && _hasVirtualNodes ) {
        ++_vnIt;
        if ( _vnIt != _grid._gDofVNode.end() ) {
          VNType *vNode = _vnIt->second;
          _index = -1 * vNode->constraint ( 0 );
        }
      }
      return _index;
    }

    /** Return the global coordinates of the active node
     */
    aol::Vec3<RealType>& coord() {
      qc::CoordType idpos;
      _grid.getIndexMapperRef().splitGlobalIndex ( _id, idpos );
      if ( _id < _grid._domainSize && _grid.isInsideDomain ( idpos ) ) {
        _coord = aol::Vec3<RealType> ( idpos );
      } else {
        VNType *vNode;
        vNode = _vnIt->second;
        _coord = vNode->_coord;
      }
      return _coord;
    }

    bool isVirtualNode() const {
      return ( _id >= _grid._domainSize );
    }

    void begin() {
      _vnIt = _grid._gDofVNode.begin();
      _id = -1;
      operator++();
    }

    bool end() {
      bool result;
      if ( _hasVirtualNodes ) result = ( _vnIt == _grid._gDofVNode.end() );
      else result = ( _id >= _grid._domainSize - 1 );
      return result;
    }

    // end of class DomainNodeIterator
  };

};

}
#endif
