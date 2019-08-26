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

#ifndef __TPCFEVIRTUALNODE_H
#define __TPCFEVIRTUALNODE_H

#include <scalarArray.h>
#include <multiArray.h>
#include <vectorExtensions.h>

#include <tpCFEBasics.h>
#include <tpCFEElement.h>

namespace tpcfe {

class ParentTetraIdentifier {
  const CFETopoElement& _el;
  const unsigned char     _tnum;

public:
  ParentTetraIdentifier ( const CFETopoElement &el, const unsigned char tnum ) : _el ( el ), _tnum ( tnum ) {}

public:
  inline bool operator< ( const ParentTetraIdentifier &other ) const {
    return ( ( _el < other._el ) || ( ( _el == other._el ) && ( _tnum < other._tnum ) ) );
  }
};

/** Storage class for virtual nodes.
 *  Each virtual node has a map< CFETopoElement, ParentTetraList > storing the indices of regular tetrahedra inside the element that constrain the virtual node (at most six)
 *  \author Schwen
 */
class ParentTetraList {
  unsigned char numPT;
  unsigned char parentTetra[6];

public:
  ParentTetraList ( ) : numPT ( 0 ) {
    for ( short i = 0; i < 6; ++i )
      parentTetra[i] = tpcfe::NODE_NOT_VIRTUAL; // identifier for unset parent tetra
  }

  // default destructor, default copy constructor and assignment operator do the right thing

  //! add parent tetrahedron
  void addPT ( const unsigned char ptn ) {
    parentTetra[ numPT ] = ptn;
    ++numPT;
  }

  unsigned char getNumPT ( ) const {
    return ( numPT );
  }

  unsigned char operator[] ( const unsigned char i ) const {
    return ( parentTetra[i] );
  }
};


template< typename RealType >
class CFEGeomVirtualNode {
public:

  // To avoid conflicts with the ISO C++ standard we need the
  // corresponding substitute uint64_t for "unsigned long long int", a 64 bit integer.
  typedef uint64_t IndexType;

  /** This type stores tetrahedra which contain this virtual node, since the virtual node can be present
   *  in multiple elements and multiple tetrahedra.
   */
  typedef map< CFETopoElement, ParentTetraList > PTMapType;

  static RealType _maximalInversionError;   //!< when inverting a matrix for constructing extrapolation weights, any matrix times its inverse with bigger Frobenius difference from identity is rejected.

  bool _isDirichlet;
  RealType _dirichletValue;

  aol::Vec3<RealType> _coord;               //!< Global coordinates of this node
  aol::Vec3<RealType> _normal;              //!< Approximated normal to interface at this node

  int _numOfConstraints;                    //!< number of constraining vertices; Length of _weight and _constraint vectors

  aol::Vector<int>  _constraints;           //!< Indices of the nodes that constrain this one.

  short _local0, _local1;                   //!< Local indices of geometrically interpolating nodes
  int   _global0, _global1;                 //!< Global indices of geometrically interpolating nodes

  PTMapType   _parentTetra;                 //!< Stores elements and parent tetrahedra in which this node lies


  //! Constructor that should be used
  CFEGeomVirtualNode ( const short local0, const short local1,
                       const int global0, const int global1,
                       const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra ) :
      _isDirichlet ( false ), _dirichletValue ( 0 ), _numOfConstraints ( 0 ), _constraints ( 0 ),
      _local0 ( local0 ), _local1 ( local1 ),
      _global0 ( global0 ), _global1 ( global1 ),
      _parentTetra() {

    ParentTetraList pt;
    pt.addPT ( ParentTetra );

    _parentTetra.insert ( pair<CFETopoElement, ParentTetraList> ( CFETopoElement ( inEl ), pt ) );

    determineGlobalCoordinates ( inEl );
  }

  CFEGeomVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEGeomVirtualNode ( const CFEGeomVirtualNode &other );
  CFEGeomVirtualNode& operator= ( const CFEGeomVirtualNode& );

public:
  /** Generate an index for the stl-map for a node being interpolated by node0 and node1 that generates indices not used for regular nodes (virtual node [0,a] is mapped to a * 2 power 32)
   */
  static IndexType mapIndex ( const int node0, const int node1 );

  /** Split index to two integers (node0 > node 1). index was concatenated either from (node0 and node1) or (node1 and node0)
   */
  static void splitIndex ( const IndexType index, int &node0, int &node1 );

  /** Return the index for the stl-map for this node
   */
  IndexType mapIndex() const {
    return mapIndex ( _global0, _global1 );
  }

  /** Return how much memory this virtual node needs (hopefully correct)
   */
  int memoryUsed() const {
    return ( sizeof ( CFEGeomVirtualNode< RealType > ) +
             _constraints.size() * sizeof ( int ) +
             _parentTetra.size() * sizeof ( ParentTetraList ) );
  }

  //! Mark virtual node as Dirichlet virtual node. Probably untested!
  void setDirichlet ( bool dirichlet, RealType value ) { _isDirichlet = dirichlet; _dirichletValue = value; }

  RealType getDirichlet ( ) const { return _dirichletValue; }

  bool isDirichlet() const { return _isDirichlet; }
  /** Add another tetrahedron and element in which this virtual node lies
   */
  void addInElement ( const tpcfe::CFEElement<RealType> &inEl, const int ParentTetra );

  /** Return the first element to which this node belongs
   */
  const CFETopoElement &getInElement ( void ) const {
    return _parentTetra.begin()->first;
  }

  /** Return the parent tetrahedron in the first element to which this node belongs
   */
  int getParentTetra ( void ) const {
    return ( _parentTetra.begin()->second ) [0];
  }

  /** Return the number of elements that carry this virtual node
   */
  int numInElements() const { return _parentTetra.size(); }

  /** Return the number of constraints for this node
   */
  int numConstraints ( ) const { return _numOfConstraints; }

  /** Return the i-th constraint
   */
  int constraint ( const int i ) const { return _constraints[i]; }

  /** Print the constraining parent-tetrahedra that have been found
   */
  void printParentTetra ( ostream &out ) const {
    for ( typename PTMapType::const_iterator vnIt = _parentTetra.begin(); vnIt != _parentTetra.end(); ++vnIt ) {
      out << vnIt->first << endl << "\t";
      for ( int i = 0; i < vnIt->second.getNumPT(); ++i ) {
        out << static_cast<int> ( ( vnIt->second ) [i] ) << ", ";
      }
      out << endl;
    }
    out << "---------------" << endl;
  }

protected:
  /** Compute grid coordinates of a node which is constrained by two other nodes.
   *  The method takes an element and computes the coordinates of the constraining nodes as well as the
   *  coordinates of the constrained nodes.
   */
  void determineGlobalCoordinates ( const CFEElement<RealType> &el );

};


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
class CFEVirtualNodeBase : public CFEGeomVirtualNode < RealType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;
  typedef typename tpcfe::ExtrapolationWeightTypeTrait< RealType, CT >::WeightType ExtrapolationWeightType;

  aol::RandomAccessContainer< ExtrapolationWeightType > _weights;                                                      //!< Weights for the nodes this virtual one is constrained by

  CFEVirtualNodeBase ( const short local0, const short local1,
                       const int global0, const int global1,
                       const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEGeomVirtualNode<RealType> ( local0, local1, global0, global1, inEl, ParentTetra ), _weights ( 0 ) {}

  CFEVirtualNodeBase ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEGeomVirtualNode<RealType> ( tn, no ) {}

  virtual ~CFEVirtualNodeBase ( ) {}

  /** Return how much memory this virtual node needs (hopefully correct)
   */
  int memoryUsed() const { return ( CFEGeomVirtualNode<RealType>::memoryUsed() + _weights.size() * sizeof ( ExtrapolationWeightType ) ); }

  /** Computes normalized gradient of the interface function.
   *
   *  The method uses the pre-computed gradients of the basis
   *  functions on the std-tetra.  An averaging is done over all
   *  tetrahedra that contain this virtual node.
   */
  void determineApproximateNormal ( const CFEStructure<RealType> &s );

  /** Return the i-th weight
   */
  ExtrapolationWeightType weight ( const int i ) const { return _weights[i]; }

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNodeBase ( const CFEVirtualNodeBase &other );
  CFEVirtualNodeBase& operator= ( const CFEVirtualNodeBase& );

protected:
  void determineConstraintsCD ( const CFEGridBase< RealType, CT, NodalCoeffType > &grid );

  //! Check whether we have at least one constraint for the virtual node.
  void postprocessDetermineConstraints ( ) {
    if ( this->_numOfConstraints == 0 ) {
      cerr << this->_global0 << " " << this->_global1 << " " << this->_coord << endl;
      throw aol::Exception ( "No constraints for virtual node found. You're in trouble here", __FILE__, __LINE__ );
    }

#ifdef VERBOSE
    {
      ExtrapolationWeightType sum = aol::ZOTrait<ExtrapolationWeightType>::zero;
      const int _numOfConstraints = this->_constraints.size();
      cerr << aol::color::green << "Number of constraints = " << aol::color::reset << this->_numOfConstraints << endl;
      for ( int i = 0; i < _numOfConstraints; ++i ) {
        cerr << aol::detailedFormat ( _weights[i] ) << " " << this->_constraints[i] << endl;
        sum += this->_weights[i];
      }
      sum -= aol::ZOTrait<ExtrapolationWeightType>::one;
      cerr << "Sum of weights minus 1 = " << sum << endl;
    }
#endif

  }
};


template < typename RealType, tpcfe::ConstraintType CT >
class CFEVirtualNodeScalar : public CFEVirtualNodeBase < RealType, CT, RealType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNodeScalar ( const short local0, const short local1,
                         const int global0, const int global1,
                         const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeBase< RealType, CT, RealType > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNodeScalar ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeBase< RealType, CT, RealType > ( tn, no ) {}

  /** Perform the extrapolation of the given data based on the weights and constraints that
   *  have been determined by determineConstraints
   */
  RealType extrapolate ( const qc::ScalarArray< RealType, qc::QC_3D > &s ) const {
    RealType value = aol::NumberTrait<RealType>::zero;
    for ( int i = 0; i < this->_numOfConstraints; ++i ) {
      value += s.get ( this->constraint ( i ) ) * this->weight ( i );
    }
    return ( value );
  }

  virtual void extrapolate ( const qc::MultiArray< RealType, qc::QC_3D > &data, aol::Vec3< RealType > &dest ) const {
    for ( short c = 0; c < 3; ++c ) {
      dest[c] = extrapolate ( data[c] );
    }
  }

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNodeScalar ( const CFEVirtualNodeScalar &other );
  CFEVirtualNodeScalar& operator= ( const CFEVirtualNodeScalar& );

protected:
  void determineConstraintsLIEHR ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeff, const CFEGridBase<RealType, CT, RealType > &grid );
  void determineConstraintsTPOS ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase<RealType, CT, RealType > &grid );
  void getTPOSLocalMatrix ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFETopoElement &sparseEl, const int parentTetra, aol::Vec<4, int> &localConstraint, aol::Matrix44<RealType> &locMat );
  void computeTPOSLocalWeightsAndConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFETopoElement &sparseEl, const int parentTetra, aol::Vec<4, RealType> &localWeight, aol::Vec<4, int> &localConstraint, bool &inversionReliable );
};


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
class CFEVirtualNode;

template < typename RealType >
class CFEVirtualNode< RealType, CFE_NONE, RealType > : public CFEVirtualNodeScalar < RealType, CFE_NONE > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_NONE > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_NONE > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_NONE, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_NONE, RealType > ( const CFEVirtualNode< RealType, CFE_NONE, RealType > &other );
  CFEVirtualNode< RealType, CFE_NONE, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_NONE, RealType >& );
};


template < typename RealType, typename NodalCoeffType >
class CFEVirtualNode< RealType, CFE_CD, NodalCoeffType > : public CFEVirtualNodeBase < RealType, CFE_CD, NodalCoeffType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeBase < RealType, CFE_CD, NodalCoeffType > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeBase < RealType, CFE_CD, NodalCoeffType >  ( tn, no ) {}

  void determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &/*nodalCoeffs*/, const CFEGridBase< RealType, CFE_CD, NodalCoeffType > &grid );

  /** Perform the extrapolation of the given data based on the weights and constraints that
   *  have been determined by determineConstraints
   */
  RealType extrapolate ( const qc::ScalarArray< RealType, qc::QC_3D > &s ) const {
    RealType value = aol::NumberTrait<RealType>::zero;
    for ( int i = 0; i < this->_numOfConstraints; ++i ) {
      value += s.get ( this->constraint ( i ) ) * this->weight ( i );
    }
    return ( value );
  }

  virtual void extrapolate ( const qc::MultiArray< RealType, qc::QC_3D > &data, aol::Vec3< RealType > &dest ) const {
    for ( short c = 0; c < 3; ++c ) {
      dest[c] = extrapolate ( data[c] );
    }
  }

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_CD, NodalCoeffType > ( const CFEVirtualNode< RealType, CFE_CD, NodalCoeffType > &other );
  CFEVirtualNode< RealType, CFE_CD, NodalCoeffType >& operator= ( const CFEVirtualNode< RealType, CFE_CD, NodalCoeffType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_LIEHR, RealType > : public CFEVirtualNodeScalar < RealType, CFE_LIEHR > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_LIEHR > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_LIEHR > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_LIEHR, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_LIEHR, RealType > ( const CFEVirtualNode< RealType, CFE_LIEHR, RealType > &other );
  CFEVirtualNode< RealType, CFE_LIEHR, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_LIEHR, RealType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_TPOS, RealType > : public CFEVirtualNodeScalar < RealType, CFE_TPOS > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_TPOS > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_TPOS > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase<RealType, CFE_TPOS, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_TPOS, RealType > ( const CFEVirtualNode< RealType, CFE_TPOS, RealType > &other );
  CFEVirtualNode< RealType, CFE_TPOS, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_TPOS, RealType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_TP, RealType > : public CFEVirtualNodeScalar < RealType, CFE_TP > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_TP > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_TP > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_TP, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_TP, RealType > ( const CFEVirtualNode< RealType, CFE_TP, RealType > &other );
  CFEVirtualNode< RealType, CFE_TP, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_TP, RealType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_DOMAIN, RealType > : public CFEVirtualNodeScalar < RealType, CFE_DOMAIN > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_DOMAIN > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_DOMAIN > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_DOMAIN, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_DOMAIN, RealType > ( const CFEVirtualNode< RealType, CFE_DOMAIN, RealType > &other );
  CFEVirtualNode< RealType, CFE_DOMAIN, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_DOMAIN, RealType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_CDWI_LIEHR, RealType > : public CFEVirtualNodeScalar < RealType, CFE_CDWI_LIEHR > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_CDWI_LIEHR > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_CDWI_LIEHR > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_CDWI_LIEHR, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_CDWI_LIEHR, RealType > ( const CFEVirtualNode< RealType, CFE_CDWI_LIEHR, RealType > &other );
  CFEVirtualNode< RealType, CFE_CDWI_LIEHR, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_CDWI_LIEHR, RealType >& );
};


template < typename RealType >
class CFEVirtualNode< RealType, CFE_CDWI_TPOS, RealType > : public CFEVirtualNodeScalar < RealType, CFE_CDWI_TPOS > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeScalar < RealType, CFE_CDWI_TPOS > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeScalar< RealType, CFE_CDWI_TPOS > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<RealType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_CDWI_TPOS, RealType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_CDWI_TPOS, RealType > ( const CFEVirtualNode< RealType, CFE_CDWI_TPOS, RealType > &other );
  CFEVirtualNode< RealType, CFE_CDWI_TPOS, RealType >& operator= ( const CFEVirtualNode< RealType, CFE_CDWI_TPOS, RealType >& );
};


template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
class CFEVirtualNodeVector: public CFEVirtualNodeBase < RealType, CT, NodalCoeffType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;
  typedef typename CFEVirtualNodeBase< RealType, CT, NodalCoeffType >::ExtrapolationWeightType ExtrapolationWeightType;

  CFEVirtualNodeVector ( const short local0, const short local1, const int global0, const int global1,
                         const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeBase < RealType, CT, NodalCoeffType > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNodeVector ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeBase < RealType, CT, NodalCoeffType > ( tn, no ) {}

  virtual void extrapolate ( const qc::MultiArray< RealType, qc::QC_3D > &data, aol::Vec3< RealType > &dest ) const;

  void determineConstraintsTPOSELAST ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CT, NodalCoeffType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNodeVector< RealType, CT, NodalCoeffType > ( const CFEVirtualNodeVector< RealType, CT, NodalCoeffType > &other );
  CFEVirtualNodeVector< RealType, CT, NodalCoeffType >& operator= ( const CFEVirtualNodeVector< RealType, CT, NodalCoeffType >& );
};


template < typename RealType, typename NodalCoeffType >
class CFEVirtualNode< RealType, CFE_TPOSELAST, NodalCoeffType > : public CFEVirtualNodeVector < RealType, CFE_TPOSELAST, NodalCoeffType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;
  typedef typename CFEVirtualNodeBase< RealType, CFE_TPOSELAST, NodalCoeffType >::ExtrapolationWeightType ExtrapolationWeightType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeVector < RealType, CFE_TPOSELAST, NodalCoeffType > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeVector < RealType, CFE_TPOSELAST, NodalCoeffType > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_TPOSELAST, NodalCoeffType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_TPOSELAST, NodalCoeffType > ( const CFEVirtualNode< RealType, CFE_TPOSELAST, NodalCoeffType > &other );
  CFEVirtualNode< RealType, CFE_TPOSELAST, NodalCoeffType >& operator= ( const CFEVirtualNode< RealType, CFE_TPOSELAST, NodalCoeffType >& );
};


template < typename RealType, typename NodalCoeffType >
class CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType > : public CFEVirtualNodeVector < RealType, CFE_CDWI_ELAST, NodalCoeffType > {
public:
  typedef typename CFEGeomVirtualNode< RealType >::IndexType IndexType;
  typedef typename CFEGeomVirtualNode< RealType >::PTMapType PTMapType;
  typedef typename CFEVirtualNodeBase< RealType, CFE_CDWI_ELAST, NodalCoeffType >::ExtrapolationWeightType ExtrapolationWeightType;

  CFEVirtualNode ( const short local0, const short local1, const int global0, const int global1,
                   const tpcfe::CFEElement<RealType> &inEl, const unsigned int ParentTetra = tpcfe::NODE_NOT_VIRTUAL )
      : CFEVirtualNodeVector < RealType, CFE_CDWI_ELAST, NodalCoeffType > ( local0, local1, global0, global1, inEl, ParentTetra ) {}

  CFEVirtualNode ( const aol::Vec3<RealType> &tn, const aol::Vec3<RealType> &no ) : CFEVirtualNodeVector < RealType, CFE_CDWI_ELAST, NodalCoeffType > ( tn, no ) {}

  void determineConstraints ( const qc::AArray<NodalCoeffType, qc::QC_3D> &nodalCoeffs, const CFEGridBase< RealType, CFE_CDWI_ELAST, NodalCoeffType > &grid );

private:
  /** Copy constructor and assignment operator not implemented */
  CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType > ( const CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType > &other );
  CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType >& operator= ( const CFEVirtualNode< RealType, CFE_CDWI_ELAST, NodalCoeffType >& );
};


/** Output this virtual node to the given stream
 */
template < typename RealType, tpcfe::ConstraintType CT, typename NodalCoeffType >
inline ostream& operator<< ( ostream &out, const CFEVirtualNode< RealType, CT, NodalCoeffType > &node ) {
  out << "CFEVirtualNode" << endl;
  out << "   coord " << node._coord << endl;
  out << "   normal " << node._normal << endl;
  for ( int i = 0; i < node._numOfConstraints; ++i ) {
    out << "   constraint " << node.constraint ( i ) << " weight " << node.weight ( i ) << endl;
  }
  return out;
}
}
#endif
