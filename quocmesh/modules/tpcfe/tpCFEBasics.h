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

#ifndef __TPCFEBASICS_H
#define __TPCFEBASICS_H

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>
#include <indexMapper.h>
#include <gridBase.h>

namespace tpcfe {

/** Enum defining how virtual nodes are constrained, thus describing which type of CFE are used.
 */
enum ConstraintType { CFE_NONE,       //!< CFE are used for finding the virtual nodes only, no extrapolation is used.
                      CFE_CD,         //!< CFE for complicated domains in scalar and vector-valued problems. Levelset defines domain of computation. No dofs outside, virtual nodes are non-dofs but slave nodes. See Liehr et al., Comput Vis Sci 12(4):171-188, 2009, doi: 10.1007/s00791-008-0093-1
                      CFE_DOMAIN,     //!< CFE used to describe a complicated domain, using non-hexahedral mesh with possibly bad tetrahedra. Virtual nodes have a global dof and their weight is 1.
                      CFE_LIEHR,      //!< CFE for discontinuous coefficients in scalar problems with non-optimal approximation order (according to the Diplom thesis of F. Liehr).
                      CFE_TPOS,       //!< CFE for discontinuous coefficients in scalar problems. See Preusser et al., SIAM J Sci Comput (5):2115-2143, 2011, doi: 10.1137/100791750
                      CFE_TP,         //!< CFE for discontinuous coefficients in scalar problems (not recommended).
                      CFE_CDWI_LIEHR, //!< CFE for complicated domains and discontinuous coefficients using the LIEHR scheme
                      CFE_CDWI_TPOS,  //!< CFE for complicated domains and discontinuous coefficients using the TPOS scheme
                      CFE_TPOSELAST,  //!< CFE for discontinuous coefficients in vector-valued problems (linearized elasticity), extension of CFE_TPOS. See Preusser et al., SIAM J Sci Comput (5):2115-2143, 2011, doi: 10.1137/100791750
                      CFE_CDWI_ELAST  //!< CFE for complicated domains and discontinuous coefficients using the TPOSELAST scheme (linearized elasticity; vector-valued problem)
};

static const std::string ConstraintTypeNames[10] = { "CFE_NONE",
                                                     "CFE_CD",
                                                     "CFE_LIEHR",
                                                     "CFE_TPOS",
                                                     "CFE_TP",
                                                     "CFE_DOMAIN",
                                                     "CFE_CDWI_LIEHR",
                                                     "CFE_CDWI_TPOS",
                                                     "CFE_TPOSELAST",
                                                     "CFE_CDWI_ELAST" };

/** collection of constants:
 *  largest structure ID, as this index, the complicated domain levelset is stored,
 *  a node has 15 AFE neighbors: 14 proper neighbors and itself,
 *  the node itself has neighbor index 7
 *  node on edge is given by two local indices; second index = 11 means regular node at first index, else virtual node on edge between the two
 *  number of edges of the standard tetrahedra in a cube
 */
enum { MAX_STRUCT_ID = 32,
       NUM_NEIGHBORS = 15,
       NEIGHBOR_SELF = 7,
       // NUM_NEIGHBORS_JC = 89,
       // NEIGHBOR_SELF_JC = 44,
       NODE_NOT_VIRTUAL = 11,
       NUM_TETRA_EDGES_IN_CUBE = 19 };


/** Type of a CFE element: index of structure interfacing element and sign pattern of element vertices
 *  \author Schwen
 */
class CFEType {
public:
  typedef signed char StructNoType;
  typedef unsigned char SignPatternType;

  StructNoType   _structureNo;  //!< index of the structure by which the element is interfaced
  SignPatternType _pureType;    //!< sign pattern of levelset function on the edges of the element

  CFEType ( ) : _structureNo ( -1 ), _pureType ( 0 ) {}

  CFEType ( const StructNoType StructureNo, const SignPatternType PureType ) : _structureNo ( StructureNo ), _pureType ( PureType ) {}

  // default copy constructor and assignment operator do the right thing

  //! an element is not interfaced if all edges have the same sign
  inline bool representsInterfaced ( ) const {
    return ( ( _pureType != 0x00 ) && ( _pureType != 0xFF ) );
  }

  inline signed char signForNonInterfacedType ( ) const {
#ifdef DEBUG
    if ( _pureType == 0x00 ) {
      return ( +1 );
    } else if ( _pureType == 0xFF ) {
      return ( -1 );
    } else {
      throw aol::Exception ( "Interfaced Elements do not have a unique sign", __FILE__, __LINE__ );
    }
#else
    return ( ( _pureType == 0x00 ) ? +1 : -1 );
#endif
  }
};



/** Abstract basis class for levelset functions defining interfaces in the CFE context (interfaces between regions with different material coefficient or between material and void)
 *  \author Schwen
 */
template< typename RealType >
class CFEStructure {
public:
  CFEStructure ( ) { }

  virtual ~CFEStructure ( ) { }

  //! Return value of the levelset function at a regular node position. Note: in general, a CFEStructure does not have a vector, so getValue ( int ) does not make sense.
  virtual RealType getValue ( const qc::CoordType& /*pos*/ ) const = 0;

  //! Compute the cut relation between the two nodes. This only makes sense for neighbors connected by regular tetra edges, but can always be called.
  virtual RealType getCutRelationBetween ( const qc::CoordType& /*pos0*/, const qc::CoordType &/*pos1*/ ) const = 0;

  virtual aol::Vec3<RealType> getInterfaceNormal ( const aol::Vec3<RealType> & /*pos*/ ) const = 0;

protected:
  inline bool differentSign ( const RealType v0, const RealType v1 ) const {
    return ( ( ( v0 < 0 ) && ( v1 > 0 ) ) || ( ( v0 > 0 ) && ( v1 < 0 ) ) );
  }

private:
  //! copy constructor should not be called
  CFEStructure ( const CFEStructure<RealType>& );

  //! assignment operator should not be called
  CFEStructure<RealType>& operator= ( const CFEStructure<RealType>& );
};


/** Affine interpolation between regular nodes on 3D data
 */
template< typename RealType >
class CFEStructureAffine : public CFEStructure<RealType> {
private:
  qc::ScalarArray<RealType, qc::QC_3D> _image;

public:
  CFEStructureAffine ( ) : _image () { }

  void setImageFrom ( const qc::ScalarArray<RealType, qc::QC_3D> &arr, const RealType adjust = aol::NumberTrait<RealType>::zero ) {
    _image.reallocate ( arr );
    _image = arr;
    int adjustCounter = 0;
    for ( int i = 0; i < _image.size(); ++i ) {
      if ( fabs ( _image[i] ) < adjust ) {
#ifdef VERBOSE
        qc::OTFILexMapper<qc::QC_3D> imap ( arr );
        cerr << "Levelset value " << _image[i] << " at position " << imap.splitGlobalIndex ( i ) << ", shifted to " << ( _image[i] < aol::NumberTrait<RealType>::zero  ?  -adjust  :  adjust ) << endl;
#endif
        _image[i] = ( _image[i] < aol::NumberTrait<RealType>::zero  ?  -adjust  :  adjust );
        ++adjustCounter;
      }
    }
    if ( adjustCounter != 0 ) {
      cerr << "Levelset function adjusted (shifted away from zero) at " << adjustCounter << " nodes" << endl;
    }
  }

public:

  RealType getValue ( const qc::CoordType& pos ) const {
    return ( _image.get ( pos ) );
  }


  RealType getCutRelationBetween ( const qc::CoordType& pos0, const qc::CoordType &pos1 ) const {
    const RealType val0 = getValue ( pos0 ), val1 = getValue ( pos1 );
    return ( val0 / ( val0 - val1 ) );
  }

  virtual aol::Vec3<RealType> getInterfaceNormal ( const aol::Vec3<RealType> &/*pos*/ ) const {
    throw aol::UnimplementedCodeException ( "tpcfe::CFEStructureAffine::getInterfaceNormal not implemented yet.", __FILE__, __LINE__ );
    aol::Vec3<RealType> n;
    return ( n );
  }
};


/** Trilinear interpolation between regular nodes on 3D data
 */
template< typename RealType >
class CFEStructureTrilin : public CFEStructure<RealType> {
private:
  qc::ScalarArray<RealType, qc::QC_3D> _image;
  static const int _numSteps = 15;
  static const RealType _threshold;

public:
  CFEStructureTrilin ( ) : _image () { }

  void setImageFrom ( const qc::ScalarArray<RealType, qc::QC_3D> &arr, const RealType adjust = aol::NumberTrait<RealType>::zero ) {
    _image.reallocate ( arr );
    _image = arr;
    int adjustCounter = 0;
    for ( int i = 0; i < _image.size(); ++i ) {
      if ( fabs ( _image[i] ) < adjust ) {
#ifdef VERBOSE
        qc::OTFILexMapper<qc::QC_3D> imap ( arr );
        cerr << "Levelset value " << _image[i] << " at position " << imap.splitGlobalIndex ( i ) << ", shifted to " << ( _image[i] < aol::NumberTrait<RealType>::zero  ?  -adjust  :  adjust ) << endl;
#endif
        _image[i] = ( _image[i] < aol::NumberTrait<RealType>::zero  ?  -adjust  :  adjust );
        ++adjustCounter;
      }
    }
    if ( adjustCounter != 0 ) {
      cerr << "Levelset function adjusted (shifted away from zero) at " << adjustCounter << " nodes" << endl;
    }
  }

public:

  RealType getValue ( const qc::CoordType& pos ) const {
    return ( _image.get ( pos ) );
  }

  RealType getCutRelationBetween ( const qc::CoordType& pos0, const qc::CoordType &pos1 ) const {
    // via secant search
    RealType b0 = aol::NumberTrait<RealType>::zero, b1 = aol::NumberTrait<RealType>::one;
    const aol::Vec3<RealType> psr0 ( pos0 ), psr1 ( pos1 );
    RealType v0 = _image.interpolate ( psr0 ), v1 = _image.interpolate ( psr1 );

    if ( differentSign ( v0, v1 ) ) {
      int i = 0;
      do {
        const RealType b2 = ( b0 * v1 - b1 * v0 ) / ( v1 - v0 );
        b0 = b1;
        v0 = v1;
        b1 = b2;
        {
          const aol::Vec3<RealType> p1 = psr0 + b1 * ( psr1 - psr0 );
          v1 = _image.interpolate ( p1 );
        }
        ++i;
      } while ( ( i < _numSteps ) && ( fabs ( v1 ) > _threshold ) );

      return ( b1 );
    } else { // no zero on edge
      return ( aol::NumberTrait<RealType>::NaN );
    }

  }

  aol::Vec3<RealType> getInterfaceNormal ( const aol::Vec3<RealType> &/*pos*/ ) const {
    throw aol::UnimplementedCodeException ( "tpcfe::CFEStructureAffine::getInterfaceNormal not implemented yet.", __FILE__, __LINE__ );
    aol::Vec3<RealType> n;
    return ( n );
  }

};


/** check whether value is smaller than threshold and print OK or FAILED and throw exception
 *  \author Schwen
 */
template< typename RealType >
void smallOrDie ( const RealType value, const RealType threshold, const char* message, const char* file, const int line ) {
  if ( fabs ( value ) > threshold ) {
    cerr << message << " " << value << " FAILED." << endl;
    throw ( aol::Exception ( message, file, line ) );
  } else {
    cerr << message << " " << value << " OK." << endl;
  }
}


/** Material parameters for isotropic linear elasticity returning local E, nu, lambda, mu, anisotropic tensor
 *  \todo think about whether this averaging makes sense
 *  Can be implemented for more general local anisotropy.
 */
template < typename RealType >
class IsotropicElasticityCoefficient {
protected:
  RealType _E;
  RealType _nu;

public:
  IsotropicElasticityCoefficient ( ) : _E ( 0 ), _nu ( 0 ) { }

  IsotropicElasticityCoefficient ( const RealType E, const RealType nu ) : _E ( E ), _nu ( nu ) { }

  inline RealType getE ( ) const { return ( _E ); }
  inline RealType getNu ( ) const { return ( _nu ); }

  RealType getLambda ( ) const;
  RealType getMu ( ) const;
  void getAnisotropicTensor ( RealType Tensor[3][3][3][3] ) const;

  void print ( ostream &os ) const;

  //! addition of E and nu do not make sense in general
  IsotropicElasticityCoefficient<RealType>& operator+= ( const IsotropicElasticityCoefficient<RealType> &other );
  //! scaling of E and nu does not make sense in general
  IsotropicElasticityCoefficient<RealType>& operator*= ( const RealType value );
  //! scaling of E and nu does not make sense in general
  IsotropicElasticityCoefficient<RealType>& operator/= ( const RealType value );
};

template< typename RealType >
ostream& operator<< ( ostream &os, const IsotropicElasticityCoefficient<RealType> &coeff ) {
  coeff.print ( os );
  return ( os );
}


template < typename RealType >
class VoigtElasticityCoefficient {
protected:
  aol::Mat<6,6,RealType> _VoigtTensor;
  static const short indexM[3][3];

public:
  VoigtElasticityCoefficient ( ) : _VoigtTensor ( ) {
  }

  VoigtElasticityCoefficient ( const aol::Mat<6,6,RealType> &VoigtTensor ) : _VoigtTensor ( VoigtTensor ) { }

  //! write this tensor in fourth-order tensor form to Tensor
  void getAnisotropicTensor ( RealType Tensor[3][3][3][3] ) const;

  //! convert Tensor to Voigt representation, assuming that it satisfies the required symmetries (average values that should be numerically equal)
  void averageFromAnisotropicTensor ( const RealType Tensor[3][3][3][3] );

  //! apply rotation to this tensor
  void rotateTensorBy ( const aol::Matrix33<RealType> rotMat );

  void print ( ostream &os ) const;

  void setFromENu ( const RealType E, const RealType nu );

  VoigtElasticityCoefficient<RealType>& operator+= ( const VoigtElasticityCoefficient<RealType> &other );
  VoigtElasticityCoefficient<RealType>& operator*= ( const RealType value );
  VoigtElasticityCoefficient<RealType>& operator/= ( const RealType value );
};

template< typename RealType >
ostream& operator<< ( ostream &os, const VoigtElasticityCoefficient<RealType> &coeff ) {
  coeff.print ( os );
  return ( os );
}


template < typename RealType, ConstraintType CT >
struct ExtrapolationWeightTypeTrait {
  typedef RealType WeightType;
};

template < typename RealType >
struct ExtrapolationWeightTypeTrait < RealType, CFE_TPOSELAST > {
  typedef aol::Matrix33<RealType> WeightType;
};

template < typename RealType >
struct ExtrapolationWeightTypeTrait < RealType, CFE_CDWI_ELAST > {
  typedef aol::Matrix33<RealType> WeightType;
};

// forward declarations

class CFEGridDefinition;
template < typename RealType, ConstraintType CT, typename NodalCoeffType > class CFEGridBase;
template < typename RealType, ConstraintType CT, typename NodalCoeffType > class CFEGrid;

template < typename RealType > class CFEElement;

class CFECube;
class CFETopoTetra;
template < typename RealType > class CFETetra;

class CFETopoTetraIterator;
template < typename RealType > class CFETetraInElementIterator;

template < typename RealType > class CFELookup;

template < typename RealType > class CFEStructure;

}

#endif
