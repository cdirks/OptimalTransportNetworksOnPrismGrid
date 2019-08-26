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

#ifndef __OPENMESHOPS_H
#define __OPENMESHOPS_H

#include <FEOpInterface.h>

namespace om {

template <typename ConfType>
class NormalProjMassOp : public aol::LumpedMassOpInterface < ConfType,NormalProjMassOp<ConfType> > {
public:
  typedef typename ConfType::RealType RealType;
protected:
  int _i, _j;
public:
  NormalProjMassOp ( const typename ConfType::InitType &Initializer,  int i, int j )
      : aol::LumpedMassOpInterface<ConfType, NormalProjMassOp<ConfType> > ( Initializer, false ),  _i ( i ), _j ( j ) {}

  inline RealType getCoeff ( const typename ConfType::ElementType &El,
                             int, const typename ConfType::DomVecType& ) const {
    aol::Vec3<RealType> n;
    El.normal ( n );
    return n[_i] * n[_j];
  }
};

template <typename ConfType>
class NormalProjection : public aol::Op<aol::MultiVector<typename ConfType::RealType> > {
public:
  typedef typename ConfType::RealType RealType;
protected:
  NormalProjMassOp<ConfType> _massProjOp00;
  NormalProjMassOp<ConfType> _massProjOp01;
  NormalProjMassOp<ConfType> _massProjOp02;
  NormalProjMassOp<ConfType> _massProjOp11;
  NormalProjMassOp<ConfType> _massProjOp12;
  NormalProjMassOp<ConfType> _massProjOp22;
  aol::LumpedMassOp<ConfType> _lumpedMassInv;

public:
  NormalProjection ( const typename ConfType::InitType &Initializer ) :
      _massProjOp00 ( Initializer, 0, 0 ),
      _massProjOp01 ( Initializer, 0, 1 ),
      _massProjOp02 ( Initializer, 0, 2 ),
      _massProjOp11 ( Initializer, 1, 1 ),
      _massProjOp12 ( Initializer, 1, 2 ),
      _massProjOp22 ( Initializer, 2, 2 ),
      _lumpedMassInv ( Initializer, true ) { }

  void applyAdd ( const aol::MultiVector<RealType> &arg, aol::MultiVector<RealType> &dest ) const {
    aol::MultiVector<RealType> tmp ( dest );

    _massProjOp00.apply ( arg[0], tmp[0] );
    _massProjOp01.applyAdd ( arg[1], tmp[0] );
    _massProjOp02.applyAdd ( arg[2], tmp[0] );

    _massProjOp01.apply ( arg[0], tmp[1] );
    _massProjOp11.applyAdd ( arg[1], tmp[1] );
    _massProjOp12.applyAdd ( arg[2], tmp[1] );

    _massProjOp02.apply ( arg[0], tmp[2] );
    _massProjOp12.applyAdd ( arg[1], tmp[2] );
    _massProjOp22.applyAdd ( arg[2], tmp[2] );

    _lumpedMassInv.applyAdd ( tmp[0], dest[0] );
    _lumpedMassInv.applyAdd ( tmp[1], dest[1] );
    _lumpedMassInv.applyAdd ( tmp[2], dest[2] );
  }
};


/*! \brief Special class for variations of TriMeshes
 * Similar to aol::FENonlinVectorDiffOpInterface< ConfType, 3, 3, Imp >
 * but - in order to reduce unnecessary conversions - the DomainType is TriMeshType (instead of MV)
 */
template <typename TriMeshConfType, typename Imp>
class VariationOpInterface : public aol::FEOpInterface<TriMeshConfType, typename TriMeshConfType::InitType, aol::MultiVector<typename TriMeshConfType::RealType> > {

public:
  typedef typename TriMeshConfType::RealType            RealType;
  typedef typename TriMeshConfType::ElementIteratorType IteratorType;
  typedef typename TriMeshConfType::InitType            TriMeshType;

  typedef aol::Mat<3, 3, RealType> NLType;

protected:
  const TriMeshType &_initializer;

public:
  explicit VariationOpInterface ( const TriMeshType &Initializer ) :
      aol::FEOpInterface< TriMeshConfType, TriMeshType, aol::MultiVector<RealType> > ( Initializer ), _initializer ( Initializer ) {
  }

  VariationOpInterface ( const TriMeshConfType & Config, const TriMeshType &Initializer ) :
      aol::FEOpInterface< TriMeshConfType, TriMeshType, aol::MultiVector<RealType> > ( Config ), _initializer ( Initializer ) {
  }

  virtual ~VariationOpInterface( ) {}
  
  void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.setZero();
    applyAdd( Arg, Dest );
  }
  
  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    TriMeshType ArgMesh( _initializer );
    ArgMesh.fromVector( Arg );
    applyAdd( ArgMesh, Dest );
  }
  
  void apply ( const TriMeshType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.setZero();
    applyAdd( Arg, Dest );
  }

  void applyAdd ( const TriMeshType &Arg, aol::MultiVector<RealType> &Dest ) const {

    int numOfVertices = _initializer.getNumberOfNodes();
    if ( Arg.getNumberOfNodes() != numOfVertices || Dest.numComponents() != 3 || Dest[0].size() != numOfVertices )
      throw aol::Exception ( "om::VariationOpInterface::applyAdd(): Mismatching input.", __FILE__, __LINE__ );

    NLType *nl_cache = new NLType[ this->getConfigurator().maxNumQuadPoints() ];

    const typename IteratorType::EndType end_it = this->getConfigurator().end();

    for ( IteratorType it = this->getConfigurator().begin(); it != end_it; ++it ) {
      const int numLocalDofs = this->getConfigurator().getNumLocalDofs ( *it );

      const typename TriMeshConfType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      //  Actually, a basis function $\phi_i$ located at the i-th node is a vector-valued function $(\phi_i^j)_{j=1,2,3}$,
      //  defined by $\phi_i^j(X_l) = \delta_{il}e_j$ and the linear interpolation between the nodes $(X_l)_l$.
      //  The Jacobi matrix $D\phi_i^j \in \R^{3,3}$ is zero besides the j-th row, which is $\nabla\tilde\phi_i$,
      //  where $\tilde\phi_i$ is the piecewise linear, real-valued function defined by $\tilde\phi_i(X_l) = \delta_{il}$.
      //  This means the j-th row of $D\phi_i^j$ does not depend on j!
      //  Hence we make use of this in the following by calculting the last term of the following expression:
      //! $A : D\phi_i^j = A_{j,} \cdot \nabla\tilde\phi_i = (A\nabla\tilde\phi_i)_j$.
      //  Note, that the evaluation of $D\phi_i^j$ will return the vector $\nabla\tilde\phi_i \in \R^3$ !!
      
      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( Arg, *it, q, bfs.getRefCoord ( q ), nl_cache[q] );

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {
        aol::Vec<3, RealType> aux;
        aux.setZero();

        for ( int q = 0; q < numQuadPoints; ++q ) {
          NLType nl;
          typename TriMeshConfType::VecType grad, tmp;

          grad = bfs.evaluateGradient ( dof, q );
          nl = nl_cache[q];
          nl *= bfs.getWeight ( q );
          nl.mult ( grad, tmp );
          aux += tmp;
        }

        aux *= this->getConfigurator().vol ( *it );

        for ( int d = 0; d < 3; d++ )
          Dest[d][ this->getConfigurator().localToGlobal ( *it, dof ) ] += aux[d];
      }
    }
    delete[] nl_cache;
  }


  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const TriMeshType &Arg,
                         const typename TriMeshConfType::ElementType &El,
                         int QuadPoint,
                         const typename TriMeshConfType::DomVecType &RefCoord,
                         NLType &NL ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().getNonlinearity ( Arg, El, QuadPoint, RefCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};

} // namespace om

#endif // __OPENMESH_OPS_H
