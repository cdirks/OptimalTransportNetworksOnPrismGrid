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

#ifndef __ELASTICOPS_H
#define __ELASTICOPS_H

#include <bemOps.h>
#include <elasticFunSol.h>

namespace bm {

//! Specifies certain properties of local solution operators
//! doublelayer: is a double layer operator
//! jumping:     ansatz space is continuous with some jumps
template <class LocalOperatorType> class OperatorTrait {};

template <class ParticleType> class OperatorTrait<ElasticSolution<ParticleType> > {
public:
 typedef IntegralOperator<ElasticSolution<ParticleType> > Operator;
 const static bool nodecentered = false;
 const static bool doublelayer  = false;
 const static bool jumping      = false;
};

template <class ParticleType> class OperatorTrait<ElasticSolutionDiff<ParticleType> > {
public:
 typedef IntegralOperator<ElasticSolutionDiff<ParticleType> > Operator;
 const static bool nodecentered = false;
 const static bool doublelayer  = true;
 const static bool jumping      = false;
};

template <class ParticleType> class OperatorTrait<LinElasticSolution<ParticleType> > {
public:
 typedef LinIntegralOperator<LinElasticSolution<ParticleType> > Operator;
 const static bool nodecentered = true;
 const static bool doublelayer  = false;
 const static bool jumping      = false;
};

template <class ParticleType> class OperatorTrait<LinElasticSolutionDiff<ParticleType> > {
public:
 typedef LinIntegralOperator<LinElasticSolutionDiff<ParticleType> > Operator;
 const static bool nodecentered = true;
 const static bool doublelayer  = true;
 const static bool jumping      = false;
};

template <class ParticleType> class OperatorTrait<LinJumpElasticSolution<ParticleType> > {
public:
 typedef LinJumpIntegralOperator<LinJumpElasticSolution<ParticleType> > Operator;
 const static bool doublelayer = false;
 const static bool jumping     = true;
};


//! Utility function for ElasticOperator
template <class ParticleType, template <class> class LocalOperatorType, class GreenType >
void setBlockInMatrix ( const bm::Boundary<ParticleType>& bnd,
                        const bm::GenElasticGreen<typename ParticleType::DataType, GreenType >& green,
                        aol::FullMatrix<typename ParticleType::DataType>& matrix, int row, int col, int rskip, int cskip,
                        typename ParticleType::DataType factor = 1, bool subtractconst = false ) {

  for ( int j = 0; j < 2; ++j ) {
    for ( int i = 0; i < 2; ++i ) {
      LocalOperatorType<ParticleType> loc ( green.asImp() , i, j );
      typename bm::OperatorTrait<LocalOperatorType<ParticleType> >::Operator op ( loc, bnd );
      int n = op.getNumRows (), m = op.getNumCols ();
      // Subtract free term coefficients from double layer operator
      if ( subtractconst ) {
        aol::Vector<typename ParticleType::DataType> cons ( m ), ft ( n );
        for ( int k = 0; k < n; ++k ) cons [k] = 1;
        op.apply ( cons, ft );
        // factor > 0 here means doublelayer (-1) and exterior problem (-1)
        for ( int k = 0; k < n; ++k ) op.add ( k, k, - ft [k] + ( ( factor > 0 && i == j ) ? 1 : 0 ) );
      }
      if ( factor != 1 ) op *= factor;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int k = 0; k < n; ++k )
        for ( int l = 0; l < m; ++l )
          matrix.ref ( row + i * rskip + k, col + j * cskip + l ) = op.ref ( k, l );
    }
  }
}

//! Utility function for ElasticEigenstress
template <class ParticleType>
void setRHS ( const bm::Boundary<ParticleType>& bnd,
              const bm::ElasticGreen<typename ParticleType::DataType>& green,
              const aol::Matrix22<typename ParticleType::DataType>& barepsilon,
              aol::Vector<typename ParticleType::DataType>& rhs, int row, bool nodecentered,
              typename ParticleType::DataType factor = 1 ) {
  int n = bnd.getNumberOfSegments ();
  aol::Matrix22<typename ParticleType::DataType> cbarepsilon; green.C.apply ( barepsilon, cbarepsilon );

  typename bm::Boundary<ParticleType>::const_iterator partit = bnd.begin (), partend = bnd.end ();
  for ( int i = 0; partit != partend; ++partit ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment ();
    typename ParticleType::ConstSegmentIteratorType end = partit->endSegment ();
    typename ParticleType::ConstSegmentIteratorType pit = end; --pit;

    for ( ; it != end; pit = it, ++it, ++i ) {
      // Normal inwards for matrix phase
      aol::Vec2<typename ParticleType::DataType> normal = nodecentered ? cornerNormal ( *pit, *it ) : it->getNormal ();
      aol::Vec2<typename ParticleType::DataType> e = cbarepsilon * normal;
      rhs [row + i] = e [0] * factor;
      rhs [row + n + i] = e [1] * factor;
    }
  }
}

//! Utility function for ElasticEigenstresses
template <class ParticleType>
void setRHS ( const bm::Boundary<ParticleType>& bnd,
              const bm::ElasticGreen<typename ParticleType::DataType>& green,
              const aol::Matrix22<typename ParticleType::DataType>& barepsilon,
              aol::MultiVector<typename ParticleType::DataType>& rhs, int row, bool nodecentered,
              typename ParticleType::DataType factor = 1 ) {
  aol::Matrix22<typename ParticleType::DataType> cbarepsilon; green.C.apply ( barepsilon, cbarepsilon );

  typename bm::Boundary<ParticleType>::const_iterator partit = bnd.begin (), partend = bnd.end ();
  for ( int i = 0; partit != partend; ++partit ) {
    typename ParticleType::ConstSegmentIteratorType it = partit->beginSegment ();
    typename ParticleType::ConstSegmentIteratorType end = partit->endSegment ();
    typename ParticleType::ConstSegmentIteratorType pit = end; --pit;

    for ( ; it != end; pit = it, ++it, ++i ) {
      // Normal inwards for matrix phase
      aol::Vec2<typename ParticleType::DataType> normal = nodecentered ? cornerNormal ( *pit, *it ) : it->getNormal ();
      aol::Vec2<typename ParticleType::DataType> e = cbarepsilon * normal;
      rhs [0] [row + i] = e [0] * factor;
      rhs [1] [row + i] = e [1] * factor;
    }
  }
}

//! Integral operator for vector valued problems (i.e. the fundamental solution is matrix valued)
/** The full matrix is assembled from 2x2 large blocks, each for one component <br>
 *  \note This will give you             \f$ + U[t] \f$ or
 *                                       \f$ - V[u] + F u \f$ for the inner problem.
 *        For the outer one you will get \f$ - U[t] \f$ or
 *                                       \f$ + V[u] + ( I - F ) u \f$.
 */
template <class ParticleType, template <class> class LocalOperatorType, class GreenType = ElasticGreen<typename ParticleType::DataType> >
class ElasticOperator : public aol::FullMatrix<typename ParticleType::DataType> {
public:
  ElasticOperator ( const bm::Boundary<ParticleType>& bnd,
                    const bm::GenElasticGreen< typename ParticleType::DataType, GreenType >& green,
                    bool interior = false )
      : aol::FullMatrix<typename ParticleType::DataType> ( ( bnd.getNumberOfSegments ()
                                                             - ( OperatorTrait<LocalOperatorType<ParticleType> >::jumping ? bnd.size () : 0 ) ) * 2,
                                                             bnd.getNumberOfSegments () * 2 )  {
    // Multiply by -1 if doublelayer operator, multiply by -1 again if solving exterior problem
    typename ParticleType::DataType factor = ( OperatorTrait<LocalOperatorType<ParticleType> >::doublelayer ? -1 : 1 ) * ( !interior ? -1 : 1 );

    setBlockInMatrix<ParticleType, LocalOperatorType> ( bnd, green, *this, 0, 0,
                                                        bnd.getNumberOfSegments () - ( OperatorTrait<LocalOperatorType<ParticleType> >::jumping ? bnd.size () : 0 ),
                                                        bnd.getNumberOfSegments (),
                                                        factor, OperatorTrait<LocalOperatorType<ParticleType> >::doublelayer );
  }

  //! Coordinate transformation
  void transform ( aol::Matrix22<typename ParticleType::DataType> trf ) {
    int n = this->getNumRows () / 2, m = this->getNumCols () / 2;

    for ( int i = 0; i < n; ++i )
      for ( int j = 0; j < m; ++j ) {
        aol::Matrix22<typename ParticleType::DataType> integral ( this->get ( i,     j ), this->get ( n + i,     j ),   // Indices of matrix kernel
                                                                  this->get ( i, m + j ), this->get ( n + i, m + j ) ); // are transposed
        aol::Matrix22<typename ParticleType::DataType> trfintegral ( 0, 0, 0, 0 );

        for ( int a = 0; a < 2; ++a )
          for ( int b = 0; b < 2; ++b )
            for ( int c = 0; c < 2; ++c )
              for ( int d = 0; d < 2; ++d )
                trfintegral [a][b] += trf [a][c] * trf [b][d] * integral [c][d];

        this->set ( i,         j, trfintegral [0][0] ); this->set ( n + i,     j, trfintegral [0][1] );
        this->set (     i, m + j, trfintegral [1][0] ); this->set ( n + i, m + j, trfintegral [1][1] );
      }
  }

  //! Rotation of coordinate system
  void rotate ( typename ParticleType::DataType angle ) {
    aol::Matrix22<typename ParticleType::DataType> rot (   cos ( angle ), sin ( angle ),
                                                         - sin ( angle ), cos ( angle ) );
    transform ( rot );
  }
};

//! Compute vector containing jump of tractions on the interface due to elastic misfit
template <class ParticleType>
class ElasticEigenstress : public aol::Vector<typename ParticleType::DataType> {
public:
  //! Constructor for case without misfit
  ElasticEigenstress ( const bm::Boundary<ParticleType>& bnd, const bm::ElasticGreen<typename ParticleType::DataType>& )
      : aol::Vector<typename ParticleType::DataType> ( 2 * bnd.getNumberOfSegments () ) {
    this->setZero ();
  };
  ElasticEigenstress ( const bm::Boundary<ParticleType>& bnd, const bm::ElasticGreen<typename ParticleType::DataType>& green,
                       const aol::Matrix22<typename ParticleType::DataType>& barepsilon, bool nodecentered, bool interior = false )
      : aol::Vector<typename ParticleType::DataType> ( 2 * bnd.getNumberOfSegments () ) {
    setRHS ( bnd, green, barepsilon, *this, 0, nodecentered, interior ? -1 : 1 );
  };
};

//! Compute multivector containing jump of tractions on the interface due to elastic misfit
template <class ParticleType>
class ElasticEigenstresses : public aol::MultiVector<typename ParticleType::DataType> {
public:
  ElasticEigenstresses ( const bm::Boundary<ParticleType>& bnd, const bm::ElasticGreen<typename ParticleType::DataType>& green,
                         const aol::Matrix22<typename ParticleType::DataType>& barepsilon, bool nodecentered, bool interior = false )
      : aol::MultiVector<typename ParticleType::DataType> ( 2, bnd.getNumberOfSegments () ) {
    setRHS ( bnd, green, barepsilon, *this, 0, nodecentered, interior ? -1 : 1 );
  };
};

template <class ParticleType>
class DisplacementNormalization : public aol::FullMatrix<typename ParticleType::DataType> {
public:
  DisplacementNormalization ( const bm::Boundary<ParticleType>& bnd, bool nodecentered )
      : aol::FullMatrix<typename ParticleType::DataType> ( 2, bnd.getNumberOfSegments () * 2 ) {
    int n = bnd.getNumberOfSegments ();
    aol::Vector<double> lens ( n );
    bnd.getLengths ( lens, nodecentered );
    aol::FullMatrix<double> lensLine ( lens, true );
    this->setBlock ( 0, 0, lensLine );
    this->setBlock ( 1, n, lensLine );
  }
};

}

#endif
