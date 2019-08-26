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

#ifndef __OSAPPROXIMATOR_H
#define __OSAPPROXIMATOR_H

#include "osTestfunction.h"

namespace tpcfe {

//! Abstract base class for FE approximation of functions.
//  Note that not all child classes require all function parameters to be set correctly.
template< class GridType, class RealType >
class Approximator {
protected:
  const GridType &_grid;
  const qc::ScalarArray<RealType, qc::QC_3D> &_data;

public:
  Approximator ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &data ) : _grid ( grid ), _data ( data ) {}

  virtual ~Approximator() {}
public:
  virtual void getApproximateVertexFunctionValues ( const tpcfe::CFEElement<RealType>         &/*el*/,
                                                    const tpcfe::CFETetra<RealType>           &/*tetra*/,
                                                    const std::vector< aol::Vec3<RealType> >  &/*vertex*/,
                                                    aol::Vec<4, RealType>                     &/*vertexValue*/ ) const = 0;

  virtual RealType getApproximateFunctionValue ( const std::vector< aol::Vec3<RealType> >     &/*vertex*/,
                                                 const aol::Vec<4, RealType>                  &/*vertexValue*/,
                                                 const aol::Vec<4, RealType>                  &/*baryCoords*/,
                                                 const aol::Vec3<RealType>                    &/*globalCoords*/ ) const = 0;

  virtual aol::Vec3<RealType> getApproximateGradientValue ( const std::vector< aol::Vec3<RealType> >     &/*vertex*/,
                                                            const aol::Vec<4, RealType>                  &/*vertexValue*/,
                                                            const aol::Vec3<RealType>                    &/*globalCoords*/ ) const = 0;

};


//! CFE extrapolation using correct extrapolation weights for vertices. Interpolation inside tetrahedra is done via barycentric coordinates.
template< class GridType >
class CFEApproximator : public Approximator<GridType, typename GridType::RealType> {
  typedef typename GridType::RealType RealType;
public:
  CFEApproximator ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &data ) : Approximator<GridType, typename GridType::RealType> ( grid, data ) {}

  virtual ~CFEApproximator() {}

public:
  virtual void getApproximateVertexFunctionValues ( const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETetra<RealType> &tetra, const std::vector< aol::Vec3<RealType> > &/*vertex*/, aol::Vec<4, RealType> &vertexValue ) const {
    for ( int i = 0; i < 4; i++ ) {
      int gIdx0, gIdx1;

      gIdx0 = el.globalIndex( tetra ( i,0 ) );
      if ( tetra ( i, 1 ) == 11 ) {
        gIdx1 = -1;
      } else {
        gIdx1 = el.globalIndex( tetra ( i,1 ) );
      }

      if ( gIdx1 == -1 ) { // non-interpolated node
        vertexValue[i] = this->_data.get ( gIdx0 );
      } else {           // interpolated node
        const typename GridType::VNType &vn = this->_grid.getVirtualNodeRef ( gIdx0, gIdx1 );
        if(vn._isDirichlet){
          vertexValue[i] = vn.getDirichlet();
        }
        else{
        vertexValue[i] = vn.extrapolate ( this->_data ); // this does not work. Maybe the constraints are set incorrectly?
        }
      }
    }
  }


  virtual RealType getApproximateFunctionValue ( const std::vector< aol::Vec3<RealType> > &/*vertex*/,
                                                 const aol::Vec<4, RealType> &vertexValue,
                                                 const aol::Vec<4, RealType> &baryCoords,
                                                 const aol::Vec3<RealType> &/*globalCoords*/ ) const {

    return ( baryCoords[0] * vertexValue[0] + baryCoords[1] * vertexValue[1] + baryCoords[2] * vertexValue[2] + baryCoords[3] * vertexValue[3] );
  }

  virtual aol::Vec3<RealType> getApproximateGradientValue ( const std::vector< aol::Vec3<RealType> >     &vertex,
                                                            const aol::Vec<4, RealType>                  &vertexValue,
                                                            const aol::Vec3<RealType>                    &/*globalCoords*/ ) const {
    aol::Vec3<RealType> approxGrad, differences;

    aol::Matrix33<RealType> dirMat;
    for ( int i = 0; i < 3; ++i ) {
      dirMat.setRow ( i, vertex[i+1] - vertex[0] );
      differences[i] = vertexValue[i+1] - vertexValue[0];
    }

    approxGrad = ( dirMat.inverse() ) * differences;

    approxGrad /= this->_grid.H();

    return ( approxGrad );
  }

};


//! In case of affine FE, interpolation for virtual nodes only uses the cutRelations. Interpolation inside tetrahedra is done via barycentric coordinates.
template< class GridType >
  class AffineFEApproximator : public Approximator<GridType, typename GridType::RealType> {
  typedef typename GridType::RealType RealType;
public:
  AffineFEApproximator ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &data ) : Approximator<GridType, typename GridType::RealType> ( grid, data ) {}

  virtual ~AffineFEApproximator() {}

public:
  virtual void getApproximateVertexFunctionValues ( const tpcfe::CFEElement<RealType> &el, const tpcfe::CFETetra<RealType> &tetra, const std::vector< aol::Vec3<RealType> > &/*vertex*/, aol::Vec<4, RealType> &vertexValue ) const {
    for ( int i = 0; i < 4; i++ ) {
      int gIdx0, gIdx1;

      gIdx0 = el.globalIndex( tetra ( i,0 ) );
      if ( tetra ( i, 1 ) == 11 ) {
        gIdx1 = -1;
      } else {
        gIdx1 = el.globalIndex( tetra ( i,1 ) );
      }

      if ( gIdx1 == -1 ) { // non-interpolated node
        vertexValue[i] = this->_data.get ( gIdx0 );
      } else {           // interpolated node

        // in contrast to CFE, we now interpolate hexahedral grid data according to the cut ratios:
        vertexValue[i] = el.getCutRelation( tetra ( i,1 ),  tetra ( i,0 ) ) * this->_data.get ( gIdx0 ) + el.getCutRelation( tetra ( i,0 ),  tetra ( i,1 ) ) * this->_data.get ( gIdx1 );

      }
    }
  }


  virtual RealType getApproximateFunctionValue ( const std::vector< aol::Vec3<RealType> > &/*vertex*/,
                                                 const aol::Vec<4, RealType> &vertexValue,
                                                 const aol::Vec<4, RealType> &baryCoords,
                                                 const aol::Vec3<RealType> &/*globalCoords*/ ) const {

    return ( baryCoords[0] * vertexValue[0] + baryCoords[1] * vertexValue[1] + baryCoords[2] * vertexValue[2] + baryCoords[3] * vertexValue[3] );
  }

  virtual aol::Vec3<RealType> getApproximateGradientValue ( const std::vector< aol::Vec3<RealType> >     &vertex,
                                                            const aol::Vec<4, RealType>                  &vertexValue,
                                                            const aol::Vec3<RealType>                    &/*globalCoords*/ ) const {
    aol::Vec3<RealType> approxGrad, differences;

    aol::Matrix33<RealType> dirMat;
    for ( int i = 0; i < 3; ++i ) {
      dirMat.setRow ( i, vertex[i+1] - vertex[0] );
      differences[i] = vertexValue[i+1] - vertexValue[0];
    }

    approxGrad = ( dirMat.inverse() ) * differences;

    approxGrad /= this->_grid.H();

    return ( approxGrad );
  }

};


//! In case of Multilinear FE, interpolation happens on each element. This is done via ScalarArray<QC_3D>::interpolate ( global coordinates )
template< class GridType, class RealType >
  class MultilinFEApproximator : public Approximator<GridType, RealType> {

public:
  MultilinFEApproximator ( const GridType &grid, const qc::ScalarArray<RealType, qc::QC_3D> &data ) : Approximator<GridType, RealType> ( grid, data ) {}
  virtual ~MultilinFEApproximator() {}

public:
  virtual void getApproximateVertexFunctionValues ( const tpcfe::CFEElement<RealType> &/*el*/, const tpcfe::CFETetra<RealType> &/*tetra*/, const std::vector< aol::Vec3<RealType> > &vertex, aol::Vec<4, RealType> &vertexValue ) const {
    for ( int i = 0; i < 4; i++ ) {
      vertexValue[i] = this->_data.interpolate ( vertex[i] );
    }
  }


  virtual RealType getApproximateFunctionValue ( const std::vector< aol::Vec3<RealType> > &/*vertex*/,
                                                 const aol::Vec<4, RealType> &/*vertexValue*/,
                                                 const aol::Vec<4, RealType> &/*baryCoords*/,
                                                 const aol::Vec3<RealType> &globalCoords ) const {

    return ( this->_data.interpolate ( globalCoords ) );
  }

  virtual aol::Vec3<RealType> getApproximateGradientValue ( const std::vector< aol::Vec3<RealType> >     &/*vertex*/,
                                                            const aol::Vec<4, RealType>                  &/*vertexValue*/,
                                                            const aol::Vec3<RealType>                    &globalCoords     ) const {
    aol::Vec3<RealType> approxGrad_f;
    aol::Vec3<RealType> approxGrad;

    this->_data.gradient ( globalCoords, approxGrad_f );

    for ( int i = 0; i < 3; ++i ) {
      approxGrad[i] = approxGrad_f[i] / this->_grid.H();
    }

    return ( approxGrad );
  }

};




template<class GridType>
void ComputeApproximationError ( const GridType &grid,
                                 const InterfaceTestFunction< typename GridType::RealType > &testFct,
                                 const Approximator<GridType, typename GridType::RealType> &approximator,
                                 typename GridType::RealType &LinfError, typename GridType::RealType &L2Error, typename GridType::RealType &H1Error, typename GridType::RealType &wH1Error,
                                 typename GridType::RealType &LinfNormFct, typename GridType::RealType &L2NormFct, typename GridType::RealType &H1NormFct, typename GridType::RealType &wH1NormFct,
                                 qc::ScalarArray<typename GridType::RealType, qc::QC_3D> *LinfDifferenceField = NULL,
                                 const int quadOrder = 1, const typename GridType::RealType tetraVolumeThreshold = 0.0 ) {

  // ATTENTION: use global coordinates here

  if ( GridType::CT != tpcfe::CFE_TPOS && GridType::CT != tpcfe::CFE_CDWI_TPOS && GridType::CT != tpcfe::CFE_CDWI_LIEHR) {
    throw aol::Exception ( "Illegal constraint type", __FILE__, __LINE__ );
  }

  typedef typename GridType::RealType RealType;

  RealType
    _Linf_difference    = 0.0,
    _L2_2_difference    = 0.0,
    _L2_2_sum           = 0.0,
    _H1x_2_difference   = 0.0,
    _H1y_2_difference   = 0.0,
    _H1z_2_difference   = 0.0,
    _H1_2_sum           = 0.0,
    _wH1x_2_difference  = 0.0,
    _wH1y_2_difference  = 0.0,
    _wH1z_2_difference  = 0.0,
    _wH1_2_sum          = 0.0,
    _Linf_fct           = 0.0;

  const RealType h = grid.H();

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    el.computeAssembleData ( grid );

    for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {

      const tpcfe::CFETetra<RealType> &tetra = *tit;

      // do not integrate outside of domain
      if( (GridType::CT == tpcfe::CFE_CDWI_TPOS || GridType::CT == tpcfe::CFE_CDWI_LIEHR) && grid.hasDomain() ){
        const int structureNo =  el.cfeType()._structureNo;
        if(el.pureCFEType() == 0){ continue; }
        if(structureNo == tpcfe::MAX_STRUCT_ID && (tetra.getSign() != -1)){ continue; }
      }

      std::vector< aol::Vec3<RealType> >   vertex ( 4 );                        // global coordinates of the vertices
      aol::Vec<4, RealType>                vertexValue;                         // approximate function values at the vertices
      unsigned int numQuadPoints = 0;
      if ( quadOrder == 1 )
        numQuadPoints = 1;
      else if ( quadOrder == 2 )
        numQuadPoints = 8;

      std::vector< aol::Vec<4, RealType> > quadPointBary ( numQuadPoints );     // barycentric coordinates of quadrature points wrt tetra
      std::vector< aol::Vec3<RealType> >   quadPoint ( numQuadPoints );         // global coordinates quadrature points
      std::vector< RealType >              quadWeight ( numQuadPoints );        // quadrature weights
      std::vector< RealType >              quadValue ( numQuadPoints );         // approximate values at quadrature points
      std::vector< aol::Vec3<RealType> >   quadGradValue ( numQuadPoints );     // approximate values of the gradient at quadrature points

      aol::Vec<4, RealType>   exactValueAtVertex;
      std::vector< RealType > exactValueAtQuadPoint ( numQuadPoints ), weightValueAtQuadPoint ( numQuadPoints );

      std::vector< aol::Vec3<RealType> > exactGradValueAtQuadPoint ( numQuadPoints );

      // determine vertices
      {
        for ( int i = 0; i < 4; i++ ) {
          tetra.computeGlobalCoordinate ( vertex[i], el, i );
        }
      }

      // barycentric coordinates of quadPoints
      if ( quadOrder == 1 ) {
        quadPointBary[0][0] = quadPointBary[0][1] = quadPointBary[0][2] = quadPointBary[0][3] = 0.25;
      } else if ( quadOrder == 2 ) {
        quadPointBary[0][0] = 0.5441518;  quadPointBary[0][1] = 0.2939988;  quadPointBary[0][2] = 0.0342028;    quadPointBary[0][3] = 1.0 - quadPointBary[0][0] - quadPointBary[0][1] - quadPointBary[0][2];
        quadPointBary[1][0] = 0.5441518;  quadPointBary[1][1] = 0.0706797;  quadPointBary[1][2] = 0.0813957;    quadPointBary[1][3] = 1.0 - quadPointBary[1][0] - quadPointBary[1][1] - quadPointBary[1][2];
        quadPointBary[2][0] = 0.1225148;  quadPointBary[2][1] = 0.5659332;  quadPointBary[2][2] = 0.0658387;    quadPointBary[2][3] = 1.0 - quadPointBary[2][0] - quadPointBary[2][1] - quadPointBary[2][2];
        quadPointBary[3][0] = 0.1225148;  quadPointBary[3][1] = 0.1360550;  quadPointBary[3][2] = 0.1566826;    quadPointBary[3][3] = 1.0 - quadPointBary[3][0] - quadPointBary[3][1] - quadPointBary[3][2];
        quadPointBary[4][0] = 0.5441518;  quadPointBary[4][1] = 0.2939988;  quadPointBary[4][2] = 0.1276466;    quadPointBary[4][3] = 1.0 - quadPointBary[4][0] - quadPointBary[4][1] - quadPointBary[4][2];
        quadPointBary[5][0] = 0.5441518;  quadPointBary[5][1] = 0.0706797;  quadPointBary[5][2] = 0.3037728;    quadPointBary[5][3] = 1.0 - quadPointBary[5][0] - quadPointBary[5][1] - quadPointBary[5][2];
        quadPointBary[6][0] = 0.1225148;  quadPointBary[6][1] = 0.5659332;  quadPointBary[6][2] = 0.2457133;    quadPointBary[6][3] = 1.0 - quadPointBary[6][0] - quadPointBary[6][1] - quadPointBary[6][2];
        quadPointBary[7][0] = 0.1225148;  quadPointBary[7][1] = 0.1360550;  quadPointBary[7][2] = 0.5847476;    quadPointBary[7][3] = 1.0 - quadPointBary[7][0] - quadPointBary[7][1] - quadPointBary[7][2];
      }

      // global coordinates of quadPoints
      for ( unsigned int q = 0; q < numQuadPoints; ++q ) {
        quadPoint[q] = quadPointBary[q][0] * vertex[0] + quadPointBary[q][1] * vertex[1] + quadPointBary[q][2] * vertex[2] + quadPointBary[q][3] * vertex[3] ;
      }

      // quadrature weights
      if ( quadOrder == 1 ) {
        quadWeight[0] = 1.0;
      } else if ( quadOrder == 2 ) {
        quadWeight[0] = 6 * 0.0091694; // the quadPoints and quadWeights for 2nd order
        quadWeight[1] = 6 * 0.0160270;
        quadWeight[2] = 6 * 0.0211570;
        quadWeight[3] = 6 * 0.0369799;
        quadWeight[4] = 6 * 0.0091694;
        quadWeight[5] = 6 * 0.0160270;
        quadWeight[6] = 6 * 0.0211570;
        quadWeight[7] = 6 * 0.0369799;
      }

      // get vertex values
      approximator.getApproximateVertexFunctionValues ( el, tetra, vertex, vertexValue );

      for ( unsigned int q = 0; q < numQuadPoints; ++q ) {
        quadValue[q]     = approximator.getApproximateFunctionValue ( vertex, vertexValue, quadPointBary[q], quadPoint[q] );
        quadGradValue[q] = approximator.getApproximateGradientValue ( vertex, vertexValue, quadPoint[q] );
      }


      // get exact values
      for ( unsigned int v = 0; v < 4; ++v ) {
        exactValueAtVertex[v] = testFct.valueAt_global ( vertex[v] );
      }

      for ( unsigned int q = 0; q < numQuadPoints; ++q ) {
        exactValueAtQuadPoint[q] = testFct.valueAt_global ( quadPoint[q] );
        exactGradValueAtQuadPoint[q] = testFct.gradient_valueAt_global ( quadPoint[q] );
        weightValueAtQuadPoint[q] = testFct.coefficientValueAtGlobal ( quadPoint[q] );
      }

#ifdef VERBOSE
      for ( int v = 0; v < 4; ++v ) {
        if ( aol::Abs ( vertexValue[v] - exactValueAtVertex[v] ) > ( 1.0e-2 * aol::Abs ( exactValueAtVertex[v] ) ) ) {
          cerr << "FUNCTIONVALUETEST Exact " << aol::mixedFormat ( exactValueAtVertex[v] ) << " approximate " << aol::mixedFormat ( vertexValue[v] )
               << " difference " << aol::mixedFormat ( exactValueAtVertex[v] - vertexValue[v] )
               << " position " << vertex[v] << " volume " << aol::mixedFormat ( tetra.getVolume() ) << " interface " << aol::mixedFormat ( testFct.interfaceValueAt_global ( vertex[v] ) ) << " " << static_cast<int>( tetra.getSign() ) << endl;
        } else {
          if ( tetra.getVolume() < 1.0 / 12.0 )
            cerr << "OKAY volume " << aol::mixedFormat ( tetra.getVolume() ) << " interface " << aol::mixedFormat ( testFct.interfaceValueAt_global ( vertex[v] ) ) << endl;
        }
      }
#endif
#ifdef VERBOSE
      for ( unsigned int q = 0; q < numQuadPoints; ++q ) {
        if ( true ) { // if ( ( exactGradValueAtQuadPoint[q] - quadGradValue[q] ).norm() / ( exactGradValueAtQuadPoint[q] ).norm() > 1e-2 ) { // relative error bigger than 1 %
          cerr << "Exact gradient = " << exactGradValueAtQuadPoint[q]  << ", approximate gradient = " << quadGradValue[q] << " at position " << quadPoint[q]
               << " reldifference " << aol::mixedFormat ( ( exactGradValueAtQuadPoint[q] - quadGradValue[q] ).norm() / exactGradValueAtQuadPoint[q].norm() ) << endl;
        }
      }
      cerr << endl;
#endif

      // update Linf_difference
      for ( int i = 0; i < 4; ++i ) { // vertices
        const RealType difference = aol::Abs ( vertexValue[i] - exactValueAtVertex[i] );
        if ( difference > _Linf_difference ) {
          _Linf_difference = difference;
        }

        if ( aol::Abs ( exactValueAtVertex[i] ) > _Linf_fct ) {
          _Linf_fct = aol::Abs ( exactValueAtVertex[i] );
        }
      }

      for ( unsigned int q = 0; q < numQuadPoints; ++q ) { // quadrature points
        const RealType difference = aol::Abs ( quadValue[q] - exactValueAtQuadPoint[q] );
        if ( difference > _Linf_difference ) {
          _Linf_difference = difference;
        }
      }

      if ( LinfDifferenceField ) {

        for ( int i = 0; i < 4; ++i ) { // vertices
          qc::CoordType vPos ( static_cast<short>( vertex[i][0] ), static_cast<short>( vertex[i][1] ), static_cast<short>( vertex[i][2] ) );
          const RealType difference = aol::Abs ( vertexValue[i] - exactValueAtVertex[i] );

#ifdef VERBOSE
          if ( difference > 0.1 )
            cerr << "FUNCTIONVALUETEST Exact " << aol::mixedFormat ( exactValueAtVertex[i] ) << " approximate " << aol::mixedFormat ( vertexValue[i] )
                 << " difference " << aol::mixedFormat ( exactValueAtVertex[i] - vertexValue[i] )
                 << " position " << vertex[i] << " volume " << aol::mixedFormat ( tetra.getVolume() ) << " interface " << aol::mixedFormat ( testFct.interfaceValueAt_global ( vertex[i] ) ) << " " << static_cast<int>( tetra.getSign() ) << endl;
#endif

          if ( difference > LinfDifferenceField->get ( vPos ) )
            LinfDifferenceField->set ( vPos, difference );
        }

        for ( unsigned int q = 0; q < numQuadPoints; ++q ) { // quadrature points
          qc::CoordType qPos ( static_cast<short>( quadPoint[q][0] ), static_cast<short>( quadPoint[q][1] ), static_cast<short>( quadPoint[q][2] ) );
          const RealType difference = aol::Abs ( quadValue[q] - exactValueAtQuadPoint[q] );

#ifdef VERBOSE
          if ( difference > 0.1 )
            cerr << "FUNCTIONVALUETEST Exact " << aol::mixedFormat ( exactValueAtQuadPoint[q] ) << " approximate " << aol::mixedFormat ( quadValue[q] )
                 << " difference " << aol::mixedFormat ( exactValueAtQuadPoint[q] - quadValue[q] )
                 << " position " << quadPoint[q] << " volume " << aol::mixedFormat ( tetra.getVolume() ) << " interface " << aol::mixedFormat ( testFct.interfaceValueAt_global ( quadPoint[q] ) ) << " " << static_cast<int>( tetra.getSign() ) << endl;
#endif

          if ( difference > LinfDifferenceField->get ( qPos ) )
            LinfDifferenceField->set ( qPos, difference );
        }

      }

      const RealType tetraVol = tetra.getVolume(), tvh3 = tetraVol * aol::Cub ( h );
      if ( tetraVol > tetraVolumeThreshold ) {

        for ( unsigned int q = 0; q < numQuadPoints; ++q ) { // quadrature points
          // sum up contributions to L2 difference
          _L2_2_difference += quadWeight[q] * tvh3 * aol::Sqr ( quadValue[q] - exactValueAtQuadPoint[q] );
          _L2_2_sum        += quadWeight[q] * tvh3 * aol::Sqr ( exactValueAtQuadPoint[q] );

          // sum up contributions to H1 difference
          const RealType
            h1xContrib   = quadWeight[q] * tvh3 * aol::Sqr ( quadGradValue[q][0] - exactGradValueAtQuadPoint[q][0] ),
            h1yContrib   = quadWeight[q] * tvh3 * aol::Sqr ( quadGradValue[q][1] - exactGradValueAtQuadPoint[q][1] ),
            h1zContrib   = quadWeight[q] * tvh3 * aol::Sqr ( quadGradValue[q][2] - exactGradValueAtQuadPoint[q][2] ),
            h1SumContrib = quadWeight[q] * tvh3 * ( aol::Sqr ( exactValueAtQuadPoint[q] ) + aol::Sqr ( quadGradValue[q][0] ) + aol::Sqr ( quadGradValue[q][1] ) + aol::Sqr ( quadGradValue[q][2] ) ),
            weight = weightValueAtQuadPoint[q];

          _H1x_2_difference  += h1xContrib;
          _H1y_2_difference  += h1yContrib;
          _H1z_2_difference  += h1zContrib;
          _H1_2_sum          += h1SumContrib;

          _wH1x_2_difference  += aol::Sqr ( weight ) * h1xContrib;
          _wH1y_2_difference  += aol::Sqr ( weight ) * h1yContrib;
          _wH1z_2_difference  += aol::Sqr ( weight ) * h1zContrib;
          _wH1_2_sum          += aol::Sqr ( weight ) * h1SumContrib;

          if ( aol::Abs ( exactValueAtQuadPoint[q] ) > _Linf_fct ) {
            _Linf_fct = aol::Abs ( exactValueAtQuadPoint[q] );
          }

        }
      }

    }
  }

  // "return" values
  LinfError   = _Linf_difference;
  L2Error     = sqrt ( _L2_2_difference );
  H1Error     = sqrt ( _L2_2_difference  + _H1x_2_difference + _H1y_2_difference + _H1z_2_difference ) ;
  wH1Error    = sqrt ( _L2_2_difference  + _wH1x_2_difference + _wH1y_2_difference + _wH1z_2_difference ) ;

  LinfNormFct = _Linf_fct;
  L2NormFct   = sqrt ( _L2_2_sum );
  H1NormFct   = sqrt ( _H1_2_sum );
  wH1NormFct  = sqrt ( _wH1_2_sum );

}

}

#endif
