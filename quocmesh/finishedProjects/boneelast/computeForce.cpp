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

/* this program loads a geometry and deformations and computes the force on a planar surface
*/

#include <parameterParser.h>
#include <preconditioner.h>
#include <solver.h>

#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEUtils.h>

//#define VERBOSE
#include "computeForce.h"

typedef double RealType;

// TODO: think about what to do concerning ConfiguratorType
template< class GridType >
void computeForceViaOp ( GridType &grid,
                         qc::MultiArray< RealType, 3> &deformations,
                         const unsigned int surfaceDirection, const int surfacePosition,
                         const RealType aleph, const RealType lambda, const RealType mu,
                         ostream &out ) {

  const RealType tolerance = 1.0e-6;
  { RealType bla = tolerance; bla *= aleph; } // to prevent unused warning ...

  typedef tpcfe::CFEConfigurator < GridType, tpcfe::CFEHybridMatrix<GridType> > ConfiguratorType;
  tpcfe::CFEElastOp< ConfiguratorType > elastop ( grid, lambda, mu );
  tpcfe::CFEMassOp< ConfiguratorType > massOp ( grid, aol::ASSEMBLED );
  massOp.assembleMatrix();
  // role of aleph??

  qc::MultiArray< RealType, 3 > tmp ( deformations, aol::STRUCT_COPY ), forces ( deformations, aol::STRUCT_COPY );
  elastop.apply ( deformations, tmp );

  aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( massOp.getMatrixRef() );
  aol::PCGInverse< aol::Vector<RealType > > massOpInv( massOp, prec, 1.0e-16, 1000 );
  massOpInv.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
  //   tpcfe::CFEMultigrid< tpcfe::CFEMassOp< ConfiguratorType > > massOpInv ( grid, massOp.getMatrixRef(), 2, 3, 3, aol::GAUSS_SEIDEL_FORWARD, 1.0, mg::MULTIGRID_V_CYCLE, 1.0e-16, 1000 );
  aol::DiagonalBlockOp< RealType > blockMassOpInv ( massOpInv );

  blockMassOpInv.apply ( tmp, forces );

  qc::ScalarArray<RealType, qc::QC_3D> forces_save ( forces[0], aol::STRUCT_COPY );
  for ( int i = 0; i < forces_save.getNumX(); ++i )
    for ( int j = 0; j < forces_save.getNumX(); ++j )
      for ( int k = 0; k < forces_save.getNumX(); ++k )
        forces_save.set ( i, j, k,  ( forces.get( aol::Vec3<short>(i, j, k) ).norm() == 0.0 ? 0.0 : 1.0 ) );
  forces_save.saveRaw ( "out/forces_struct.raw", qc::PGM_UNSIGNED_CHAR_BINARY, 0.0, 1.0 );

  forces.saveRaw ( "out/forces%d.raw", qc::PGM_UNSIGNED_CHAR_BINARY, 0, forces.getMaxAbsValue() );

  const int beginPos = ( surfacePosition == -1 ? 0 : surfacePosition );
  const int endPos =   ( surfacePosition == -1 ? grid.getSize()[surfaceDirection] : surfacePosition+1 );

  std::map< int, aol::Vec3<RealType> > vecForcesStorage; // ( endPos ); // no big waste; want to start indexing at beginPos.

#if 1

  for ( qc::RectangularIterator<qc::QC_3D> rit ( grid ); rit.notAtEnd(); ++rit ) {
    if ( grid.getDomainRef().get ( *rit ) < 0 ) {
      // we are only interested in regular nodes here, the virtual nodes are treated below.
      vecForcesStorage[ (*rit)[ surfaceDirection ] ] += aol::Sqr( grid.H() ) * forces.get( *rit );
    }
  }

  // \todo use the virtual triangles in the cutting / intersection plane
  for ( typename GridType::VMapTypeConstIt it = grid.virtualNodeMapRef().begin(); it != grid.virtualNodeMapRef().end(); ++it ) {
    //    cerr << it->second->coord << " " << it->second->extrapolate ( deformations[0] ) << " " << it->second->extrapolate ( deformations[1] ) << " " << it->second->extrapolate ( deformations[2] ) << endl;
  }

#else

  for ( int surfacePos = beginPos; surfacePos < endPos; ++surfacePos ) {
    tpcfe::CFEElementIterator<RealType> elit = grid.elementIterator();
    for ( elit.start(); !(elit.atEnd()); ++elit ) { // loop over all elements = cubes
      const short type = elit->cfeType(); // type characterizes how the cube is interfaced by the interface, i. e. signature of eight nodes (inside/outside)
      if ( type != 0 && type != 255 )
        elit->computeCutRelations ( grid ); // this fails for types 0 and 255 (non-interfaced cubes), but there are no cut relations anyway

      tpcfe::CFETetrahedronIterator<RealType> tit;
      for ( tit.start ( type, -1 ); !(tit.atEnd()); ++tit ) { // loop over all negative (i. e. inner) tetrahedra

        const tpcfe::CFETetra<RealType> *tetra = *tit;
        aol::Vec3<RealType> VertexCoords[4]; // these are in image coordinates / quoc units
        int boundary_type = 0, boundary_points = 0;

        for ( int i = 0; i < 4; i++ ) { // loop over vertices of tetra to find out where they lie in relation to the plane
          tetra->computeGlobalCoordinate ( VertexCoords[i], ( *elit ), i );
          if( aol::Abs ( VertexCoords[i][surfaceDirection] - surfacePos ) > tolerance )
            boundary_type += ( 1 << i );
          else
            ++boundary_points;
        }

        if ( boundary_points == 3 ) { // these have one in-plane face

          aol::Vec3<int> mi; // mapping of indices such that vertices 0, 1, 2 lie in the plane and 3 does not
          int nonbnd_vertex;

          switch( boundary_type ) {
          case 1:
            mi.set( 1, 2, 3 );
            nonbnd_vertex = 0;
            break;
          case 2:
            mi.set( 0, 2, 3 );
            nonbnd_vertex = 1;
            break;
          case 4:
            mi.set( 0, 1, 3 );
            nonbnd_vertex = 2;
            break;
          case 8:
            mi.set( 0, 1, 2 );
            nonbnd_vertex = 3;
              break;
          default:
            mi.set( -1,-1,-1);
            nonbnd_vertex = -1;
          }

          aol::Vec3<RealType> force_at_vertex[4];

          for ( int i = 0; i < 4; ++i ) { // i = vertex of tetra

            const int
              lIdx0 = ( *tetra ) ( i, 0 ),
              lIdx1 = ( *tetra ) ( i, 1 ),
              gIdx0 = elit->_globalIndex[lIdx0];

            if ( lIdx1 == 11 ) { // at regular nodes, we have data
              force_at_vertex[i][0] = forces[0][gIdx0];
              force_at_vertex[i][1] = forces[1][gIdx0];
              force_at_vertex[i][2] = forces[2][gIdx0];
            } else {             // at virtual nodes, we need to interpolate
              const int gIdx1 = elit->_globalIndex[lIdx1];
              const tpcfe::CFEVirtualNode<RealType, tpcfe::CFE_CD> &vn =  grid.getVirtualNodeRef ( gIdx0, gIdx1 );
              force_at_vertex[i][0] = vn.extrapolate ( forces[0] );
              force_at_vertex[i][1] = vn.extrapolate ( forces[1] );
              force_at_vertex[i][2] = vn.extrapolate ( forces[2] );
            }
          }

          aol::Vec3<RealType> vec_force = ( 1.0 / 3.0 ) * ( force_at_vertex[0] + force_at_vertex[1] + force_at_vertex[2] );

          aol::Vec3<RealType> vertex[3];
          for ( int i = 0; i < 3; ++i ) {
            vertex[i] = aleph * grid.H() * VertexCoords[ mi[i] ]; // usage of aleph correct??
          }

          // here, we compute the area of the in-plane triangle:
          aol::Vec3<RealType> crossProduct;
          crossProduct.makeCrossProduct( vertex[2] - vertex[0], vertex[2] - vertex[1] );

          const RealType area_contribution = 0.5 * crossProduct.norm();
          // 1/2 for area of triangle instead of parallelogram;

          // compute and add the contributions to total force
          aol::Vec3<RealType> vec_force_contrib = area_contribution * vec_force;

          vecForcesStorage[ surfacePos ] += vec_force_contrib;
        }
      }
    }
  }
#endif

  for ( int surfacePos = beginPos; surfacePos < endPos; ++surfacePos ) {
    aol::Vec3<RealType> total_force;
    out << "at surfacePos = " << surfacePos;
    out << ": total force = " << vecForcesStorage[surfacePos] << endl;
  }

}


int main ( int argc, char** argv ) {
  try {

    aol::ParameterParser params ( ( argc == 2 ? argv[1] : "computeForce.par" ) ); // read parameters from file

    //    params.dump();

    const int depth = params.getInt ( "depth" );
    const RealType aleph = static_cast<RealType>( params.getDouble("aleph") );    // unit is m / 1 by def

    const RealType E = static_cast<RealType>( params.getDouble("Emodulus") );     // unit is N m^{-2} by def
    const RealType nu = static_cast<RealType>( params.getDouble("nu") );          // unit is 1 by def


    cerr << "Creating domain ... ";

    tpcfe::CFEGrid< RealType, tpcfe::CFE_CD, tpcfe::IsotropicElasticityCoefficient<RealType> > grid ( depth );
    qc::ScalarArray<RealType, qc::QC_3D> levelset ( grid );

    char levelset_filename[1024];
    params.getString ( "geometry_filename", levelset_filename );
    cerr << levelset_filename << endl;
    levelset.load ( levelset_filename );

    grid.setDomainFrom ( levelset );
    grid.detectAndInitVirtualNodes();

    qc::BitArray<qc::QC_3D> dummy_dirichlet ( qc::GridSize<qc::QC_3D>::createFrom ( grid ) );
    dummy_dirichlet.setAll ( false );
    grid.setDirichletMask ( dummy_dirichlet );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    cerr << " done." << endl;

    cerr << "Total inner volume in SI units is " << grid.getTotalInnerVolume() * aol::Cub( aleph ) << endl;

    cerr << "Loading deformations ";

    qc::MultiArray< RealType, 3> deformation ( grid );                           // unit is m -- this means correct boundary conditions are necessary in the elasticity computation
    char deformations_filenamemask[1024];
    params.getString ( "deformations_fnmask", deformations_filenamemask );
    cerr << deformations_filenamemask << endl;
    deformation.load ( deformations_filenamemask ); // these are assumed to contain deformation in SI units (i. e. meters)
    cerr << "x min = " << deformation[0].getMinValue() << ", x max = " << deformation[0].getMaxValue() << endl
         << "y min = " << deformation[1].getMinValue() << ", y max = " << deformation[1].getMaxValue() << endl
         << "z min = " << deformation[2].getMinValue() << ", z max = " << deformation[2].getMaxValue() << endl << endl;

    cerr << " done." << endl;

#if 0
    // necessary if we deal with outputs from the solutionStepSaver
    cerr << "Loading and adding boundary conditions";
    qc::MultiArray< RealType, 3 > BCs ( grid );
    char dirichletBC_fnmask[1024];
    params.getString ( "dirichletBCs_fnmask", dirichletBC_fnmask );
    BCs.load(dirichletBC_fnmask);
    deformation += BCs;
    cerr << " done. " << endl;
#endif

    const unsigned int surface_direction = params.getInt( "surface_direction" );   // unsigned direction of the surface normal

    const RealType tolerance = 1e-8; // for in-plane points, not in SI units

    aol::Vec3<RealType> surface_outer_normal;
    surface_outer_normal [ surface_direction ] = 1.0;

    const unsigned int surface_position  = params.getInt( "surface_position" );    // index of the layer, no unit
    unsigned int surface_pos = surface_position;

    aol::Vector<double> fst_dummy;

    tpcfe::IsotropicElasticityCoefficient<RealType> ENuMinus ( E, nu ), ENuPlus ( aol::NumberTrait<RealType>::NaN, aol::NumberTrait<RealType>::NaN );
    qc::AArray< tpcfe::IsotropicElasticityCoefficient<RealType>, qc::QC_3D > isoCoeff ( grid );
    tpcfe::setCoeffForLevelset ( isoCoeff, levelset, ENuMinus, ENuPlus );

    tpcfe::computeForce ( grid, deformation, surface_direction, surface_outer_normal, tolerance, surface_pos, aleph, isoCoeff, cerr, fst_dummy );

  } catch ( aol::Exception e ) {
    e.dump();
    return ( EXIT_FAILURE );
  }

  return ( EXIT_SUCCESS );
}
