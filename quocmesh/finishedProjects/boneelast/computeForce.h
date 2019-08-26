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

#ifndef __COMPUTEFORCE_H
#define __COMPUTEFORCE_H

namespace tpcfe {

template< class GridType >
void computeForce ( const GridType &grid,
                    const qc::MultiArray< typename GridType::RealType, 3> &deformations,
                    const unsigned int surfaceDirection, const aol::Vec3< typename GridType::RealType > &surfaceOuterNormal, const typename GridType::RealType tolerance, const int surfacePosition,
                    const typename GridType::RealType aleph, const qc::AArray< typename GridType::NodalCoeffType, qc::QC_3D > &nodalCoeff,
                    ostream &out, aol::Vector< typename GridType::RealType > &forces_storage ) {

  // outer normal is not redundant because of direction

  typedef typename GridType::RealType RealType;

  const int
    beginPos = ( surfacePosition == -1 ? 0 : surfacePosition ),
    endPos =   ( surfacePosition == -1 ? grid.getSize()[surfaceDirection] : surfacePosition + 1 );

  forces_storage.resize ( endPos - beginPos );

  for ( int surfacePos = beginPos; surfacePos < endPos; ++surfacePos ) {

    //      out << "Position: " << surfacePos << " ";

    bool triang_parallel_normal = false, triang_antipara_normal = false;
    RealType  total_area = 0.0; // this sums up the area, unit is m^{2}

    aol::Vec3<RealType> total_force; // set zero automatically, this sums up the force contributions, unit is N

    const qc::GridSize<qc::QC_3D> gridSize ( grid );
    for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {

      if ( aol::Abs ( (*it)[surfaceDirection] - surfacePos ) < 2 ) {

        tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

        el.computeAssembleData ( grid );

        tpcfe::CFEWeightProvider< RealType, typename GridType::NodalCoeffType > weightProvider ( nodalCoeff, el );

        signed char which = -2;
        if ( GridType::CT == tpcfe::CFE_CD )
          which = -1;
        else if ( GridType::CT == tpcfe::CFE_TPOSELAST )
          which = 0;
        else
          throw aol::Exception ( "tpcfe::computeForce: illegal constraint type", __FILE__, __LINE__ );

        for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, which ); tit.notAtEnd(); ++tit ) { // loop over all negative (i. e. inner) tetrahedra

          const tpcfe::CFETetra<RealType>& tetra = *tit;
          aol::Vec3<RealType> VertexCoords[4]; // these are in image coordinates / quoc units
          int boundary_type = 0, boundary_points = 0;

          for ( int i = 0; i < 4; i++ ) { // loop over vertices of tetra to find out where they lie in relation to the plane
            tetra.computeGlobalCoordinate ( VertexCoords[i], el, i );
            if ( aol::Abs ( VertexCoords[i][surfaceDirection] - surfacePos ) > tolerance )
              boundary_type += ( 1 << i );
            else
              ++boundary_points;
          }

          if ( boundary_points == 3 ) { // these have one in-plane face

            aol::Vec3<int> mi; // mapping of indices such that vertices 0, 1, 2 lie in the plane and 3 does not
            int nonbnd_vertex;

            switch ( boundary_type ) {
            case 1:
              mi.set ( 1, 2, 3 );
              nonbnd_vertex = 0;
              break;
            case 2:
              mi.set ( 0, 2, 3 );
              nonbnd_vertex = 1;
              break;
            case 4:
              mi.set ( 0, 1, 3 );
              nonbnd_vertex = 2;
              break;
            case 8:
              mi.set ( 0, 1, 2 );
              nonbnd_vertex = 3;
              break;
            default:
              mi.set ( -1, -1, -1 );
              nonbnd_vertex = -1;
            }

#ifdef VERBOSE
            cerr << mi[0] << " " << mi[1] << " " << mi[2] << " " << nonbnd_vertex << " " << endl;
            cerr << "Vertex coords: " << VertexCoords[mi[0]] << "     " << VertexCoords[mi[1]] << "     " << VertexCoords[mi[2]] << "     " << VertexCoords[nonbnd_vertex] << endl;
#endif

            aol::Vec3<RealType> directions[3]; // directed edges of the tetra with adjacent out-of-plane vertex, unit is m
            for ( int i = 0; i < 3; ++i ) {
              directions[i] = aleph * grid.H() * ( VertexCoords[ nonbnd_vertex ] - VertexCoords[ mi[i] ] );  // conversion of quoc image coordinates to world coordinates in SI units
            }

#ifdef VERBOSE
            cerr << "Directions: " << directions[0] << "         " << directions[1] << "         " << directions[2] << "         " << endl;
#endif

            triang_parallel_normal |= ( ( directions[0] * surfaceOuterNormal ) < 0 ); // detected tetra in opposite normal direction so far?
            triang_antipara_normal |= ( ( directions[0] * surfaceOuterNormal ) > 0 ); // detected tetea in normal direction so far?

            aol::Matrix33<RealType> directions_matrix;
            for ( int i = 0; i < 3; ++i )
              directions_matrix.setRow ( i, directions[i] );                    // unit is m

#ifdef VERBOSE
            cerr << "Directions matrix: " << endl << directions_matrix << endl;
#endif

            aol::Matrix33<RealType> TetraMatrix = directions_matrix.inverse();  // unit is m^{-1}
            // for computing the gradient from the difference quotients in direction of the edges above -- this is exact for tetra-wise linear functions

#ifdef VERBOSE
            cerr << "Tetra matrix: " << endl << TetraMatrix << endl;
#endif

            aol::Vec3<RealType> deformation_at_vertex[4];                       // unit is m

            for ( int i = 0; i < 4; ++i ) { // i = vertex of tetra

              const int
                lIdx0 = tetra ( i, 0 ),
                lIdx1 = tetra ( i, 1 ),
                gIdx0 = el.globalIndex(lIdx0);

              if ( lIdx1 == 11 ) { // at regular nodes, we have data
                deformation_at_vertex[i][0] = deformations[0][gIdx0];            // unit of deformation is m
                deformation_at_vertex[i][1] = deformations[1][gIdx0];
                deformation_at_vertex[i][2] = deformations[2][gIdx0];
              } else {             // at virtual nodes, we need to interpolate
                const int gIdx1 = el.globalIndex(lIdx1);
                const typename GridType::VNType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
                vn.extrapolate ( deformations, deformation_at_vertex[i] ); // extrapolation preserves unit
              }

            }

            for ( int i = 0; i < 4; ++i ) {
#ifdef VERBOSE
              cerr << "Deformations at vertex " << i << ": " << deformation_at_vertex[i][0] << "  " << deformation_at_vertex[i][1] << "  " << deformation_at_vertex[i][2] << endl;
#endif
#ifdef GNUPLOT_OUTPUT
              cerr << VertexCoords[i][0] << " " << VertexCoords[i][1] << " " << VertexCoords[i][2] << " "
                   << deformation_at_vertex[i][0] << "  " << deformation_at_vertex[i][1] << "  " << deformation_at_vertex[i][2] << "   plot_deform" << endl;
#endif
            }

            aol::Vec3<RealType>
              difference_deformations[3],                                       // unit is m
              gradient_deformations[3];                                         // unit is 1

            for ( int i = 0; i < 3; ++i ) {
              for ( int j = 0; j < 3; ++j ) { // difference between the deformations at out-of-plane vertex and other vertices

                difference_deformations[j][i] = deformation_at_vertex[ nonbnd_vertex ][j] - deformation_at_vertex[ mi[i] ][j];       // unit is m
                // first index in difference = spatial direction of deformation, this is second index in deformation_at_vertex

              }
            }

#ifdef VERBOSE
            cerr << "Differences: " << difference_deformations[0] << "                " << difference_deformations[1] << "                " << difference_deformations[2] << endl;
#endif

            for ( int i = 0; i < 3; ++i ) { // here, we compute the gradient
              gradient_deformations[i] = TetraMatrix * difference_deformations[i];   // unit is m m^{-1} = 1
            }

#ifdef VERBOSE
            cerr << "Gradients: " << gradient_deformations[0] << "                " << gradient_deformations[1] << "                " << gradient_deformations[2] << endl;
#endif
#ifdef GNUPLOT_OUTPUT
            cerr << 0.25 * ( VertexCoords[0][0] + VertexCoords[1][0] + VertexCoords[2][0] + VertexCoords[3][0] ) << " "
                 << 0.25 * ( VertexCoords[0][1] + VertexCoords[1][1] + VertexCoords[2][1] + VertexCoords[3][1] ) << " "
                 << 0.25 * ( VertexCoords[0][2] + VertexCoords[1][2] + VertexCoords[2][2] + VertexCoords[3][2] ) << " "
                 << gradient_deformations[0][0] << " " << gradient_deformations[0][1] << " " << gradient_deformations[0][2] << " "
                 << gradient_deformations[1][0] << " " << gradient_deformations[1][1] << " " << gradient_deformations[1][2] << " "
                 << gradient_deformations[2][0] << " " << gradient_deformations[2][1] << " " << gradient_deformations[2][2] << "   plot_gradient " << endl;
#endif
            aol::Matrix33<RealType> sigma; // this will become the sigma tensor, unit is N m^{-2}

            RealType CTensor[3][3][3][3];
            weightProvider.meanWeight ( tetra.getSign() ).getAnisotropicTensor ( CTensor );

            for ( short i = 0; i < 3; ++i ) {
              for ( short j = 0; j < 3; ++j ) {
                for ( short k = 0; k < 3; ++k ) {
                  for ( short l = 0; l < 3; ++l ) {
                    sigma.add ( i, j, CTensor[i][j][k][l] * ( gradient_deformations[k][l] + gradient_deformations[l][k] ) / 2. );
                  }
                }
              }
            }

            aol::Vec3<RealType> sigma_times_normal = sigma * surfaceOuterNormal;   // unit is N m^{-2} * 1 = N m^{-2}

            aol::Vec3<RealType> vertex[3]; // in-plane vertex coordinates, unit is m
            for ( int i = 0; i < 3; ++i ) {
              vertex[i] = aleph * grid.H() * VertexCoords[ mi[i] ];  // conversion to SI units and global coordinates
            }

            // here, we compute the area of the in-plane triangle:
            aol::Vec3<RealType> crossProduct;                                   // unit is m^{2}
            crossProduct.makeCrossProduct ( vertex[2] - vertex[0], vertex[2] - vertex[1] );

            const RealType area_contribution = 0.5 * crossProduct.norm();       // unit is m^{2}
            // 1/2 for area of triangle instead of parallelogram;

#ifdef VERBOSE
            cerr << "Area: " << area_contribution << endl;
#endif

            total_area += area_contribution;

            // compute and add the contributions to total force
            aol::Vec3<RealType> vec_force_contrib = area_contribution * sigma_times_normal;  // unit is m^{2} * N m^{-2} = N
            //            cerr << vec_force_contrib << endl;
            total_force += vec_force_contrib;                                   // unit is N

#ifdef VERBOSE
            cerr << "Force: " << vec_force_contrib << endl;
            cerr << endl;
#endif

          }
        }
      }
    }

    // if we have a cutting plane for which we sum up upper and lower tetrahedra, we consider two times the area. So we need to divide by 2 ...
    if ( triang_parallel_normal && triang_antipara_normal ) {
      total_area  /= 2.0;
      total_force /= 2.0;
    }

    out << "area = " << total_area << " m^2." << endl;
    out << "force (vector) = " << total_force << "   ";
    out << "force norm     = " << total_force.norm() << " N." << endl;

    forces_storage[ surfacePos - beginPos ] = total_force.norm();

  }

  cerr << "Average force = " << aol::longScientificFormat ( forces_storage.getMeanValue() )
       << ", standard deviation = " << aol::longScientificFormat ( forces_storage.getStdDev() )
       << ", relative std dev = " << aol::longScientificFormat ( forces_storage.getStdDev() / forces_storage.getMeanValue() ) << endl;

}

template< class GridType >
void getAverageSigmasEpsilons ( GridType &grid,
                                qc::MultiArray< typename GridType::RealType, 3> &deformations,
                                const unsigned int surfaceDirection, const aol::Vec3<typename GridType::RealType> &surfaceOuterNormal, const typename GridType::RealType tolerance,
                                const typename GridType::RealType aleph, const typename GridType::RealType lambda, const typename GridType::RealType mu,
                                std::vector< aol::Matrix33<typename GridType::RealType> > &sigmas, std::vector< aol::Matrix33<typename GridType::RealType> > &epsilons ) {
  typedef typename GridType::RealType RealType;

  sigmas.clear();   sigmas.resize ( grid.getSize()[surfaceDirection] );
  epsilons.clear(); epsilons.resize ( grid.getSize()[surfaceDirection] );

  for ( int surfacePos = 0; surfacePos < grid.getSize()[surfaceDirection] ; ++surfacePos ) {

    bool triang_parallel_normal = false, triang_antipara_normal = false;

    const qc::GridSize<qc::QC_3D> gridSize ( grid );
    for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
      tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

      const short type = el.cfeType(); // type characterizes how the cube is interfaced by the interface, i. e. signature of eight nodes (inside/outside)
      if ( type != 0 && type != 255 )
        el.computeCutRelations ( grid ); // this fails for types 0 and 255 (non-interfaced cubes), but there are no cut relations anyway

      for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); ! ( tit.atEnd() ); ++tit ) { // loop over all negative (i. e. inner) tetrahedra

        const tpcfe::CFETetra<RealType> *tetra = *tit;
        aol::Vec3<RealType> VertexCoords[4]; // these are in image coordinates / quoc units
        int boundary_type = 0, boundary_points = 0;

        for ( int i = 0; i < 4; i++ ) { // loop over vertices of tetra to find out where they lie in relation to the plane
          tetra->computeGlobalCoordinate ( VertexCoords[i], el, i );
          if ( aol::Abs ( VertexCoords[i][surfaceDirection] - surfacePos ) > tolerance )
            boundary_type += ( 1 << i );
          else
            ++boundary_points;
        }

        if ( boundary_points == 3 ) { // these have one in-plane face

          aol::Vec3<int> mi; // mapping of indices such that vertices 0, 1, 2 lie in the plane and 3 does not
          int nonbnd_vertex;

          switch ( boundary_type ) {
          case 1:
            mi.set ( 1, 2, 3 );
            nonbnd_vertex = 0;
            break;
          case 2:
            mi.set ( 0, 2, 3 );
            nonbnd_vertex = 1;
            break;
          case 4:
            mi.set ( 0, 1, 3 );
            nonbnd_vertex = 2;
            break;
          case 8:
            mi.set ( 0, 1, 2 );
            nonbnd_vertex = 3;
            break;
          default:
            mi.set ( -1, -1, -1 );
            nonbnd_vertex = -1;
          }

          aol::Vec3<RealType> directions[3]; // directed edges of the tetra with adjacent out-of-plane vertex, unit is m
          for ( int i = 0; i < 3; ++i ) {
            directions[i] = aleph * grid.H() * ( VertexCoords[ nonbnd_vertex ] - VertexCoords[ mi[i] ] );  // conversion of quoc image coordinates to world coordinates in SI units
          }

          triang_parallel_normal |= ( ( directions[0] * surfaceOuterNormal ) < 0 ); // detected tetra in opposite normal direction so far?
          triang_antipara_normal |= ( ( directions[0] * surfaceOuterNormal ) > 0 ); // detected tetea in normal direction so far?

          aol::Matrix33<RealType> directions_matrix;
          for ( int i = 0; i < 3; ++i )
            directions_matrix.setRow ( i, directions[i] );                    // unit is m

          aol::Matrix33<RealType> TetraMatrix = directions_matrix.inverse();  // unit is m^{-1}
          // for computing the gradient from the difference quotients in direction of the edges above -- this is exact for tetra-wise linear functions

          aol::Vec3<RealType> deformation_at_vertex[4];                       // unit is m

          for ( int i = 0; i < 4; ++i ) { // i = vertex of tetra

            const int
              lIdx0 = ( *tetra ) ( i, 0 ),
              lIdx1 = ( *tetra ) ( i, 1 ),
              gIdx0 = el.globalIndex(lIdx0);

            if ( lIdx1 == 11 ) { // at regular nodes, we have data
              deformation_at_vertex[i][0] = deformations[0][gIdx0];            // unit of deformation is m
              deformation_at_vertex[i][1] = deformations[1][gIdx0];
              deformation_at_vertex[i][2] = deformations[2][gIdx0];
            } else {             // at virtual nodes, we need to interpolate
              const int gIdx1 = el.globalIndex(lIdx1);
              const tpcfe::CFEVirtualNode<RealType, tpcfe::CFE_CD, RealType> &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
              deformation_at_vertex[i][0] = vn.extrapolate ( deformations[0] );  // extrapolation preserves unit
              deformation_at_vertex[i][1] = vn.extrapolate ( deformations[1] );
              deformation_at_vertex[i][2] = vn.extrapolate ( deformations[2] );
            }

          }

          aol::Vec3<RealType>
            difference_deformations[3],                                       // unit is m
            gradient_deformations[3];                                         // unit is 1

          for ( int i = 0; i < 3; ++i ) {
            for ( int j = 0; j < 3; ++j ) { // difference between the deformations at out-of-plane vertex and other vertices

              difference_deformations[j][i] = deformation_at_vertex[ nonbnd_vertex ][j] - deformation_at_vertex[ mi[i] ][j];       // unit is m
              // first index in difference = spatial direction of deformation, this is second index in deformation_at_vertex

            }
          }

          for ( int i = 0; i < 3; ++i ) { // here, we compute the gradient
            gradient_deformations[i] = TetraMatrix * difference_deformations[i];   // unit is m m^{-1} = 1
          }

          aol::Matrix33<RealType> epsilon;
          for ( int i = 0; i < 3; ++i ) {
            for ( int j = 0; j < 3; ++j ) {
              epsilon.set ( i, j,  0.5 * ( gradient_deformations[i][j] + gradient_deformations[j][i] )  );
            }
          }
          aol::Matrix33<RealType> sigma; // this will become the sigma tensor, unit is N m^{-2}

          for ( int i = 0; i < 3; ++i ) {
            for ( int j = 0; j < 3; ++j ) {
              sigma.set ( i, j, mu * ( gradient_deformations[i][j] + gradient_deformations[j][i] ) );
            }
            sigma.add ( i, i,  lambda * ( gradient_deformations[0][0] + gradient_deformations[1][1] + gradient_deformations[2][2] ) );
          }

          aol::Vec3<RealType> vertex[3]; // in-plane vertex coordinates, unit is m
          for ( int i = 0; i < 3; ++i ) {
            vertex[i] = aleph * grid.H() * VertexCoords[ mi[i] ];  // conversion to SI units and global coordinates
          }

          // here, we compute the area of the in-plane triangle:
          aol::Vec3<RealType> crossProduct;                                   // unit is m^{2}
          crossProduct.makeCrossProduct ( vertex[2] - vertex[0], vertex[2] - vertex[1] );

          const RealType area_contribution = 0.5 * crossProduct.norm();       // unit is m^{2}
          // 1/2 for area of triangle instead of parallelogram;

          epsilon *= area_contribution;
          sigma *= area_contribution;
          epsilons[ surfacePos ] += epsilon;
          sigmas[ surfacePos ] += sigma;

        }
      }
    }

    // if we have a cutting plane for which we sum up upper and lower tetrahedra, we consider two times the area. So we need to divide by 2 ...
    if ( triang_parallel_normal && triang_antipara_normal ) {
      epsilons[ surfacePos ] *= 0.5;
      sigmas[ surfacePos ]   *= 0.5;
    }
  }
}



//! this method is for complicated domains
template< class GridType >
void getSigmaEpsilonViaPartialTetTraversal ( const GridType &grid,
                                             const qc::MultiArray< typename GridType::RealType, 3> &deformations,
                                             const typename GridType::RealType aleph, const typename GridType::RealType lambda, const typename GridType::RealType mu,
                                             const aol::Vec3<int> &lowerBnd, const aol::Vec3<int> &upperBnd,
                                             aol::Matrix33<typename GridType::RealType> &AverageSigma, aol::Matrix33<typename GridType::RealType> &AverageEpsilon ) {
  typedef typename GridType::RealType RealType;

  RealType vol = aol::NumberTrait<RealType>::zero;

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    // outside region of interest
    if ( ( el[ 0 ] <  lowerBnd[0] ) || ( el[ 1 ] <  lowerBnd[1] ) || ( el[ 2 ] <  lowerBnd[2] ) ||
         ( el[ 0 ] >= upperBnd[0] ) || ( el[ 1 ] >= upperBnd[1] ) || ( el[ 2 ] >= upperBnd[2] ) ) {
      continue;
    }

    // const CFEType type = el.cfeType();

    el.computeAssembleData ( grid );

    for ( tpcfe::CFETetraInElementIterator<RealType> tit ( el, -1 ); tit.notAtEnd(); ++tit ) {

      const tpcfe::CFETetra<RealType> &tetra = *tit;
      aol::Vec3<RealType> VertexCoords[4];

      for ( short vtx = 0; vtx < 4; vtx++ )
        tetra.computeGlobalCoordinate ( VertexCoords[ vtx ], el, vtx );

      aol::Vec3<RealType> directions[3];
      for ( short i = 0; i < 3; ++i ) {
        directions[i] = aleph * grid.H() * ( VertexCoords[ i ] - VertexCoords[ 3 ] );
      }

      aol::Matrix33<RealType> directionsMatrix;
      for ( short i = 0; i < 3; ++i )
        directionsMatrix.setRow ( i, directions[i] );

      aol::Matrix33<RealType> TetraMatrix = directionsMatrix.inverse();  // unit is m^{-1}

      aol::Vec3<RealType> deformationAtVertex[4];

      for ( short vtx = 0; vtx < 4; ++vtx ) {

        const int
          lIdx0 = tetra( vtx, 0 ),
          lIdx1 = tetra( vtx, 1 ),
          gIdx0 = el.globalIndex(lIdx0);

        if ( lIdx1 == 11 ) { // at regular nodes, we have data
          deformationAtVertex[vtx] = deformations.get ( gIdx0 );
        } else {             // at virtual nodes, we need to interpolate
          const int gIdx1 = el.globalIndex(lIdx1);
          const typename GridType::VNType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
          vn.extrapolate ( deformations, deformationAtVertex[vtx] );  // extrapolation preserves unit
        }
      }

      aol::Vec3<RealType>
        differenceDeformations[3],
        gradientDeformations[3];

      for ( short i = 0; i < 3; ++i )
        for ( short j = 0; j < 3; ++j )
          differenceDeformations[j][i] = deformationAtVertex[ i ][j] - deformationAtVertex[ 3 ][j];

      for ( short i = 0; i < 3; ++i )
        gradientDeformations[i] = TetraMatrix * differenceDeformations[i];

      aol::Matrix33<RealType> epsilon, sigma;

      for ( short i = 0; i < 3; ++i )
        for ( short j = 0; j < 3; ++j )
          epsilon.set ( i, j,  0.5 * ( gradientDeformations[i][j] + gradientDeformations[j][i] )  );

      for ( short i = 0; i < 3; ++i ) {
        for ( short j = 0; j < 3; ++j ) {
          sigma.set ( i, j, mu * ( gradientDeformations[i][j] + gradientDeformations[j][i] ) );
        }
        sigma.add ( i, i,  lambda * ( gradientDeformations[0][0] + gradientDeformations[1][1] + gradientDeformations[2][2] ) );
      }

      RealType Vol = tetra.getVolume() * aol::Cub ( grid.H() );
      vol += Vol;
      sigma *= Vol;
      epsilon *= Vol;

      AverageSigma += sigma;
      AverageEpsilon += epsilon;

    } // end tetra in element loop
  } // end element loop

  cerr << "Volume = " << vol << endl;

  const RealType bboxVolume = aol::Cub ( aleph ) * ( grid.H() * ( upperBnd[0] - lowerBnd[0] ) *
                                                     grid.H() * ( upperBnd[1] - lowerBnd[1] ) *
                                                     grid.H() * ( upperBnd[2] - lowerBnd[2] ) );

  cerr << "bboxVolume = " << bboxVolume << ", vol = " << vol << endl;

  AverageSigma /= bboxVolume; // average over total volume of bounding box
  AverageEpsilon /= vol;      // average over inner volume treated
}

//! this method is for complicated domains
template< class GridType >
void getSigmaEpsilonViaFullTetTraversal ( const GridType &grid,
                                          const qc::MultiArray< typename GridType::RealType, 3> &deformations,
                                          const typename GridType::RealType aleph, const typename GridType::RealType lambda, const typename GridType::RealType mu,
                                          aol::Matrix33<typename GridType::RealType> &AverageSigma, aol::Matrix33<typename GridType::RealType> &AverageEpsilon ) {
  getSigmaEpsilonViaPartialTetTraversal<GridType> ( grid, deformations, aleph, lambda, mu, aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( grid.getNumX() - 1, grid.getNumY() - 1, grid.getNumZ() - 1 ), AverageSigma, AverageEpsilon );
}


//! this method is for jumping coefficients
template< class GridType >
void getSigmaEpsilonViaPartialTetTraversal ( GridType &grid,
                                             qc::MultiArray< typename GridType::RealType, 3> &deformations,
                                             const qc::AArray< typename GridType::NodalCoeffType, qc::QC_3D > &nodalCoeff,
                                             const typename GridType::RealType aleph,
                                             const aol::Vec3<int> &lowerBnd, const aol::Vec3<int> &upperBnd,
                                             aol::Matrix33<typename GridType::RealType> &AverageSigma, aol::Matrix33<typename GridType::RealType> &AverageEpsilon ) {

  typedef typename GridType::RealType RealType;

  RealType vol = aol::NumberTrait<RealType>::zero;

  const qc::GridSize<qc::QC_3D> gridSize ( grid );
  for ( typename GridType::FullElementIterator it ( grid ); it.notAtEnd(); ++it ) {
    tpcfe::CFEElement<RealType> el ( *it, gridSize, grid.getElType ( *it ) );

    // outside region of interest
    if ( ( el[ 0 ] <  lowerBnd[0] ) || ( el[ 1 ] <  lowerBnd[1] ) || ( el[ 2 ] <  lowerBnd[2] ) ||
         ( el[ 0 ] >= upperBnd[0] ) || ( el[ 1 ] >= upperBnd[1] ) || ( el[ 2 ] >= upperBnd[2] ) ) {
      continue;
    }

    // const CFEType type = el.cfeType();
    el.computeAssembleData ( grid );

    tpcfe::CFEWeightProvider < RealType, typename GridType::NodalCoeffType > weightProvider ( nodalCoeff, el );

    for ( CFETetraInElementIterator<RealType> tit ( el ); tit.notAtEnd(); ++tit ) {

      const CFETetra<RealType> &tetra = *tit;
      aol::Vec3<RealType> VertexCoords[4];

      for ( short vtx = 0; vtx < 4; vtx++ )
        tetra.computeGlobalCoordinate ( VertexCoords[ vtx ], el, vtx );

      aol::Vec3<RealType> directions[3];
      for ( short i = 0; i < 3; ++i ) {
        directions[i] = aleph * grid.H() * ( VertexCoords[ i ] - VertexCoords[ 3 ] );
      }

      aol::Matrix33<RealType> directionsMatrix;
      for ( short i = 0; i < 3; ++i )
        directionsMatrix.setRow ( i, directions[i] );

      aol::Matrix33<RealType> TetraMatrix = directionsMatrix.inverse();  // unit is m^{-1}

      aol::Vec3<RealType> deformationAtVertex[4];

      for ( short vtx = 0; vtx < 4; ++vtx ) {

        const int
          lIdx0 = tetra( vtx, 0 ),
          lIdx1 = tetra( vtx, 1 ),
          gIdx0 = el.globalIndex(lIdx0);

        if ( lIdx1 == 11 ) { // at regular nodes, we have data
          deformationAtVertex[vtx] = deformations.get ( gIdx0 );
          // cerr << gIdx0 << " " << deformationAtVertex[vtx] << endl;
        } else {             // at virtual nodes, we need to interpolate
          const int gIdx1 = el.globalIndex(lIdx1);
          const typename GridType::VNType &vn = grid.getVirtualNodeRef ( gIdx0, gIdx1 );
          vn.extrapolate ( deformations, deformationAtVertex[vtx] );  // extrapolation preserves unit
          // cerr << gIdx0 << " " << gIdx1 << " " << deformationAtVertex[vtx] << endl;
        }
      }

      aol::Vec3<RealType>
        differenceDeformations[3],
        gradientDeformations[3];

      for ( short i = 0; i < 3; ++i )
        for ( short j = 0; j < 3; ++j )
          differenceDeformations[j][i] = deformationAtVertex[ i ][j] - deformationAtVertex[ 3 ][j];

      for ( short i = 0; i < 3; ++i )
        gradientDeformations[i] = TetraMatrix * differenceDeformations[i];

      aol::Matrix33<RealType> epsilon, sigma;

      for ( short i = 0; i < 3; ++i )
        for ( short j = 0; j < 3; ++j )
          epsilon.set ( i, j,  0.5 * ( gradientDeformations[i][j] + gradientDeformations[j][i] )  );

      RealType localTensor[3][3][3][3];
      weightProvider.meanWeight ( (*tit).getSign() ).getAnisotropicTensor ( localTensor );

      for ( short i = 0; i < 3; ++i ) {
        for ( short j = 0; j < 3; ++j ) {
          for ( short k = 0; k < 3; ++k ) {
            for ( short l = 0; l < 3; ++l ) {
              sigma.add ( i, j,  localTensor[i][j][k][l] * epsilon.get ( k, l ) );
            }
          }
        }
      }

      RealType Vol = tetra.getVolume() * aol::Cub ( grid.H() );
      vol += Vol;
      sigma *= Vol;
      epsilon *= Vol;

      AverageSigma += sigma;
      AverageEpsilon += epsilon;

    } // end tetra in element loop
  } // end element loop

  cerr << "Volume = " << vol << endl;

  const RealType bboxVolume = aol::Cub ( aleph ) * ( grid.H() * ( upperBnd[0] - lowerBnd[0] ) *
                                                     grid.H() * ( upperBnd[1] - lowerBnd[1] ) *
                                                     grid.H() * ( upperBnd[2] - lowerBnd[2] ) );

  cerr << "bboxVolume = " << bboxVolume << ", vol = " << vol << endl;


  AverageSigma /= bboxVolume ; // volume of bounding box
  AverageEpsilon /= vol;       // divide by volume treated
}

//! this method is for jumping coefficients
template< class GridType >
void getSigmaEpsilonViaFullTetTraversal ( GridType &grid,
                                          qc::MultiArray< typename GridType::RealType, 3> &deformations,
                                          const qc::AArray< typename GridType::NodalCoeffType, qc::QC_3D > &nodalCoeff,
                                          const typename GridType::RealType aleph,
                                          aol::Matrix33<typename GridType::RealType> &AverageSigma, aol::Matrix33<typename GridType::RealType> &AverageEpsilon ) {
  getSigmaEpsilonViaPartialTetTraversal<GridType> ( grid, deformations, nodalCoeff, aleph, aol::Vec3<int> ( 0, 0, 0 ), aol::Vec3<int> ( grid.getNumX() - 1, grid.getNumY() - 1, grid.getNumZ() - 1 ), AverageSigma, AverageEpsilon );
}



enum analysis_face { FACE_FRONT, FACE_BACK };
#ifdef VERBOSE
#define PRINT_FIT_DATA
#endif

template < typename GridType >
typename GridType::RealType computePreNuValueViaOctave ( const GridType &grid, const qc::MultiArray< typename GridType::RealType, qc::QC_3D> &u, const short compression_dir, const short analysis_dir, const analysis_face Face ) {
  char filename[1024], command[1024];
  static int number = 0;

  {
    const int d0 = analysis_dir, d1 = ( analysis_dir + 1 ) % 3, d2 = ( analysis_dir + 2 ) % 3;

    sprintf ( filename, "out/fit%03d.dat", number );
    ofstream fitout ( filename );
    qc::GridDefinition::OldFullBoundaryNodeIterator bnit;
    const typename GridType::RealType h = grid.H();

    aol::Vec3<typename GridType::RealType> normal;
    normal [ analysis_dir ] = ( Face == FACE_FRONT ? -1.0 : 1.0 );

#ifdef PRINT_FIT_DATA
    cerr << "Data for fitting:" << endl;
#endif
    for ( bnit = grid.begin(); bnit != grid.end(); ++bnit ) {
      qc::CoordType pos = *bnit;
      const int fac = ( Face == FACE_FRONT ? 0 : grid.getSize()[analysis_dir] - 1 );
      if ( ( pos[ d0 ] == fac ) && ( aol::Sqr ( h * pos[ d1 ] - 0.5 ) + aol::Sqr ( h * pos[d2] - 0.5 ) < ( 1.0 / 9.0 ) ) && grid.isInsideDomain ( pos ) ) {
        fitout << aol::longScientificFormat ( pos[ d1 ] * h ) << " " << aol::longScientificFormat ( pos[ d2 ] * h ) << " " << aol::longScientificFormat ( u.get ( pos ) * normal ) << endl;
#ifdef PRINT_FIT_DATA
        cerr << aol::longScientificFormat ( pos[ d1 ] * h ) << " " << aol::longScientificFormat ( pos[ d2 ] * h ) << " " << aol::longScientificFormat ( u.get ( pos ) * normal ) << endl;
#endif
      }
    }
    sprintf ( filename, "out/fit%03d.input", number );
    ofstream fitplot ( filename );
    sprintf( command, "fit f(x,y) 'out/fit%03d.dat' using 1:2:3:(1) via a, b, c, m, n", number );
    fitplot << "f ( x, y ) = a * (x-m)*(x-m) + b * (y-n)*(y-n) + c" << endl
    << "a = -0.01; b = -0.01; c = 0.001; m = 0.5; n = 0.5;" << endl
    << command << endl;
  } // close ofstreams

  sprintf ( command, "gnuplot < out/fit%03d.input", number );
  system ( command );
  sprintf ( command, "tail -n 15 fit.log | head -n 5 > out/fit%03d.logtmp", number );
  system ( command );
  sprintf ( filename, "out/fit%03d.logtmp", number );
  ifstream fitlog ( filename );

  aol::Vector<float> coeff ( 5 );

  for ( int i = 0; i < 5; ++i ) {
    char line [1024];
    fitlog.getline ( line, 1024 );
    cerr << "found " << line << endl;
    sscanf ( line, "%*s %*s %f", &coeff[i] );
    cerr << "Coeff[i] = " << coeff[i] << endl;
  }

  sprintf ( filename, "out/plot%03d.input", number );
  ofstream plotout ( filename );
  sprintf ( filename, "out/fit%03d_", number );
  sprintf ( command, "splot \"out/fit%03d.dat\" w p, ", number );
  plotout << "set terminal postscript eps color enhanced" << endl
  << "set output \"" << filename << compression_dir << "_" << analysis_dir << "_" <<  ( Face == FACE_FRONT ? "f" : "b" ) << ".eps\"" << endl
  << "set samples 100,100" << endl
  << "set isosamples 40,40" << endl
  << command << coeff[0] << " * (x-" << coeff[3] << ")*(x-" << coeff[3] << " ) + " << coeff[1] << " * (y-" << coeff[4] << ")*(y-" << coeff[4] << ") + " << coeff[2] << endl;

  sprintf ( command, "gnuplot < out/plot%03d.input", number );
  system ( command );

  cerr << "Function:" << endl;
  cerr << "f ( x, y ) = " << coeff[0] << " * (x-" << coeff[3] << ")*(x-" << coeff[3] << " ) + " << coeff[1] << " * (y-" << coeff[4] << ")*(y-" << coeff[4] << ") + " << coeff[2] << endl;

  number++;

  return ( coeff[2] );
}



// the following method does not really belong here and requires documentation ...
template< typename GridType >
void setVaryingTransverselyIsotropicTensorCoeffs ( const GridType &grid,
                                                   const qc::ScalarArray<typename GridType::RealType, qc::QC_3D> &levelset,
                                                   const int nRods,
                                                   const aol::Mat<6,6,typename GridType::RealType> VoigtTensorMinus[3],
                                                   const aol::Mat<6,6,typename GridType::RealType> &VoigtTensorPlus,
                                                   qc::AArray< typename GridType::NodalCoeffType, qc::QC_3D > &coeff,
                                                   const typename GridType::RealType xRotAngle = 0 ) {
  typedef typename GridType::RealType RealType;

  qc::RectangularContainer < aol::Mat<6,6,RealType>,  qc::AArray< aol::Mat<6,6,RealType>, qc::QC_3D >,  qc::QC_3D > dualCoeff ( aol::Vec3<int> ( -nRods, -nRods, -nRods ), aol::Vec3<int> ( 2 * nRods, 2 * nRods, 2 * nRods ) );

  for ( qc::RectangularIterator<qc::QC_3D> rit ( dualCoeff ); rit.notAtEnd(); ++rit ) {
    if ( ( aol::Abs ( (*rit)[0] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 0 ) ) {
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[0];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[1];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[2];
      dualCoeff.getRef( *rit ) /= 3.0;
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 1 ) ) {
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[0];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[1];
      dualCoeff.getRef( *rit ) /= 2.0;
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 0 ) ) {
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[0];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[2];
      dualCoeff.getRef( *rit ) /= 2.0;
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 0 ) ) {
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[1];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[2];
      dualCoeff.getRef( *rit ) /= 2.0;
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 1 ) ) {
      dualCoeff.getRef( *rit ) = VoigtTensorMinus[0];
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 0 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 1 ) ) {
      dualCoeff.getRef( *rit ) = VoigtTensorMinus[1];
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 0 ) ) {
      dualCoeff.getRef( *rit ) = VoigtTensorMinus[2];
    } else if ( ( aol::Abs ( (*rit)[0] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[1] % 2 ) == 1 ) && ( aol::Abs ( (*rit)[2] % 2 ) == 1 ) ) {
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[0];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[1];
      dualCoeff.getRef( *rit ) += VoigtTensorMinus[2];
      dualCoeff.getRef( *rit ) /= 3.0;
    } else {
      cerr << "This should not happen." << endl;
      abort();
    }
  }

  // set interpolated tensor values
  aol::Matrix33<RealType> rotMat, rotMatScaled;
  rotMat.setRotationAboutX ( xRotAngle );
  rotMatScaled = rotMat;
  rotMatScaled /= cos ( xRotAngle );

  for ( typename GridType::FullNodeIterator fnit ( grid ); fnit.notAtEnd(); ++fnit ) {
    if ( levelset.get ( *fnit ) < 0 ) {
      aol::Mat<6,6,RealType> LocalVoigtTensor, dwarf;

      // poor man's multilinear interpolation

      aol::Vec3<RealType> pos( *fnit );
      pos *= ( grid.H() * ( nRods + 1 ) );

      if ( xRotAngle != 0 ) {
        pos -= aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );
        aol::Vec3<RealType> rotPos = rotMatScaled * pos;
        rotPos += aol::Vec3<RealType> ( 0.5, 0.5, 0.5 );

        pos = rotPos;
      }

      // multilinear interpolation
      aol::Vec3<short> posInt ( static_cast<int>( pos[0] ), static_cast<int>( pos[1] ), static_cast<int>( pos[2] ) ), fromPos;
      aol::Vec3<RealType> posDiff ( - posInt );
      posDiff += pos;

      fromPos = posInt + aol::Vec3<short> ( 0, 0, 0 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= ( 1-posDiff[0] ) * ( 1-posDiff[1] ) * ( 1-posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 0, 0, 1 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= ( 1-posDiff[0] ) * ( 1-posDiff[1] ) * (   posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 0, 1, 0 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= ( 1-posDiff[0] ) * (   posDiff[1] ) * ( 1-posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 1, 0, 0 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= (   posDiff[0] ) * (1- posDiff[1] ) * ( 1-posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 0, 1, 1 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= ( 1-posDiff[0] ) * (   posDiff[1] ) * (   posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 1, 0, 1 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= (   posDiff[0] ) * ( 1-posDiff[1] ) * (   posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 1, 1, 0 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= (   posDiff[0] ) * (   posDiff[1] ) * ( 1-posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      fromPos = posInt + aol::Vec3<short> ( 1, 1, 1 );
      if ( fromPos[0] < dualCoeff.getUpper()[0] && fromPos[1] < dualCoeff.getUpper()[1] && fromPos[2] < dualCoeff.getUpper()[2] ) {
        dwarf = dualCoeff.getRef ( fromPos );
        dwarf *= (   posDiff[0] ) * (   posDiff[1] ) * (   posDiff[2] );
        LocalVoigtTensor += dwarf;
      }

      typename GridType::NodalCoeffType LocalVoigtCoeff ( LocalVoigtTensor );
      if ( xRotAngle != 0 ) {
        LocalVoigtCoeff.rotateTensorBy ( rotMat );
      }

      coeff.getRef ( *fnit ) = LocalVoigtCoeff;

    } else {
      coeff.getRef ( *fnit ) = typename GridType::NodalCoeffType( VoigtTensorPlus );
    }

  }
}


template< typename RealType >
void convertTensor ( const aol::Matrix33<RealType> sigmas[3][3], aol::FullMatrix<RealType> &AllSigmas ) { // the first should be a reference, but we do not really care
  static const int imap[3][3] = { { 0, 5, 4 },
                                  { 8, 1, 3 },
                                  { 7, 6, 2 } };

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      for ( int k = 0; k < 3; ++k ) {
        for ( int l = 0; l < 3; ++l ) {
          AllSigmas.set ( imap[i][j], imap[k][l], sigmas[i][j][k][l] );
        }
      }
    }
  }
}


template< typename RealType >
void averageAndDumpTensor ( const aol::FullMatrix<RealType> &AllSigmas, aol::FullMatrix<RealType> &VoigtSigmas ) {

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      VoigtSigmas.set ( i, j, AllSigmas.get ( i, j ) );
    }
  }

  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 3; j < 6; ++j ) {
      cerr << "Averageing "
           << aol::detailedFormat ( AllSigmas.get ( i, j ) ) << " and " << aol::detailedFormat ( AllSigmas.get ( i, j+3 ) )
           << " for entry " << i << " " << j << endl;
      VoigtSigmas.set ( i, j, ( AllSigmas.get ( i, j ) + AllSigmas.get ( i, j+3 ) ) / 2 );
    }
  }

  for ( int i = 3; i < 6; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      cerr << "Averageing "
           << aol::detailedFormat ( AllSigmas.get ( i, j ) ) << " and " << aol::detailedFormat ( AllSigmas.get ( i+3, j ) )
           << " for entry " << i << " " << j << endl;
      VoigtSigmas.set ( i, j, ( AllSigmas.get ( i, j ) + AllSigmas.get ( i+3, j ) ) / 2 );
    }
  }

  for ( int i = 3; i < 6; ++i ) {
    for ( int j = 3; j < 6; ++j ) {
      cerr << "Averageing "
           << aol::detailedFormat ( AllSigmas.get( i, j ) ) << " " << aol::detailedFormat ( AllSigmas.get( i, j + 3 ) ) << " "
           << aol::detailedFormat ( AllSigmas.get( i + 3, j ) ) << " " << aol::detailedFormat ( AllSigmas.get( i + 3, j + 3 ) )
           << " for entry " << i << " " << j << endl;
      VoigtSigmas.set ( i, j,  ( AllSigmas.get( i, j ) +  AllSigmas.get( i, j+3 ) +  AllSigmas.get( i+3, j ) +  AllSigmas.get( i+3, j+3 ) ) / 4 );
    }
  }

  cerr << AllSigmas << endl << endl << VoigtSigmas << endl << endl;

  cerr << endl << endl << "finally: the TENSOR:" << endl;
  for ( int i = 0; i < 6; ++i ) {
    for ( int j = 0; j < 6; ++j ) {
      cerr << aol::detailedFormat ( VoigtSigmas.get( i, j ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  aol::FullMatrix<RealType> VoigtSigmasSymmetrized ( 6, 6 );
  for ( int i = 0; i < 6; ++i ) {
    for ( int j = 0; j < 6; ++j ) {
      VoigtSigmasSymmetrized.set ( i, j, 0.5 * ( VoigtSigmas.get( i, j ) + VoigtSigmas.get( j, i ) ) );
    }
  }

  cerr << "and its symmetrized version: " << endl;
  for ( int i = 0; i < 6; ++i ) {
    for ( int j = 0; j < 6; ++j ) {
      cerr << aol::detailedFormat ( VoigtSigmasSymmetrized.get( i, j ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

}


template< typename RealType >
void convertAverageAndDumpSigmaTensor ( const aol::Matrix33<RealType> sigmas[3][3] ) { // this should be a reference, but we do not really care

  aol::FullMatrix<RealType> AllSigmas( 9, 9 );

  convertTensor ( sigmas, AllSigmas );

  aol::FullMatrix<RealType> VoigtSigmas( 6, 6 );

  averageAndDumpTensor ( AllSigmas, VoigtSigmas );
}


template< typename RealType >
void convertAverageAndDumpTensors ( const aol::Matrix33<RealType> sigmas[3][3], const aol::Matrix33<RealType> epsilons[3][3] ) { // these should be references, but we do not really care

  aol::FullMatrix<RealType> AllSigmas( 9, 9 ), AllEpsilons ( 9, 9 );

  convertTensor ( epsilons, AllEpsilons );
  convertTensor ( sigmas, AllSigmas );

  aol::FullMatrix<RealType> VoigtSigmas( 6, 6 ), VoigtEpsilons ( 6, 6 );

  averageAndDumpTensor ( AllEpsilons, VoigtEpsilons );
  averageAndDumpTensor ( AllSigmas, VoigtSigmas );
}
}

#endif
