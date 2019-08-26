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

#include <multiArray.h>
#include <parameterParser.h>
#include <tpCFEElastOp.h>
#include <tpCFELevelsets.h>
#include <tpCFEMultigrid.h>
#include <tpCFEUtils.h>
#include <shapeLevelsetGenerator.h>


typedef double RealType;
static const tpcfe::ConstraintType CT = tpcfe::CFE_CDWI_ELAST;
//typedef tpcfe::CFEGrid < RealType, CT, tpcfe::IsotropicElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEGrid < RealType, CT, tpcfe::VoigtElasticityCoefficient<RealType> > GridType;
typedef tpcfe::CFEHybridMatrix<GridType> MatrixType;
typedef tpcfe::CFEConfigurator < GridType, MatrixType > ConfiguratorType;

typedef tpcfe::CFEJCElastOp< ConfiguratorType > ElastOpType;


int main ( const int, const char**  ) {

  try {

    aol::StopWatch timer;

    const int level = 5;
    GridType grid ( level );

    qc::ScalarArray<RealType, qc::QC_3D> lsDomain ( grid ), lsSpongiosa ( grid );

    const RealType
      outerRadX = 0.45,
      outerRadY = 0.35,
      innerRadX = 0.39,
      innerRadY = 0.29;

    qc::ShapeLevelsetGenerator<RealType>::generateEllipticColumnLevelset ( lsDomain,    outerRadX, outerRadY );
    qc::ShapeLevelsetGenerator<RealType>::generateEllipticColumnLevelset ( lsSpongiosa, innerRadX, innerRadY );

    grid.setDomainFrom ( lsDomain );
    grid.addStructureFrom ( lsSpongiosa ); // if both domain and structure are to be set, domain must be set first.

    grid.setCoarsenedFlag ( true ); // so that the hybrid matrix works properly.

#if 0
    for ( typename GridType::FullElementIterator elit ( grid ); elit.notAtEnd(); ++elit ) {
      if ( (*elit)[0] == grid.getNumX() / 2 && (*elit)[2] == grid.getNumZ() / 2) {
        cerr << *elit << " " << static_cast<int> ( grid.getElType ( *elit )._structureNo ) << " " << static_cast<int> ( grid.getElType ( *elit )._pureType ) << endl;
      }
    }
    // abort();
#endif

    qc::AArray< GridType::NodalCoeffType, qc::QC_3D > coeff ( grid );

#if 0
    const RealType
      EMinus  = 2.0, // use something like the spongiosa tensor
      nuMinus = 0.38,
      EPlus   = 13.0, // use compacta tensor
      nuPlus  = 0.35;

    cerr << EMinus << " " << nuMinus << " " << EPlus << " " << nuPlus << endl;
    const GridType::NodalCoeffType ENuMinus ( EMinus, nuMinus ), ENuPlus ( EPlus, nuPlus );
    tpcfe::setCoeffForLevelset ( coeff, lsSpongiosa, ENuMinus, ENuPlus );
#else

    aol::Mat<6,6,RealType> VoigtTensorSpongiosaMiddle, VoigtTensorSpongiosaBoundary, VoigtTensorCompacta;

    {
      const RealType
        E  = 13.0,
        nu = 0.32,
        Gmix = E / ( 2 * ( 1 + nu ) );

      aol::Mat<6,6,RealType> VoigtTensorInverse;

      VoigtTensorInverse[0][0] = 1.  / E   ;   VoigtTensorInverse[0][1] = -nu / E   ;   VoigtTensorInverse[0][2] = -nu / E   ;
      VoigtTensorInverse[1][0] = -nu / E   ;   VoigtTensorInverse[1][1] = 1   / E   ;   VoigtTensorInverse[1][2] = -nu / E   ;
      VoigtTensorInverse[2][0] = -nu / E   ;   VoigtTensorInverse[2][1] = -nu / E   ;   VoigtTensorInverse[2][2] = 1   / E   ;
      VoigtTensorInverse[3][3] = 1.  / Gmix;
      VoigtTensorInverse[4][4] = 1.  / Gmix;
      VoigtTensorInverse[5][5] = 1.  / Gmix;

      VoigtTensorCompacta.makeInverse ( VoigtTensorInverse );

      cerr << "VoigtTensorCompacta" << endl << VoigtTensorCompacta << endl << endl;
    }

    {
      aol::Mat<6,6,RealType> &vt = VoigtTensorSpongiosaBoundary; // 0400_0800_0450
      vt[0][0] =  0.950989295280;   vt[0][1] =  0.354433368999;   vt[0][2] =  0.338035685375;   vt[0][3] = -0.015816983336;   vt[0][4] =  0.018706780716;   vt[0][5] =  0.040829802527;
      vt[1][0] =  0.354433368999;   vt[1][1] =  0.963793581034;   vt[1][2] =  0.329892340833;   vt[1][3] = -0.041330967114;   vt[1][4] =  0.026974899395;   vt[1][5] = -0.053671224538;
      vt[2][0] =  0.338035685375;   vt[2][1] =  0.329892340833;   vt[2][2] =  1.731158996245;   vt[2][3] =  0.059518009600;   vt[2][4] =  0.060780128290;   vt[2][5] =  0.010453529553;
      vt[3][0] = -0.015816983336;   vt[3][1] = -0.041330967114;   vt[3][2] =  0.059518009600;   vt[3][3] =  0.386777671345;   vt[3][4] = -0.005495389659;   vt[3][5] =  0.022630274088;
      vt[4][0] =  0.018706780716;   vt[4][1] =  0.026974899395;   vt[4][2] =  0.060780128290;   vt[4][3] = -0.005495389659;   vt[4][4] =  0.402203799543;   vt[4][5] =  0.000289253999;
      vt[5][0] =  0.040829802527;   vt[5][1] = -0.053671224538;   vt[5][2] =  0.010453529553;   vt[5][3] =  0.022630274088;   vt[5][4] =  0.000289253999;   vt[5][5] =  0.332245222527;
      cerr << "VoigtTensorSpongiosaBoundary" << endl << VoigtTensorSpongiosaBoundary << endl << endl;
    }

    {
      aol::Mat<6,6,RealType> &vt = VoigtTensorSpongiosaMiddle; // 0800_0600_0450
      vt[0][0] =  0.367089282167;   vt[0][1] =  0.142321867856;   vt[0][2] =  0.178939478737;   vt[0][3] = -0.012769168907;   vt[0][4] =  0.009366334168;   vt[0][5] = -0.001421463851;
      vt[1][0] =  0.142321867856;   vt[1][1] =  0.479144529612;   vt[1][2] =  0.183489372411;   vt[1][3] = -0.025595494896;   vt[1][4] = -0.008276986035;   vt[1][5] =  0.013439979443;
      vt[2][0] =  0.178939478737;   vt[2][1] =  0.183489372411;   vt[2][2] =  1.380548553295;   vt[2][3] =  0.026997232253;   vt[2][4] =  0.045736275258;   vt[2][5] =  0.002927310558;
      vt[3][0] = -0.012769168907;   vt[3][1] = -0.025595494896;   vt[3][2] =  0.026997232253;   vt[3][3] =  0.211675465462;   vt[3][4] =  0.004821881749;   vt[3][5] =  0.000253858815;
      vt[4][0] =  0.009366334168;   vt[4][1] = -0.008276986035;   vt[4][2] =  0.045736275258;   vt[4][3] =  0.004821881749;   vt[4][4] =  0.193078825136;   vt[4][5] = -0.009867626549;
      vt[5][0] = -0.001421463851;   vt[5][1] =  0.013439979443;   vt[5][2] =  0.002927310558;   vt[5][3] =  0.000253858815;   vt[5][4] = -0.009867626549;   vt[5][5] =  0.110468654468;
      cerr << "VoigtTensorSpongiosaMiddle" << endl << VoigtTensorSpongiosaMiddle << endl << endl;
    }

    for ( qc::RectangularIterator<qc::QC_3D> rit ( grid ); rit.notAtEnd(); ++rit ) {
      if ( lsDomain.get ( *rit ) > 0 ) {
        coeff.getRef ( *rit ) = GridType::NodalCoeffType ( aol::Mat<6,6,RealType> () );
      } else {
        if ( lsSpongiosa.get ( *rit ) > 0 ) {
          coeff.getRef ( *rit ) = GridType::NodalCoeffType ( VoigtTensorCompacta );
        } else {

          const RealType lambda = sqrt ( aol::Sqr ( grid.H() * ( (*rit)[0] - grid.getNumX()/2 ) / innerRadX ) + aol::Sqr ( grid.H() * ( (*rit)[1] - grid.getNumY()/2 ) / innerRadY ) );
          aol::Mat<6,6,RealType> tmp1 = VoigtTensorSpongiosaBoundary, tmp2 = VoigtTensorSpongiosaMiddle;
          tmp1 *= lambda;
          tmp2 *= ( 1.0 - lambda );
          tmp1 += tmp2;
          coeff.getRef ( *rit ) = GridType::NodalCoeffType ( tmp1 );

          // if ( (*rit)[0] == grid.getNumX() / 2 && (*rit)[2] == grid.getNumZ() / 2) {
          //  cerr << *rit << " " << lambda << " " << endl << tmp1 << endl;
          // }
        }
      }
    }
#endif


    grid.relaxedDetectAndInitVirtualNodes ( coeff, 1.0e-14, 1.0e-14 );

    qc::MultiArray< RealType, qc::QC_3D > dirichletBCs ( grid ), rhs ( grid ), soln ( grid ), dwarf ( grid );

    qc::BitArray<qc::QC_3D> DirichletMask ( grid );

    tpcfe::setZShift ( grid, dirichletBCs, DirichletMask, -0.2 );

    grid.setDirichletMask ( DirichletMask );
    grid.setDOFMaskFromDirichletAndDomainNodeMask();

    timer.start();
    ElastOpType elastOp ( grid, coeff );
    timer.stop();
    cerr << "Setting up elasticity matrix took " << timer.elapsedWallClockTime() << " seconds" << endl;

    dwarf -= dirichletBCs;

    elastOp.apply ( dwarf, rhs );

    elastOp.restrictDirichletEntries();
    grid.restrictDirichletNodes ( rhs );

    timer.start();

    aol::SSORPreconditioner< aol::MultiVector<RealType>, MatrixType > prec ( elastOp.getBlockMatrixRef() );
    aol::PCGInverse< aol::MultiVector<RealType> > solver ( elastOp.getBlockMatrixRef(), prec, 1.0e-16, 200000 );
    solver.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );

    solver.apply ( rhs, soln );

    cerr << "Memusage = " << ( aol::memusage() >> 20 ) << " MiB" << endl;

    timer.stop();

    cerr << "Solver (including setup) took " << timer.elapsedWallClockTime() << " seconds" << endl;

    soln += dirichletBCs;

    soln.save ( "soln.dat.bz2", qc::PGM_DOUBLE_BINARY );

    aol::TriangMesh<float> tMesh;

    { // semi-useful
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );

      itg.determineInterfaceTriangulation();
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "InterfaceU.pov" );
      // tMesh.saveAsUDPLY ( "InterfaceU.ply.bz2" );

      itg.deformVertices ( soln, 1.0 );
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "InterfaceD.pov" );
      tMesh.saveAsUDPLY ( "InterfaceD.ply.bz2" );
    }

    {
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_X );
      itg.determineSliceTriangulation ( qc::QC_Y, grid.getNumY() / 2 );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "SliceU.ply.bz2" );

      itg.deformVertices ( soln, 1.0 );
      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "SliceD.ply.bz2" );
    }

    {
      tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > pseudoGrid ( level );
      pseudoGrid.setDomainFrom ( lsDomain );
      pseudoGrid.detectAndInitVirtualNodes();

      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator< tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > > itg ( pseudoGrid, soln, qc::QC_Z );
      itg.determineInterfaceAndBoundaryTriangulation();

      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "DomainU.ply.bz2" );
    }

    {
      tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > pseudoGrid ( level );
      pseudoGrid.setDomainFrom ( lsSpongiosa );
      pseudoGrid.detectAndInitVirtualNodes();

      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator< tpcfe::CFEGrid < RealType, tpcfe::CFE_CD > > itg ( pseudoGrid, soln, qc::QC_Z );
      itg.determineInterfaceAndBoundaryTriangulation();

      itg.writeToTriangMesh ( tMesh );
      tMesh.saveAsUDPLY ( "SpongiosaU.ply.bz2" );
    }

    if ( false ) { // not useful
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineBoundaryTriangulation( +1 );
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "BoundaryPU.pov" );
      tMesh.saveAsUDPLY ( "BoundaryPU.ply.bz2" );

      itg.deformVertices ( soln, 1.0 );
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "BoundaryPD.pov" );
      tMesh.saveAsUDPLY ( "BoundaryPD.ply.bz2" );
    }

    if ( false ) { // not useful
      tpcfe::CFEInterfaceTriangulationWithComponentDataGenerator<GridType> itg ( grid, soln, qc::QC_Z );
      itg.determineBoundaryTriangulation( -1 );
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "BoundaryMU.pov" );
      tMesh.saveAsUDPLY ( "BoundaryMU.ply.bz2" );

      itg.deformVertices ( soln, 1.0 );
      itg.writeToTriangMesh ( tMesh );
      // tMesh.saveAsPov ( "BoundaryMD.pov" );
      tMesh.saveAsUDPLY ( "BoundaryMD.ply.bz2" );
    }

  } catch ( aol::Exception &exc ) {
    exc.dump();
    return ( EXIT_FAILURE );
  }
  return ( EXIT_SUCCESS );
}
