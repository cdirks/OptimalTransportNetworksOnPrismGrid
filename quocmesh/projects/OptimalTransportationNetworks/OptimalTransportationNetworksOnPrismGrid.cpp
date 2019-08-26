#include<aol.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<initializer_list>
#include<scalarArray.h>
#include<FEOpInterface.h>
#include<Newton.h>
#include<omp.h>
#include<smallMat.h>
#include<parameterParser.h>
#include<unordered_map>
#include<chrono>
#include<sys/stat.h>

#include<Timer.h>
#include<Example.h>
#include<adaptiveFEPrismMesh.h>
#include<adaptiveFEPrismMeshConfigurator.h>
#include<adaptiveFEPrismMeshRPOp.h>
#include<SemiMassMatrixAssembler.h>
#include<MixedGaussFDFEQuadrature3D.h>
#include<SpecialAdaptiveFEPrismMeshRPOp.h>
#include<PDMinimizationSolver.h>
#include<adaptivePrismMeshValueMap.h>
#include<AuxFunctions.h>



typedef double RealType;


//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
// Main 
int main () {

  // Create results folder
  //string results = "Results"; 
  mkdir("Results",0777);  
  #pragma GCC diagnostic push 
  #pragma GCC diagnostic ignored "-Wunused-result"
  system("rm -r Results/*");
  #pragma GCC diagnostic pop
  
  Timer timerAll;
  timerAll.reset();

  //------------------------------------
  // Get parameters
  aol::ParameterParser parser ( "../../../quocmesh/projects/OptimalTransportationNetworks/Parameters.par" );
  const string problemType = parser.getString ( "problemtype" );
  const int INITIALGRIDLEVELXY = parser.getInt ( "initialgridlevelxy" );
  const int INITIALGRIDLEVELZ = parser.getInt ( "initialgridlevelz" );
  const int NUMREFINEMENTS = parser.getInt ( "numrefinements" );
  const int MAXIMALGRIDLEVELXY = parser.getInt ( "maximalgridlevelxy" );
  const int MAXIMALGRIDLEVELZ = parser.getInt ( "maximalgridlevelz" );
  int maxIter = parser.getInt ( "maxiter" );
  int maxIterFinalRun = parser.getInt ( "maxiterfinalrun" );
  const int maxIterProject = parser.getInt ( "maxiterproject" );
  const RealType tol = parser.getReal < RealType > ( "tol" );  
  const RealType tolProject = parser.getReal < RealType > ( "tolproject" );
  const RealType epsilon = parser.getReal < RealType > ( "epsilon" );
  const RealType a = parser.getReal < RealType > ( "a" );
  const RealType tauFactor = parser.getReal < RealType > ( "taufactor" );
  const RealType refinementTol = parser.getReal < RealType > ( "refinementtol" );
  const string exampleType = parser.getString ( "example" );
  const string stoppingCriterion = parser.getString ( "stoppingcriterion" );
  const string refinementCriterion = parser.getString ( "refinementcriterion" );
  const int precon = parser.getInt ( "preconditioning" );
  const int MAXLEVELXY = MAXIMALGRIDLEVELXY + 1;
  const int MAXLEVELZ = MAXIMALGRIDLEVELZ + 1;
  cout << "ooooooooooooooooooooooooooooooooo" << endl;
  if ( problemType == "BranchedTransport" ) cout << "BRANCHED TRANSPORT: epsilon = " << epsilon << endl;
  if ( problemType == "UrbanPlanning" ) cout << "URBAN PLANNING: a = " << a << ", epsilon = " << epsilon << endl;
  cout << "ooooooooooooooooooooooooooooooooo" << endl;  
  //------------------------------------
  typedef MixedGaussFDFEQuadrature3D < RealType > QuadRule;
  typedef AdaptiveFEPrismMesh < RealType > GridType; 
  typedef AdaptiveFEPrismMeshConfigurator < RealType, GridType, QuadRule > ConfiguratorType;
  typedef SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > RPOpType; 
  typedef aol::SparseMatrix < RealType > MatrixType;
  typedef ConfiguratorType::RealType RealType;
  typedef aol::Vector < RealType > VectorType;
  typedef aol::CSRMatrix < RealType > CSRMatrixType;
  typedef aol::TriangleIntegration < RealType, qc::QC_2D, 1 > QuadRule2D; 
  typedef AdaptiveFETriangMesh < RealType > GridType2D; 
  typedef aol::TriangMeshConfigurator < RealType, GridType2D, QuadRule2D > ConfiguratorType2D;
  typedef Example < ConfiguratorType2D, ConfiguratorType > ExampleType; 
  //------------------------------------
  // Get grid  
  GridType grid ( INITIALGRIDLEVELXY, INITIALGRIDLEVELZ );
  //------------------------------------
  // Initialize variables
  VectorType v, phi1, phi2, phi3, vLastRun, phi1LastRun, phi2LastRun, phi3LastRun;
  int gridDofs, gridNodes, numHangingNodes, valueMapSize = 0;
  //------------------------------------
  // Get example
  ExampleType example ( exampleType, INITIALGRIDLEVELXY, problemType ); 
  //------------------------------------
  // Set power lookuptable 
  VectorType powerLookuptable;
  if ( problemType == "BranchedTransport" ) setPowerLookuptable ( powerLookuptable, MAXLEVELZ, epsilon ); 
  //------------------------------------
  // Define interpolation helper 
  AdaptivePrismMeshValueMap < RealType > interpolator;
  //------------------------------------
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Repeat iteration 
  bool STOP = false;
  for ( int run = 0; run <= NUMREFINEMENTS; ++run ) {
  
    cout << "---------------> Iteration round " << run << " started." << endl;
  
    // Set configurator and rp for current grid 
    ConfiguratorType configurator ( grid );   
    RPOpType rp ( grid, configurator, MAXLEVELXY, MAXLEVELZ ); 
    gridDofs = grid.getNumberOfDofs();
    gridNodes = grid.getNumberOfNodes(); 
    numHangingNodes = grid.getNumberOfHangingNodes();
      
    // Set finite element matrix assembler 
    SemiMassMatrixAssemblerX < ConfiguratorType > xAssembler ( configurator, grid ); 
    SemiMassMatrixAssemblerY < ConfiguratorType > yAssembler ( configurator, grid );
    SemiMassMatrixAssemblerZ < ConfiguratorType > zAssembler ( configurator, grid );
      
    // Set primal dual solver 
    PDMinimizationSolver < ConfiguratorType, CSRMatrixType, ExampleType > solver ( rp, tol, maxIterProject, tolProject, epsilon, a, tauFactor, stoppingCriterion, precon );  
      
    MatrixType matrixX ( gridNodes, gridNodes ); MatrixType matrixY ( gridNodes, gridNodes ); MatrixType matrixZ ( gridNodes, gridNodes );
    MatrixType matrixXT ( gridDofs, gridDofs ); MatrixType matrixYT ( gridDofs, gridDofs ); MatrixType matrixZT ( gridDofs, gridDofs );

    if ( run == 0 ) {
        
      // Set initial variables 
      v.resize ( gridDofs );
      example.setInitialImage ( rp, v ); example.setBoundaryMask ( rp, run ); 
      phi1.resize ( gridDofs ); phi2.resize ( gridDofs ); phi3.resize ( gridDofs );
      
      // Set lookuptables
      rp.setLineSegmentList(); 
      rp.setConstraintList(); 

      // Define finite element matrices
      xAssembler.assembleAddMatrix ( matrixX ); yAssembler.assembleAddMatrix ( matrixY ); zAssembler.assembleAddMatrix ( matrixZ );
      matrixX.transposeTo ( matrixXT ); matrixY.transposeTo ( matrixYT ); matrixZ.transposeTo ( matrixZT );
      
      // Set step size
      if ( precon == 0 ) {
        solver.setSigma ( matrixX, matrixY, matrixZ );   
        solver.setTau ( matrixXT, matrixYT, matrixZT );   
      }
      
    }  
    else {
        
      // Resize variables
      v.resize ( gridNodes );
      phi1.resize ( gridNodes ); phi2.resize ( gridNodes ); phi3.resize ( gridNodes );
      example.setBoundaryMask ( rp, run );   
      
      // Interpolate results
      MatrixType interpolMat ( gridNodes, valueMapSize );   
      interpolator.makeInterpolationMatrixZConst ( rp, interpolMat );      
      interpolMat.apply ( vLastRun, v );    
      interpolMat.apply ( phi1LastRun, phi1 ); interpolMat.apply ( phi2LastRun, phi2 ); interpolMat.apply ( phi3LastRun, phi3 );
      rp.restrictVectorZConst ( v ); 
      rp.restrictVectorZConst ( phi1 ); rp.restrictVectorZConst ( phi2 ); rp.restrictVectorZConst ( phi3 ); 
      solver.renewBoundaries ( v, example );  
      
      // Set line segment list and constraint list
      rp.setLineSegmentList();     
      if ( problemType == "BranchedTransport" ) {
          rp.updateConstraintListBT ( powerLookuptable, phi1, phi2, 1e-3, true );
          solver.projectionOntoIntegralConstraintBT ( powerLookuptable, rp, phi1, phi2 ); 
      }
      if ( problemType == "UrbanPlanning" ) {
          rp.updateConstraintListUP ( epsilon, a, phi1, phi2, 1e-3, true );           
          solver.projectionOntoIntegralConstraintUP ( rp, phi1, phi2 );   
      }
      
      // Define finite element matrices
      xAssembler.assembleAddMatrix ( matrixX ); yAssembler.assembleAddMatrix ( matrixY ); zAssembler.assembleAddMatrix ( matrixZ );
      if ( numHangingNodes > 0 ) { 
        rp.applyPTPInPlaceZConstSparse ( matrixX, matrixY, matrixZ ); 
        rp.fastRestrictSparseMatrix ( matrixX ); rp.fastRestrictSparseMatrix ( matrixY ); rp.fastRestrictSparseMatrix ( matrixZ );
      }
      matrixX.transposeTo ( matrixXT ); matrixY.transposeTo ( matrixYT ); matrixZ.transposeTo ( matrixZT );
      
      // Set step size
      if ( precon == 0 ) {
        solver.setSigma ( matrixX, matrixY, matrixZ );   
        solver.setTau ( matrixXT, matrixYT, matrixZT );
      }
      
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Iteration
    if ( run == NUMREFINEMENTS ) maxIter = maxIterFinalRun;
    cout << "Start primal dual algorithm." << endl;
    RealType lastError;
    Timer timer;
    timer.reset();
    lastError = solver.minimize ( example, maxIter, powerLookuptable, grid, rp, matrixX, matrixY, matrixZ, matrixXT, matrixYT, matrixZT, v, phi1, phi2, phi3, run ); 
    timer.end();
    cout << "Elapsed time for iteration in round " << run << " : " << timer.elapsed() << " sec." << endl;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Save result
    saveImageAsVTK ( grid, rp, v, "v", run );
    saveImageAsVTK ( grid, rp, phi1, "phi1", run ); saveImageAsVTK ( grid, rp, phi2, "phi2", run ); saveImageAsVTK ( grid, rp, phi3, "phi3", run );

    /*
    // Create image from vector field (for vector field construction)
    if ( run == NUMREFINEMENTS ) {
      createImageFromVectorFieldLaplace ( grid, phi1, phi2, INITIALGRIDLEVELXY );
      saveImageAsVTKVectorField ( grid, rp, phi1, phi2, phi3, "phiVectorField", run );
    }
    */
    
    // Set lastIteration variables
    vLastRun.resize ( gridNodes );
    phi1LastRun.resize ( gridNodes ); phi2LastRun.resize ( gridNodes ); phi3LastRun.resize ( gridNodes );
    rp.extendVectorZConst ( v, vLastRun );
    rp.extendVectorZConst ( phi1, phi1LastRun ); rp.extendVectorZConst ( phi2, phi2LastRun ); rp.extendVectorZConst ( phi3, phi3LastRun );
    
    // Check if all constraints are satisfied 
    int numNotSatisfied = 0;
    if ( problemType == "BranchedTransport" ) numNotSatisfied = solver.checkConstraintUniformBT ( rp, powerLookuptable, phi1, phi2 ); 
    if ( problemType == "UrbanPlanning" ) numNotSatisfied = solver.checkConstraintUniformUP ( rp, phi1, phi2 ); 
    cout << "Run " << run << ": Constraints not satisfied: " << numNotSatisfied << endl;
    
    //oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    // REFINE GRID
    int numMarked;    
    VectorType phi1Diff ( gridDofs ); VectorType phi2Diff ( gridDofs ); VectorType phi3Diff ( gridDofs );
    
    // Threshold v and defift binary v 
    VectorType vThresh ( v );
    vThresh.threshold ( 0.5, 0.0, 1.0 );
    VectorType solution; 
    example.deliftImage ( rp, vThresh, solution, run );

    if ( run < NUMREFINEMENTS ) {
    
        // If error is too large: Emergency break 
        if ( lastError > 1e10 ) {
            cout << "EMERGENCY BREAK!!!" << endl;
            phi1.setZero();
            phi2.setZero();
            phi3.setZero();   
        }
        
        // If algorithm converged: Normal refinement, otherwise refine additionally with z gradient 
        bool strong = false;        
        if ( lastError > 1e-3 ) strong = true;
        valueMapSize = interpolator.getIndexFromNodes ( rp );
        if ( refinementCriterion == "LocalPrimalDualGap" ) {
          VectorType PDgapPrimal ( rp.getGridRef().getNumberOfElements() );
          VectorType PDgapDual ( rp.getGridRef().getNumberOfElements() );
          numMarked = solver.computeLocalPrimalDualGap ( problemType, PDgapPrimal, PDgapDual, example, powerLookuptable, grid, rp, matrixX, matrixY, matrixZ, matrixXT, matrixYT, matrixZT, v, phi1, phi2, phi3, refinementTol, run, strong );   
        }
        else if ( refinementCriterion == "LocalPrimalDualGapOnly" ) {
          VectorType PDgapPrimal ( rp.getGridRef().getNumberOfElements() );
          VectorType PDgapDual ( rp.getGridRef().getNumberOfElements() );
          numMarked = solver.computeLocalPrimalDualGapOnly ( problemType, PDgapPrimal, PDgapDual, example, powerLookuptable, grid, rp, v, phi1, phi2, phi3, refinementTol );      
        }
        else if ( refinementCriterion == "LocalGradient" ) {
          numMarked = solver.computeLocalGradient ( grid, rp, v, refinementTol );
        }
        
        if ( numMarked == 0 ) STOP = true;
        
    }
    else STOP = true; 
    
    // If no more elements were refined or maxNumRefinements is reached: Stop iteration
    if ( STOP ) break;
    
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  timerAll.end();
  cout << "Overall time: " << timerAll.elapsed() << " sec." << endl;
  RealType numElementsUniform = static_cast<RealType>( 2 * (1<<MAXIMALGRIDLEVELXY) * (1<<MAXIMALGRIDLEVELXY) * (1<<MAXIMALGRIDLEVELZ) );
  RealType numNodesUniform = static_cast<RealType>( ((1<<MAXIMALGRIDLEVELXY)+1) * ((1<<MAXIMALGRIDLEVELXY)+1) * ((1<<MAXIMALGRIDLEVELZ)+1) );
  cout << "Total number of elements in final round: " << grid.getNumberOfElements() << " (" << 100.0*static_cast<RealType>(grid.getNumberOfElements())/numElementsUniform << " % of uniform resolution grid)" << endl;
  cout << "Total number of nodes in final round: " << grid.getNumberOfNodes() << " (" << 100.0*static_cast<RealType>(grid.getNumberOfNodes())/numNodesUniform << " % of uniform resolution grid)" << endl;
  
  return 0;
  
}

