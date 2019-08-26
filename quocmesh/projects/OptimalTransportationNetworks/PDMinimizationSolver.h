#ifndef __PDMINIMIZATIONSOLVER_H
#define __PDMINIMIZATIONSOLVER_H

#include<aol.h>
#include<MixedGaussFDFEQuadrature3D.h>
#include<SemiMassMatrixAssembler.h>
#include<adaptivePrismMeshValueMap.h>
#include<AuxFunctions.h>

// MOSEK
#ifdef USE_MOSEK
#include<fusion.h>
#include<monty.h>
using namespace mosek::fusion;
using namespace monty;
#endif

template < typename ConfiguratorType, typename MatrixType, typename ExampleType > 
class PDMinimizationSolver {

public:
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef aol::Vector < RealType > VectorType;
    typedef SpecialAdaptiveFEPrismMeshRPOp < ConfiguratorType > RPOpType;
    typedef typename ConfiguratorType::InitType GridType;
    typedef aol::SparseMatrix < RealType > SparseMatrixType;
    typedef MixedGaussFDFEQuadrature3D < RealType > QuadRule;
    typedef typename AdaptiveFEPrismMeshRPOp < ConfiguratorType >::IndexMapType IndexMapType;
    
    PDMinimizationSolver ( RPOpType & rp, const RealType tol, const int maxIterProject, const RealType tolProject, const RealType epsilon, const RealType a, const RealType tauFactor, const string stoppingCrit, const int precon ) {
        _tol = tol;
        _maxIterProject = maxIterProject;
        _tolProject = tolProject;
        _size = rp.getGridRef().getNumberOfDofs();
        _epsilon = epsilon;
        _a = a;
        _tauFactor = tauFactor;
        _stoppingCrit = stoppingCrit;
        _precon = precon; 
    }
    
    
public: 
    
    RealType _tol;
    int _maxIterProject;
    RealType _tolProject;
    RealType _sigma;
    RealType _tau;
    int _size;
    RealType _epsilon;
    RealType _a;
    RealType _tauFactor;
    string _stoppingCrit;
    int _precon; 
    
    RealType getSigma () const {
        return _sigma;
    }
    
    RealType getTau () const {
        return _tau;
    }
    
    void setSigma ( SparseMatrixType &matrixX, SparseMatrixType &matrixY, SparseMatrixType &matrixZ ) { 
        RealType maxRowSumX, maxRowSumY, maxRowSumZ, sum;
        int numRows = matrixX.getNumRows();
        maxRowSumX = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixX, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumX ) {
                maxRowSumX = sum;
            }
        }
        maxRowSumY = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixY, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumY ) {
                maxRowSumY = sum;
            }
        }
        maxRowSumZ = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixZ, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumZ ) {
                maxRowSumZ = sum;
            }
        }
        _sigma = 1.0 / ( aol::Max ( maxRowSumX, maxRowSumY, maxRowSumZ ));
        
    }
     
    void setTau ( SparseMatrixType &matrixX, SparseMatrixType &matrixY, SparseMatrixType &matrixZ ) { 
        RealType maxRowSumX, maxRowSumY, maxRowSumZ, sum;
        int numRows = matrixX.getNumRows();
        maxRowSumX = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixX, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumX ) {
                maxRowSumX = sum;
            }
        }
        maxRowSumY = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixY, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumY ) {
                maxRowSumY = sum;
            }
        }
        
        maxRowSumZ = 0.0;
        for ( int i = 0; i < numRows; ++i ) {
            sum = 0.0;
            for ( aol::SparseMatrixRowIterator < RealType > rowit ( matrixZ, i ); rowit.notAtEnd(); ++rowit ) {
                sum += aol::Abs ( rowit->value );
            }
            if ( sum > maxRowSumZ ) {
                maxRowSumZ = sum;
            }
            
        }
        _tau = 1.0 / ( aol::Max ( maxRowSumX, maxRowSumY, maxRowSumZ ));
        
    }
    
    void setSigmaTauSpectralNorm ( SparseMatrixType &matrixX, SparseMatrixType &matrixY, SparseMatrixType &matrixZ, SparseMatrixType &matrixXT, SparseMatrixType &matrixYT, SparseMatrixType &matrixZT ) {
      int gridDofs = matrixX.getNumRows();
      SparseMatrixType MTM ( gridDofs, gridDofs );
      MTM.addMatrixProduct ( matrixXT, matrixX );
      MTM.addMatrixProduct ( matrixYT, matrixY );
      MTM.addMatrixProduct ( matrixZT, matrixZ );
      RealType maxEigenValue = MTM.getMaxEigenValue(); 
      _sigma = 1.0/sqrt(aol::Abs(maxEigenValue)) - 1e-3;
      _tau = _sigma;  
    }  
    
    void setSigmaTauFrobeniusNorm ( SparseMatrixType &matrixX, SparseMatrixType &matrixY, SparseMatrixType &matrixZ ) { 
      int gridDofs = matrixX.getNumRows();
      RealType norm = matrixX.getFrobeniusNormSqr();
      norm += matrixY.getFrobeniusNormSqr() + matrixZ.getFrobeniusNormSqr();
      _sigma = 1.0/sqrt(norm);
      _tau = _sigma; 
    }
    
    void projectionOntoIntegralConstraintBT ( const VectorType & powerLookuptable, RPOpType &rp, VectorType & phi1, VectorType & phi2 ) { 
        
        RealType alpha = 1.0-_epsilon;
        int maxLevelGrid = rp.getGridRef().getMaxLevelz();
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) );
        
        int constraintListSize = rp.getConstraintList().size();
    
        // Go through list of line segments
        #pragma omp parallel for 
        for ( int i = 0; i < constraintListSize; ++i ) { 
            
            if ( rp.getConstraintList()[i].size() != 0 ) {  
                
                int lineLength = rp.getLineSegmentList()[i].size();
                int numberOfConstraints = rp.getConstraintList()[i].size();
                VectorType phi1Vec ( lineLength );
                VectorType phi2Vec ( lineLength );
                VectorType zCoordsVec ( lineLength );
                std::vector < int > interpolPoints;
                aol::Vec2 < RealType > refCoords;
                int iter = 0;
                int k, t1, t2, totalIntervalLength;
                RealType error, hz2, A, B, beta, intervalLength, phi1Old, phi2Old, power, normPhi, value;
                 
                // Get phi1Vec and phi2Vec
                k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1Vec[k] = phi1[lineit->second]; 
                    phi2Vec[k] = phi2[lineit->second];
                    zCoordsVec[k] = lineit->first;
                }
                
                // Set lastDiffPhi1 and lastDiffPhi2
                VectorType lastDiffPhi1 ( numberOfConstraints + 1 );
                VectorType lastDiffPhi2 ( numberOfConstraints + 1 );
                lastDiffPhi1.setZero();
                lastDiffPhi2.setZero();
                
                error = 1.0;
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Dykstra projection algorithm
                while ( iter < _maxIterProject && error > _tolProject ) {   
                
                    // Iterate over single constraints
                    for ( int j = 0; j < numberOfConstraints; ++j ) {
                        
                        t1 = rp.getConstraintList()[i][j].first;
                        t2 = rp.getConstraintList()[i][j].second;
                        totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1]; 
                        
                        A = - static_cast < RealType > ( totalIntervalLength ) * hz * lastDiffPhi1[j];
                        B = - static_cast < RealType > ( totalIntervalLength ) * hz * lastDiffPhi2[j];
                        
                        for ( k = t1; k < t2; ++k ) {
                            intervalLength = static_cast < RealType > ( zCoordsVec[k+1] - zCoordsVec[k] ) * hz;
                            A += intervalLength * phi1Vec[k];
                            B += intervalLength * phi2Vec[k];
                        }
                        beta = aol::Min ( 0.0, powerLookuptable[totalIntervalLength] / ( sqrt ( A*A+B*B ) ) - 1.0 / ( static_cast < RealType > ( totalIntervalLength ) * hz ) ); 
                        
                        for ( k = t1; k < t2; ++k ) {
                            phi1Vec[k] += beta*A - lastDiffPhi1[j];
                            phi2Vec[k] += beta*B - lastDiffPhi2[j];
                        }
                        
                        lastDiffPhi1[j] = beta*A;
                        lastDiffPhi2[j] = beta*B;
                        
                    }
      
                    // Additional constraint 
                    t2 = lineLength - 1;
                    
                    phi1Old = phi1Vec[t2] - lastDiffPhi1[numberOfConstraints];
                    phi2Old = phi2Vec[t2] - lastDiffPhi2[numberOfConstraints];
                
                    hz2 = 1.0 / ( static_cast < RealType > ( 1 << maxLevelGrid ) );
                    power = pow ( hz2, alpha );
                    normPhi = sqrt ( phi1Old*phi1Old + phi2Old*phi2Old );
                    value = aol::Min ( 1.0, power / normPhi ); 
                
                    phi1Vec[t2] = phi1Old * value;
                    phi2Vec[t2] = phi2Old * value;
                    
                    lastDiffPhi1[numberOfConstraints] = phi1Vec[t2] - phi1Old;
                    lastDiffPhi2[numberOfConstraints] = phi2Vec[t2] - phi2Old;
                    
                    ++iter;
      
                    // Compute error 
                    error = 0.0;
                    for ( int j = 0; j < numberOfConstraints; ++j ) {
                        
                        int t1 = rp.getConstraintList()[i][j].first;
                        int t2 = rp.getConstraintList()[i][j].second;
                        int totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1];   
                        
                        A = 0.0;
                        B = 0.0;
                        for ( int kk = t1; kk < t2; ++kk ) {
                            intervalLength = static_cast < RealType > ( zCoordsVec[kk+1] - zCoordsVec[kk] );
                            A += phi1Vec[kk] * intervalLength; 
                            B += phi2Vec[kk] * intervalLength; 
                        }
                        A *= hz;
                        B *= hz;
                        
                        if ( sqrt ( A*A + B*B ) > powerLookuptable[totalIntervalLength]*static_cast<RealType>(totalIntervalLength)*hz + 1e-10 ) {
                            error = 1.0;
                            break;   
                        }
                            
                    }
                
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
                // Write new values to phi1 and phi2Out
                k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1[lineit->second] = phi1Vec[k];
                    phi2[lineit->second] = phi2Vec[k];
                }
                
            }
            
        }
        
    }
    
    void projectionOntoIntegralConstraintUP ( RPOpType &rp, VectorType & phi1, VectorType & phi2 ) { 

        int maxLevelGrid = rp.getGridRef().getMaxLevelz();
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) ); 
        
        int constraintListSize = rp.getConstraintList().size();
    
        // Go through list of line segments
        #pragma omp parallel for 
        for ( int i = 0; i < constraintListSize; ++i ) { 
            
            if ( rp.getConstraintList()[i].size() != 0 ) {  
                
                int lineLength = rp.getLineSegmentList()[i].size();
                int numberOfConstraints = rp.getConstraintList()[i].size();
                VectorType phi1Vec ( lineLength );
                VectorType phi2Vec ( lineLength );
                VectorType phi1VecTilde ( lineLength );
                VectorType phi2VecTilde ( lineLength );
                VectorType zCoordsVec ( lineLength );
                std::vector < int > interpolPoints;
                aol::Vec2 < RealType > refCoords;
                int iter = 0;
                int k, t1, t2, totalIntervalLength;
                RealType error, hz2, A, B, Aold, Bold, beta, intervalLength, phi1Old, phi2Old, power, normPhi, value, minimum, factor;
                 
                // Get phi1Vec and phi2Vec
                k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1Vec[k] = phi1[lineit->second]; 
                    phi2Vec[k] = phi2[lineit->second];
                    zCoordsVec[k] = lineit->first;
                }
                
                // Set lastDiffPhi1 and lastDiffPhi2
                VectorType lastDiffPhi1 ( numberOfConstraints + lineLength + 1 );
                VectorType lastDiffPhi2 ( numberOfConstraints + lineLength + 1 );
                lastDiffPhi1.setZero();
                lastDiffPhi2.setZero();
                
                error = 1.0;
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Dykstra projection algorithm
                while ( iter < _maxIterProject && error > _tolProject ) {   
                
                    // Iterate over single constraints
                    for ( int j = 0; j < numberOfConstraints; ++j ) {
                    
                        t1 = rp.getConstraintList()[i][j].first;
                        t2 = rp.getConstraintList()[i][j].second;
                        totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1]; 
                        
                        A = - static_cast < RealType > ( totalIntervalLength ) * hz * lastDiffPhi1[j];
                        B = - static_cast < RealType > ( totalIntervalLength ) * hz * lastDiffPhi2[j];
                        
                        for ( k = t1; k < t2; ++k ) {
                            intervalLength = static_cast < RealType > ( zCoordsVec[k+1] - zCoordsVec[k] ) * hz;
                            A += intervalLength * phi1Vec[k];
                            B += intervalLength * phi2Vec[k];
                        }
                        
                        
                        minimum = aol::Min ( hz * static_cast < RealType > ( totalIntervalLength ) + _epsilon, hz * static_cast < RealType > ( totalIntervalLength ) * _a ); 
                        beta = aol::Min ( 0.0, minimum / ( sqrt ( A*A+B*B ) * hz * static_cast < RealType > ( totalIntervalLength ) ) - 1.0 / ( hz * static_cast < RealType > ( totalIntervalLength ) ) );
                        
                        for ( k = t1; k < t2; ++k ) { 
                            phi1Vec[k] += beta*A - lastDiffPhi1[j];
                            phi2Vec[k] += beta*B - lastDiffPhi2[j];
                        }
                        
                        lastDiffPhi1[j] = beta*A;
                        lastDiffPhi2[j] = beta*B;
                        
                    }
                    
                    // Additional constraint: |phi_x|<=a 
                    for ( int k = 0; k < lineLength; ++k ) { 
                        
                        phi1VecTilde[k] = phi1Vec[k] - lastDiffPhi1[numberOfConstraints+k];
                        phi2VecTilde[k] = phi2Vec[k] - lastDiffPhi2[numberOfConstraints+k];
                        
                        factor = aol::Min ( 1.0, _a / sqrt ( phi1VecTilde[k]*phi1VecTilde[k] + phi2VecTilde[k]*phi2VecTilde[k] + 1e-10 ) );
    
                        phi1Vec[k] = factor * phi1VecTilde[k];
                        phi2Vec[k] = factor * phi2VecTilde[k];
                            
                        lastDiffPhi1[numberOfConstraints+k] = phi1Vec[k] - phi1VecTilde[k];
                        lastDiffPhi2[numberOfConstraints+k] = phi2Vec[k] - phi2VecTilde[k];
                            
                    }
      
                    // Additional constraint 
                    t2 = lineLength - 1;
                    phi1Old = phi1Vec[t2] - lastDiffPhi1[numberOfConstraints];
                    phi2Old = phi2Vec[t2] - lastDiffPhi2[numberOfConstraints];
                    hz2 = 1.0 / ( static_cast < RealType > ( 1 << maxLevelGrid ) ); 
                    power = aol::Min ( hz2 * _a, hz2 + _epsilon );
                    normPhi = sqrt ( phi1Old*phi1Old + phi2Old*phi2Old );
                    value = aol::Min ( 1.0, power / normPhi ); 
                    phi1Vec[t2] = phi1Old * value;
                    phi2Vec[t2] = phi2Old * value;
                    lastDiffPhi1[numberOfConstraints] = phi1Vec[t2] - phi1Old;
                    lastDiffPhi2[numberOfConstraints] = phi2Vec[t2] - phi2Old;
                    
                    ++iter;
      
                    // Compute error 
                    error = 0.0;
                    for ( int j = 0; j < numberOfConstraints; ++j ) {
                        
                        int t1 = rp.getConstraintList()[i][j].first;
                        int t2 = rp.getConstraintList()[i][j].second;
                        int totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1];   
                        
                        A = 0.0;
                        B = 0.0;
                        for ( int kk = t1; kk < t2; ++kk ) {
                            intervalLength = static_cast < RealType > ( zCoordsVec[kk+1] - zCoordsVec[kk] );
                            A += phi1Vec[kk] * intervalLength; 
                            B += phi2Vec[kk] * intervalLength; 
                        }
                        A *= hz;
                        B *= hz;
                        
                        if ( sqrt ( A*A + B*B ) > aol::Min ( _a * static_cast<RealType>(totalIntervalLength)*hz, static_cast<RealType>(totalIntervalLength)*hz + _epsilon ) + 1e-10 ) {
                            error = 1.0;
                            break;   
                        }
                            
                    }
                
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
                // Write new values to phi1 and phi2Out
                k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1[lineit->second] = phi1Vec[k];
                    phi2[lineit->second] = phi2Vec[k];                
                }
                
            }
            
        }
        
    }
    
    RealType minimize ( ExampleType & example, const int maxIter, VectorType & powerLookuptable, GridType & grid, RPOpType &rp, const SparseMatrixType &matrixX, const SparseMatrixType &matrixY, 
                        const SparseMatrixType &matrixZ, const SparseMatrixType &matrixXT, const SparseMatrixType &matrixYT, const SparseMatrixType &matrixZT, VectorType &v, VectorType &phi1, VectorType &phi2, VectorType &phi3, const int run ) {
        
        VectorType vBar ( v ); VectorType vOld ( v.size() ); VectorType vDiff ( v.size() );
        VectorType K1v ( v.size() ); VectorType K2v ( v.size() ); VectorType K3v ( v.size() );
        VectorType KTphi1 ( v.size() ); VectorType KTphi2 ( v.size() ); VectorType KTphi3 ( v.size() );
        int gridDofs = rp.getGridRef().getNumberOfDofs();
        int iter = 0;
        RealType error = 1.0;
        int stopVal = 1000;
        
        // Define finite element matrices in CSR format 
        MatrixType matrixXCSR; matrixXCSR = matrixX;
        MatrixType matrixYCSR; matrixYCSR = matrixY;
        MatrixType matrixZCSR; matrixZCSR = matrixZ; 
        MatrixType matrixXTCSR; matrixXTCSR = matrixXT; 
        MatrixType matrixYTCSR; matrixYTCSR = matrixYT; 
        MatrixType matrixZTCSR; matrixZTCSR = matrixZT; 
        
        // Include step size directly in CSR matrices 
        if ( _precon == 0 ) {
            _sigma *= 1.0/_tauFactor*0.9;
            _tau *= _tauFactor*0.9;
            matrixXCSR *= _sigma; matrixYCSR *= _sigma; matrixZCSR *= _sigma; 
            matrixXTCSR *= _tau; matrixYTCSR *= _tau; matrixZTCSR *= _tau;
        }
        else {
            this->setPreconditioner ( matrixXCSR, matrixYCSR, matrixZCSR, matrixXTCSR, matrixYTCSR, matrixZTCSR ); 
            matrixXCSR *= 1.0/_tauFactor; matrixYCSR *= 1.0/_tauFactor; matrixZCSR *= 1.0/_tauFactor; 
            matrixXTCSR *= _tauFactor; matrixYTCSR *= _tauFactor; matrixZTCSR *= _tauFactor; 
        }
        

        // ITERATION 
        while ( iter <= maxIter && aol::Abs(error) > _tol ) {
            
            if ( ( iter % stopVal ) == 0 ) {
                if ( _stoppingCrit == "PrimalDualGap" ) {
                    if ( iter == maxIter ) {
                        matrixX.apply ( v, K1v ); matrixY.apply ( v, K2v ); matrixZ.apply ( v, K3v ); 
                        matrixXT.apply ( phi1, KTphi1 ); matrixYT.apply ( phi2, KTphi2 ); matrixZT.apply ( phi3, KTphi3 );   
                        if ( example._problemType == "BranchedTransport" ) error = computePrimalDualGapBT ( example, powerLookuptable, grid, rp, v, phi1, phi2, phi3, K1v, K2v, K3v, KTphi1, KTphi2, KTphi3 );
                        if ( example._problemType == "UrbanPlanning" ) error = computePrimalDualGapUP ( example, grid, rp, v, phi1, phi2, phi3, K1v, K2v, K3v, KTphi1, KTphi2, KTphi3 ); 
                        cout << "Iteration " << iter << ", primal dual gap = " << error << endl; 
                    }
                    else {
                        cout << "Iteration " << iter << endl;    
                    }
                }
                else if ( _stoppingCrit == "PrimalDifference" ) {
                    if ( iter == 0 ) error = 1.0;
                    else {
                        vDiff = v;
                        vDiff -= vOld;
                        error = vDiff.norm();
                    }
                    cout << "Iteration " << iter << ", primal difference = " << error << endl; 
                }
                else {
                    VectorType vDiff ( v );
                    cout << "Stopping criterion not implemented" << endl;
                }
                if ( aol::Abs(error) < _tol ) {
                    cout << "Iteration stopped in " << iter << ", error = " << error << endl;
                    break;
                }
                if ( aol::Abs(error) > 1e10 ) {
                    cout << "EMERGENCY BREAK! Iteration stopped in " << iter << ", error = " << error << endl;
                    phi1.setZero();
                    phi2.setZero();
                    phi3.setZero();
                    break;
                }
            }
            
            // Save old iteration step
            vOld = v;
            
            // Update constraint list 
            if ( ( iter % 100 ) == 0 && iter <= static_cast < int > ( maxIter * 0.5 ) && run != 0 ) {
                if ( example._problemType == "BranchedTransport" ) rp.updateConstraintListBT ( powerLookuptable, phi1, phi2, 1e-5, false );  
                if ( example._problemType == "UrbanPlanning" ) rp.updateConstraintListUP ( _epsilon, _a, phi1, phi2, 1e-5, false ); 
            } 
            
            // Dual prox 
            matrixXCSR.apply ( vBar, K1v ); matrixYCSR.apply ( vBar, K2v ); matrixZCSR.apply ( vBar, K3v );
            phi1 += K1v; phi2 += K2v; phi3 += K3v; 
            if ( example._problemType == "BranchedTransport" ) this->projectionOntoIntegralConstraintBT ( powerLookuptable, rp, phi1, phi2 );
            if ( example._problemType == "UrbanPlanning" ) this->projectionOntoIntegralConstraintUP ( rp, phi1, phi2 );
            phi3.thresholdFromBelow ( 0.0, 0.0 );
            
            // Primal prox 
            matrixXTCSR.apply ( phi1, KTphi1 );
            matrixYTCSR.apply ( phi2, KTphi2 );   
            matrixZTCSR.apply ( phi3, KTphi3 );
            for ( int i = 0; i < v.size(); ++i ) {
                if ( !example._boundaryMask[i] ) {
                    v[i] = aol::Min ( aol::Max ( 0.0, v[i] - ( KTphi1[i] + KTphi2[i] + KTphi3[i] ) ), 1.0 ); 
                }
                else {
                    v[i] = example._boundaryVals[i];   
                }
            }
            
            // Overrelaxation 
            this->overrelaxation ( &vBar[0], &v[0], &vOld[0], gridDofs ); 
            
            ++iter;
            
        }
        
        return aol::Abs(error);
        
    }
    
    
    void setPreconditioner ( MatrixType &mX, MatrixType &mY, MatrixType &mZ, MatrixType &mXT, MatrixType &mYT, MatrixType &mZT ) { 
        // Sigma1 
        for ( int row = 0; row < mX.getNumRows(); ++row ) {
          RealType sigma1Val = mX.getRowSumAbs ( row );
          if ( sigma1Val > 1e-10 ) 
            mX.multiplyRowByScalar ( row, 1.0/sigma1Val ); 
        }
        // Sigma2
        for ( int row = 0; row < mY.getNumRows(); ++row ) {
          RealType sigma2Val = mY.getRowSumAbs ( row );
          if ( sigma2Val > 1e-10 ) 
            mY.multiplyRowByScalar ( row, 1.0/sigma2Val ); 
        }
        // Sigma3
        for ( int row = 0; row < mZ.getNumRows(); ++row ) {
          RealType sigma3Val = mZ.getRowSumAbs ( row );
          if ( sigma3Val > 1e-10 )
            mZ.multiplyRowByScalar ( row, 1.0/sigma3Val ); 
        }
        // Tau 
        for ( int row = 0; row < mXT.getNumRows(); ++row ) {
          RealType tauVal = mXT.getRowSumAbs ( row );
          tauVal += mYT.getRowSumAbs ( row ) + mZT.getRowSumAbs ( row );
          if ( tauVal > 1e-10 ) {
            mXT.multiplyRowByScalar ( row, 1.0/tauVal );
            mYT.multiplyRowByScalar ( row, 1.0/tauVal );
            mZT.multiplyRowByScalar ( row, 1.0/tauVal );
          }
        }
    }
  
    void renewBoundaries ( VectorType & v, ExampleType & example ) { 
        for ( int i = 0; i < v.size(); ++i ) {
            if ( example._boundaryMask[i] ) {
                v[i] = example._boundaryVals[i];   
            }
        }            
    }
    
    template < typename PointerType >
    void overrelaxation ( PointerType *vBar, const PointerType *v, const PointerType *vOld, const int &size ) {
        for ( int i = 0; i < size; ++i, ++vBar, ++v, ++vOld )
            *vBar = ( *v - *vOld ) + *v;
    } 
    
    RealType computePrimalDualResidual ( const MatrixType & matrixX, const MatrixType & matrixY, const MatrixType & matrixZ, const MatrixType & matrixXT, const MatrixType & matrixYT, const MatrixType & matrixZT, const VectorType & v, 
                                         const VectorType & phi1, const VectorType & phi2, const VectorType & phi3, const VectorType & vOld, const VectorType & phi1Old, const VectorType & phi2Old, const VectorType & phi3Old ) {
        
        int gridDofs = v.size();
        VectorType vDiff ( gridDofs );
        VectorType phi1Diff ( gridDofs ); VectorType phi2Diff ( gridDofs ); VectorType phi3Diff ( gridDofs );
        VectorType K1vdiff ( gridDofs ); VectorType K2vdiff ( gridDofs ); VectorType K3vdiff ( gridDofs );
        VectorType K1Tphi1diff ( gridDofs ); VectorType K2Tphi2diff ( gridDofs ); VectorType K3Tphi3diff ( gridDofs );
        VectorType primalRes ( gridDofs ); VectorType dualRes1 ( gridDofs ); VectorType dualRes2 ( gridDofs ); VectorType dualRes3 ( gridDofs );
    
        // Difference between iterates
        for ( int i = 0; i < gridDofs; ++i ) {
            vDiff[i]= v[i] - vOld[i];   
            phi1Diff[i] = phi1[i] - phi1Old[i]; phi2Diff[i] = phi2[i] - phi2Old[i]; phi3Diff[i] = phi3[i] - phi3Old[i];
        }
  
        // Primal residual
        matrixXT.apply ( phi1Diff, K1Tphi1diff ); matrixYT.apply ( phi2Diff, K2Tphi2diff ); matrixZT.apply ( phi3Diff, K3Tphi3diff );
        for ( int i = 0; i < gridDofs; ++i ) {
            primalRes[i] = vDiff[i] / _tau - ( K1Tphi1diff[i] + K2Tphi2diff[i] + K3Tphi3diff[i] );   
        }
  
        // Dual residual
        matrixX.apply ( vDiff, K1vdiff ); matrixY.apply ( vDiff, K2vdiff ); matrixZ.apply ( vDiff, K3vdiff );
        for ( int i = 0; i < gridDofs; ++i ) {
            dualRes1[i] = phi1Diff[i] / _sigma - K1vdiff[i];
            dualRes2[i] = phi2Diff[i] / _sigma - K2vdiff[i];
            dualRes3[i] = phi3Diff[i] / _sigma - K3vdiff[i];
        }
  
        // Norm of residuals
        RealType primalResVal = 0.0;
        RealType dualResVal = 0.0;
        for ( int i = 0; i < gridDofs; ++i ) {
            primalResVal += abs(primalRes[i]);
            dualResVal += abs(dualRes1[i]) + abs(dualRes2[i]) + abs(dualRes3[i]);
        }
        
        return ( primalResVal + dualResVal );
        
    }
    
    #ifdef USE_MOSEK
    void computePhiOpt ( RPOpType & rp, const VectorType & v, const VectorType & K1v, const VectorType & K2v, const VectorType & K3v, VectorType & phi1Opt, VectorType & phi2Opt, VectorType & phi3Opt ) {
     
      // Parameters 
      int gridDofs = rp.getGridRef().getNumberOfDofs();
      rp.setConstraintListFull(); 
      int constraintListSize = rp.getConstraintListFull().size();  
      RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) ); 
      
      // Go through list of line segments and optimize independently 
      #pragma omp parallel for 
      for ( int i = 0; i < constraintListSize; ++i ) { 
          
        int numberOfConstraints = rp.getConstraintListFull()[i].size();

        // Extract corresponding values of Kv 
        int lineLength = rp.getLineSegmentList()[i].size();
        int k = 0, t1, t2, totalIntervalLength;
        RealType RHS; 
        VectorType zCoordsLoc ( lineLength );
        std::vector < RealType > K1vLoc ( lineLength ); std::vector < RealType > K2vLoc ( lineLength ); std::vector < RealType > K3vLoc ( lineLength ); 
        for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
        K1vLoc[k] = K1v[lineit->second];
        K2vLoc[k] = K2v[lineit->second];
        K3vLoc[k] = K3v[lineit->second];
        zCoordsLoc[k] = lineit->first;
        }
        std::shared_ptr < ndarray < RealType, 1 > > K1vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K2vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K3vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::copy ( K1vLoc.begin(), K1vLoc.end(), K1vLocPtr->begin() );
        std::copy ( K2vLoc.begin(), K2vLoc.end(), K2vLocPtr->begin() );
        std::copy ( K3vLoc.begin(), K3vLoc.end(), K3vLocPtr->begin() );
        
        // Model 
        Model::t M = new Model("primalEnergyLoc"); 
    
        // Variables 
        mosek::fusion::Variable::t phi1OptLoc = M->variable ( "phi1OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi2OptLoc = M->variable ( "phi2OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi3OptLoc = M->variable ( "phi3OptLoc", lineLength, Domain::inRange(0.0,1e5f) ); // ATTENTION: Needs to be bounded from above 
        
        // Constraints and additional variables 
        for ( int j = 0; j < numberOfConstraints; ++j ) {
            t1 = rp.getConstraintListFull()[i][j].first;
            t2 = rp.getConstraintListFull()[i][j].second;
            totalIntervalLength = zCoordsLoc[t2] - zCoordsLoc[t1];  
            RHS = pow ( hz * static_cast<RealType>(totalIntervalLength), 1.0-_epsilon ); 
            // Compute length vector 
            std::vector < RealType > lengthVec ( lineLength ); 
            for ( int kk = 0; kk < lineLength; ++kk ) {
                if ( kk >= t1 && kk < t2 ) 
                    lengthVec[kk] = static_cast < RealType > ( zCoordsLoc[kk+1] - zCoordsLoc[kk] ) * hz;
                else 
                    lengthVec[kk] = 0.0;
            }
            std::shared_ptr < ndarray < RealType, 1 > > lengthVecPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
            std::copy ( lengthVec.begin(), lengthVec.end(), lengthVecPtr->begin() );
            // Add cone constraint 
            M->constraint ( Expr::vstack ( RHS, Expr::dot(lengthVecPtr,phi1OptLoc), Expr::dot(lengthVecPtr,phi2OptLoc) ), Domain::inQCone(3) ); 
        }    
        
        // Objective 
        M->objective ("obj", ObjectiveSense::Maximize, Expr::add ( Expr::add( Expr::dot(K1vLocPtr,phi1OptLoc), Expr::dot(K2vLocPtr,phi2OptLoc) ), Expr::dot(K3vLocPtr,phi3OptLoc) ) );
        
        // Solve the problem and get solution  
        M->solve();   
        M->acceptedSolutionStatus ( AccSolutionStatus::Optimal ); 
        ndarray < RealType, 1 > phi1OptLocSol = *(phi1OptLoc->level()); 
        ndarray < RealType, 1 > phi2OptLocSol = *(phi2OptLoc->level());
        ndarray < RealType, 1 > phi3OptLocSol = *(phi3OptLoc->level());

        // Write result to phiOpt (global)
        k = 0;
        for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            phi1Opt[lineit->second] = phi1OptLocSol[k];
            phi2Opt[lineit->second] = phi2OptLocSol[k];
            phi3Opt[lineit->second] = phi3OptLocSol[k];
        } 
        
        M->dispose(); 
            
      } 
        
    }
    
    RealType computePrimalDualGapBT ( ExampleType & example, VectorType & powerLookuptable, GridType & grid, RPOpType & rp, const VectorType & v, const VectorType & phi1, const VectorType & phi2, const VectorType & phi3, 
                                      const VectorType & K1v, const VectorType & K2v, const VectorType & K3v, const VectorType & K1Tphi1, const VectorType & K2Tphi2, const VectorType & K3Tphi3 ) {
        
      cout << "USE MOSEK VERSION" << endl;  
        
      int gridDofs = rp.getGridRef().getNumberOfDofs();
      rp.setConstraintListFull(); 
      int constraintListSize = rp.getConstraintListFull().size();
      VectorType phi1Opt ( gridDofs ); VectorType phi2Opt ( gridDofs ); VectorType phi3Opt ( gridDofs );   
      RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) );
      
      //------------------------------------------------------------------------------------------------------
      // Compute phiOpt using MOSEK 

      // Go through list of line segments and optimize independently 
      #pragma omp parallel for 
      for ( int i = 0; i < constraintListSize; ++i ) { 
          
        int numberOfConstraints = rp.getConstraintListFull()[i].size();

        // Extract corresponding values of Kv 
        int lineLength = rp.getLineSegmentList()[i].size();
        int k = 0, t1, t2, totalIntervalLength;
        RealType RHS; 
        VectorType zCoordsLoc ( lineLength );
        std::vector < RealType > K1vLoc ( lineLength ); std::vector < RealType > K2vLoc ( lineLength ); std::vector < RealType > K3vLoc ( lineLength ); 
        for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
        K1vLoc[k] = K1v[lineit->second];
        K2vLoc[k] = K2v[lineit->second];
        K3vLoc[k] = K3v[lineit->second];
        zCoordsLoc[k] = lineit->first;
        }
        std::shared_ptr < ndarray < RealType, 1 > > K1vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K2vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K3vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::copy ( K1vLoc.begin(), K1vLoc.end(), K1vLocPtr->begin() );
        std::copy ( K2vLoc.begin(), K2vLoc.end(), K2vLocPtr->begin() );
        std::copy ( K3vLoc.begin(), K3vLoc.end(), K3vLocPtr->begin() );
        
        // Model 
        Model::t M = new Model("primalEnergyLoc"); 
    
        // Variables 
        mosek::fusion::Variable::t phi1OptLoc = M->variable ( "phi1OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi2OptLoc = M->variable ( "phi2OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi3OptLoc = M->variable ( "phi3OptLoc", lineLength, Domain::inRange(0.0,1e5f) ); // ATTENTION: Needs to be bounded from above 
        
        // Constraints and additional variables 
        for ( int j = 0; j < numberOfConstraints; ++j ) {
            t1 = rp.getConstraintListFull()[i][j].first;
            t2 = rp.getConstraintListFull()[i][j].second;
            totalIntervalLength = zCoordsLoc[t2] - zCoordsLoc[t1];  
            RHS = pow ( hz * static_cast<RealType>(totalIntervalLength), 1.0-_epsilon ); 
            // Compute length vector 
            std::vector < RealType > lengthVec ( lineLength ); 
            for ( int kk = 0; kk < lineLength; ++kk ) {
                if ( kk >= t1 && kk < t2 ) 
                    lengthVec[kk] = static_cast < RealType > ( zCoordsLoc[kk+1] - zCoordsLoc[kk] ) * hz;
                else 
                    lengthVec[kk] = 0.0;
            }
            std::shared_ptr < ndarray < RealType, 1 > > lengthVecPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
            std::copy ( lengthVec.begin(), lengthVec.end(), lengthVecPtr->begin() );
            // Add cone constraint 
            M->constraint ( Expr::vstack ( RHS, Expr::dot(lengthVecPtr,phi1OptLoc), Expr::dot(lengthVecPtr,phi2OptLoc) ), Domain::inQCone(3) ); 
        }    
        
        // Objective 
        M->objective ("obj", ObjectiveSense::Maximize, Expr::add ( Expr::add( Expr::dot(K1vLocPtr,phi1OptLoc), Expr::dot(K2vLocPtr,phi2OptLoc) ), Expr::dot(K3vLocPtr,phi3OptLoc) ) );
        
        // Solve the problem and get solution  
        try {
          M->solve();   
          M->acceptedSolutionStatus ( AccSolutionStatus::Optimal ); 
          ndarray < RealType, 1 > phi1OptLocSol = *(phi1OptLoc->level()); 
          ndarray < RealType, 1 > phi2OptLocSol = *(phi2OptLoc->level());
          ndarray < RealType, 1 > phi3OptLocSol = *(phi3OptLoc->level());

          // Write result to phiOpt (global)
          k = 0;
          for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            phi1Opt[lineit->second] = phi1OptLocSol[k];
            phi2Opt[lineit->second] = phi2OptLocSol[k];
            phi3Opt[lineit->second] = phi3OptLocSol[k];
          } 
            
          M->dispose(); 
        }
        catch ( const SolutionError &e ) {
           
          // Write result to phiOpt (global)
          k = 0;
          for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            phi1Opt[lineit->second] = phi1[lineit->second];
            phi2Opt[lineit->second] = phi2[lineit->second];
            phi3Opt[lineit->second] = phi3[lineit->second];
          }   
            
        }
            
      }
      //------------------------------------------------------------------------------------------------------
      
      // Get vOpt and compute PD gap 
      VectorType vOpt ( v );
      RealType primal = 0.0;
      RealType dual = 0.0;
      RealType saddlePoint = 0.0;
      VectorType divPhiVec ( gridDofs );
      for ( int i = 0; i < gridDofs; ++i ) {
        divPhiVec[i] = K1Tphi1[i] + K2Tphi2[i] + K3Tphi3[i];
        if ( !example._boundaryMask[i] ) {
          if ( divPhiVec[i] >= 0.0 ) { 
            vOpt[i] = 0.0;
          }
          else if ( divPhiVec[i] < 0.0 ) {
            vOpt[i] = 1.0;
          }
        }
        else {
          vOpt[i] = example._boundaryVals[i];   
        }
        // Compute primal dual gap 
        primal += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
        dual += vOpt[i]*divPhiVec[i];
        saddlePoint += phi1[i]*K1v[i] + phi2[i]*K2v[i] + phi3[i]*K3v[i];            
      }   

      cout << "Primal = " << primal << ", dual = " << dual << ", saddle point = " << saddlePoint << endl;
      
      return ( primal - dual );

    }
    
    RealType computePrimalDualGapUP ( ExampleType & example, GridType & grid, RPOpType & rp, const VectorType & v, const VectorType & phi1, const VectorType & phi2, const VectorType & phi3, 
                                      const VectorType & K1v, const VectorType & K2v, const VectorType & K3v, const VectorType & K1Tphi1, const VectorType & K2Tphi2, const VectorType & K3Tphi3 ) {
        
      int gridDofs = rp.getGridRef().getNumberOfDofs();
      rp.setConstraintListFull(); 
      int constraintListSize = rp.getConstraintListFull().size();
      VectorType phi1Opt ( gridDofs ); VectorType phi2Opt ( gridDofs ); VectorType phi3Opt ( gridDofs );   
      RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) );
      
      //------------------------------------------------------------------------------------------------------
      // Compute phiOpt using MOSEK 

      // Go through list of line segments and optimize independently 
      #pragma omp parallel for 
      for ( int i = 0; i < constraintListSize; ++i ) { 
          
        int numberOfConstraints = rp.getConstraintListFull()[i].size();

        // Extract corresponding values of Kv 
        int lineLength = rp.getLineSegmentList()[i].size();
        int k = 0, t1, t2, totalIntervalLength;
        RealType RHS; 
        VectorType zCoordsLoc ( lineLength );
        std::vector < RealType > K1vLoc ( lineLength ); std::vector < RealType > K2vLoc ( lineLength ); std::vector < RealType > K3vLoc ( lineLength ); 
        for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
        K1vLoc[k] = K1v[lineit->second];
        K2vLoc[k] = K2v[lineit->second];
        K3vLoc[k] = K3v[lineit->second];
        zCoordsLoc[k] = lineit->first;
        }
        std::shared_ptr < ndarray < RealType, 1 > > K1vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K2vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::shared_ptr < ndarray < RealType, 1 > > K3vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
        std::copy ( K1vLoc.begin(), K1vLoc.end(), K1vLocPtr->begin() );
        std::copy ( K2vLoc.begin(), K2vLoc.end(), K2vLocPtr->begin() );
        std::copy ( K3vLoc.begin(), K3vLoc.end(), K3vLocPtr->begin() );
        
        // Model 
        Model::t M = new Model("primalEnergyLoc"); 
    
        // Variables 
        mosek::fusion::Variable::t phi1OptLoc = M->variable ( "phi1OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi2OptLoc = M->variable ( "phi2OptLoc", lineLength, Domain::unbounded() );
        mosek::fusion::Variable::t phi3OptLoc = M->variable ( "phi3OptLoc", lineLength, Domain::inRange(0.0,1e8f) ); // ATTENTION: Needs to be bounded from above 
        
        // Constraints and additional variables 
        for ( int j = 0; j < numberOfConstraints; ++j ) {
            t1 = rp.getConstraintListFull()[i][j].first;
            t2 = rp.getConstraintListFull()[i][j].second;
            totalIntervalLength = zCoordsLoc[t2] - zCoordsLoc[t1];  
            RHS = aol::Min ( hz*static_cast<RealType>(totalIntervalLength)+_epsilon, hz*static_cast<RealType>(totalIntervalLength)*_a ); 
            // Compute length vector 
            std::vector < RealType > lengthVec ( lineLength ); 
            for ( int kk = 0; kk < lineLength; ++kk ) {
                if ( kk >= t1 && kk < t2 ) 
                    lengthVec[kk] = static_cast < RealType > ( zCoordsLoc[kk+1] - zCoordsLoc[kk] ) * hz;
                else 
                    lengthVec[kk] = 0.0;
            }
            std::shared_ptr < ndarray < RealType, 1 > > lengthVecPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
            std::copy ( lengthVec.begin(), lengthVec.end(), lengthVecPtr->begin() );
            // Add cone constraint 
            M->constraint ( Expr::vstack ( RHS, Expr::dot(lengthVecPtr,phi1OptLoc), Expr::dot(lengthVecPtr,phi2OptLoc) ), Domain::inQCone(3) ); 
        }    
        
        // Constraint |phi|<=a 
        for ( int j = 0; j < lineLength; ++j ) {
            M->constraint ( Expr::vstack ( _a, phi1OptLoc->index(j), phi2OptLoc->index(j) ), Domain::inQCone(3) );
        }
        
        // Objective 
        M->objective ("obj", ObjectiveSense::Maximize, Expr::add ( Expr::add( Expr::dot(K1vLocPtr,phi1OptLoc), Expr::dot(K2vLocPtr,phi2OptLoc) ), Expr::dot(K3vLocPtr,phi3OptLoc) ) );
        
        // Solve the problem and get solution  
        M->solve();   
        M->acceptedSolutionStatus ( AccSolutionStatus::Optimal ); 
        ndarray < RealType, 1 > phi1OptLocSol = *(phi1OptLoc->level()); 
        ndarray < RealType, 1 > phi2OptLocSol = *(phi2OptLoc->level());
        ndarray < RealType, 1 > phi3OptLocSol = *(phi3OptLoc->level());

        // Write result to phiOpt (global)
        k = 0;
        for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            phi1Opt[lineit->second] = phi1OptLocSol[k];
            phi2Opt[lineit->second] = phi2OptLocSol[k];
            phi3Opt[lineit->second] = phi3OptLocSol[k];
        } 
        
        M->dispose(); 
            
      }
      //------------------------------------------------------------------------------------------------------
      
      // Get vOpt and compute PD gap 
      VectorType vOpt ( v );
      RealType primal = 0.0;
      RealType dual = 0.0;
      RealType saddlePoint = 0.0;
      VectorType divPhiVec ( gridDofs );
      for ( int i = 0; i < gridDofs; ++i ) {
        divPhiVec[i] = K1Tphi1[i] + K2Tphi2[i] + K3Tphi3[i];
        if ( !example._boundaryMask[i] ) {
          if ( divPhiVec[i] >= 0.0 ) {   
            vOpt[i] = 0.0;
          }
          else if ( divPhiVec[i] < 0.0 ) {
            vOpt[i] = 1.0;
          }
        }
        else {
          vOpt[i] = example._boundaryVals[i];   
        }
        // Compute primal dual gap 
        primal += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
        dual += vOpt[i]*divPhiVec[i];
        saddlePoint += phi1[i]*K1v[i] + phi2[i]*K2v[i] + phi3[i]*K3v[i];            
      }   

      cout << "Primal = " << primal << ", dual = " << dual << ", saddle point = " << saddlePoint << endl;
      
      return ( primal - dual );

    }
    #else 
    RealType computePrimalDualGapBT ( ExampleType & example, VectorType & powerLookuptable, GridType & /*grid*/, RPOpType & rp, const VectorType & v, const VectorType & phi1, const VectorType & phi2, const VectorType & phi3, 
                                      const VectorType & K1v, const VectorType & K2v, const VectorType & K3v, const VectorType & K1Tphi1, const VectorType & K2Tphi2, const VectorType & K3Tphi3 ) {
        
        int gridDofs = rp.getGridRef().getNumberOfDofs();
        
        // Set variables
        VectorType vOpt ( v );
        VectorType phi1Opt ( phi1 );
        VectorType phi2Opt ( phi2 );
        VectorType phi3Opt ( phi3 );
        VectorType phi1OptOld ( gridDofs );
        VectorType phi2OptOld ( gridDofs );
        VectorType phi3OptOld ( gridDofs );
        VectorType phi1OptCheck ( gridDofs );
        VectorType phi2OptCheck ( gridDofs );
        VectorType phi3OptCheck ( gridDofs );
        
        // Get phiOpt 
        int iter = 0;
        RealType error = 1.0;
        RealType dt = 1.0;
        RealType energyOld, energyNew, energy;
        bool check = false;
        while ( iter < 1000 && error > 1e-5 ) {
            phi1OptOld = phi1Opt;
            phi2OptOld = phi2Opt;
            phi3OptOld = phi3Opt;
            energyOld = 0.0;
            for ( int i = 0; i < gridDofs; ++i ) {
                energyOld += phi1Opt[i] * K1v[i] + phi2Opt[i] * K2v[i] + phi3Opt[i] * K3v[i];
            }
            check = false;
            dt = 10.0;
            while (check == false) {
                for ( int i = 0; i < gridDofs; ++i ) {
                    phi1OptCheck[i] = phi1Opt[i] + dt * K1v[i];
                    phi2OptCheck[i] = phi2Opt[i] + dt * K2v[i];
                    phi3OptCheck[i] = phi3Opt[i] + dt * K3v[i];
                }
                this->projectionOntoIntegralConstraintBT ( powerLookuptable, rp, phi1OptCheck, phi2OptCheck );
                phi3OptCheck.thresholdFromBelow ( 0.0, 0.0 );
                energyNew = 0.0;
                for ( int i = 0; i < gridDofs; ++i ) {
                    energyNew += phi1OptCheck[i] * K1v[i] + phi2OptCheck[i] * K2v[i] + phi3OptCheck[i] * K3v[i];
                }
                if ( energyNew > energyOld ) {
                    check = true;
                    phi1Opt = phi1OptCheck;
                    phi2Opt = phi2OptCheck;
                    phi3Opt = phi3OptCheck;
                }
                else {
                    dt = dt / 2.0;
                }
                if (dt < 1e-10) {
                    check = true;
                    iter = 100000;
                }
            }
            error = 0.0;
            energy = 0.0;
            for ( int i = 0; i < gridDofs; ++i ) {
                error += aol::Abs ( phi1Opt[i] - phi1OptOld[i] ) + aol::Abs ( phi2Opt[i] - phi2OptOld[i] ) + aol::Abs ( phi3Opt[i] - phi3OptOld[i] );
                energy += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
            }
            error = error / gridDofs;
            ++iter;
        }  
        
        // Get vOpt and compute PD gap 
        RealType primal = 0.0;
        RealType dual = 0.0;
        RealType saddlePoint = 0.0;
        VectorType divPhiVec ( gridDofs );
        for ( int i = 0; i < gridDofs; ++i ) {
            divPhiVec[i] = K1Tphi1[i] + K2Tphi2[i] + K3Tphi3[i];
            if ( !example._boundaryMask[i] ) {
                if ( divPhiVec[i] >= (1e-8) ) { 
                    vOpt[i] = 0.0;
                }
                else if ( divPhiVec[i] < -(1e-8) ) {
                    vOpt[i] = 1.0;
                }
            }
            else {
                vOpt[i] = example._boundaryVals[i];   
            }
            // Compute primal dual gap 
            primal += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
            dual += vOpt[i]*divPhiVec[i];
            saddlePoint += phi1[i]*K1v[i] + phi2[i]*K2v[i] + phi3[i]*K3v[i];            
        }   
        
        return ( primal - dual );
        
    }   
    
    RealType computePrimalDualGapUP ( ExampleType & example, GridType & /*grid*/, RPOpType & rp, const VectorType & v, const VectorType & phi1, const VectorType & phi2, const VectorType & phi3, 
                                      const VectorType & K1v, const VectorType & K2v, const VectorType & K3v, const VectorType & K1Tphi1, const VectorType & K2Tphi2, const VectorType & K3Tphi3 ) {
        
        int gridDofs = rp.getGridRef().getNumberOfDofs();
        
        // Set variables
        VectorType vOpt ( v );
        VectorType phi1Opt ( phi1 );
        VectorType phi2Opt ( phi2 );
        VectorType phi3Opt ( phi3 );
        VectorType phi1OptOld ( gridDofs );
        VectorType phi2OptOld ( gridDofs );
        VectorType phi3OptOld ( gridDofs );
        VectorType phi1OptCheck ( gridDofs );
        VectorType phi2OptCheck ( gridDofs );
        VectorType phi3OptCheck ( gridDofs );
            
        // Get phiOpt 
        int iter = 0;
        RealType error = 1.0;
        RealType dt = 1.0;
        RealType energyOld, energyNew, energy;
        bool check = false;
        while ( iter < 1000 && error > 1e-5 ) {
            phi1OptOld = phi1Opt;
            phi2OptOld = phi2Opt;
            phi3OptOld = phi3Opt;
            energyOld = 0.0;
            for ( int i = 0; i < gridDofs; ++i ) {
                energyOld += phi1Opt[i] * K1v[i] + phi2Opt[i] * K2v[i] + phi3Opt[i] * K3v[i];
            }
            check = false;
            dt = 10.0;
            while (check == false) {
                for ( int i = 0; i < gridDofs; ++i ) {
                    phi1OptCheck[i] = phi1Opt[i] + dt * K1v[i];
                    phi2OptCheck[i] = phi2Opt[i] + dt * K2v[i];
                    phi3OptCheck[i] = phi3Opt[i] + dt * K3v[i];
                }
                this->projectionOntoIntegralConstraintUP ( rp, phi1OptCheck, phi2OptCheck );
                phi3OptCheck.thresholdFromBelow ( 0.0, 0.0 );
                energyNew = 0.0;
                for ( int i = 0; i < gridDofs; ++i ) {
                    energyNew += phi1OptCheck[i] * K1v[i] + phi2OptCheck[i] * K2v[i] + phi3OptCheck[i] * K3v[i];
                }
                if ( energyNew > energyOld ) {
                    check = true;
                    phi1Opt = phi1OptCheck;
                    phi2Opt = phi2OptCheck;
                    phi3Opt = phi3OptCheck;
                }
                else {
                    dt = dt / 2.0;
                }
                if (dt < 1e-10) {
                    check = true;
                    iter = 100000;
                }
            }
            error = 0.0;
            energy = 0.0;
            for ( int i = 0; i < gridDofs; ++i ) {
                error += aol::Abs ( phi1Opt[i] - phi1OptOld[i] ) + aol::Abs ( phi2Opt[i] - phi2OptOld[i] ) + aol::Abs ( phi3Opt[i] - phi3OptOld[i] );
                energy += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
            }
            error = error / gridDofs;
            ++iter;
        }        
        
        // Get vOpt 
        RealType primal = 0.0;
        RealType dual = 0.0;
        VectorType divPhiVec ( gridDofs );
        for ( int i = 0; i < gridDofs; ++i ) {
            divPhiVec[i] = K1Tphi1[i] + K2Tphi2[i] + K3Tphi3[i];
            if ( !example._boundaryMask[i] ) {
                if ( divPhiVec[i] > 1e-10 ) {
                    vOpt[i] = 0.0;
                }
                else if ( divPhiVec[i] < -(1e-10) ) {
                    vOpt[i] = 1.0;
                }
            }
            else {
                vOpt[i] = example._boundaryVals[i];   
            }
            // Compute primal dual gap 
            primal += phi1Opt[i]*K1v[i] + phi2Opt[i]*K2v[i] + phi3Opt[i]*K3v[i];
            dual += vOpt[i]*divPhiVec[i];            
        }     
        
        return ( primal - dual );
        
    }   
    
    #endif
    
    int checkConstraintUniformBT ( RPOpType & rp, VectorType & powerLookuptable, VectorType & phi1, VectorType & phi2 ) {
        
        int numNotSatisfied = 0;
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) );
        RealType alpha = 1.0-_epsilon;
        
        // Update constraint list to all possible constraints
        rp.updateConstraintListBT ( powerLookuptable, phi1, phi2, 1000.0, false );
        
        // Go through all constraints
        int constraintListSize = rp.getConstraintList().size();
        for ( int i = 0; i < constraintListSize; ++i ) {
            
            if ( rp.getConstraintList()[i].size() != 0 ) { 
                
                int lineLength = rp.getLineSegmentList()[i].size();
                int numberOfConstraints = rp.getConstraintList()[i].size();
                VectorType phi1Vec ( lineLength );
                VectorType phi2Vec ( lineLength );
                VectorType zCoordsVec ( lineLength );
                std::vector < int > interpolPoints;
                aol::Vec2 < RealType > refCoords;
                int t1, t2, totalIntervalLength;
                RealType A, B, intervalLength;
                
                // Get phi1Vec and phi2Vec
                int k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1Vec[k] = phi1[lineit->second];
                    phi2Vec[k] = phi2[lineit->second];
                    zCoordsVec[k] = lineit->first;
                }
                
                // Go through all constraints
                for ( int j = 0; j < numberOfConstraints; ++j ) { 
                    
                    t1 = rp.getConstraintList()[i][j].first;
                    t2 = rp.getConstraintList()[i][j].second;
                    totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1];
                    
                    // Compute integral
                    A = 0.0;
                    B = 0.0;
                    for ( k = t1; k < t2; ++k ) {
                        intervalLength = static_cast < RealType > ( ( zCoordsVec[k+1] - zCoordsVec[k] ) ) * hz;
                        A += intervalLength * phi1Vec[k];
                        B += intervalLength * phi2Vec[k];  
                    }
                    
                    // Check constraint 
                    if ( sqrt ( A*A + B*B ) > pow ( totalIntervalLength, alpha ) + 1e-15 ) {  
                        ++numNotSatisfied;   
                    }
                    
                }
                
            }
            
        }
        
        return numNotSatisfied;
        
    }
    
    int checkConstraintUniformUP ( RPOpType & rp, VectorType & phi1, VectorType & phi2 ) { 
        
        int numNotSatisfied = 0;
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rp.getStopLevelz() ) ); 
        
        // Update constraint list to all possible constraints
        rp.updateConstraintListUP ( _epsilon, _a, phi1, phi2, 1000.0, false );
        
        rp.extendVectorZConst ( phi1 );
        rp.extendVectorZConst ( phi2 );
        
        // Go through all constraints
        int constraintListSize = rp.getConstraintList().size();
        for ( int i = 0; i < constraintListSize; ++i ) {
            
            if ( rp.getConstraintList()[i].size() != 0 ) { 
                
                int lineLength = rp.getLineSegmentList()[i].size();
                int numberOfConstraints = rp.getConstraintList()[i].size();
                VectorType phi1Vec ( lineLength );
                VectorType phi2Vec ( lineLength );
                VectorType zCoordsVec ( lineLength );
                std::vector < bool > hangingVec ( lineLength );
                std::vector < int > interpolPoints;
                aol::Vec2 < RealType > refCoords;
                int t1, t2, totalIntervalLength;
                RealType A, B, intervalLength;
                
                // Get phi1Vec and phi2Vec
                int k = 0;
                for ( std::map <int,int>::const_iterator lineit = rp.getLineSegmentList()[i].begin(); lineit != rp.getLineSegmentList()[i].end(); ++lineit, ++k ) {
                    phi1Vec[k] = phi1[lineit->second];
                    phi2Vec[k] = phi2[lineit->second];
                    zCoordsVec[k] = lineit->first;
                    hangingVec[k] = rp.getGridRef().isHangingNode ( lineit->second );
                }
                
                // Go through all constraints
                for ( int j = 0; j < numberOfConstraints; ++j ) { 
                    
                    t1 = rp.getConstraintList()[i][j].first;
                    t2 = rp.getConstraintList()[i][j].second;
                    totalIntervalLength = zCoordsVec[t2] - zCoordsVec[t1];
                    
                    // Compute integral
                    A = 0.0;
                    B = 0.0;
                    for ( k = t1; k < t2; ++k ) {
                        intervalLength = static_cast < RealType > ( ( zCoordsVec[k+1] - zCoordsVec[k] ) ) * hz;
                        A += intervalLength * phi1Vec[k];
                        B += intervalLength * phi2Vec[k];  
                    }
                    
                    // Check constraint 
                    if ( sqrt ( A*A + B*B ) > aol::Min ( _a * static_cast<RealType>(totalIntervalLength)*hz, static_cast<RealType>(totalIntervalLength)*hz + _epsilon ) + 1e-15 ) { 
                        ++numNotSatisfied;   
                    }
                    
                }
                
            }
            
        }
        
        rp.restrictVectorZConst ( phi1 );
        rp.restrictVectorZConst ( phi2 );
        
        return numNotSatisfied;
        
    }
    
    RealType checkBoundaryConditions ( RPOpType & rp, typename GridType::ElementIterator & subelement, VectorType & vBoundaryDiff, const RealType zCoordMotherEl ) {         
      // Check if element is the lower one in mother element 
      RealType valZ = 0.0;
      aol::Vec < 6, int > nodeInds = subelement.getNodeIndices();
      aol::Vec3 < RealType > coordsNode3 = rp.getGridRef().getNodeCoords ( nodeInds[3] );
      if ( coordsNode3[2] < zCoordMotherEl ) {
        aol::Vec3 < RealType > coordsNode0 = rp.getGridRef().getNodeCoords ( nodeInds[0] ); 
        aol::Vec3 < RealType > coordsNode1 = rp.getGridRef().getNodeCoords ( nodeInds[1] );
        aol::Vec3 < RealType > coordsNode2 = rp.getGridRef().getNodeCoords ( nodeInds[2] );
        RealType node0OnBound = static_cast < RealType > ( coordsNode0[1] == 0.0 || coordsNode0[1] == 1.0 || coordsNode0[0] == 0.0 || coordsNode0[0] == 1.0 );
        RealType node1OnBound = static_cast < RealType > ( coordsNode1[1] == 0.0 || coordsNode1[1] == 1.0 || coordsNode1[0] == 0.0 || coordsNode1[0] == 1.0 );
        RealType node2OnBound = static_cast < RealType > ( coordsNode2[1] == 0.0 || coordsNode2[1] == 1.0 || coordsNode2[0] == 0.0 || coordsNode2[0] == 1.0 );
        valZ = aol::Abs ( vBoundaryDiff[nodeInds[3]]-vBoundaryDiff[nodeInds[0]] ) * node0OnBound;
        valZ += aol::Abs ( vBoundaryDiff[nodeInds[4]]-vBoundaryDiff[nodeInds[1]] ) * node1OnBound;
        valZ += aol::Abs ( vBoundaryDiff[nodeInds[5]]-vBoundaryDiff[nodeInds[2]] ) * node2OnBound;
      }
      return valZ; 
    }
    
    bool elementOnBoundary ( RPOpType & rp, typename GridType::ElementIterator & element ) {
      aol::Vec < 6, int > nodeInds = element.getNodeIndices();
      aol::Vec3 < RealType > coordsNode0 = rp.getGridRef().getNodeCoords ( nodeInds[0] ); 
      aol::Vec3 < RealType > coordsNode1 = rp.getGridRef().getNodeCoords ( nodeInds[1] ); 
      aol::Vec3 < RealType > coordsNode2 = rp.getGridRef().getNodeCoords ( nodeInds[2] ); 
      bool check1 = ( coordsNode0[1] == 0.0 || coordsNode0[1] == 1.0 || coordsNode0[0] == 0.0 || coordsNode0[0] == 1.0 );
      bool check2 = ( coordsNode1[1] == 0.0 || coordsNode1[1] == 1.0 || coordsNode1[0] == 0.0 || coordsNode1[0] == 1.0 ); 
      bool check3 = ( coordsNode2[1] == 0.0 || coordsNode2[1] == 1.0 || coordsNode2[0] == 0.0 || coordsNode2[0] == 1.0 ); 
      return ( check1 || check2 || check3 );
    }
    
    bool elementOnFrontBackBoundary ( RPOpType & rp, typename GridType::ElementIterator & element ) {
      aol::Vec < 6, int > nodeInds = element.getNodeIndices();
      aol::Vec3 < RealType > coordsNode0 = rp.getGridRef().getNodeCoords ( nodeInds[0] ); 
      aol::Vec3 < RealType > coordsNode1 = rp.getGridRef().getNodeCoords ( nodeInds[1] ); 
      aol::Vec3 < RealType > coordsNode2 = rp.getGridRef().getNodeCoords ( nodeInds[2] ); 
      bool check1 = ( coordsNode0[1] == 0.0 || coordsNode0[1] == 1.0 );
      bool check2 = ( coordsNode1[1] == 0.0 || coordsNode1[1] == 1.0 ); 
      bool check3 = ( coordsNode2[1] == 0.0 || coordsNode2[1] == 1.0 ); 
      return ( check1 || check2 || check3 );
    }
    
    int computeLocalGradient ( GridType & grid, RPOpType & rp, VectorType & v, const RealType gradientTol ) {
     
        rp.extendVectorZConst ( v );
        QuadRule quadrule;
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFunc ( rp.getConfiguratorRef(), v );
        int numElements = grid.getNumberOfElements();
        VectorType gradientXVal ( numElements );
        VectorType gradientYVal ( numElements );
        VectorType gradientZVal ( numElements );
        VectorType boundaryXYVal ( numElements );
        VectorType boundaryZVal ( numElements );
        typename ConfiguratorType::VecType gradVVal;
        RealType minBaseArea = 0.5 / static_cast < RealType > ( 1 << ( 2 * ( rp.getStopLevelxy() - 1 ) ) ); 
        
        // Go through all elements and compute gradient
        int index = 0;
        int numMarkedXY = 0; int numMarkedZ = 0;
        RealType boundVal; 
        typename ConfiguratorType::DomVecType refCoords;
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          // Case 1: Element is not on boundaries  
          for ( int q = 0; q < QuadRule::numQuadPoints; ++q ) {
            refCoords = quadrule.getRefCoord ( q ); 
            vFunc.evaluateGradient ( *element, refCoords, gradVVal );
            if ( q == 0 ) {
              gradientXVal[index] = aol::Abs(gradVVal[0]);
              gradientYVal[index] = aol::Abs(gradVVal[1]);
            }
            gradientZVal[index] += aol::Abs(gradVVal[2]);            
          }
          gradientZVal[index] *= 1.0/static_cast<RealType>(QuadRule::numQuadPoints);
          // Case 2: Element is on boundaries 
          if ( elementOnBoundary ( rp, element ) ) { 
            if ( gradientXVal[index] > 1e-3 || gradientYVal[index] > 1e-3 ) {
              gradientXVal[index] = 0.0;
              gradientYVal[index] = 0.0;
              boundaryXYVal[index] = 1.0;   
            }
            if ( gradientZVal[index] > 1e-3 ) {
              gradientZVal[index] = 0.0;
              boundaryZVal[index] = 1.0;
            }
          }
        }
        
        // Mark elements for refinement 
        RealType maxGradX = gradientXVal.getMaxAbsValue();
        RealType maxGradY = gradientYVal.getMaxAbsValue();
        RealType maxGradZ = gradientZVal.getMaxAbsValue();
        if ( maxGradX > 1e-10 ) gradientXVal *= 1.0/maxGradX;
        if ( maxGradY > 1e-10 ) gradientYVal *= 1.0/maxGradY;
        if ( maxGradZ > 1e-10 ) gradientZVal *= 1.0/maxGradZ;
        index = 0;
        bool elementWasZMarked;
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {  
          elementWasZMarked = false;
          // Check if xy gradient is large 
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( gradientXVal[index] + gradientYVal[index] > gradientTol ) {   
              if ( element.getBaseArea() > minBaseArea ) { 
                grid.markXY ( index ); 
                ++numMarkedXY;   
              }
            }
            // Element is on boundary 
            else {
              if ( boundaryXYVal[index] > 0.5 ) {
                if ( element.getBaseArea() > minBaseArea ) {
                  grid.markXY ( index );
                  ++numMarkedXY;
                }
              }
            }
          }
          // Check if z gradient is large 
          if ( !grid.isMarkedForZRefinement ( index ) ) {
            if ( gradientZVal[index] > 0.6 ) {  
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index ); 
                ++numMarkedZ; 
                elementWasZMarked = true; 
              }
            }
            // Element is on boundary 
            else {
              if ( boundaryZVal[index] > 0.5 ) {
                if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                  grid.markZ ( index );
                  ++numMarkedZ; 
                  elementWasZMarked = true; 
                }
              }
            }
          }

          // ADDITIONAL REFINEMENT: NEIGHBORS OF REFINED ELEMENT 
          if ( elementWasZMarked ) { 
            // Find neighbors 
            aol::Vec < 8, GlobalIndex > neighbors = grid.getNeighborsOfElement ( index );
            // Mark neighbors for refinement 
            for ( int j = 0; j < 6; ++j ) {
              if ( neighbors[j] >= 0 ) {
                typename RPOpType::GridType::ElementIterator neighborIt ( grid, neighbors[j] );
                if ( !grid.isMarkedForZRefinement ( neighbors[j] ) && neighborIt.getLevelz() < rp.getStopLevelz() -1 ) {
                  grid.markZ ( neighbors[j] );
                  ++numMarkedZ;
                }
                // Check neighbor's neighbours 
                aol::Vec < 8, GlobalIndex > neighborsNeighbors = grid.getNeighborsOfElement ( neighbors[j] );
                for ( int k = 0; k < 6; ++k ) {
                  if ( neighborsNeighbors[k] >= 0 ) {
                    typename RPOpType::GridType::ElementIterator neighborsNeighborIt ( grid, neighborsNeighbors[k] );
                    if ( !grid.isMarkedForZRefinement ( neighborsNeighbors[k] ) && neighborsNeighborIt.getLevelz() < rp.getStopLevelz()-1 ) {
                      grid.markZ ( neighborsNeighbors[k] );
                      ++numMarkedZ;
                    }
                    // Check level 2... 
                    aol::Vec < 8, GlobalIndex > neighborsNeighbors2 = grid.getNeighborsOfElement ( neighborsNeighbors[k] );
                    for ( int l = 0; l < 6; ++l ) {
                      if ( neighborsNeighbors2[l] >= 0 ) {
                        typename RPOpType::GridType::ElementIterator neighborsNeighborIt2 ( grid, neighborsNeighbors2[l] );
                        if ( !grid.isMarkedForZRefinement ( neighborsNeighbors2[l] ) && neighborsNeighborIt2.getLevelz() < rp.getStopLevelz()-1 ) {
                          grid.markZ ( neighborsNeighbors2[l] );
                          ++numMarkedZ;  
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          
          
  
        }
        rp.restrictVectorZConst ( v );
          
        cout << "Marked " << numMarkedXY << " elements (" << 100.0*static_cast<RealType>(numMarkedXY)/static_cast<RealType>(numElements) << "%) for xy refinement" << endl; 
        cout << "Marked " << numMarkedZ << " elements (" << 100.0*static_cast<RealType>(numMarkedZ)/static_cast<RealType>(numElements) << "%) for z refinement" << endl; 
        
        grid.refineMarkedElements();   
        
        int numRefined = grid.getNumberOfElements() - numElements;
        cout << "Refined " << numRefined << " elements (" << 100.0*static_cast<RealType>(numRefined)/static_cast<RealType>(numElements) << "%) in total" << endl; 

        return numMarkedXY + numMarkedZ;

    }
    
    #ifdef USE_MOSEK
    int computeLocalPrimalDualGap ( const string & problemType, VectorType & PDgapPrimal, VectorType & PDgapDual, ExampleType &example, VectorType & powerLookuptable, 
                                    GridType & grid, RPOpType & rp, const SparseMatrixType & matrixX, const SparseMatrixType & matrixY, const SparseMatrixType & matrixZ, const SparseMatrixType & matrixXT, const SparseMatrixType & matrixYT, const SparseMatrixType & matrixZT, VectorType & v, VectorType & phi1, VectorType & phi2, VectorType & phi3, const RealType refinementTol, const int run, const bool strong ) { 
        
        QuadRule quadrule;    
        int gridDofs = rp.getGridRef().getNumberOfDofs(); int gridNodes = grid.getNumberOfNodes();
        
        // Create finer grid and corresponding interpolator 
        AdaptiveFEPrismMesh < RealType > gridFiner ( grid );        
        gridFiner.refineAll();      
        AdaptivePrismMeshValueMap < RealType > interpolator;
        interpolator.getIndexFromNodes ( rp );
        ConfiguratorType configuratorFiner ( gridFiner );
        RPOpType rpFiner ( gridFiner, configuratorFiner, rp.getStopLevelxy(), rp.getStopLevelz() );
        int finerGridDofs = rpFiner.getGridRef().getNumberOfDofs(); int finerGridNodes = rpFiner.getGridRef().getNumberOfNodes();      
        SparseMatrixType interpolMat ( finerGridNodes, gridNodes ); 
        interpolator.makeInterpolationMatrixZConst ( rpFiner, interpolMat ); 
        
        // Define variables
        VectorType K1v ( gridDofs ); VectorType K2v ( gridDofs ); VectorType K3v ( gridDofs );
        VectorType K1vFiner ( finerGridNodes ); VectorType K2vFiner ( finerGridNodes ); VectorType K3vFiner ( finerGridNodes );
        matrixX.apply ( v, K1v ); matrixY.apply ( v, K2v ); matrixZ.apply ( v, K3v );
        VectorType vFiner ( finerGridNodes ); VectorType phi1Finer ( finerGridNodes ); VectorType phi2Finer ( finerGridNodes ); VectorType phi3Finer ( finerGridNodes );
        VectorType PDgap ( grid.getNumberOfElements() );
        
        // Interpolate variables to finer grid 
        rp.extendVectorZConst ( v );
        rp.extendVectorZConst ( phi1 ); rp.extendVectorZConst ( phi2 ); rp.extendVectorZConst ( phi3 ); 
        rp.extendVectorZConst ( K1v ); rp.extendVectorZConst ( K2v ); rp.extendVectorZConst ( K3v );
        interpolMat.apply ( v, vFiner ); interpolMat.apply ( phi1, phi1Finer ); interpolMat.apply ( phi2, phi2Finer ); interpolMat.apply ( phi3, phi3Finer );
        interpolMat.apply ( K1v, K1vFiner ); interpolMat.apply ( K2v, K2vFiner ); interpolMat.apply ( K3v, K3vFiner );
        rp.restrictVectorZConst ( v ); rp.restrictVectorZConst ( phi1 ); rp.restrictVectorZConst ( phi2 ); rp.restrictVectorZConst ( phi3 ); 
        rp.restrictVectorZConst ( K1v ); rp.restrictVectorZConst ( K2v ); rp.restrictVectorZConst ( K3v );
        rpFiner.restrictVectorZConst ( vFiner ); rpFiner.restrictVectorZConst ( phi1Finer ); rpFiner.restrictVectorZConst ( phi2Finer ); rpFiner.restrictVectorZConst ( phi3Finer ); 
        rpFiner.restrictVectorZConst ( K1vFiner ); rpFiner.restrictVectorZConst ( K2vFiner ); rpFiner.restrictVectorZConst ( K3vFiner ); 
        VectorType vFinerCorrectBounds ( vFiner );
        
        // Set finer grid matrix assemblers 
        SemiMassMatrixAssemblerX < ConfiguratorType > xAssemblerFiner ( configuratorFiner, gridFiner );
        SemiMassMatrixAssemblerY < ConfiguratorType > yAssemblerFiner ( configuratorFiner, gridFiner );
        SemiMassMatrixAssemblerZ < ConfiguratorType > zAssemblerFiner ( configuratorFiner, gridFiner );        
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get phiOpt  
        
        // Set example and lookuptables
        ExampleType exampleFiner ( example ); 
        exampleFiner.setBoundaryMask ( rpFiner );
        rpFiner.setLineSegmentList();
        rpFiner.setConstraintListFull(); 
        int constraintListSize = rpFiner.getConstraintListFull().size();
        VectorType phi1Opt ( phi1Finer ); VectorType phi2Opt ( phi2Finer ); VectorType phi3Opt ( phi3Finer );
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rpFiner.getStopLevelz() ) );
        
        // Go through list of line segments and optimize independently 
        #pragma omp parallel for 
        for ( int i = 0; i < constraintListSize; ++i ) { 
            
          int numberOfConstraints = rpFiner.getConstraintListFull()[i].size();
        
          // Extract corresponding values of Kv 
          int lineLength = rp.getLineSegmentList()[i].size();
          int k = 0, t1, t2, totalIntervalLength;
          RealType RHS; 
          VectorType zCoordsLoc ( lineLength );
          std::vector < RealType > K1vLoc ( lineLength ); std::vector < RealType > K2vLoc ( lineLength ); std::vector < RealType > K3vLoc ( lineLength ); 
          for ( std::map <int,int>::const_iterator lineit = rpFiner.getLineSegmentList()[i].begin(); lineit != rpFiner.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            K1vLoc[k] = K1vFiner[lineit->second];
            K2vLoc[k] = K2vFiner[lineit->second];
            K3vLoc[k] = K3vFiner[lineit->second];
            zCoordsLoc[k] = lineit->first;
          }
          std::shared_ptr < ndarray < RealType, 1 > > K1vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
          std::shared_ptr < ndarray < RealType, 1 > > K2vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
          std::shared_ptr < ndarray < RealType, 1 > > K3vLocPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
          std::copy ( K1vLoc.begin(), K1vLoc.end(), K1vLocPtr->begin() );
          std::copy ( K2vLoc.begin(), K2vLoc.end(), K2vLocPtr->begin() );
          std::copy ( K3vLoc.begin(), K3vLoc.end(), K3vLocPtr->begin() );
          
          // Model 
          Model::t M = new Model("primalEnergyLoc"); 
          
          // Variables 
          mosek::fusion::Variable::t phi1OptLoc = M->variable ( "phi1OptLoc", lineLength, Domain::unbounded() );
          mosek::fusion::Variable::t phi2OptLoc = M->variable ( "phi2OptLoc", lineLength, Domain::unbounded() );
          mosek::fusion::Variable::t phi3OptLoc = M->variable ( "phi3OptLoc", lineLength, Domain::inRange(0.0,1e8f) ); // ATTENTION: Needs to be bounded from above 
          
          // Constraints and additional variables 
          for ( int j = 0; j < numberOfConstraints; ++j ) {
            t1 = rpFiner.getConstraintListFull()[i][j].first;
            t2 = rpFiner.getConstraintListFull()[i][j].second;
            totalIntervalLength = zCoordsLoc[t2] - zCoordsLoc[t1];  
            RHS = pow ( hz * static_cast<RealType>(totalIntervalLength), 1.0-_epsilon ); 
            // Compute length vector 
            std::vector < RealType > lengthVec ( lineLength ); 
            for ( int kk = 0; kk < lineLength; ++kk ) {
                if ( kk >= t1 && kk < t2 ) 
                    lengthVec[kk] = static_cast < RealType > ( zCoordsLoc[kk+1] - zCoordsLoc[kk] ) * hz;
                else 
                    lengthVec[kk] = 0.0;
            }
            std::shared_ptr < ndarray < RealType, 1 > > lengthVecPtr ( new ndarray < RealType, 1 >(shape(lineLength)) );
            std::copy ( lengthVec.begin(), lengthVec.end(), lengthVecPtr->begin() );
            // Add cone constraint 
            M->constraint ( Expr::vstack ( RHS, Expr::dot(lengthVecPtr,phi1OptLoc), Expr::dot(lengthVecPtr,phi2OptLoc) ), Domain::inQCone(3) ); 
          }
          
          // Objective 
          M->objective ("obj", ObjectiveSense::Maximize, Expr::add ( Expr::add( Expr::dot(K1vLocPtr,phi1OptLoc), Expr::dot(K2vLocPtr,phi2OptLoc) ), Expr::dot(K3vLocPtr,phi3OptLoc) ) );
            
          // Solve the problem and get solution  
          M->solve();   
          M->acceptedSolutionStatus ( AccSolutionStatus::Optimal ); 
          ndarray < RealType, 1 > phi1OptLocSol = *(phi1OptLoc->level()); 
          ndarray < RealType, 1 > phi2OptLocSol = *(phi2OptLoc->level());
          ndarray < RealType, 1 > phi3OptLocSol = *(phi3OptLoc->level());

          // Write result to phiOpt (global)
          k = 0;
          for ( std::map <int,int>::const_iterator lineit = rpFiner.getLineSegmentList()[i].begin(); lineit != rpFiner.getLineSegmentList()[i].end(); ++lineit, ++k ) {
            phi1Opt[lineit->second] = phi1OptLocSol[k];
            phi2Opt[lineit->second] = phi2OptLocSol[k];
            phi3Opt[lineit->second] = phi3OptLocSol[k];
          } 
            
          // Destroy model   
          M->dispose(); 
            
        }


        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get vOpt
        
        // Set variables 
        VectorType K1Tphi1 ( gridDofs ); VectorType K2Tphi2 ( gridDofs ); VectorType K3Tphi3 ( gridDofs );
        matrixXT.apply ( phi1, K1Tphi1 ); matrixYT.apply ( phi2, K2Tphi2 ); matrixZT.apply ( phi3, K3Tphi3 );
        VectorType K1Tphi1Finer ( finerGridNodes ); VectorType K2Tphi2Finer ( finerGridNodes ); VectorType K3Tphi3Finer ( finerGridNodes );
        rp.extendVectorZConst ( K1Tphi1 ); rp.extendVectorZConst ( K2Tphi2 ); rp.extendVectorZConst ( K3Tphi3 ); 
        interpolMat.apply ( K1Tphi1, K1Tphi1Finer ); interpolMat.apply ( K2Tphi2, K2Tphi2Finer ); interpolMat.apply ( K3Tphi3, K3Tphi3Finer );
        rpFiner.restrictVectorZConst ( K1Tphi1Finer ); rpFiner.restrictVectorZConst ( K2Tphi2Finer ); rpFiner.restrictVectorZConst ( K3Tphi3Finer ); 
        VectorType vOpt ( vFiner ); VectorType vOptCheck ( vFiner ); VectorType divPhiFiner ( finerGridDofs );
        
        // Update v
        RealType energyNew, energyOld; 
        for ( int i = 0; i < finerGridDofs; ++i ) {
          divPhiFiner[i] = K1Tphi1Finer[i] + K2Tphi2Finer[i] + K3Tphi3Finer[i];
          energyOld = vFiner[i] * divPhiFiner[i];
          if ( divPhiFiner[i] > 1e-10 ) vOptCheck[i] = 0.0;
          else if ( divPhiFiner[i] < -(1e-10) ) vOptCheck[i] = 1.0;
          energyNew = vOptCheck[i] * divPhiFiner[i];
          if ( energyNew > energyOld + 1e-8 ) vOpt[i] = vOptCheck[i];
        }
        this->renewBoundaries ( vOpt, exampleFiner );
        
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get coefficient vectors 
        VectorType vBoundaryDiff ( finerGridDofs ); VectorType vDiff ( finerGridDofs );
        VectorType phi1Diff ( finerGridDofs ); VectorType phi2Diff ( finerGridDofs ); VectorType phi3Diff ( finerGridDofs );
        for ( int i = 0; i < finerGridDofs; ++i ) {
            vBoundaryDiff[i] = vFiner[i] - vFinerCorrectBounds[i];
            phi1Diff[i] = phi1Opt[i] - phi1Finer[i]; phi2Diff[i] = phi2Opt[i] - phi2Finer[i]; phi3Diff[i] = phi3Opt[i] - phi3Finer[i];
            vDiff[i] = vOpt[i] - vFiner[i];
        }
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        
        // Get discrete functions 
        rp.extendVectorZConst ( v ); rpFiner.extendVectorZConst ( vBoundaryDiff ); 
        rpFiner.extendVectorZConst ( phi1Diff ); rpFiner.extendVectorZConst ( phi2Diff ); rpFiner.extendVectorZConst ( phi3Diff );
        rpFiner.extendVectorZConst ( K1vFiner ); rpFiner.extendVectorZConst ( K2vFiner ); rpFiner.extendVectorZConst ( K3vFiner );
        rpFiner.extendVectorZConst ( vDiff ); rpFiner.extendVectorZConst ( divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFunc ( rp.getConfiguratorRef(), v );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vBoundaryDiffFunc ( rpFiner.getConfiguratorRef(), vBoundaryDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1DiffFunc ( rpFiner.getConfiguratorRef(), phi1Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2DiffFunc ( rpFiner.getConfiguratorRef(), phi2Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3DiffFunc ( rpFiner.getConfiguratorRef(), phi3Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K1vFunc ( rpFiner.getConfiguratorRef(), K1vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K2vFunc ( rpFiner.getConfiguratorRef(), K2vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K3vFunc ( rpFiner.getConfiguratorRef(), K3vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vDiffFunc ( rpFiner.getConfiguratorRef(), vDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > divPhiFunc ( rpFiner.getConfiguratorRef(), divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFinerFunc ( rpFiner.getConfiguratorRef(), vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1FinerFunc ( rpFiner.getConfiguratorRef(), phi1Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2FinerFunc ( rpFiner.getConfiguratorRef(), phi2Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3FinerFunc ( rpFiner.getConfiguratorRef(), phi3Finer );
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Compute PD gap and mark elements for refinement
        
        RealType volElement, primalGap, dualGap, valueElementPrimal, valueElementDual, valueSubelementPrimal, valueSubelementDual, valPhi1Diff, valPhi2Diff, valPhi3Diff, divPhi, valVOpt, gradVX, gradVY, gradVZ, valVDiff; 
        RealType boundaryCheckValXY, boundaryCheckValZ, gradientXVal, gradientYVal, gradientZVal, gradientTol = 0.0005, zCoordMotherEl; 
        if ( strong ) gradientTol *= 0.1;
        int index = 0, numChildren, J, L; 
        typename ConfiguratorType::DomVecType refCoords;
        VectorType boundaryCheckXY ( grid.getNumberOfElements() );
        VectorType boundaryCheckZ ( grid.getNumberOfElements() );
        VectorType gradientX ( grid.getNumberOfElements() );
        VectorType gradientY ( grid.getNumberOfElements() );
        VectorType gradientZ ( grid.getNumberOfElements() ); 
        
        typename ConfiguratorType::VecType gradVVal, gradBoundaryDiff, gradV, gradPhi1, gradPhi2, gradPhi3;
        
        // Go through all elements and compute primal dual gap 
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          valueElementPrimal = 0.0; valueElementDual = 0.0; boundaryCheckValXY = 0.0, boundaryCheckValZ = 0.0;
          volElement = rp.getConfiguratorRef().vol ( *element );
          zCoordMotherEl = rp.getGridRef().getNodeCoords ( element.getNodeIndices()[3] )[2];
            
          // Go through all subelements   
          auto gotSubelements = rpFiner.getGridRef().getSubelementList().find ( element.getIndex() );
          numChildren = gotSubelements->second.size();
          for ( int i = 0; i < numChildren; ++i ) {
              
            typename RPOpType::GridType::ElementIterator subelement ( rpFiner.getGridRef(), gotSubelements->second[i] );  
            aol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices();
            valueSubelementPrimal = 0.0; valueSubelementDual = 0.0;  
            gradientXVal = 0.0; gradientYVal = 0.0; gradientZVal = 0.0;
            
            // Dual part 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixX; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixY; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixZ; 
            xAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixX );
            yAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixY );
            zAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixZ );
            for ( int j = 0; j < 6; ++j ) {
              J = nodeIndsSubelement[j];  
              for ( int l = 0; l < 6; ++l ) {
                L = nodeIndsSubelement[l];
                valueSubelementDual += vDiff[J] * ( phi1Finer[L]*localMatrixX[j][l] + phi2Finer[L]*localMatrixY[j][l] + phi3Finer[L]*localMatrixZ[j][l] );
              }
            }
            
            // Primal part, boundaryCheck and gradient 
            for ( int q = 0; q < QuadRule::numQuadPoints; ++q ) {
              refCoords = quadrule.getRefCoord ( q ); 
              // Primal part 
              vFinerFunc.evaluateGradient ( *subelement, refCoords, gradV ); 
              valPhi1Diff = phi1DiffFunc.evaluate ( *subelement, refCoords );
              valPhi2Diff = phi2DiffFunc.evaluate ( *subelement, refCoords );
              valPhi3Diff = phi3DiffFunc.evaluate ( *subelement, refCoords ); 
              valueSubelementPrimal += ( valPhi1Diff*gradV[0] + valPhi2Diff*gradV[1] + valPhi3Diff*gradV[2] ) * quadrule.getWeight ( q );
              // Boundary check 
              vBoundaryDiffFunc.evaluateGradient ( *subelement, refCoords, gradBoundaryDiff ); 
              boundaryCheckValXY += aol::Abs ( gradBoundaryDiff[0] ) + aol::Abs ( gradBoundaryDiff[1] ); 
              boundaryCheckValZ += this->checkBoundaryConditions ( rpFiner, subelement, vBoundaryDiff, zCoordMotherEl ); 
              // Gradient 
              vFunc.evaluateGradient ( *element, refCoords, gradVVal );
              gradientXVal += gradVVal[0];
              gradientYVal += gradVVal[1];
              gradientZVal += gradVVal[2]; 
            }
            
            // Add subelement values to element values 
            valueElementPrimal += valueSubelementPrimal;
            valueElementDual += valueSubelementDual;  
            gradientX[index] = gradientXVal;
            gradientY[index] = gradientYVal; 
            gradientZ[index] = gradientZVal; 
                
          }
          
          gradientX[index] *= volElement;
          gradientY[index] *= volElement;
          gradientZ[index] *= volElement; 
          
          // Set value of primal dual gap 
          valueElementPrimal *= volElement;
          valueElementDual *= volElement;
          PDgapPrimal[index] = valueElementPrimal;  
          PDgapDual[index] = valueElementDual;
          PDgap[index] = PDgapPrimal[index] - PDgapDual[index];
          
          // Set value of boundaryCheck 
          boundaryCheckXY[index] = boundaryCheckValXY;
          boundaryCheckZ[index] = boundaryCheckValZ;
        
        }
        
        // Save pd gap 
        saveImageAsVTKCell ( grid, PDgapPrimal, "pdgapprimal", run );
        saveImageAsVTKCell ( grid, PDgapDual, "pdgapdual", run );
        saveImageAsVTKCell ( grid, PDgap, "pdgap", run );
         
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        // Refine where pd gapaol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices(); or gradient is large 
        // Get maximum pd gap 
        RealType maxPDgap = PDgap.getMaxValue();
        index = 0;
        int numMarkedXY = 0, numMarkedZ = 0;
        RealType minBaseArea = 0.5 / static_cast < RealType > ( 1 << ( 2 * ( rp.getStopLevelxy() - 1 ) ) );  
        
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          // Check pd gap   
          if ( aol::Abs(PDgap[index]) >= refinementTol*maxPDgap && maxPDgap > 1e-5 ) { 
            if ( element.getBaseArea() > minBaseArea ) {
              grid.markXY ( index );
              ++numMarkedXY; 
            }
            if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {  
              grid.markZ ( index );
              ++numMarkedZ;
            }
          }

          // Check boundary conditions: If vFiner does not equal vFinerCorrectBounds, element needs to be refined 
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( boundaryCheckXY[index] >= 1e-7 ) { 
              if ( element.getBaseArea() > minBaseArea ) {   
                grid.markXY ( index ); 
                ++numMarkedXY; 
              }
            }
            if ( boundaryCheckZ[index] >= 1e-7 ) {
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index );
                ++numMarkedZ;
              }
            }
          }
          
          // Check gradient: If gradient is large, refine element  
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( aol::Abs(gradientX[index]) + aol::Abs(gradientY[index]) > gradientTol ) {
              if ( element.getBaseArea() > minBaseArea ) {  
                grid.markXY ( index ); 
                ++numMarkedXY;
              }
            }
          }
          
          if ( strong && !grid.isMarkedForZRefinement ( index ) ) {
            if ( aol::Abs(gradientZ[index]) > gradientTol ) {
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index ); 
                ++numMarkedZ; 
              }
            }
          }
          
        }
        
        cout << "Marked " << numMarkedXY << " elements for xy refinement" << endl; 
        cout << "Marked " << numMarkedZ << " elements for z refinement" << endl; 
        
        grid.refineMarkedElements();
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        
        return numMarkedXY + numMarkedZ;
        
    }    
    #else
    int computeLocalPrimalDualGap ( const string & problemType, VectorType & PDgapPrimal, VectorType & PDgapDual, ExampleType &example, VectorType & powerLookuptable, 
                                    GridType & grid, RPOpType & rp, const SparseMatrixType & matrixX, const SparseMatrixType & matrixY, const SparseMatrixType & matrixZ, const SparseMatrixType & matrixXT, const SparseMatrixType & matrixYT, const SparseMatrixType & matrixZT, VectorType & v, VectorType & phi1, VectorType & phi2, VectorType & phi3, const RealType refinementTol, const int run, const bool strong ) { 
        
        QuadRule quadrule;    
        int gridDofs = rp.getGridRef().getNumberOfDofs(); int gridNodes = grid.getNumberOfNodes();
        
        // Create finer grid and corresponding interpolator 
        AdaptiveFEPrismMesh < RealType > gridFiner ( grid );        
        gridFiner.refineAll();      
        AdaptivePrismMeshValueMap < RealType > interpolator;
        interpolator.getIndexFromNodes ( rp );
        ConfiguratorType configuratorFiner ( gridFiner );
        RPOpType rpFiner ( gridFiner, configuratorFiner, rp.getStopLevelxy(), rp.getStopLevelz() );
        int finerGridDofs = rpFiner.getGridRef().getNumberOfDofs(); int finerGridNodes = rpFiner.getGridRef().getNumberOfNodes();      
        SparseMatrixType interpolMat ( finerGridNodes, gridNodes ); 
        interpolator.makeInterpolationMatrixZConst ( rpFiner, interpolMat ); 
        
        // Define variables
        VectorType K1v ( gridDofs ); VectorType K2v ( gridDofs ); VectorType K3v ( gridDofs );
        VectorType K1vFiner ( finerGridNodes ); VectorType K2vFiner ( finerGridNodes ); VectorType K3vFiner ( finerGridNodes );
        matrixX.apply ( v, K1v ); matrixY.apply ( v, K2v ); matrixZ.apply ( v, K3v );
        VectorType vFiner ( finerGridNodes ); VectorType phi1Finer ( finerGridNodes ); VectorType phi2Finer ( finerGridNodes ); VectorType phi3Finer ( finerGridNodes );
        VectorType PDgap ( grid.getNumberOfElements() );
        
        // Interpolate variables to finer grid 
        rp.extendVectorZConst ( v );
        rp.extendVectorZConst ( phi1 ); rp.extendVectorZConst ( phi2 ); rp.extendVectorZConst ( phi3 ); 
        rp.extendVectorZConst ( K1v ); rp.extendVectorZConst ( K2v ); rp.extendVectorZConst ( K3v );
        interpolMat.apply ( v, vFiner ); interpolMat.apply ( phi1, phi1Finer ); interpolMat.apply ( phi2, phi2Finer ); interpolMat.apply ( phi3, phi3Finer );
        interpolMat.apply ( K1v, K1vFiner ); interpolMat.apply ( K2v, K2vFiner ); interpolMat.apply ( K3v, K3vFiner );
        rp.restrictVectorZConst ( v ); rp.restrictVectorZConst ( phi1 ); rp.restrictVectorZConst ( phi2 ); rp.restrictVectorZConst ( phi3 ); 
        rp.restrictVectorZConst ( K1v ); rp.restrictVectorZConst ( K2v ); rp.restrictVectorZConst ( K3v );
        rpFiner.restrictVectorZConst ( vFiner ); rpFiner.restrictVectorZConst ( phi1Finer ); rpFiner.restrictVectorZConst ( phi2Finer ); rpFiner.restrictVectorZConst ( phi3Finer ); 
        rpFiner.restrictVectorZConst ( K1vFiner ); rpFiner.restrictVectorZConst ( K2vFiner ); rpFiner.restrictVectorZConst ( K3vFiner ); 
        VectorType vFinerCorrectBounds ( vFiner );
        
        // Set finer grid matrix assemblers 
        SemiMassMatrixAssemblerX < ConfiguratorType > xAssemblerFiner ( configuratorFiner, gridFiner );
        SemiMassMatrixAssemblerY < ConfiguratorType > yAssemblerFiner ( configuratorFiner, gridFiner );
        SemiMassMatrixAssemblerZ < ConfiguratorType > zAssemblerFiner ( configuratorFiner, gridFiner );        
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get phiOpt  
        
        // Set example and lookuptables
        ExampleType exampleFiner ( example ); 
        exampleFiner.setBoundaryMask ( rpFiner );
        rpFiner.setLineSegmentList();
        rpFiner.setConstraintListFull(); 
        VectorType phi1Opt ( phi1Finer ); VectorType phi2Opt ( phi2Finer ); VectorType phi3Opt ( phi3Finer );
        VectorType phi1OptOld ( phi1Finer ); VectorType phi2OptOld ( phi1Finer ); VectorType phi3OptOld ( phi1Finer ); 
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rpFiner.getStopLevelz() ) );
        int iter = 0; RealType error = 1.0; RealType energyOld, energy;
        RealType dt = 1.0;
        
        // Iteration 
        while ( iter < 100 && error > 1e-8 ) { 
            phi1OptOld = phi1Opt;
            phi2OptOld = phi2Opt;
            phi3OptOld = phi3Opt; 
            energyOld = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                energyOld += phi1Opt[i] * K1vFiner[i] + phi2Opt[i] * K2vFiner[i] + phi3Opt[i] * K3vFiner[i];
            }
            for ( int i = 0; i < finerGridDofs; ++i ) {
                phi1Opt[i] = phi1Opt[i] + dt * K1vFiner[i];
                phi2Opt[i] = phi2Opt[i] + dt * K2vFiner[i];
                phi3Opt[i] = phi3Opt[i] + dt * K3vFiner[i];
            }
            if ( problemType == "BranchedTransport" ) this->projectionOntoIntegralConstraintBT ( powerLookuptable, rp, phi1Opt, phi2Opt );
            else if ( problemType == "UrbanPlanning" ) this->projectionOntoIntegralConstraintUP ( rp, phi1Opt, phi2Opt );
            phi3Opt.thresholdFromBelow ( 0.0, 0.0) ;
            energy = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                energy += phi1Opt[i] * K1v[i] + phi2Opt[i] * K2v[i] + phi3Opt[i] * K3v[i];
            }
            error = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                error += aol::Abs ( phi1Opt[i] - phi1OptOld[i] ) + aol::Abs ( phi2Opt[i] - phi2OptOld[i] ) + aol::Abs ( phi3Opt[i] - phi3OptOld[i] );
            }
            error = error / static_cast<RealType>(finerGridDofs);
            ++iter;
        }    
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get vOpt
        
        // Set variables 
        VectorType K1Tphi1 ( gridDofs ); VectorType K2Tphi2 ( gridDofs ); VectorType K3Tphi3 ( gridDofs );
        matrixXT.apply ( phi1, K1Tphi1 ); matrixYT.apply ( phi2, K2Tphi2 ); matrixZT.apply ( phi3, K3Tphi3 );
        VectorType K1Tphi1Finer ( finerGridNodes ); VectorType K2Tphi2Finer ( finerGridNodes ); VectorType K3Tphi3Finer ( finerGridNodes );
        rp.extendVectorZConst ( K1Tphi1 ); rp.extendVectorZConst ( K2Tphi2 ); rp.extendVectorZConst ( K3Tphi3 ); 
        interpolMat.apply ( K1Tphi1, K1Tphi1Finer ); interpolMat.apply ( K2Tphi2, K2Tphi2Finer ); interpolMat.apply ( K3Tphi3, K3Tphi3Finer );
        rpFiner.restrictVectorZConst ( K1Tphi1Finer ); rpFiner.restrictVectorZConst ( K2Tphi2Finer ); rpFiner.restrictVectorZConst ( K3Tphi3Finer ); 
        VectorType vOpt ( vFiner ); VectorType vOptCheck ( vFiner ); VectorType divPhiFiner ( finerGridDofs );
        
        // Update v
        RealType energyNew;
        energyOld = 0.0;
        for ( int i = 0; i < finerGridDofs; ++i ) {
          divPhiFiner[i] = K1Tphi1Finer[i] + K2Tphi2Finer[i] + K3Tphi3Finer[i];
          energyOld = vFiner[i] * divPhiFiner[i];
          if ( divPhiFiner[i] > 1e-10 ) vOptCheck[i] = 0.0;
          else if ( divPhiFiner[i] < -(1e-10) ) vOptCheck[i] = 1.0;
          energyNew = vOptCheck[i] * divPhiFiner[i];
          if ( energyNew > energyOld + 1e-8 ) vOpt[i] = vOptCheck[i];
        }
        this->renewBoundaries ( vOpt, exampleFiner );
        
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get coefficient vectors 
        VectorType vBoundaryDiff ( finerGridDofs ); VectorType vDiff ( finerGridDofs );
        VectorType phi1Diff ( finerGridDofs ); VectorType phi2Diff ( finerGridDofs ); VectorType phi3Diff ( finerGridDofs );
        for ( int i = 0; i < finerGridDofs; ++i ) {
            vBoundaryDiff[i] = vFiner[i] - vFinerCorrectBounds[i];
            phi1Diff[i] = phi1Opt[i] - phi1Finer[i]; phi2Diff[i] = phi2Opt[i] - phi2Finer[i]; phi3Diff[i] = phi3Opt[i] - phi3Finer[i];
            vDiff[i] = vOpt[i] - vFiner[i];
        }
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        
        // Get discrete functions 
        rp.extendVectorZConst ( v ); rpFiner.extendVectorZConst ( vBoundaryDiff ); 
        rpFiner.extendVectorZConst ( phi1Diff ); rpFiner.extendVectorZConst ( phi2Diff ); rpFiner.extendVectorZConst ( phi3Diff );
        rpFiner.extendVectorZConst ( K1vFiner ); rpFiner.extendVectorZConst ( K2vFiner ); rpFiner.extendVectorZConst ( K3vFiner );
        rpFiner.extendVectorZConst ( vDiff ); rpFiner.extendVectorZConst ( divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFunc ( rp.getConfiguratorRef(), v );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vBoundaryDiffFunc ( rpFiner.getConfiguratorRef(), vBoundaryDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1DiffFunc ( rpFiner.getConfiguratorRef(), phi1Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2DiffFunc ( rpFiner.getConfiguratorRef(), phi2Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3DiffFunc ( rpFiner.getConfiguratorRef(), phi3Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K1vFunc ( rpFiner.getConfiguratorRef(), K1vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K2vFunc ( rpFiner.getConfiguratorRef(), K2vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K3vFunc ( rpFiner.getConfiguratorRef(), K3vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vDiffFunc ( rpFiner.getConfiguratorRef(), vDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > divPhiFunc ( rpFiner.getConfiguratorRef(), divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFinerFunc ( rpFiner.getConfiguratorRef(), vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1FinerFunc ( rpFiner.getConfiguratorRef(), phi1Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2FinerFunc ( rpFiner.getConfiguratorRef(), phi2Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3FinerFunc ( rpFiner.getConfiguratorRef(), phi3Finer );
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Compute PD gap and mark elements for refinement
        
        RealType volElement, primalGap, dualGap, valueElementPrimal, valueElementDual, valueSubelementPrimal, valueSubelementDual, valPhi1Diff, valPhi2Diff, valPhi3Diff, divPhi, valVOpt, gradVX, gradVY, gradVZ, valVDiff; 
        RealType boundaryCheckValXY, boundaryCheckValZ, gradientXVal, gradientYVal, gradientZVal, gradientTol = 0.0005, zCoordMotherEl; 
        if ( strong ) gradientTol *= 0.1;
        int index = 0, numChildren, J, L; 
        typename ConfiguratorType::DomVecType refCoords;
        VectorType boundaryCheckXY ( grid.getNumberOfElements() );
        VectorType boundaryCheckZ ( grid.getNumberOfElements() );
        VectorType gradientX ( grid.getNumberOfElements() );
        VectorType gradientY ( grid.getNumberOfElements() );
        VectorType gradientZ ( grid.getNumberOfElements() ); 
        
        typename ConfiguratorType::VecType gradVVal, gradBoundaryDiff, gradV, gradPhi1, gradPhi2, gradPhi3;
        
        // Go through all elements and compute primal dual gap 
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          valueElementPrimal = 0.0; valueElementDual = 0.0; boundaryCheckValXY = 0.0, boundaryCheckValZ = 0.0;
          volElement = rp.getConfiguratorRef().vol ( *element );
          zCoordMotherEl = rp.getGridRef().getNodeCoords ( element.getNodeIndices()[3] )[2];
            
          // Go through all subelements   
          auto gotSubelements = rpFiner.getGridRef().getSubelementList().find ( element.getIndex() );
          numChildren = gotSubelements->second.size();
          for ( int i = 0; i < numChildren; ++i ) {
              
            typename RPOpType::GridType::ElementIterator subelement ( rpFiner.getGridRef(), gotSubelements->second[i] );  
            aol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices();
            valueSubelementPrimal = 0.0; valueSubelementDual = 0.0;  
            gradientXVal = 0.0; gradientYVal = 0.0; gradientZVal = 0.0;
            
            // Dual part 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixX; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixY; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixZ; 
            xAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixX );
            yAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixY );
            zAssemblerFiner.prepareLocalMatrix ( *subelement, localMatrixZ );
            for ( int j = 0; j < 6; ++j ) {
              J = nodeIndsSubelement[j];  
              for ( int l = 0; l < 6; ++l ) {
                L = nodeIndsSubelement[l];
                valueSubelementDual += vDiff[J] * ( phi1Finer[L]*localMatrixX[j][l] + phi2Finer[L]*localMatrixY[j][l] + phi3Finer[L]*localMatrixZ[j][l] );
              }
            }
            
            // Primal part, boundaryCheck and gradient 
            for ( int q = 0; q < QuadRule::numQuadPoints; ++q ) {
              refCoords = quadrule.getRefCoord ( q ); 
              // Primal part 
              vFinerFunc.evaluateGradient ( *subelement, refCoords, gradV ); 
              valPhi1Diff = phi1DiffFunc.evaluate ( *subelement, refCoords );
              valPhi2Diff = phi2DiffFunc.evaluate ( *subelement, refCoords );
              valPhi3Diff = phi3DiffFunc.evaluate ( *subelement, refCoords ); 
              valueSubelementPrimal += ( valPhi1Diff*gradV[0] + valPhi2Diff*gradV[1] + valPhi3Diff*gradV[2] ) * quadrule.getWeight ( q );
              // Boundary check 
              vBoundaryDiffFunc.evaluateGradient ( *subelement, refCoords, gradBoundaryDiff ); 
              boundaryCheckValXY += aol::Abs ( gradBoundaryDiff[0] ) + aol::Abs ( gradBoundaryDiff[1] ); 
              boundaryCheckValZ += this->checkBoundaryConditions ( rpFiner, subelement, vBoundaryDiff, zCoordMotherEl ); 
              // Gradient 
              vFunc.evaluateGradient ( *element, refCoords, gradVVal );
              gradientXVal += gradVVal[0];
              gradientYVal += gradVVal[1];
              gradientZVal += gradVVal[2]; 
            }
            
            // Add subelement values to element values 
            valueElementPrimal += valueSubelementPrimal;
            valueElementDual += valueSubelementDual;  
            gradientX[index] = gradientXVal;
            gradientY[index] = gradientYVal; 
            gradientZ[index] = gradientZVal; 
                
          }
          
          gradientX[index] *= volElement;
          gradientY[index] *= volElement;
          gradientZ[index] *= volElement; 
          
          // Set value of primal dual gap 
          valueElementPrimal *= volElement;
          valueElementDual *= volElement;
          PDgapPrimal[index] = valueElementPrimal;  
          PDgapDual[index] = valueElementDual;
          PDgap[index] = PDgapPrimal[index] - PDgapDual[index];
          
          // Set value of boundaryCheck 
          boundaryCheckXY[index] = boundaryCheckValXY;
          boundaryCheckZ[index] = boundaryCheckValZ;
        
        }
        
        // Save pd gap 
        saveImageAsVTKCell ( grid, PDgapPrimal, "pdgapprimal", run );
        saveImageAsVTKCell ( grid, PDgapDual, "pdgapdual", run );
        saveImageAsVTKCell ( grid, PDgap, "pdgap", run );
         
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        // Refine where pd gapaol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices(); or gradient is large 
        // Get maximum pd gap 
        RealType maxPDgap = PDgap.getMaxValue();
        index = 0;
        int numMarkedXY = 0, numMarkedZ = 0;
        RealType minBaseArea = 0.5 / static_cast < RealType > ( 1 << ( 2 * ( rp.getStopLevelxy() - 1 ) ) );  
        
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          // Check pd gap   
          if ( aol::Abs(PDgap[index]) >= refinementTol*maxPDgap && maxPDgap > 1e-5 ) { 
            if ( element.getBaseArea() > minBaseArea ) {
              grid.markXY ( index );
              ++numMarkedXY; 
            }
            if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {  
              grid.markZ ( index );
              ++numMarkedZ;
            }
          }

          // Check boundary conditions: If vFiner does not equal vFinerCorrectBounds, element needs to be refined 
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( boundaryCheckXY[index] >= 1e-7 ) { 
              if ( element.getBaseArea() > minBaseArea ) {   
                grid.markXY ( index ); 
                ++numMarkedXY; 
              }
            }
            if ( boundaryCheckZ[index] >= 1e-7 ) {
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index );
                ++numMarkedZ;
              }
            }
          }
          
          // Check gradient: If gradient is large, refine element  
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( aol::Abs(gradientX[index]) + aol::Abs(gradientY[index]) > gradientTol ) {
              if ( element.getBaseArea() > minBaseArea ) {  
                grid.markXY ( index ); 
                ++numMarkedXY;
              }
            }
          }
          
          if ( strong && !grid.isMarkedForZRefinement ( index ) ) {
            if ( aol::Abs(gradientZ[index]) > gradientTol ) {
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index ); 
                ++numMarkedZ; 
              }
            }
          }
          
        }
        
        cout << "Marked " << numMarkedXY << " elements for xy refinement" << endl; 
        cout << "Marked " << numMarkedZ << " elements for z refinement" << endl; 
        
        grid.refineMarkedElements();
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        
        return numMarkedXY + numMarkedZ;
        
    }    
    #endif 
    
    int computeLocalPrimalDualGapOnly ( const string & problemType, VectorType & PDgapPrimal, VectorType & PDgapDual, ExampleType &example, VectorType & powerLookuptable, 
                                        GridType & grid, RPOpType & rp, VectorType & v, VectorType & phi1, VectorType & phi2, VectorType & phi3, const RealType refinementTol ) { 
        
        QuadRule quadrule;    
        
        // Create finer grid and corresponding interpolator 
        AdaptiveFEPrismMesh < RealType > gridFiner ( grid );  
        if ( gridFiner.getMaxLevelz() == rp.getStopLevelz() - 1 ) gridFiner.refineAllXY();
        else gridFiner.refineAll();  
        AdaptivePrismMeshValueMap < RealType > interpolator;
        interpolator.getIndexFromNodes ( rp );
        ConfiguratorType configuratorFiner ( gridFiner );
        RPOpType rpFiner ( gridFiner, configuratorFiner, rp.getStopLevelxy(), rp.getStopLevelz() );
        int finerGridDofs = rpFiner.getGridRef().getNumberOfDofs(); int finerGridNodes = rpFiner.getGridRef().getNumberOfNodes(); int gridNodes = rp.getGridRef().getNumberOfNodes(); 
        SparseMatrixType interpolMat ( finerGridNodes, gridNodes ); 
        interpolator.makeInterpolationMatrixZConst ( rpFiner, interpolMat ); 
        
        // Set finer FE matrices 
        SemiMassMatrixAssemblerX < ConfiguratorType > xAssembler ( configuratorFiner, gridFiner ); 
        SemiMassMatrixAssemblerY < ConfiguratorType > yAssembler ( configuratorFiner, gridFiner );
        SemiMassMatrixAssemblerZ < ConfiguratorType > zAssembler ( configuratorFiner, gridFiner );
        SparseMatrixType matrixX ( finerGridNodes, finerGridNodes );
        SparseMatrixType matrixY ( finerGridNodes, finerGridNodes );
        SparseMatrixType matrixZ ( finerGridNodes, finerGridNodes );
        SparseMatrixType matrixXT ( finerGridDofs, finerGridDofs );
        SparseMatrixType matrixYT ( finerGridDofs, finerGridDofs );
        SparseMatrixType matrixZT ( finerGridDofs, finerGridDofs );
        xAssembler.assembleAddMatrix ( matrixX );
        yAssembler.assembleAddMatrix ( matrixY );
        zAssembler.assembleAddMatrix ( matrixZ );
        rpFiner.applyPTPInPlaceZConstSparse ( matrixX, matrixY, matrixZ ); 
        rpFiner.fastRestrictSparseMatrix ( matrixX );
        rpFiner.fastRestrictSparseMatrix ( matrixY );
        rpFiner.fastRestrictSparseMatrix ( matrixZ );
        matrixX.transposeTo ( matrixXT );
        matrixY.transposeTo ( matrixYT );
        matrixZ.transposeTo ( matrixZT );
        
        // Define variables: phiOpt 
        VectorType vFiner ( finerGridNodes ); 
        VectorType K1vFiner ( finerGridDofs ); VectorType K2vFiner ( finerGridDofs ); VectorType K3vFiner ( finerGridDofs );
        VectorType PDgap ( grid.getNumberOfElements() );
        rp.extendVectorZConst ( v );
        interpolMat.apply ( v, vFiner );
        rp.restrictVectorZConst ( v );
        rpFiner.restrictVectorZConst ( vFiner );
        matrixX.apply ( vFiner, K1vFiner );
        matrixY.apply ( vFiner, K2vFiner );
        matrixZ.apply ( vFiner, K3vFiner );
        
        // Define variables: vOpt 
        VectorType phi1Finer ( finerGridNodes ); VectorType phi2Finer ( finerGridNodes ); VectorType phi3Finer ( finerGridNodes );
        VectorType K1Tphi1Finer ( finerGridDofs ); VectorType K2Tphi2Finer ( finerGridDofs ); VectorType K3Tphi3Finer ( finerGridDofs );
        VectorType divPhiFiner ( finerGridDofs );
        rp.extendVectorZConst ( phi1 ); rp.extendVectorZConst ( phi2 ); rp.extendVectorZConst ( phi3 );
        interpolMat.apply ( phi1, phi1Finer ); interpolMat.apply ( phi2, phi2Finer ); interpolMat.apply ( phi3, phi3Finer );
        rp.restrictVectorZConst ( phi1 ); rp.restrictVectorZConst ( phi2 ); rp.restrictVectorZConst ( phi3 );
        rpFiner.restrictVectorZConst ( phi1Finer ); rpFiner.restrictVectorZConst ( phi2Finer ); rpFiner.restrictVectorZConst ( phi3Finer );
        matrixXT.apply ( phi1Finer, K1Tphi1Finer );
        matrixYT.apply ( phi2Finer, K2Tphi2Finer );
        matrixZT.apply ( phi3Finer, K3Tphi3Finer );
        VectorType vFinerCorrectBounds ( vFiner );
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get phiOpt  
        
        // Set example and lookuptables
        ExampleType exampleFiner ( example ); 
        exampleFiner.setBoundaryMask ( rpFiner );
        rpFiner.setLineSegmentList();
        rpFiner.setConstraintListFull(); 
        VectorType phi1Opt ( phi1Finer ); VectorType phi2Opt ( phi2Finer ); VectorType phi3Opt ( phi3Finer );
        VectorType phi1OptOld ( phi1Finer ); VectorType phi2OptOld ( phi1Finer ); VectorType phi3OptOld ( phi1Finer ); 
        RealType hz = 1.0 / ( static_cast < RealType > ( 1 << rpFiner.getStopLevelz() ) );
        int iter = 0; RealType error = 1.0; RealType energyOld, energy;
        RealType dt = 1.0;
        
        // Iteration 
        while ( iter < 100 && error > 1e-8 ) { 
            phi1OptOld = phi1Opt;
            phi2OptOld = phi2Opt;
            phi3OptOld = phi3Opt; 
            energyOld = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                energyOld += phi1Opt[i] * K1vFiner[i] + phi2Opt[i] * K2vFiner[i] + phi3Opt[i] * K3vFiner[i];
            }
            for ( int i = 0; i < finerGridDofs; ++i ) {
                phi1Opt[i] = phi1Opt[i] + dt * K1vFiner[i];
                phi2Opt[i] = phi2Opt[i] + dt * K2vFiner[i];
                phi3Opt[i] = phi3Opt[i] + dt * K3vFiner[i];
            }
            if ( problemType == "BranchedTransport" ) this->projectionOntoIntegralConstraintBT ( powerLookuptable, rp, phi1Opt, phi2Opt );
            if ( problemType == "UrbanPlanning" ) this->projectionOntoIntegralConstraintUP ( rp, phi1Opt, phi2Opt );
            phi3Opt.thresholdFromBelow ( 0.0, 0.0) ;
            energy = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                energy += phi1Opt[i] * K1vFiner[i] + phi2Opt[i] * K2vFiner[i] + phi3Opt[i] * K3vFiner[i];
            }
            error = 0.0;
            for ( int i = 0; i < finerGridDofs; ++i ) {
                error += aol::Abs ( phi1Opt[i] - phi1OptOld[i] ) + aol::Abs ( phi2Opt[i] - phi2OptOld[i] ) + aol::Abs ( phi3Opt[i] - phi3OptOld[i] );
            }
            error = error / static_cast<RealType>(finerGridDofs);
            ++iter;
        }  
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get vOpt
        
        // Set variables 
        VectorType vOpt ( vFiner ); VectorType vOptCheck ( vFiner ); 
        
        // Update v
        RealType energyNew; 
        energyOld = 0.0;
        for ( int i = 0; i < finerGridDofs; ++i ) {
          divPhiFiner[i] = K1Tphi1Finer[i] + K2Tphi2Finer[i] + K3Tphi3Finer[i];
          energyOld = vFiner[i] * divPhiFiner[i];
          if ( divPhiFiner[i] > 1e-10 ) vOptCheck[i] = 0.0;
          else if ( divPhiFiner[i] < -(1e-10) ) vOptCheck[i] = 1.0;
          energyNew = vOptCheck[i] * divPhiFiner[i];
          if ( energyNew > energyOld + 1e-8 ) vOpt[i] = vOptCheck[i];
        }
        this->renewBoundaries ( vOpt, exampleFiner );
        
        //saveImageAsVTK ( gridFiner, rpFiner, vOpt, "vOpt", run );
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Get coefficient vectors 
        VectorType vBoundaryDiff ( finerGridDofs ); VectorType vDiff ( finerGridDofs );
        VectorType phi1Diff ( finerGridDofs ); VectorType phi2Diff ( finerGridDofs ); VectorType phi3Diff ( finerGridDofs );
        this->renewBoundaries ( vFinerCorrectBounds, exampleFiner );
        for ( int i = 0; i < finerGridDofs; ++i ) {
            vBoundaryDiff[i] = vFiner[i] - vFinerCorrectBounds[i];
            phi1Diff[i] = phi1Opt[i] - phi1Finer[i]; phi2Diff[i] = phi2Opt[i] - phi2Finer[i]; phi3Diff[i] = phi3Opt[i] - phi3Finer[i];
            vDiff[i] = vOpt[i] - vFiner[i];
        }
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        
        // Get discrete functions 
        rp.extendVectorZConst ( v ); rpFiner.extendVectorZConst ( vBoundaryDiff ); 
        rpFiner.extendVectorZConst ( phi1Diff ); rpFiner.extendVectorZConst ( phi2Diff ); rpFiner.extendVectorZConst ( phi3Diff );
        rpFiner.extendVectorZConst ( K1vFiner ); rpFiner.extendVectorZConst ( K2vFiner ); rpFiner.extendVectorZConst ( K3vFiner );
        rpFiner.extendVectorZConst ( vDiff ); rpFiner.extendVectorZConst ( divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFunc ( rp.getConfiguratorRef(), v );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vBoundaryDiffFunc ( rpFiner.getConfiguratorRef(), vBoundaryDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1DiffFunc ( rpFiner.getConfiguratorRef(), phi1Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2DiffFunc ( rpFiner.getConfiguratorRef(), phi2Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3DiffFunc ( rpFiner.getConfiguratorRef(), phi3Diff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K1vFunc ( rpFiner.getConfiguratorRef(), K1vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K2vFunc ( rpFiner.getConfiguratorRef(), K2vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > K3vFunc ( rpFiner.getConfiguratorRef(), K3vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vDiffFunc ( rpFiner.getConfiguratorRef(), vDiff );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > divPhiFunc ( rpFiner.getConfiguratorRef(), divPhiFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > vFinerFunc ( rpFiner.getConfiguratorRef(), vFiner );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi1FinerFunc ( rpFiner.getConfiguratorRef(), phi1Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi2FinerFunc ( rpFiner.getConfiguratorRef(), phi2Finer );
        aol::DiscreteFunctionDefault < ConfiguratorType, VectorType > phi3FinerFunc ( rpFiner.getConfiguratorRef(), phi3Finer );
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo  
        // Compute PD gap and mark elements for refinement
        RealType volElement, primalGap, dualGap, valueElementPrimal, valueElementDual, valueSubelementPrimal, valueSubelementDual, valPhi1Diff, valPhi2Diff, valPhi3Diff, divPhi, valVOpt, valVDiff; 
        RealType boundaryCheckValXY, boundaryCheckValZ, zCoordMotherEl; 
        int index = 0, numChildren, J, L; 
        typename ConfiguratorType::DomVecType refCoords;
        VectorType boundaryCheckXY ( grid.getNumberOfElements() );
        VectorType boundaryCheckZ ( grid.getNumberOfElements() );
        
        typename ConfiguratorType::VecType gradVVal, gradBoundaryDiff, gradV, gradPhi1, gradPhi2, gradPhi3;
        
        // Go through all elements and compute primal dual gap 
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          valueElementPrimal = 0.0; valueElementDual = 0.0; boundaryCheckValXY = 0.0, boundaryCheckValZ = 0.0;
          volElement = rp.getConfiguratorRef().vol ( *element );
          zCoordMotherEl = rp.getGridRef().getNodeCoords ( element.getNodeIndices()[3] )[2];
            
          // Go through all subelements   
          auto gotSubelements = rpFiner.getGridRef().getSubelementList().find ( element.getIndex() );
          numChildren = gotSubelements->second.size();
          for ( int i = 0; i < numChildren; ++i ) {
              
            typename RPOpType::GridType::ElementIterator subelement ( rpFiner.getGridRef(), gotSubelements->second[i] );  
            aol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices();
            valueSubelementPrimal = 0.0; valueSubelementDual = 0.0;  
            
            // Dual part 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixX; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixY; 
            aol::Mat <ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrixZ; 
            xAssembler.prepareLocalMatrix ( *subelement, localMatrixX );
            yAssembler.prepareLocalMatrix ( *subelement, localMatrixY );
            zAssembler.prepareLocalMatrix ( *subelement, localMatrixZ );
            for ( int j = 0; j < 6; ++j ) {
              J = nodeIndsSubelement[j];  
              for ( int l = 0; l < 6; ++l ) {
                L = nodeIndsSubelement[l];
                valueSubelementDual += vDiff[J] * ( phi1Finer[L]*localMatrixX[j][l] + phi2Finer[L]*localMatrixY[j][l] + phi3Finer[L]*localMatrixZ[j][l] );
              }
            }
            
            // Primal part, boundaryCheck and gradient 
            for ( int q = 0; q < QuadRule::numQuadPoints; ++q ) {
              refCoords = quadrule.getRefCoord ( q ); 
              // Primal part 
              vFinerFunc.evaluateGradient ( *subelement, refCoords, gradV ); 
              valPhi1Diff = phi1DiffFunc.evaluate ( *subelement, refCoords );
              valPhi2Diff = phi2DiffFunc.evaluate ( *subelement, refCoords );
              valPhi3Diff = phi3DiffFunc.evaluate ( *subelement, refCoords ); 
              valueSubelementPrimal += ( valPhi1Diff*gradV[0] + valPhi2Diff*gradV[1] + valPhi3Diff*gradV[2] ) * quadrule.getWeight ( q );
              // Boundary check 
              vBoundaryDiffFunc.evaluateGradient ( *subelement, refCoords, gradBoundaryDiff ); 
              boundaryCheckValXY += aol::Abs ( gradBoundaryDiff[0] ) + aol::Abs ( gradBoundaryDiff[1] ); 
              boundaryCheckValZ += this->checkBoundaryConditions ( rpFiner, subelement, vBoundaryDiff, zCoordMotherEl ); 
            }
            
            // Add subelement values to element values 
            valueElementPrimal += valueSubelementPrimal;
            valueElementDual += valueSubelementDual;  
    
          }
          
          // Set value of primal dual gap 
          valueElementPrimal *= volElement;
          valueElementDual *= volElement;
          PDgapPrimal[index] = valueElementPrimal;  
          PDgapDual[index] = valueElementDual;
          PDgap[index] = aol::Abs(PDgapPrimal[index] - PDgapDual[index]);
          
          // Set value of boundaryCheck 
          boundaryCheckXY[index] = boundaryCheckValXY;
          boundaryCheckZ[index] = boundaryCheckValZ;
        
          // Check if element is on lower or upper boundary 
          if ( elementOnFrontBackBoundary ( rp, element ) ) {
            boundaryCheckXY[index] += 1.0;   
          }
          
        }
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        // Refine where pd gapaol::Vec < 6, int > nodeIndsSubelement = subelement.getNodeIndices(); or gradient is large 
        // Get maximum pd gap 
        RealType maxPDgap = PDgap.getMaxValue();
        if ( maxPDgap > 1e-10 ) PDgap *= 1.0/maxPDgap;
        index = 0;
        int numMarkedXY = 0, numMarkedZ = 0;
        RealType minBaseArea = 0.5 / static_cast < RealType > ( 1 << ( 2 * ( rp.getStopLevelxy() - 1 ) ) );  
        
        for ( typename RPOpType::GridType::ElementIterator element ( rp.getGridRef() ); element.notAtEnd(); ++element, ++index ) {
          
          // Check pd gap   
          if ( PDgap[index] >= refinementTol && maxPDgap > 1e-8 ) { 
            if ( element.getBaseArea() > minBaseArea ) {
              grid.markXY ( index );
              ++numMarkedXY; 
            }
            if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {  
              grid.markZ ( index );
              ++numMarkedZ;
            }
          }
          
          // Check boundary conditions: If vFiner does not equal vFinerCorrectBounds, element needs to be refined 
          if ( !grid.isMarkedForXYRefinement ( index ) ) {
            if ( boundaryCheckXY[index] >= 1e-7 ) { 
              if ( element.getBaseArea() > minBaseArea ) {   
                grid.markXY ( index ); 
                ++numMarkedXY; 
              }
            }
            if ( boundaryCheckZ[index] >= 1e-7 ) {
              if ( element.getLevelz() < rp.getStopLevelz() - 1 ) {
                grid.markZ ( index );
                ++numMarkedZ;
              }
            }
          }
          
        }
        
        cout << "Marked " << numMarkedXY << " elements for xy refinement" << endl; 
        cout << "Marked " << numMarkedZ << " elements for z refinement" << endl; 
        
        grid.refineMarkedElements();
        
        //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo 
        
        return numMarkedXY + numMarkedZ;
        
    }   
    
};

#endif
