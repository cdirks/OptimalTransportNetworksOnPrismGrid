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

#ifndef __AHMED_H
#define __AHMED_H

/**************************************************************/
/***   Interface to the AHMED hierarchical matrix library   ***/
/**************************************************************/
#ifdef USE_EXTERNAL_AHMED

#if AHMED_INTERFACE_VERSION != 0 && AHMED_INTERFACE_VERSION != 1
#error "AHMED_INTERFACE_VERSION needs to be set to 0 or 1!"
#endif


// includes for old interface
#if AHMED_INTERFACE_VERSION == 0

#include <bemesh.h>
#include <bemOps.h>
#include <elasticFunSol.h>
#include <laplaceFunSol.h>

WARNING_OFF ( old-style-cast )
WARNING_OFF ( unused-parameter )
WARNING_OFF ( redundant-decls )

#include <H.h>
#include <dof.h>
#include <cluster_pca.h>
#include <blcluster.h>
#include <Matrix.h>

#include <aca.h>

WARNING_ON ( old-style-cast )
WARNING_ON ( unused-parameter )
WARNING_ON ( redundant-decls )

#ifdef _OPENMP
#include <omp.h>
#endif

#endif // AHMED_INTERFACE_VERSION == 0


// includes for new interface
#if AHMED_INTERFACE_VERSION == 1

#include <unistd.h>      // getpid
#include <elasticOps.h>  // includes boundary.h, bemOps.h, elasticFunSol.h

WARNING_OFF ( old-style-cast )
#include <ahmedIncludes.h>
WARNING_ON  ( old-style-cast )

#endif // AHMED_INTERFACE_VERSION == 1


/****************************************************************/
/***  Parts common to old and new interface                   ***/
/****************************************************************/
namespace bm  {

  //! \brief Class collecting all relevant parameters for AHMED
  struct AHMEDParams  {
    // all values are directly accessible
  public:

    // HMatrix approximation options
    double  clusterEta;       //!< Cluster parameter
    int     minCluster;       //!< Minimal dimension of blocks
    double  lowRankEps;       //!< Maximal relative error caused by low rank approximation of system matrix
    int     maxRank;          //!< Maximal rank of low rank approximation before using dense representation
    bool    agglomerate;      //!< Perform agglomeration of HMatrix blocks

    // Preconditioning options
    bool    precondition;     //!< Perform HLU preconditioning when using AHMED's built-in solvers
    double  precondEps;       //!< Maximal relative error caused by low rank approximation of preconditioner
    bool    preAgglomerate;   //!< Perform agglomeration for preconditioner

    // Solver options
    double  solverEps;        //!< Stopping tolerance (not squared!) when using AHMED's built-in solvers
    int     maxIter;          //!< Maximal number of iteration when using AHMED's built-in solvers

    // Output options
    bool    verbose;          //!< Verbose output and timing
    bool    analyze;          //!< Analyze the approximation quality, warning: ruins performance


    //! \brief Default constructor, assigns default values
    AHMEDParams()
    : clusterEta     ( 0.8 )
    , minCluster     ( 10 )
    , lowRankEps     ( 1e-5 )
    , maxRank        ( 1000 )
    , agglomerate    ( true )
    , precondition   ( true )
    , precondEps     ( 1e-2 )
    , preAgglomerate ( true )
    , solverEps      ( 1e-8 )
    , maxIter        ( 1000 )
    , verbose        ( false )
    , analyze        ( false )
    {}

    //! \brief Constructor reading parameters using TypedParameterParser object
    //! If certain parameters are not specified default values will be used
    explicit AHMEDParams( aol::TypedParameterParser & parser )
    : clusterEta     ( 0.8 )
    , minCluster     ( 10 )
    , lowRankEps     ( 1e-5 )
    , maxRank        ( 1000 )
    , agglomerate    ( true )
    , precondition   ( true )
    , precondEps     ( 1e-2 )
    , preAgglomerate ( true )
    , solverEps      ( 1e-8 )
    , maxIter        ( 1000 )
    , verbose        ( false )
    , analyze        ( false )
    {
      if ( parser.hasVariable( "clusterEta" ) )      parser.get ( "clusterEta",     clusterEta     );
      if ( parser.hasVariable( "minCluster" ) )      parser.get ( "minCluster",     minCluster     );
      if ( parser.hasVariable( "HEpsilon" ) )        parser.get ( "HEpsilon",       lowRankEps     );
      if ( parser.hasVariable( "maxRank" ) )         parser.get ( "maxRank",        maxRank        );
      if ( parser.hasVariable( "Hagglomerate" ) )    parser.get ( "Hagglomerate",   agglomerate    );
      if ( parser.hasVariable( "Hprecondition" ) )   parser.get ( "Hprecondition",  precondition   );
      if ( parser.hasVariable( "HLUEpsilon" ) )      parser.get ( "HLUEpsilon",     precondEps     );
      if ( parser.hasVariable( "HLUagglomerate" ) )  parser.get ( "HLUagglomerate", preAgglomerate );
      if ( parser.hasVariable( "SolverEpsilon" ) )   parser.get ( "SolverEpsilon",  solverEps      );
      if ( parser.hasVariable( "maxIter" ) )         parser.get ( "maxIter",        maxIter        );
      if ( parser.hasVariable( "Hverbose" ) )        parser.get ( "Hverbose",       verbose        );
      if ( parser.hasVariable( "Hanalyze" ) )        parser.get ( "Hanalyze",       analyze        );
    }

    //! \brief Print current options
    void print () const  {
      cout.setf( ios::boolalpha );
      cout << "Current AHMED parameters" << endl
           << "  clusterEta:      " << clusterEta << endl
           << "  minCluster:      " << minCluster << endl
           << "  HEpsilon:        " << lowRankEps << endl
           << "  maxRank:         " << maxRank << endl
           << "  Hagglomerate:    " << agglomerate << endl
           << "  Hprecondition:   " << precondition << endl
           << "  HLUEpsilon:      " << precondEps << endl
           << "  HLUagglomerate:  " << preAgglomerate << endl
           << "  SolverEpsilon:   " << solverEps << endl
           << "  maxIter:         " << maxIter << endl
           << "  Hverbose:        " << verbose << endl
           << "  Hanalyze:        " << analyze << endl << endl;
      cout.unsetf( ios::boolalpha );
    }
  };


  typedef enum {
    RENORMALIZE_DIAGONAL_NOT,
    RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT
  } RenormalizeDiagonalType;


  //! \brief Interface class for bm::Segment subclasses with collocation node in the center
  template <class SegmentType>
  class centralCollNode
#if AHMED_INTERFACE_VERSION == 0
  : public dof
#endif
  {
    private:
      const SegmentType seg;
    public:
      centralCollNode ( const SegmentType& s ) : seg ( s ) {}
      centralCollNode& operator= ( const centralCollNode & other )  {
        seg = other.seg;
      }
      unsigned getdim () const { return 2; }
      double getcenter ( const unsigned i ) const { return 0.5 * ( seg.getStart () [i] + seg.getEnd () [i] ); }
      double getradius2 () const { return 0.25 * seg.getDirection ().normSqr (); }
      const SegmentType& getSegment () const { return seg; }
#if AHMED_INTERFACE_VERSION == 0
      virtual ~centralCollNode () {}
#endif
  };


  //! \brief Interface class for bm::Segment subclasses with collocation node at the node between two neighbouring segments
  template <class SegmentType>
  class nodalCollNode
#if AHMED_INTERFACE_VERSION == 0
  : public dof
#endif
  {
    private:
      const SegmentType pseg, nseg;
    public:
      nodalCollNode ( const SegmentType& ps, const SegmentType& ns ) : pseg ( ps ), nseg ( ns ) {}
      nodalCollNode& operator= ( const nodalCollNode & other )  {
        pseg = other.pseg;
        nseg = other.nseg;
      }
      unsigned getdim () const { return 2; }
      double getcenter ( const unsigned i ) const { return nseg.getStart () [i]; }
      double getradius2 () const { return std::max ( pseg.getDirection ().normSqr (), nseg.getDirection ().normSqr () ); }
      const SegmentType& getSegment ( int i ) const { return i ? nseg : pseg; } // 0 -> previous, 1 -> next
#if AHMED_INTERFACE_VERSION == 0
      virtual ~nodalCollNode () {}
#endif
  };

}

/****************************************************************/
/***  New interface, i.e. for AHMED from git repository       ***/
/****************************************************************/
#if AHMED_INTERFACE_VERSION == 1
namespace bm  {

  //! \brief Helper class to generate dof list
  template <class ParticleType, class DofType>
  class DofHelper {};

  //! \brief Specialization for central collocation points on segments
  template <class ParticleType>
  class DofHelper<ParticleType, centralCollNode<typename ParticleType::ConstSegmentType> >  {
    private:
      typedef typename ParticleType::ConstSegmentType  SegmentType;
      typedef centralCollNode<SegmentType>             DofType;
    public:
      static int getNumberOfDofs( const bm::Boundary<ParticleType>& bnd )  {
        return bnd.getNumberOfSegments();
      }
      static DofType allocateDof( const SegmentType & prevseg, const SegmentType & nextseg )  {
        return new DofType ( nextseg );
      }
  };

  //! \brief Specialization for nodal collocation points between segments
  template <class ParticleType>
  class DofHelper<ParticleType, nodalCollNode<typename ParticleType::ConstSegmentType> >  {
    private:
      typedef typename ParticleType::ConstSegmentType  SegmentType;
      typedef nodalCollNode<SegmentType>               DofType;
    public:
      static int getNumberOfDofs( const bm::Boundary<ParticleType>& bnd )  {
        return bnd.getNumberOfPoints();
      }
      static DofType * allocateDof( const SegmentType & prevseg, const SegmentType & nextseg )  {
        return new DofType ( prevseg, nextseg );
      }
  };


  //! \brief Interface to AHMED's GeoInfoForClustering
  template <class DofType>
  class DofInfo : public AHMED::GeoInfoForClustering  {
    protected:
      DofType ** _dofs;
    public:
      DofInfo( DofType ** dofs )
      : _dofs( dofs ) {};

      //! \brief implementation for GeoInfoForClustering: read j-th coordinate of i-th degree of freedom.
      double readCoordinate(const unsigned i, const unsigned j) const
      {
        return _dofs[i]->getcenter(j);
      }
  };


  //! \brief Implementation of AHMED's MATGEN interface
  template <class SLOperatorType, class DLOperatorType, class ParticleType, class DofType >
  class EntryEvaluator {};

  //! \brief Specialization for nodal collocation points between segments
  template <class SLOperatorType, class DLOperatorType, class ParticleType>
  class EntryEvaluator< SLOperatorType, DLOperatorType, ParticleType, nodalCollNode<typename ParticleType::ConstSegmentType> >
  {
    public:
      typedef typename ParticleType::DataType  value_type;

    private:
      typedef typename ParticleType::ConstSegmentType SegmentType;
      typedef typename ParticleType::DataType         DataType;
      typedef nodalCollNode<SegmentType>              DofType;

    protected:
      const SLOperatorType   _slOp;
      const DLOperatorType   _dlOp;
      DataType               _tau;
      DofType **             _dofs;
      const aol::BitVector & _neumannIndicator;
      const AHMED::Perm &    _perm;
      bool                   _invert;
      bool                   _interior;

    public:
      EntryEvaluator( const SLOperatorType & slOp, const DLOperatorType & dlOp, DofType ** dofs,
                      const aol::BitVector & neumannIndicator, const AHMED::Perm & perm, const bool interior )
      : _slOp             ( slOp )
      , _dlOp             ( dlOp )
      , _tau              ( 0. )
      , _dofs             ( dofs )
      , _neumannIndicator ( neumannIndicator )
      , _perm             ( perm )
      , _invert           ( false )
      , _interior         ( interior )
      {}

      //! \brief Set invert flag for the neumannIndicator
      void setInvert ( bool inv )  {
        _invert = inv;
      }

      //! \brief Implementation of AHMED's MATGEN evaluation routine
      //!
      //! copy from nodalMatrixEntry, adapted to new interface
      void cmpbl( const unsigned b1, const unsigned n1, const unsigned b2, const unsigned n2, value_type* p ) const {

        DataType val = 0.;
        for ( unsigned int j = 0; j < n2; ++j ) {
          SegmentType seg2a = static_cast<DofType*> ( _dofs[ _perm.origIdx(j+b2) ] ) -> getSegment ( 0 ),
                      seg2b = static_cast<DofType*> ( _dofs[ _perm.origIdx(j+b2) ] ) -> getSegment ( 1 );
          for ( unsigned int i = 0; i < n1; ++i ) {
            // these segments (second argument to evaluateLocally) prescribe the evaluation point of the fundamental solution
            SegmentType seg1a = static_cast<DofType*> ( _dofs[ _perm.origIdx(i+b1) ] ) -> getSegment ( 0 ),
                        seg1b = static_cast<DofType*> ( _dofs[ _perm.origIdx(i+b1) ] ) -> getSegment ( 1 );

            // do XOR with invert flag
            if ( _neumannIndicator[ _perm.origIdx(j+b2) ] != _invert )  {
              val  = -_slOp.evaluateLocally ( seg2b, seg1b, false );
              val += -_slOp.evaluateLocally ( seg2a, seg1a, true );
            }
            else  {
              val  = _dlOp.evaluateLocally ( seg2b, seg1b, false );
              val += _dlOp.evaluateLocally ( seg2a, seg1a, true );
            }

            // fac not used anymore
            if ( !_interior ) val *= -1.;

            // TODO: is tau always 0 ?!
            // *** Evil! Correct only for parametric
            if ( (_tau != 0) && (seg1a.sameParticle ( seg2a ) ) ) {

              if ( seg1a == seg2a ) {
                // diagonal
                DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
                val += 2 * _tau / ( hip * him );
              } else if ( seg1b == seg2a ) {
                // 1 before 2
                DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
                aol::Vec2<DataType> ni = cornerNormal ( seg1a, seg1b ), nip = cornerNormal ( seg2a, seg2b );
                val -= 2 * _tau * ni * nip / ( hip * ( hip + him ) );
              } else if ( seg2b == seg1a ) {
                // 2 before 1
                DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
                aol::Vec2<DataType> ni = cornerNormal ( seg1a, seg1b ), nim = cornerNormal ( seg2a, seg2b );
                val -= 2 * _tau * ni * nim / ( him * ( hip + him ) );
              }
            }
            *p++ = val;
          }
        }
      }

      //! \brief Returns expected size of the matrix entries
      double scale(const unsigned /*b1*/, const unsigned /*n1*/, const unsigned /*b2*/, const unsigned /*n2*/) const {
        return 1.0;

      }
  };


  //! \brief  An H-matrix for the approximation of a BEM operator with H-LU preconditoner (copied from AHMED examples)
  //! \tparam T1  Value type for the H-matrix. May be double, float, dcomp or scomp.
  //! \tparam T2  Value type of the preconditoner. May be T1 or typename num_traits<T1>::coarse_type (lower precision type).
  //!
  //! Matrix is derived from GeHMat which provides methods for approximation etc. and Matrix which provides solvers etc.
  template <typename T1, typename T2 = T1>
  class BEMMatrix final : public AHMED::GeHmat<T1>, public AHMED::Matrix<T1>
  {
  public:
    typedef T1 value_type;
    typedef T2 precond_type;

    /// constructor
    BEMMatrix( AHMED::BlockClusterTree_sptr blTree)
    : AHMED::GeHmat<value_type>(blTree), AHMED::Matrix<value_type>(), L_(blTree), U_(blTree)
    {
      L_ = precond_type(1.0);
      U_ = precond_type(1.0);
    }

    /// additional constructor for using GeHmat's constructor from four blocks
    BEMMatrix( AHMED::GeHmat<value_type> A11, AHMED::GeHmat<value_type> A12,
               AHMED::GeHmat<value_type> A21, AHMED::GeHmat<value_type> A22 )
    : AHMED::GeHmat<value_type>(A11, A12, A21, A22)
    , AHMED::Matrix<value_type>()
    , workingCopy_ ( *this )  // convert to precond_type. needed for genLUDecomp and useTree below
    , L_( AHMED::useTree<precond_type>( workingCopy_ ) )
    , U_( AHMED::useTree<precond_type>( workingCopy_ ) )
    {
      L_ = precond_type(1.0);
      U_ = precond_type(1.0);
    }

    ~BEMMatrix() = default;

    /// initialize preconditioner with pre Coarsening
    void initPrecond(const AHMED::SvdParam &params, const bool preCoarse)
    {
      const bool succ =
      workingCopy_.genLUDecomp_noCopy(L_, U_, params, preCoarse);
      if (!succ) {
        std::cerr << "HLU decomposition not successfull." << std::endl;
        L_ = precond_type(1.0);
        U_ = precond_type(1.0);
      }
      workingCopy_.clear();  // free memory
    }

    size_t memoryInBytes_LU() const
    {
      return L_.memoryInBytes() + U_.memoryInBytes();
    }

    unsigned numRows() const override
    {
      return AHMED::GeHmat<value_type>::numRows();
    }
    unsigned numCols() const override
    {
      return AHMED::GeHmat<value_type>::numColumns();
    }

    void amux(const value_type &d, const value_type *x, value_type *y) const
    override
    {
      this->mltaVec(d, x, y);
    }

    void precond_apply(value_type *x) const override
    {
      precond_apply(x, std::is_same<value_type, precond_type>());
    }

  private:
    /// wrapper function for precond_apply with value_type == precond_type
    void precond_apply(value_type *x, const std::true_type) const
    {
      HLU_solve(L_, U_, x);
    }

    /// wrapper function for precond_apply with value_type != precond_type
    void precond_apply(value_type *x, const std::false_type) const
    {
      auxVec_.assign(x, x + numRows());
      HLU_solve(L_, U_, auxVec_.data());
      std::copy_n(auxVec_.data(), numRows(), x);
    }

    //----------------------------------------------------------------
    AHMED::GeHmat<precond_type> workingCopy_; //!< working copy for genLUDecomp_noCopy
    AHMED::LtHmat<precond_type> L_;           //!< lower triangular H-matrix of lu decomp of this
    AHMED::UtHmat<precond_type> U_;           //!< upper triangular H-matrix of lu decomp of this

    /**
     *  @brief Auxilliary vector for precond_apply() if
     *    precond_type != value_type.
     */
    mutable std::vector<precond_type> auxVec_;
  };


  //! \brief Base class for next generation HMatrix Ops, doing everything needed to use AHMED for quocmesh BEM problems
  //! \tparam  ParticleType  Type of particle used in bm::boundary (only tested with ParaParticle so far)
  //! \tparam  DofType       Either nodal collocation nodes between segments (nodalCollNode) or central collocation nodes on segments (centralCollNode)
  //! \tparam  SLOpType      Single layer operator acting on Neumann values
  //! \tparam  DLOpType      Double layer operator acting on Dirichlet values
  template <class ParticleType, class DofType, class SLOpType, class DLOpType>
  class HMatrixOpBase : public aol::Op< aol::Vector <typename ParticleType::DataType> >
  {
    protected:
      typedef typename ParticleType::DataType  RealType;
      typedef BEMMatrix<RealType,float>        LhsHMatrixType; //!< float used for preconditioning
      typedef AHMED::GeHmat<RealType>          RhsHMatrixType;

      const AHMEDParams             _par;
      const int                     _dim;
      int                           _numDofs;
      DofType **                    _dofs;
      AHMED::Perm                   _indexPerm;
      AHMED::cluster_sptr           _clTreePtr;
      AHMED::BlockClusterTree_sptr  _blclTreePtr;
      LhsHMatrixType *              _lhsHMatrix;
      RhsHMatrixType *              _rhsHMatrix;

      //! \brief Constructor (protected, derived objects have to be used)
      HMatrixOpBase( const int dim, const bm::Boundary<ParticleType>& bnd, const aol::BitVector & neumannIndicator, const AHMEDParams & par )
      : _par         ( par )
      , _dim         ( dim )
      , _numDofs     ( 0 )
      , _dofs        ( NULL )
      , _indexPerm   ()
      , _clTreePtr   ()
      , _blclTreePtr ()
      , _lhsHMatrix  ( NULL )
      , _rhsHMatrix  ( NULL )
      {
        _numDofs = DofHelper<ParticleType,DofType>::getNumberOfDofs( bnd );
        _dofs = new DofType * [_numDofs];
        _indexPerm.initIdentity( _numDofs );

        // construct dofs and list them in _dofs
        int i = 0, n = 0;
        typename std::list<ParticleType>::const_iterator it = bnd.begin (), end = bnd.end ();
        for ( ; it != end; ++it ) {

          typename ParticleType::ConstSegmentIteratorType sit = it->beginSegment (), send = it-> endSegment (), psit = send; --psit;
          for ( ; sit != send; psit = sit, ++sit, i++ ) {

            // list dofs in sequential order
            _dofs[i] = DofHelper<ParticleType,DofType>::allocateDof ( *psit, *sit );

            // we only need to permute the indicies to have Neumann dofs in the first section and Dirichlet dofs in the second
            if ( neumannIndicator[i] )  {
              _indexPerm.swap_index( i, n );
              n++;
            }
          }
        }

        // if numNeumann = 0 or neumannIndicator.size() the cluster tree should not be separated a priori
        // the case 0 is handled correctly by genClusterTree2d_pca, so use 0 also in the other case
        int numNeumann = neumannIndicator.numTrue();
        if ( numNeumann == neumannIndicator.size() ) numNeumann = 0;

        // NOTE the sepaerated cluster tree leads to dense rectangular blocks in the H matrix which make HLU preconditioning fail
        //      it is possible to use non separated cluster trees but then blocks with mixed operators have to be approximated
//         numNeumann = 0;

        // let AHMED generate the cluster tree
        _clTreePtr = AHMED::genClusterTree2d_pca( std::make_shared< DofInfo<DofType> > ( _dofs ), _numDofs, _par.minCluster, _indexPerm, /*sep*/ numNeumann );
        if ( _par.verbose )  {
          cout << "Cluster tree:        ";
          AHMED::displayClusterInfo( _clTreePtr.get() );
        }
        if ( _par.analyze )  {
          cout << "Cluster leaves:" << endl;
          printClusterLeaves( _clTreePtr.get() );
          cout << endl;
        }

        // let AHMED generate the block cluster tree
        // the same cluster tree is used for dofs and collocation points
        // however the resulting matrix will not be symmetric, so do not use the symmetric version
        _blclTreePtr = AHMED::genBlockClusterTree( _clTreePtr, _clTreePtr, _par.clusterEta  /*, optional max depth set to numeric limits max by default */ );
        if ( _par.verbose )  {
          cout << "Block cluster tree:  ";
          AHMED::displayInfo( *_blclTreePtr );
        }
      }


      //! \brief Destructor
      ~HMatrixOpBase () {
        for ( int i = 0; i < _numDofs; ++i ) delete _dofs [i];
        delete [] _dofs;
        delete _lhsHMatrix;
        delete _rhsHMatrix;
      }


    public:
      //! \brief Print the partition of indices in generated cluster tree
      void printClusterLeaves( const AHMED::cluster* cl ) {
        cout << "Index " << cl->getidx() << ",  Size " << cl->size() << endl;
        if ( cl->isleaf() ) {
          cout << "  Leaf with indicies ";
          for ( unsigned int i = cl->getnbeg(); i < cl->getnend(); ++i )
            cout << _indexPerm.origIdx(i) << " ";
          cout << endl;
        }
        else
          for ( int i = 0; i < 2; ++i )
            printClusterLeaves( cl->getson(i) );
      }

      //! \brief Map quocmesh vector to permuted AHMED vector
      void mapBtoH ( const aol::Vector<typename ParticleType::DataType>& bvec, aol::Vector<typename ParticleType::DataType>& hvec ) const {
        int i = 0, j, qns = _numDofs * _dim;
        for ( j = 0; j < qns; j += _numDofs )
          for ( i = 0; i < _numDofs; ++i )
            hvec [ j + _indexPerm.permIdx(i) ] = bvec [j + i];
        // remaining entries
        for ( i = j; i < bvec.size (); ++i ) hvec [i] = bvec [i];
      }

      //! \brief Map permuted AHMED vector to quocmesh vector
      void mapHtoB ( const aol::Vector<typename ParticleType::DataType>& hvec, aol::Vector<typename ParticleType::DataType>& bvec ) const {
        int i = 0, j, qns = _numDofs * _dim;
        for ( j = 0; j < qns; j += _numDofs )
          for ( i = 0; i < _numDofs; ++i )
            bvec [ j + _indexPerm.origIdx(i) ] = hvec [j + i];
        // remaining entries
        for ( i = j; i < bvec.size (); ++i ) bvec [i] = hvec [i];
      }

      //! \brief Apply add right hand side op
      void applyAddRhsOp( const aol::Vector<RealType> & arg, aol::Vector<RealType> & dest, RealType factor = 1. ) const  {
        _rhsHMatrix->mltaVec( factor, arg.getData(), dest.getData() );
      }

      //! \brief Apply right hand side op
      void applyRhsOp( const aol::Vector<RealType> & arg, aol::Vector<RealType> & dest, RealType factor = 1. ) const  {
        dest.setZero();
        applyAddRhsOp( arg, dest, factor );
      }

      //! \brief Apply add left hand side op
      void applyAddLhsOp( const aol::Vector<RealType> & arg, aol::Vector<RealType> & dest, RealType factor = 1. ) const  {
        _lhsHMatrix->mltaVec( factor, arg.getData(), dest.getData() );
      }

      //! \brief Apply left hand side op
      void applyLhsOp( const aol::Vector<RealType> & arg, aol::Vector<RealType> & dest, RealType factor = 1. ) const  {
        dest.setZero();
        applyAddLhsOp( arg, dest, factor );
      }

      //! \brief Implements aol::Op interface, to use this for solving a linear system of equations the lhs op is used
      void applyAdd( const aol::Vector<RealType> & arg, aol::Vector<RealType> & dest ) const  {
        applyAddLhsOp( arg, dest, 1. );
      }


      //! \brief Solve using AHMED's GMRes with preconditioning (if requested)
      //! \param[in]  rhs      right hand side
      //! \param[out] sol      solution
      //! \param[in]  restart  number of iterations after which a restart in the GMRes algorithm is triggered
      void GMResSolve( const aol::Vector<RealType> & rhs, aol::Vector<RealType> & sol, const unsigned restart = 100 )
      {
        AHMED::RealTimer * timer = NULL;
        if ( this->_par.verbose )  {
          timer = new AHMED::RealTimer;
          timer->start();
          cout << "\nSolving system with GMRes ... " << flush;
        }

        const auto info = AHMED::GMRes( *_lhsHMatrix, this->_par.solverEps, this->_par.maxIter, restart, rhs.getData(), sol.getData() );

        if ( this->_par.verbose )  {
          double time = timer->current();
          delete timer;
          if ( info.first > this->_par.solverEps )
            cout << "did not converge!\n";
          else
            cout << "done.\n";
          cout << "Accuracy: " << info.first << ",  Steps: " << info.second << endl;
          cout << "Time:    " << time << " s -- " << 1000. * time/(2*this->_numDofs) << " ms/N." << endl;
        }
        else  {
          if ( info.first > this->_par.solverEps )  {
            cerr << "GMRes did not converge!\n";
            cerr << "Accuracy: " << info.first << ",  Steps: " << info.second << endl;
          }
        }
      }


      //! \brief Solve using AHMED's BiCGStab with preconditioning (if requested)
      //! \param[in]  rhs      right hand side
      //! \param[out] sol      solution
      void BiCGStabSolve( const aol::Vector<RealType> & rhs, aol::Vector<RealType> & sol )
      {
        AHMED::RealTimer * timer = NULL;
        if ( this->_par.verbose )  {
          timer = new AHMED::RealTimer;
          timer->start();
          cout << "\nSolving system with BiCGStab ... " << flush;
        }

        RealType outeps = this->_par.solverEps;
        unsigned int outiter = this->_par.maxIter;
        const unsigned int status = AHMED::BiCGStab( *_lhsHMatrix, rhs.getData(), sol.getData(), outeps, outiter );

        if ( this->_par.verbose )  {
          double time = timer->current();
          delete timer;
          if ( status != 0 )
            cout << "did not converge!\n";
          else
            cout << "done.\n";
          cout << "Accuracy: " << outeps << ",  Steps: " << outiter << endl;
          cout << "Time:    " << time << " s -- " << 1000. * time/(2*this->_numDofs) << " ms/N." << endl;
        }
        else  {
          if ( status != 0 )  {
            cerr << "BiCGStab did not converge!\n";
            cerr << "Accuracy: " << outeps << ",  Steps: " << outiter << endl;
          }
        }
      }
  };

  //! \brief Derived next generation HMatrix Op for problems in elasticity
  //! \tparam  ParticleType  Type of particle used in bm::boundary (tested with ParaParticle)
  //! \tparam  DofType       Either nodal collocation nodes between segments (nodalCollNode) or central collocation nodes on segments (centralCollNode)
  //! \tparam  SLOpType      Single layer operator acting on Neumann values
  //! \tparam  DLOpType      Double layer operator acting on Dirichlet values
  //!
  //! This operator employs the AHMED library to solve boundary integral equations of the type
  //! \f[ u = U[t] - V[u] + F u, \text{ or } u = -U[t] + V[u] + ( \Id - F ) u \f]
  //! for the inner or the outer problem respectively. Here u is the displacement and t the normal tension on the boundary.<br/>
  //! The boundary mesh is given by a \a bm::Boundary object, U is a single layer operator and V a double layer operator
  //! given by the templates SLOpType and DLOpType. The double layer operator is supposed to also compute and add F.
  //! If this cannot be done (correctly) \a RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT can be used to correct the diagonal entries.
  //! (note that this will generate an additional \f$ \Id \f$ later applied to u)<br/>
  //! A vector of bools \a neumannIndicator is passed to identify degrees of freedom where tensions (Neumann values)
  //! are to be computed. It is assumed that displacements (Dirichlet values) are to be computed on the remaining dofs.
  //! According to this the integral equation is rearranged resulting in a right hand side operator acting on known values
  //! which can be applied using \a HMatrixOpBase::applyRhsOp and a left hand side operator acting on unknown values which
  //! can be applied using \a HMatrixOpBase::applyLhsOp or used to solve the system of equations.
  template <class ParticleType, class DofType, class SLOpType, class DLOpType>
  class ElasticHMatrixOp : public HMatrixOpBase< ParticleType, DofType, SLOpType, DLOpType >
  {
    private:
      typedef typename HMatrixOpBase< ParticleType, DofType, SLOpType, DLOpType >::RealType        RealType;
      typedef typename HMatrixOpBase< ParticleType, DofType, SLOpType, DLOpType >::LhsHMatrixType  LhsHMatrixType;
      typedef typename HMatrixOpBase< ParticleType, DofType, SLOpType, DLOpType >::RhsHMatrixType  RhsHMatrixType;
      typedef EntryEvaluator< SLOpType, DLOpType, ParticleType, DofType >                          EvalType;

    protected:
      const int _dim;             //!< Dimension, currently always 2
      ElasticHMatrixOp() {};      //!< Default constructor forbidden

    public:
      //! \brief Constructor to be used
      //! \param[in] bnd               Quocmesh's bm boundary mesh
      //! \param[in] neumannIndicator  Vector of bools, true iff DOF i is an unknown Neumann value
      //! \param[in] green             bm::ElasticGreen fundamental solution
      //! \param[in] par               AHMEDParams struct containing all relevant parameters for AHMED
      //! \param[in] renorm            (optional) Toggle correction of double layer operator's diagonal entries
      //! \param[in] interior          (optional) Toggle interior problem
      ElasticHMatrixOp( const bm::Boundary<ParticleType>& bnd, const aol::BitVector & neumannIndicator, const bm::ElasticGreen<RealType> & green, const AHMEDParams & par,
                        RenormalizeDiagonalType renorm = RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT, bool interior = true
                      )
      : HMatrixOpBase< ParticleType, DofType, SLOpType, DLOpType > ( 2, bnd, neumannIndicator, par )
      , _dim          ( 2 )
      {
        AHMED::RealTimer timer;
        // allocate blocks of left and right hand side matrices
        RhsHMatrixType *** lhsBlocks = new RhsHMatrixType ** [_dim];
        RhsHMatrixType *** rhsBlocks = new RhsHMatrixType ** [_dim];
        for ( int i = 0; i < _dim; ++i )  {
          lhsBlocks[i] = new RhsHMatrixType * [_dim];
          rhsBlocks[i] = new RhsHMatrixType * [_dim];
          for ( int j = 0; j < _dim; ++j )  {
            lhsBlocks[i][j] = new RhsHMatrixType ( this->_blclTreePtr );  // not a typo!
            rhsBlocks[i][j] = new RhsHMatrixType ( this->_blclTreePtr );
          }
        }

        // approximate blocks
        const AHMED::SvdParam svdParam ( this->_par.lowRankEps, this->_par.maxRank, AHMED::relative, AHMED::Euclidean );
        const AHMED::AcaRowApproximator< EvalType > approx ( svdParam );
        // template parameter seems to be used for value_type only, so use the same for all blocks

        for ( int i = 0; i < _dim; ++i )  {
          for ( int j = 0; j < _dim; ++j )  {
            const SLOpType slOp ( green, i, j );
            const DLOpType dlOp ( green, i, j );
            EvalType eval ( std::move( slOp ), std::move( dlOp ), this->_dofs, neumannIndicator, this->_indexPerm, interior );

            // do not invert the neumannIndicator
            eval.setInvert( false );

            // approximate lhs block
            if ( this->_par.verbose ) timer.restart();
            lhsBlocks[i][j]->approximate( eval, approx );
            cerr << endl;
            if ( this->_par.verbose ) AHMED::io::displayInfo( lhsBlocks[i][j]->memoryInBytes(), this->_numDofs, timer.current(), sizeof(RealType) );
            if ( this->_par.analyze ) lhsBlocks[i][j]->displayApproximationError( eval );

            // agglomerate lhs block
            if ( this->_par.agglomerate )  {
              if ( this->_par.verbose ) timer.restart();
              lhsBlocks[i][j]->agglomerate( svdParam );
              if ( this->_par.verbose )  {
                cout << "Agglomerated:" << endl;
                AHMED::io::displayInfo( lhsBlocks[i][j]->memoryInBytes(), this->_numDofs, timer.current(), sizeof(RealType) );
              }
              if ( this->_par.analyze ) lhsBlocks[i][j]->displayApproximationError( eval );
            }


            // invert the neumannIndicator and start over
            eval.setInvert( true );

            // approximate rhs block
            if ( this->_par.verbose ) timer.restart();
            rhsBlocks[i][j]->approximate( eval, approx );
            cerr << endl;
            if ( this->_par.verbose ) AHMED::io::displayInfo( rhsBlocks[i][j]->memoryInBytes(), this->_numDofs, timer.current(), sizeof(RealType) );
            if ( this->_par.analyze ) rhsBlocks[i][j]->displayApproximationError( eval );

            // agglomerate rhs block
            if ( this->_par.agglomerate )  {
              if ( this->_par.verbose ) timer.restart();
              rhsBlocks[i][j]->agglomerate( svdParam );
              if ( this->_par.verbose )  {
                cout << "Agglomerated:" << endl;
                AHMED::io::displayInfo( lhsBlocks[i][j]->memoryInBytes(), this->_numDofs, timer.current(), sizeof(RealType) );
              }
              if ( this->_par.analyze ) rhsBlocks[i][j]->displayApproximationError( eval );
            }
          }
        }

        // correct diagonal entries for double layer operator
        if ( renorm == RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT ) {
          if ( this->_par.verbose ) timer.restart();

          // for diagonal CRSMatrix: will be used for both column and row indices, therefore size n+1
          std::vector<unsigned int> diagInd ( this->_numDofs+1 );
          std::iota( diagInd.begin(), diagInd.end(), 0 );

          const int numNeumann = neumannIndicator.numTrue();
          aol::Vector<RealType> onesNeumann ( this->_numDofs ), onesDirichlet ( this->_numDofs );
          // should do the same as querying neumannIndicator[ _indexPerm.origIdx(i) ]
          for ( int i = 0; i < numNeumann; ++i )              onesNeumann[i]   = 1.;
          for ( int i = numNeumann; i < this->_numDofs; ++i ) onesDirichlet[i] = 1.;

          for ( int i = 0; i < _dim; ++i )  {
            for ( int j = 0; j < _dim; ++j )  {
              // problem is that the double layer operator is partly included in rhsBlocks and lhsBlocks
              aol::Vector<RealType> resN ( this->_numDofs );
              // res should contain the result of the double layer operator being applied to ones (times -1)
              rhsBlocks[i][j]->mltaVec( -1., onesNeumann.getData(),   resN.getData() );
              lhsBlocks[i][j]->mltaVec( -1., onesDirichlet.getData(), resN.getData() );

              // if we are on the diagonal and face an exterior problem subtract the identity
              if ( ( !interior ) && ( i == j ) )  {
                resN -= onesDirichlet;
                resN -= onesNeumann;
              }

              // split the entries
              aol::Vector<RealType> resD ( resN );
              resN *= onesNeumann;
              resD *= onesDirichlet;

              // now subtract, svdParam controls truncation which should not occur for diagonal blocks anyway
              rhsBlocks[i][j]->addCRS ( resN.getData(), diagInd.data(), diagInd.data(), svdParam );
              cerr << endl;
              lhsBlocks[i][j]->addCRS ( resD.getData(), diagInd.data(), diagInd.data(), svdParam );
              cerr << endl;
            }
          }

          if ( this->_par.verbose ) {
            double time = timer.current();
            cout << "Renormalize diagonal:" << endl
                 << "Time:    " << time << " s -- " << 1000. * time/(2*this->_numDofs) << " ms/N." << endl;
          }
        }

        // merge blocks
        this->_lhsHMatrix = new LhsHMatrixType ( std::move( *lhsBlocks[0][0] ), std::move( *lhsBlocks[0][1] ), std::move( *lhsBlocks[1][0] ), std::move( *lhsBlocks[1][1] ) );
        this->_rhsHMatrix = new RhsHMatrixType ( std::move( *rhsBlocks[0][0] ), std::move( *rhsBlocks[0][1] ), std::move( *rhsBlocks[1][0] ), std::move( *rhsBlocks[1][1] ) );

        if ( this->_par.verbose )  {
          stringstream pid;
          pid << getpid();
          string filename = "lhsHMatrix_" + pid.str() + ".ps";
          cout << "Writing visualization to " << filename << endl;
          this->_lhsHMatrix->visualizeWithPostscript( filename );
          filename = "rhsHMatrix_" + pid.str() + ".ps";
          cout << "Writing visualization to " << filename << endl;
          this->_rhsHMatrix->visualizeWithPostscript( filename );
        }

        // preconditioning
        const AHMED::SvdParam svdParamPrecond ( this->_par.precondEps, this->_par.maxRank, AHMED::relative, AHMED::Euclidean );
        if ( this->_par.precondition )  {
          if ( this->_par.verbose ) timer.restart();
          this->_lhsHMatrix->initPrecond( svdParamPrecond, this->_par.preAgglomerate );
          if ( this->_par.verbose ) { cout << "Preconditioner:" << endl; AHMED::io::displayInfo( this->_lhsHMatrix->memoryInBytes_LU(), 2*this->_numDofs, timer.current(), sizeof(float) ); }
        }

        // delete blocks
        for ( int i = 0; i < _dim; ++i ) {
          for ( int j = 0; j < _dim; ++j ) {
            // the internet says that the moved-from object is left in a valid but unspecified state
            // and that it should be destroyed
            delete lhsBlocks[i][j];
            delete rhsBlocks[i][j];
          }
          delete [] lhsBlocks [i];
          delete [] rhsBlocks [i];
        }
        delete [] lhsBlocks;
        delete [] rhsBlocks;
      }
  };
}
#endif


/****************************************************************/
/***  Old interface, i.e. for AHMED subversion revision 174   ***/
/****************************************************************/
#if AHMED_INTERFACE_VERSION == 0

#ifdef USE_PTHREADS
#include <pthread.h>

namespace bm {
  template <class EntryType> class HMatrixOp;
}

namespace {
  template <class DataType>
  struct HMatrixVectorMultiplicationParameter {
    blcluster*         tree1;
    blcluster*         tree2;
    mblock<DataType>** blocks;
    DataType*          arg1;
    DataType*          arg2;
    DataType*          dest;
  };
  template <class DataType>
  void* HMatrixVectorMultiplicationThread ( void* par ) {
    HMatrixVectorMultiplicationParameter<DataType>* p = reinterpret_cast<HMatrixVectorMultiplicationParameter<DataType>*> (par);
    multa_H_vec ( 1, p->tree1, p->blocks, p->arg1, p->dest );
    multa_H_vec ( 1, p->tree2, p->blocks, p->arg2, p->dest );
    return NULL;
  }
  template <class EntryType>
  struct HMatrixApproximationParameter {
    int                       numblock;
    int                       start;
    int                       step;
    bm::HMatrixOp<EntryType>* op;
    blcluster**               blocklist;
    const EntryType*          entries;
    dof**                     segments;
    double                    epsilon;
    int                       maxrank;
  };
  template <class EntryType>
  void* HMatrixApproximationThread ( void* par ) {
    HMatrixApproximationParameter<EntryType>* p = reinterpret_cast<HMatrixApproximationParameter<EntryType>*> (par);
    for ( int i = p->start; i < p->numblock; i += p->step ) {
      //cerr << pthread_self () << " -> " << i << endl;
      blcluster* block = p->blocklist [i];
      mblock<typename EntryType::DataType>*& matrixblock = p->op->_matrix.blcks [ block -> idx ];
      p->op->thread ( * ( p->entries ), matrixblock, block, p->segments + block -> b1, p->segments + block -> b2, p->epsilon, p->maxrank );
    }
    return NULL;
  }
  template <class EntryType>
  struct HMatrixApproximationParameterMulti {
    int                       numblock;
    int                       n;
    int                       start;
    int                       step;
    bm::HMatrixOp<EntryType>* op;
    blcluster**               blocklist;
    const EntryType***        entries;
    dof**                     segments;
    double                    epsilon;
    int                       maxrank;
  };
  template <class EntryType>
  void* HMatrixApproximationThreadMulti ( void* par ) {
    HMatrixApproximationParameterMulti<EntryType>* p = reinterpret_cast<HMatrixApproximationParameterMulti<EntryType>*> (par);
    // cerr << "[";
    for ( int i = p->start; i < 4 * p->numblock; i += p->step ) {
      blcluster* block = p->blocklist [i];
      mblock<typename EntryType::DataType>*& matrixblock = p->op->_matrix.blcks [ block -> idx ];
      int off1 = block -> b1, off2 = block -> b2;
      int o1 = off1 / p->n, o2 = off2 / p->n;
      off1 %= p->n; off2 %= p->n;
      p->op->thread ( * ( p->entries [o1] [o2] ), matrixblock, block, p->segments + off1, p->segments + off2, p->epsilon, p->maxrank );
    }
    // cerr << "]";
    return NULL;
  }
}
#endif

namespace bm {

template <class ParticleType>
class centralClusterTree {
private:
  dof** segs;
  cluster2d_pca* tree;
  blcluster* blocktree;
  int numsegs;
  unsigned int numclus, numblock;

  int getDoF ( int i ) const { return segs [i] -> idx; }

public:
  centralClusterTree ( const bm::Boundary<ParticleType>& bnd, double clustereta = 1.1, int minclussize = 4,
                       bool constraint = false ) {

    numsegs = bnd.getNumberOfSegments ();
    segs = new dof* [numsegs];

    bm::ConstAllSegmentIterator<ParticleType> it = bnd.beginSegment ();
    for ( int i = 0; i < numsegs; ++i, ++it ) {
      segs [i] = new bm::centralCollNode<typename ParticleType::ConstSegmentType> ( *it );
      segs [i]->idx = segs [i]->vidx = i;
    }

    tree = new cluster2d_pca ( segs, 0, numsegs );

    tree->subdivide ( segs, minclussize, numclus );
    blocktree = genblcltree ( tree, tree, clustereta, numblock );

    // Add one level on top for constraint
    if ( constraint ) {

      blcluster* consrow = new blcluster;
      blcluster* lagrcol = new blcluster;
      blcluster* brentry = new blcluster;
      blcluster* topnode = new blcluster;

      // Constraint row
      consrow->b1 = numsegs;
      consrow->b2 = 0;
      consrow->n1 = 1;
      consrow->n2 = numsegs;
      consrow->pos  = -1;
      consrow->ccl  = tree;
      consrow->idx  = numblock++;
      consrow->adm  = true;

      // Lagrange multiplier column
      lagrcol->b1 = 0;
      lagrcol->b2 = numsegs;
      lagrcol->n1 = numsegs;
      lagrcol->n2 = 1;
      lagrcol->pos  = 1;
      lagrcol->rcl  = tree;
      lagrcol->idx  = numblock++;
      lagrcol->adm  = true;

      // Bottom right entry
      brentry->b1 = numsegs;
      brentry->b2 = numsegs;
      brentry->n1 = 1;
      brentry->n2 = 1;
      brentry->pos  = 0;
      brentry->idx  = numblock++;
      brentry->adm  = true;

      // New top node
      topnode->b1 = 0;
      topnode->b2 = 0;
      topnode->n1 = numsegs + 1;
      topnode->n2 = numsegs + 1;
      topnode->pos = 0;
      topnode->adm = false;
      topnode->son1 = blocktree;
      topnode->son2 = lagrcol;
      topnode->son3 = consrow;
      topnode->son4 = brentry;

      blocktree = topnode;
    }
  }

  ~centralClusterTree () {
    for ( int i = 0; i < numsegs; ++i ) delete segs [i];
    delete [] segs;
    delete tree;
    delete blocktree;
  }

dof** getDoFs () const { return segs; }
  int getNumDoFs () const { return numsegs; }
  blcluster* getBlockTree () const { return blocktree; }
  int getNumBlocks () const { return numblock; }
  void mapBtoH ( const aol::Vector<typename ParticleType::DataType>& bvec,
                 aol::Vector<typename ParticleType::DataType>& hvec ) const {
    int i;
      for ( i = 0; i < numsegs; ++i )
  hvec [i] = bvec [ getDoF ( i ) ];
    for ( i = i; i < bvec.size (); ++i ) hvec [i] = bvec [i];
  }
  void mapHtoB ( const aol::Vector<typename ParticleType::DataType>& hvec,
                 aol::Vector<typename ParticleType::DataType>& bvec ) const {
    int i;
      for ( i = 0; i < numsegs; ++i )
  bvec [ getDoF ( i ) ] = hvec [i];
    for ( i = i; i < bvec.size (); ++i ) bvec [i] = hvec [i];
  }
};

template <class ParticleType>
class nodalClusterTree {
private:
  dof** segs;
  cluster2d_pca* tree;
  blcluster* blocktree;
  int numsegs;
  unsigned int numclus, numblock;
  int quadruplications;

  int getDoF ( int i ) const { return segs [i] -> idx; }

public:
  nodalClusterTree ( const bm::Boundary<ParticleType>& bnd, double clustereta = 1.1, int minclussize = 4,
         bool constraint = false ) {

    quadruplications = 0;

    numsegs = bnd.getNumberOfSegments ();
    segs = new dof* [numsegs];

    int i = 0;
    typename std::list<ParticleType>::const_iterator it = bnd.begin (), end = bnd.end ();
    for ( ; it != end; ++it ) {

      typename ParticleType::ConstSegmentIteratorType sit = it->beginSegment (), send = it-> endSegment (), psit = send; --psit;
      for ( ; sit != send; psit = sit, ++sit, ++i ) {

        segs [i] = new bm::nodalCollNode<typename ParticleType::ConstSegmentType> ( *psit, *sit );
        segs [i]->idx = segs [i]->vidx = i;
      }
    }

    tree = new cluster2d_pca ( segs, 0, numsegs );

    numclus = 1;
    tree->subdivide ( segs, minclussize, numclus );
    blocktree = genblcltree ( tree, tree, clustereta, numblock );

    // Add one level on top for constraint
    if ( constraint ) {

      blcluster* consrow = new blcluster;
      blcluster* lagrcol = new blcluster;
      blcluster* brentry = new blcluster;
      blcluster* topnode = new blcluster;

      // Constraint row
      consrow->b1 = numsegs;
      consrow->b2 = 0;
      consrow->n1 = 1;
      consrow->n2 = numsegs;
      consrow->pos  = -1;
      consrow->ccl  = tree;
      consrow->idx  = numblock++;
      consrow->adm  = true;

      // Lagrange multiplier column
      lagrcol->b1 = 0;
      lagrcol->b2 = numsegs;
      lagrcol->n1 = numsegs;
      lagrcol->n2 = 1;
      lagrcol->pos  = 1;
      lagrcol->rcl  = tree;
      lagrcol->idx  = numblock++;
      lagrcol->adm  = true;

      // Bottom right entry
      brentry->b1 = numsegs;
      brentry->b2 = numsegs;
      brentry->n1 = 1;
      brentry->n2 = 1;
      brentry->pos  = 0;
      brentry->idx  = numblock++;
      brentry->adm  = true;

      // New top node
      topnode->b1 = 0;
      topnode->b2 = 0;
      topnode->n1 = numsegs + 1;
      topnode->n2 = numsegs + 1;
      topnode->pos = 0;
      topnode->adm = false;
      topnode->son1 = blocktree;
      topnode->son2 = lagrcol;
      topnode->son3 = consrow;
      topnode->son4 = brentry;

      blocktree = topnode;
    }
  }

  void increaseIndices ( blcluster* subtree, int off1, int off2, int idxoff, char pos = 0 ) {
    subtree->b1 += off1;
    subtree->b2 += off2;
    if ( pos ) subtree->pos = pos; // Will never need to set block on diagonal after quadruplication
  if ( subtree->isleaf () ) { subtree->idx += idxoff; return; }
    if ( subtree->son1 ) increaseIndices ( subtree->son1, off1, off2, idxoff, pos );
    if ( subtree->son2 ) increaseIndices ( subtree->son2, off1, off2, idxoff, pos );
    if ( subtree->son3 ) increaseIndices ( subtree->son3, off1, off2, idxoff, pos );
    if ( subtree->son4 ) increaseIndices ( subtree->son4, off1, off2, idxoff, pos );
  }

  void quadruplicate () {

    ++quadruplications;

    // Generate new top node
    blcluster* topnode = new blcluster;
    topnode->b1   = 0;
    topnode->b2   = 0;
    topnode->n1   = 2 * blocktree->n1;
    topnode->n2   = 2 * blocktree->n2;
    topnode->pos  = 0;
    topnode->adm  = false;

    // Make three copies of old tree and connect
    topnode->son1 = blocktree;
    topnode->son2 = new blcluster ( blocktree );
    topnode->son3 = new blcluster ( blocktree );
    topnode->son4 = new blcluster ( blocktree );
    blocktree     = topnode;

    // Adapt indices
    topnode->son2->pos = + 1;
    topnode->son3->pos = -1;
    increaseIndices ( blocktree->son2, blocktree->son1->n1, 0, numblock, 1 );
    increaseIndices ( blocktree->son3, 0, blocktree->son1->n2, 2*numblock, -1 );
    increaseIndices ( blocktree->son4, blocktree->son1->n1, blocktree->son1->n2, 3*numblock );
  }

  ~nodalClusterTree () {
    for ( int i = 0; i < numsegs; ++i ) delete segs [i];
    delete [] segs;
    delete tree;
    delete blocktree;
  }

dof** getDoFs () const { return segs; }
  int getNumDoFs () const { return numsegs; }
  blcluster* getBlockTree () const { return blocktree; }
  int getNumBlocks () const { return numblock; }
  void mapBtoH ( const aol::Vector<typename ParticleType::DataType>& bvec,
                 aol::Vector<typename ParticleType::DataType>& hvec ) const {
    int i = 0, j, qns = numsegs * 1 << quadruplications;
    for ( j = 0; j < qns; j += numsegs )
      for ( i = 0; i < numsegs; ++i )
        hvec [j + i] = bvec [j + getDoF ( i ) ];
    for ( i = j; i < bvec.size (); ++i ) hvec [i] = bvec [i];
  }
  void mapHtoB ( const aol::Vector<typename ParticleType::DataType>& hvec,
                 aol::Vector<typename ParticleType::DataType>& bvec ) const {
    int i = 0, j, qns = numsegs * 1 << quadruplications;
    for ( j = 0; j < qns; j += numsegs )
      for ( i = 0; i < numsegs; ++i )
        bvec [ j + getDoF ( i ) ] = hvec [j + i];
    for ( i = j; i < bvec.size (); ++i ) bvec [i] = hvec [i];
  }
  int getQuadruplications () const {
    return quadruplications;
  }
};

//! Cluster Tree for mixed BV problems and nodal basis functions
//! nIndicator[i] is supposed to be true iff DOF i is an unknown Neumann value
template <class ParticleType>
class nodalMixedClusterTree {
private:
  dof** segs;
  cluster2d_pca* dTree;
  cluster2d_pca* nTree;
  blcluster* blocktree;
  int numsegs, numdSegs, numnSegs;
  unsigned int numclus, numblock;
  unsigned int numblock_son[4];
  int quadruplications;

  int getDoF ( int i ) const { return segs [i] -> idx; }

public:
  nodalMixedClusterTree ( const bm::Boundary<ParticleType>& bnd,
                          const aol::BitVector& nIndicator,
                          int numNeumannSegs,
                          double clustereta = 1.1, int minclussize = 4 ) {

    quadruplications = 0;
    numsegs  = bnd.getNumberOfSegments();
    numnSegs = numNeumannSegs;
    numdSegs = numsegs - numnSegs;
    if ( (numnSegs == 0) || (numdSegs == 0) )
      throw ( aol::ParameterException ( "Mixed boundary values version: both Dirichlet and Neumann Data must not be empty!", __FILE__, __LINE__ ) );
    segs  = new dof* [numsegs];

    int i = 0, d = 0, n = 0;
    typename std::list<ParticleType>::const_iterator it = bnd.begin (), end = bnd.end ();
    for ( ; it != end; ++it ) {

      typename ParticleType::ConstSegmentIteratorType sit = it->beginSegment (), send = it-> endSegment (), psit = send; --psit;
      for ( ; sit != send; psit = sit, ++sit, i++ ) {

        if ( nIndicator[i] )  {
          segs [n] = new bm::nodalCollNode<typename ParticleType::ConstSegmentType> ( *psit, *sit );
          segs [n]->idx = segs [n]->vidx = i;
          n++;
        }
        else  {
          segs [numnSegs+d] = new bm::nodalCollNode<typename ParticleType::ConstSegmentType> ( *psit, *sit );
          segs [numnSegs+d]->idx = segs [numnSegs+d]->vidx = i;
          d++;
        }
      }
    }
    if ( n != numnSegs ) throw ( aol::ParameterException ( "Mismatch of BitVector and number of neumann segments!", __FILE__, __LINE__ ) );

//     for ( int i = 0; i < numsegs; i++ ) cerr << segs[i]->idx << " ";

    // Create one tree for Neumann DOFs and another for Dirichlet DOFs
    // Second argument is beginning DOF-Index in this cluser, thrid argument is beginning DOF-Index in next cluser
    nTree = new cluster2d_pca ( segs, 0,        numnSegs );
    dTree = new cluster2d_pca ( segs, numnSegs, numsegs  );

    numclus = 0;
    // Use global DOF-list! cluser2d_pca knows where his DOFs are!
    unsigned int addclus = 1; nTree->subdivide ( segs, minclussize, addclus ); numclus += addclus;
                 addclus = 1; dTree->subdivide ( segs, minclussize, addclus ); numclus += addclus;

    // Create 4 blocks
    numblock = 0;
    numblock_son[0] = 0; blcluster* ul = genblcltree ( nTree, nTree, clustereta, numblock_son[0] ); numblock += numblock_son[0];
    numblock_son[1] = 0; blcluster* ur = genblcltree ( nTree, dTree, clustereta, numblock_son[1] ); numblock += numblock_son[1];
    numblock_son[2] = 0; blcluster* ll = genblcltree ( dTree, nTree, clustereta, numblock_son[2] ); numblock += numblock_son[2];
    numblock_son[3] = 0; blcluster* lr = genblcltree ( dTree, dTree, clustereta, numblock_son[3] ); numblock += numblock_son[3];

//     cerr << endl << "Check this:" << endl;
//     cerr << "ul: " << ul->b1 << " " << ul->b2 << " " << ul->n1 << " "  << ul->n2 << " "  << ( ul->pos == 0 ? "diag" : "nondiag" ) << endl;
//     cerr << "ur: " << ur->b1 << " " << ur->b2 << " " << ur->n1 << " "  << ur->n2 << " "  << ( ur->pos == 0 ? "diag" : "nondiag" )  << endl;
//     cerr << "ll: " << ll->b1 << " " << ll->b2 << " " << ll->n1 << " "  << ll->n2 << " "  << ( ll->pos == 0 ? "diag" : "nondiag" )  << endl;
//     cerr << "lr: " << lr->b1 << " " << lr->b2 << " " << lr->n1 << " "  << lr->n2 << " "  << ( lr->pos == 0 ? "diag" : "nondiag" )  << endl;

    // New top node
    blcluster* top = new blcluster;
    top->b1 = 0;
    top->b2 = 0;
    top->n1 = numsegs;
    top->n2 = numsegs;
    top->pos = 0;
    top->adm = false;
    top->son1 = ul;
    top->son2 = ur;
    top->son3 = ll;
    top->son4 = lr;

    blocktree = top;

    // Modify non upper left blocks but don't set an offset here! genblctree is hell of intelligent!
    increaseIndices ( blocktree->son2, 0, 0, numblock_son[0],                                     1 );
    increaseIndices ( blocktree->son3, 0, 0, numblock_son[0] + numblock_son[1],                  -1 );
    increaseIndices ( blocktree->son4, 0, 0, numblock_son[0] + numblock_son[1] + numblock_son[2]    );

//     cerr << "MixedTree generated " << numclus << " clusters and " << numblock << " blocks" << endl;
  }

  void increaseIndices ( blcluster* subtree, int off1, int off2, int idxoff, char pos = 0 ) {
    subtree->b1 += off1;
    subtree->b2 += off2;
    if ( pos ) subtree->pos = pos; // An off-diagonal block will never move to the diagonal
    if ( subtree->isleaf () ) { subtree->idx += idxoff; return; }
    if ( subtree->son1 ) increaseIndices ( subtree->son1, off1, off2, idxoff, pos );
    if ( subtree->son2 ) increaseIndices ( subtree->son2, off1, off2, idxoff, pos );
    if ( subtree->son3 ) increaseIndices ( subtree->son3, off1, off2, idxoff, pos );
    if ( subtree->son4 ) increaseIndices ( subtree->son4, off1, off2, idxoff, pos );
  }

  void quadruplicate () {

    ++quadruplications;

    // Generate new top node
    blcluster* topnode = new blcluster;
    topnode->b1   = 0;
    topnode->b2   = 0;
    topnode->n1   = 2 * blocktree->n1;
    topnode->n2   = 2 * blocktree->n2;
    topnode->pos  = 0;
    topnode->adm  = false;

    // Make three copies of old tree and connect
    topnode->son1 = blocktree;
    topnode->son2 = new blcluster ( blocktree );
    topnode->son3 = new blcluster ( blocktree );
    topnode->son4 = new blcluster ( blocktree );
    blocktree     = topnode;

    // Adapt indices
    topnode->son2->pos = + 1; // Unnecessary?
    topnode->son3->pos = -1;
    increaseIndices ( blocktree->son2, blocktree->son1->n1, 0, numblock, 1 );
    increaseIndices ( blocktree->son3, 0, blocktree->son1->n2, 2*numblock, -1 );
    increaseIndices ( blocktree->son4, blocktree->son1->n1, blocktree->son1->n2, 3*numblock );
  }

  ~nodalMixedClusterTree () {
    for ( int i = 0; i < numsegs; ++i )  delete segs [i];
    delete [] segs;
    delete dTree;
    delete nTree;
    delete blocktree;
  }

  dof** getDoFs () const { return segs; }
  int getNumDoFs () const { return numsegs; }
  blcluster* getBlockTree () const { return blocktree; }
  int getNumBlocks () const { return numblock; }

  void mapBtoH ( const aol::Vector<typename ParticleType::DataType>& bvec,
                 aol::Vector<typename ParticleType::DataType>& hvec ) const {
    int i = 0, j, qns = numsegs * 1 << quadruplications;
    for ( j = 0; j < qns; j += numsegs )
      for ( i = 0; i < numsegs; ++i )
        hvec [j + i] = bvec [j + getDoF ( i ) ];
    for ( i = j; i < bvec.size (); ++i ) hvec [i] = bvec [i];
  }

  void mapHtoB ( const aol::Vector<typename ParticleType::DataType>& hvec,
                 aol::Vector<typename ParticleType::DataType>& bvec ) const {
    int i = 0, j, qns = numsegs * 1 << quadruplications;
    for ( j = 0; j < qns; j += numsegs )
      for ( i = 0; i < numsegs; ++i )
        bvec [ j + getDoF ( i ) ] = hvec [j + i];
    for ( i = j; i < bvec.size (); ++i ) bvec [i] = hvec [i];
  }

  int getQuadruplications () const {
    return quadruplications;
  }

  unsigned int getQuadrant ( unsigned int idx ) const {
    idx %= numblock;
    unsigned int barrier = numblock, ret = 4;
    do  {
      ret--;
      barrier -= numblock_son[ret];
    } while (idx < barrier);
    return ret+1;
  }
};

template <class _LocalOperatorType>
class centralMatrixEntry {
public:
  typedef _LocalOperatorType LocalOperatorType;
  typedef typename LocalOperatorType::ParticleType ParticleType;
  typedef typename LocalOperatorType::ParticleType::DataType DataType;
  typedef typename LocalOperatorType::ParticleType::ConstSegmentType SegmentType;
  typedef centralClusterTree<ParticleType> TreeType;
private:
  const LocalOperatorType op;
  DataType tau;
public:
  centralMatrixEntry ( const LocalOperatorType& o, DataType t = 0 ) : op ( o ), tau ( t ) {}
  void operator () ( unsigned int n1, dof** segs1, unsigned int n2, dof** segs2, DataType* res ) const {
    for ( unsigned int j = 0; j < n2; ++j ) {
      SegmentType seg2 = static_cast <centralCollNode <SegmentType>* > ( segs2 [j] ) -> getSegment ();
      for ( unsigned int i = 0; i < n1; ++i ) {
        SegmentType seg1 = static_cast <centralCollNode <SegmentType>* > ( segs1 [i] ) -> getSegment ();

  // *** Evil! Correct only for rectangle
        DataType add = 0;
        if ( seg1.sameParticle ( seg2 ) )
          if ( seg2.isFollowedBy ( seg1 ) || seg1.isFollowedBy ( seg2 ) ) {
            aol::Vec2<double> dir = seg1.getDirection ();
            add = 2 * tau / dir.normSqr ();
          }

        *res++ = op.evaluateLocally ( seg2, seg1 ) - add;
      }
    }
  }
};

//! Class for computing matrix entries when using linear nodal basis functions
//! @warning    _TreeType has to be based on nodal basis functions!
template <class _LocalOperatorType, class _TreeType = nodalClusterTree<typename _LocalOperatorType::ParticleType> >
class nodalMatrixEntry {
public:
  typedef _LocalOperatorType LocalOperatorType;
  typedef typename LocalOperatorType::ParticleType ParticleType;
  typedef typename LocalOperatorType::ParticleType::DataType DataType;
  typedef typename LocalOperatorType::ParticleType::ConstSegmentType SegmentType;
  typedef _TreeType TreeType;
private:
  const LocalOperatorType op;
  DataType tau;
  DataType fac;
public:
  nodalMatrixEntry ( const LocalOperatorType& o, DataType t = 0, DataType f = 1 ) : op ( o ), tau ( t ), fac ( f ) {}
  void operator () ( unsigned int n1, dof** segs1, unsigned int n2, dof** segs2, DataType* res ) const {
    for ( unsigned int j = 0; j < n2; ++j ) {
      SegmentType seg2a = static_cast <nodalCollNode <SegmentType>* > ( segs2 [j] ) -> getSegment ( 0 ),
                  seg2b = static_cast <nodalCollNode <SegmentType>* > ( segs2 [j] ) -> getSegment ( 1 );
      for ( unsigned int i = 0; i < n1; ++i ) {
        SegmentType seg1a = static_cast <nodalCollNode <SegmentType>* > ( segs1 [i] ) -> getSegment ( 0 ),
                    seg1b = static_cast <nodalCollNode <SegmentType>* > ( segs1 [i] ) -> getSegment ( 1 );

        DataType val = op.evaluateLocally ( seg2b, seg1b, false );
        val += op.evaluateLocally ( seg2a, seg1a, true );

        val *= fac;


        // *** Evil! Correct only for parametric
        if ( (tau != 0) && (seg1a.sameParticle ( seg2a ) ) ) {

          if ( seg1a == seg2a ) {
            // diagonal
            DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
            val += 2 * tau / ( hip * him );
          } else if ( seg1b == seg2a ) {
            // 1 before 2
            DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
            aol::Vec2<DataType> ni = cornerNormal ( seg1a, seg1b ), nip = cornerNormal ( seg2a, seg2b );
            val -= 2 * tau * ni * nip / ( hip * ( hip + him ) );
          } else if ( seg2b == seg1a ) {
            // 2 before 1
            DataType him = seg1a.getDirection () . norm (), hip = seg1b.getDirection () . norm ();
            aol::Vec2<DataType> ni = cornerNormal ( seg1a, seg1b ), nim = cornerNormal ( seg2a, seg2b );
            val -= 2 * tau * ni * nim / ( him * ( hip + him ) );
          }
        }
        *res++ = val;
      }
    }
  }
};

//! Trait for entry types
/** DoubleLayerEntryType specifies the corresponding EntryType when mixing single and double layer operator.
  * Currently only nodalMixedClusterTree versions are relevant, all other will fail to call getQuadrant
  */
template <class EntryType> class EntryTypeTrait {};

template <class ParticleType> class EntryTypeTrait<nodalMatrixEntry<LinElasticSolution<ParticleType>, nodalMixedClusterTree<ParticleType> > > {
public:
 typedef nodalMatrixEntry<LinElasticSolutionDiff<ParticleType>, nodalMixedClusterTree<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<nodalMatrixEntry<LinElasticSolutionDiff<ParticleType>, nodalMixedClusterTree<ParticleType> > > {
public:
 typedef nodalMatrixEntry<LinElasticSolutionDiff<ParticleType>, nodalMixedClusterTree<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<nodalMatrixEntry<LinElasticSolution<ParticleType> > > {
public:
 typedef nodalMatrixEntry<LinElasticSolutionDiff<ParticleType>, nodalClusterTree<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<nodalMatrixEntry<LinElasticSolutionDiff<ParticleType> > > {
public:
 typedef nodalMatrixEntry<LinElasticSolutionDiff<ParticleType>, nodalClusterTree<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<centralMatrixEntry<ElasticSolution<ParticleType> > > {
public:
 typedef centralMatrixEntry<ElasticSolutionDiff<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<centralMatrixEntry<ElasticSolutionDiff<ParticleType> > > {
public:
 typedef centralMatrixEntry<ElasticSolutionDiff<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<centralMatrixEntry<ElasticSolutionGrad<ParticleType> > > {
public:
 typedef centralMatrixEntry<ElasticSolutionGrad<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<centralMatrixEntry<OffCenterElasticSolutionGrad<ParticleType> > > {
public:
 typedef centralMatrixEntry<OffCenterElasticSolutionGrad<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<centralMatrixEntry<SingleLayerPotential<ParticleType> > > {
public:
 typedef centralMatrixEntry<SingleLayerPotential<ParticleType> > DoubleLayerEntryType;
};

template <class ParticleType> class EntryTypeTrait<nodalMatrixEntry<LinSingleLayerPotential<ParticleType> > > {
public:
 typedef nodalMatrixEntry<LinSingleLayerPotential<ParticleType> > DoubleLayerEntryType;
};


template <class DataType>
class HLUPreconditioner : public aol::Op< aol::Vector <DataType> > {
private:
  mblock<DataType>** lBlocks;
  mblock<DataType>** uBlocks;
  blcluster* blockcluster;

public:
  HLUPreconditioner ( const HMatrix<DataType>& mat, double delta = 0.1, int maxrank = 100 )
      : lBlocks ( NULL ), uBlocks ( NULL ), blockcluster ( NULL ) {
    genLUprecond ( mat.blclTree, mat.blcks, delta, maxrank, blockcluster, lBlocks, uBlocks, true );
  }
  ~HLUPreconditioner () {
    freembls ( blockcluster, lBlocks );
    freembls ( blockcluster, uBlocks );
    delete blockcluster;
  }

  void apply ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    dest = arg;
    HLU_solve ( blockcluster, lBlocks, uBlocks, & ( dest [0] ) );
  }
  void applyAdd ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    aol::Vector<DataType> temp ( dest );
    apply ( arg, dest );
    dest += temp;
  }
  void transposedApply ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    dest = arg;
    HLUT_solve ( blockcluster, lBlocks, uBlocks, & ( dest [0] ) );
  }
  void transposedApplyAdd ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    aol::Vector<DataType> temp ( dest );
    transposedApply ( arg, dest );
    dest += temp;
  }
};

template <class DataType>
class HDiagPreconditioner : public aol::DiagonalMatrix<DataType> {
private:
  mblock<DataType>** blocks;

public:
  HDiagPreconditioner ( const HMatrix<DataType>& mat )
      : aol::DiagonalMatrix<DataType>( mat.m, mat.n ), blocks ( mat.blcks ) {

#ifdef _OPENMP
    if ( mat.blclTree->isleaf() ) visitDiag( mat.blclTree );
    else  {
      if ( mat.blclTree->son1->isleaf() || mat.blclTree->son4->isleaf() )   {
#pragma omp parallel sections
        {
#pragma omp section
          { visitDiag( mat.blclTree->son1 ); }
#pragma omp section
          { visitDiag( mat.blclTree->son4 ); }
        }
      }
      else  {
        if ( mat.blclTree->son1->son1->isleaf() || mat.blclTree->son1->son4->isleaf() || mat.blclTree->son4->son1->isleaf() || mat.blclTree->son4->son4->isleaf() ) {
#pragma omp parallel sections
          {
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4 ); }
          }
        }
        else  {
#pragma omp parallel sections
          {
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4->son4 ); }
          }
        }
      }
    }
#else
    visitDiag( mat.blclTree );
#endif
  }

  void visitDiag( blcluster* tree )  {
    if ( tree->isleaf() )  {
      double * arg =  new double [tree->n1];
      double * dest = new double [tree->n1];
      int offset = tree->b1;
      for ( unsigned int i = 0; i < tree->n1; i++ )  {
        for ( unsigned int j = 0; j < tree->n1; j++ ) { arg[j]=0.; dest[j]=0.; }
        arg[i] = 1.;

        blocks [tree->idx] -> multa_mbl_vec ( 1, arg, dest );
        this->set( offset+i, offset+i, 1./dest[i] );
      }
      delete [] arg;
      delete [] dest;
    }
    else  {
      visitDiag( tree->son1 );
      visitDiag( tree->son4 );
    }
  }

  // apply, applyAdd provided by DiagonalMatrix
};

template <class DataType>
class HDiagBlockPreconditioner : public aol::SparseMatrix<DataType> {
private:
  mblock<DataType>** blocks;

public:
  HDiagBlockPreconditioner ( const HMatrix<DataType>& mat )
      : aol::SparseMatrix<DataType>( mat.m, mat.n ), blocks ( mat.blcks ) {

#ifdef _OPENMP
    if ( mat.blclTree->isleaf() ) visitDiag( mat.blclTree );
    else  {
      if ( mat.blclTree->son1->isleaf() || mat.blclTree->son4->isleaf() )   {
#pragma omp parallel sections
        {
#pragma omp section
          { visitDiag( mat.blclTree->son1 ); }
#pragma omp section
          { visitDiag( mat.blclTree->son4 ); }
        }
      }
      else  {
        if ( mat.blclTree->son1->son1->isleaf() || mat.blclTree->son1->son4->isleaf() || mat.blclTree->son4->son1->isleaf() || mat.blclTree->son4->son4->isleaf() ) {
#pragma omp parallel sections
          {
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4 ); }
          }
        }
        else  {
#pragma omp parallel sections
          {
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son1->son4->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son1->son4 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4->son1 ); }
#pragma omp section
            { visitDiag( mat.blclTree->son4->son4->son4 ); }
          }
        }
      }
    }
#else
    visitDiag( mat.blclTree );
#endif
  }

  void visitDiag( blcluster* tree )  {
    if ( tree->isleaf() )  {
      aol::FullMatrix<double> diagBlock    ( tree->n1, tree->n1 );
      aol::FullMatrix<double> inverseBlock ( tree->n1, tree->n1 );
      double * arg =  new double [tree->n1];
      double * dest = new double [tree->n1];
      for ( unsigned int j = 0; j < tree->n1; j++ ) { arg[j]=0.; dest[j]=0.; }
      int offset = tree->b1;
      for ( unsigned int i = 0; i < tree->n1; i++ )  {
        arg[i] = 1.;
        blocks [tree->idx] -> multa_mbl_vec ( 1, arg, dest );
        for ( unsigned int j = 0; j < tree->n1; j++ ) {
          diagBlock.ref( j, i ) = dest[j];
          arg[j]=0.;
          dest[j]=0.;
        }
      }
      inverseBlock.makeInverse( diagBlock );
      for ( unsigned int i = 0; i < tree->n1; i++ )  {
        for ( unsigned int j = 0; j < tree->n1; j++ )  {
          this->set( offset+i, offset+j, inverseBlock.ref( i, j ) );
        }
      }
      delete [] arg;
      delete [] dest;
    }
    else  {
      visitDiag( tree->son1 );
      visitDiag( tree->son4 );
    }
  }

  // apply, applyAdd provided by SparseMatrix
};

typedef enum {
  RESORT_VECTOR_NOT,
  RESORT_VECTOR_BY_SEGMENT_TREE
} ResortVectorType;

/** @todo ML:
 * Move blocktree extensions to appropriate subclasses of HMatrixOp ?
 * blocktree for vector operators and transmission systems
 * elasticity operators, subtract image of constant from diagonal
 *
 * Try new version of parallel ACA calls
 */

template <class EntryType>
class HMatrixOp : public aol::Op< aol::Vector <typename EntryType::DataType> > {

private:

#ifdef USE_PTHREADS
  friend void* HMatrixApproximationThread<EntryType> (void*);
  friend void* HMatrixApproximationThreadMulti<EntryType> (void*);
#endif

  typedef typename EntryType::DataType DataType;

  HMatrix <DataType> _matrix;
  const typename EntryType::TreeType& _tree;
  ResortVectorType _resort;

  void makeBlockList ( blcluster* node, blcluster**& blocklist ) {
    if ( node->isleaf() ) *blocklist++ = node;
    else {
      makeBlockList ( node->son1, blocklist );
      makeBlockList ( node->son2, blocklist );
      if ( node->son3 ) makeBlockList ( node->son3, blocklist );
      makeBlockList ( node->son4, blocklist );
    }
  }

  template<class FunctorType>
  void approximate ( const FunctorType& f, mblock<DataType>* block,
                     dof** segs1, dof** segs2, blcluster* blockcluster, double epsilon, unsigned maxrank ) {
    unsigned n1 = blockcluster->n1, n2 = blockcluster->n2;
    unsigned maxk = n1 * n2 / ( n1 + n2 ) + 1;
    if ( maxk > maxrank ) maxk = maxrank;

    bool succ = false;

    if ( blockcluster->isadm () ) {

      // Bottom right entry
      if ( blockcluster->rcl == NULL && blockcluster->ccl == NULL ) {
        block->setdense ();
        DataType* data = block->getdata ();
        data [0] = 0;
        return;
      }

      // Lagrange multiplier column
      if ( blockcluster->rcl != NULL && blockcluster->ccl == NULL ) {
        block->setdense ();
        DataType* data = block->getdata ();
        for ( unsigned int i = 0; i < n1; ++i ) data [i] = 1;
        return;
      }

      // Constraint row
      if ( blockcluster->rcl == NULL && blockcluster->ccl != NULL ) {
        block->setdense ();
        DataType* data = block->getdata ();
        for ( unsigned int i = 0; i < n2; ++i ) data [i] = sqrt ( segs2 [i]->getradius2 () );
        return;
      }

      // Normal block
      if ( maxk ) {
        succ = ACAn ( f, n1, segs1, n2, segs2, blockcluster->icom, epsilon, maxk, block );
      }
    }

    if ( !succ ) {
      block->setdense ();
      f ( n1, segs1, n2, segs2, block->getdata() );
    }
  }

  template<class FunctorType>
  void thread ( const FunctorType& f, mblock<DataType>*& matrix , blcluster* block,
                dof** segs1, dof** segs2, double epsilon, unsigned maxrank ) {

    matrix = new mblock<DataType> ( block->n1, block->n2 );
    QUOC_ASSERT ( matrix != NULL );
    approximate ( f, matrix, segs1, segs2, block, epsilon, maxrank );
  }

  void subtractFromDiagonal ( const aol::Vector<DataType>& diag, const blcluster* block, int off, int modulus ) {

    if ( ( static_cast<int>( block->b1 ) - static_cast<int>( block->b2 ) ) % modulus ) return; // non-diagonal with respect to modulus

    if ( block->isleaf () ) {
      mblock<DataType>* data = _matrix.blcks [ block->idx ];
      if ( ! data->isdense () )
        throw ( aol::ParameterException ( "HMatrixOp::subtractFromDiagonal, diagonal admissible block must be dense",
                                         __FILE__, __LINE__ ) );
      if ( block->n1 != block->n2 )
        throw ( aol::ParameterException ( "HMatrixOp::subtractFromDiagonal, diagonal block must be square",
                                         __FILE__, __LINE__ ) );

      int n = block->n1;
      for ( int i = 0; i < n; ++i ) {
        data->getdata () [ i * n + i ] -= diag [ off + i ];
      }
    } else {
      subtractFromDiagonal ( diag, block->son1, off,                   modulus );
      subtractFromDiagonal ( diag, block->son4, off + block->son1->n1, modulus );
    }
  }

public:

  HMatrixOp ( const EntryType& entries, const typename EntryType::TreeType& tree,
              double epsilon = 1E-4, int maxrank = 100 )
    : _matrix ( tree.getNumDoFs (), tree.getNumDoFs (), tree.getBlockTree () ), _tree ( tree ), _resort ( RESORT_VECTOR_NOT ) {

    _matrix.blcks = allocmbls ( tree.getBlockTree () );

    dof** segments = tree.getDoFs ();

    int numblock = tree.getNumBlocks ();
    blcluster** blocklist = new blcluster* [numblock];
    { blcluster** tempblocklist = blocklist; makeBlockList ( tree.getBlockTree (), tempblocklist ); }

#if defined ( USE_PTHREADS ) && ! defined ( _OPENMP )
    HMatrixApproximationParameter<EntryType> par1, par2;
    par1.numblock = numblock;
    par1.step = 2;
    par1.op = this;
    par1.blocklist = blocklist;
    par1.entries = &entries;
    par1.segments = segments;
    par1.epsilon = epsilon;
    par1.maxrank = maxrank;

    par2 = par1;

    par1.start = 0;
    par2.start = 1;

    pthread_t thr1, thr2;
    pthread_create (&thr1, NULL, HMatrixApproximationThread<EntryType>, reinterpret_cast<void*> (&par1));
    pthread_create (&thr2, NULL, HMatrixApproximationThread<EntryType>, reinterpret_cast<void*> (&par2));
    pthread_join (thr1, NULL);
    pthread_join (thr2, NULL);
#else
#ifdef _OPENMP
#pragma omp parallel for schedule ( guided )
#endif
    for ( int i = 0; i < numblock; i++ ) {

      blcluster* block = blocklist [i];
      mblock<DataType>*& matrixblock = _matrix.blcks [ block -> idx ];
      thread ( entries, matrixblock, block, segments + block -> b1, segments + block -> b2, epsilon, maxrank );
    }
#endif

    delete [] blocklist;

    if ( aol::debugging::matrix ) {
      std::ofstream os ( "matrix.ps" );
      psoutputH ( os, tree.getBlockTree (), tree.getNumDoFs (), _matrix.blcks, false );
      os.close();
    }
  }
  template <class ParameterType>
  HMatrixOp ( int dim, const ParameterType& parameter, const typename EntryType::TreeType& tree,
              double epsilon = 1E-4, int maxrank = 100,
        RenormalizeDiagonalType renorm = RENORMALIZE_DIAGONAL_NOT,
        ResortVectorType resort = RESORT_VECTOR_NOT,
        double factor = 1 )
    : _matrix ( 2 * tree.getNumDoFs (), 2 * tree.getNumDoFs (), tree.getBlockTree () ), _tree ( tree ), _resort ( resort ) {

    const EntryType*** entries = new const EntryType** [ dim ];
    for ( int i = 0; i < dim; ++i ) {
      entries [i] = new const EntryType* [ dim ];
      for ( int j = 0; j < dim; ++j ) {
  entries [i][j] = new const EntryType ( typename EntryType::LocalOperatorType ( parameter, j, i ), 0, factor );
      }
    }

    if ( tree.getQuadruplications () != 1 )
      throw ( aol::ParameterException ( "HMatrixOp, need quadruplicated tree", __FILE__, __LINE__ ) );

    _matrix.blcks = allocmbls ( tree.getBlockTree () );

    dof** segments = tree.getDoFs ();

    int numblock = tree.getNumBlocks ();
    int n = tree.getNumDoFs ();

    blcluster** blocklist = new blcluster* [4 * numblock];
  { blcluster** tempblocklist = blocklist; makeBlockList ( tree.getBlockTree (), tempblocklist ); }

#if defined ( USE_PTHREADS ) && ! defined ( _OPENMP )
    HMatrixApproximationParameterMulti<EntryType> par1, par2;
    par1.numblock = numblock;
    par1.n = n;
    par1.step = 2;
    par1.op = this;
    par1.blocklist = blocklist;
    par1.entries = entries;
    par1.segments = segments;
    par1.epsilon = epsilon;
    par1.maxrank = maxrank;

    par2 = par1;

    par1.start = 0;
    par2.start = 1;

    pthread_t thr1, thr2;
    pthread_create (&thr1, NULL, HMatrixApproximationThreadMulti<EntryType>, reinterpret_cast<void*> (&par1));
    pthread_create (&thr2, NULL, HMatrixApproximationThreadMulti<EntryType>, reinterpret_cast<void*> (&par2));
    pthread_join (thr1, NULL);
    pthread_join (thr2, NULL);
#else
#ifdef _OPENMP
#pragma omp parallel for schedule ( guided )
#endif
    for ( int i = 0; i < 4 * numblock; i++ ) {
      blcluster* block = blocklist [i];
      mblock<DataType>*& matrixblock = _matrix.blcks [ block -> idx ];
      int off1 = block -> b1, off2 = block -> b2;
      int o1 = off1 / n, o2 = off2 / n;
      off1 %= n; off2 %= n;
      thread ( * ( entries [o1] [o2] ), matrixblock, block, segments + off1, segments + off2, epsilon, maxrank );
    }
#endif

    delete [] blocklist;

    for ( int i = 0; i < dim; ++i ) {
      for ( int j = 0; j < dim; ++j ) {
  delete entries [i][j];
      }
      delete [] entries [i];
    }
    delete [] entries;

    if ( aol::debugging::matrix ) {
      std::ofstream os ( "matrix.ps" );
      psoutputH ( os, tree.getBlockTree (), 4 * numblock, _matrix.blcks, false );
      os.close();
    }

    if ( renorm == RENORMALIZE_DIAGONAL_RIGID_BODY_ARGUMENT ) {
      aol::Vector<DataType> one ( 2 * n ), temp ( 2 * n ), res1 ( 2 * n ), res2 ( 2 * n );
      for ( int i = 0; i < 2*n; ++i ) one [i] = ( i < n ? 1 : 0 );
      this->apply ( one, temp );
      if ( _resort == RESORT_VECTOR_BY_SEGMENT_TREE ) _tree.mapBtoH ( temp, res1 );
      else res1 = temp;
      if ( factor < 0 ) res1 += one;
      for ( int i = 0; i < 2*n; ++i ) one [i] = ( i < n ? 0 : 1 );
      this->apply ( one, temp );
      if ( _resort == RESORT_VECTOR_BY_SEGMENT_TREE ) _tree.mapBtoH ( temp, res2 );
      else res2 = temp;
      if ( factor < 0 ) res2 += one;
      subtractFromDiagonal ( res1, _matrix.blclTree->son1, 0, n );
      subtractFromDiagonal ( res1, _matrix.blclTree->son3, n, n );
      subtractFromDiagonal ( res2, _matrix.blclTree->son2, 0, n );
      subtractFromDiagonal ( res2, _matrix.blclTree->son4, n, n );
    }
  }

  ~HMatrixOp () {}

  //! Mixes a single layer operator with a corresponding double layer operator
  /** The HMatrixOp calling this function has to be initialized with the single layer operator.
   *  The argument has to be initialized with the corresponding double layer operator.
   *  This routine will swap entries that act on the searched Neumann values.
   *  The calling HMatrixOp will contain the right hand side operator and the argument will contain
   *  the left hand side operator afterwards.
   *  @warning  It is assumed that both operators are based on EXACTLY THE SAME nodalMixedCluserTree!
   */
  void mixWithDoubleLayerOp( HMatrixOp<typename EntryTypeTrait<EntryType>::DoubleLayerEntryType >& dop, bool print = false, bool suppressWarning = false )
  {
    if ( (!suppressWarning) && (&_tree != &dop.getTree()) )
      cerr << "*** WARNING ***  mixWithDoubleLayerOp: Trees are not stored at the same location!" << endl;

    HMatrix<DataType> & smat = getHMatrix();
    HMatrix<DataType> & dmat = dop.getHMatrix();

    int numblock = _tree.getNumBlocks ();
    blcluster** blocklist = new blcluster* [4 * numblock];
  { blcluster** tempblocklist = blocklist; makeBlockList ( _tree.getBlockTree (), tempblocklist ); } // don't know why but SIGSEGV otherwise

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < 4 * numblock; i++ ) {
      blcluster* block = blocklist [i];
      unsigned int quadrant = _tree.getQuadrant( block->idx );
      if ( (quadrant == 1) || (quadrant == 3) )  {  // swap
        int idx = block->idx;
        mblock<DataType>* temp = smat.blcks [ idx ];
        smat.blcks [ idx ] = dmat.blcks [ idx ];
        dmat.blcks [ idx ] = temp;
      }
    }

    delete [] blocklist;
    blocklist = NULL;

    if ( print || aol::debugging::matrix ) {
      std::ofstream os ( "RHSOp.ps" );
      psoutputH ( os, _tree.getBlockTree (), 4 * numblock, smat.blcks, false );
      os.close();
      std::ofstream os2 ( "LHSOp.ps" );
      psoutputH ( os2, _tree.getBlockTree (), 4 * numblock, dmat.blcks, false );
      os2.close();
    }
  }

  void applyAdd ( const aol::Vector <DataType>& arg,
                  aol::Vector <DataType>& dest ) const {

    //! @todo ML: Ahem!
    double* x1 = const_cast <DataType*> ( & ( arg [0] ) );
    double* y1 = & ( dest [0] );

    aol::Vector<DataType> tempx ( arg.size () );
    aol::Vector<DataType> tempy ( dest.size () );

    if ( _resort == RESORT_VECTOR_BY_SEGMENT_TREE ) {
      _tree.mapBtoH ( arg, tempx );
      x1 = const_cast <DataType*> ( & ( tempx [0] ) );
      y1 = const_cast <DataType*> ( & ( tempy [0] ) );
    }

#ifdef _OPENMP
    blcluster* tree = _matrix.blclTree;
    mblock<DataType>** blocks = _matrix.blcks;

    if ( tree->isleaf() )
      blocks [tree->idx] -> multa_mbl_vec ( 1, x1, y1 );
    else {

      double* x2 = x1 + tree->son1->n2;
      double* y2 = y1 + tree->son1->n1;

      if ( tree->son1->isleaf () || tree->son2->isleaf () || tree->son3->isleaf () || tree->son4->isleaf () ) {

#pragma omp parallel sections
        {
#pragma omp section
          {
            multa_H_vec ( 1, tree->son1, blocks, x1, y1 );
            multa_H_vec ( 1, tree->son2, blocks, x2, y1 );
          }
#pragma omp section
          {
            multa_H_vec ( 1, tree->son3, blocks, x1, y2 );
            multa_H_vec ( 1, tree->son4, blocks, x2, y2 );
          }
        }
      } else {

        double* x11 = x1;
        double* x12 = x1 + tree->son1->son1->n2;
        double* x21 = x2;
        double* x22 = x2 + tree->son2->son1->n2;

        double* y11 = y1;
        double* y12 = y1 + tree->son1->son1->n1;
        double* y21 = y2;
        double* y22 = y2 + tree->son3->son1->n1;

#pragma omp parallel sections
        {
#pragma omp section
          {
            multa_H_vec ( 1, tree->son1->son1, blocks, x11, y11 );
            multa_H_vec ( 1, tree->son1->son2, blocks, x12, y11 );
            multa_H_vec ( 1, tree->son2->son1, blocks, x21, y11 );
            multa_H_vec ( 1, tree->son2->son2, blocks, x22, y11 );
          }
#pragma omp section
          {
            multa_H_vec ( 1, tree->son1->son3, blocks, x11, y12 );
            multa_H_vec ( 1, tree->son1->son4, blocks, x12, y12 );
            multa_H_vec ( 1, tree->son2->son3, blocks, x21, y12 );
            multa_H_vec ( 1, tree->son2->son4, blocks, x22, y12 );
          }
#pragma omp section
          {
            multa_H_vec ( 1, tree->son3->son1, blocks, x11, y21 );
            multa_H_vec ( 1, tree->son3->son2, blocks, x12, y21 );
            multa_H_vec ( 1, tree->son4->son1, blocks, x21, y21 );
            multa_H_vec ( 1, tree->son4->son2, blocks, x22, y21 );
          }
#pragma omp section
          {
            multa_H_vec ( 1, tree->son3->son3, blocks, x11, y22 );
            multa_H_vec ( 1, tree->son3->son4, blocks, x12, y22 );
            multa_H_vec ( 1, tree->son4->son3, blocks, x21, y22 );
            multa_H_vec ( 1, tree->son4->son4, blocks, x22, y22 );
          }
        }
      }
    }
#else
#ifdef USE_PTHREADS
    blcluster* tree = _matrix.blclTree;
    mblock<DataType>** blocks = _matrix.blcks;

    if ( tree->isleaf() )
      blocks [tree->idx] -> multa_mbl_vec ( 1, x1, y1 );
    else {

      double* x2 = x1 + tree->son1->n2;
      double* y2 = y1 + tree->son1->n1;

      HMatrixVectorMultiplicationParameter<DataType> par1, par2;
      par1.blocks = par2.blocks = blocks;
      par1.arg1 = par2.arg1 = x1;
      par1.arg2 = par2.arg2 = x2;

      par1.tree1 = tree->son1;
      par1.tree2 = tree->son2;
      par1.dest = y1;
      par2.tree1 = tree->son3;
      par2.tree2 = tree->son4;
      par2.dest = y2;

      pthread_t thr1, thr2;
      pthread_create (&thr1, NULL, HMatrixVectorMultiplicationThread<DataType>, reinterpret_cast<void*> (&par1));
      pthread_create (&thr2, NULL, HMatrixVectorMultiplicationThread<DataType>, reinterpret_cast<void*> (&par2));
      pthread_join (thr1, NULL);
      pthread_join (thr2, NULL);
    }
#else
    _matrix.amux ( 1, x1, y1 );
#endif
#endif
    if ( _resort == RESORT_VECTOR_BY_SEGMENT_TREE ) {
      tempx.reallocate ( dest.size () );
      _tree.mapHtoB ( tempy, tempx );
      dest += tempx;
    }
  }

  const HMatrix <DataType>& getHMatrix () const { return _matrix; }
  ResortVectorType getResortVectorSetting () const { return _resort; }
  HMatrix <DataType>& getHMatrix () { return _matrix; }
  const typename EntryType::TreeType& getTree () const { return _tree; }
};

//! \warning untested!
template <class EntryType>
class HMatrixMVOp : public aol::Op< aol::MultiVector <typename EntryType::DataType> > {

public:

  typedef typename EntryType::DataType DataType;

private:

  HMatrix <DataType> _matrix;
  const typename EntryType::TreeType& _tree;
  ResortVectorType _resort;

  void multiplyAdd ( const aol::MultiVector <DataType>& arg, aol::MultiVector <DataType>& dest, int maxLev, int lev, int i, int j, const blcluster* subtree, const mblock<DataType>** blocks) const {

    if ( lev == maxLev ) {
      if ( _resort == RESORT_VECTOR_BY_SEGMENT_TREE ) {
  aol::Vector<DataType> tempx ( arg [i].size () );
  aol::Vector<DataType> tempy ( dest [j].size () );
  _tree.mapBtoH ( arg [i], tempx );
  multa_H_vec ( 1, subtree, blocks, tempx, tempy );
  tempx.reallocate ( tempy.size () );
  _tree.mapHtoB ( tempy, tempx );
  dest [j] += tempx;
      } else {
  multa_H_vec ( 1, subtree, blocks, arg [i], dest [j] );
      }
    } else {
      multiplyAdd (arg, dest, maxLev, lev + 1, 2*i + 0, 2*j + 0, subtree->son1, blocks);
      multiplyAdd (arg, dest, maxLev, lev + 1, 2*i + 1, 2*j + 0, subtree->son2, blocks);
      multiplyAdd (arg, dest, maxLev, lev + 1, 2*i + 0, 2*j + 1, subtree->son3, blocks);
      multiplyAdd (arg, dest, maxLev, lev + 1, 2*i + 1, 2*j + 1, subtree->son4, blocks);
    }
  }

public:

  // For this to work, the matrix must have a block structure compatible to the structure of the two multivectors
  HMatrixMVOp (const HMatrixOp<EntryType>& mat ) : _matrix (mat.getHMatrix ()), _tree (mat.getTree ()), _resort (mat.getResortVectorSetting ()) {}

  void applyAdd ( const aol::MultiVector <DataType>& arg,
                  aol::MultiVector <DataType>& dest ) const {

    int num = arg.numComponents (), lev = aol::ld ( num );

    QUOC_ASSERT ( num == dest.numComponents () && 1 << lev == num );

    multiplyAdd ( arg, dest, lev, 0, 0, 0, _matrix.blclTree );
  }
};


}

#endif // AHMED_INTERFACE_VERSION == 0
#endif // USE_EXTERNAL_AHMED
#endif // __AHMED_H
