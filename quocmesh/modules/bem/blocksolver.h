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

#ifndef __BLOCKSOLVER_H
#define __BLOCKSOLVER_H

#include <aol.h>
#include <vec.h>
#include <matrix.h>
#include <scalarArray.h>
#include <matrixInverse.h>

#include <bemesh.h>

namespace aol {

namespace {
template <class DataType, class IndexType>
struct Ref {
  IndexType col; DataType fac; DataType rhs;
  Ref () : col ( 0 ), fac ( 0 ), rhs ( 0 ) {}
  Ref ( IndexType c, DataType f, DataType r ) : col ( c ), fac ( f ), rhs ( r ) {}
  bool operator == ( const Ref& ref ) const {
    return col == ref.col && fac == ref.fac && rhs == ref.rhs;
  }
  bool operator != ( const Ref& ref ) const {
    return col != ref.col || fac != ref.fac || rhs != ref.rhs;
  }
};
}

//! This class gets a number of underdetermined linear equation systems
//! and additionally a number of lines with only up to two non-zero entries
//! so that the full system is quadratic
//! It stores references to matrices that will not be destroyed
template <class DataType = double, class IndexType = int>
class BlockSolver : public aol::Op<aol::Vector<DataType> > {

public:

  struct Line {
    struct Entry {
      IndexType block, var;
      DataType fac;
      Entry ( IndexType b, IndexType v, DataType f ) : block ( b ), var ( v ), fac ( f ) {};
    };
    Line ( const Entry& aa, const Entry& bb, DataType r ) : a ( aa ), b ( bb ), rhs ( r ) {};
    Entry a, b;
    DataType rhs;
  };

private:

  std::vector<const aol::FullMatrix<DataType>* > _mata;
  std::vector<const aol::FullMatrix<DataType>* > _matb;
  std::vector<const aol::Vector<DataType>* > _rhs;
  std::vector<bool> _delete;

  std::vector<Line> _line;

  mutable aol::Vector<DataType> _rhstmp;

  // Used to speed up search in modifyConstBlock
  IndexType _nextPositionInLineVector;
  // Caching of decomposed matrix
  mutable aol::QRInverse<DataType> * _QRinv;


public:

  BlockSolver ()
  : _nextPositionInLineVector(0), _QRinv( NULL )
  {}

  ~BlockSolver () {
    for ( int i = 0; i < static_cast<int> ( _delete.size () ); ++i )
      if ( _delete [i] ) {
        if ( _mata [i] ) delete _mata [i];
        if ( _matb [i] ) delete _matb [i];
        if ( _rhs [i] ) delete _rhs [i];
      }
    if ( _QRinv ) delete _QRinv;
  }

  void appendBlockRef ( const aol::FullMatrix<DataType>& mat, const aol::Vector<DataType>& rhs, bool del = false ) {
    _mata.push_back ( &mat );
    _matb.push_back ( NULL );
    _rhs.push_back ( &rhs );
    _delete.push_back ( del );
  };

  void appendBlockRef ( const aol::FullMatrix<DataType>& mata, const aol::FullMatrix<DataType>& matb, const aol::Vector<DataType>& rhs,
                        bool del = false ) {
    _mata.push_back ( &mata );
    _matb.push_back ( &matb );
    _rhs.push_back ( &rhs );
    _delete.push_back ( del );
  };

  void appendLine ( const Line& line ) {
    _line.push_back ( line );

    if ( aol::debugging::lowlevel )
      cerr << aol::mixedFormat ( line.a.fac ) << " *  x [ " << aol::intFormat ( line.a.var ) << " ] in block " << aol::intFormat ( line.a.block ) << "   +   "
     << aol::mixedFormat ( line.b.fac ) << " *  x [ " << aol::intFormat ( line.b.var ) << " ] in block " << aol::intFormat ( line.b.block ) << "   =   " << aol::mixedFormat ( line.rhs ) << endl;

    //if (line.a.block == line.b.block && line.a.fac != 0 && line.b.fac != 0)
    //throw aol::ParameterException ("aol::BlockSolver::appendLine,x crossreference within block forbidden", __FILE__, __LINE__);
  };

  void appendZeroBlock ( IndexType block, IndexType from, IndexType to ) {
    Line line ( typename Line::Entry ( block, from, 1 ), typename Line::Entry ( 0, 0, 0 ), 0 );

    for ( int i = from; i < to; ++i ) {

      line.a.var = i;
      appendLine ( line );
    }
  };

  void appendConstBlock ( IndexType block, IndexType from, IndexType to, DataType rhs ) {
    Line line ( typename Line::Entry ( block, from, 1 ), typename Line::Entry ( 0, 0, 0 ), rhs );

    for ( int i = from; i < to; ++i ) {

      line.a.var = i;
      appendLine ( line );
    }
  };

  //! Modifies an already set block.
  //! Searches for specific lines.
  void modifyConstBlock ( IndexType block, IndexType from, IndexType to, DataType rhs ) {
    for ( int i = from; i < to; ++i ) {
      // search for this line
      IndexType j = _nextPositionInLineVector;
      bool found = false;
      do  {
        found = ( _line[j].a.block == block ) && ( _line[j].a.var == i );
        j = (j+1) % _line.size();
      } while ( ( !found ) && ( j != _nextPositionInLineVector ) );

      if ( !found )  {
        stringstream s; s << "aol::BlockSolver::modifyConstBlock, line (block " << block << ", var " << i << ") does not exist";
        throw aol::ParameterException ( s.str() , __FILE__, __LINE__ );
      }
      else  {
        _line[ (j-1) % _line.size() ].rhs = rhs;
        _nextPositionInLineVector = j;
      }
    }
  };

  void appendVectorBlock ( IndexType block, IndexType from, IndexType to, const aol::Vector<DataType>& rhs ) {
    Line line ( typename Line::Entry ( block, from, 1 ), typename Line::Entry ( 0, 0, 0 ), 0 );

    for ( int i = from; i < to; ++i ) {

      line.a.var = i;
      line.rhs = rhs [i-from];
      appendLine ( line );
    }
  };

  template <class ConfigType>
  void appendFunctionBlock ( IndexType block, IndexType from, IndexType to,
                             const ConfigType& config, DataType ( *function ) ( const ConfigType&, IndexType ) ) {
    Line line ( typename Line::Entry ( block, from, 1 ), typename Line::Entry ( 0, 0, 0 ), 0 );

    for ( int i = from; i < to; ++i ) {

      line.a.var = i;
      line.rhs = ( *function ) ( config, i - from );

      appendLine ( line );
    }
  }

  void appendEqual ( IndexType blocka, IndexType linea, IndexType blockb, IndexType lineb, bool negative = false ) {
    Line line ( typename Line::Entry ( blocka, linea, 1 ), typename Line::Entry ( blockb, lineb, negative ? 1 : -1 ), 0 );
    appendLine ( line );
  }

  void appendEqualOffset ( IndexType blocka, IndexType linea, IndexType blockb, IndexType lineb, DataType offset, bool negative = false ) {
    Line line ( typename Line::Entry ( blocka, linea, 1 ), typename Line::Entry ( blockb, lineb, negative ? 1 : -1 ), offset );
    appendLine ( line );
  }

  void appendEqualBlock ( IndexType blocka, IndexType froma, IndexType toa, IndexType blockb, IndexType fromb, bool negative = false ) {
    Line line ( typename Line::Entry ( blocka, froma, 1 ), typename Line::Entry ( blockb, fromb, negative ? 1 : -1 ), 0 );

    for ( int a = froma, b = fromb; a < toa; ++a, ++b ) {

      line.a.var = a;
      line.b.var = b;
      appendLine ( line );
    }
  }

  void appendEqualBlockOffset ( IndexType blocka, IndexType froma, IndexType toa, IndexType blockb, IndexType fromb, DataType offset, bool negative = false ) {
    Line line ( typename Line::Entry ( blocka, froma, 1 ), typename Line::Entry ( blockb, fromb, negative ? 1 : -1 ), offset );

    for ( int a = froma, b = fromb; a < toa; ++a, ++b ) {

      line.a.var = a;
      line.b.var = b;
      appendLine ( line );
    }
  }

  void appendEqualBlockOffsetFactor ( IndexType blocka, IndexType froma, IndexType toa, IndexType blockb, IndexType fromb, DataType offset, DataType factor ) {
    Line line ( typename Line::Entry ( blocka, froma, factor ), typename Line::Entry ( blockb, fromb, -1 ), offset );

    for ( int a = froma, b = fromb; a < toa; ++a, ++b ) {

      line.a.var = a;
      line.b.var = b;
      appendLine ( line );
    }
  }

  void appendEqualBlockI ( IndexType blocka, IndexType froma, IndexType toa, IndexType blockb, IndexType fromb,
                           bool negative = false, DataType rhs = 0 ) {
    Line line ( typename Line::Entry ( blocka, froma, 1 ), typename Line::Entry ( blockb, fromb, negative ? 1 : -1 ), rhs );

    for ( int a = froma, b = fromb + toa - froma - 1; a < toa; ++a, --b ) {

      line.a.var = a;
      line.b.var = b;
      appendLine ( line );
    }
  }

  void appendEqualBlockILinear ( IndexType blocka, IndexType froma, IndexType toa, IndexType blockb, IndexType fromb,
                                 DataType fromval = 0., DataType toval = 0., bool negative = false ) {
    Line line ( typename Line::Entry ( blocka, froma, 1 ), typename Line::Entry ( blockb, fromb, negative ? 1 : -1 ), 0. );

    DataType valstep = (toval-fromval) / (toa-froma-1);
    DataType val = fromval;
    for ( int a = froma, b = fromb + toa - froma - 1; a < toa; ++a, --b, val+=valstep ) {
      line.a.var = a;
      line.b.var = b;
      line.rhs   = val;
      appendLine ( line );
    }
  }

  void appendPeriodicSquare ( IndexType block, IndexType start, IndexType sidelen,
                              bool negative = false, bool corner = false ) {
    int o = corner ? 1 : 0;

    appendEqualBlockI ( block, start + o, start + sidelen, block, start + 2 * sidelen + o, negative ); // T-B
    appendEqualBlockI ( block, start + sidelen + o, start + 2 * sidelen, block, start + 3 * sidelen + o, negative ); // L-R

    if ( corner ) {
      appendEqualBlockI ( block, start, start + 1, block, start + 2 * sidelen, negative ); // TR-BL
      appendEqualBlockI ( block, start + sidelen, start + sidelen + 1, block, start + 3 * sidelen, negative ); // TL-BR
    }
  }

  void appendPeriodicSquareJ ( IndexType block, IndexType start, IndexType sidelen,
                               bool negative = false ) {
    ++sidelen;
    // intermediate points
    appendEqualBlockI ( block, start + 1, start + sidelen - 1, block, start + 2 * sidelen + 1, negative ); // T-B
    appendEqualBlockI ( block, start + sidelen + 1, start + 2 * sidelen - 1, block, start + 3 * sidelen + 1, negative ); // L-R
    // endpoints of top and bottom line
    appendEqualBlockI ( block, start, start + 1, block, start + sidelen - 1, !negative ); // BR-BL
    appendEqualBlockI ( block, start + sidelen - 1, start + sidelen, block, start + 2 * sidelen, negative ); // TR-BR
    appendEqualBlockI ( block, start + 2 * sidelen, start + 2 * sidelen + 1, block, start + 3 * sidelen - 1, !negative ); // TL-TR
    // endpoints of left and right line
    appendEqualBlockI ( block, start + sidelen, start + sidelen + 1, block, start + 2 * sidelen - 1, !negative ); // BR-TR
    appendEqualBlockI ( block, start + 2 * sidelen - 1, start + 2 * sidelen, block, start + 3 * sidelen, negative ); // TR-TL
    appendEqualBlockI ( block, start + 3 * sidelen, start + 3 * sidelen + 1, block, start + 4 * sidelen - 1, !negative ); // TL-TR
  }

  void appendPeriodicShiftedSquareX ( IndexType block, IndexType start, IndexType sidelen, DataType shift,
                                      bool negative = false, bool corner = false ) {
    int o = corner ? 1 : 0;

    appendEqualBlockI ( block, start + o, start + sidelen, block, start + 2 * sidelen + o, negative, 0 ); // T-B
    appendEqualBlockI ( block, start + sidelen + o, start + 2 * sidelen, block, start + 3 * sidelen + o, negative, - shift ); // L-R

    if ( corner ) {
      appendEqualBlockI ( block, start, start + 1, block, start + sidelen, negative, shift ); // BR-BL
      appendEqualBlockI ( block, start + sidelen, start + sidelen + 1, block, start + 2 * sidelen, negative, 0 ); // TR-BR
      appendEqualBlockI ( block, start + 2 * sidelen, start + 2 * sidelen + 1, block, start + 3 * sidelen, negative, - shift ); // TL-TR
    }
  }

  void appendPeriodicShiftedSquareY ( IndexType block, IndexType start, IndexType sidelen, DataType shift,
                                      bool negative = false, bool corner = false ) {
    int o = corner ? 1 : 0;

    appendEqualBlockI ( block, start + o, start + sidelen, block, start + 2 * sidelen + o, negative, shift ); // T-B
    appendEqualBlockI ( block, start + sidelen + o, start + 2 * sidelen, block, start + 3 * sidelen + o, negative, 0 ); // L-R

    if ( corner ) {
      appendEqualBlockI ( block, start, start + 1, block, start + sidelen, negative, 0 ); // BR-BL
      appendEqualBlockI ( block, start + sidelen, start + sidelen + 1, block, start + 2 * sidelen, negative, shift ); // TR-BR
      appendEqualBlockI ( block, start + 2 * sidelen, start + 2 * sidelen + 1, block, start + 3 * sidelen, negative, 0 ); // TL-TR
    }
  }

  //! Append boundary conditions for an affine periodic square
  /**
   *  @note    negative values in tilt correspond to pulling loads (BG: counter-intuitive?)
   *  @warning when supplying positive values make sure your domain is large enough to avoid intersections!
   */
  void appendPeriodicTiltedSquare ( IndexType block, IndexType start, IndexType sidelen, aol::Vec2<DataType> tilt,
                                    bool negative = false, bool corner = false ) {
    int o = corner ? 1 : 0;

    appendEqualBlockI ( block, start + o, start + sidelen, block, start + 2 * sidelen + o, negative, tilt [1] );                 // bottom - top = tilt[1]
    appendEqualBlockI ( block, start + sidelen + o, start + 2 * sidelen, block, start + 3 * sidelen + o, negative, - tilt [0] ); // right - left = -tilt[0]

    if ( corner ) {
      appendEqualBlockI ( block, start, start + 1, block, start + sidelen, negative, tilt [0] ); // BR-BL
      appendEqualBlockI ( block, start + sidelen, start + sidelen + 1, block, start + 2 * sidelen, negative, tilt [1] ); // TR-BR
      appendEqualBlockI ( block, start + 2 * sidelen, start + 2 * sidelen + 1, block, start + 3 * sidelen, negative, - tilt [0] ); // TL-TR
    }
  }

  void apply_fast ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest, bool iterRefine = false ) const {
    // Faster, solves _line by column operations

    IndexType w = 0, h = 0;
    std::vector<IndexType> blockstart; blockstart.push_back ( 0 );
    for ( typename std::vector<const aol::FullMatrix<DataType>* >::const_iterator ait = _mata.begin (), bit = _matb.begin ();
          ait != _mata.end (); ++ait, ++bit ) {
      IndexType blockw = ( *ait )->getNumCols () + ( ( *bit ) ? ( *bit )->getNumCols () : 0 );
      IndexType blockh = ( *ait )->getNumRows ();
      w += blockw;
      h += blockh;
      blockstart.push_back ( w );
    }
    blockstart.pop_back ();

    IndexType fullw = w, fullh = h;

    fullh += _line.size ();
    w -= _line.size ();

    // Store crossreferences in vector
    // ref [i] = {j,f,r} means col [i] = f * col [j] + r
    std::vector<Ref<DataType, IndexType> > ref ( fullw );
    for ( int i = 0; i < fullw; ++i ) ref [i] = Ref<DataType, IndexType> ( i, 1, 0 );

    for ( typename std::vector<Line>::const_iterator lit = _line.begin (); lit != _line.end (); ++lit ) {

      // Set column to explicit value
      if ( lit->a.fac != 0 && lit->b.fac == 0 ) {
        IndexType i = blockstart [lit->a.block] + lit->a.var;
        Ref<DataType, IndexType> self ( i, 1, 0 );
        if ( ref [i] != self )
          throw aol::ParameterException ( "aol::BlockSolver::apply, cannot set value, already related", __FILE__, __LINE__ );

        ref [i] = Ref<DataType, IndexType> ( 0, 0, lit->rhs / lit->a.fac );
      }

      // Corellate columns
      if ( lit->a.fac != 0 && lit->b.fac != 0 ) {
        IndexType i = blockstart [lit->a.block] + lit->a.var, j = blockstart [lit->b.block] + lit->b.var;

        if ( aol::debugging::lowlevel ) {
          cerr << lit->b.block << " (" << lit->b.var << ")    --    " << lit->a.block << " (" << lit->a.var << ")    --    " << endl;
          cerr << "col [" << j << "] = " <<  - lit->a.fac / lit->b.fac << " * col [" << i << "] + " << lit->rhs / lit->b.fac << endl;
        }

        Ref<DataType, IndexType> self ( j, 1, 0 );
        if ( ref [j] != self )
          throw aol::ParameterException ( "aol::BlockSolver::apply, cannot relate, already related", __FILE__, __LINE__ );
        if ( i > j )
          throw aol::ParameterException ( "aol::BlockSolver::apply, wrong relation direction", __FILE__, __LINE__ );

        ref [j] = Ref<DataType, IndexType> ( i, - lit->a.fac / lit->b.fac, lit->rhs / lit->b.fac );
      }
    }

    // Chain down
    for ( IndexType i = 0; i < fullw; ++i ) {
      IndexType col = ref [i].col;
      DataType fac = ref [i].fac, rhs = ref [i].rhs, oldfac = fac;
      while ( ref [col].col != col ) {
        rhs += oldfac * ref [col].rhs;
        fac *= ref [col].fac;
        oldfac = ref [col].fac;
        col = ref [col].col;
      }
      ref [i].col = col;
      ref [i].fac = fac;
      ref [i].rhs = rhs;
    }

    // Shift together
    aol::Vector<IndexType> unused ( fullw );
    unused.setZero ();

    unused [0] = ( ref [0].col == 0 && ref [0].fac == 1 ) ? 0 : 1;
    for ( int i = 1; i < fullw; ++i ) {
      if ( ref [i].col == i && ref [i].fac == 1 )
        unused [i] = unused [i-1];
      else
        unused [i] = unused [i-1] + 1;
    }

    for ( int i = 1; i < fullw; ++i )
      ref [i].col -= unused [ref [i].col];

    // If DOF i is assigned a fixed value fac will be 0 and this dof will not exist in the condensed system.
    // Therefore col is irrelevant and can also be invalid; we correct this here.
    for ( int i = 1; i < fullw; ++i )
      if ( ref [i].fac == 0 ) ref [i].col = 0;

    // Build matrix
    aol::FullMatrix<DataType> mat ( h, w );
    aol::Vector<DataType> rhs ( arg );

    // Column operations
    IndexType r = 0, c = 0, i = 0;
    for ( typename std::vector<const aol::FullMatrix<DataType>* >::const_iterator ait = _mata.begin (), bit = _matb.begin ();
          ait != _mata.end (); ++ait, ++bit ) {
      for ( i = 0; i < ( *ait )->getNumCols (); ++i ) {
        if ( ref [c+i].fac == 0 ) {
          aol::Vector<DataType> col ( ( *ait )->getNumRows () );
          ( *ait )->getSubColumn ( 0, i, col );
          col *= -ref [c+i].rhs;
          rhs.addBlock ( r, col );
        } else {
          IndexType j = ref [c+i].col;
          aol::Vector<DataType> col ( ( *ait )->getNumRows () );
          ( *ait )->getSubColumn ( 0, i, col );
          aol::Vector<DataType> col2 ( col );
          col *= -ref [c+i].rhs;
          rhs.addBlock ( r, col );
          col2 *= ref [c+i].fac;
          mat.addSubColumn ( r, j, col2 );
        }
      }
      c += i;
      if ( *bit ) {
        for ( i = 0; i < ( *bit )->getNumCols (); ++i ) {
          if ( ref [c+i].fac == 0 ) {
            aol::Vector<DataType> col ( ( *bit )->getNumRows () );
            ( *bit )->getSubColumn ( 0, i, col );
            col *= -ref [c+i].rhs;
            rhs.addBlock ( r, col );
          } else {
            IndexType j = ref [c+i].col;
            aol::Vector<DataType> col ( ( *bit )->getNumRows () );
            ( *bit )->getSubColumn ( 0, i, col );
            aol::Vector<DataType> col2 ( col );
            col *= -ref [c+i].rhs;
            rhs.addBlock ( r, col );
            col2 *= ref [c+i].fac;
            mat.addSubColumn ( r, j, col2 );
          }
        }
        c += i;
      }
      r += ( *ait )->getNumRows ();
    }

    if ( aol::debugging::matrix ) {
      cerr << mat;

      qc::ScalarArray<DataType, qc::QC_2D> arr;
      arr.visualizeMatrix( mat );
      arr.setOverflowHandling ( aol::CLIP_THEN_SCALE, -0.2, 0.4 );
      arr.save ( "mat.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );
    }

    // Solve
    aol::Vector<DataType> tempdest ( w );
    if ( !_QRinv ) _QRinv = new aol::QRInverse<DataType> ( mat );

    _QRinv->apply ( rhs, tempdest );

    // Iterative refinement
    if ( iterRefine )
    {
      int iter = 0;
      aol::Vector<DataType> res ( h );
      do {
        mat.apply( tempdest, res );
        res.addMultiple( rhs, -1. );
        res*=-1.;
        _QRinv->applyAdd ( res, tempdest );
        iter++;
      } while ( (res.norm() > 1e-15) && (iter < 10) );
    }

    // Restore full vector
    for ( int i = 0; i < fullw; ++i ) dest [i] = ref [i].fac * tempdest [ref [i].col] + ref [i].rhs;
  }

  void apply ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    IndexType w = 0, h = 0;
    std::vector<IndexType> blockstart; blockstart.push_back ( 0 );
    for ( typename std::vector<const aol::FullMatrix<DataType>* >::const_iterator ait = _mata.begin (), bit = _matb.begin ();
          ait != _mata.end (); ++ait, ++bit ) {
      IndexType blockw = ( *ait )->getNumCols () + ( ( *bit ) ? ( *bit )->getNumCols () : 0 );
      IndexType blockh = ( *ait )->getNumRows ();
      w += blockw;
      h += blockh;
      blockstart.push_back ( w );
    }
    blockstart.pop_back ();

    h += _line.size ();

    aol::FullMatrix<DataType> fullmat ( h, w );

    IndexType r = 0, c = 0;
    for ( typename std::vector<const aol::FullMatrix<DataType>* >::const_iterator ait = _mata.begin (), bit = _matb.begin ();
          ait != _mata.end (); ++ait, ++bit ) {
      fullmat.setBlock ( r, c, **ait );
      c += ( *ait )->getNumCols ();
      if ( *bit ) {
        fullmat.setBlock ( r, c, **bit );
        c += ( *bit )->getNumCols ();
      }
      r += ( *ait )->getNumRows ();
    }

    for ( typename std::vector<Line>::const_iterator lit = _line.begin (); lit != _line.end (); ++lit, ++r ) {

      if ( lit->a.fac ) fullmat.set ( r, blockstart [lit->a.block] + lit->a.var, lit->a.fac );
      if ( lit->b.fac ) fullmat.set ( r, blockstart [lit->b.block] + lit->b.var, lit->b.fac );
    }

    if ( aol::debugging::matrix ) {
      cerr << fullmat;

      qc::ScalarArray<DataType, qc::QC_2D> arr;
      arr.visualizeMatrix( fullmat );
      arr.setOverflowHandling ( aol::CLIP_THEN_SCALE, -0.2, 0.4 );
      arr.save ( "fullmat.pgm", qc::PGM_UNSIGNED_CHAR_BINARY );
    }

    aol::QRInverse<DataType> fullinv ( fullmat );
    fullinv.apply ( arg, dest );

    if ( aol::debugging::matrix ) {
      aol::Vector<DataType> res ( arg );
      fullmat.apply ( dest, res );
      res -= arg;
      cerr << "Residual : " << aol::mixedFormat ( res.norm () ) << endl;
    }
  };

  const aol::Vector<DataType>& getRHS_fast () const {
    int i = 0;
    for ( typename std::vector<const aol::Vector<DataType>* >::const_iterator it = _rhs.begin (); it != _rhs.end (); ++it ) {
      i += ( *it )->size ();
    }
    _rhstmp.resize ( i );

    i = 0;
    for ( typename std::vector<const aol::Vector<DataType>* >::const_iterator it = _rhs.begin (); it != _rhs.end ();
          i += ( *it )->size (), ++it ) {
      _rhstmp.setBlock ( i, **it );
    }

    return _rhstmp;
  };

  const aol::Vector<DataType>& getRHS () const {
    int i = 0;
    for ( typename std::vector<const aol::Vector<DataType>* >::const_iterator it = _rhs.begin (); it != _rhs.end (); ++it ) {
      i += ( *it )->size ();
    }
    _rhstmp.resize ( i + _line.size () );

    i = 0;
    for ( typename std::vector<const aol::Vector<DataType>* >::const_iterator it = _rhs.begin (); it != _rhs.end ();
          i += ( *it )->size (), ++it ) {
      _rhstmp.setBlock ( i, **it );
    }

    for ( typename std::vector<Line>::const_iterator lit = _line.begin (); lit != _line.end (); ++lit, ++i ) {

      _rhstmp [i] = lit->rhs;
    }

    return _rhstmp;
  };

  void applyAdd ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {
    aol::Vector<DataType> temp ( dest );
    apply ( arg, dest );
    dest += temp;
  };
};
}

#endif
