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

#ifndef __BZIPIOSTREAM_H
#define __BZIPIOSTREAM_H

#include "aol.h"

#ifdef USE_LIB_BZ2
#include <bzlib.h>
// Unfortunately bzlib.h defines a lot of crap under some platforms.
// We need to kill the defines here.
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef PACKED
#undef PACKED
#endif
#ifdef rad1
#undef rad1
#endif
#ifdef rad2
#undef rad2
#endif
#ifdef rad3
#undef rad3
#endif
#endif

namespace aol {

// VC++ 2012 always produces warning C4250 when deriving from the std stream classes, but this warning should be ignored.
#if defined ( _MSC_VER ) && ( _MSC_VER == 1700 )
#pragma warning ( push )
#pragma warning ( disable : 4250 )
#endif

/**
 * \brief Can be used like an ifstream, but supports on the fly bzip2 decompression using libbz2.
 *
 * Automatically decides whether to decompress the input file or not based on the suffix of
 * the constructor argument FileName.
 *
 * \author Berkels
 */
class Bzipifstream : public stringstream{
public:
  Bzipifstream ( const char *FileName ){
    FILE *tempf;
    char buf[200000];
    tempf = fopen ( FileName, "rb" );
    if ( !tempf ){
      string errorMessage = "Could not open \"";
      errorMessage += FileName;
      errorMessage += "\" for reading.\n";
      throw Exception ( errorMessage.c_str(), __FILE__, __LINE__ );
    }

    if ( hasBzipSuffx ( FileName ) ) {
#ifdef USE_LIB_BZ2
      BZFILE *tempb;
      int bzerror;

      tempb = BZ2_bzReadOpen ( &bzerror, tempf, 0, 0, NULL, 0 );
      while ( bzerror == BZ_OK ) {
        const int bytesRead = BZ2_bzRead ( &bzerror, tempb, buf, sizeof ( buf ) );
        (*this).write ( buf, bytesRead );
      }
      BZ2_bzReadClose ( &bzerror, tempb );
#else
      throw Exception ( "Reading bz2 compressed files with Bzipifstream without using bzlib is impossible.\nDefine USE_LIB_BZ2 and link bzlib, e.g by using CFLAG += -DUSE_LIB_BZ2 and LFLAGS += -lbz2, to remedy this.\n", __FILE__, __LINE__ );
#endif // USE_LIB_BZ2
    } else {
      while ( !feof(tempf) ) {
        const int bytesRead = static_cast<int>(fread ( buf, 1, sizeof ( buf ), tempf ));
        (*this).write ( buf, bytesRead );
      }
    }
    fclose ( tempf );
  }

private:
  template< typename AnyThing >
  Bzipifstream& operator<< ( const AnyThing& ); // do not implement
};

/**
 * \brief Can be used like an ofstream, but supports on the fly bzip2 compression using libbz2.
 *
 * \author Berkels
 */
class Bzipofstream : public stringstream{
  FILE *outFile;
  const bool compressOutput;
public:
  Bzipofstream ( const char *FileName )
    : outFile ( NULL ),
      compressOutput ( hasBzipSuffx ( FileName ) )
  {
    outFile = fopen ( FileName, "wb" );
    if ( !outFile ){
      string errorMessage = "Could not open \"";
      errorMessage += FileName;
      errorMessage += "\" for output.\n";
      throw Exception ( errorMessage.c_str(), __FILE__, __LINE__ );
    }
  }
  ~Bzipofstream(){
    close();
  }
  void close() {
    char buf[200000];
    // We can't do anything when there is no outfile (probably close was already called).
    if ( outFile == NULL )
      return;

    if ( compressOutput ){
#ifdef USE_LIB_BZ2
      BZFILE *b;
      int bzerror;
      b = BZ2_bzWriteOpen ( &bzerror, outFile, 9, 0, 0 );
      while ( !this->eof() ) {
        (*this).read ( buf, sizeof ( buf ) );
        const int bytesRead = static_cast<int>(this->gcount());
        BZ2_bzWrite ( &bzerror, b, buf, bytesRead );
      }
      unsigned int temp;
      BZ2_bzWriteClose ( &bzerror, b, 0, &temp, &temp );
#else
      throw Exception ( "Writing bz2 compressed files with Bzipofstream without using bzlib is impossible.\nDefine USE_LIB_BZ2 and link bzlib, e.g by using CFLAG += -DUSE_LIB_BZ2 and LFLAGS += -lbz2, to remedy this.\n", __FILE__, __LINE__ );
#endif
    }
    else{
      while ( !this->eof() ) {
        (*this).read ( buf, sizeof ( buf ) );
        const int bytesRead = static_cast<int>(this->gcount());
        fwrite( buf, 1, bytesRead, outFile );
      }
    }
    fclose( outFile );
    outFile = NULL;
  }

private:
  template< typename AnyThing >
  Bzipofstream& operator>> ( const AnyThing& ); // do not implement
};

// Turn on the warning C4250 again.
#if defined ( _MSC_VER ) && ( _MSC_VER == 1700 )
#pragma warning ( pop )
#endif

}// end namespace aol

#endif //__BZIPIOSTREAM_H

