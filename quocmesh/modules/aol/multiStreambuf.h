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

#ifndef __MULTISTREAMBUF_H
#define __MULTISTREAMBUF_H

#include <aol.h>

namespace aol {

/** \brief Output stream buffer that sends its
 *  content to multiple other stream buffers.
 *
 *  \author von Deylen
 */
class MultiStreambuf : public streambuf {
public:
  MultiStreambuf ();

  size_t addStreambuf    ( streambuf * buffer );
  size_t removeStreambuf ( streambuf * buffer );

protected:
  typedef std::char_traits<char> traits_type;
  typedef traits_type::int_type  int_type;

  typedef list<streambuf*> StreambufList;
  StreambufList _streambufPtrList;

  int_type overflow ( int_type c );
  int_type sync ();
};

/**
 *  \brief Print cout and cerr not only to console
 *  (or whereever you have redirected it), but
 *  also to a file.
 *
 *  When you create an instance of this class, the
 *  standard cout / cerr buffers are replaced by
 *  "output buffer collections" that distribute every
 *  write attemp to the previously used cout / cerr buffer
 *  as well as to a file output buffer.
 *
 *  It is possible to iterate this procedure by
 *  creating multiple instances of AdditionalOutputToFile
 *  (but they have to be destroyed in reverse order).
 *
 *  \author von Deylen
 */
class AdditionalOutputToFile {
public:
  AdditionalOutputToFile ( string Filename );
  ~AdditionalOutputToFile ();

protected:
  string      _filename;
  ofstream    _filestream;
  streambuf * _previousCoutStreambuf;
  streambuf * _previousCerrStreambuf;
  streambuf * _previousClogStreambuf;
  MultiStreambuf _multiCoutStreambuf;
  MultiStreambuf _multiCerrStreambuf;
  MultiStreambuf _multiClogStreambuf;
};

} // end of namespace aol.

#endif
