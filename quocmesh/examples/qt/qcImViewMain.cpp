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

#include "qcImViewApp.hpp"

int main ( int argc, char *argv[] ) {
  // Set some basic information about our application. This influences how Qt stores our settings.
  QCoreApplication::setOrganizationName ( "AICES" );
  QCoreApplication::setOrganizationDomain ( "aices.rwth-aachen.de" );
  QCoreApplication::setApplicationName ( "Quoc Image Viewer" );

  // Check if a path is stored in the application settings. If so, add it to the search path.
  // Such a setting can be added manually with
  // settings.setValue( "path", "/usr/local/macports108/bin" );
  // and is stored persistently in the settings file.
  QSettings settings;
  QString path = settings.value( "path" ).toString();
  if ( path.isEmpty() == false )
    aol::appendSearchPath ( path.toLatin1().constData() );

  QuocViewerApp app ( argc, argv );
  return app.exec();
}
