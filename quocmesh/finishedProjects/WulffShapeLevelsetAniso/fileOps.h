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

#include <configurators.h>
#include <fstream>
#include <stdlib.h>

using namespace std;


/**************************************************************
 * some useful methods for managing files
 * @author Nemitz
 */



// ------------------------------------------------------------------
// Methode öffnet eine Datei und fügt die Variableninhalte hinten
// an, wenn die Datei noch nicht existiert wird sie angelegt.
// ------------------------------------------------------------------

void appendAreaAndVolume( char *filename, int timeStep, double value )
{
  // Datei wenn möglich zum anhängen öffnen
  ifstream fin;
  fin.open( filename );

  if ( fin.good() )      // d.h. Datei existiert bereits
  {
    ofstream fout( filename, ios::app );
    if ( fout.fail() )   // Datei konnte nicht zum erweitern geöffnet werden
      cerr << "ACHTUNG: KONNTE DATEI NICHT FÜR AUSGABE ÖFFNEN!";
    else
    {                    // Werte anhängen
      fout << timeStep << " " << value << endl;
      fout.close();
    }
  }
  else                   // d.h. eine neue Datei muss erzeugt werden
  {
    ofstream fout( filename );
    if ( fout.fail() )   // Datei konnte nicht zum erweitern geöffnet werden
      cerr << "ACHTUNG: KONNTE NEUE DATEI NICHT ERZEUGEN!";
    else
    {                    // Werte anhängen
      fout << timeStep << " " << value << endl;
      fout.close();
    }
  }
}
