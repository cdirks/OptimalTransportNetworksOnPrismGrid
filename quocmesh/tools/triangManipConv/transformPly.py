#!/usr/bin/python
# -*- coding: utf-8 -*-

# scales, shifts, and rotates all points in a .ply-file

from sys  import *
from os   import *
from math import *

if len( argv ) < 10:
  print r"scales, shifts, and rotates all points in a .ply (or ud-ply) file"
  print r"usage: "+argv[0]+" <input file> <output file> <scale factor> <x-shift> <y-shift> <z-shift> <x-rot> <y-rot> <z-rot>\n"

else:
  # get number of vertices
  inFile = file( argv[1], 'r' )
  inMesh = inFile.read()
  inFile.close()
  wordlist = inMesh.split()
  vertexNum = int(wordlist[wordlist.index("vertex")+1])

  # copy off header
  body = inMesh.split( "end_header" )
  outMesh = file( argv[2], 'w' )
  outMesh.write( body[0] )
  outMesh.write( "end_header\n" )

  # transform all points
  lines = body[1].split( "\n" )
  for i in range( 1, vertexNum + 1 ) :
    coords = lines[i].split()
    x = float(coords[0]) * float(argv[3]) + float(argv[4])
    y = float(coords[1]) * float(argv[3]) + float(argv[5])
    z = float(coords[2]) * float(argv[3]) + float(argv[6])
    # rotation about x-axis
    yold = y
    y = cos(float(argv[7])) * yold - sin(float(argv[7])) * z
    z = sin(float(argv[7])) * yold + cos(float(argv[7])) * z
    # rotation about y-axis
    zold = z
    z = cos(float(argv[8])) * zold - sin(float(argv[8])) * x
    x = sin(float(argv[8])) * zold + cos(float(argv[8])) * x
    # rotation about z-axis
    xold = x
    x = cos(float(argv[9])) * xold - sin(float(argv[9])) * y
    y = sin(float(argv[9])) * xold + cos(float(argv[9])) * y
    outMesh.write( str(x) + " " + str(y) + " " + str(z) )
    for j in range( 3, len(coords) ) :
      outMesh.write( " " + coords[j] )
    outMesh.write( "\n" )

  # copy polygons
  for i in range( vertexNum + 1, len( lines ) ) :
    outMesh.write( lines[i] + "\n" )
  outMesh.flush()
