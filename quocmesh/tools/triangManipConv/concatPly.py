#!/usr/bin/python
# -*- coding: utf-8 -*-

# concatenates multiple .ply-files into a single one

from sys import *
from os  import *

if len( argv ) < 3:
  print r"concatenates multiple .ply-files into a single one"
  print r"usage: "+argv[0]+" <input file list> <output file>\n"

else:
  # count vertices and faces
  vertexNums = []
  faceNums = []
  vertexCount = 0
  faceCount = 0
  for i in range( 1, len( argv ) - 1 ):
    inFile = file( argv[i], 'r' )
    inMesh = inFile.read()
    inFile.close()
    wordlist = inMesh.split()
    vertexNums.append( int(wordlist[wordlist.index("vertex")+1]) )
    faceNums.append( int(wordlist[wordlist.index("face")+1]) )
    vertexCount = vertexCount + vertexNums[i-1]
    faceCount = faceCount + faceNums[i-1]

  # write .ply-header
  inFile = file( argv[1], 'r' )
  inHeader = inFile.read().split( "end_header" )
  inFile.close()
  outMesh = file( argv[len( argv )-1], 'w' )
  outMesh.write ( "ply\n" )
  outMesh.write ( "format ascii 1.0\n" )
  outMesh.write ( "comment: concatenated from" )
  for i in range( 1, len( argv ) - 1 ):
    outMesh.write( " "+argv[i] )
  outMesh.write( "\n" )
  outMesh.write ( "element vertex "+str(vertexCount)+"\n" )
  props = (inHeader[0].split( "element" )[1]).split( "property" )
  for i in range( 1, len( props ) ):
    outMesh.write ( "property" + props[i] )
  outMesh.write ( "element face "+str(faceCount)+"\n" )
  props = (inHeader[0].split( "element" )[2]).split( "property" )
  for i in range( 1, len( props ) ):
    outMesh.write ( "property" + props[i] )
  outMesh.write ( "end_header\n" )

  # copy vertices
  for i in range( 1, len( argv ) - 1 ):
    inFile = file( argv[i], 'r' )
    inMesh = inFile.read()
    inFile.close()
    lines = (inMesh.split( "end_header" )[1]).split( "\n" )
    for j in range( 1, vertexNums[i-1] + 1 ):
      outMesh.write( lines[j]+"\n" )

  # copy faces
  vertexCount = 0
  for i in range( 1, len( argv ) - 1 ):
    inFile = file( argv[i], 'r' )
    inMesh = inFile.read()
    inFile.close()
    lines = (inMesh.split( "end_header" )[1]).split( "\n" )
    for j in range( vertexNums[i-1] + 1, len( lines ) ):
      words = lines[j].split()
      if len( words ) > 1:
        outMesh.write( words[0] )
        for k in range( 1, len( words ) ):
          outMesh.write( " " + str( int( words[k] ) + vertexCount ) )
        outMesh.write( "\n" )
    vertexCount = vertexCount + vertexNums[i-1]

  outMesh.flush()
