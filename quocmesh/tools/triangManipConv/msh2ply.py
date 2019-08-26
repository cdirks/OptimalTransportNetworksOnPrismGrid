#!/usr/bin/python
# -*- coding: utf-8 -*-

# converts a .msh-file into a list of .ply-files, one for each mesh

from sys import *
from os  import *

if len( argv ) != 3:
  print r"converts a .msh-file into a list of .ply-files, one for each mesh"
  print r"usage: "+argv[0]+" <input file> <output file basename>\n"

else:
  inFile = file( argv[1], 'r' )
  inData = inFile.read()
  inMeshes = inData.split("Mesh: ")

  for i in range( 1, len( inMeshes ) ):
    meshParts = inMeshes[i].split ( "_Count: " )
    vertexLines = meshParts[1].split( "\n" )
    triangLines = meshParts[3].split( "\n" )

    outMesh = file( argv[2]+str(i)+".ply", 'w' )
    outMesh.write ( "ply\n" )
    outMesh.write ( "format ascii 1.0\n" )
    outMesh.write ( "comment: converted from "+argv[1]+" via QuocMesh msh2ply\n" )
    outMesh.write ( "element vertex "+vertexLines[0].split(" ")[0]+"\n" )
    outMesh.write ( "property float x\n" )
    outMesh.write ( "property float y\n" )
    outMesh.write ( "property float z\n" )
    outMesh.write ( "property float I\n" )
    outMesh.write ( "element face "+triangLines[0]+"\n" )
    outMesh.write ( "property list uchar int vertex_index\n" )
    outMesh.write ( "end_header\n" )

    for l in range ( 1, len(vertexLines) - 1 ):
      outMesh.write( vertexLines[l]+"\n" )

    for l in range ( 1, len(triangLines) ):
      if triangLines[l].strip():
        outMesh.write( "3 "+triangLines[l]+"\n" )

    outMesh.flush()
