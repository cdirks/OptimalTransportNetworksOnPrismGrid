# This script produces a tex-file, which itself encodes a ps or pdf, in which
# four columns of images are displayed. The i-th row corresponds to the i-th image.
# The first column depicts the image itself, the second to fourth column the
# length-, volume- and surface-deformation energy.
# The number of images as well as the wordstem of the actual image have to be
# passed by the command line.

import os
import sys

if len(sys.argv) != 3:
  print r"Correct usage would be: python geodesicTex.py <parameter file name> <number of levelsets>"
else:
  # read in parameter file
  f = file(sys.argv[1], 'r')
  text = f.read()
  f.close()
  numLevelsets = int(sys.argv[2]);

  # read in parameters
  wordlist = text.split()
  coarsestLevel  = int(wordlist[wordlist.index("coarsestLevel")+1])
  finestLevel    = int(wordlist[wordlist.index("GridDepth")+1])
  numberOfImages = int(wordlist[wordlist.index("numberOfPoints")+1])
  iterations     = int(wordlist[wordlist.index("maxIterations")+1])
  gamma          = float(wordlist[wordlist.index("gamma")+1])
  directory      = wordlist[wordlist.index("destDirectory")+1]

  showFirstEnergy = False

  # produce energy graphs
  for level in range( coarsestLevel, finestLevel+1 ):
    elasticEnergy  = [ float(0) ]
    mismatchEnergy = [ float(0) ]
    for iteration in range(iterations):
      elasticEnergy  = elasticEnergy + [ float(0) ]
      mismatchEnergy = mismatchEnergy + [ float(0) ]

    dataFile = file(directory+"/energyValuesOnLevel"+str(level)+".txt",'r')
    data = dataFile.read()
    dataFile.close()
    dataList = data.split()
    pos = int(2)

    # produce deformation energy plots
    # produce temporary gnuplot data file
    for image in range(numberOfImages):
      dat = file("tmpFile."+str(image),"w")
      for iteration in range(iterations+1):
        if ( iteration != int(0) ) | ( showFirstEnergy ):
          dat.write(str(iteration)+" "+dataList[pos]+"\n")
        elasticEnergy[iteration] = elasticEnergy[iteration] + float(dataList[pos])
        pos = pos + 1
      dat.close()
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"hyperelastic energy\"\n")
    f.write("set output \""+directory+"/elasticEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile."+str(0)+"\" w l")
    for image in range(1,numberOfImages-1):
      f.write(", \""+"tmpFile."+str(image)+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    for image in range(numberOfImages):
      os.system(r"rm tmpFile."+str(image))

    # produce level set energy plots
    # produce temporary gnuplot data file
    pos = pos + 2
    for image in range(numLevelsets*numberOfImages):
      dat = file("tmpFile."+str(image),"w")
      for iteration in range(iterations+1):
        if ( iteration != int(0) ) | ( showFirstEnergy ):
          dat.write(str(iteration)+" "+dataList[pos]+"\n")
        mismatchEnergy[iteration] = mismatchEnergy[iteration] + ( gamma * float(dataList[pos]) )
        pos = pos + 1
      dat.close()
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"mismatch energy\"\n")
    f.write("set output \""+directory+"/mismatchEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile."+str(0)+"\" w l")
    for k in range(numLevelsets):
      for image in range(k*numLevelsets+1,k*numLevelsets+numberOfImages-1):
        f.write(", \""+"tmpFile."+str(image)+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    for image in range(numLevelsets*numberOfImages):
      os.system(r"rm tmpFile."+str(image))

    # produce total energy plots
    # produce temporary gnuplot data file
    pos = pos + 2
    dat = file("tmpFile.0","w")
    for iteration in range(iterations+1):
      if ( iteration != int(0) ) | ( showFirstEnergy ):
        dat.write(str(iteration)+" "+dataList[pos]+"\n")
      pos = pos + 1
    dat.close()
    dat = file("tmpFile.elastic","w")
    for iteration in range(iterations+1):
      if ( iteration != int(0) ) | ( showFirstEnergy ):
        dat.write(str(iteration)+" "+str(elasticEnergy[iteration])+"\n")
    dat.close()
    dat = file("tmpFile.mismatch","w")
    for iteration in range(iterations+1):
      if ( iteration != int(0) ) | ( showFirstEnergy ):
        dat.write(str(iteration)+" "+str(mismatchEnergy[iteration])+"\n")
    dat.close()
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"total energy\"\n")
    f.write("set output \""+directory+"/allEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile."+str(0)+"\" w l")
    f.write(", \""+"tmpFile.elastic"+"\" w l")
    f.write(", \""+"tmpFile.mismatch"+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"total energy\"\n")
    f.write("set output \""+directory+"/totalEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile."+str(0)+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"hyperelastic energy\"\n")
    f.write("set output \""+directory+"/totalElasticEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile.elastic"+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"energy decrease\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"iteration\"\n")
    f.write("set ylabel \"mismatch energy\"\n")
    f.write("set output \""+directory+"/totalMismatchEnergy"+str(level)+".png\"\n")
    f.write("plot \""+"tmpFile.mismatch"+"\" w l")
    f.close()
    # produce energy graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
    os.system(r"rm tmpFile.0")
    os.system(r"rm tmpFile.elastic")
    os.system(r"rm tmpFile.mismatch")

  # open new tex-file
  texFile = file(directory+r"/geodesic.tex", 'w')

  # write header
  texFile.write(r"\documentclass[a4paper,12pt]{article}"+"\n")
  texFile.write(r"\usepackage{graphicx}"+"\n")
  texFile.write(r"%\usepackage{amsmath}"+"\n\n")

  # change layout
  texFile.write(r"\setlength{\textwidth}{19cm}"+"\n")
  texFile.write(r"\setlength{\textheight}{29cm}"+"\n")
  texFile.write(r"\setlength{\oddsidemargin}{-2.3cm}"+"\n")
  texFile.write(r"\setlength{\evensidemargin}{-2.3cm}"+"\n")
  texFile.write(r"\setlength{\topmargin}{-2.3cm}"+"\n")
  texFile.write(r"\setlength{\linewidth}{19cm}"+"\n")
  texFile.write(r"\setlength{\headheight}{0cm}"+"\n")
  texFile.write(r"\setlength{\headsep}{0cm}"+"\n")
  texFile.write(r"\setlength{\parindent}{0cm}"+"\n")
  texFile.write(r"\setlength{\parskip}{0.1cm}"+"\n")
  texFile.write(r"\setlength{\linewidth}{19cm}"+"\n\n")

  # write tex-body
  texFile.write(r"\begin{document}"+"\n")

  # include images of the geodesics, level sets, deformations, and energies
  for level in range( coarsestLevel, finestLevel+1 ):

    # depict complete geodesic
    for levelset in range( numLevelsets ):
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for i in range(numberOfImages):
        ind = level * 10000 + iterations * 100 + i
        texFile.write(r"  \includegraphics{"+directory+r"/object"+str(levelset)+"_"+str(ind)+r".png}%"+"\n")
      texFile.write(r"}"+"\n\n")
      # depict pullback of shapes for complete geodesic
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for i in range(1,numberOfImages):
        ind = level * 10000 + iterations * 100 + i
        texFile.write(r"  \includegraphics{"+directory+r"/backDeformed"+str(levelset)+"Image"+str(ind)+r".png}%"+"\n")
      ind = level * 10000 + iterations * 100 + numberOfImages-1
      texFile.write(r"  \includegraphics{"+directory+r"/object"+str(levelset)+"_"+str(ind)+r".png}%"+"\n")
      texFile.write(r"}"+"\n\n")

      # depict iterates of the objects and level set functions
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for iteration in range(iterations+1):
        ind = level * 10000 + iteration * 100 + 0
        texFile.write(r"  \includegraphics{"+directory+r"/mismatch"+str(levelset)+"Energy"+str(ind)+r".png}%"+"\n")
      texFile.write(r"}\\"+"\n")
      for i in range(1,numberOfImages-1):
        # show evolution of object
        texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
        for iteration in range(iterations+1):
          ind = level * 10000 + iteration * 100 + i
          texFile.write(r"  \includegraphics{"+directory+r"/object"+str(levelset)+"_"+str(ind)+r".png}%"+"\n")
        texFile.write(r"}\\"+"\n")
        # show evolution of pullbacks of objects
        texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
        for iteration in range(iterations+1):
          ind = level * 10000 + iteration * 100 + i
          texFile.write(r"  \includegraphics{"+directory+r"/backDeformed"+str(levelset)+"Image"+str(ind)+r".png}%"+"\n")
        texFile.write(r"}\\"+"\n")
        # show evolution of level set function and/or mismatch energy
        texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
        for iteration in range(iterations+1):
          ind = level * 10000 + iteration * 100 + i
          texFile.write(r"  \includegraphics{"+directory+r"/mismatch"+str(levelset)+"Energy"+str(ind)+r".png}%"+"\n")
          # texFile.write(r"  \includegraphics{"+directory+r"/levelsetFunc"+str(ind)+r".png}%"+"\n")
        texFile.write(r"}\\"+"\n")
        # show derivative wrt level set function
        texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
        for iteration in range(iterations+1):
          ind = level * 10000 + iteration * 100 + i
          texFile.write(r"  \includegraphics{"+directory+r"/levelset"+str(levelset)+"Grad"+str(ind)+r".png}%"+"\n")
        texFile.write(r"}"+"\n\n")

      texFile.write(r"\newpage"+"\n\n")

    # depict iterates of the deformations and the local energy
    for i in range(numberOfImages-1):
      # show evolution of object
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for iteration in range(iterations+1):
        ind = level * 10000 + iteration * 100 + i
        texFile.write(r"  \includegraphics{"+directory+r"/dispChessPat"+str(ind+1)+r".png}%"+"\n")
      texFile.write(r"}\\"+"\n")
      # show evolution of level set function
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for iteration in range(iterations+1):
        ind = level * 10000 + iteration * 100 + i
        texFile.write(r"  \includegraphics{"+directory+r"/elasticEnergy"+str(ind)+r".png}%"+"\n")
      texFile.write(r"}\\"+"\n")
      # show derivative wrt level set function
      texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
      for iteration in range(iterations+1):
        ind = level * 10000 + iteration * 100 + i
        texFile.write(r"  \includegraphics{"+directory+r"/displGrad"+str(ind)+r".png}%"+"\n")
      texFile.write(r"}"+"\n\n")

    texFile.write(r"\newpage"+"\n\n")

  # depict evolution of energy
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/mismatchEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}\\"+"\n")
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/elasticEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}\\"+"\n")
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/allEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}\\"+"\n")
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/totalEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}\\"+"\n")
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/totalElasticEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}\\"+"\n")
  texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
  for level in range( coarsestLevel, finestLevel+1 ):
    texFile.write(r"  \includegraphics{"+directory+r"/totalMismatchEnergy"+str(level)+r".png}%"+"\n")
  texFile.write(r"}"+"\n\n")

  # append parameter file
  f = file(sys.argv[1], 'r')
  text = f.read()
  f.close()
  texFile.write(r"\newpage"+"\n\n")
  texFile.write(r"\begin{verbatim}"+"\n")
  texFile.write(text)
  texFile.write("\n"+r"\end{verbatim}"+"\n\n")

  texFile.write(r"\end{document}"+"\n")

  # save and close the file
  texFile.close()

  # convert images
  #cmd = r"for FILE in "+directory+r"/*.pgm; do convert $FILE ${FILE%pgm}png; done"
  #print cmd
  #os.system(cmd)

  # produce pdf
  cmd = 'pdflatex %s/geodesic.tex' %(directory)
  print cmd
  os.system(cmd)
