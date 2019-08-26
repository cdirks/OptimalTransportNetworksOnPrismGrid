# This script produces a tex-file, which itself encodes a ps or pdf, in which
# four columns of images are displayed. The i-th row corresponds to the i-th image.
# The first column depicts the image itself, the second to fourth column the
# length-, volume- and surface-deformation energy.
# The number of images as well as the wordstem of the actual image have to be
# passed by the command line.

import os
import sys

if len(sys.argv) != 2:
  print r"Correct usage would be: python geodesicTex.py <parameter file name>"
else:
  # read in parameter file
  f = file(sys.argv[1], 'r')
  text = f.read()
  f.close()
  wordlist = text.split()
  directory  = wordlist[wordlist.index("destDirectory")+1]
  iterations = int(wordlist[wordlist.index("maxIterations")+1])

  numberList = [int(2),int(3),int(5),int(9),int(17),int(25),int(33)]
  for index in range(1,len(numberList)):
    numberOfImages = numberList[index]
    # compute geodesic with numberOfImages points
    f = file(sys.argv[1], 'r')
    text = f.read()
    f.close()
    text = text.replace(r"numberOfPoints",r"numberOfPoints "+str(numberOfImages)+r" #")
    text = text.replace(r"InitNumberOfShapes",r"InitNumberOfShapes "+str(numberList[index-1])+r" #")
    f = file(sys.argv[1], 'w')
    f.write(text)
    f.close()
    os.system(r"averagePhaseField "+sys.argv[1])

    # compute geodesic length
    dataFile = file(directory+"/energyValuesOnLevel7.txt",'r')
    data = dataFile.read()
    dataFile.close()
    dataList = data.split()
    geodesicLength = float(0)
    for image in range(numberOfImages):
      geodesicLength = geodesicLength + float(dataList[image*(iterations+1)+iterations])
    geodesicLength = geodesicLength*(numberOfImages-float(1))

    # add geodesic length to plot file
    if not os.path.exists(directory+r"/geodesicLength.txt"):
      f = file(directory+r"/geodesicLength.txt", 'w')
      f.close()
    f = file(directory+r"/geodesicLength.txt", 'r')
    text = f.read()
    f.close()
    text = text+str(numberOfImages)+" "+str(geodesicLength)+"\n"
    f = file(directory+r"/geodesicLength.txt", 'w')
    f.write(text)
    f.close()

    # plot geodesic length
    # produce temporary gnuplot command file
    f = file(r"tmpFile.tmpFile",'w')
    f.write("set title \"geodesic length\"\n")
    f.write("set terminal png\n")
    f.write("\0")
    f.write("set xlabel \"points on discrete geodesic\"\n")
    f.write("set ylabel \"geodesic length\"\n")
    f.write("set output \""+directory+"/geodesicLength.png\"\n")
    f.write("plot \""+directory+"/geodesicLength.txt\" w l")
    f.close()
    # produce graph
    os.system(r"gnuplot tmpFile.tmpFile")
    os.system(r"rm tmpFile.tmpFile")
