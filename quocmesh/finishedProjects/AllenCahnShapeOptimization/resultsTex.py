# This script produces a tex-file, which itself encodes a ps or pdf, in which
# four columns of images are displayed. The i-th row corresponds to the i-th image.
# The first column depicts the image itself, the second to fourth column the
# length-, volume- and surface-deformation energy.
# The number of images as well as the wordstem of the actual image have to be
# passed by the command line.

import os
import sys

if len(sys.argv) != 2:
  print r"Correct usage would be: python resultTex.py <parameter file name>"
else:
  # read in parameter file
  f = file(sys.argv[1], 'r')
  text = f.read()
  f.close()

  # read in parameters
  wordlist = text.split()
  coarsestLevel  = int(wordlist[wordlist.index("coarsestLevel")+1])
  finestLevel    = int(wordlist[wordlist.index("GridDepth")+1])
  iterations     = int(wordlist[wordlist.index("maxRepetitions")+1])
  directory      = wordlist[wordlist.index("destDirectory")+1]

  # open new tex-file
  texFile = file(directory+r"/result.tex", 'w')

  # write header
  texFile.write(r"\documentclass[a4paper,12pt]{article}"+"\n")
  texFile.write(r"\usepackage{graphicx}"+"\n")
  texFile.write(r"%\usepackage{amsmath}"+"\n\n")

  # change layout
  texFile.write(r"\setlength{\textwidth}{19cm}"+"\n")
  texFile.write(r"\setlength{\textheight}{29cm}"+"\n")
  texFile.write(r"\setlength{\oddsidemargin}{-2cm}"+"\n")
  texFile.write(r"\setlength{\evensidemargin}{-2cm}"+"\n")
  texFile.write(r"\setlength{\topmargin}{-2cm}"+"\n")
  texFile.write(r"\setlength{\linewidth}{19cm}"+"\n")
  texFile.write(r"\setlength{\headheight}{0cm}"+"\n")
  texFile.write(r"\setlength{\headsep}{0cm}"+"\n")
  texFile.write(r"\setlength{\parindent}{0cm}"+"\n")
  texFile.write(r"\setlength{\parskip}{0.1cm}"+"\n")
  texFile.write(r"\setlength{\linewidth}{19cm}"+"\n\n")

  # write tex-body
  texFile.write(r"\begin{document}"+"\n")

  # include images of the iterations
  for level in range( coarsestLevel, finestLevel+1 ):

    # depict iterations on that level
    texFile.write(r"\resizebox{\linewidth}{!}{%"+"\n")
    for i in range(iterations):
      ind = level * 10 + i
      texFile.write(r"  \includegraphics[width=100pt]{"+directory+r"/phasefield"+str(ind)+r".png}\hspace{10pt}%"+"\n")
    texFile.write(r"}"+"\n\n")

  # depict deformed geometry
  if os.path.exists(directory+r"/deformation.png"):
    texFile.write(r"\reflectbox{\includegraphics[width=.2\linewidth,angle=180]{"+directory+r"/deformation.png}}%"+"\n\n")

  # append parameter file
  f = file(sys.argv[1], 'r')
  text = f.read()
  f.close()
  texFile.write(r"%\newpage"+"\n\n")
  texFile.write(r"{\small%"+"\n")
  texFile.write(r"\begin{verbatim}"+"\n")
  texFile.write(text)
  texFile.write("\n"+r"\end{verbatim}"+"\n")
  texFile.write(r"}"+"\n\n")

  texFile.write(r"\end{document}"+"\n")

  # save and close the file
  texFile.close()

  # convert images
  cmd = r"for FILE in "+directory+r"/phasefield*.pgm; do convert $FILE ${FILE%pgm}png; done"
  print cmd
  os.system(cmd)

  # produce pdf
  cmd = 'pdflatex %s/result.tex' %(directory)
  print cmd
  os.system(cmd)
