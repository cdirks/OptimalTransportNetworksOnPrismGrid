#!/usr/bin/perl
#
# Converts an obj-file into a ply-file (both file formats representing triangular, surfaces in 3D).
#
# Usage: perl obj2ply sourceFile.obj destinationFile.ply
#
# authors:
# Stefan Paulus, paulus@igg.uni-bonn.de
# Institute of Geodesy and Geoinformation
# Professorship of Geodesy
#
# Benedikt Wirth
# Institute for Numerical Simulation
#
# 12.04.2011


use Fcntl; # for reading pointer settings


# how many files have to be processed
$ARGV = @ARGV;
#print   "$ARGV \n";

##do for the whole list of files in the argument list:
for ($i=0;$i<$ARGV;$i= $i+2)
{
      #read outputfilename & open file
      #perl [x] name [x] filename[0] filename[1] ...
      my $inputfilename = @ARGV[$i];
      my $outputname = @ARGV[$i];
      my $j= $i+1;
      $outputname = @ARGV[$j];
      open(IN, "<$inputfilename") or die $! ; # open for reading
      open(PLY, ">$outputname") or die $! ;   # open for writing

      # $| = 1; means write text directly $| = 0; means wait with this
      $| = 1;
      print "IN : $inputfilename ";
      print "OUT: $outputname ";
      $| = 0;
      my $cou;            # read variable for lines
      my $amountVertices = 0; # number of vertices in file
      my $amountFaces    = 0; # number of faces in file

      # Amount of Vertices and Faces calculation
      # read first, because this must be written to the output header
      while (<IN>){
              $cou=$_;
              # count all lines with vertices
              if ((index($cou,"v") == 0)) {
                      $amountVertices = $amountVertices +1;
              }
              # count all lines with faces
              if ((index($cou,"f") == 0)) {
                      $amountFaces = $amountFaces +1;
              }
      }
      #print $amountVertices." vertices counted! \n" ;
      #print $amountFaces." faces  counted!" ;

      # reading pointer set to zero
      # we have to copy from the beginning
      $adr = 0;
      seek(IN,$adr,0);

      #write header of output file
      print PLY "ply\n";
      print PLY "format ascii 1.0\n";
      print PLY "comment author: Stefan Paulus\n";
      print PLY "comment IGG institut for geodesy and geoinformation\n";
      print PLY "comment university of bonn\n";
      print PLY "element vertex $amountVertices\n";
      print PLY "property float x\n";
      print PLY "property float y\n";
      print PLY "property float z\n";
      print PLY "element face $amountFaces\n";
      print PLY "property list uchar int vertex_indices\n";
      print PLY "end_header\n";
      #end of header writing

      # kills the process when order of vertices and faces is mixed!
      my $orderBit = 0;
      # counts all the written lines
      my $counter = 0;
      # sign for error
      my $error = 0;

      # read Input file, process and write to output file
      while (<IN>){
              $cou=$_;
              $counter++;
              ## no vertice, no face, not interesting!
              if ((index($cou,"g") == 0) || (index($cou,"#") == 0))  {
                 next;
              }
              ## copy all lines with vertices
              if ((index($cou,"v") == 0) && ($orderBit == 0)) {
                      @parts = split(/ /, $cou);

                      print PLY "@parts[1] ";
                      print PLY "@parts[2] ";
                      print PLY "@parts[3]";
              }

              ## copy all lines with faces
              if ((index($cou,"f") == 0)) {
                      $orderBit = 1;
                      @parts = split(/ /, $cou);

                      ## THE  BIG TRICK! -> OBJ starts the indexes with 1 to max,
                      ## but PLY with 0 to max-1
                      $first  = @parts[1]-1;# + $headerOffSet;
                      $second = @parts[2]-1;# + $headerOffSet;
                      $third  = @parts[3]-1;# + $headerOffSet;

                      print PLY "3 ";
                      print PLY "$first ";
                      print PLY "$second ";
                      print PLY "$third\n";

                      if (($first >= $amountVertices)||($first<0)){
                              print "Error : Index out of Range 1! LINE $counter. $first \n";
                              $error++;
                      }
                      if (($second >= $amountVertices)||($second<0)){
                              print "Error : Index out of Range 2! LINE $counter. $second \n";
                              $error++;
                      }
                      if (($third >= $amountVertices)||($third <0)){
                              print "Error : Index out of Range 3! LINE $counter. $third \n";
                              $error++;
                      }
                }
                ## close everything, throw exception and die,
                # if order of vertices and faces is damaged!
                if (((index($cou,"v") == 0)&&($orderBit != 0))||(error >0)) {
                     print "order of input file damaged!";
                     close(IN);
                     close(PLY);
                     unlink($outputname);
                     die;
               }
      }
      close(IN);
      close(PLY);
      $| = 1;
      print "done! \n";
      $| = 0;
}