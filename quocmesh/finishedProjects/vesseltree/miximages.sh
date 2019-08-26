#!/bin/tcsh
foreach i($*)
    ls $i
    convert -crop 1024x1024 $i.eps $i:r:r.ve.ppm
    convert -geometry 50% $i:r:r.ve.ppm $i:r.ves.ppm
    convert -crop 128x128+1+1 $i.ppm $i:r:r.ti.ppm
    convert -geometry 400% $i:r:r.ti.ppm $i:r:r.tis.ppm
    pnmcomp -alpha /tmp/schwen/alpha_1024_full.pgm $i:r:r.ves.ppm $i:r:r.tis.ppm $i:r:r.mix.ppm
end
