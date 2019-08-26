# plot tau
# Author:	Jingfeng Han

set time
set title "plot tau"
set style data linespoints
plot "./results/TauLog_R2T.txt" using 1:2 title "SSD ||tau_R2T||", \
     "./results/TauLog_T2R.txt" using 1:2 title "SSD ||taz_T2R||"