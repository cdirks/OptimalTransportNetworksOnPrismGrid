# plot SSD
# Author:	Jingfeng Han

set time
set title "plot SSD_{REG}"
set style data linespoints
plot "./results/Diff_R.txt" using 1:2 title "SSD ||R-T_u||", \
     "./results/Diff_T.txt" using 1:2 title "SSD ||R_u-T||"