# plot E_{REG}, E_{EXT}, E_{INT}, E_{CON}
# Author:	Jingfeng Han

set time
set title "plot E_{REG} R2T"
set style data linespoints
plot "./results/EnergyLog_R2T.txt" using 1:2 title "E_{REG}", \
     "./results/EnergyLog_R2T.txt" using 1:3 title "E_{EXT}", \
     "./results/EnergyLog_R2T.txt" using 1:4 title "E_{INT}", \
     "./results/EnergyLog_R2T.txt" using 1:5 title "E_{CON}"