set data style linespoints
set grid
plot 'ro.dat' using 1:2 title 'anal�tica'
replot 'ro.dat' using 1:3 title 'num�rica'
replot 'ro.dat' using 1:4 title 'geometria'
set xlabel 'x (m)'
set ylabel 'ro (km/m^3)'
set title ''                              
replot
