set data style linespoints
set grid
plot 'fm.dat' using 1:2 title 'anal�tico'
replot 'fm.dat' using 1:3 title 'num�rico'
set xlabel 'x (m)'
set ylabel 'empuxo (kg/s)'
set title ''                              
replot
