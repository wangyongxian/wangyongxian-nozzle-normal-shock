set data style linespoints
set grid
plot 'fm.dat' using 1:2 title 'analítico'
replot 'fm.dat' using 1:3 title 'numérico'
set xlabel 'x (m)'
set ylabel 'coeficiente de descarga (kg/s)'
set title ''                              
replot
