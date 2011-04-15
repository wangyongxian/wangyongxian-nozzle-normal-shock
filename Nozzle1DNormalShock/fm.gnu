set data style linespoints
set grid
plot 'fm.dat' using 1:2 title 'fluxo de massa'
replot 'fm.dat' using 1:3 title 'fluxo de massa analitico'
set xlabel 'x (m)'
set ylabel 'fluxo de massa (kg/s)'
set title ''                              
replot
