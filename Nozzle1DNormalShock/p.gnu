set data style linespoints
set grid
plot 'p.dat' using 2:3 title 'pressão'
set xlabel 'x (m)'
set ylabel 'p (Pa)'
set title 'Resultado para N=12 (18/4/2011)'                              
replot
