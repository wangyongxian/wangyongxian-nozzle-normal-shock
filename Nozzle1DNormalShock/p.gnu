set data style linespoints
set grid
plot 'p.dat' using 1:2 title 'anal�tica'
replot 'p.dat' using 1:3 title 'num�rica'
replot 'p.dat' using 1:4 title 'geometria'
set xlabel 'x (m)'
set ylabel 'p (Pa)'
set title 'Resultado para N='                              
replot
