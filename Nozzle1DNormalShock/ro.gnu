set data style linespoints
set grid
plot                                           'RO.dat' using 1:2 title 'anal�tica 2'
replot                                           'RO.dat' using 1:3 title 'num�rica'
replot                                           'RO.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel 'Ro(km/m^3
set title                                     'teste titulo'
replot
