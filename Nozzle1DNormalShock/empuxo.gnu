set data style linespoints
set grid
plot                                           'fm.dat' using 1:2 title 'anal�tica 2'
replot                                           'fm.dat' using 1:3 title 'num�rica'
replot                                           'fm.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel 'Empuxo (k
set title                                     'teste titulo'
replot
