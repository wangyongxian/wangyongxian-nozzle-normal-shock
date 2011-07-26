set data style linespoints
set grid
plot                                            'T.dat' using 1:2 title 'analítica 2'
replot                                            'T.dat' using 1:3 title 'numérica'
replot                                            'T.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel     'T(K)'
set title                                     'teste titulo'
replot
