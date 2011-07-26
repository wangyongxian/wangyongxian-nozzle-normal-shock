set data style linespoints
set grid
plot                                            'P.dat' using 1:2 title 'analítica 2'
replot                                            'P.dat' using 1:3 title 'numérica'
replot                                            'P.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel    'P(Pa)'
set title                                     'teste titulo'
replot
