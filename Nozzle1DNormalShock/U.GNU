set data style linespoints
set grid
plot                                            'u.dat' using 1:2 title 'analítica 2'
replot                                            'u.dat' using 1:3 title 'numérica'
replot                                            'u.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel   'u(m/s)'
set title                                     'teste titulo'
replot
