set data style linespoints
set grid
plot                                         'Mach.dat' using 1:2 title 'analítica 2'
replot                                         'Mach.dat' using 1:3 title 'numérica'
replot                                         'Mach.dat' using 1:4 title 'geometria'
set xlabel     'x(m)'
set ylabel     'Mach'
set title                                     'teste titulo'
replot
