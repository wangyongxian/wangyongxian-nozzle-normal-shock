set data style linespoints
set grid
plot 'mach.dat' using 1:2 title 'analítica' 
replot 'mach.dat' using 1:3 title 'numérica' 
replot 'mach.dat' using 1:4 title 'geometria'
set xlabel 'x (m)'
set ylabel 'Mach'
set title 'Resultado para N='                              
replot
