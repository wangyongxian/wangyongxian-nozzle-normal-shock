set data style linespoints
set grid
plot 'u.dat' using 1:2 title 'numerica' 
replo 'u.dat' using 1:3 title 'analitica' 
set xlabel 'x (m)'
set ylabel 'u (m/s)'
set title 'Resultado para N=12 (18/4/2011)'                              
replot
