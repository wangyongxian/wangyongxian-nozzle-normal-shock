set data style linespoints
set grid
plot 'u.dat' using 2:3 title 'velocidades' 
set xlabel 'x (m)'
set ylabel 'u (m/s)'
set title 'PROG10_CFD1, MÉTODO SIMPLEC, N=12 (5/10/2010)'                
replot
