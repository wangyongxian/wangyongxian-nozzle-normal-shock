set data style linespoints
set grid
plot 'p.dat' using 2:3 title 'pressão'
set xlabel 'x (m)'
set ylabel 'p (Pa)'
set title 'PROG10_CFD1, MÉTODO SIMPLEC, N=12 (5/10/2010)'                
replot
