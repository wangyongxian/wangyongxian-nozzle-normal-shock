set data style linespoints
set grid
plot 'fm.dat' using 2:3 title 'fluxo de massa'
set xlabel 'x (m)'
set ylabel 'fluxo de massa (kg/s)'
set title 'PROG10_CFD1, MÉTODO SIMPLEC, N=12 (5/10/2010)'                
replot
