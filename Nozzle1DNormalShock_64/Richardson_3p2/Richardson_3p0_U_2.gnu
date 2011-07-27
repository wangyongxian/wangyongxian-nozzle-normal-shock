set grid
set time
set data style linespoints
set key right bottom
set xlabel 'log (h)'
set ylabel 'logaritmo do módulo de erros estimados'
  plot 'U_pL.dat'   using 1:2 title 'E(pL)'
replot 'U_p2.dat'   using 1:2 title 'E(p2)'
replot 'U_p3.dat'   using 1:2 title 'E(p3)'
replot 'U_12.dat'   using 1:3 title 'Uri_12(Th)'
replot 'U_12.dat'   using 1:4 title 'GCI_pL(Th)'
replot 'Uri_pU.dat' using 1:3 title 'Uri_p(Th)'
replot 'Uri_pU.dat' using 1:4 title 'GCI_p(Th)'
replot 'Uc_12.dat'  using 1:2 title 'Uc(Tc_12)'
replot 'Uc_13.dat'  using 1:2 title 'Uc(Tc_13)'
set title 'Resultado para N=12 (27/7/2011)'
replot
