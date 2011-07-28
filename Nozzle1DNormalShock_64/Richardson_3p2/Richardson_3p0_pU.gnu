set grid
set time
set data style linespoints
set key right top
set yrange[-2:8]
set xlabel 'log (h)'
set ylabel 'ordem aparente'
  plot 'pL.dat'     using 1:2 title 'pL'
replot 'p2.dat'     using 1:2 title 'pV(2)'
replot 'p3.dat'     using 1:2 title 'pV(3)'
replot 'pU_h.dat'   using 1:2 title 'pU(Th)'
replot 'pU_i.dat'   using 1:2 title 'pU(Ti_pU)'
replot 'pU_c12.dat' using 1:2 title 'pU(Tc_12)'
replot 'pU_c13.dat' using 1:2 title 'pU(Tc_13)'
replot 'pU_bi.dat'  using 1:2 title 'pU(Tbi_pU)'
set title 'Resultado para N=12 (28/7/2011)'
replot
