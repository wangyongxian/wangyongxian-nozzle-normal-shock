set grid
set time
set data style linespoints
set key right bottom
set xlabel 'log (h)'
set ylabel 'logaritmo do módulo de erros verdadeiros'
  plot 'E_pL_p2.dat' using 1:2 title 'E(pL)'
replot 'E_pL_p2.dat' using 1:3 title 'E(p2)'
replot 'E_p3.dat'  using 1:2 title 'E(p3)'
replot 'Eh.dat'     using 1:2 title 'E(Th)'
replot 'Ei_12.dat'  using 1:2 title 'E(Ti_12)'
replot 'Ei_pU.dat'  using 1:2 title 'E(Ti_pU)'
replot 'Ec_12.dat'  using 1:2 title 'E(Tc_12)'
replot 'Ec_13.dat'  using 1:2 title 'E(Tc_13)'
replot 'Ei_13.dat'  using 1:2 title 'E(Ti_13)'
replot 'Ei_23.dat'  using 1:2 title 'E(Ti_23)'
replot 'Ei_bi.dat'  using 1:2 title 'E(Ti_bi)'
replot 'Ei_tri.dat' using 1:2 title 'E(Ti_tri)'
replot 'Ebi.dat'    using 1:2 title 'E(Tbi_pU)'
set title 'T(1/2), Peclet 1Dp 1.3, CDS, S(0), Pe=10, Richardson 3.2 (10/7/2009)'
replot
