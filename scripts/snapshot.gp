set term png color
set size 0.75

set out '1.png'
set title 'Reynolds number'
set xlabel 'time'
set ylabel 'Re_lambda'
pl 'stat1.gp' u 2:7 w lp

set out '2.png'
set title 'eta*kmax'
set ylabel 'eta * kmax'
pl 'stat2.gp' u 2:7 w lp


