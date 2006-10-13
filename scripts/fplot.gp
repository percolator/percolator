#!/usr/bin/gnuplot
#set terminal postscript landscape color 22
set terminal postscript eps color 22
set output "fplot.ps"
set title ''
set xlabel "FDR"
set ylabel "# true positives"
plot [0:0.15] 'fplot.out' u 2:1 t 'Only XCorr' w lines lt 1,'fplotItr.out' u 2:1 t 'LMS LD' w lines lt 2
