#!/usr/bin/gnuplot
#set terminal postscript landscape color 22
set terminal postscript portrait color 22
set output "cloudItr.ps"
set title 'LMS Linear Discrimination'
set xlabel "deltaCN"
set ylabel "Xcorr"
plot 'randomItr.out' t 'Random' w points pt 4,'forwardItr.out' t 'Forward' w points pt 1
