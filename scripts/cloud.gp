#!/usr/bin/gnuplot
set terminal postscript eps color 22
#set terminal postscript landscape color 22
#set terminal postscript portrait color 22
set output "cloud.ps"
set title 'Raw scores'
set xlabel "deltaCN"
set ylabel "Xcorr"
set arrow from 0,2.4116 to 0.9,2.4116 nohead lt 3
set arrow from 0,3.816 to 0.9,0.365 nohead lt 4
plot 'random.out' t 'Random' w points pt 4,'forward.out' t 'Forward' w points pt 1
