#!/bin/tcsh
percolator -r Percolator.res -w per.w --gist-out gist $1 $2 $3
percolator -m 5 -r 'Percolator (re-ranked).res' --gist-out m5_gist -o rerank.norm.sqt -s rerank.shuffled.sqt $1 $2 $3
pepproph.py gist.data gist.label Y >! 'PeptideProphet (tryptic).res'
pepproph.py gist.data gist.label >! 'PeptideProphet.res'
washburn.py gist.data gist.label >! Washburn.pnt
dta.py gist.data gist.label Y >! 'DTASelect (tryptic).pnt'
dta.py gist.data gist.label >! DTASelect.pnt
echo "0 0 0 1 0" >! xcorr.w
score.py gist.data gist.label xcorr.w >! SEQUEST.res
score.py gist.data gist.label xcorr.w  Y >! 'SEQUEST (tryptic).res'
plotROCs.py
mv roc.eps main_tp_fdr.eps

foreach n (1 2 3)
  mkdir $n
  score.py gist.data gist.label per.w N $n >! $n/Percolator.res
  sqtScore.pl rerank.norm.sqt rerank.shuffled.sqt $n >! $n/'Percolator (re-ranked).res'
  pepproph.py gist.data gist.label Y $n >! $n/'PeptideProphet (tryptic).res'
  pepproph.py gist.data gist.label N $n >! $n/PeptideProphet.res
  washburn.py gist.data gist.label $n >! $n/Washburn.pnt
  dta.py gist.data gist.label Y $n >! $n/'DTASelect (tryptic).pnt'
  dta.py gist.data gist.label N $n >! $n/DTASelect.pnt
  score.py gist.data gist.label xcorr.w N $n >! $n/SEQUEST.res
  score.py gist.data gist.label xcorr.w Y $n >! $n/'SEQUEST (tryptic).res'
  cd $n
    plotROCs.py
    mv roc.eps ../charge"$n"_tp_fdr.eps
  cd ..
end
