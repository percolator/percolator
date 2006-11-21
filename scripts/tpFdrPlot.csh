#!/bin/tcsh
percolator -r percolator.res --gist-out gist -o per.sqt -s shuf.sqt $1 $2 $3
bestScore.py per.sqt shuf.sqt >! percolator_rerank.res
pepproph.py gist.data gist.label xx >! prohet_tryp.res
pepproph.py gist.data gist.label >! prohet.res
washburn.py gist.data gist.label >! washburn.pnt
dta.py gist.data gist.label xx >! DTAselect_tryp.pnt
dta.py gist.data gist.label >! DTAselect.pnt
echo "0 0 0 1 0" >! xcorr.w
echo "0 0 0 1 0 0 0 0 0 0 0 100 100 0" >! xcorr_tryp.w
score.py gist.data gist.label xcorr.w >! xcorr.res
score.py gist.data gist.label xcorr_tryp.w >! xcorr_tryp.res
plotROCs.py
