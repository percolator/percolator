#!/bin/tcsh
#$ -S /bin/tcsh
cd /gs/home/lukall/cvs/percolator/scan
/gs/home/lukall/cvs/percolator/src/Debug/percolator -F $FDR -p $CPOS -n $CNEG -r $FDR-$CPOS-$CNEG.res /nfs/noble2/lukall/percolator/normal-01.sqt /nfs/noble2/lukall/percolator/shuffled-01.sqt /nfs/noble2/lukall/percolator/shuffled2-01.sqt
