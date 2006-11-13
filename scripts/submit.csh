#!/bin/tcsh
set opt="-j y -m n"
foreach pos (0.01 0.03 0.1 0.3 1 3 10 30)
  foreach neg (0.1 0.3 1 3 10 30 100 300 1000)
    foreach fdr (0.03)
      set vars="-v CPOS=$pos -v CNEG=$neg -v FDR=$fdr -o $cwd/out-$fdr-$pos-$neg"
      set cmd="qsub $opt $vars -N l$fdr-$pos-$neg /gs/home/lukall/cvs/percolator/scripts/wrap.csh"
      echo $cmd
      ssh -A clic1 $cmd
    end
  end
end
