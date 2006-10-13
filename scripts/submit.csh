#!/bin/tcsh
set opt="-j y -m n"
foreach pos (-2 -1.5 -1 -0.5 0 0.5 1 1.5 2)
  foreach neg (-2 -1.5 -1 -0.5 0 0.5 1 1.5 2)
    set vars="-v CPOS=$pos -v CNEG=$neg -v PYTHONPATH=$PYTHONPATH"
    set cmd="qsub $opt $vars -N six$pos-$neg $cwd/wrap.csh"
    echo $cmd
    ssh -A clic1 $cmd
  end
end
